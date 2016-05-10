/* Copyright 2013 Eric Chang and Rutwik Parikh */

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

#include "./eba.h"

/* The number of vertices in the graph. */
extern int n;
solution_t *cur_solution;

bool modified = 1;

// Free every element of an ilist
void list_free(ilist *current);

// Unset all locks for all nodes in the path up to but not including "end" or
// until the end of the (NULL-terminated) list. Note that if you set end
// to be -1 then we just unset all nodes in the list.
void unget_locks(ilist *current, int end) {
  while (current != NULL && current->vertex != end) {
    cur_solution->owner[current->vertex] = -1;
    omp_unset_lock(&(cur_solution->lock[current->vertex]));
    current = current->next;
  }
}

// Try to acquire locks for all the nodes in the path. Return true on success
// or false on failure.
bool get_locks(ilist *path) {
  ilist *current = path;
  int procId = omp_get_thread_num();
  while (current != NULL) {
    int u = current->vertex;
    // Try to acquire the current vertex's lock
    if (omp_test_lock(&(cur_solution->lock[u]))) {
      if (cur_solution->removed[u]) {
        cur_solution->owner[u] = -1;
        omp_unset_lock(&(cur_solution->lock[u]));
        unget_locks(path, u);
        return false;
      }
      cur_solution->owner[u] = procId;
      current = current->next;
    } else {
      // Busy loop to try and acquire the lock if we're a higher processor
      while (procId > cur_solution->owner[u]) {
        if (omp_test_lock(&(cur_solution->lock[u]))) {
          cur_solution->owner[u] = procId;
          current = current->next;
        }
      }
      // If we're a lower processor, give up.
      if (procId < cur_solution->owner[u]) {
        unget_locks(path, u);
        return false;
      }
    }
  }
  return true;
}

// Build an alternating path backward from tail through even to an exposed node
ilist *find_head_exposed(int even, ilist *tail, int avoid) {
  // Make sure that we avoid the last element of the tail (or we'd form a loop)
  if (even == avoid) {
    return NULL;
  }
  // Check if the even node is exposed
  if (cur_solution->matching[even] == -1) {
    return tail;
  }
  // Since the current node is even and not exposed, it must be matched
  // with an odd node if it is in the path
  int a = cur_solution->matching[even];
  // Add the next piece of the path
  ilist *head = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  ilist *next1 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  ilist *next2 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  *next2 = ilist(a, tail);
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is an even i
    if (cur_solution->forest[a][i] && cur_solution->even[i] && even != i) {
      *next1 = ilist(i, next2);
      if ((head = find_head_exposed(i, next1, avoid)) != NULL) {
        return head;
      }
    }
  }
  // Return NULL if we don't find anything
  free(next1);
  free(next2);
  free(head);
  return NULL;
}

// Finds an alternating path from prev to an exposed node through even
ilist *find_tail_exposed(int prev, int even) {
  // Check if the even node is exposed
  if (cur_solution->matching[even] == -1) {
    ilist *last = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    ilist *tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    *last = ilist(even, NULL);
    *tail = ilist(prev, last);
    return tail;
  }

  // Since the current node is even and not exposed, it must be matched
  // with an odd node if it is in the path
  int a = cur_solution->matching[even];
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is an even i
    if (cur_solution->forest[a][i] && cur_solution->even[i] && even != i) {
      ilist *tail;
      if ((tail = find_tail_exposed(a, i)) != NULL) {
        ilist *tail1 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
        ilist *tail2 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
        *tail1 = ilist(even, tail);
        *tail2 = ilist(prev, tail1);
        return tail2;
      }
    }
  }
  // Return NULL if we don't find anything
  return NULL;
}

// Checks path to see if it's an augmenting path. Returns true if it is,
// false if it is not.
bool check_path(ilist *path) {
  // First vertex of path must be exposed
  if (cur_solution->matching[path->vertex] != -1) {
    return false;
  }
  path = path->next;
  // Every two elements of the path must be paired
  int old = -1;
  while (path->next != NULL) {
    old = path->vertex;
    path = path->next;
    // Every two elements must be in the matching together
    if (cur_solution->matching[old] != path->vertex ||
        cur_solution->matching[path->vertex] != old) {
       return false;
    }
    // If next is NULL, then the path is of even length.
    if (path->next == NULL) {
      return false;
    }
    path = path->next;
  }
  // Last vertex of path must be exposed
  if (cur_solution->matching[path->vertex] != -1) {
    return false;
  }
  return true;
}

// Create an alternating path from sections from an exposed node to u,
// u to v, and v to an exposed node
ilist *find_path_exposed(int u, int v) {
  ilist *tail = find_tail_exposed(u, v);
  if (tail == NULL) {
    return NULL;
  }
  // Find the exposed node in the tail and make sure it isn't the first node
  // of the head path
  ilist *path = tail;
  while (path->next != NULL) {
    path = path->next;
  }
  int last = path->vertex;
  path = find_head_exposed(u, tail, last);
  if (path == NULL) {
    list_free(path);
    return NULL;
  }
  if (!get_locks(path)) {
    list_free(path);
    return NULL;
  }
  // If the path is legal, return it
  if (check_path(path)) {
    return path;
  } else {
    unget_locks(path, -1);
    list_free(path);
    return NULL;
  }
}

// Finds a path to the target and set the last element to be first,
// later forming a cycle
ilist *find_cycle(int prev, int current, int target, ilist *first) {
  // Check if the current node is the target
  if (current == target) {
    ilist *tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    *tail = ilist(prev, first);
    return tail;
  }
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is current
    if (i != prev && cur_solution->forest[current][i]) {
      ilist *tail;
      if ((tail = find_cycle(current, i, target, first)) != NULL) {
        ilist *new_tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
        *new_tail = ilist(prev, tail);
        return new_tail;
      }
    }
  }
  // Return NULL if we don't find anything
  return NULL;
}

// Checks cycle to see if it is a blossom and if so returns the cycle starting
// from the root. Otherwise returns NULL.
ilist *find_root(ilist *cycle) {
  int old = cycle->vertex;
  cycle = cycle->next;
  // This loop must terminate since there must exist some vertex with
  // a matching not in the blossom, or not in a matching (the root)
  while (cur_solution->matching[cycle->vertex] == old ||
      cur_solution->matching[cycle->vertex] == cycle->next->vertex) {
    old = cycle->vertex;
    cycle = cycle->next;
  }
  int root = cycle->vertex;
  cycle = cycle->next;
  // Now we check if the blossom rooted at root is an odd-length alternating
  // cycle. Iterate a pair at a time and check that they're in the matching.
  // Since we know that cycle is a cycle and root is an element of it, this
  // loop must terminate one way or another.
  while (cycle->vertex != root) {
    old = cycle->vertex;
    cycle = cycle->next;
    if (cycle->vertex == root || cur_solution->matching[old] != cycle->vertex
        || cur_solution->matching[cycle->vertex] != old) {
      return NULL;
    }
    cycle = cycle->next;
  }
  return cycle;
}

// Finds an odd-lengthed alternating cycle. Note that all cycles are
// guaranteed to be odd-length. Returns null if no such cycle exists.
blossom *find_blossom(int a, int b) {
  // Create the last element of the cycle
  ilist *first = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  *first = ilist(a, NULL);
  // Find the path that begins and ends on the same node
  ilist *path;
  if ((path = find_cycle(a, b, a, first)) == NULL) {
    free(first);
    return NULL;
  }
  // Allocate memory for the blossom we know must exist
  blossom *bloss = reinterpret_cast<blossom *>(malloc(sizeof(blossom)));
  // Try to acquire all the locks, if we fail, then free the path and
  // return a "fake" blossom to know not to look for an augmenting path
  if (!get_locks(path->next)) {
    list_free(path);
    *bloss = blossom(-1, NULL);
    return bloss;
  }
  // Close the cycle by removing the first element of the path and replacing
  // it with the last element
  first->next = path->next;
  free(path);
  // Find the root of the cycle
  ilist *root_cycle = find_root(first);
  // Returns NULL if the blossom is invalid
  if (root_cycle == NULL) {
    // Unset the locks and free the cycle
    ilist *second = first->next;
    first->next = NULL;
    // -1 is a dummy value to unset locks for the whole list
    unget_locks(second, -1);
    free(first);
    list_free(second);
    *bloss = blossom(-1, NULL);
    return bloss;
  }
  *bloss = blossom(root_cycle->vertex, root_cycle);
  return bloss;
}

// Removes all non-root elements of a blossom's cycle from the graph and forest
void remove_blossom(blossom *b) {
  int r = b->vertex;
  ilist *current = b->first->next;
  // Iterate over the elements in the list, removing them from the graph,
  // the forest, and the matching
  #pragma omp critical
  {
    while (current->vertex != r) {
      int v = current->vertex;
      // Find all of v's neighbors and remove them. Keep track of what
      // edges we change in the blossom's modified array
      for (int i = 0; i < n; i++) {
        if (cur_solution->graph[v][i] && i != r) {
          // Remove v -> i and i -> v from the graph.
          // The 1 in modified means we need to add it back in later
          b->modified[v][i] = 1;
          b->modified[i][v] = 1;
          cur_solution->graph[v][i] = 0;
          cur_solution->graph[i][v] = 0;
          // Remove the edges of v from the forest as well
          if (cur_solution->forest[v][i]) {
            cur_solution->forest[r][i] = 1;
            cur_solution->forest[i][r] = 1;
            cur_solution->forest[v][i] = 0;
            cur_solution->forest[i][v] = 0;
          }
          // If there is no edge between r and i in the graph, we need to
          // add it in. The -1 in modified means we need to remove it later
          if (!(cur_solution->graph[r][i])) {
            b->modified[r][i] = -1;
            b->modified[i][r] = -1;
            cur_solution->graph[r][i] = 1;
            cur_solution->graph[i][r] = 1;
          }
        }
      }
      // Disconnect v from r if they were connected. The 1 in modified
      // means that we have to add the edge back later
      if (cur_solution->graph[v][r]) {
        b->modified[v][r] = 1;
        b->modified[r][v] = 1;
        cur_solution->graph[v][r] = 0;
        cur_solution->graph[r][v] = 0;
        cur_solution->forest[v][r] = 0;
        cur_solution->forest[r][v] = 0;
      }
      // Every node in the blossom except for the root must be in a matching
      cur_solution->matching[v] = -1;
      // The removed nodes are neither even nor odd
      cur_solution->even[v] = 0;
      cur_solution->odd[v] = 0;
      cur_solution->removed[v] = 1;
      current = current->next;
      cur_solution->owner[v] = -1;
      cur_solution->owner[v] = -1;
      omp_unset_lock(&(cur_solution->lock[v]));
    }
  }
  // If we remove the matched element of r, r isn't in a matching anymore
  if (cur_solution->removed[cur_solution->matching[r]]) {
    cur_solution->matching[r] = -1;
  }
    cur_solution->owner[r] = -1;
  cur_solution->owner[r] = -1;
  omp_unset_lock(&(cur_solution->lock[r]));
}

// Given a path, augment it. Returns true if the graph was fully augmented
bool augment(ilist *path) {
  int old = path->vertex;
  ilist *current = path->next;
  free(path);
  ilist *temp;
  // Add the first edge to the matching
  cur_solution->matching[old] = current->vertex;
  cur_solution->even[old] = 0;
  cur_solution->matching[current->vertex] = old;
  cur_solution->odd[current->vertex] = 0;
  cur_solution->owner[old] = -1;
  cur_solution->owner[current->vertex] = -1;
  omp_unset_lock(&(cur_solution->lock[old]));
  omp_unset_lock(&(cur_solution->lock[current->vertex]));
  // Iterate two elements at a time and set their matchings to each other
  while (current->next != NULL) {
    temp = current->next;
    free(current);
    current = temp;
    // Since paths are odd length, we never seg fault
    old = current->vertex;
    temp = current->next;
    free(current);
    current = temp;
    cur_solution->matching[old] = current->vertex;
    cur_solution->even[old] = 0;
    cur_solution->matching[current->vertex] = old;
    cur_solution->odd[current->vertex] = 0;
    cur_solution->owner[old] = -1;
    cur_solution->owner[current->vertex] = -1;
    omp_unset_lock(&(cur_solution->lock[old]));
    omp_unset_lock(&(cur_solution->lock[current->vertex]));
  }
  return true;
}

// Lift all the blossoms back into the graph and add half the edges in
// their cycles to the matching
void lift_blossom(blossom *b) {
  ilist *current = b->first;
  ilist *temp;
  // If the blossom is not in a matching (-1) then it doesn't matter what the
  // root is, just leave it as the old root.
  // However, if the blossom is in a matching, we have to identify
  // the new root that the matching is attached to.
  // If an edge was never added from the old root to what it is matched
  // to now (value of -1 in modified), then the old root is still the
  // root, but otherwise we have to look for the new root.
  if (b->modified[b->vertex][cur_solution->matching[b->vertex]] == -1) {
    // Since the old root (first element of the cycle) is not the new root,
    // keep looking
    current = current->next;
    // Iterate along the cycle until we hit the new root, which is any
    // vertex in the cycle that can be in the blossom's current matching,
    // meaning that an edge exists between the new root and the old root's
    // matched vertex
    while (b->modified[current->vertex][cur_solution->matching[b->vertex]]
        != 1) {
      current = current->next;
    }
    // Move the matched edge from the old root to the new root
    cur_solution->matching[current->vertex] =
        cur_solution->matching[b->vertex];
    cur_solution->matching[cur_solution->matching[b->vertex]] =
        current->vertex;
    cur_solution->matching[b->vertex] = -1;
  }
  int new_root = current->vertex;
  // This is guaranteed to exist since an odd-length cycle has at
  // least 3 nodes/edges
  current = current->next;
  // Iterate along the cycle 2 at a time until we hit the old root, which
  // is the end of the path
  while (current->vertex != new_root) {
    int old = current->vertex;
    temp = current->next;
    free(current);
    current = temp;
    cur_solution->matching[old] = current->vertex;
    cur_solution->matching[current->vertex] = old;
    temp = current->next;
    free(current);
    current = temp;
  }
  free(current);
  free(b);
}

/* Checks an edge u-v to see if it completes an alternating path.
 * Augments the path if found. If a blossom
 * is discovered, contract it. Otherwise add the edge to the forest.
 */
void find_alternating_path(int u, int v) {
  // Found an edge u-v, now case on it
  if (!(cur_solution->even[v] || cur_solution->odd[v])) {
    // If v is not in even, then it must be in a matching. Mark v
    // as odd, and mark its matched vertex as even. Then, add u-v
    // and v-M(v) to F
    int m_v = cur_solution->matching[v];
    // Try to get the locks for v and M(v). Lock the smaller of the two first
    // so we don't deadlock.
    int low = v;
    if (m_v < v) {
      low = m_v;
    }
    // If we fail, then we know someone else already added them to the forest.
    if (omp_test_lock(&(cur_solution->lock[low]))) {
      // If we get the lock after someone else, v will already be even or odd,
      // as will M(v)
      if (!(cur_solution->even[v] || cur_solution->odd[v])) {
        int m_v = cur_solution->matching[v];
        cur_solution->forest[u][v] = 1;
        cur_solution->forest[v][u] = 1;
        cur_solution->odd[v] = 1;
        cur_solution->even[m_v] = 1;
        cur_solution->forest[v][m_v] = 1;
        cur_solution->forest[m_v][v] = 1;
      }
      omp_unset_lock(&(cur_solution->lock[low]));
    }
  } else if (cur_solution->even[v]) {
    // We know v is in even, so now we case on whether u and v are
    // connected in the forest
    blossom *b;
    // If u and v are connected and v is even, we form a blossom.
    // Run DFS to find the odd-length cycle with alternating
    // membership in the matching.
    // If it doesn't return NULL, then every edge in the blossom must be
    // locked.
    if ((b = find_blossom(u, v)) != NULL) {
      // If the vertex is set to -1, then we failed to find a blossom due
      // to relinquishing locks. This is to avoid checking for an
      // alternating path when we just return NULL
      if (b->vertex == -1) {
        free(b);
      } else {
        // Create a new blist element and add it to the existing list
        blist *l = reinterpret_cast<blist *>(malloc(sizeof(blist)));
        // Note that we lock the stack before adding to it to avoid
        // concurrency issues. Doesn't slow us since shouldn't have too many
        // blossoms at once.
        omp_set_lock(&(cur_solution->blist_lock));
        *l = blist(b, cur_solution->blossoms);
        cur_solution->blossoms = l;
        omp_unset_lock(&(cur_solution->blist_lock));
        // Remove all vertices but the root of the blossom from the graph
        // and forest, and remove all vertices of the blossom from the
        // matching.
        // At the end, we unlock all the nodes in the blossom.
        remove_blossom(b);
      }
    } else {
      ilist *path = NULL;
      // Since u and v are not connected, there must be an augmenting
      // path through u and v. Run DFS to find a path through u and v
      // to two exposed nodes. Note that we have all the locks if
      // we don't return NULL.
      if ((path = find_path_exposed(u, v)) != NULL) {
        // Augment the path; should always work. Note that it also
        // unlocks all the nodes in the path.
        modified = augment(path);
      }
    }
  }
}

/*
 * Runs Edmonds' Blossom Algorithm parallelizing over even nodes u
 */
void solve_eba(solution_t *solution) {
  cur_solution = solution;
  modified = 1;
  // Do some basic matchings first
  for (int i = 0; i < n; i++) {
    if (solution->matching[i] == -1) {
      for (int j = 0; j < n; j++) {
        if (i != j && solution->graph[i][j] && solution->matching[j] == -1) {
          solution->matching[i] = j;
          solution->matching[j] = i;
          break;
        }
      }
    }
  }
  // Look for augmentable paths, break out when none remain
  while (modified) {
    // Place every exposed vertex in even
    for (int i = 0; i < n; i++) {
      // For all non-removed vertices, if it is not in the matching,
      // mark it as even
      if (!(solution->removed[i]) && solution->matching[i] == -1) {
        solution->even[i] = 1;
      } else {
        solution->even[i] = 0;
      }
    }
    memset(solution->odd, 0, sizeof(solution->odd));
    memset(solution->forest, 0, sizeof(solution->forest));
    // If modified is set, then solution has been modified
    modified = 0;
    // If flag is set, then we changed the forest
    bool flag = 1;
    // Keep going until we find a path, give up if there's no new data being
    // added (flag stays 0)
    while (flag && !modified) {
      flag = 0;
      // Try every pair of u and v
      #pragma omp parallel for schedule(dynamic)
      for (int u = 0; u < n; u++) {
        if (!modified && !(solution->removed[u]) && solution->even[u]) {
          // Iterate over non-odd neighbors of u
          for (int v = 0; v < n; v++) {
            // Check all edges from an even u to a not odd v
            if (!modified && solution->graph[u][v] && !(solution->odd[v])) {
              // Add u-v to the forest and see if it completes a blossom
              // or alternating path
              find_alternating_path(u, v);
              flag = 1;
            }
          }
        }
      }
    }
  }
  // Lift the contracted blossoms and add new edges to M
  blist *current = solution->blossoms;
  blist *temp;
  while (current != NULL) {
    lift_blossom(current->b);
    temp = current->next;
    free(current);
    current = temp;
  }
  return;
}
