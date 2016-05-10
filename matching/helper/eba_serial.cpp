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

solution_t *seq_solution;

extern bool modified;

// Free every element of an ilist
void list_free(ilist *current) {
  ilist *temp = current;
  while (current != NULL) {
    temp = current->next;
    free(current);
    current = temp;
  }
}

// Build an alternating path backward from tail through even to an exposed node
ilist *find_head_exposed_serial(int even, ilist *tail) {
  // Check if the even node is exposed
  if (seq_solution->matching[even] == -1) {
    return tail;
  }

  // Since the current node is even and not exposed, it must be matched
  // with an odd node if it is in the path
  int a = seq_solution->matching[even];
  // Add the next piece of the path
  ilist *head = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  ilist *next1 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  ilist *next2 = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  *next2 = ilist(a, tail);
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is an even i
    if (seq_solution->forest[a][i] && seq_solution->even[i] && even != i) {
      *next1 = ilist(i, next2);
      if ((head = find_head_exposed_serial(i, next1)) != NULL) {
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
ilist *find_tail_exposed_serial(int prev, int even) {
  // Check if the even node is exposed
  if (seq_solution->matching[even] == -1) {
    ilist *last = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    ilist *tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    *last = ilist(even, NULL);
    *tail = ilist(prev, last);
    return tail;
  }

  // Since the current node is even and not exposed, it must be matched
  // with an odd node if it is in the path
  int a = seq_solution->matching[even];
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is an even i
    if (seq_solution->forest[a][i] && seq_solution->even[i] && even != i) {
      ilist *tail;
      if ((tail = find_tail_exposed_serial(a, i)) != NULL) {
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

// Create an alternating path from sections from an exposed node to u,
// u to v, and v to an exposed node
ilist *find_path_exposed_serial(int u, int v) {
  ilist *tail = find_tail_exposed_serial(u, v);
  ilist *path = find_head_exposed_serial(u, tail);
  return path;
}

// Finds a path to the target and set the last element to be first,
// later forming a cycle
ilist *find_cycle_serial(int prev, int current, int target, ilist *first) {
  // Check if the current node is the target
  if (current == target) {
    ilist *tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
    *tail = ilist(prev, first);
    return tail;
  }
  // Iterate over the neighbors of a
  for (int i = 0; i < n; i++) {
    // Check this neighbor of a to see if it is current
    if (i != prev && seq_solution->forest[current][i]) {
      ilist *tail;
      if ((tail = find_cycle_serial(current, i, target, first)) != NULL) {
        ilist *new_tail = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
        *new_tail = ilist(prev, tail);
        return new_tail;
      }
    }
  }
  // Return NULL if we don't find anything
  return NULL;
}

// Finds an odd-lengthed alternating cycle. Note that all cycles are
// guaranteed to be odd-length. Returns null if no such cycle exists.
blossom *find_blossom_serial(int a, int b) {
  // Create the last element of the cycle
  ilist *first = reinterpret_cast<ilist *>(malloc(sizeof(ilist)));
  *first = ilist(a, NULL);
  // Find the path that begins and ends on the same node
  ilist *path;
  if ((path = find_cycle_serial(a, b, a, first)) == NULL) {
    free(first);
    return NULL;
  }
  // Close the cycle by removing the first element of the path and replacing
  // it with the last element
  first->next = path->next;
  free(path);
  // Find the root of the cycle
  int old = first->vertex;
  first = first->next;
  // This loop must terminate since there must exist some vertex with
  // a matching not in the blossom, or not in a matching
  while (seq_solution->matching[first->vertex] == old ||
      seq_solution->matching[first->vertex] == first->next->vertex) {
    old = first->vertex;
    first = first->next;
  }
  blossom *bloss = reinterpret_cast<blossom *>(malloc(sizeof(blossom)));
  *bloss = blossom(first->vertex, first);
  return bloss;
}

// Removes all non-root elements of a blossom's cycle from the graph and forest
void remove_blossom_serial(blossom *b) {
  int r = b->vertex;
  ilist *current = b->first->next;
  // Iterate over the elements in the list, removing them from the graph,
  // the forest, and the matching
  while (current->vertex != r) {
    int v = current->vertex;
    // Find all of v's neighbors and remove them. Keep track of what
    // edges we change in the blossom's modified array
    for (int i = 0; i < n; i++) {
      if (seq_solution->graph[v][i] && i != r) {
        // Remove v -> i and i -> v from the graph.
        // The 1 in modified means we need to add it back in later
        b->modified[v][i] = 1;
        b->modified[i][v] = 1;
        seq_solution->graph[v][i] = 0;
        seq_solution->graph[i][v] = 0;
        // Remove the edges of v from the forest as well
        if (seq_solution->forest[v][i]) {
          seq_solution->forest[r][i] = 1;
          seq_solution->forest[i][r] = 1;
          seq_solution->forest[v][i] = 0;
          seq_solution->forest[i][v] = 0;
        }
        // If there is no edge between r and i in the graph, we need to
        // add it in. The -1 in modified means we need to remove it later
        if (!(seq_solution->graph[r][i])) {
          b->modified[r][i] = -1;
          b->modified[i][r] = -1;
          seq_solution->graph[r][i] = 1;
          seq_solution->graph[i][r] = 1;
        }
      }
    }
    // Disconnect v from r if they were connected. The 1 in modified
    // means that we have to add the edge back later
    if (seq_solution->graph[v][r]) {
      b->modified[v][r] = 1;
      b->modified[r][v] = 1;
      seq_solution->graph[v][r] = 0;
      seq_solution->graph[r][v] = 0;
      seq_solution->forest[v][r] = 0;
      seq_solution->forest[r][v] = 0;
    }
    // Every node in the blossom except for the root must be in a matching
    seq_solution->matching[v] = -1;
    // The removed nodes are neither even nor odd
    seq_solution->even[v] = 0;
    seq_solution->odd[v] = 0;
    seq_solution->removed[v] = 1;
    current = current->next;
  }
  // If we remove the matched element of r, r isn't in a matching anymore
  if (seq_solution->removed[seq_solution->matching[r]]) {
    seq_solution->matching[r] = -1;
  }
}

// Given a path, augment it. Returns true if the graph was fully augmented
bool augment_serial(ilist *path) {
  int old = path->vertex;
  ilist *current = path->next;
  free(path);
  ilist *temp;
  // Add the first edge to the matching
  seq_solution->matching[old] = current->vertex;
  seq_solution->matching[current->vertex] = old;
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
    seq_solution->matching[old] = current->vertex;
    seq_solution->even[old] = 0;
    seq_solution->matching[current->vertex] = old;
    seq_solution->odd[current->vertex] = 0;
  }
  return true;
}

// Lift all the blossoms back into the graph and add half the edges in
// their cycles to the matching
void lift_blossom_serial(blossom *b) {
  ilist *current = b->first;
  ilist *temp;
  // If the blossom is not in a matching (-1) then it doesn't matter what the
  // root is, just leave it as the old root.
  // However, if the blossom is in a matching, we have to identify
  // the new root that the matching is attached to.
  // If an edge was never added from the old root to what it is matched
  // to now (value of -1 in modified), then the old root is still the
  // root, but otherwise we have to look for the new root.
  if (b->modified[b->vertex][seq_solution->matching[b->vertex]] == -1) {
    // Since the old root (first element of the cycle) is not the new root,
    // keep looking
    current = current->next;
    // Iterate along the cycle until we hit the new root, which is any
    // vertex in the cycle that can be in the blossom's current matching,
    // meaning that an edge exists between the new root and the old root's
    // matched vertex
    while (b->modified[current->vertex][seq_solution->matching[b->vertex]]
        != 1) {
      current = current->next;
    }
    // Move the matched edge from the old root to the new root
    seq_solution->matching[current->vertex] =
        seq_solution->matching[b->vertex];
    seq_solution->matching[seq_solution->matching[b->vertex]] =
        current->vertex;
    seq_solution->matching[b->vertex] = -1;
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
    seq_solution->matching[old] = current->vertex;
    seq_solution->matching[current->vertex] = old;
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
void find_alternating_path_serial(int u, int v) {
  // Found an edge u-v, now case on it
  if (!(seq_solution->even[v] || seq_solution->odd[v])) {
    int m_v = seq_solution->matching[v];
    seq_solution->forest[u][v] = 1;
    seq_solution->forest[v][u] = 1;
    seq_solution->odd[v] = 1;
    seq_solution->even[m_v] = 1;
    seq_solution->forest[v][m_v] = 1;
    seq_solution->forest[m_v][v] = 1;
  } else if (seq_solution->even[v]) {
    // We know v is in even, so now we case on whether u and v are
    // connected in the forest
    blossom *b;
    // If u and v are connected and v is even, we form a blossom.
    // Run DFS to find the odd-length cycle with alternating
    // membership in the matching
    if ((b = find_blossom_serial(u, v)) != NULL) {
      // Create a new blist element and add it to the existing list
      blist *l = reinterpret_cast<blist *>(malloc(sizeof(blist)));
      *l = blist(b, seq_solution->blossoms);
      seq_solution->blossoms = l;
      // Remove all vertices but the root of the blossom from the graph
      // and forest, and remove all vertices of the blossom from the
      // matching
      remove_blossom_serial(b);
    } else {
      ilist *path = NULL;
      // Since u and v are not connected, there must be an augmenting
      // path through u and v. Run DFS
      // to find the path to any exposed node
      // Get the path between the exposed nodes
      path = find_path_exposed_serial(u, v);
      // Exchange membership of edges in the path in the matching if we
      // find a path
      if (path != NULL) {
        // Augment the path; should always work.
        modified = augment_serial(path);
        if (!modified) {
          printf("Augmenting failed!\n");  //----------------
          exit(-1);
        }
      }
      if (path == NULL) {  //--------------------
        printf("Found a NULL alternating path.\n");
        exit(-1);
      }
    }
  }
}

/*
 * Runs Edmonds' Blossom Algorithm sequentially
 */
void solve_eba_serial(solution_t *solution) {
  seq_solution = solution;
  modified = 1;
  // Do some basic matchings first
  for (int i = 0; i < n; i++) {
    if (solution->matching[i] == -1) {
      for (int j = 0; j < n; j++) {
        if (solution->graph[i][j] && solution->matching[j] == -1) {
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
    // Keep going until we find a path
    while (flag && !modified) {
      flag = 0;
      // Try every pair of u and v
      for (int u = 0; u < n; u++) {
        if (modified) {
          break;
        }
        if (!(solution->removed[u]) && solution->even[u]) {
          for (int v = 0; v < n; v++) {
            if (modified) {
              break;
            }
            // Check all edges from an even u to a not-odd v
            if (solution->graph[u][v] && !solution->odd[v]) {
              // Add u-v to the forest and see if it completes a blossom
              // or alternating path
              find_alternating_path_serial(u, v);
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
    lift_blossom_serial(current->b);
    temp = current->next;
    free(current);
    current = temp;
  }
  return;
}
