#ifndef _EBA_H_
#define _EBA_H_

/* The maximum number of vertices and maximum number of edges in the
 * matching.
 */
#define MAX_N 500
#define NTHREADS 6

// List of integers for blossom, forms a cycle at root of blossom
struct ilist {
  int vertex;
  ilist *next;
  
  ilist(int v, ilist *n) {
    vertex = v;
    next = n;
  }
};

// Contracted blossom with original root
struct blossom {
  int modified[MAX_N][MAX_N];
  int vertex;
  ilist *first;
  
  blossom(int v, ilist *path) {
    vertex = v;
    first = path;
  }
};

// List of blossoms
struct blist {
  blossom *b;
  blist *next;
  
  blist(blossom *new_b, blist *new_l) {
    b = new_b;
    next = new_l;
  }
};

typedef struct solution_s {
  // Locks for every vertex
  omp_lock_t lock[MAX_N];
  short owner[MAX_N];
  
  // These are shrunk with the blossoms but persist between augmentations
  short matching[MAX_N];
  bool graph[MAX_N][MAX_N];
  bool forest[MAX_N][MAX_N];
  bool removed[MAX_N];
  
  // These are reset after every augmentation
  bool even[MAX_N];
  bool odd[MAX_N];
  
  // These are added to as we shrink blossoms and then lifted at the end
  blist *blossoms;
  omp_lock_t blist_lock;
} solution_t;

void solve_eba(solution_t* solution);

#endif
