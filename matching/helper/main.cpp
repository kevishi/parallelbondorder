/* Copyright 2013 Eric Chang and Rutwik Parikh */

#include <fcntl.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>  // For STDIN_FILENO.

#include "./CycleTimer.h"
#include "./eba.h"

/* The number of vertices in the graph. */
int n = 1;

/* The adjacency matrix of the graph. adj[i][j] is 1 if an edge
 * exists between vertices i and j, 0 otherwise. */
bool adj[MAX_N][MAX_N];

void solve_eba(solution_t *solution);
void solve_eba_serial(solution_t *solution);

// ============================================================================
// Timer: returns time in seconds
// ============================================================================

void usage(char *program) {
  printf("Usage: %s [options]\n"
         "  Solves the Maximal Matching problem using the edges provided \n"
         "  on stdin.\n"
         "\n"
         "Program Options:\n"
         "  -b  --bench        Run benchmarking\n"
         "  -i <file>          Use <file> instead of standard in\n"
         "  -?  --help         This message\n", program);
}

void parse_args(int argc, char **argv) {
  int opt;

  static struct option long_opts[] = {
    {"benchmark", 0, 0, 'b'},
    {"input", 1, 0, 'i'},
    {"help", 0, 0, '?'},
    {0, 0, 0, 0}
  };

  while ((opt = getopt_long(argc, argv, "i:?h", long_opts, NULL)) != EOF) {
    switch (opt) {
    case 'i':
      if (dup2(open(optarg, O_RDONLY), STDIN_FILENO) < 0) {
        perror(optarg);
        exit(EXIT_FAILURE);
      }
      break;
    case 'h':                  /* Explicit fall through */
    case '?':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
}

/* Parses standard in for an adjaceny matrix, loading into the global array. */
void parse_adj_matrix() {
  scanf("%d", &n);
  if ((unsigned) n > MAX_N) {
    printf("Invalid number of vertices: saw %d but max is %d\n", n, MAX_N);
    exit(EXIT_FAILURE);
  }

  int i, j;
  int cur;

  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      scanf("%d", &cur);
      if (cur < 0 || cur > 1) {
        printf("Invalid input; inputs must be 0 or 1.\n");
      }
      adj[i][j] = cur;
      adj[j][i] = adj[i][j];
    }
  }
}

int main(int argc, char **argv) {
  parse_args(argc, argv);
  parse_adj_matrix();
  omp_set_num_threads(NTHREADS);
  //---------------------SEQUENTIAL SOLUTION------------------
  // Initialize solution
  solution_t solution;
  memset(&solution, 0, sizeof(solution));
  memcpy(solution.graph, adj, sizeof(adj));
  memset(solution.matching, -1, sizeof(solution.matching));

  double startTime = CycleTimer::currentSeconds();
  solve_eba_serial(&solution);
  double endTime = CycleTimer::currentSeconds();

  printf("\nSequential solution found!\n");
  int matched_vertices[MAX_N];
  memset(matched_vertices, 0, sizeof(matched_vertices));
  int matching_size = 0;
  printf("Solution took %.4fs sequentially\n",
      endTime - startTime);
  //---------------------SEQUENTIAL SOLUTION AGAIN------------------
  // Note that we run the sequential solution again to take into
  // account the effects of locality of the prior run
  // Initialize solution
  memset(&solution, 0, sizeof(solution));
  memcpy(solution.graph, adj, sizeof(adj));
  memset(solution.matching, -1, sizeof(solution.matching));

  startTime = CycleTimer::currentSeconds();
  solve_eba_serial(&solution);
  endTime = CycleTimer::currentSeconds();

  printf("\nSequential solution found (again)!\n");
  memset(matched_vertices, 0, sizeof(matched_vertices));
  matching_size = 0;
  // We print out every edge in the matching with the smaller vertex first.
  // Note that i = n-1 either has already been matched or is not in the
  // matching
  for (int i = 0; i < n-1; i++) {
    // Nothing to check if i is not in the matching
    if (solution.matching[i] == -1) {
      continue;
    }
    if (!(adj[i][solution.matching[i]])) {
      printf("Error: Edge in matching %d -> %d not in graph.\n",
          i, solution.matching[i]);
    }
    // If i is not matched, then matching[i] will be -1
    if (i < solution.matching[i]) {
      // Check if the vertices are already touching an edge
      // in the matching. 0 means the vertex is not in the matching
      // yet, 1 means it is in the matching
      if (matched_vertices[i] == 0) {
        matched_vertices[i] = 1;
        if (matched_vertices[solution.matching[i]] == 0) {
          matched_vertices[solution.matching[i]] = 1;
          // If the edge is valid, increment the matching size
          matching_size++;
        } else {
          printf("Error: Edge %d -> %d is invalid, vertex %d already in "
              "the matching.\n",
              i, solution.matching[i], solution.matching[i]);
        }
      } else {
          printf("Error: Edge %d -> %d is invalid, vertex %d already in "
              "the matching.\n", i, solution.matching[i], i);
      }
    }
    if (solution.matching[solution.matching[i]] != i) {
      printf("Error: Matching %d -> %d not mutual; %d -> %d.\n",
            i, solution.matching[i], solution.matching[i],
            solution.matching[solution.matching[i]]);
    }
  }
  printf("Matching contained %d edges.\n", matching_size);
  printf("Solution took %.4fs sequentially\n",
      endTime - startTime);

  /*for (int i =  0; i < n; i++) {
    if (i % 5 == 0) {
      printf("\n");
    }
    printf("%d: %d\t", i, solution.matching[i]);
  }
  printf("\n");*/
  //--------------------OMP SOLUTION--------------------
  // Initialize solution
  memset(&solution, 0, sizeof(solution));
  memcpy(solution.graph, adj, sizeof(adj));
  memset(solution.matching, -1, sizeof(solution.matching));
  memset(solution.owner, -1, sizeof(solution.owner));
  // Initialize locks
  for (int i = 0; i < n; i++) {
    omp_init_lock(&(solution.lock[i]));
  }
  omp_init_lock(&(solution.blist_lock));

  startTime = CycleTimer::currentSeconds();
  solve_eba(&solution);
  endTime = CycleTimer::currentSeconds();

  printf("\nOMP solution found!\n");
  memset(matched_vertices, 0, sizeof(matched_vertices));
  matching_size = 0;
  // We print out every edge in the matching with the smaller vertex first.
  // Note that i = n-1 either has already been matched or is not in the
  // matching
  for (int i = 0; i < n-1; i++) {
    // Nothing to check if i is not in the matching
    if (solution.matching[i] == -1) {
      continue;
    }
    if (!(adj[i][solution.matching[i]])) {
      printf("Error: Edge in matching %d -> %d not in graph.\n",
          i, solution.matching[i]);
    }
    // If i is not matched, then matching[i] will be -1
    if (i < solution.matching[i]) {
      // Check if the vertices are already touching an edge
      // in the matching. 0 means the vertex is not in the matching
      // yet, 1 means it is in the matching
      if (matched_vertices[i] == 0) {
        matched_vertices[i] = 1;
        if (matched_vertices[solution.matching[i]] == 0) {
          matched_vertices[solution.matching[i]] = 1;
          // If the edge is valid, increment the matching size
          matching_size++;
        } else {
          printf("Error: Edge %d -> %d is invalid, vertex %d already in "
              "the matching.\n",
              i, solution.matching[i], solution.matching[i]);
        }
      } else {
          printf("Error: Edge %d -> %d is invalid, vertex %d already in "
              "the matching.\n", i, solution.matching[i], i);
      }
    }
    if (solution.matching[solution.matching[i]] != i) {
      printf("Error: Matching %d -> %d not mutual; %d -> %d.\n",
            i, solution.matching[i], solution.matching[i],
            solution.matching[solution.matching[i]]);
    }
  }
  printf("Matching contained %d edges.\n", matching_size);
  printf("Solution took %.4fs on %d processors\n",
      endTime - startTime, omp_get_max_threads());

  /*for (int i =  0; i < n; i++) {
    if (i % 5 == 0) {
      printf("\n");
    }
    printf("%d: %d\t", i, solution.matching[i]);
  }
  printf("\n");*/
  // Destroy locks
  for (int i = 0; i < n; i++) {
    omp_destroy_lock(&(solution.lock[i]));
  }

  return 0;
}
