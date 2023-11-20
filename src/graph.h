#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <stdbool.h>
#include "header.h"
#include "basic.h"

void dfs(int node, bool *visited, int **adjMatrix, int numNodes, int *connectedComponent, int* componentSize);
void findConnectedComponents(int **adjMatrix, int numNodes, int* numComponents, int *componentSizes, int **connectedComponent);