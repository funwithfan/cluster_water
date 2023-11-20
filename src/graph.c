#include "graph.h"

// Depth-First Search
void dfs(int node, bool *visited, int **adjMatrix, int numNodes, int *connectedComponent, int* componentSize) {
    int i;
    visited[node] = true;
    connectedComponent[(*componentSize)++] = node;
    
    for (i = 0; i < numNodes; i++) {
        if (adjMatrix[node][i] == 1 && !visited[i]) {
            dfs(i, visited, adjMatrix, numNodes, connectedComponent, componentSize);
        }
    }
}


// find connected components and their sizes
void findConnectedComponents(int **adjMatrix, int numNodes, int* numComponents, int *componentSizes, int **connectedComponent) {
    int i, j;
    bool visited[numNodes];

    *numComponents = 0;

    // Initialize visited array
    initializeBoolArray(numNodes, visited, false);
    
    // Traverse through each node
    for (i = 0; i < numNodes; i++) {
        if (!visited[i]) {
            (*numComponents)++;
            int componentSize = 0;
            dfs(i, visited, adjMatrix, numNodes, connectedComponent[i], &componentSize);
            componentSizes[*numComponents - 1] = componentSize;
            //printf("there are %d nodes in component #%d\n", componentSize, *numComponents);
        }
    }
}

// Calculate the degree of each node
void calculateNodesDegree(int **adjMatrix, int numNodes, int* nodesDegree) {
    int i, j, degree;
    for(i = 0;i < numNodes; i++) {
        degree = 0;
        for(j = 0; j < numNodes; j++) {
            degree += adjMatrix[i][j];
        }
        nodesDegree[i] = degree;
    }
}
