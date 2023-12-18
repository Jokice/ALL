#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_VERTICES 100
#define INF INT_MAX
/*// Node to represent an edge in the adjacency list
typedef struct Node{
    int vertex;
    int weight;
    struct Node* next;
} Node;

// List to represent the adjacency list for each vertex
typedef struct List {
    Node* head;
} List;

// Graph structure with an array of adjacency lists
typedef struct {
    int vertices;
    List* array;
} Graphadj;

// Function to create a new node
Node* createNode(int vertex, int weight) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->vertex = vertex;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}
// Function to create a graph with a given number of vertices
Graphadj* createGraph(int vertices)
 {
    Graphadj* graph = (Graphadj*)malloc(sizeof(Graph));
    graph->vertices = vertices;
    graph->array = (List*)malloc(vertices * sizeof(List));

    // Initialize each adjacency list as empty
    for (int i = 0; i < vertices; ++i) {
        graph->array[i].head = NULL;
    }

    return graph;
}

// Function to add an edge to an undirected graph
void addEdge(Graphadj* graph, int src, int dest, int weight) {
    // Add an edge from src to dest
    Node* newNode = createNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;

    // Since the graph is undirected, add an edge from dest to src as well
    newNode = createNode(src, weight);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}

// Function to print the adjacency list representation of the graph
void printGraph(Graphadj* graph) {
    for (int i = 0; i < graph->vertices; ++i) {
        printf("%d -> ", i);
        Node* current = graph->array[i].head;
        while (current != NULL) {
            printf("(%d, %d) -> ", current->vertex, current->weight);
            current = current->next;
        }
        printf("NULL\n");
    }
}

*/
// 定义图的结构体
typedef struct 
{
    int vertices; // 顶点数
    int edges;    // 边数
    int matrix[MAX_VERTICES][MAX_VERTICES]; // 邻接矩阵
} Graph;

// 创建无向加权图的邻接矩阵
void createWeightedGraph(Graph *G) 
{
    printf("Enter the number of vertices: ");
    scanf("%d",&G->vertices);

    printf("Enter the number of edges: ");
    scanf("%d", &G->edges);

    // 初始化邻接矩阵
    for (int i = 0; i < G->vertices; ++i) {
        for (int j = 0; j < G->vertices; ++j) {
            G->matrix[i][j] = INF;
        }
    }

    // 输入边的权值
    for (int k = 0; k < G->edges; ++k) {
        int i, j, weight;
        printf("Enter edge (i, j) and its weight: ");
        scanf("%d %d %d", &i, &j, &weight);
        G->matrix[i][j] = weight;
        G->matrix[j][i] = weight; // 无向图的邻接矩阵是对称的
    }
}

// 非递归深度优先搜索
void DFS_non_recursive(Graph *G, int start) 
{
    int stack[MAX_VERTICES];
    int top = -1;
    int visited[MAX_VERTICES] = {0};

    printf("DFS Non-Recursive: ");
    printf("%d ", start);
    visited[start] = 1;
    stack[++top] = start;

    while (top != -1)
	{
        int current = stack[top];

        int found = 0;
        for (int i = 0; i < G->vertices; ++i) 
		{
            if (G->matrix[current][i] != INF && !visited[i])
			 {
                printf("%d ", i);
                visited[i] = 1;
                stack[++top] = i;
                found = 1;
                break;
            }
        }

        if (!found)
		{
            --top;
        }
    }

    printf("\n");
}

// 广度优先搜索
void BFS(Graph *G, int start)
{
    int queue[MAX_VERTICES];
    int front = 0, rear = 0;
    int visited[MAX_VERTICES] = {0};

    printf("BFS: ");
    printf("%d ", start);
    visited[start] = 1;
    queue[rear++] = start;

    while (front != rear) {
        int current = queue[front++];

        for (int i = 0; i < G->vertices; ++i) {
            if (G->matrix[current][i] != INF && !visited[i]) {
                printf("%d ", i);
                visited[i] = 1;
                queue[rear++] = i;
            }
        }
    }

    printf("\n");
}

// Prim 算法构建最小生成树
void primMST(Graph*G) 
{
    int parent[MAX_VERTICES]; // 保存最小生成树的父节点
    int key[MAX_VERTICES];    // 保存顶点的权值
    int inMST[MAX_VERTICES];  // 记录顶点是否在最小生成树中

    // 初始化 key 和 inMST 数组
    for (int i = 0; i < G->vertices; ++i) {
        key[i] = INF;
        inMST[i] = 0;
    }

    // 从顶点 0 开始构建最小生成树
    key[0] = 0;
    parent[0] = -1;

    for (int count = 0; count < G->vertices - 1; ++count)
	 {
        int minKey = INF;
        int u;

        // 选取权值最小的顶点
        for (int v = 0; v < G->vertices; ++v) 
		{
            if (!inMST[v] && key[v] < minKey)
			{
                minKey = key[v];
                u = v;
            }
        }

        inMST[u] = 1;

        // 更新与 u 相邻的顶点的权值
        for (int v = 0; v < G->vertices; ++v) 
		{
            if (G->matrix[u][v] != INF && !inMST[v] && G->matrix[u][v] < key[v]) 
			{
                parent[v] = u;
                key[v] = G->matrix[u][v];
            }
        }
    }

    // 输出最小生成树的边
    printf("Prim's MST Edges:\n");
    for (int i = 1; i < G->vertices; ++i) 
	{
        printf("(%d, %d) - Weight: %d\n", parent[i], i, G->matrix[i][parent[i]]);
    }
}

// Dijkstra 算法求最短路径
void dijkstra(Graph *G, int start) 
{
    int distance[MAX_VERTICES]; // 保存从 start 到各顶点的最短距离
    int visited[MAX_VERTICES];  // 记录顶点是否被访问过

    // 初始化 distance 和 visited 数组
    for (int i = 0; i < G->vertices; ++i) 
	{
        distance[i] = INF;
        visited[i] = 0;
    }

    // 从起始顶点开始
    distance[start] = 0;

    for (int count = 0; count < G->vertices - 1; ++count) 
	{
        int minDist = INF;
        int u;
        // 选取未访问顶点中距离最小的顶点
        for (int v = 0; v < G->vertices; ++v)
		 {
            if (!visited[v] && distance[v] < minDist)
			 {
                minDist = distance[v];
                u = v;
            }
        }

        visited[u] = 1;

        // 更新与 u 相邻的顶点的距离
        for (int v = 0; v < G->vertices; ++v) 
		{
                if (!visited[v] && G->matrix[u][v] != INF && distance[u] != INF && 
                distance[u] + G->matrix[u][v] < distance[v])
				 {
                distance[v] = distance[u] + G->matrix[u][v];
                }
        }
    }

    // 输出从起始顶点到每个顶点的最短路径
    printf("Dijkstra's Shortest Paths from vertex %d:\n", start);
    for (int i = 0; i < G->vertices; ++i) 
	{
        printf("To vertex %d: ",i);
        if (distance[i] == INF)
		{
            printf("No path\n");
        } 
		else 
		{
            // 输出路径中的顶点信息
            int path[MAX_VERTICES];
            int pathLength = 0;
            // 回溯构建路径
            int j = i;
            while (j != start) 
			{
                path[pathLength++] = j;
                j = path[j];
            }
            // 输出路径
            printf("%d ", start);
            for (int k = pathLength - 1; k >= 0; --k) 
			{
                printf("%d ", path[k]);
            }
            printf(" - Distance: %d\n", distance[i]);
        }
    }
}

int main() {
    Graph G;
    createWeightedGraph(&G);
     // Create a graph with 4 vertices
    /*Graphadj* graph = createGraph(4);

    // Add edges to the graph
    addEdge(graph, 0, 1, 2);
    addEdge(graph, 0, 2, 4);
    addEdge(graph, 1, 2, 1);
    addEdge(graph, 1, 3, 7);
    addEdge(graph, 2, 3, 3);

    // Print the adjacency list representation of the graph
    printf("Adjacency List Representation of the Graph:\n");
    printGraph(graph);
    */
    // 非递归深度优先搜索
    int startDFS = 0;
    DFS_non_recursive(&G, startDFS);

    // 广度优先搜索
    int startBFS = 0;
    BFS(&G, startBFS);

    // Prim算法构建最小生成树
    primMST(&G);

    // Dijkstra算法求最短路径
    int startDijkstra = 0;
    dijkstra(&G, startDijkstra);

    return 0;
}

