#ifndef GRAPH_H
#define GRAPH_H

#define GRAPH_INCREMENT 16

typedef struct vertex
{
	/**
	 * Whatever you want.
	 */
	void *data;

	/**
	 * Size of the data stored in data.
	 */
	int alloc;

	/**
	 * The index identifying the vertex.
	 */
	int id;

	/**
	 * The number of neighbors.
	 */
	int degree;

	/**
	 * Flag set to 1 if the back pointer is set, 0 otherwise.
	 */
	char bflag;

	/**
	 * Use when doing DFS or BFS; indicates the path taken to current vertex.
	 */
	struct vertex *back;

	/**
	 * List of neighboring vertices.
	 */
	struct vertex *neighbors;
} vertex_t;

typedef struct graph
{
	/**
	 * The total number of vertices.
	 */
	int nvert;

	/**
	 * The total number of edges.
	 */
	int nedge;

	/**
	 * The number of bytes to allocate for each vertex's data.
	 */
	int step;

	/**
	 * Pointer to the first vertex.
	 */
	vertex_t *vertex; /* the pointer to the first vertex */
} graph_t;

#endif

