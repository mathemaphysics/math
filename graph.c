#include <stdio.h>
#include <stdlib.h>
#include "util_config.h"

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

/**
 * Initializes the vertex_t object with the given identification
 * number and the specified data size.
 * @param obj_in Input vertex_t object
 * @param id_in The id number to assign to the vertex
 * @param alloc_in The size of the data stored in data
 */
int vertex_init( vertex_t *obj_in, int id_in, int alloc_in )
{
	obj_in->data = (void*) malloc( alloc_in );
	if( obj_in->data == NULL )
		return -1;
	obj_in->alloc = alloc_in;
	obj_in->id = 0;
	obj_in->degree = 0;
	obj_in->bflag = 0;
	obj_in->back = NULL;
	obj_in->neighbors = NULL;

	return 0;
}

int vertex_reset( vertex_t *obj_in )
{
	obj_in->bflag = 0;
	obj_in->back = NULL;

	return 0;
}

int vertex_free( vertex_t *obj_in )
{
	free( obj_in->data );

	return 0;
}

/**
 * Initializes the graph object to the specified number of vertices
 * using assumed data size of alloc_in bytes.
 * @param obj_in This is the graph_t object
 * @param nvtx_in The number of vertices to initialize
 * @param alloc_in The number of bytes to allocate for data in each vertex
 */
int graph_init( graph_t *obj_in, int nvtx_in, int alloc_in )
{
	int i;

	obj_in->vertex = (vertex_t*) malloc( nvtx_in * sizeof(vertex_t) );
	obj_in->nvert = nvtx_in;
	obj_in->nedge = 0;
	obj_in->step = alloc_in;
	for(i=0;i<nvtx_in;i++)
		vertex_init( &(obj_in->vertex[i]), i, alloc_in );
}

int graph_add_vertex( graph_t *obj_in, void *data_in )
{
	
}

int graph_add_edge(  )
{
	
}

int graph_free( graph_t *obj_in )
{
	int i;

	for(i=0;i<obj_in->nvert;i++)
	{
		
	}
}

