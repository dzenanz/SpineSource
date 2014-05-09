
/* ============================================================================
 * p2/cell/vertex.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include <vec3.h>

#include "cell.hh"
#include "edge.hh"
#include "vertex.hh"

/* ----------------------------------------------------------------------------
 * Vertex
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Vertex *Vertex::make(Cell *cell)
{
  assert(cell!=0);

  return new Vertex(cell);
}

void Vertex::kill(Vertex *vertex)
{
  assert(vertex!=0);

  delete vertex;
}

/* -- public instance methods ---------------------------------------------- */

void Vertex::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}

void Vertex::addEdge(Edge *edge)
{
  assert(edge!=0);

  // only keep track of one edge in the orbit--this one is as good as any

  this->edge = edge;
}

void Vertex::removeEdge(Edge *edge)
{
  assert(edge!=0);

  // replace the arbitrary edge with another edge in the orbit
  // use null if this is the only edge
  // assumes that the edge hasn't been actually removed yet

  Edge *next = edge->Onext();

  this->edge = next!=edge ? next : 0;
}

unsigned Vertex::realIndex()
{
	if (_realIndex==-1) //calculate it if invalid
	{
		Vertex **p = std::find(&cell->vertices[0], &cell->vertices[cell->vertexCount-1], this);
		_realIndex = p - &cell->vertices[0];
	}
	return _realIndex;
}

/* -- protected instance methods ------------------------------------------- */

Vertex::Vertex(Cell *cell)
{
  assert(cell!=0);

  this->pos[0]			= 0.0;
  this->pos[1]			= 0.0;
  this->pos[2]			= 0.0;
  this->subdivPos[0]	= 0.0;
  this->subdivPos[1]	= 0.0;
  this->subdivPos[2]	= 0.0;
  this->cell			= cell;
  this->id				= cell->makeVertexID();
  this->curv			= 0;
  this->dist            = 0;
  this->edge			= 0;
  this->level			= 0;
  this->stepSize		= 0.0;
  this->lastH			= 0;
  this->lastIntensity	= 0;
  this->valence			= 0;
  this->_realIndex      = -1;

  cell->addVertex(this);
}

Vertex::~Vertex()
{
  cell->removeVertex(this);
}

