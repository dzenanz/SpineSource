
/* ============================================================================
 * p2/cell/face.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell.hh"
#include "edge.hh"
#include "face.hh"

/* ----------------------------------------------------------------------------
 * Face
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Face *Face::make(Cell *cell)
{
  assert(cell!=0);

  return new Face(cell);
}

void Face::kill(Face *face)
{
  assert(face!=0);

  delete face;
}

/* -- public instance methods ---------------------------------------------- */

void Face::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}

void Face::addEdge(Edge *edge)
{
  assert(edge!=0);

  // only keep track of one edge in the orbit--this one is as good as any

  this->edge = edge;
}

void Face::removeEdge(Edge *edge)
{
  assert(edge!=0);

  // replace the arbitrary edge with another edge in the orbit
  // use null if this is the only edge
  // assumes that the edge hasn't been actually removed yet

  Edge *next = edge->Lnext();

  this->edge = next!=edge ? next : 0;
}

/* -- protected instance methods ------------------------------------------- */

Face::Face(Cell *cell)
{
  assert(cell!=0);

  this->cell = cell;
  this->id   = cell->makeFaceID();
  this->edge = 0;

  cell->addFace(this);
}

Face::~Face()
{
  cell->removeFace(this);
}

