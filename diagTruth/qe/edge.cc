
/* ============================================================================
 * p2/cell/edge.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

/* ----------------------------------------------------------------------------
 * QuadEdge
 * ------------------------------------------------------------------------- */

class QuadEdge
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize the edges of this quad edge with no connections.
   */
  QuadEdge()
  {
    edges[0].index = 0;
    edges[1].index = 1;
    edges[2].index = 2;
    edges[3].index = 3;

    edges[0].next = edges+0;
    edges[1].next = edges+3;
    edges[2].next = edges+2;
    edges[3].next = edges+1;

    unsigned int id = Edge::nextID;

    edges[0].id = id+0;
    edges[1].id = id+1;
    edges[2].id = id+2;
    edges[3].id = id+3;

    Edge::nextID = id+4;
  }

  /* -- public instance variables ------------------------------------------ */

  public:

  /*
   * The edges of this quad edge.
   */
  Edge edges[4];

};

/* ----------------------------------------------------------------------------
 * Edge
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Edge *Edge::make()
{
  return (new QuadEdge())->edges;
}

void Edge::kill(Edge *edge)
{
  assert(edge!=0);

  // detach the edge from its cell
  splice(edge, edge->Oprev());
  splice(edge->Sym(), edge->Sym()->Oprev());

  // free the quad edge that the edge belongs to
  delete (QuadEdge*)(edge-edge->index);
}

void Edge::splice(Edge *a, Edge *b)
{
  assert(a!=0);
  assert(b!=0);

  // see Guibas and Stolfi

  Edge* alpha = a->Onext()->Rot();
  Edge* beta  = b->Onext()->Rot();

  Edge* t1 = b->Onext();
  Edge* t2 = a->Onext();
  Edge* t3 = beta->Onext();
  Edge* t4 = alpha->Onext();

  a->next     = t1;
  b->next     = t2;
  alpha->next = t3;
  beta->next  = t4;
}

/* -- public instance methods ---------------------------------------------- */

void Edge::setID(unsigned int id)
{
  assert(id>0);

  this->id = id;
}

void Edge::setOrg(Vertex *org)
{
  // add this edge to the (vertex) orbit of _org_

  vertex = org;

  org->addEdge(this);
}

void Edge::setDest(Vertex *dest)
{
  // add this edge to the (vertex) orbit of _dest_

  Sym()->vertex = dest;

  dest->addEdge(Sym());
}

void Edge::setLeft(Face *left)
{
  // add this edge to the (face) orbit of _left_

  Rot()->face = left;

  left->addEdge(this);
}

void Edge::setRight(Face *right)
{
  // add this edge to the (face) orbit of _right_

  InvRot()->face = right;

  right->addEdge(Sym());
}

Real Edge::length()
{
	return (Org()->pos - Dest()->pos).length();
}

vec3 Edge::asVector()
{
	return Dest()->pos-Org()->pos;
}

/* -- protected instance methods ------------------------------------------- */

Edge::Edge()
{
  // _index_ is initialized by QuadEdge
  // _next_ is initialized by QuadEdge
  // _id_ is initialized by QuadEdge

  vertex = 0;
  face   = 0;
  midpoint=vec3(std::numeric_limits<float>::quiet_NaN(),0,0);
}

Edge::~Edge()
{
}

/* -- private class variables ---------------------------------------------- */

unsigned int Edge::nextID = 4;

