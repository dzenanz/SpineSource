
/* ============================================================================
 * p2/cell/cell.cc
 * ========================================================================= */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell.hh"
#include "face.hh"
#include "vertex.hh"
#include "obj.hh"

/* ----------------------------------------------------------------------------
 * Cell
 * ------------------------------------------------------------------------- */

/* -- public class methods ------------------------------------------------- */

Cell *Cell::make()
{
  // create a looping edge that connects to itself at a single vertex
  // the edge delimits two faces
  // this is the smallest cell that is consistent with our invariants

  Cell *cell = new Cell();

  Vertex *vertex = Vertex::make(cell);
  Face   *left   = Face::make(cell);
  Face   *right  = Face::make(cell);
  Edge   *edge   = Edge::make()->InvRot();

  edge->setOrg(vertex);
  edge->setDest(vertex);
  edge->setLeft(left);
  edge->setRight(right);

  return cell;
}

Cell *Cell::makeTetrahedron()
{
  // use the Euler operators to construct a tetrahedron

  Cell *cell = make();

  // grab the initial vertex

  Vertex *vertex1;

  {
    CellVertexIterator iterator(cell);

    vertex1 = iterator.next();

    assert(vertex1!=0);
  }

  // grab the initial edge and the initial faces

  Edge *edge1 = vertex1->getEdge();
  Face *left  = edge1->Left();
  Face *right = edge1->Right();

  // drop in four vertices along the initial edge

  Vertex *vertex2 = cell->makeVertexEdge(vertex1, left, right)->Dest();
  Vertex *vertex3 = cell->makeVertexEdge(vertex2, left, right)->Dest();
  Vertex *vertex4 = cell->makeVertexEdge(vertex3, left, right)->Dest();

  // cut each of the faces in half from complementary vertices

  Face *front  = cell->makeFaceEdge(left, vertex2, vertex4)->Right();
  Face *bottom = cell->makeFaceEdge(right, vertex1, vertex3)->Right();

  return cell;
}

void Cell::kill(Cell *cell)
{
  assert(cell!=0);

  delete cell;
}

/* -- public instance methods (Euler operators) ---------------------------- */

Edge *Cell::makeVertexEdge(Vertex *vertex, Face *left, Face *right)
{
  assert(vertex!=0);
  assert(left!=0);
  assert(right!=0);

  // locate the edges to the right of each of the faces in the orbit of the
  // vertex

  Edge *edge  = vertex->getEdge();
  Edge *edge1 = getOrbitLeft(edge, right);
  Edge *edge2 = getOrbitLeft(edge, left);

  if (edge1==0)
  {
    (void)fprintf(stderr, "Cell::makeVertexEdge: unable to locate right face %u on vertex %u\n",
		 right->getID(), vertex->getID());
    abort();
  }

  if (edge2==0)
  {
    (void)fprintf(stderr, "Cell::makeVertexEdge: unable to locate left face %u on vertex %u\n",
		 left->getID(), vertex->getID());
    abort();
  }

  // create a new vertex and copy the position of the vertex of origin

  Vertex *vertexNew = Vertex::make(this);

  vertexNew->pos = vertex->pos;

  // create a new edge and rotate it to make a clockwise loop

  Edge *edgeNew = Edge::make()->Rot();

  // connect the origin (and destination) of the new edge to _vertex_ so that
  // the left face of the edge is _left_
  // this makes a loop on the inside of _left_

  Edge::splice(edge2, edgeNew);

  // split the origin and destination of the loop so that the right face of the
  // edge is now _right_
  // this results in a non-loop edge dividing _left_ from _right_

  Edge::splice(edge1, edgeNew->Sym());

  // initialize the secondary attributes of the new edge

  edgeNew->setOrg(edge1->Org());
  edgeNew->setLeft(edge2->Left());
  edgeNew->setRight(edge1->Left());

  // all edges leaving the destination orbit of the new edge now have the new
  // vertex as their vertex of origin

  setOrbitOrg(edgeNew->Sym(), vertexNew);

  return edgeNew;
}

void Cell::killVertexEdge(Edge *edge)
{
  assert(edge!=0);

  // locate _edge1_ and _edge2_ as in _makeVertexEdge_

  Edge *edge1 = edge->Oprev();
  Edge *edge2 = edge->Lnext();

  // use _edge1_ for _edge2_ if the destination vertex is isolated

  if (edge2==edge->Sym())
    edge2 = edge1;

  // inverse of _makeVertexEdge_

  Edge::splice(edge1, edge->Sym());
  Edge::splice(edge2, edge);

  // all edges leaving the destination orbit of the deleted edge now have its
  // origin vertex as their vertex of origin

  setOrbitOrg(edge2, edge1->Org());

  // don't use the deleted edge as a reference edge any more

  edge1->Org()->addEdge(edge1);
  edge1->Left()->addEdge(edge1);
  edge2->Left()->addEdge(edge2);

  // reclaim the vertex and the edge

  Vertex::kill(edge->Dest());
  Edge::kill(edge);
}

Edge *Cell::makeFaceEdge(Face *face, Vertex *org, Vertex *dest)
{
  assert(face!=0);
  assert(org!=0);
  assert(dest!=0);

  // locate the edges leaving each of the vertices in the orbit of the face

  Edge *edge  = face->getEdge();
  Edge *edge1 = getOrbitOrg(edge, org);
  Edge *edge2 = getOrbitOrg(edge, dest);

  if (edge1==0)
  {
    (void)fprintf(stderr, "Cell::makeFaceEdge: unable to locate origin vertex %u on face %u\n",
		 org->getID(), face->getID());
    abort();
  }

  if (edge2==0)
  {
    (void)fprintf(stderr, "Cell::makeFaceEdge: unable to locate destination vertex %u on face %u\n",
		 dest->getID(), face->getID());
    abort();
  }

  // create a new face

  Face *faceNew = Face::make(this);

  // create a new (non-loop) edge

  Edge *edgeNew = Edge::make();

  // connect the destination of the new edge to the origin of _edge2_
  // both faces of the edge are now _face_

  Edge::splice(edge2, edgeNew->Sym());

  // connect the origin of the new edge to _edge1_
  // _face_ is split in half along the new edge, with the new face introduced
  // on the right

  Edge::splice(edge1, edgeNew);

  // initialize the secondary attributes of the new edge

  edgeNew->setOrg(edge1->Org());
  edgeNew->setDest(edge2->Org());
  edgeNew->setLeft(edge2->Left());

  // all edges in the right orbit of the new edge (i.e. the left orbit of its
  // Sym) now have the new face as their left face

  setOrbitLeft(edgeNew->Sym(), faceNew);

  return edgeNew;
}

void Cell::killFaceEdge(Edge *edge)
{
  assert(edge!=0);

  // locate _edge1_ and _edge2_ as in _makeFaceEdge_

  Edge *edge1 = edge->Oprev();
  Edge *edge2 = edge->Lnext();

  // use _edge2_ for _edge1_ if the right face is inside a loop

  if (edge1==edge->Sym())
    edge1 = edge2;

  // inverse of _makeFaceEdge_

  Edge::splice(edge2, edge->Sym());
  Edge::splice(edge1, edge);

  // all edges in the right orbit of the deleted edge now have its left face
  // as their left face

  setOrbitLeft(edge1, edge2->Left());

  // don't use the deleted edge as a reference edge any more

  edge1->Org()->addEdge(edge1);
  edge2->Org()->addEdge(edge2);
  edge2->Left()->addEdge(edge2);

  // reclaim the face and the edge

  Face::kill(edge->Right());
  Edge::kill(edge);
}

/* -- public instance methods ---------------------------------------------- */

void Cell::addVertex(Vertex *vertex)
{
  assert(vertex!=0);

  // expand the vertex array, if necessary

  if (vertexCount>=vertexSize)
  {
    unsigned int vertexSizeNew = vertexSize*2;
    Vertex     **verticesNew   = new Vertex*[vertexSizeNew];

    for (unsigned int i = 0; i<vertexCount; i++)
      verticesNew[i] = vertices[i];

    delete[] vertices;

    vertices   = verticesNew;
    vertexSize = vertexSizeNew;
  }

  // add the vertex to the end of the array

  vertices[vertexCount++] = vertex;
}

void Cell::removeVertex(Vertex *vertex)
{
  assert(vertex!=0);

  // locate the vertex in the array and replace it with the current last vertex
  // if already the last vertex, just overwrite it
  // slow: should make this a doubly-linked list ???

  for (unsigned int i = vertexCount; i>0; i--)
    if (vertices[i-1]==vertex)
    {
      vertices[i-1] = vertices[--vertexCount];
      return;
    }

  //assert(0);
}

void Cell::addFace(Face *face)
{
  assert(face!=0);

  // expand the face array, if necessary

  if (faceCount>=faceSize)
  {
    unsigned int faceSizeNew = faceSize*2;
    Face       **facesNew    = new Face*[faceSizeNew];

    for (unsigned int i = 0; i<faceCount; i++)
      facesNew[i] = faces[i];

    delete[] faces;

    faces    = facesNew;
    faceSize = faceSizeNew;
  }

  // add the face to the end of the array

  faces[faceCount++] = face;
}

void Cell::removeFace(Face *face)
{
  assert(face!=0);

  // locate the face in the array and replace it with the current last face
  // if already the last face, just overwrite it
  // slow: should make this a doubly-linked list ???

  for (unsigned int i = faceCount; i>0; i--)
    if (faces[i-1]==face)
    {
      faces[i-1] = faces[--faceCount];
      return;
    }

  //assert(0);
}

Cell* Cell::deepCopy()
{
	return objCloneCell(this);
}

/* -- protected instance methods ------------------------------------------- */

Cell::Cell()
{
  // preallocate enough elements for a cube

  vertices    = new Vertex*[8];
  vertexCount = 0;
  vertexSize  = 8;
  vertexID    = 1;

  faces     = new Face*[6];
  faceCount = 0;
  faceSize  = 6;
  faceID    = 1;
}

Cell::~Cell()
{
  // reclaim each of the vertices and faces still owned by the cell
  // go in backwards order so that when the elements try to remove themselves,
  // it will be linear time

  {
    for (unsigned int i = vertexCount; i>0; i--)
      Vertex::kill(vertices[i-1]);
  }

  {
    for (unsigned int i = faceCount; i>0; i--)
      Face::kill(faces[i-1]);
  }

  // reclaim the vertex and face arrays

  delete[] vertices;
  delete[] faces;
}

/* -- private instance methods --------------------------------------------- */

Edge *Cell::getOrbitOrg(Edge *edge, Vertex *org)
{
  assert(edge!=0);
  assert(org!=0);

  // traverse the Lnext orbit of _edge_ looking for an edge whose origin is
  // _org_

  Edge *scan = edge;

  do
  {
    if (scan->Org()==org)
      return scan;

    scan = scan->Lnext();
  }
  while (scan!=edge);

  return 0;
}

void Cell::setOrbitOrg(Edge *edge, Vertex *org)
{
  assert(edge!=0);
  assert(org!=0);

  // traverse the Onext orbit of _edge_, setting the origin of each edge to
  // _org_

  Edge *scan = edge;

  do
  {
    scan->setOrg(org);

    scan = scan->Onext();
  }
  while (scan!=edge);
}

Edge *Cell::getOrbitLeft(Edge *edge, Face *left)
{
  assert(edge!=0);
  assert(left!=0);

  // traverse the Onext orbit of _edge_ looking for an edge whose left face is
  // _left

  Edge *scan = edge;

  do
  {
    if (scan->Left()==left)
      return scan;

    scan = scan->Onext();
  }
  while (scan!=edge);

  return 0;
}

void Cell::setOrbitLeft(Edge *edge, Face *left)
{
  assert(edge!=0);
  assert(left!=0);

  // traverse the Lnext orbit of _edge_, setting the left face of each edge to
  // _left_

  Edge *scan = edge;

  do
  {
    scan->setLeft(left);

    scan = scan->Lnext();
  }
  while (scan!=edge);
}

unsigned int Cell::countEdges()
{
    unsigned int count=0;
    CellFaceIterator cfi(this);
    Face *f;
    while ((f=cfi.next())!=0)
    {
        FaceEdgeIterator fei(f);
        while (fei.next()!=0)
            count++;
    }
    return count;
}