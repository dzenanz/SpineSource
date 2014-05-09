
/* ============================================================================
 * p2/cell/cell.hh
 * ========================================================================= */

#ifndef cellINCLUDED
#define cellINCLUDED

#include "edge.hh"
#include "face.hh"
#include "vertex.hh"

class CellVertexIterator;
class CellFaceIterator;
class CellForwardVertexIterator;

/* ----------------------------------------------------------------------------
 * Cell
 * ------------------------------------------------------------------------- */

/*
 * An enclosed volume, bounded by a set of vertices and faces.
 *
 * Vertices   the vertices of the cell;
 *            all are nonnull
 * VertexIDs  an increasing sequence of positive integers used to number
 *            distinct vertices;
 *            all are positive
 * Faces      the faces of the cell;
 *            all are nonnull
 * FaceIDs    an increasing sequence of positive integers used to number
 *            distinct faces;
 *            all are positive
 */
class Cell
{

  /* -- public class methods ----------------------------------------------- */

  public:

  /*
   * Return a new, degenerate cell consisting of a single closed edge (loop),
   * a single vertex at the origin, and a pair of faces.
   * <- the new cell;
   *    will be nonnull
   */
  static Cell *make();

  /*
   * Return a new cell with the topology of a tetrahedron and all vertices at
   * the origin.
   * <- the new tetrahedron;
   *    will be nonnull
   */
  static Cell *makeTetrahedron();

  /*
   * Release the storage occupied by a given cell.
   * cell -> the cell to kill;
   *         must be nonnull
   */
  static void kill(Cell *cell);

  /* -- public instance methods (Euler operators) -------------------------- */

  /*
   * Use these methods to construct cells with guaranteed consistent topology.
   * Other means of modifying a cell can potentially produce bizarre results.
   */

  public:

  /*
   * Return a new edge formed by splitting a given vertex between a given pair
   * of faces.
   * A new vertex is introduced at the destination of the new edge.
   * The new edge has _left_ along its left and _right_ along its right.
   * vertex      -> the vertex to split to make the new edge;
   *                must be nonnull;
   *                must share an edge with both _left_ and _right_
   * left, right -> the faces adjacent to the new edge;
   *                must be nonnull;
   *                must share an edge with _vertex_
   * <- the new edge;
   *    will be nonnull
   */
  Edge *makeVertexEdge(Vertex *vertex, Face *left, Face *right);

  /*
   * Delete a given edge from this cell, along with its destination vertex.
   * edge -> the edge to kill;
   *         must be nonnull
   */
  void killVertexEdge(Edge *edge);

  /*
   * Return a new edge formed by splitting a given face through a given pair
   * of vertices.
   * A new face is introduced to the right of the new edge.
   * The new edge has _org_ as its origin and _dest_ as its destination.
   * face      -> the face to divide to make the new edge;
   *              must be nonnull;
   *              must have both _org_ and _dest_ on its perimiter
   * org, dest -> the vertices for the endpoints of the new edge;
   *              must be nonnull;
   *              must be located on the perimiter of _face_
   * <- the new edge;
   *    will be nonnull
   */
  Edge *makeFaceEdge(Face *face, Vertex *org, Vertex *dest);

  /*
   * Delete a given edge from this cell, along with its right face.
   * edge -> the edge to kill;
   *         must be nonnull
   */
  void killFaceEdge(Edge *edge);

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Return the number of vertices in this cell.
   * <- the number of vertices
   */
  unsigned int countVertices();

  /*
   * Add a given vertex to this cell.
   * vertex -> the vertex to add;
   *           must be nonnull;
   *           must not be in the cell
   */
  void addVertex(Vertex *vertex);

  /*
   * Remove a given vertex from this cell.
   * vertex -> the vertex to remove;
   *           must be nonnull;
   *           must be in the cell
   */
  void removeVertex(Vertex *vertex);

  /*
   * Return a new vertex ID.
   * <- a new vertex ID;
   *    will be positive
   */
  unsigned int makeVertexID();

  /*
   * Return the number of faces in this cell.
   * <- the number of faces
   */
  unsigned int countFaces();

  /*
   * Add a given face to this cell.
   * face -> the face to add;
   *         must be nonnull
   *         must not be in the cell
   */
  void addFace(Face *face);

  /*
   * Remove a given face from this cell.
   * face -> the face to remove;
   *         must be nonnull;
   *         must be in the cell
   */
  void removeFace(Face *face);

  /*
   * Return a new face ID.
   * <- a new face ID;
   *    will be positive
   */
  unsigned int makeFaceID();

  /*
   * Performs a deep copy. Very inefficient for large meshes!
   * <- copy of this cell;
   */
  Cell* deepCopy();

  /* -- protected instance methods ----------------------------------------- */

  protected:

  /*
   * Initialize this cell consisting of no vertices and no faces.
   */
  Cell();

  /*
   * Release the storage occupied by the contents of this cell.
   */
  ~Cell();

  /* -- private instance methods ------------------------------------------- */

  private:

  /*
   * Return the edge with a given origin vertex in the face orbit of a given
   * edge.
   * edge -> an edge of the orbit to look for the vertex in;
   *         must be nonnull
   * org  -> the origin vertex to look for;
   *         must be nonnull
   * <- the edge in the same face orbit as _edge_ with origin vertex _org_;
   *    null if not found
   */
  Edge *getOrbitOrg(Edge *edge, Vertex *org);

  /*
   * Set the origin of the vertex orbit of a given edge to a given vertex.
   * edge -> an edge of the orbit to set the origin vertex of;
   *         must be nonnull
   * org  -> the new origin vertex;
   *         must be nonnull
   */
  void setOrbitOrg(Edge *edge, Vertex *org);

  /*
   * Return the edge with a given left face in the vertex orbit of a given
   * edge.
   * edge -> an edge of the orbit to look for the face in;
   *         must be nonnull
   * left -> the left face to look for;
   *         must be nonnull
   * <- the edge in the same vertex orbit as _edge_ with left face _left_;
   *    null if not found
   */
  Edge *getOrbitLeft(Edge *edge, Face *left);

  /*
   * Set the left face of the face orbit of a given edge to a given face.
   * edge -> an edge of the orbit to set the left face of;
   *         must be nonnull
   * left -> the new left face;
   *         must be nonnull
   */
  void setOrbitLeft(Edge *edge, Face *left);

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The vertices in this cell.
   * Nonnull.
   */
  Vertex **vertices;

  /*
   * The number of vertices in this cell.
   * Less than or equal to _vertexSize_.
   */
  unsigned int vertexCount;

  /*
   * The number of vertices allocated for this cell.
   * Greater than or equal to _vertexCount_.
   */
  unsigned int vertexSize;

  /*
   * The next unused vertex ID.
   */
  unsigned int vertexID;

  /*
   * The faces in this cell.
   * Nonnull.
   */
  Face **faces;

  /*
   * The number of faces in this cell.
   * Less than or equal to _faceSize_.
   */
  unsigned int faceCount;

  /*
   * The number of faces allocated for this cell.
   * Greater than or equal to _faceCount_.
   */
  unsigned int faceSize;

  /*
   * The next unused face ID.
   */
  unsigned int faceID;

  /* -- friend classes ----------------------------------------------------- */

  friend class CellVertexIterator;
  friend class CellFaceIterator;
  friend class CellForwardVertexIterator;

  friend unsigned Vertex::realIndex();

};

/* -- inline instance methods ---------------------------------------------- */

inline unsigned int Cell::countVertices()
{
  return vertexCount;
}

inline unsigned int Cell::makeVertexID()
{
  return vertexID++;
}

inline unsigned int Cell::countFaces()
{
  return faceCount;
}

inline unsigned int Cell::makeFaceID()
{
  return faceID++;
}

/* ----------------------------------------------------------------------------
 * CellVertexIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the vertices of a given cell in arbitrary order.
 */
class CellVertexIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this vertex iterator over a given cell.
   * cell -> the cell to iterate the vertices of;
   *         must be nonnull
   */
  CellVertexIterator(Cell *cell)
  {
    this->cell  = cell;
    this->count = cell->vertexCount;
  }

  /*
   * Release the storage occupied by this vertex iterator.
   */
  ~CellVertexIterator()
  {
  }

  /*
   * Return the next vertex of this vertex iterator, if any.
   * <- the next vertex of this vertex iterator;
   *    null if none
   */
  Vertex *next()
  {
    // iterate the array in reverse order so that the current vertex can be
    // removed during iteration

    if (count<1)
      return 0;

    return cell->vertices[--count];
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell whose vertices are being iterated.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The number of vertices left to iterate.
   */
  unsigned int count;

};

/* ----------------------------------------------------------------------------
 * CellFaceIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the faces of a given cell in arbitrary order.
 */
class CellFaceIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this face iterator over a given cell.
   * cell -> the cell to iterate the faces of;
   *         must be nonnull
   */
  CellFaceIterator(Cell *cell)
  {
    this->cell  = cell;
    this->count = cell->faceCount;
  }

  /*
   * Release the storage occupied by this face iterator.
   */
  ~CellFaceIterator()
  {
  }

  /*
   * Return the next face of this face iterator, if any.
   * <- the next face of this face iterator;
   *    null if none
   */
  Face *next()
  {
    // iterate the array in reverse order so that the current face can be
    // removed during iteration

    if (count<1)
      return 0;

    return cell->faces[--count];
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell whose faces are being iterated.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The number of faces left to iterate.
   */
  unsigned int count;

};

/* ----------------------------------------------------------------------------
 * CellForwardVertexIterator
 * ------------------------------------------------------------------------- */

/*
 * Enumerates the vertices of a given cell in order in which they were created.
 */
class CellForwardVertexIterator
{

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Initialize this vertex iterator over a given cell.
   * cell -> the cell to iterate the vertices of;
   *         must be nonnull
   */
  CellForwardVertexIterator(Cell *cell)
  {
    this->cell  = cell;
	this->currentIndex = 0;
  }

  /*
   * Release the storage occupied by this vertex iterator.
   */
  ~CellForwardVertexIterator()
  {
  }

  /*
   * Return the next vertex of this vertex iterator, if any.
   * <- the next vertex of this vertex iterator;
   *    null if none
   */
  Vertex *next()
  {
	  if (currentIndex>=cell->vertexCount)
		return 0;
    return cell->vertices[currentIndex++];
  }

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The cell whose vertices are being iterated.
   * Nonnull.
   */
  Cell *cell;

  /*
   * The number of vertices left to iterate.
   */
  unsigned int currentIndex;

};

#endif /* #ifndef cellINCLUDED */

