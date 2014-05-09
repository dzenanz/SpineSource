
/* ============================================================================
 * p2/cell/edge.hh
 * ========================================================================= */

#ifndef edgeINCLUDED
#define edgeINCLUDED

#include <vec3.h>
#include <limits>

class Face;
class QuadEdge;
class Vertex;

/* ----------------------------------------------------------------------------
 * Edge
 * ------------------------------------------------------------------------- */

/*
 * A directed edge from one vertex to another, adjacent to two faces.
 * Based on Dani Lischinski's code from Graphics Gems IV.
 * Original quad-edge data structure due to Guibas and Stolfi (1985).
 *
 * ID     the ID number assigned to the edge;
 *        positive
 * data   generic data attached to the edge by the client
 * Org    the vertex of origin for the edge;
 *        null if currently unknown
 * Dest   the destination vertex for the edge;
 *        null if currently unknown
 * Left   the left face of the edge;
 *        null if currently unknown
 * Right  the right face of the edge;
 *        null if currently unknown
 */
class Edge
{

  /* -- public class methods ----------------------------------------------- */

  public:

  /*
   * Return a new, unconnected edge.
   * <- the new edge;
   *    will be nonnull
   */
  static Edge *make();

  /*
   * Release the storage occupied by a given edge.
   * edge -> the edge to kill;
   *         must be nonnull
   */
  static void kill(Edge *edge);

  /*
   * Splice a given pair of edges.
   * a, b -> the edges to splice;
   *         must be nonnull
   *
   * This operator affects the two edge rings around the origins of a and b,
   * and, independently, the two edge rings around the left faces of a and b.
   * In each case, (i) if the two rings are distinct, Splice will combine
   * them into one; (ii) if the two are the same ring, Splice will break it
   * into two separate pieces.
   * Thus, Splice can be used both to attach the two edges together, and
   * to break them apart. See Guibas and Stolfi (1985) p.96 for more details
   * and illustrations.
   */
  static void splice(Edge *a, Edge *b);

  /* -- public instance variables ------------------------------------------ */

  /*
   * The midpoint location for usage in subdivision. Contains NaN if unintialized.
   */
  vec3 midpoint;

  /* -- public instance methods -------------------------------------------- */

  public:

  /*
   * Return the ID of this edge.
   * <- the ID of this edge;
   *    will be positive
   */
  unsigned int getID();

  /*
   * Change the ID of this edge.
   * id -> the new id for this edge;
   *       must be positive
   */
  void setID(unsigned int id);

  /*
   * Return the origin vertex of this edge.
   * <- the origin of this edge;
   *    null if currently unknown
   */
  Vertex *Org();

  /*
   * Return the destination vertex of this edge.
   * <- the destination of this edge;
   *    null if currently unknown
   */
  Vertex *Dest();

  /*
   * Change the origin vertex of this edge to a given vertex.
   * org  -> the new origin vertex of this edge;
   *         null if currently unknown
   */
  void setOrg(Vertex *org);

  /*
   * Change the destination vertex of this edge to a given vertex.
   * dest -> the new destination vertex of this edge;
   *         null if currently unknown
   */
  void setDest(Vertex *dest);

  /*
   * Return the left face of this edge.
   * <- the left face of this edge;
   *    null if currently unknown
   */
  Face *Left();

  /*
   * Return the right face of this edge.
   * <- the right face of this edge;
   *    null if currently unknown
   */
  Face *Right();

  /*
   * Change the left face of this edge to a given face.
   * left  -> the new left face of this edge;
   *          null if currently unknown
   */
  void setLeft(Face *left);

  /*
   * Change the right face of this edge to a given face.
   * right -> the new right face of this edge;
   *          null if currently unknown
   */
  void setRight(Face *right);

  /*
   * Return the dual of this edge, directed from its right to its left.
   * <- the right to left dual of this edge;
   *    will be nonnull
   */
  Edge *Rot();

  /*
   * Return the dual of this edge, directed from its left to its right.
   * <- the left to right dual of this edge;
   *    will be nonnull
   */
  Edge *InvRot();

  /*
   * Return the edge from the destination to the origin of this edge.
   * <- the symmetric of this edge;
   *    will be nonnull
   */
  Edge *Sym();

  /*
   * Return the next ccw edge around (from) the origin of this edge.
   * <- the next edge from the origin;
   *    will be nonnull
   */
  Edge *Onext();

  /*
   * Return the next cw edge around (from) the origin of this edge.
   * <- the previous edge from the origin;
   *    will be nonnull
   */
  Edge *Oprev();

  /*
   * Return the next ccw edge around (into) the destination of this edge.
   * <- the next edge to the destination;
   *    will be nonnull
   */
  Edge *Dnext();

  /*
   * Return the next cw edge around (into) the destination of this edge.
   * <- the previous edge to the destination;
   *    will be nonnull
   */
  Edge* Dprev();

  /*
   * Return the ccw edge around the left face following this edge.
   * <- the next left face edge;
   *    will be nonnull
   */
  Edge* Lnext();

  /*
   * Return the ccw edge around the left face before this edge.
   * <- the previous left face edge;
   *    will be nonnull
   */
  Edge* Lprev();

  /*
   * Return the edge around the right face ccw following this edge.
   * <- the next right face edge;
   *    will be nonnull
   */
  Edge* Rnext();

  /*
   * Return the edge around the right face ccw before this edge.
   * <- the previous right face edge;
   *    will be nonnull
   */
  Edge* Rprev();

  /*
   * Return the length of this edge.
   * <- length;
   */
  Real length();
  
  /*
   * Returns this edge as vector (Dest-Org)
   */
  vec3 asVector();

  /* -- protected instance methods ----------------------------------------- */

  protected:

  /*
   * Initialize this edge with no connections.
   */
  Edge();

  /*
   * Release the storage occupied by the contents of this edge.
   */
  ~Edge();

  /* -- private class variables -------------------------------------------- */

  private:

  /*
   * The next unused edge ID number.
   * Positive.
   */
  static unsigned int nextID;

  /* -- private instance variables ----------------------------------------- */

  private:

  /*
   * The index of this edge in its quad-edge structure.
   * Between 0 and 3 inclusive.
   */
  unsigned int index;

  /*
   * The next ccw edge around (from) the origin of this edge.
   * Nonnull.
   */
  Edge *next;

  /*
   * The ID of this edge.
   * Positive.
   */
  unsigned int id;

  /*
   * The origin vertex of this edge, if prime.
   * Null if not prime.
   */
  Vertex *vertex;

  /*
   * The target face of this edge, if dual.
   * Null if not dual.
   */
  Face *face;

  /* -- friend classes ----------------------------------------------------- */

  friend class QuadEdge;

};

/* -- inline instance methods ---------------------------------------------- */

inline unsigned int Edge::getID()
{
  return id;
}

inline Edge* Edge::Rot()
{
  return index<3 ? this+1 : this-3;
}

inline Edge* Edge::InvRot()
{
  return index>0 ? this-1 : this+3;
}

inline Edge* Edge::Sym()
{
  return index<2 ? this+2 : this-2;
}

inline Edge* Edge::Onext()
{
  return next;
}

inline Edge* Edge::Oprev()
{
  return Rot()->Onext()->Rot();
}

inline Edge* Edge::Dnext()
{
  return Sym()->Onext()->Sym();
}

inline Edge* Edge::Dprev()
{
  return InvRot()->Onext()->InvRot();
}

inline Edge* Edge::Lnext()
{
  return InvRot()->Onext()->Rot();
}

inline Edge* Edge::Lprev()
{
  return Onext()->Sym();
}

inline Edge* Edge::Rnext()
{
  return Rot()->Onext()->InvRot();
}

inline Edge* Edge::Rprev()
{
  return Sym()->Onext();
}

inline Vertex* Edge::Org()
{
  return vertex;
}

inline Vertex* Edge::Dest()
{
  return Sym()->vertex;
}

inline Face* Edge::Left()
{
  return Rot()->face;
}

inline Face* Edge::Right()
{
  return InvRot()->face;
}

#endif /* #ifndef edgeINCLUDED */

