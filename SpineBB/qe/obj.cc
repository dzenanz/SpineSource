// obj.cc - read Wavefront "obj" file format for 3-D models
// This understands only a subset of the file format, namely
//	# comment
//	v <X> <Y> <Z>
//	f <V1> <V2> ... <VN>
//
// Paul Heckbert	10 Feb 1999


// Most of the code in this file is only relevant
// if you need to modify objReadCell

// If you can think of a more elegant way to do this, that avoids the creation
// of all these temporary data structures, and going through the model twice,
// please let me know.

#include <cstdlib>
#include <cstring>	// for strcmp
#include <fstream>	// for file I/O
#include <iomanip>	// for file input
#include <sstream>

#include <cassert>
#include <cstdio>

#include <vec3.h>
using namespace std;

#include "list.hh"
#include "array.hh"
#include "cell.hh"

// ---------------------------- some data structures used by objReadCell only

struct Tvert;

typedef List<Tvert> Vlist;	// linked list of pointers to vertices

struct Tface {			// a (temporary) face
    Vlist vlist;		// the vertices of this face, in ccw order
    int   no;                   // face number
    Face *face;                 // final face in cell, null if not inst. yet
    // need anything else in here??
};

ostream &operator<<(ostream &s, const Tface &f)	// print a Tface
    { return s << 'f' << f.no << " with vertices " << f.vlist; }

class Tsector {
public:

  Tvert *p;     // first ccw vertex
  Tface *f;     // intervening face
  Tvert *q;     // second ccw vertex

  Tsector(Tvert *p, Tface *f, Tvert *q)
    { this->p = p; this->f = f, this->q = q; }

};

typedef List<Tsector> Arc;	// arc of consec. edges emanating from a vertex
				// in counterclockwise order
				// (a linked list of pointers to other vertex)

typedef List<Arc> Arclist;	// unordered collection of arcs about a vertex
				// (a linked list of pointers to Arcs)
				// when done, this (linear) list contains
				// the ccw cycle of edges about a vertex

// For example, for vertex v below,
//
//       c------ b------i
//      / \     /      /
//     /   \   /      /
//    /     \ /      /
//   d------ v -----a--h
//    \     / \        |
//     \   /   \       |
//      \ /     \      |
//       e-------f-----g
//
// some valid Arcs are the lists (a,b), (a,b,c), (b,c), (c,d),
// (f,a), (e,f,a,b), etc. because those are the other endpoints of
// edges emanating from v, in counterclockwise (ccw) order.
// An arc always consists of at least two vertices.
// A valid Arclist is any set of disjoint arcs, in arbitrary order.
// When done, the Arclist for this vertex would be a single Arc.
// It would be a cyclic permutation of (a,b,c,d,e,f).

struct Tvert {			// a (temporary) vertex
    int no;			// ??for debugging
    int done;			// is topology fully set & arclist complete?
    vec3 p;			// position
    Arclist arclist;		// info about the vertices adjacent to this one
    Vertex *vertex;             // final vertex in cell, null if not id. yet
    int instantiated;           // true if identified and instantiated
};

ostream &operator<<(ostream &s, const Tvert &t)	// print a Tvert
    { return s << t.no; }

ostream &operator<<(ostream &s, const Tsector &sector)	// print a Tsector
    { return s << sector.p->no << "-" << sector.q->no; }

// ---------------------------- obj file input

static void merge_arc(Tvert *v, Tvert *p, Tvert *q, Tface *f) {
    // Merge the arc (p,q) into the list of arcs around vertex v.
    // Cases:
    //  1. ( bef &&  aft) it connects two existing arcs
    //  2. ( bef && !aft) it goes on the end of an existing arc
    //  3. (!bef &&  aft) it goes on the beginning of an existing arc
    //  4. (!bef && !aft) it does not connect with an existing arc
    // cout << "merge_arc " << *v << " " << *p << " " << *q << endl;
    // cout << "before, arclist=" << v->arclist;
    List_item<Arc> *a, *aft_item;
    Arc *bef = 0, *aft = 0;
    Tsector *sector = new Tsector(p, f, q);
    for (a=v->arclist.first(); a; a=a->next()) {
	// a->obj is an Arc
	if (a->obj->last()->obj->q==p) bef = a->obj;
	if (a->obj->first()->obj->p==q) {aft = a->obj; aft_item = a;}
    }
    // cout << "  bef=" << *bef << "  aft=" << *aft;
    // now concatenate the three arcs bef, (p,q), and aft
    // where bef and aft might be null
    if (bef) {
	if (aft) {	// 1. ( bef &&  aft) it connects two existing arcs
	    bef->append(sector);		// insert new sector
	    if (bef==aft) {
		// done with vertex! connecting these would make arc circular
		// cout << v->arclist << " done" << endl;
		v->done = 1;
		return;
	    }
	    // now we'll merge two arcs in the arclist
	    v->arclist.remove(aft_item);	// remove following arc
	    bef->concat(aft);			// and concat it into previous
	}
	else		// 2. ( bef && !aft) it goes on the end of existing arc
	    bef->append(sector);
    }
    else {
	if (aft)	// 3. (!bef &&  aft) it goes on beg. of existing arc
	    aft->prepend(sector);
	else {		// 4. (!bef && !aft) it doesn't connect w. existing arc
	    Arc *arc = new Arc;
	    assert(arc);
	    arc->append(sector);
	    v->arclist.append(arc);
	}
    }
    // cout << "after, arclist=" << v->arclist;
}

static void add_arcs(Vlist &vlist, Tface *f) {
    // cout << "add_arcs " << vlist;
    // vlist is not a circular list, but we need to step through all
    // consecutive triples as if it were circular
    List_item<Tvert> *u, *v, *w;
    for (u=vlist.last(), v=vlist.first(), w=v->next(); w; u=v, v=w, w=w->next())
	merge_arc(v->obj, w->obj, u->obj, f);
    merge_arc(v->obj, vlist.first()->obj, u->obj, f);  // one more that we missed
}

/*
 * identified   <=> Tvert has been associated with a particular Vertex
 * instantiated <=> Tface has been associated with a particular Face AND
 *                  all vertices of the face have been identified
 * instantiated <=> Tvert has been identified AND
 *                  all adjacent Tfaces have been instantiated
 */

/*
 * Return true if a given pair of vertices is connected directly by an edge
 * along a given left face.
 * vertex1, vertex2 -> the vertices to check;
 *                     must be nonnull
 * left             -> the left face to check for;
 *                     must be nonnull
 * <- true if there is an edge from _vertex1_ to _vertex2_ with left face
 *    _left_
 */
static int isConnected(Vertex *vertex1, Vertex *vertex2, Face *left)
{
  assert(vertex1!=0);
  assert(vertex2!=0);
  assert(left!=0);

  // check the orbit of vertex1 for an edge to vertex2

  VertexEdgeIterator edges(vertex1);

  Edge *edge;

  while ((edge = edges.next())!=0)
    if (edge->Dest()==vertex2 && edge->Left()==left)
      return 1;

  return 0;
}

/*
 * Return the face to the right of a given face around a given vertex.
 * vertex -> the vertex to look for the face around;
 *           must be nonnull
 * left   -> the left face to return the right face of;
 *           must be nonnull
 * <- the face to the right of _left_ around _vertex_;
 *    null if none
 */
static Face *RightFace(Vertex *vertex, Face *left)
{
  assert(vertex!=0);
  assert(left!=0);

  // check the left face of each edge in the orbit of the vertex

  Edge *start = vertex->getEdge();
  Edge *scan  = start;

  do
  {
    if (scan->Left()==left)
      return scan->Right();

    scan = scan->Onext();
  }
  while (scan!=start);

  return 0;
}

/*
 * Return true if a given vertex is adjacent to a given face.
 * face   -> the face to look for the vertex in;
 *           must be nonnull
 * vertex -> the vertex to look for;
 *           must be nonnull
 * <- true if _vertex_ is on _face_
 */
static int hasVertex(Face *face, Vertex *vertex)
{
  assert(face!=0);
  assert(vertex!=0);

  // check the origin vertex of each edge on the face

  FaceEdgeIterator edges(face);

  Edge *edge;

  while ((edge = edges.next())!=0)
    if (edge->Org()==vertex)
      return 1;

  return 0;
}

/*
 * Return true if a given face includes all the identified vertices on a given
 * Tvert list.
 * face  -> the face to check;
 *          must be nonnull
 * vlist -> the vertex list to check against;
 *          must be nonnull
 * <- true if _face_ is adjacent to all the vertices on _vlist_
 */
static int hasVertices(Face *face, Vlist *vlist)
{
  assert(face!=0);
  assert(vlist!=0);

  // check each vertex on the list

  for (List_item<Tvert> *vi = vlist->first(); vi!=0; vi = vi->next())
  {
    Vertex *vertex = vi->obj->vertex;

    if (vertex!=0 && !hasVertex(face, vertex))
      return 0;
  }

  return 1;
}

/*
 * Return a face that can be used to instantiate a given Tface.
 * cell -> the cell to get the face from;
 *         must be nonnull
 * f    -> Tface to get the face for;
 *         must be nonnull
 * <- a face that can be used to instantiate _f_;
 *    null if none are available
 */
static Face *getFace(Cell *cell, Tface *f)
{
  assert(cell!=0);
  assert(f!=0);

  // locate all the unused faces in the cell

  Face       **faces = new Face*[cell->countFaces()];
  unsigned int count = 0;

  {
    CellFaceIterator iterator(cell);

    Face *face;

    while ((face = iterator.next())!=0)
      if (face->data==0)
	faces[count++] = face;
  }

  // discard any faces that don't include all the identified vertices of the
  // Tface

  {
    unsigned int i = 0;

    while (i<count)
    {
      Face *face = faces[i];

      if (hasVertices(face, &f->vlist))
	i++;
      else
	faces[i] = faces[--count];
    }
  }

  Face *face = count>0 ? faces[0] : 0;

  delete[] faces;

  return face;
}

/*
 * Instantiate a given Tface in a given cell by identifying its vertices.
 * cell -> the cell to instantiate the face in;
 *         must be nonnull
 * f    -> the Tface to instantiate;
 *         must be nonnull
 */
static void makeFace(Cell *cell, Tface *f)
{
  assert(cell!=0);
  assert(f!=0);

  // get the face to use for the Tface

  Face *face = getFace(cell, f);

  assert(face!=0);

  // connect all pairs of identified vertices on the face, as necessary

  {
    for (List_item<Tvert> *vi = f->vlist.first(); vi!=0; vi = vi->next())
    {
      Vertex *vertex1 = vi->obj->vertex;
      Vertex *vertex2;

      if (vertex1!=0)
      {
	// find the next identified vertex, even if just itself

	List_item<Tvert> *vj = vi;

	for (;;)
	{
	  vj = vj->next();

	  if (vj==0)
	    vj = f->vlist.first();

	  vertex2 = vj->obj->vertex;

	  if (vertex2!=0)
	    break;
	}

	// connect the vertices, if necessary

	if (!isConnected(vertex1, vertex2, face))
	  (void)cell->makeFaceEdge(face, vertex1, vertex2)->Right();
      }
    }
  }

  // find the first identified vertex

  List_item<Tvert> *vi0 = f->vlist.first();

  while (vi0->obj->vertex==0)
    vi0 = vi0->next();

  // identify all the following and preceding vertices

  List_item<Tvert> *vi     = vi0;
  Vertex           *vertex = vi0->obj->vertex;

  for (;;)
  {
    vi = vi->next();

    if (vi==0)
      vi = f->vlist.first();

    if (vi==vi0)
      break;

    Tvert *v = vi->obj;

    if (v->vertex==0)
    {
      Face *right = RightFace(vertex, face);

      assert(right!=0);

      v->vertex      = cell->makeVertexEdge(vertex, face, right)->Dest();
      v->vertex->pos = v->p;

      v->vertex->setID(v->no);
    }

    vertex = v->vertex;
  }

  // the face is now instantiated

  f->face = face;

  face->setID(f->no);
  face->data = f;
}

/*
 * Instantiate a given identified Tvert in a given cell by instantiating its
 * adjacent faces.
 * cell -> the cell to instantiate the Tvert in;
 *         must be nonnull
 * v    -> the Tvert to instantiate;
 *         must be nonnull
 */
static void makeVertex(Cell *cell, Tvert *v)
{
  assert(cell!=0);
  assert(v!=0);

  // find the first sector with an identified p vertex

  List_item<Tsector> *wi0 = v->arclist.first()->obj->first();

  while (wi0->obj->p->vertex==0)
    wi0 = wi0->next();

  // instantiate all following sectors of the vertex in counterclockwise order

  List_item<Tsector> *wi = wi0;

  do
  {
    Tface *f = wi->obj->f;

    if (f->face==0)
      makeFace(cell, f);

    wi = wi->next();

    if (wi==0)
      wi = v->arclist.first()->obj->first();
  }
  while (wi!=wi0);

  // the vertex is now instantiated

  v->instantiated = 1;
}

static void print_quadedge(Array<Tvert> verts, List<Tface> faces) {
    // print vertices around each face and vertex currently

    cout << "VERTICES OF EACH FACE:" << endl;
    List_item<Tface> *fi;
    for (fi=faces.first(); fi; fi=fi->next()) {
	cout << "face:";
	List_item<Tvert> *vi;
	for (vi=fi->obj->vlist.first(); vi; vi=vi->next())
	    cout << " " << vi->obj->no;
	cout << endl;
    }
    cout << endl;

    cout << "VERTICES AROUND EACH VERTEX:" << endl;
    int i;
    for (i=0; i<verts.num(); i++) {
	Tvert *v;
	v = &verts[i];
	cout << "around vertex " << v->no << ":";
	assert(v->done);
	assert(v->arclist.length()==1);
	// step through the Tverts in the first (and only) arc of arclist
	List_item<Tsector> *wi;
	for (wi=v->arclist.first()->obj->first(); wi; wi=wi->next()) {
	    Tsector *sector = wi->obj;
	    cout << " " << *sector;
	}
	cout << endl;
    }
}

void check_closed(Array<Tvert> verts) {
    // check to see if polyhedron is closed (else we'd crash soon anyway)
    int i;
    for (i=0; i<verts.num(); i++) {
	Arclist &al = verts[i].arclist;
	if (!verts[i].done || al.length()!=1) {
	    if (al.length()==0)
		cerr << "\nERROR in OBJ file: unused vertex "
		    << verts[i].no << endl;
	    else if (!verts[i].done)
		cerr << "\nERROR in OBJ file: vertex " << verts[i].no
		    << " is not surrounded by polygons" << endl;
	    else
		cerr << "\nERROR in OBJ file: repeated face: " <<
			*al.first()->next()->obj->first()->obj->f << endl;
	    exit(1);
	}
    }
}

static Cell *build_quadedge(Array<Tvert> verts, List<Tface> faces) {

  check_closed(verts);

  // create a cell and fetch its initial vertex

  Cell *cell = Cell::make();

  Vertex *vertex1;

  {
    CellVertexIterator vertices(cell);

    vertex1 = vertices.next();
  }

  // instantiate a face of the initial vertex

  {
    Tvert *v = &verts[0];

    v->vertex      = vertex1;
    v->vertex->pos = v->p;

    v->vertex->setID(v->no);

    makeFace(cell, v->arclist.first()->obj->first()->obj->f);
  }

  // instantiate identified vertices until all are instantiated

  for (;;)
  {
    int instantiated = 1;

    for (int i = 0; i<verts.num(); i++)
    {
      Tvert *v = &verts[i];

      if (v->vertex!=0 && !v->instantiated)
	makeVertex(cell, v);

      instantiated &= v->instantiated;
    }

    if (instantiated)
      break;
  }

  // reset the data pointers of all faces

  {
    CellFaceIterator iterator(cell);

    Face *face;

    while ((face = iterator.next())!=0)
      face->data = 0;
  }

  return cell;
}

static Cell *objReadCell(istream &s, const char *streamname) {

    // warning: this routine does a lousy job of checking for errors
    // (e.g. spaces at the end of input lines)
    // that should be fixed!

    char tok[20];

    Array<Tvert> verts;		// all the vertices
    int nvert = 0;		// current vertex number (counts up)
    int nface = 0;		// current face number (counts up)

    List<Tface> faces;		// all the faces

    while (s >> setw(sizeof tok) >> tok) {
        // cout << "(" << tok << ")" << endl;
	if (!strcmp(tok, "v")) {		// vertex command
	    nvert++;
	    float x, y, z;
	    s >> x >> y >> z;
	    // cout << "verts[" << nvert << "]=" << Vec3(x, y, z) << endl;
	    verts[nvert-1].p = vec3(x, y, z);
	    verts[nvert-1].no = nvert;
	    verts[nvert-1].done = 0;
	    verts[nvert-1].vertex = 0;
	    verts[nvert-1].instantiated = 0;
	    // OBJ vertex numbering starts at 1, but we shift them
	    // all to start at 0 for indexing into vert[]
	}
	else if (!strcmp(tok, "f")) {	// face command
	    nface++;
	    Tface *f = new Tface;
	    f->face = 0;
	    f->no   = nface;
	    assert(f);
	    faces.append(f);
	    int n;
	    for (n=0; s.peek()!='\n'; n++) {
		int v;
		// cout << "  peek='" << (char)s.peek() << "' ";
		s >> v;
		// cout << "  got " << v << endl;
		v--;		// we start numbering at 0, not 1
		f->vlist.append(&verts[v]);
	    }
	    // cout << "done gobbling" << endl;
	    // cout << "f " << f->vlist;
	    add_arcs(f->vlist, f);		// add the topological info in face f
	}
	else if (!strcmp(tok, "#")) {
	    // gobble comment
	    s.ignore(1000, '\n');
	}
	else {
	    cerr << "objReadCell: I can't parse this OBJ file, hit token (" << tok
		<< ")" << endl;
	    exit(1);
	}
    }
    /*
    cout << "obj file " << streamname << " contained "
	<< nvert << " vertices, "
	<< faces.length() << " faces "
	<< endl;
    */

    // print_quadedge(verts, faces);

    return build_quadedge(verts, faces);
}

Cell *objReadCell(const char *file) {
    ifstream s(file, ios::in);
    if (!s) {
	cerr << "objReadCell: can't read " << file << endl;
	return 0;
    }
    return objReadCell(s, file);
}

Cell *objReadCellFromString(string obj_file_contents)
{
	stringstream ss(ios_base::in);
	ss.str(obj_file_contents);
	return objReadCell(ss, "string");
}

#ifdef OBJ_MAIN	// compile with -DOBJ_MAIN to test objReadCell

void main(int argc, char **argv) {
    if (argc!=2) exit(1);
    objReadCell(argv[1]);
}

#endif

static void objWriteCell(Cell *cell, ostream &s, const char *streamname)
{
  // renumber vertices in current order
  {
    CellVertexIterator vertices(cell);

    Vertex      *vertex;
    unsigned int id = 1;

    while ((vertex = vertices.next())!=0)
      vertex->setID(id++);
  }

  // write vertices in same order: explicit ids become implicit ids

  s << "# " << cell->countVertices() << " vertices\n";

  {
    CellVertexIterator vertices(cell);

    Vertex *vertex;

    while ((vertex = vertices.next())!=0)
      s << "v " << vertex->pos[0] << " "
	        << vertex->pos[1] << " "
	        << vertex->pos[2] << "\n";
  }

  // write faces in any order

  s << "# " << cell->countFaces() << " faces\n";

  {
    CellFaceIterator faces(cell);

    Face *face;

    while ((face = faces.next())!=0)
    {
      s << "f";

      FaceEdgeIterator edges(face);

      Edge *edge;

      while ((edge = edges.next())!=0)
		s << " " << edge->Org()->getID();

      s <<endl;
    }
  }
}

void objWriteCell(Cell *cell, const char *file)
{
    ofstream s(file, ios::out);
    if (!s) {
	cerr << "objWriteCell: can't write " << file << endl;
	return;
    }
    objWriteCell(cell, s, file);
}

Cell *objCopyCell(Cell *cell)
{
  // renumber vertices in current order
  // yuk: should really leave ids intact ???

  {
    CellVertexIterator vertices(cell);

    Vertex      *vertex;
    unsigned int id = 1;

    while ((vertex = vertices.next())!=0)
      vertex->setID(id++);
  }

  Array<Tvert> verts;		// all the vertices
  int nvert = 0;		// current vertex number (counts up)

  List<Tface> faces;		// all the faces
  int nface = 0;		// current face number (counts up)

  // make the Tverts

  {
    CellVertexIterator vertices(cell);

    Vertex *vertex;

    while ((vertex = vertices.next())!=0)
    {
      nvert++;
      verts[nvert-1].p = vertex->pos;
      verts[nvert-1].no = nvert;
      verts[nvert-1].done = 0;
      verts[nvert-1].vertex = 0;
      verts[nvert-1].instantiated = 0;
      // OBJ vertex numbering starts at 1, but we shift them
      // all to start at 0 for indexing into vert[]
    }
  }

  // make the Tfaces

  {
    CellFaceIterator iterator(cell);

    Face *face;

    while ((face = iterator.next())!=0)
    {
      nface++;
      Tface *f = new Tface;
      f->face = 0;
      f->no   = nface;
      assert(f);
      faces.append(f);
      {
	FaceEdgeIterator edges(face);

	Edge *edge;

	while ((edge = edges.next())!=0)
	  f->vlist.append(&verts[edge->Org()->getID()-1]);
      }
      add_arcs(f->vlist, f);		// add the topological info in face f
    }
  }

  // go to town

  return build_quadedge(verts, faces);
}

Cell *objCloneCell(Cell *cell)
{
  return objCopyCell(cell); //implement properly later
}
