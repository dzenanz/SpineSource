#include "declarations.h"
#include "itkTriangleCell.h"
#include <fstream>
#include <assert.h>
#include <qmessagebox.h>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTransform.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSTLWriter.h>
#include <vtkSmartPointer.h>

typedef MeshType::CellType CellType;
typedef CellType::CellAutoPointer CellAutoPointer;
typedef itk::TriangleCell< CellType > TriangleType;
typedef MeshType::PointsContainerIterator  PointIterator;
typedef MeshType::CellsContainerIterator  CellIterator;

unsigned int exp2i(int exponent) //2^exp
{
	unsigned int r=1;
	for (int i=1; i<=exponent; i++)
		r*=2;
	return r;
}

Edge *split(Edge *e)
{
    // split edge e, putting new vertex at midpoint

    // get Cell pointer from vertex (Edges don't have one)
    Cell *c = e->Org()->getCell();

    // split, creating new edge and vertex (sets topology only)
    Edge *enew = c->makeVertexEdge(e->Org(), e->Left(), e->Right());

    // At this point enew->Dest()==e->Org(),
    // and enew->Dest(), the new vertex, is between enew and e.
    // You might want to check the definition of makeVertexEdge to
    // convince yourself of this.

    // position new vertex at midpoint (note use of Vec3::operator+)
    enew->Dest()->pos =(enew->Org()->pos + e->Dest()->pos)*0.5;
	enew->Dest()->subdivPos =(enew->Org()->subdivPos + e->Dest()->subdivPos)*0.5;

    return enew;	// return new edge
}

void translate(Cell *qe, const vec3& vector)
{
	CellVertexIterator cellVertices(qe);
	Vertex *v;
	while ((v = cellVertices.next()) != 0)
	{
		v->pos+=vector;
	}
}

MeshType::Pointer qe2itkMesh(Cell *cell)
//copied from objWriteCell
{
	MeshType::Pointer mesh=MeshType::New();
	MeshType::PointType p;

	// renumber vertices in current order
	CellVertexIterator vertices(cell);
	Vertex *vertex;
	unsigned int id = 0;

	while ((vertex = vertices.next())!=0)
	{
		p[0]=vertex->pos[0];
		p[1]=vertex->pos[1];
		p[2]=vertex->pos[2];
		mesh->SetPoint( id, p );
		vertex->setID(++id);
	}
	
    CellFaceIterator faces(cell);
    Face *face;
	id=0;

    while ((face = faces.next())!=0)
    {
		CellAutoPointer tr;
		tr.TakeOwnership(new TriangleType);
		FaceEdgeIterator edges(face);
		Edge *edge;
		
		unsigned int pid=0;
		while ((edge = edges.next())!=0)
		{
			tr->SetPointId(pid++, edge->Org()->getID()-1);
		}
		mesh->SetCell(id, tr);
		mesh->SetCellData(id, 0.0);
		face->setID(++id);
    }

	return mesh;
}

//converts vertex positions from physical to index space of an image
void qePhysicalToIndex(Cell *qe, VisualizingImageType::Pointer image)
{
	CellVertexIterator vertices(qe);
	Vertex *vertex;
	VisualizingImageType::PointType p;
	itk::ContinuousIndex<double, 3> ind;

	while ((vertex = vertices.next())!=0)
	{
		p[0]=vertex->pos[0];
		p[1]=vertex->pos[1];
		p[2]=vertex->pos[2];
		if (image->TransformPhysicalPointToContinuousIndex(p, ind))
		{
			vertex->pos[0]=ind[0];
			vertex->pos[1]=ind[1];
			vertex->pos[2]=ind[2];
		}
		else
			throw itk::ExceptionObject(__FILE__,__LINE__,"Point outside of volume!",__FUNCTION__);
	}
}

//converts vertex positions from index to physical space of an image
void qeIndexToPhysical(Cell *qe, VisualizingImageType::Pointer image)
{
	CellVertexIterator vertices(qe);
	Vertex *vertex;
	VisualizingImageType::PointType p;
	itk::ContinuousIndex<double, 3> ind;

	while ((vertex = vertices.next())!=0)
	{
		ind[0]=vertex->pos[0];
		ind[1]=vertex->pos[1];
		ind[2]=vertex->pos[2];
		image->TransformContinuousIndexToPhysicalPoint(ind, p);
		vertex->pos[0]=p[0];
		vertex->pos[1]=p[1];
		vertex->pos[2]=p[2];
	}
}

//converts vertex positions from index to physical space of an image
void itkMeshIndexToPhysical(MeshType::Pointer mesh, VisualizingImageType::Pointer image)
{
	PointIterator pointIterator = mesh->GetPoints()->Begin();
	PointIterator pointEnd      = mesh->GetPoints()->End();
	VisualizingImageType::PointType p;
	itk::ContinuousIndex<double, 3> ind;

	while( pointIterator != pointEnd )
	{
		ind[0]=pointIterator.Value()[0];
		ind[1]=pointIterator.Value()[1];
		ind[2]=pointIterator.Value()[2];
		image->TransformContinuousIndexToPhysicalPoint(ind, p);
		pointIterator.Value()=p;
		++pointIterator;
	}
}

void writeMeshObj(MeshType::Pointer mesh, const char * fileName)
{
	std::ofstream f(fileName);

	f<<"# "<<mesh->GetNumberOfPoints()<<" vertices\n";
	PointIterator pointIterator = mesh->GetPoints()->Begin();
	PointIterator pointEnd      = mesh->GetPoints()->End();
	while( pointIterator != pointEnd )
	{
		f<<"v "<<pointIterator.Value()[0]<<' '<<pointIterator.Value()[1]<<' '<<pointIterator.Value()[2]<< std::endl;
		++pointIterator;
	}

	f<<"# "<<mesh->GetNumberOfCells()<<" faces\n";	
	CellIterator cellIterator = mesh->GetCells()->Begin();
	CellIterator cellEnd      = mesh->GetCells()->End();
	while( cellIterator != cellEnd )
	{
		CellType * cell = cellIterator.Value();
		if ( cell->GetType() == CellType::TRIANGLE_CELL)
		{
			TriangleType * t=(TriangleType *)cell;
			TriangleType::PointIdIterator pit = t->PointIdsBegin();
			f<<"f "<<(*pit+1);
			pit++;
			f<<' '<<(*pit+1);
			pit++;
			f<<' '<<(*pit+1)<<'\n';
		}
		++cellIterator;
	}

	f.close();
}

void calcFaceNormals(Cell *s)
{
	CellFaceIterator fi(s);
	Face *f;
	while ((f = fi.next()) != 0)
	{
		Edge *e=f->getEdge();
		vec3 a=e->asVector(), b=e->Lnext()->asVector();
		f->normal=a^b*0.5; //normal length equal to surface area of this triangle face
	}
}

vtkMatrix4x4* rotate(vec3 from, vec3 to)
{
    vtkMatrix4x4 *m=vtkMatrix4x4::New();
    vtkSmartPointer<vtkTransform> rot=vtkSmartPointer<vtkTransform>::New();
    vec3 rotAxis=from^to; //cross product is rotation axis
    double angle=acos(from*to/(from.length()*to.length())); //angle from dot product
    rot->RotateWXYZ(angle*180/M_PI, rotAxis.ptr());
    m->DeepCopy(rot->GetMatrix());
    return m;
}

void rotate(Cell *qe, vec3 from, vec3 to)
{
    vtkSmartPointer<vtkMatrix4x4> rot=rotate(from, to);
    applyMatrix(qe, rot);
}

//returns average edge length (also updates average edge length for each vertex)
float getAverageEdgeLength(Cell *qe, float &minLen, float &maxLen, float &stdDev)
{
	double sum=0;
	maxLen=0;
	minLen=std::numeric_limits<float>::max();
	float len, sum1, avg;
	unsigned count=0;
	Vertex *v;
	Edge *e;

	CellVertexIterator cvi(qe);
	while ((v = cvi.next())!=0)
	{
		VertexEdgeIterator vei(v);
		Edge *e;
        sum1=0;
		while ((e = vei.next())!=0)
		{
			len=e->length();
			sum1+=len;
			if (len<minLen)
				minLen=len;
			if (len>maxLen)
				maxLen=len;
		}
        sum+=sum1;
        count+=v->valence;
        //v->dist=sum1/v->valence; //average edge length in this vertex's neigborhood
	}
	avg=sum/count;

	sum=0;
    count/=2; //we now count every edge only from 1 direction
	cvi=CellVertexIterator(qe);
	while ((v = cvi.next())!=0)
	{
		VertexEdgeIterator vei(v);
		while ((e = vei.next())!=0)
			if (e->Dest()<v) //account for every edge only once
				sum+=(e->length()-avg)*(e->length()-avg);
	}
	stdDev=sqrt(sum/count);
	return avg;
}

void calcVertexValence(Cell *qe, bool showResult)
{
	assert(qe);
	CellVertexIterator cvi(qe);
	Vertex *v;
	typedef std::map<unsigned char, unsigned int> vmap;
	vmap valence;

	while ( (v=cvi.next())!=0 )
	{
		Edge *orig=v->getEdge();
		unsigned int count=1;
		for (Edge *e=orig->Onext(); e!=orig; e=e->Onext())
			count++;
		v->valence=count;
		valence[count]++;
	}
	if (showResult)
	{
		QString s;
		for (vmap::iterator it=valence.begin(); it!=valence.end(); it++)
			s+="valence["+QString::number(it->first)+"]="+QString::number(it->second)+"\n";
		QMessageBox::warning(0, "Vertex valence", s);
	}
}

vtkActor *poly2actor(vtkPolyData *poly, bool wireframe, vec3 color)
{
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);
    vtkActor *actor=vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color[0],color[1],color[2]);
    if (wireframe)
        actor->GetProperty()->SetRepresentationToWireframe();
    return actor;
}

void writeSTL(vtkPolyData* poly, const char * fileName)
{
    vtkSmartPointer<vtkSTLWriter> writer=vtkSmartPointer<vtkSTLWriter>::New();
    writer->SetInputData(poly);
    writer->SetFileName(fileName);
    writer->Write();
}

void applyMatrix(Cell *qe, vtkMatrix4x4 *matrix)
{
    CellVertexIterator it3(qe);
    Vertex *v;
    while ((v = it3.next()) != 0)
        v->pos=matrix*v->pos;
}