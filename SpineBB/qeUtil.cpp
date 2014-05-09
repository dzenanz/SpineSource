#include "declarations.h"
#include "itkTriangleCell.h"
#include <fstream>
#include <qmessagebox.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTransform.h>
#include <vtkSTLWriter.h>

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
	enew->Dest()->stepSize=(enew->Org()->stepSize + e->Dest()->stepSize)*0.5;
	enew->Dest()->lastIntensity=(enew->Org()->lastIntensity + e->Dest()->lastIntensity)/2;
	enew->Dest()->lastH=(enew->Org()->lastH + e->Dest()->lastH)/2;

    return enew;	// return new edge
}

vec3 normal(Vertex *v)
{
	Edge *e=v->getEdge();
	vec3 n(0,0,0),n0;
	do
	{
		//n=sum(n0*angle/(|a|+|b|)),  n0=a^b/(|a|*|b|)
        vec3 a=e->asVector(), b=e->Onext()->asVector();
		
		double cosFi=a*b/(a.length()*b.length());
		if (cosFi>=1)
			cosFi=0.9999;
		else if (cosFi<=-1)
			cosFi=-0.9999;

		n0=a^b/(a.length()*b.length());
		n+=n0*acos(cosFi)/(a.length()+b.length());

		e=e->Onext();
	} while (e!=v->getEdge());
	n.normalize();
	return n;
}

float curvature(Vertex *v)
{
	//sum of all angles around this vertex are equal to 180 degrees on plain surfaces
	//that sum is lower on protrusions and dents, and higher on rims of those
	Edge *e=v->getEdge(), *next=e->Onext();
	double angle=0;
	vec3 n, c=e->asVector();
	do
	{
		n=next->asVector();
		
		double cosFi=(c*n)/(c.length()*n.length());
		if (cosFi>1)
			cosFi=1;
		else if (cosFi<-1)
			cosFi=-1;
		angle+=acos(cosFi);

		e=next;
		c=n;
		next=next->Onext();
	} while (e!=v->getEdge());
     
	return abs(6.283185307179586476925286766559-angle); //2pi-angle
}

void updateNormalsAndCurvatures(Cell *qe)
{
	CellVertexIterator it(qe);
	Vertex *v;
	while ((v = it.next()) != 0)
	{
		v->normal=normal(v);
	}

	CellVertexIterator it2(qe);
	while ((v = it2.next()) != 0)
	{
		v->curv=curvature(v);
	}
}

void updateVertexDistances(Cell *qe)
{
	CellVertexIterator it(qe);
	Vertex *v;
	while ((v = it.next()) != 0)
		v->dist=0;

    CellVertexIterator it2(qe);
	while ((v = it2.next()) != 0)
	{
		VertexEdgeIterator vei(v);
        Edge *e;
        v->dist*=-1; //any contributions to this vertex by processed vertices has to be inverted
        while ((e = vei.next()) != 0)
            if (e->Dest()->dist<=0) //vertex unprocessed
	        {
                float len=e->asVector().length();
                v->dist+=len;
                e->Dest()->dist-=len; //opposite vertices (unprocessed) get negative contributions
	        }
        v->dist/=v->valence;
	}
}

vec3 supposedMoveDirection(Vertex *v, float avgEdgeLen)
//calculated direction in which vertex is supposed to move
//in order to have a more uniform surface distribution
{
    vec3 smd(0,0,0);
    Edge *e;
    VertexEdgeIterator vei(v);
    while ((e = vei.next()) != 0)
    {
        vec3 edgeVec=e->asVector();
        float len=edgeVec.normalize();
        smd+=edgeVec*(len-avgEdgeLen);
    }
    smd/=v->valence;
    return smd;
}

void rearrangeVertices(Cell *qe, float avgEdgeLen, const vec3& center)
{
    vec3 smd(0,0,0);
    Vertex *v;
    CellVertexIterator it2(qe);
	while ((v = it2.next()) != 0)
        if (abs(v->dist-avgEdgeLen)/avgEdgeLen>0.05) //difference from average >5%
        {
            vec3 cvv=v->pos-center;
            vec3 smd=supposedMoveDirection(v, avgEdgeLen);
            //now we need to get rid of component along center-vertex vector (inflation direction)
            //move vertex only perpendicular to this vector
            vec3 planeNormal=cvv^smd;
            vec3 allowedMoveDirection=planeNormal^cvv;
            //vector perpendicular to cvv and in the same plane as supposedMoveDirection
            allowedMoveDirection.normalize();
            v->pos+=allowedMoveDirection*(allowedMoveDirection*smd);
        }
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

float getAverageEdgeLength(Cell *qe, float &minLen, float &maxLen, float &stdDev)
//returns average edge length
{
	double sum=0;
	maxLen=0;
	minLen=std::numeric_limits<float>::max();
	float len, avg;
	unsigned count=0;
	Vertex *v;
	Edge *e;

	CellVertexIterator cvi(qe);
	while ((v = cvi.next())!=0)
	{
		VertexEdgeIterator vei(v);
		Edge *e;
		while ((e = vei.next())!=0)
			if (e->Dest()<v) //account for every edge only once
			{
				len=e->length();
				sum+=len;
				count++;
				if (len<minLen)
					minLen=len;
				if (len>maxLen)
					maxLen=len;
			}
	}
	avg=sum/count;

	sum=0;
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

void writeSTL(Cell *qe, const char * fileName)
{
    vtkSTLWriter *writer=vtkSTLWriter::New();
    vtkPolyData *poly=qe2vtk(qe);
    writer->SetInputData(poly);
    writer->SetFileName(fileName);
    writer->Write();
    poly->Delete();
}

void applyMatrix(Cell *qe, vtkMatrix4x4 *matrix)
{
    CellVertexIterator it3(qe);
    Vertex *v;
    while ((v = it3.next()) != 0)
        v->pos=matrix*v->pos;
}