#include "MainLogic.h"

#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "itksys/SystemTools.hxx"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkLineIterator.h"
#include "itkGradientImageFilter.h"
#include "itkTriangleMeshToBinaryImageFilter.h"

#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkOBJReader.h>
#include "vtkOBJWriter.h"
#include <vtkRendererSource.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkArrowSource.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkClipPolyData.h>
#include <vtkFillHolesFilter.h>
#include <vtkClipClosedSurface.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkBox.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkCleanPolyData.h>

#include <vec3.h>

#include <QPainter>
#include <QFile>
#include <QDir>
#include <qmessagebox.h>

#include <Eigen/Eigen>

using namespace std;

MainLogic::MainLogic(MainWindow& mainWindow)
:mainForm(mainWindow)
{
	connect(&mainForm, SIGNAL(volumeOpened(std::string)), SLOT(volumeOpen(std::string)));
	connect(mainForm.painter, SIGNAL(majorUpdate(int, int)), SLOT(majorSliceUpdate(int, int)));
	current=0;
	visualizing=0;

    if (QDir("D:/Temp").exists())
        tempDir="D:/Temp/";
    else if (QDir("C:/Temp").exists())
        tempDir="C:/Temp/";
    else
    {
        tempDir=QDir::tempPath().toStdString();
        //add trailing slash if needed
        if (strcmp(&tempDir.at(tempDir.size() - 1),"\\") && strcmp(&tempDir.at(tempDir.size() - 1),"/"))
            tempDir+='/';
    }
}

std::string MainLogic::tempDir;
//disease thresholds:
double crushedVertebraThreshold=0.2;
double significantSpondylolisthesis=0.25;

template <class T>
inline bool withinEpsilonEnvironment(T point, T x, float eps)
{
	return x>=point-point*eps && x<=point+point*eps;
}

void MainLogic::volumeOpen(std::string filename)
{
	try {
	mainForm.statusbar->showMessage("Reading dataset: "+QString::fromStdString(filename));
	typedef  itk::ImageFileReader< InternalImageType >  ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( filename.c_str() );
	reader->Update();
	current=reader->GetOutput();

    typedef itk::MinimumMaximumImageCalculator< InternalImageType >  CalculatorType;
	CalculatorType::Pointer calculatorI = CalculatorType::New();
	calculatorI->SetImage( current );
	calculatorI->ComputeMaximum();
	maxValue = calculatorI->GetMaximum();

    int i=1;
    while (pow(2.0f,i)<maxValue)
        i++;
    maxValue2=pow(2.0f,i);

	lImage=InternalImageType::New();
	hImage=InternalImageType::New();
	std::string fnNoExt=itksys::SystemTools::GetFilenamePath(filename)+'/'
		+itksys::SystemTools::GetFilenameWithoutLastExtension(filename);

    mainForm.statusbar->showMessage("Calculating LH values");
    calcLHvalues(current, maxValue/100, false, lImage, hImage); //epsilon = 1%

    if (mainForm.actionSave_LH_images_and_histogram->isChecked())
    {
	    CalculatorType::Pointer calculatorL = CalculatorType::New();
	    calculatorL->SetImage( lImage );
	    calculatorL->ComputeMaximum();
	    maxValue = max(calculatorL->GetMaximum(), maxValue);
    	
	    CalculatorType::Pointer calculatorH = CalculatorType::New();
	    calculatorH->SetImage( hImage );
	    calculatorH->ComputeMaximum();
	    maxValue = max(calculatorH->GetMaximum(), maxValue);
        while (pow(2.0f,i)<maxValue)
            i++;
        maxValue2=pow(2.0f,i);
    	
        VisualizingImageType::Pointer lVis, hVis;

	    typedef itk::ShiftScaleImageFilter < InternalImageType, VisualizingImageType > RescaleImageFilterType;
	    RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
	    rescale->SetInput( lImage );
	    rescale->SetScale(256/maxValue2);
        rescale->Update();
	    lVis=rescale->GetOutput();

	    RescaleImageFilterType::Pointer rescale2 = RescaleImageFilterType::New();
	    rescale2->SetInput( hImage );
	    rescale2->SetScale(256/maxValue2);
        rescale2->Update();
	    hVis=rescale2->GetOutput();
    	
	    //write LH files
	    typedef itk::ImageFileWriter<InternalImageType> WriterType;
	    WriterType::Pointer writer1=WriterType::New();
	    writer1->SetFileName(fnNoExt+"_L.mha");
	    writer1->SetInput(lImage);
	    writer1->Update();

	    writer1->SetFileName(fnNoExt+"_H.mha");
        writer1->SetInput(hImage);
	    writer1->Update();

	    calc2DJointHistogram(lVis, hVis, fnNoExt+"_LH.png");
    }

	mainForm.updateVisualization();
	mainForm.updateOverlay();
	volume_filename=filename;
	mainForm.setWindowTitle(QString("Spine Analyzer - ")+QString::fromStdString(filename));
    QApplication::processEvents();
	}   catch (std::exception& e)
	{
		QMessageBox mb(QMessageBox::Critical, QString("Error opening file"),
			QString("Error opening file ")+QString::fromStdString(filename)+"\nException content:\n"+e.what(),
            QMessageBox::Ok);
		mb.exec();
	}
}

void writeImage(InternalImageType::Pointer image, std::string fileName, bool compressed)
{
    typedef itk::ImageFileWriter<InternalImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    writer->SetUseCompression(compressed);
    writer->Update();
}

void writeImage(VisualizingImageType::Pointer image, std::string fileName, bool compressed)
{
    typedef itk::ImageFileWriter<VisualizingImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    writer->SetUseCompression(compressed);
    writer->Update();
}

inline double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void savePath(QPainterPath &path, ofstream& f)
{
	for (int i=0; i<path.elementCount(); i++)
		f<<path.elementAt(i).x<<' '<<path.elementAt(i).y<<'\n';
	f.close();
}

bool subdivideMesh(Cell *qe, float maxEdge)
//splits edges whose length exceeds parameter maxEdge
{
	bool was=false;
	// first, set the splitme bits on longest edges of triangles, if the edges is longer than allowed
	CellFaceIterator cellFaces(qe);
	Face *f;
	while ((f = cellFaces.next()) != 0)
	{
	    // visit each face of cell qe
	    FaceEdgeIterator faceEdges(f);
	    Edge *edge;
	    while ((edge = faceEdges.next()) != 0)
		{
			int splitme = edge->Org() < edge->Dest() && edge->length()>maxEdge; // my Sym's bit will be the complement of mine
			// edge must be longest in either left or right triangle
			splitme=splitme && (edge->Lnext()->length()<=edge->length() && edge->Onext()->length()<=edge->length() || 
				edge->Dnext()->length()<=edge->length() && edge->Rnext()->length()<=edge->length() );
			//improve this for cases of isoscales and equilateral triangles (see if any other edge was already marked)
			edge->midpoint[2] = splitme? 1.0 : 0.0; // set bit
			if (splitme)
				was=true;
	    }
	}
	if (!was)
		return false;

	//now subdivide marked edges
	CellFaceIterator cellFaceIt(qe);
	while ((f = cellFaceIt.next()) != 0)
	{
	    // visit each face of cell qe
	    FaceEdgeIterator faceEdges(f);
	    Edge *edge;
	    while ((edge = faceEdges.next()) != 0)
		{
			// visit each edge of face f
			// if its "splitme" bit set then split it
			if (edge->midpoint[2]>0.0)
			{
				Edge *enew = split(edge);
				// clear splitme bits on two sub-edges and
				// their Syms to avoid recursive splitting
				edge->midpoint[2] = 0.0;
				edge->Sym()->midpoint[2] = 0.0;
				enew->midpoint[2] = 0.0;
				enew->Sym()->midpoint[2] = 0.0;

				qe->makeFaceEdge(enew->Left(), enew->Dest(), enew->Onext()->Dest());
				qe->makeFaceEdge(enew->Right(), enew->Dest(), enew->Rnext()->Org());
				break;
			}
	    }
	}
	return true;
}

bool pos2val(InternalImageType::Pointer image, vec3& pos, InternalImageType::PixelType &val)
//returns true if given position is within image (in that case also sets val to that voxel's value)
{
	itk::Point<float, 3> p(pos.ptr());
	typedef itk::LinearInterpolateImageFunction<InternalImageType,float> InterpolatorType;
	static InterpolatorType::Pointer interp=InterpolatorType::New();
	interp->SetInputImage(image);
	//p[0]=pos.x();
	//p[1]=pos.y();
	//p[2]=pos.z();
	InterpolatorType::ContinuousIndexType ind;
	if (image->TransformPhysicalPointToContinuousIndex<float>(p, ind)) //coordinates within image
	{			
		val=interp->EvaluateAtContinuousIndex( ind );
		return true;
	}
	else
		return false;
}

void smooth(Cell *qe, double smoothingFactor=0.2)
{
	CellVertexIterator cvi(qe);
	Vertex *v;
	while ( (v=cvi.next())!=0 )
	{
		VertexEdgeIterator vei(v);
		Edge *e;
		vec3 u(0,0,0); //position update
		unsigned count=0;
		while ( (e=vei.next())!=0 )
		{
			u+=(e->Dest()->pos-v->pos);
			count++;
		}
		v->normal=u*smoothingFactor/count; //keep new position in normal
	}
	
	CellVertexIterator cvi2(qe);
	while ( (v=cvi2.next())!=0 )
		v->pos+=v->normal; //copy smoothed position from normal
}

void MainLogic::inflate(Cell *qe, vec3& center, float avgEdgeLen)
{	
	InternalImageType::IndexType ind;
	itk::Point<float, 3> p;
	InternalImageType::PixelType pixel;
	InternalImageType::SpacingType sp = current->GetSpacing();
	vec3 pos;
	CellVertexIterator it3(qe);
	Vertex *v;
	while ((v = it3.next()) != 0)
	{
		vec3 n=v->pos-center; //center-vertex vector
		n.normalize();

		p[0]=v->pos[0]; p[1]=v->pos[1]; p[2]=v->pos[2];
		if (!current->TransformPhysicalPointToIndex(p, ind))
			//due to numerical instabilities during edge splitting, a point can end up outside
			continue; //continue while loop

		InternalImageType::PixelType lval=lImage->GetPixel(ind), hval=hImage->GetPixel(ind);
		
		if (lval>=lowThreshold && hval<=highThreshold) //inner area
		{
			pos=v->pos+n*v->stepSize;
			if (pos2val(current, pos, pixel)) //if within image
			{
				v->pos=pos;
				v->lastH=hval;
				v->lastIntensity=pixel;
			}
			else
				v->stepSize*=0.75;
		}
		else
		{
			float factor, cosPhi;
			cosPhi=(v->normal*n)/(v->normal.length()*n.length());
			//cosine of the angle between normal and center-vertex vector
			if (v->curv>0.5)
				factor=0.5/v->curv;
			else
				factor=1;

            pos=v->pos+n*factor*cosPhi*v->stepSize;

            if (avgEdgeLen>0) //if inflating for subdivision method
            {
                vec3 smd=supposedMoveDirection(v, avgEdgeLen);
                if (n*smd>0)
                pos+=smd;
            }
            
			if (pos2val(current, pos, pixel)) //if within image
			{				
				if (pixel>=lowThreshold && pixel<=highThreshold && withinEpsilonEnvironment(v->lastH,hval,0.1)
					&& (v->lastIntensity>=pixel||withinEpsilonEnvironment(v->lastIntensity, pixel, 0.05)) )
				{
					v->pos=pos;
					v->lastH=hval;
					v->lastIntensity=pixel;						
				}
				else
					v->stepSize*=0.75;
			}
			else
				v->stepSize*=0.75;
		}
	}
}

double calcAvgDistance(QPainterPath &path, QPointF center)
{
	unsigned count=0;
	double dist=0;
	for (int i=0; i<path.elementCount(); i++)
	{
		count++;
		dist+=distance(center.x(), center.y(), path.elementAt(i).x, path.elementAt(i).y);
	}
	return dist/count;
}

void saveInitPicture(const QPixmap& slice, const QPainterPath& path, QString filename)
{
    QImage pm(slice.size(), QImage::Format_RGB888);
	QPainter pnt(&pm);
    pnt.drawPixmap(0,0,slice);
    QPen pen(QBrush(Qt::yellow), 2);
	pnt.setPen(pen);
	pnt.drawPath(path);
	pnt.end();
	pm.save(filename);
}

void calcThresholds(InternalImageType::Pointer image, QPainterPath &path, double xc, double yc, int zc, InternalImageType::PixelType maxVal,
					InternalImageType::PixelType& highThreshold, InternalImageType::PixelType& lowThreshold)
{
	InternalImageType::IndexType ind;
	typedef std::map<InternalImageType::PixelType, unsigned> pCounter;
	pCounter vox;
	unsigned count=0;
    InternalImageType::PixelType noiseUnit=maxVal/1024;
	QRectF bb=path.boundingRect(); //path's bounding box
	for (int i=bb.top(); i<bb.bottom(); i++) 
		for (int k=bb.left(); k<bb.right(); k++) 
			if (path.contains(QPoint(k,i)))
			{
				ind[0]=k;
				ind[1]=i;
				ind[2]=zc;
                InternalImageType::PixelType pix=image->GetPixel(ind);
                pix=floor(pix/noiseUnit)*noiseUnit;
				vox[pix]++;
				count++;
			}
	
	for (pCounter::iterator it=vox.begin(); it!=vox.end(); )
		if (it->second<(count/1000.0)+1) //this voxel value is due to noise (less than 0.1%)
			vox.erase(it++);
		else
			it++;

	highThreshold=vox.rbegin()->first+noiseUnit; //add 1 noiseUnit to round up
	lowThreshold=vox.begin()->first; //already rounded down
}

void clipWrite(Cell *qe, InternalImageType::Pointer image, const char *filename)
{
    vtkPolyData *pd=qe2vtk(qe);

    VisualizingImageType::DirectionType d=image->GetDirection();
    vtkSmartPointer<vtkMatrix4x4> mat=vtkSmartPointer<vtkMatrix4x4>::New();
    for (int i=0; i<3; i++)
        for (int k=0; k<3; k++)
            mat->SetElement(i,k, d(i,k));
    
    vtkTransform *t=vtkTransform::New();
    t->PostMultiply();
    t->Scale(image->GetSpacing().GetDataPointer());
    t->Concatenate(mat); //rotation
    t->Translate(image->GetOrigin().GetDataPointer());
    vtkTransform *it=vtkTransform::New();
    it->SetMatrix(t->GetMatrix());
    it->Inverse();

    vtkTransformPolyDataFilter *itp=vtkTransformPolyDataFilter::New();
    itp->SetTransform(it);
    itp->SetInputData(pd);

    int xs=image->GetLargestPossibleRegion().GetSize(0);
    int ys=image->GetLargestPossibleRegion().GetSize(1);
    int zs=image->GetLargestPossibleRegion().GetSize(2);

    vtkPlaneCollection *pc=vtkPlaneCollection::New();
    {
        vtkPlane *x0=vtkPlane::New();
        x0->SetOrigin(-0.5,0,0);
        x0->SetNormal(1,0,0);
        pc->AddItem(x0);

        vtkPlane *y0=vtkPlane::New();
        y0->SetOrigin(0,-0.5,0);
        y0->SetNormal(0,1,0);
        pc->AddItem(y0);

        vtkPlane *z0=vtkPlane::New();
        z0->SetOrigin(0,0,-0.5);
        z0->SetNormal(0,0,1);
        pc->AddItem(z0);

        vtkPlane *x1=vtkPlane::New();
        x1->SetOrigin(xs-0.5,0,0);
        x1->SetNormal(-1,0,0);
        pc->AddItem(x1);

        vtkPlane *y1=vtkPlane::New();
        y1->SetOrigin(0,ys-0.5,0);
        y1->SetNormal(0,-1,0);
        pc->AddItem(y1);

        vtkPlane *z1=vtkPlane::New();
        z1->SetOrigin(0,0,zs-0.5);
        z1->SetNormal(0,0,-1);
        pc->AddItem(z1);
    }

    vtkClipClosedSurface *ccs=vtkClipClosedSurface::New();
    ccs->SetInputConnection(itp->GetOutputPort());
    ccs->SetClippingPlanes(pc);
    ccs->Update();

    vtkCleanPolyData *cpd=vtkCleanPolyData::New();
    cpd->SetInputConnection(ccs->GetOutputPort());
    cpd->Update();

    vtkTransformPolyDataFilter *tp=vtkTransformPolyDataFilter::New();
    tp->SetTransform(t);
    tp->SetInputConnection(cpd->GetOutputPort());
    tp->Update();

    vtkOBJWriter *w=vtkOBJWriter::New();
    w->SetFileName(filename);
    w->SetInputConnection(tp->GetOutputPort());
    w->Update();
}

float mergability(Face *f1, Face *f2, float maxAngle) //returns mergability for 2 adjecant faces
{
	double cosFi=f1->normal*f2->normal/(f1->normal.length()*f2->normal.length()); //cosine of the angle
	if (acos(cosFi)>maxAngle) //fi>maxAngle
		return 0;
	return cosFi*(f1->normal.length()+f2->normal.length())
		/(1+abs(f1->normal.length()-f2->normal.length())); //prefer merging same sized faces
}

void removeDanglingEdges(Cell *s)
{
	CellVertexIterator cvi(s);
	Vertex *v;
	bool work_exists;
	do
	{
		work_exists=false;
		while ((v = cvi.next()) != 0)
		{
			Edge *e=v->getEdge();
			unsigned count=0;
			while (e==e->Onext()) //this vertex is on the end of a dangling edge
			{
				Edge *t=e;
				e=e->Lnext();
				v=t->Org();
				if (t->Left()->getEdge()==t || t->Left()->getEdge()==t->Sym())
					t->Left()->addEdge(e);
				e->Org()->addEdge(e);
				Edge::kill(t);
				s->removeVertex(v);
				count++;
				work_exists=true;
			}
			if (count>1) //we need to reinitialize the iterator, or we might fetch a deleted vertex
				cvi=CellVertexIterator(s);
		}
	} while (work_exists);

	//make copies of duplicate vertices in face vertex list
	CellFaceIterator cellFaceIt(s);
	Face *f;
	while ((f = cellFaceIt.next()) != 0)
	{
	    // visit each face
	    FaceEdgeIterator faceEdges(f);
		std::vector<Edge*> edges;
		Edge *edge;
	    while ((edge = faceEdges.next()) != 0)
			edges.push_back(edge);

		for (int i=0; i<edges.size()-1; i++)
			for (int k=i+1; k<edges.size(); k++)
				if (edges[i]->Org()==edges[k]->Org())
				{
					edge=s->makeVertexEdge(edges[i]->Org(), edges[i]->Right(), edges[k]->Right());
					//maybe nudge position a bit?
				}
	}
}

bool mergeStep(Cell *s, float maxAngle) //merge faces by removing edges, returns false when merges are not possible
{
	bool merged=false;
	CellFaceIterator cellFaceIt(s);
	Face *f;
	while ((f = cellFaceIt.next()) != 0)
	{
	    // visit each face
	    FaceEdgeIterator faceEdges(f);
	    Edge *edge;
		float maxMerge=0;
		Edge *eMerge;
	    
		while ((edge = faceEdges.next()) != 0)
		{
			Face *ftm=edge->Right();
			if (ftm==f) //can happen due to a redundant edge
				continue;
			float m=mergability(f, ftm, maxAngle);
			if (m>maxMerge)
			{
				maxMerge=m;
				eMerge=edge;
			}
	    }

		if (maxMerge>0)
		{
			Face *right=eMerge->Right();
			//preserve face size and orientation through normals
			f->normal+=right->normal;
			s->killFaceEdge(eMerge);
			merged=true;
		}
	}	
	return merged;
}

inline bool vecLength(const vec3& a, const vec3& b)
{
	return a.length2()>b.length2();
}

bool Axis(vec3& axis, const std::vector<vec3>& normalsOfMergedFaces)
//first 2 normals are averaged to extract axis
{
	vec3 a=normalsOfMergedFaces.at(0), b=normalsOfMergedFaces.at(1);
	bool retVal=true;
	if (a*b>=0) //top and bottom pointing normals are not nearly opposite
		retVal=false;
	else
		b=-b;
	axis=a+b;
	axis.normalize();
	return retVal;
}
std::vector<vec3> adjecantVertebrae(Cell *qe, float avgDist, vec3& center, 
										 InternalImageType::Pointer image, vec3 lastAdjecant=vec3(0,0,0))
{
	std::vector<vec3> f, fNew; //list of faces to be returned
    Cell *s=qe->deepCopy();
	
	calcFaceNormals(s); //calculate initial face normals

	//0.1rad = 5.8degrees; 0.78rad = 45deg
	static int iter=0;
	//objWriteCell(s, ("D:\\Temp\\merge"+QString::number(++iter)+"_before.obj").toStdString().c_str());
	for (float angle=0.1; angle<=0.35; angle+=0.1) //due to numerical instabilities the angle ends up being 4.000001
	{
		while (mergeStep(s, angle)) //do mergesteps while merging is possible
		{
			//objWriteCell(s, ("D:\\Temp\\step"+QString::number(++iter).rightJustified(2,'0')+".obj").toStdString().c_str());
		}
	//removeDanglingEdges(s); removeDanglingEdges(s); removeDanglingEdges(s); //once is sometimes not enough
	//objWriteCell(s, ("D:\\Temp\\merge"+QString::number(iter)+"angle"+QString::number(angle)+".obj").toStdString().c_str());
	}
	//objWriteCell(s, ("D:\\Temp\\merge"+QString::number(iter)+"after.obj").toStdString().c_str());

	CellFaceIterator cellFaceIt(s);
	Face *face;
	while ((face = cellFaceIt.next()) != 0)
		f.push_back(face->normal); //put all resulting normals into an array

	//sort(f.begin(), f.end(), vecLength); //largest vec not interesting any more, using parameter lastAdjecant
	if (lastAdjecant!=vec3(0,0,0))
	{
		vec3 p;
		fNew.reserve(f.size()/2);
		for (std::vector<vec3>::iterator it=f.begin(); it!=f.end(); ++it)
		{
			float cosFi=(*it)*lastAdjecant/(it->length()*lastAdjecant.length());
			if (cosFi>0.2) //in the general wanted direction
			//cosFi>0.2 instead of cosFi>0 excludes directions normal to wanted direction
			{
				//InternalImageType::PixelType pix;
				//p=*it;
				//p.normalize();
				//p*=avgDist*2;
				//p+=center; //candidate position of adjecant vertebra
				//if (pos2val(image, p, pix)) //coordinates within image
				//above makes orientation of the last vertebra incorrect
				{
					if (cosFi>0.7) //direction similar to last adjecant vector
						(*it)*=10*(cosFi-0.6); //increase importance of that normal up to 4x
					fNew.push_back(*it);
				}
			}
		}
		f.swap(fNew);
		fNew.clear();
	}

	sort(f.begin(), f.end(), vecLength);
	f.resize(4);
	return f;
}

vtkActor *poly2actor(vtkPolyData *poly, bool wireframe)
{
    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputData(poly);
    vtkActor *actor=vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1,1,0);
    if (wireframe)
        actor->GetProperty()->SetRepresentationToWireframe();
    return actor;
}

vtkActor* symVert2vtk(Vertebra *vertebra)
{
    vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();
    arrowSource->Update(); //pointing in 1,0,0 direction

    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(arrowSource->GetOutputPort());
    vtkActor *actor =vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.8, 0.0, 0.2);

    vec3 arrow(1,0,0), rotAxis=arrow^vertebra->axis; //cross product is rotation axis
    double angle=acos(arrow*vertebra->axis/(arrow.length()*vertebra->axis.length())); //angle from dot product

    actor->RotateWXYZ(angle*180/M_PI, rotAxis[0], rotAxis[1], rotAxis[2]);
    actor->SetScale(vertebra->axis.length()*1.2);
    actor->SetPosition(vertebra->center[0], vertebra->center[1], vertebra->center[2]);
    return actor;
}

double MainLogic::segmentOne(Vertebra *vertebra)
//returns average distance from center to surface
{
	InternalImageType::SpacingType sp = current->GetSpacing();
	double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);
	double minSpacing=min(min(sp[0], sp[1]), sp[2]);
	string tricube="# Wavefront OBJ format, triangulated cube: 8 vertices, 18 edges, 12 faces\n\
\n\
# vertices\n\
v -1.0 -1.0 -1.0\n\
v -1.0 -1.0  1.0\n\
v -1.0  1.0 -1.0\n\
v -1.0  1.0  1.0\n\
v  1.0 -1.0 -1.0\n\
v  1.0 -1.0  1.0\n\
v  1.0  1.0 -1.0\n\
v  1.0  1.0  1.0\n\
\n\
# faces\n\
f 8 2 6\n\
f 2 8 4\n\
f 5 3 7\n\
f 7 6 5\n\
f 8 3 4\n\
f 3 8 7\n\
f 4 1 2\n\
f 1 4 3\n\
f 3 5 1\n\
f 6 7 8\n\
f 5 2 1\n\
f 2 5 6\n\
";
	vertebra->qe=objReadCellFromString(tricube);
	CellVertexIterator it3(vertebra->qe);
	Vertex *v;
	while ((v = it3.next()) != 0)
	{
		v->stepSize=minSpacing;
		v->lastH=highThreshold;
		v->lastIntensity=(highThreshold+lowThreshold)/2;
	}
    translate(vertebra->qe, vertebra->center);
    
    if (mainForm.actionInteractive_inflation->isChecked())
    {
        vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
    }
	updateNormalsAndCurvatures(vertebra->qe);

    double oldDist, dist2d=vertebra->radius;
    vertebra->radius=1.4;
	int iter=0;
	do
	{
		inflate(vertebra->qe, vertebra->center, 0);
		if (iter>15) //required when maxEdge given to subdivideMesh=2.95*avgSpacing (not needed for 1.95)
			smooth(vertebra->qe, 0.2);
		subdivideMesh(vertebra->qe, 2.95*avgSpacing);
		updateNormalsAndCurvatures(vertebra->qe);

		if (mainForm.actionInteractive_inflation->isChecked())
		{
			mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra->actor);
            vertebra->actor->Delete();
            vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
			mainForm.vis->GetRenderWindow()->Render();
            QApplication::processEvents();
		}

        vertebra->calcNewCenter();
		oldDist=vertebra->radius;
        vertebra->calcAvgDistance(); //update radius
		iter++;
	} while (iter<200 //max iteration count
		&& vertebra->radius<dist2d*1.5 //up to 50% greater than user's approximate selection
		&& (vertebra->radius-oldDist)>0.02*avgSpacing); //and the growth rate is not too slow
	
	updateNormalsAndCurvatures(vertebra->qe);
    if (mainForm.actionInteractive_inflation->isChecked())
    {
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra->actor);
        vertebra->actor->Delete();
    }
    vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
    mainForm.vis->GetRenderWindow()->Render();
    QApplication::processEvents();

	return vertebra->radius;
}

template<typename _MatrixType> class myLU: public Eigen::PartialPivLU<_MatrixType>
{
public:
    myLU(): Eigen::PartialPivLU<_MatrixType>() {}
    myLU(Index size): Eigen::PartialPivLU<_MatrixType>(size) {}
    myLU(const MatrixType& matrix): Eigen::PartialPivLU<_MatrixType>(matrix) {}

    void markInitialized()
    {
        m_isInitialized=true;
    }
};

void updateWeightMatrix(Eigen::MatrixXd& weights, myLU<Eigen::MatrixXd>& LU,
	Cell *qe, unsigned maxLevel, unsigned freeLevels, unsigned& freeCount)
{
	QString levelInfo=QString::number(maxLevel)+"_"+QString::number(freeLevels);
    QDir cache("./cache_LU");
    if (!cache.exists())
        cache.mkpath(".");
	QFile fw("./cache_LU/W_"+levelInfo);
	QFile flu("./cache_LU/LU_"+levelInfo);
    QFile fp("./cache_LU/P_"+levelInfo);
	if (fw.exists()&&flu.exists()&&fp.exists())
	{
		Vertex *v;
		freeCount=0;
		CellVertexIterator cvi0(qe);
		while ((v = cvi0.next())!=0)
			if (v->level<=freeLevels)
				freeCount++;
		//we now have count of free vertices
		weights.resize(freeCount, qe->countVertices());
		fw.open(QIODevice::ReadOnly);
		fw.read((char *)weights.data(), freeCount*qe->countVertices()*sizeof(double));
		fw.close();
        LU=myLU<Eigen::MatrixXd>(freeCount);
#ifdef _DEBUG
        LU.markInitialized(); //do away with failing assertion
#endif
		flu.open(QIODevice::ReadOnly);
		flu.read((char *)LU.matrixLU().data(), freeCount*freeCount*sizeof(double));
		flu.close();
		fp.open(QIODevice::ReadOnly);
		fp.read((char *)LU.permutationP().indices().data(), freeCount*sizeof(int));
		fp.close();
	}
	else
	{
		Eigen::MatrixXd t=createWeightMatrix(qe, maxLevel, freeLevels, freeCount);
		weights=t;
		t*=t.transpose();
        LU.compute(t);
		fw.open(QIODevice::WriteOnly);
		fw.write((char *)weights.data(), freeCount*qe->countVertices()*sizeof(double));
		fw.close();
		flu.open(QIODevice::WriteOnly);
		flu.write((char *)LU.matrixLU().data(), freeCount*freeCount*sizeof(double));
		flu.close();
		fp.open(QIODevice::WriteOnly);
		fp.write((char *)LU.permutationP().indices().data(), freeCount*sizeof(int));
		fp.close();
	}	
}

double MainLogic::segmentOneButterfly(Vertebra *vertebra)
//returns average distance from center to surface
{
	InternalImageType::SpacingType sp = current->GetSpacing();
	double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);
	double minSpacing=min(min(sp[0], sp[1]), sp[2]);
    double dist2d=vertebra->radius;
	vertebra->qe=objReadCell("subdivBody.obj");
	//vertebra->qe=objReadCell("icosahedron.obj");
	calcVertexValence(vertebra->qe, false);
    vertebra->calcAvgDistance();
	CellVertexIterator it3(vertebra->qe);
	Vertex *v;
	while ((v = it3.next()) != 0)
	{
		v->stepSize=minSpacing;
		v->lastH=highThreshold;
		v->lastIntensity=(highThreshold+lowThreshold)/2;
	}

	rotate(vertebra->qe, vec3(0,0,1), vertebra->axis);
    translate(vertebra->qe, vertebra->center);
	copyPos2SubdivPos(vertebra->qe);
    if (mainForm.actionInteractive_inflation->isChecked())
    {
        vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
    }

	unsigned maxLevel=1, freeLevels=0, freeCount;
	for (unsigned i=1; i<=maxLevel; i++)
		topologicallySubdivide(vertebra->qe, i);

	Eigen::Matrix<double,Eigen::Dynamic,3> P;
	P.resize((vertebra->qe)->countVertices(), 3);
	Eigen::MatrixXd W;
    myLU<Eigen::MatrixXd> LU;
	//QTime timeInv=QTime::currentTime();
	updateWeightMatrix(W, LU, vertebra->qe, maxLevel, freeLevels, freeCount);

	//QTime time0=QTime::currentTime();
	calculateSubdivisionPositions(vertebra->qe, maxLevel,freeLevels);
	//QTime time1=QTime::currentTime();
	//{
	//	CellVertexIterator cvi(vertebra->qe);
	//	while ((v = cvi.next())!=0)
	//		if (v->level>freeLevels) //freelevels' positions are not calculated
	//		{
	//			v->pos=vec3(0,0,0);
	//			for (unsigned i=0; i<freeCount; i++)
	//			{
	//				v->pos[0]+=P(i,0)*W.coeff(i,v->realIndex());
	//				v->pos[1]+=P(i,1)*W.coeff(i,v->realIndex());
	//				v->pos[2]+=P(i,2)*W.coeff(i,v->realIndex());
	//			}
	//		}
	//}
	//QTime time2=QTime::currentTime();
	//QMessageBox::warning(0, "Operation duration (ms)", QString("Standard subdivision: ")+QString::number(time0.msecsTo(time1))
	//										+"\nMatrix multiplication subdivision: "+QString::number(time1.msecsTo(time2))
	//										+"\nWeight creation and matrix inversion: "+QString::number(timeInv.msecsTo(time0)));
	swapSubdivPos(vertebra->qe); copyPos2SubdivPos(vertebra->qe); //together equivalent to copySubdivPos2Pos
	updateNormalsAndCurvatures(vertebra->qe);

	double oldDist;
    vec3 origCenter=vertebra->center;
	float edgeLen, stdDev, minEdge, maxEdge;
    edgeLen=getAverageEdgeLength(vertebra->qe, minEdge, maxEdge, stdDev);
	int iter=0, extraIter=0;
	do
	{
        if (mainForm.actionMaintain_vertex_surface_density->isChecked())
		    inflate(vertebra->qe, vertebra->center, edgeLen);
        else
            inflate(vertebra->qe, vertebra->center, 0);
		if (mainForm.actionSubdivisionOptimal->isChecked())
		{
			if (iter>5)
				smooth(vertebra->qe, 0.5);
			CellVertexIterator cvi0(vertebra->qe);
			while ((v = cvi0.next())!=0)
			{
				P(v->realIndex(),0)=v->pos[0];
				P(v->realIndex(),1)=v->pos[1];
				P(v->realIndex(),2)=v->pos[2];
			}
		}

        //objWriteCell(vertebra->qe, (QString::fromStdString(filename_base)+QString::number(vertebraNumber)
        //             +"_"+QString::number(iter)+"before.obj").toStdString().c_str());
		oldDist=vertebra->radius;
		vertebra->calcAvgDistance(); //update radius
		vertebra->calcNewCenter();
		edgeLen=getAverageEdgeLength(vertebra->qe, minEdge, maxEdge, stdDev);
        if (mainForm.actionMaintain_vertex_surface_density->isChecked())
            rearrangeVertices(vertebra->qe, edgeLen,vertebra->center );

        if (freeLevels<maxLevel)
		{
			if (mainForm.actionSubdivisionHeuristic->isChecked())
				normalizeSubdivisionHierarchy(vertebra->qe, maxLevel, freeLevels);
			else
			{
				Eigen::Matrix<double,Eigen::Dynamic,3> N;
                N=LU.solve(W*P); //N=WWTinv*(W*P);
				CellVertexIterator cvi(vertebra->qe);
				while ((v = cvi.next())!=0)
					if (v->level<=freeLevels)
					{
						v->pos[0]=N(v->realIndex(),0);
						v->pos[1]=N(v->realIndex(),1);
						v->pos[2]=N(v->realIndex(),2);
					}
				calculateSubdivisionPositions(vertebra->qe, maxLevel, freeLevels);
			}
		}
		//objWriteCell(vertebra->qe, (QString::fromStdString(filename_base)+QString::number(vertebraNumber)
  //                   +"_"+QString::number(iter)+"optimal.obj").toStdString().c_str());

		if (edgeLen>avgSpacing*1.95)
		{
            if (stdDev>0.8 && (vertebra->center-origCenter).length()>1) //large difference between shortest and longest edge
            //and center has moved at least 1mm (this condition prevents infinite loops)
            {
                if (mainForm.actionInteractive_inflation->isChecked())
                {
                    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra->actor);
                    vertebra->actor->Delete();
                }
                //objWriteCell(vertebra->qe, (filename_base+"_"+Vertebra::labels[vertebra->labelIndex]+"_wrong.obj").c_str());
                Cell::kill(vertebra->qe);
                QApplication::processEvents();
                vertebra->radius=dist2d;
                return segmentOneButterfly(vertebra); //start with larger size of initial shape to avoid wasting time?
                //restart segmentation of this vertebra with much better center estimate
            }

			topologicallySubdivide(vertebra->qe, ++maxLevel);
			calculateSubdivisionPositions(vertebra->qe, maxLevel, freeLevels);
			if (maxLevel-2>freeLevels)
				freeLevels++;
			if (mainForm.actionSubdivisionOptimal->isChecked() && freeLevels<maxLevel)
			{
				P.resize((vertebra->qe)->countVertices(), 3);
				updateWeightMatrix(W, LU, vertebra->qe, maxLevel, freeLevels, freeCount);
			}
		}
		
		updateNormalsAndCurvatures(vertebra->qe);
		if (mainForm.actionInteractive_inflation->isChecked())
		{
			mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra->actor);
            vertebra->actor->Delete();
            vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
			mainForm.vis->GetRenderWindow()->Render();
            QApplication::processEvents();
		}

		if (extraIter>0)
			extraIter--;
		if (extraIter==0 && vertebra->radius-oldDist<0.2*minSpacing)
			if (++freeLevels<=maxLevel)
			{
				extraIter=exp2i(maxLevel-freeLevels-1);
				if (mainForm.actionSubdivisionOptimal->isChecked() && freeLevels<maxLevel)
				{
					P.resize((vertebra->qe)->countVertices(), 3);
					updateWeightMatrix(W, LU, vertebra->qe, maxLevel, freeLevels, freeCount);
				}
			}
		iter++;
	} while (extraIter>0 || iter<200 //if extra iterations given OR maxIter not reached and
		&& vertebra->radius<dist2d*1.5 //up to 50% greater than user's approximate selection
		&& vertebra->radius-oldDist>0.2*minSpacing); //and the growth rate is not too slow
	
	updateNormalsAndCurvatures(vertebra->qe);

    if (mainForm.actionInteractive_inflation->isChecked())
    {
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra->actor);
        vertebra->actor->Delete();
    }
    vertebra->actor=poly2actor(qe2vtk(vertebra->qe), mainForm.actionShow_wireframe->isChecked());
    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra->actor);
    mainForm.vis->GetRenderWindow()->Render();
    QApplication::processEvents();
	return vertebra->radius;
}

void MainLogic::segmentationLoop(const vec3& last2Center, const Vertebra& last,
                                 const vec3& lastAdjecantVertebraVector, double sizeGrowth, int labelIndexDirection)
{
    Vertebra *vertebra1=new Vertebra(*this);
    if (labelIndexDirection==1)
        vertebra.push_back(vertebra1);
    else
        vertebra.insert(vertebra.begin(),vertebra1);
    vertebra1->labelIndex=last.labelIndex+labelIndexDirection;
    vertebra1->radius=last.radius;

    vec3 ccVec, ccn;
    std::vector<vec3> adjecantVertebraeList;

	ccVec=last.center-last2Center;

	ccn=ccVec;
	ccn.normalize();
	vertebra1->center=last.center+(lastAdjecantVertebraVector*0.5+ccn*0.5)*ccVec.length()*(1.0+sizeGrowth);
    vertebra1->axis=labelIndexDirection*ccVec/2;

	mainForm.statusbar->showMessage("Segmenting vertebra "+QString::fromStdString(Vertebra::labels[vertebra1->labelIndex]));
	time=QTime::currentTime();
    
	if (!mainForm.actionSubdivisionDisabled->isChecked())
        segmentOneButterfly(vertebra1);
	else
		segmentOne(vertebra1);

    adjecantVertebraeList=adjecantVertebrae(vertebra1->qe, vertebra1->radius, vertebra1->center, current, lastAdjecantVertebraVector);

    f<<time.toString("hh:mm:ss.zzz\t").toStdString()<<time.msecsTo(QTime::currentTime())
		<<"\tSegmenting vertebra "<<Vertebra::labels[vertebra1->labelIndex]<<"\n";	

    if (abs(vertebra1->radius-last.radius)<0.3*vertebra1->radius && vertebra1->volume/last.volume>0.5 && vertebra1->volume/last.volume<2 &&
        abs(vertebra1->radius-vertebra0->radius)<0.8*vertebra0->radius && vertebra1->volume/vertebra0->volume>0.08 && vertebra1->volume/vertebra0->volume<12.5 &&
        (vertebra1->center-last.center).length()>vertebra1->radius /*if center guess is too close */)
    {
        clipWrite(vertebra1->qe, current, (filename_base+"_"+Vertebra::labels[vertebra1->labelIndex]+".obj").c_str());
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(symVert2vtk(vertebra1));
        adjecantVertebraeList[0].normalize();
        ccn=vertebra1->center-last.center;
        sizeGrowth=(ccn.length()-ccVec.length())/ccn.length(); //expect the same growth in this direction
        segmentationLoop(last.center, *vertebra1, adjecantVertebraeList[0], sizeGrowth, labelIndexDirection); //recursion
    }
    else
    {
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra1->actor);
        vertebra1->actor->Delete();
        if (labelIndexDirection==1)
            vertebra.pop_back();
        else
            vertebra.erase(vertebra.begin());
    }
}

void MainLogic::segmentVertebrae(vec3& originalCenter, double dist2D, int labelIndex)
//segment first and second vertebra, then call segmentation loop
{
	InternalImageType::SpacingType sp = current->GetSpacing();
	vec3 lastAdjecantVertebraVector, firstAdjecantVertebraVector;
    vertebra.clear();
    vertebra0=new Vertebra(*this);
    vertebra.push_back(vertebra0);
    vertebra0->center=originalCenter;
    vertebra0->axis=vec3(0,0,1);
    vertebra0->radius=dist2D;
    vertebra0->labelIndex=labelIndex;

	mainForm.statusbar->showMessage("Segmenting initial vertebra...");
	QTime time=QTime::currentTime();
	if (!mainForm.actionSubdivisionDisabled->isChecked())
        segmentOneButterfly(vertebra0);
		//first vertebra, axis is just a guess - and it does not influence result that much
	else
        segmentOne(vertebra0);

	f.open((filename_base+"_times.csv").c_str());
	f<<"StartTime	Duration|ms	Description\n";
	f<<time.toString("hh:mm:ss.zzz\t").toStdString()
		<<time.msecsTo(QTime::currentTime())<<"\tInitial vertebra\n";

    clipWrite(vertebra0->qe, current, (filename_base+"_"+Vertebra::labels[vertebra0->labelIndex]+".obj").c_str());
	
    std::vector<vec3> adjecantVertebraeList=adjecantVertebrae(vertebra0->qe, vertebra0->radius, originalCenter, current);

    Vertebra *vertebra1=new Vertebra(*this);
    vertebra.push_back(vertebra1);

	//segment first adjacent vertebra
	mainForm.statusbar->showMessage("Segmenting first adjacent vertebra...");
	std::vector<vec3>::iterator it=adjecantVertebraeList.begin();
	const int numberOfFactors=4;
	double distanceFactors[numberOfFactors]={2.0,1.8,2.2,1.6};
	bool found=false;
    int labelIndexDirection=0;
	while (it!=adjecantVertebraeList.end())
	{
		for (int i=0; i<numberOfFactors; i++)
		{
			vec3 p=*it;
			p.normalize();
			vertebra1->axis=p;
            p*=vertebra0->radius*distanceFactors[i];
			p+=originalCenter; //candidate position of adjecant vertebra
            vertebra1->center=p;
            vertebra1->radius=vertebra0->radius;
            vertebra1->volume=vertebra0->volume;
			time=QTime::currentTime();
            
            if ((*it)*vec3(0,0,1)>0)//pointing up
            {
                vertebra1->labelIndex=vertebra0->labelIndex+1;
                labelIndexDirection=+1;
            }
            else
            {
                vertebra1->labelIndex=vertebra0->labelIndex-1;
                labelIndexDirection=-1;
            }            

			if (!mainForm.actionSubdivisionDisabled->isChecked())
				segmentOneButterfly(vertebra1);
			else
				segmentOne(vertebra1);
			f<<time.toString("hh:mm:ss.zzz\t").toStdString()
				<<time.msecsTo(QTime::currentTime())<<"\tFirst attempt to locate adjacent vertebra\n";	
			if (abs(vertebra1->radius-vertebra0->radius)<0.2*vertebra0->radius && vertebra1->volume/vertebra0->volume>0.5 && vertebra1->volume/vertebra0->volume<2)
			//within 20% of the first one by distance (cubed for volume)
			{
				firstAdjecantVertebraVector=(p-originalCenter);
				lastAdjecantVertebraVector=firstAdjecantVertebraVector;
				found=true;
                if (firstAdjecantVertebraVector*vec3(0,0,1)<0)//pointing down
                {
                    //we are going down the spine so places should be reversed
                    vertebra[0]=vertebra1;
                    vertebra[1]=vertebra0;
                }
				break;
			}
            else //unsuccessful, delete
            {
                mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra1->actor);
                vertebra1->actor->Delete();
                vertebra1->actor=0;
            }
		}
		if (found)
			break;
		else
			it++; //next vector
	}
	if (!found)
	{
		if (vertebra1->actor)
        {
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra1->actor);
            vertebra1->actor->Delete();
        }
        vertebra.pop_back(); //we only managed to segment the first vertebra		
	}
    else
    {
	    if (it!=adjecantVertebraeList.begin()) //first guess is not successful
	    {
		    vertebra0->axis=*it;
		    vertebra0->axis.normalize();
	    }
	    else
		    Axis(vertebra0->axis, adjecantVertebraeList);
        vertebra0->axis*=labelIndexDirection*vertebra0->radius;
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(symVert2vtk(vertebra0));

        clipWrite(vertebra1->qe, current, (filename_base+"_"+Vertebra::labels[vertebra1->labelIndex]+".obj").c_str());
        adjecantVertebraeList=adjecantVertebrae(vertebra1->qe, vertebra1->radius, vertebra1->center, current, lastAdjecantVertebraVector);

	    Axis(vertebra1->axis, adjecantVertebraeList);
	    vertebra1->axis*=labelIndexDirection*vertebra1->radius;
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(symVert2vtk(vertebra1));

        adjecantVertebraeList[0].normalize();
        segmentationLoop(vertebra0->center, *vertebra1, adjecantVertebraeList[0], 0/*initial guess*/, labelIndexDirection);

        double sizeGrowth=pow((double)(vertebra[vertebra.size()-1]->center-vertebra[vertebra.size()-2]->center).length()
                                     /(vertebra[1]->center-vertebra[0]->center).length(), 1.0/(vertebra.size()-1))-1;
        //find step of geometric progression
    
	    //now go in the opposite direction (from the first vertebra)
        firstAdjecantVertebraVector.normalize();
        segmentationLoop(vertebra1->center, *vertebra0, -firstAdjecantVertebraVector, labelIndexDirection*sizeGrowth, -labelIndexDirection);
    }

    postSegmentation();
    f.close(); //save execution times
}

void MainLogic::postSegmentation()
{
    QApplication::processEvents();
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Writing segmentation masks"<<endl;
    for (int i=0; i<vertebra.size(); i++)
        if (mainForm.actionSave_binary_masks->isChecked())
        {
            mainForm.statusbar->showMessage(QString("Saving ")+Vertebra::labels[vertebra[i]->labelIndex]+" binary mask");
            //if (!mainForm.actionShow_mask_overlays->isChecked())
                vertebra[i]->getMask();
            writeImage(vertebra[i]->mask, filename_base+"_"+Vertebra::labels[vertebra[i]->labelIndex]+".mha", true);
        }

    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Finished writing masks"<<endl;

    if (vertebra.size()>1) //=1 => error or dummy test
    {
        //fit positions and orientations to polyline. initial guess: a0=x0, b0=y0; other ai,bi =0
        Eigen::VectorXd x,y;
        polynomialFit(4, vertebra,x,y, true);

        //show polyline for debugging purposes
        double lowerZ=vertebra[0]->center[2];
        double upperZ=vertebra[vertebra.size()-1]->center[2];
        vtkActor *polyLine=centerline(x,y,lowerZ,upperZ, vec3(0,0,1)); //blue
        diagnose(vertebra,x,y, filename_base);

        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(polyLine);
        mainForm.vis->GetRenderWindow()->Render();
        QApplication::processEvents();
    
        f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Finished diagnosis"<<endl;
    }
    return;
}

unsigned int nextPowerOf2(unsigned int number)
//return next power of two which is greater or equal to the given number
{
	double r=log((double)number)/log(2.0);
	return exp2i(ceil(r));
}

void MainLogic::majorSliceUpdate(int sliceNumber, int labelIndex)
{
    mainForm.painter->path.closeSubpath();
	mainForm.painter->path.setFillRule(Qt::WindingFill);

	mainForm.statusbar->showMessage("Calculating center of region...");
	double xc=0, yc=0;
	unsigned long long count=0;
	QRectF lh(mainForm.painter->path.boundingRect());
	for (int i=lh.y(); i<lh.y()+lh.height(); i++)
		for (int k=lh.x(); k<lh.x()+lh.width() ; k++)
			if (mainForm.painter->path.contains(QPointF(k,i)))
			{
				xc+=k;
				yc+=i;
				count++;
			}
	if (count>0)
	{
		xc/=count;
		yc/=count;
	}

	mainForm.statusbar->showMessage("Calculating threshold...");
	calcThresholds(current, mainForm.painter->path, xc, yc, sliceNumber, maxValue2, highThreshold, lowThreshold);


	InternalImageType::PointType p;
	InternalImageType::IndexType ind;
	ind[0]=xc; ind[1]=yc; ind[2]=sliceNumber;
	current->TransformIndexToPhysicalPoint(ind, p);
	vec3 center(p[0], p[1], p[2]);

	filename_base=tempDir+itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);

    ofstream f((QString::fromStdString(filename_base)+".init").toStdString().c_str());
    f<<"file\n";
    f<<volume_filename<<'\n';
    f<<sliceNumber<<'\n';
    f<<Vertebra::labels[labelIndex]<<'\n';
    f<<"vertebra outline:\n";
	savePath(mainForm.painter->path, f);

	mainForm.statusbar->showMessage("Segmenting vertebrae...");
	InternalImageType::SpacingType sp = current->GetSpacing();
	double dist2d=calcAvgDistance(mainForm.painter->path, QPointF(xc, yc))*pow(sp[0]*sp[1], 0.5); //geometric mean

	segmentVertebrae(center, dist2d, labelIndex);
	
    saveInitPicture(mainForm.painter->slices[sliceNumber], mainForm.painter->path, QString::fromStdString(filename_base)+"_"+QString::number(sliceNumber)+".png");

    mainForm.vis->GetRenderWindow()->Render();

	mainForm.statusbar->showMessage("Ready");
}