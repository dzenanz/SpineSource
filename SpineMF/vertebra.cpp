#include "vertebra.h"
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkArrowSource.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkAddImageFilter.h"
#include "vtkCellLocator.h"
#include "vtkTransform.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkTransformFilter.h"
#include "MainLogic.h"
#include <algorithm>
#include <vector>

Cell* Vertebra::model=0;
double Vertebra::modelRadius;
vtkPolyData* Vertebra::polyModel;
//const char *vertebraMeshFilename="vertebraModel.obj";
const char *vertebraMeshFilename="subdivBody.obj";
const char* Vertebra::labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

vtkPolyData* qe2vtk(Cell *qe)
{
    vtkPolyData *poly=vtkPolyData::New();
    vtkSmartPointer<vtkPoints> pt = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> p = vtkSmartPointer<vtkCellArray>::New();
    
    pt->Allocate(qe->countVertices());
    p->Allocate(qe->countFaces());

	CellFaceIterator cellFaces(qe);
	Face *f;
	while ((f = cellFaces.next()) != 0)
	{			
		Edge *e=f->getEdge();

		Vertex *v=e->Org();
        vtkIdType ids[3];
        ids[0]=pt->InsertNextPoint(v->pos.ptr());

		e=e->Lnext();
		v=e->Org();
		ids[1]=pt->InsertNextPoint(v->pos.ptr());

		e=e->Lnext();
		v=e->Org();
		ids[2]=pt->InsertNextPoint(v->pos.ptr());

        p->InsertNextCell(3, ids);
	}

    poly->SetPoints(pt);
    poly->SetPolys(p);
    return poly;
}

void dumpPoints(vtkPoints *points, char *filename)
{
    std::ofstream f(filename);
    double p[3];
    for (int i=0; i<points->GetNumberOfPoints(); i++)
    {
        points->GetPoint(i,p);
        f<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
    }
    f.close();
}

//align model with current mesh. vertices correspond 100%. model scaling included
vtkMatrix4x4* Vertebra::calcRigidTransform(unsigned upToLevel)
{
    vtkSmartPointer<vtkPoints> cP=vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> mP=vtkSmartPointer<vtkPoints>::New();
    CellVertexIterator cviC(qe), cviM(model);
    Vertex *v;

    unsigned i=0;
    cP->SetNumberOfPoints(30722); //5th level of subdivision
	while ((v = cviC.next()) != 0)
        if (v->level<=upToLevel)
	    {
            cP->SetPoint(i,v->pos.ptr());
            i++;
	    }
    cP->SetNumberOfPoints(i);

    mP->SetNumberOfPoints(i);
    i=0;
	while ((v = cviM.next()) != 0)
        if (v->level<=upToLevel)
	    {
            mP->SetPoint(i,v->pos.ptr());
            i++;
	    }
    
    //dumpPoints(cP,"D:/Temp/pTarget.txt");
    //dumpPoints(mP,"D:/Temp/pSource.txt");

    vtkSmartPointer<vtkLandmarkTransform> lt =
        vtkSmartPointer<vtkLandmarkTransform>::New();
    lt->SetTargetLandmarks(cP);
    lt->SetSourceLandmarks(mP);
    vtkMatrix4x4 *mat=vtkMatrix4x4::New();
    lt->GetMatrix(mat); //calculates the transform
    return mat;
}

//change InternalUpdate to allow custom initial matrix
class vtkICP: public vtkIterativeClosestPointTransform
{
public:
  static vtkICP *New();
  vtkTypeMacro(vtkICP,vtkIterativeClosestPointTransform);
protected:
  vtkICP() :vtkIterativeClosestPointTransform() {}
  ~vtkICP() {}
//ignore centroids and use matrix from LandmarkTransform
void InternalUpdate()
{
  // Check source, target

  if (this->Source == NULL || !this->Source->GetNumberOfPoints())
    {
    vtkErrorMacro(<<"Can't execute with NULL or empty input");
    return;
    }

  if (this->Target == NULL || !this->Target->GetNumberOfPoints())
    {
    vtkErrorMacro(<<"Can't execute with NULL or empty target");
    return;
    }

  // Create locator

  this->CreateDefaultLocator();
  this->Locator->SetDataSet(this->Target);
  this->Locator->SetNumberOfCellsPerBucket(1);
  this->Locator->BuildLocator();

  // Create two sets of points to handle iteration

  int step = 1;
  if (this->Source->GetNumberOfPoints() > this->MaximumNumberOfLandmarks)
    {
    step = this->Source->GetNumberOfPoints() / this->MaximumNumberOfLandmarks;
    vtkDebugMacro(<< "Landmarks step is now : " << step);
    }

  vtkIdType nb_points = this->Source->GetNumberOfPoints() / step;

  vtkPoints *points1 = vtkPoints::New();
  points1->SetNumberOfPoints(nb_points);

  vtkPoints *closestp = vtkPoints::New();
  closestp->SetNumberOfPoints(nb_points);

  vtkPoints *points2 = vtkPoints::New();
  points2->SetNumberOfPoints(nb_points);

  // Fill with initial positions (sample dataset using step)

  vtkTransform *accumulate = vtkTransform::New();
  accumulate->PostMultiply();

  vtkIdType i;
  int j;
  double p1[3], p2[3];

  //here starts the initialization with custom matrix
  accumulate->SetMatrix(LandmarkTransform->GetMatrix());
  accumulate->Update();

  for (i = 0, j = 0; i < nb_points; i++, j += step)
  {
    double outPoint[3];
    accumulate->InternalTransformPoint(this->Source->GetPoint(j), outPoint);
    points1->SetPoint(i, outPoint);
  }


  // Go
  
  vtkIdType cell_id;
  int sub_id;
  double dist2, totaldist = 0;
  double outPoint[3];

  vtkPoints *temp, *a = points1, *b = points2;

  this->NumberOfIterations = 0;

  do 
    {
    // Fill points with the closest points to each vertex in input

    for(i = 0; i < nb_points; i++)
      {
      this->Locator->FindClosestPoint(a->GetPoint(i),
                                      outPoint,
                                      cell_id,
                                      sub_id,
                                      dist2);
      closestp->SetPoint(i, outPoint);
      }
    
    // Build the landmark transform

    this->LandmarkTransform->SetSourceLandmarks(a);
    this->LandmarkTransform->SetTargetLandmarks(closestp);
    this->LandmarkTransform->Update();

    // Concatenate (can't use this->Concatenate directly)
    
    accumulate->Concatenate(this->LandmarkTransform->GetMatrix());
  
    this->NumberOfIterations++;
    vtkDebugMacro(<< "Iteration: " << this->NumberOfIterations);
    if (this->NumberOfIterations >= this->MaximumNumberOfIterations) 
      {
      break;
      }

    // Move mesh and compute mean distance if needed

    if (this->CheckMeanDistance)
      {
      totaldist = 0.0;
      }

    for(i = 0; i < nb_points; i++)
      {
      a->GetPoint(i, p1);
      this->LandmarkTransform->InternalTransformPoint(p1, p2);
      b->SetPoint(i, p2);
      if (this->CheckMeanDistance)
        {
        if (this->MeanDistanceMode == VTK_ICP_MODE_RMS) 
          {
          totaldist += vtkMath::Distance2BetweenPoints(p1, p2);
          } else {
          totaldist += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
          }
        }
      }
    
    if (this->CheckMeanDistance)
      {
      if (this->MeanDistanceMode == VTK_ICP_MODE_RMS) 
        {
        this->MeanDistance = sqrt(totaldist / (double)nb_points);
        } else {
        this->MeanDistance = totaldist / (double)nb_points;
        }
      vtkDebugMacro("Mean distance: " << this->MeanDistance);
      if (this->MeanDistance <= this->MaximumMeanDistance)
        {
        break;
        }
      }

    temp = a;
    a = b;
    b = temp;

    } 
  while (1);

  // Now recover accumulated result

  this->Matrix->DeepCopy(accumulate->GetMatrix());

  accumulate->Delete();
  points1->Delete();
  closestp->Delete();
  points2->Delete();
}
private:
  vtkICP(const vtkICP&);  // Not implemented.
  void operator=(const vtkICP&);  // Not implemented.
};

vtkStandardNewMacro(vtkICP);

//for aligning the final result with a healthy model for comparison
//or for getting correct orientation so shape model can be used
vtkMatrix4x4* Vertebra::calcICP(bool calcScaling)
{
    vtkMatrix4x4 *temp;
    vtkSmartPointer<vtkICP> icp=vtkSmartPointer<vtkICP>::New();
    getPoly();
    
    icp->SetSource(polyModel);
    icp->SetTarget(poly);
	icp->SetCheckMeanDistance(1);

    temp=rotate(vec3(0,0,1), axis);
    vtkSmartPointer<vtkTransform> t=vtkSmartPointer<vtkTransform>::New();
    t->PostMultiply();
    t->SetMatrix(temp);
    t->Scale(radius/modelRadius, radius/modelRadius, radius/modelRadius);
    t->Translate(center.ptr());

    icp->GetLandmarkTransform()->GetMatrix()->DeepCopy(t->GetMatrix());
    if (!calcScaling)
        icp->GetLandmarkTransform()->SetModeToRigidBody();
    
    icp->GetMatrix(temp); //execute ICP
    return temp;
}

vtkActor* Vertebra::axis2vtk(bool useMRM)
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

    if (useMRM)
    {
        vtkSmartPointer<vtkTransform> t=vtkSmartPointer<vtkTransform>::New();
        t->RotateY(-90);
        t->Scale(20,20,20);
        t->PostMultiply();
        t->Concatenate(mrm);
        actor->SetUserTransform(t);
    }
    else
    {
        vec3 arrow(1,0,0), rotAxis=arrow^axis; //cross product is rotation axis
        double angle=acos(arrow*axis/(arrow.length()*axis.length())); //angle from dot product

        actor->RotateWXYZ(angle*180/M_PI, rotAxis[0], rotAxis[1], rotAxis[2]);
        actor->SetScale(axis.length()*1.2);
        actor->SetPosition(center[0], center[1], center[2]);
    }
    return actor;
}

////voxelization using ITK implementation (bugged)
//VisualizingImageType::Pointer Vertebra::getMask()
//{
//    MeshType::Pointer mesh=qe2itkMesh(qe);
//    typedef itk::TriangleMeshToBinaryImageFilter<MeshType,VisualizingImageType> MeshFilterType;
//    MeshFilterType::Pointer meshFilter = MeshFilterType::New();
//    meshFilter->SetInfoImage(logic.visualizing);
//    meshFilter->SetInput(mesh);
//    meshFilter->Update();
//    mask=meshFilter->GetOutput();
//    return mask;
//}

//voxelization using VTK implementation (workaround)
VisualizingImageType::Pointer Vertebra::getMask()
{
    VisualizingImageType::Pointer whiteITK=VisualizingImageType::New();
    //whiteITK->CopyInformation(logic.visualizing);
    whiteITK->SetRegions(logic.visualizing->GetLargestPossibleRegion());
    whiteITK->Allocate();
    whiteITK->FillBuffer(1);

    typedef itk::ImageToVTKImageFilter<VisualizingImageType> itkVtkConverter;
    itkVtkConverter::Pointer conv=itkVtkConverter::New();
    conv->SetInput(whiteITK);
    conv->Update();
    vtkSmartPointer<vtkImageData> whiteVTK=conv->GetOutput();

    vtkSmartPointer<vtkPolyData> poly=qe2vtk(qe);

    VisualizingImageType::DirectionType d=logic.visualizing->GetDirection();
    vtkSmartPointer<vtkMatrix4x4> Mitk=vtkSmartPointer<vtkMatrix4x4>::New();
    for (int i=0; i<3; i++)
        for (int k=0; k<3; k++)
            Mitk->SetElement(i,k, d(i,k));

    vtkSmartPointer<vtkMatrix4x4> Ms=vtkSmartPointer<vtkMatrix4x4>::New();
    Ms->SetElement(0,0,logic.visualizing->GetSpacing()[0]);
    Ms->SetElement(1,1,logic.visualizing->GetSpacing()[1]);
    Ms->SetElement(2,2,logic.visualizing->GetSpacing()[2]);

    vtkMatrix4x4::Multiply4x4(Mitk,Ms,Mitk);
    VisualizingImageType::PointType origin=logic.visualizing->GetOrigin();    
    for (int i=0; i<3; i++)
        Mitk->SetElement(i,3, origin[i]);
    Mitk->Invert();

    vtkSmartPointer<vtkTransform> t=vtkSmartPointer<vtkTransform>::New();
    t->SetMatrix(Mitk);
    vtkSmartPointer<vtkTransformFilter> tf=vtkSmartPointer<vtkTransformFilter>::New();
    tf->SetInputData(poly);
    tf->SetTransform(t);
    tf->Update();
    poly->SetPoints(tf->GetOutput()->GetPoints());

    vtkSmartPointer<vtkPolyDataToImageStencil> voxelizer =
        vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    voxelizer->SetInputData(poly);
    voxelizer->SetOutputWholeExtent(whiteVTK->GetExtent());
    voxelizer->Update();

    vtkImageStencil *stencil=vtkImageStencil::New(); //crashes with smart pointer
    //vtkSmartPointer<vtkImageStencil> stencil=vtkSmartPointer<vtkImageStencil>::New();
    stencil->SetInputData(conv->GetOutput());
    stencil->SetStencilConnection(voxelizer->GetOutputPort());
    stencil->SetBackgroundValue(0);
    stencil->Update();
    whiteVTK=stencil->GetOutput();

    typedef itk::VTKImageToImageFilter<VisualizingImageType> vtk2itkConverter;
    vtk2itkConverter::Pointer conv2=vtk2itkConverter::New();
    conv2->SetInput(whiteVTK);
    conv2->Update();
    mask=conv2->GetOutput();
    mask->SetDirection(logic.visualizing->GetDirection());
    mask->SetSpacing(logic.visualizing->GetSpacing());
    mask->SetOrigin(logic.visualizing->GetOrigin());
    return mask;
}

vtkPolyData* Vertebra::getPoly()
{
    if (poly)
        poly->Delete();
    poly=qe2vtk(qe);
    return poly;
}

double Vertebra::calcAvgDistance() //average distance of vertices to center
{
    CellVertexIterator cvi(qe);
    Vertex *v;
    radius=0;
    while ( (v=cvi.next())!=0 )
        radius+=(center-v->pos).length();
    radius/=qe->countVertices();
    return radius;
}

//calculates center, volume and surface area
void Vertebra::calcNewCenter()
{
    double tv;
    volume=0;
    surface=0;
    vec3 newC(0,0,0), tc, a,b,c;
    CellFaceIterator cfi(qe);
    Face *f;
    while ( (f=cfi.next()) != 0 )
    {
        FaceEdgeIterator fei(f);
        Edge *e=fei.next();
        a=e->Org()->pos;
        e=fei.next();
        b=e->Org()->pos;
        e=fei.next();
        c=e->Org()->pos;
        tv=(a-center) * ((b-center)^(c-center)) / 6.0;
        //volume of tetrahedron consisting of this face and the old center
        tc=(a+b+c+center)/4.0; //center of mass of this tetrahedron
        newC+=(tc*tv);
        volume+=tv;
        surface+=((b-a)^(c-a)).length()/2; //vector product is double surface area
    }
    center=newC/volume;
}

vec3 operator*(vtkMatrix4x4* matrix, const vec3 point)
{
    float p4[4];
    p4[0]=point[0];
    p4[1]=point[1];
    p4[2]=point[2];
    p4[3]=1;
    matrix->MultiplyPoint(p4,p4);
    return vec3(p4);
}

void Vertebra::addMask()
{
    typedef itk::AddImageFilter<VisualizingImageType> addType;
    addType::Pointer addFilter=addType::New();
    addFilter->SetInput1(logic.masks);
    addFilter->SetInput2(mask);
    addFilter->InPlaceOn();
    addFilter->Update();
    logic.masks=addFilter->GetOutput();
    logic.mainForm.updateOverlay();
}

inline bool vecLength(const vec3& a, const vec3& b)
{
    return a.length2()>b.length2();
}

//fit average vertebra shape and extract endplate vector from it
vec3 Vertebra::calcAxis(vec3 lastAxis)
{
    vtkSmartPointer<vtkTransform> t=vtkSmartPointer<vtkTransform>::New();
    axis=lastAxis; //initialize ICP with this orientation
    vtkMatrix4x4* m=calcICP(true);
    t->SetMatrix(m);
    float temp[4];
    t->GetPosition(temp);
    center[0]=temp[0];
    center[1]=temp[1];
    center[2]=temp[2];
    t->PostMultiply();
    t->Translate(-center[0], -center[1], -center[2]);
    float *tmp=t->TransformFloatPoint(0,0,1);
    axis[0]=tmp[0];
    axis[1]=tmp[1];
    axis[2]=tmp[2];
    axis.normalize();
    return axis;
}