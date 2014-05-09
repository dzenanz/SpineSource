//declarations of globally used types and entry points for top-level functions
#pragma once

#include "obj.hh" //includes all QuadEdge related stuff
#include <string>
#include <array>
#include <Eigen/Eigen>

#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkMesh.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCommand.h"
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <QProgressDialog>
#include <QApplication>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

const qreal PointRadius = 5;
extern const char *vertebraMeshFilename;

typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Image<unsigned char, 2> SliceImageType;
typedef itk::Image<float, 2> InternalSliceImageType;
typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<InternalImageType::OffsetType, 3> VectorMapType;
typedef itk::Image<itk::CovariantVector<float, 3>, 3> GradientImageType;
typedef itk::Mesh<float,3> MeshType;

class Log1Functor
{
public:
   float operator()( float input )
   {
	   return log(1.0+input);
   }
  unsigned long long operator()( unsigned long long input )
   {
	   return log(1.0+input);
   }
};
typedef itk::UnaryFunctorImageFilter<InternalImageType, InternalImageType, Log1Functor> Log1Type;

void polynomialFit(unsigned degree, const std::vector<class Vertebra *> vertebra, Eigen::VectorXd &x, Eigen::VectorXd &y, bool ignoreOrientations=true);
vtkActor * centerline(Eigen::VectorXd &x, Eigen::VectorXd &y, double fromZ, double toZ, const vec3 color);
void diagnose(const std::vector<Vertebra *> vertebra, const Eigen::VectorXd x, const Eigen::VectorXd y, std::string filename_base);

unsigned int exp2i(int exponent); //2^exp, returns at least 1, even if argument is negative
vec3 operator*(vtkMatrix4x4* matrix, const vec3 point);
void applyMatrix(Cell *qe, vtkMatrix4x4 *matrix);
void writeSTL(vtkPolyData* poly, const char * fileName);

//conversions
vtkPolyData* qe2vtk(Cell *qe);
vtkActor *poly2actor(vtkPolyData *poly, bool wireframe, vec3 color=vec3(1,1,0)); //default yellow
MeshType::Pointer qe2itkMesh(Cell *cell);
void qePhysicalToIndex(Cell *qe, VisualizingImageType::Pointer image);
void qeIndexToPhysical(Cell *qe, VisualizingImageType::Pointer image);
void itkMeshIndexToPhysical(MeshType::Pointer mesh, VisualizingImageType::Pointer image);

void translate(Cell *qe, const vec3& vector);
vtkMatrix4x4* rotate(vec3 from, vec3 to);
void rotate(Cell *qe, vec3 from, vec3 to);
void calcFaceNormals(Cell *s);
float getAverageEdgeLength(Cell *qe, float &minLen, float &maxLen, float &stdDev);
void calcVertexValence(Cell *qe, bool showResult=true);

//splits edge e at midpoint, returns enew (the rear part of original edge)
//enewOrg->enewDest==eOrg->eDest
Edge *split(Edge *e);

void calcLHvalues(InternalImageType::Pointer in, const float eps, InternalImageType::Pointer l, InternalImageType::Pointer h);
void calc2DJointHistogram(VisualizingImageType::Pointer x_aka_l, VisualizingImageType::Pointer y_aka_h, std::string savefilename);

//smooth interpolating subdivision using modified butterfly scheme
void topologicallySubdivide(Cell *qe, unsigned level);
void calculateSubdivisionPositions(Cell *qe, unsigned maxLevel, unsigned subdivFromLevel);
void normalizeHeuristic(Cell *qe, int maxLevel, int howManyFinestLevels);
void copyPos2SubdivPos(Cell *qe);
void copySubdivPos2Pos(Cell *qe);
void swapSubdivPos(Cell *qe);
Eigen::MatrixXd createWeightMatrix(Cell *qe, unsigned maxLevel, unsigned subdivFromLevel, unsigned &freeCount);

//QtProgressDialog
class FilterProgress : public itk::Command 
{
public:
  typedef  FilterProgress		   Self;
  typedef  itk::Command            Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  QProgressDialog *progress;
  static const int maximum=1000;

  FilterProgress()
  {
	progress=new QProgressDialog("Filter progress", QString(), 0, maximum);
	progress->setWindowModality(Qt::ApplicationModal);
	progress->setValue(0);
  };

  ~FilterProgress()
  {
	  progress->close();
	  delete progress;
  }
public:
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
	if( typeid( itk::ProgressEvent )   ==  typeid( event ) )
		{
		const itk::ProcessObject * process = dynamic_cast< const itk::ProcessObject * >( object );
		const int value = static_cast<int>( process->GetProgress() * maximum );
		progress->setValue( value );
		QApplication::processEvents();
		}
    }
}; //FilterProgress