#pragma once
#include "declarations.h"
#include "cell.hh"
#include <vtkMatrix4x4.h>
#include <vtkPropAssembly.h>
#include <vtkActor.h>
#include <QMutex>
#include <QWaitCondition>

class MainLogic;

class Vertebra
{
public:
    static Cell *model; //shape model
    static double modelRadius;
    static vtkPolyData *polyModel;

    QMutex renderLock;
    Cell *qe; //quad-edge mesh of this vertebral body
    vtkActor *actor; //for visualization
    double radius, surface, volume; //parameters of this vb
    InternalImageType::IndexType centerIndex; //to be set before the segmentation comences
    InternalImageType::PixelType highThreshold, lowThreshold, avgValue, sigmaValue;
    vec3 center, axis;

    vtkPolyData *poly; //vtk polygonal mesh for display
    vtkMatrix4x4 *mrm; //model registration matrix

    int labelIndex;
    static const char* labels[];

    //per-vertebra diagnosis result
    int crushedFlag, spondylFlag;
    double crushed, spondyl;

    InternalImageType::Pointer cl_dfTh, cl_value, clCombined;
    VisualizingImageType::Pointer mask;

    MainLogic &logic;
    static QMutex maskLock;
    //methods
    Vertebra(MainLogic &ml)
    :logic(ml)
    {
        axis=vec3(0,0,0);
        poly=0;
        mrm=0; //maybe identity?
        //labelIndex=-1; //invalid index
        //radius=surface=volume=0;

        if (model==0)
        {
            model=objReadCell(vertebraMeshFilename);
            calcVertexValence(model, false);
            copyPos2SubdivPos(model);
            for (unsigned i=1; i<=4; i++)
                topologicallySubdivide(model, i);
            calculateSubdivisionPositions(model, 4, 0);
            copySubdivPos2Pos(model);
            qe=model; modelRadius=calcCurrentRadius(); qe=0; //avgDist uses qe
            polyModel=qe2vtk(model);
        }
        
    }

    VisualizingImageType::Pointer getVisualizationMask();
    VisualizingImageType::Pointer getOriginalMask();
    VisualizingImageType::Pointer getMaskFromImageInformation(
        VisualizingImageType::RegionType region,
        VisualizingImageType::SpacingType spacing,
        VisualizingImageType::PointType origin,
        VisualizingImageType::DirectionType direction);
    void addMask();
    vtkActor* axis2vtk(bool useMRM);    
    vtkPolyData *getPoly();
    vtkMatrix4x4* calcRigidTransform(unsigned upToLevel);
    vtkMatrix4x4* calcICP(bool calcScaling);
    double calcCurrentRadius();
    void calcNewCenter();
    vec3 calcAxis(vec3 lastAxis=vec3(0,0,1));
};