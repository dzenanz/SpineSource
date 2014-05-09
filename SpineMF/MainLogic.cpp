#include "MainLogic.h"

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iomanip>

#include "itksys/SystemTools.hxx"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkLineIterator.h"
#include "itkGradientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkCastImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCustomCanny.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkVectorMagnitudeImageFilter.h"

#include <vtkMatrix4x4.h>
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
#include <vtkSphereSource.h>
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

#include "vec3.h"
#include "itkStructureTensorImageFilter.h"
#include "itkBinaryClosingByReconstructionImageFilter.h"
#include "itkFastBilateralImageFilter.h"
#include "voxelClassification.h"

#include <QPainter>
#include <QBitmap>
#include <QFile>
#include <QDir>
#include <qmessagebox.h>


using namespace std;

MainLogic::MainLogic(MainWindow& mainWindow)
:mainForm(mainWindow)
{
    connect(&mainForm, SIGNAL(volumeOpened(std::string)), SLOT(volumeOpen(std::string)));
	connect(mainForm.painter, SIGNAL(approximateBoundaryDrawn(int,int)), SLOT(approximateBoundaryDrawn(int,int)));
    connect(mainForm.painter, SIGNAL(singlePointPicked(int,int)), SLOT(singlePointPicked(int,int)));
	connect(mainForm.painter, SIGNAL(vertebraCentersPicked(int,int)), SLOT(vertebraCentersPicked(int,int)));
	connect(mainForm.painter, SIGNAL(threePointsPicked(int,int)), SLOT(threePointsPicked(int,int)));
    connect(&mainForm, SIGNAL(dicomOpened(std::string)), SLOT(dicomOpen(std::string)));
    connect(mainForm.painter, SIGNAL(majorUpdate(int)), SLOT(approximateBoundaryDrawn(int)));
    current=0;
    visualizing=0;
    vertebra1=0;
    masks=0;
    vi=0;
    maxValue=800; //reasonable initial guess

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
const float percentMoreInflated=0.00; //0.02=>2%
const unsigned heuristicLevels=2;
//disease thresholds:
double crushedVertebraThreshold=0.2;
double significantSpondylolisthesis=0.25;

template <class T>
inline bool withinEpsilonEnvironment(T point, T x, float eps)
{
    return x>=point-point*eps && x<=point+point*eps;
}

inline InternalImageType::IndexType vec3ToIndex(vec3 point, InternalImageType::Pointer image)
{
    InternalImageType::IndexType ind;
    InternalImageType::PointType p;
    p[0]=point[0]; p[1]=point[1]; p[2]=point[2];
    image->TransformPhysicalPointToIndex(p, ind);
    return ind;
}

void MainLogic::dicomOpen(std::string directory)
{
    try
    {
        mainForm.statusbar->showMessage("Reading dataset: "+QString::fromStdString(directory));

        typedef itk::GDCMImageIO ImageIOType;
        ImageIOType::Pointer dicomIO = ImageIOType::New();
        typedef itk::ImageSeriesReader< InternalImageType > ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetImageIO(dicomIO);

        typedef itk::GDCMSeriesFileNames NamesGeneratorType;
        NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
        nameGenerator->SetUseSeriesDetails(true);
        nameGenerator->SetDirectory(directory.c_str());

        typedef std::vector< std::string > SeriesIdContainer;
        const SeriesIdContainer& seriesUID = nameGenerator->GetSeriesUIDs();

        std::string seriesIdentifier;
        seriesIdentifier = seriesUID.begin()->c_str(); //Choose First Series Found In Directory

        typedef std::vector<std::string> FileNamesContainer;
        FileNamesContainer fileNames;
        fileNames = nameGenerator->GetFileNames(seriesIdentifier);

        reader->SetFileNames(fileNames);
        reader->Update();
        current=reader->GetOutput();

        volume_filename=directory;
        afterOpen(dicomIO->GetMetaDataDictionary());
    }
    catch (std::exception& e)
    {
        QMessageBox mb(QMessageBox::Critical, QString("Error opening DICOM series"),
            QString("Error opening DICOM series in directory ")+QString::fromStdString(directory)+"\nException content:\n"+e.what(),
            QMessageBox::Ok);
        mb.exec();
    }
}

void MainLogic::volumeOpen(std::string filename)
{
    try
    {
        mainForm.statusbar->showMessage("Reading dataset: "+QString::fromStdString(filename));
        typedef  itk::ImageFileReader< InternalImageType >  ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( filename.c_str() );
        reader->Update();
        current=reader->GetOutput();
        volume_filename=filename;
        afterOpen(current->GetMetaDataDictionary());
    }
    catch (std::exception& e)
    {
        QMessageBox mb(QMessageBox::Critical, QString("Error opening file"),
            QString("Error opening file ")+QString::fromStdString(filename)+"\nException content:\n"+e.what(),
            QMessageBox::Ok);
        mb.exec();
    }
}

void writeImage(InternalImageType::Pointer image, std::string fileName, bool compressed=false)
{
    typedef itk::ImageFileWriter<InternalImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    writer->SetUseCompression(compressed);
    writer->Update();
}

void writeImage(VisualizingImageType::Pointer image, std::string fileName, bool compressed=true)
{
    typedef itk::ImageFileWriter<VisualizingImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    writer->SetUseCompression(compressed);
    writer->Update();
}

void MainLogic::saveLHfiles(int i)
{
    typedef itk::MinimumMaximumImageCalculator< InternalImageType >  CalculatorType;
    CalculatorType::Pointer calculatorL = CalculatorType::New();
    calculatorL->SetImage( lImage );
    calculatorL->ComputeMaximum();
    maxValue = max(calculatorL->GetMaximum(), maxValue);

    CalculatorType::Pointer calculatorH = CalculatorType::New();
    calculatorH->SetImage( hImage );
    calculatorH->ComputeMaximum();
    maxValue = max(calculatorH->GetMaximum(), maxValue);
    while (pow(2.0f,i)-1<maxValue)
        i++;
    maxValue2=pow(2.0f,i)-1;

    VisualizingImageType::Pointer lVis, hVis;

    typedef itk::ShiftScaleImageFilter < InternalImageType, VisualizingImageType > RescaleImageFilterType;
    RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
    rescale->SetInput( lImage );
    rescale->SetScale(255/maxValue);
    rescale->Update();
    lVis=rescale->GetOutput();

    RescaleImageFilterType::Pointer rescale2 = RescaleImageFilterType::New();
    rescale2->SetInput( hImage );
    rescale2->SetScale(255/maxValue);
    rescale2->Update();
    hVis=rescale2->GetOutput();

    std::string fnNoExt=itksys::SystemTools::GetFilenamePath(volume_filename)+'/'
        +itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);

    writeImage(lImage, fnNoExt+"_L.mha");
    writeImage(hImage, fnNoExt+"_H.mha");
    calc2DJointHistogram(lVis, hVis, fnNoExt+"_LH.png");
}

class DicomWindowRescaler
{
#define ymax (1024)
public:
    float center, width;
    float operator()( float x )
    {
        if	(x <= center - 0.5 - (width-1)/2)
            return 0;
        else if	(x > center - 0.5 + (width-1)/2)
            return ymax;
        else
            return ((x - (center - 0.5)) / (width-1) + 0.5) * ymax;
    }
};

//apply denoising, calculate LH values and distance field (DF) from gradient magnitude
void MainLogic::afterOpen(const itk::MetaDataDictionary &metaData)
{
    filename_base=tempDir+itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);
    if (masks)
        masks=0;

    series_description="";
    double windowWidth=0, windowCenter=0;
    typedef itk::MetaDataObject<std::string> MetaDataStringType;
    itk::MetaDataDictionary::ConstIterator itr = metaData.Begin();
    itk::MetaDataDictionary::ConstIterator end = metaData.End();

    //check presence of the metadata and use it if found
    while(itr!=end)
    {
        MetaDataStringType::Pointer entryValue =
            dynamic_cast<MetaDataStringType *>(itr->second.GetPointer());
        if (entryValue)
        {
            if (itr->first == "0008|103e") //Series Description
                series_description = entryValue->GetMetaDataObjectValue();
            if (itr->first == "0028|1050") //Window Center
                windowCenter=QString::fromStdString(entryValue->GetMetaDataObjectValue()).toDouble();
            if (itr->first == "0028|1051") //Window Width
                windowWidth=QString::fromStdString(entryValue->GetMetaDataObjectValue()).toDouble();
        }
        ++itr;
    }

    if (windowCenter!=0 && windowWidth!=0 && mainForm.actionUse_DICOMs_windowing->isChecked())
    {
        typedef itk::UnaryFunctorImageFilter<InternalImageType, InternalImageType,
            DicomWindowRescaler> WindowRescalerFilter;
        WindowRescalerFilter::Pointer wr=WindowRescalerFilter::New();
        wr->SetInput(current);
        wr->GetFunctor().center=windowCenter;
        wr->GetFunctor().width=windowWidth;
        wr->Update();
        current=wr->GetOutput();
        if (mainForm.actionSave_debug_images->isChecked())
            writeImage(current, tempDir+"dicomWindow.mha");
    }
    
    typedef itk::MinimumMaximumImageCalculator< InternalImageType >  CalculatorType;
    CalculatorType::Pointer calculatorI = CalculatorType::New();
    calculatorI->SetImage( current );
    calculatorI->ComputeMaximum();
    maxValueOriginal = calculatorI->GetMaximum();

    mainForm.on_actionMedian_denoising_triggered();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(current, tempDir+"median.mha");

    calculatorI->SetImage( current );
    calculatorI->ComputeMaximum();
    maxValue = calculatorI->GetMaximum(); //needed to get parameters right for bilateral filter

    mainForm.statusbar->showMessage("Bilateral filtering");
    typedef itk::FastBilateralImageFilter<InternalImageType, InternalImageType> BilateralType;
    BilateralType::Pointer bilateral=BilateralType::New();
    bilateral->SetInput(current);
    bilateral->SetDomainSigma(10); //default 4 gives very similar result
    bilateral->SetRangeSigma(0.03*maxValue);
    //FilterProgress::Pointer progress=FilterProgress::New(); //progress bar not needed for fast bilateral
    //bilateral->AddObserver(itk::ProgressEvent(),  progress);
    bilateral->Update();
    current=bilateral->GetOutput();
    calculatorI->SetImage( current );
    calculatorI->ComputeMaximum();
    maxValue = calculatorI->GetMaximum();
    mainForm.updateVisualization();	
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(current, tempDir+"bilateral.mha");
    
    int i=1;
    while (pow(2.0f,i)-1<maxValue)
        i++;
    maxValue2=pow(2.0f,i)-1;

    mainForm.statusbar->showMessage("Calculating LH values");
    lImage=InternalImageType::New();
    hImage=InternalImageType::New();
    calcLHvalues(current, maxValue/100, lImage, hImage); //epsilon = 1%
    if (mainForm.actionSave_LH_images_and_histogram->isChecked())
        saveLHfiles(i);

    mainForm.statusbar->showMessage("Calculating structure tensors");
    typedef itk::itkStructureTensorImageFilter<InternalImageType> TensorFilterType;
    TensorFilterType::Pointer filter=TensorFilterType::New();
    filter->SetInput(current);
    filter->SetSigma(1.0);
    filter->Update();
	    
    mainForm.statusbar->showMessage("Conducting eigen-analysis");
    typedef itk::Image<TensorFilterType::CVector, 3> CovariantImageType;
    typedef itk::SymmetricEigenAnalysisImageFilter<TensorFilterType::OutputImageType,
        CovariantImageType> EigenFilterType;
    EigenFilterType::Pointer ef=EigenFilterType::New();
    ef->SetInput(filter->GetOutput());
    ef->SetDimension(3);
    ef->OrderEigenValuesBy(EigenFilterType::FunctorType::OrderByMagnitude);
    ef->Update();

    typedef itk::UnaryFunctorImageFilter<CovariantImageType, InternalImageType,
        EigenvaluesToSurfel> Eigen2TypeFilter;
    Eigen2TypeFilter::Pointer e2t=Eigen2TypeFilter::New();
    e2t->SetInput(ef->GetOutput());
    e2t->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(e2t->GetOutput(), tempDir+"eigen2surfel.mha");

    mainForm.statusbar->showMessage("Calculating gradient");
    typedef itk::GradientMagnitudeImageFilter<InternalImageType, InternalImageType> GtMtype;
    GtMtype::Pointer gtmFilter=GtMtype::New();
    gtmFilter->SetInput(current);
    gtmFilter->Update();
    InternalImageType::Pointer gm=gtmFilter->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(gm, tempDir+"gm.mha");

    CalculatorType::Pointer calculatorGM = CalculatorType::New();
    calculatorGM->SetImage(gm);
    calculatorGM->ComputeMaximum();
    maxGM=calculatorGM->GetMaximum();

    typedef itk::MultiplyImageFilter<InternalImageType> multType;
    multType::Pointer mult=multType::New();
    mult->SetInput1(e2t->GetOutput());
    mult->SetInput2(gm);
    mult->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(mult->GetOutput(), tempDir+"mult.mha");

    CalculatorType::Pointer calculatorMult = CalculatorType::New();
    calculatorMult->SetImage( mult->GetOutput() );
    calculatorMult->ComputeMaximum();

    typedef itk::BinaryThresholdImageFilter<InternalImageType, VisualizingImageType> BinarizorType;
    BinarizorType::Pointer binMult=BinarizorType::New();
    binMult->SetInput(mult->GetOutput());
    binMult->SetLowerThreshold(calculatorMult->GetMaximum()/10);
    binMult->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(binMult->GetOutput(), tempDir+"multBin.mha");

    typedef itk::ShiftScaleImageFilter < InternalImageType, InternalImageType > RescaleImageFilterType;
    RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
    rescale->SetInput( mult->GetOutput() );
    rescale->SetScale(1/calculatorMult->GetMaximum());
    rescale->Update();

    mainForm.statusbar->showMessage("Canny edge detection");
    typedef itk::CustomCanny<InternalImageType, InternalImageType> cannyType;
    cannyType::Pointer cannyFilter=cannyType::New();
    cannyFilter->SetVariance(2);
    cannyFilter->SetUpperThreshold(0.15);
    cannyFilter->SetLowerThreshold(0.05);
    cannyFilter->SetInput(current);
    cannyFilter->SetInput(1, rescale->GetOutput());
    cannyFilter->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(cannyFilter->GetOutput(), tempDir+"canny.mha");

    typedef itk::BinaryThresholdImageFilter<InternalImageType, VisualizingImageType> BinarizorType;
    BinarizorType::Pointer bin=BinarizorType::New();
    bin->SetInput(cannyFilter->GetOutput());
    bin->SetLowerThreshold(0.5); //used values are 0.0 and 1.0
    bin->Update();

    mainForm.statusbar->showMessage("Calculating distance field from edges");
    typedef itk::SignedMaurerDistanceMapImageFilter<VisualizingImageType, InternalImageType> DistanceMapType;
    DistanceMapType::Pointer dm=DistanceMapType::New();
    dm->SetInput(bin->GetOutput());
    dm->SetUseImageSpacing(true);
    dm->SquaredDistanceOff();
    dm->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(dm->GetOutput(), tempDir+"cannyDF.mha");

    DistanceMapType::Pointer dmMulti=DistanceMapType::New();
    dmMulti->SetInput(binMult->GetOutput());
    dmMulti->SetUseImageSpacing(true);
    dmMulti->SquaredDistanceOff();
    dmMulti->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(dmMulti->GetOutput(), tempDir+"multiDF.mha");

    mainForm.statusbar->showMessage("Executing classifiers");
    typedef itk::UnaryFunctorImageFilter<InternalImageType,InternalImageType,dfEdgeClassifier> dfClassifier;
    dfClassifier::Pointer dfCl=dfClassifier::New();
    dfCl->SetInput(dm->GetOutput());
    dfCl->Update();
    cl_dfCanny=dfCl->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(cl_dfCanny, tempDir+"clCanny.mha");

    dfClassifier::Pointer dfClMulti=dfClassifier::New();
    dfClMulti->SetInput(dmMulti->GetOutput());
    dfClMulti->Update();
    cl_dfMult=dfClMulti->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(cl_dfMult, tempDir+"clMulti.mha");

    typedef itk::BinaryFunctorImageFilter<InternalImageType, InternalImageType,
        InternalImageType, lhClassifierFunctor> lhClassifier;
    lhClassifier::Pointer lhCl=lhClassifier::New();
    lhCl->GetFunctor().maxVal=maxValue;
    lhCl->SetInput1(lImage);
    lhCl->SetInput2(hImage);
    lhCl->Update();
    cl_LH1=lhCl->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(cl_LH1, tempDir+"clLH1.mha");

    typedef itk::TernaryFunctorImageFilter<InternalImageType, InternalImageType,
        InternalImageType, InternalImageType, lhClassifierFunctor2> lhClassifier2;
    lhClassifier2::Pointer lhCl2=lhClassifier2::New();
    lhCl2->GetFunctor().maxVal=maxValue;
    lhCl2->SetInput1(lImage);
    lhCl2->SetInput2(hImage);
    lhCl2->SetInput3(current);
    lhCl2->Update();
    cl_LH2=lhCl2->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(cl_LH2, tempDir+"clLH2.mha");

    std::vector<InternalImageType::Pointer> classifiers(4);
    classifiers[0]=cl_LH1;
    classifiers[1]=cl_LH2;
    classifiers[2]=cl_dfCanny;
    classifiers[3]=cl_dfMult;
    InternalImageType::Pointer clCombinedDS=combineClassifiers(classifiers, mean);
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(clCombinedDS, tempDir+"clCombinedDS.mha");

    mainForm.statusbar->showMessage("Ready");
    mainForm.setWindowTitle(QString("Spine Analyzer - ")+QString::fromStdString(volume_filename));
    QApplication::processEvents();
}

inline double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//returns true if the given position is within image (in that case also sets val to that voxel's value)
bool pos2val(InternalImageType::Pointer image, vec3& pos, InternalImageType::PixelType &val)
{
    itk::Point<float, 3> p(pos.ptr());
    typedef itk::LinearInterpolateImageFunction<InternalImageType,float> InterpolatorType;
    static InterpolatorType::Pointer interp=InterpolatorType::New();
    interp->SetInputImage(image);
    InterpolatorType::ContinuousIndexType ind;
    if (image->TransformPhysicalPointToContinuousIndex<float>(p, ind)) //coordinates within image
    {
        val=interp->EvaluateAtContinuousIndex( ind );
        return true;
    }
    else
        return false;
}

void MainLogic::inflate(float stepSize)
{
    InternalImageType::PixelType pixelFw, pixelBw;
    InternalImageType::SpacingType sp = vertebra[vi]->clCombined->GetSpacing();
    //float minSpacing=min(min(sp[0], sp[1]), sp[2]);
    vec3 posFw, posBw;
    CellVertexIterator it3(vertebra[vi]->qe);
    Vertex *v;
    while ((v = it3.next()) != 0)
    {
        vec3 n=v->pos-vertebra[vi]->center; //center-vertex vector
        n.normalize();

        vec3 mp;
        if (useShape)//use model to steer inflation direction
        {
            mp=vertebra[vi]->mrm*Vertebra::model->vertices[v->realIndex()]->pos;
            //if ((mp-v->pos).length2()>avgSpacing)
            //    clamp to avgSpacing
        }

        if (v->wasDeflated)
            posFw=v->pos+n*stepSize*0.5;
        else
            posFw=v->pos+n*stepSize;
        if (v->wasDeflated)
            posBw=v->pos-n*stepSize;
        else
            posBw=v->pos-n*stepSize*0.5;

        if (useShape)
        {
            posFw=0.5*posFw+0.5*mp; //50% shape model, 50% inflation
            posBw=0.5*posBw+0.5*mp; //50% shape model, 50% inflation
        }

        if (!pos2val(vertebra[vi]->clCombined, posFw, pixelFw)) //if outside
        {
            v->pos=posBw; //deflate
            v->wasDeflated=true;
        }
        if (pixelFw>0.5)
        {
            v->pos=posFw; //inflate
            v->wasDeflated=false;
        }
        else
        {
            if (!pos2val(vertebra[vi]->clCombined, posBw, pixelBw)) //if outside
                return; //do not change wasDeflated property

            if (pixelFw<=pixelBw)
            {
                v->pos=posFw; //inflate
                v->wasDeflated=false;
            }
            else
            {
                v->pos=posBw; //deflate
                v->wasDeflated=true;
            }
        }
    }
}

double calcAvgDistance(const QPainterPath &path, QPointF center)
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

void MainLogic::calcThresholds(const QPainterPath &path)
{
    InternalImageType::IndexType ind;
    double xc=vertebra[vi]->centerIndex[0]*0.2,
           yc=vertebra[vi]->centerIndex[1]*0.2; //move 20% towards center
    int zc=vertebra[vi]->centerIndex[2];
    QPainterPath p(path.elementAt(0)*0.8+QPointF(xc,yc));
    for (int i=1; i<path.elementCount(); i++)
    {
        p.lineTo(path.elementAt(i)*0.8+QPointF(xc,yc));
    }

    QImage pm(current->GetLargestPossibleRegion().GetSize(0), current->GetLargestPossibleRegion().GetSize(1), QImage::Format_ARGB32);
    QPainter pnt(&pm);
    pnt.setPen(Qt::blue);
    pnt.drawPath(p);
    QPen pen(QBrush(Qt::yellow), 3);
    pnt.setPen(pen);
    pnt.drawPath(path);
    pnt.setPen(Qt::red);
    pnt.drawPath(path);
    pnt.end();
    pm.save(QString::fromStdString(filename_base)+"_"+QString::number(mainForm.painter->sliceSlider->value())+".png");

    typedef std::map<InternalImageType::PixelType, unsigned> pCounter;
    pCounter vox;
    unsigned long long count=0;
    InternalImageType::PixelType noiseUnit=(maxValue2+1)/1024;
    QRectF bb=path.boundingRect();
    for (int z=max(0, zc-1); z<=min(zc+1, (int)current->GetLargestPossibleRegion().GetSize(2)-1); z++) //3 slices
        for (int i=bb.top(); i<bb.bottom(); i++)
            for (int k=bb.left(); k<bb.right(); k++)
                if (p.contains(QPoint(k,i)) || z==zc&&path.contains(QPoint(k,i)))
                {
                    ind[0]=k;
                    ind[1]=i;
                    ind[2]=z;
                    InternalImageType::PixelType pix=current->GetPixel(ind);
                    pix=floor(pix/noiseUnit)*noiseUnit;
                    vox[pix]++;
                    count++;
                }

    double sum=0;
    for (pCounter::iterator it=vox.begin(); it!=vox.end();)
        if (it->second<(count/10000.0)+1) //this voxel value is due to noise (less than 0.01%)
        {
            count-=it->second; //reduce the count appropriately
            vox.erase(it++); //remove the voxel
        }
        else
        {
            sum+=it->first*it->second; //val*count
            it++;
        }
    vertebra[vi]->avgValue=sum/count;
    
    sum=0; //now calculate standard deviation
    for (pCounter::iterator it=vox.begin(); it!=vox.end();)
    {
        sum+=it->second*(it->first-vertebra[vi]->avgValue)*(it->first-vertebra[vi]->avgValue);
        it++;
    }
    vertebra[vi]->sigmaValue=sqrt(sum/count);

    vertebra[vi]->highThreshold=vox.rbegin()->first+noiseUnit; //add 1 noiseUnit to round up
    vertebra[vi]->lowThreshold=vox.begin()->first; //already rounded down
}

void MainLogic::recalcThresholds(Cell *s, vec3 center, double lastDist)
{
    typedef std::map<InternalImageType::PixelType, unsigned> pCounter;
    pCounter vox;
    unsigned long long count=0;
    InternalImageType::PixelType noiseUnit=(maxValue2+1)/1024;
    InternalImageType::PointType p;
    p[0]=center[0]; p[1]=center[1]; p[2]=center[2];
    current->TransformPhysicalPointToIndex(p, vertebra[vi]->centerIndex);

    InternalImageType::SpacingType sp = current->GetSpacing();
    InternalImageType::RegionType region;
    region.SetIndex(0, max<double>(0.0, vertebra[vi]->centerIndex[0]-lastDist/(2*sp[0])));
    region.SetIndex(1, max<double>(0.0, vertebra[vi]->centerIndex[1]-lastDist/(2*sp[1])));
    region.SetIndex(2, max<double>(0.0, vertebra[vi]->centerIndex[2]-lastDist/(2*sp[2])));
    region.SetSize(0, min<unsigned long long>(lastDist/sp[0], current->GetLargestPossibleRegion().GetSize(0)-region.GetIndex(0)));
    region.SetSize(1, min<unsigned long long>(lastDist/sp[1], current->GetLargestPossibleRegion().GetSize(1)-region.GetIndex(1)));
    region.SetSize(2, min<unsigned long long>(lastDist/sp[2], current->GetLargestPossibleRegion().GetSize(2)-region.GetIndex(2)));

    if (region.GetSize(0)<=0 || region.GetSize(0)<=0 ||region.GetSize(0)<=0
        || region.GetIndex(0)>=current->GetLargestPossibleRegion().GetSize(0)
        || region.GetIndex(1)>=current->GetLargestPossibleRegion().GetSize(1)
        || region.GetIndex(2)>=current->GetLargestPossibleRegion().GetSize(2))
        return; //region outside of image completely

    Cell *qe=s->deepCopy();
    CellVertexIterator cvi(qe);
    Vertex *v;
    center*=0.2;
    while ( (v=cvi.next())!=0 )
        v->pos=v->pos*0.8+center; //move all vertices 20% inwards
    MeshType::Pointer mesh=qe2itkMesh(qe);
    Cell::kill(qe);

    typedef itk::TriangleMeshToBinaryImageFilter<MeshType,VisualizingImageType> MeshFilterType;
    MeshFilterType::Pointer meshFilter = MeshFilterType::New();
    meshFilter->SetInfoImage(visualizing);
    meshFilter->SetInput(mesh);
    meshFilter->Update();
    binVertebra=meshFilter->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(binVertebra, filename_base+"_reduced"+QString::number(vertebra.size()).toStdString()+".mha");

    typedef itk::ImageRegionConstIterator<InternalImageType> IterType;
    IterType iter(current, region);
    while(!iter.IsAtEnd())
    {
        if (binVertebra->GetPixel(iter.GetIndex()))
        {
            InternalImageType::PixelType pix = iter.Get();
            pix=floor(pix/noiseUnit)*noiseUnit;
            vox[pix]++;
            count++;
        }
        ++iter;
    }

    double sum=0;
    for (pCounter::iterator it=vox.begin(); it!=vox.end();)
        if (it->second<(count/10000.0)+1) //this voxel value is due to noise (less than 0.01%)
        {
            count-=it->second; //reduce the count appropriately
            vox.erase(it++); //remove the voxel
        }
        else
        {
            sum+=it->first*it->second; //val*count
            it++;
        }
    vertebra[vi]->avgValue=sum/count;
    
    sum=0; //now calculate standard deviation
    for (pCounter::iterator it=vox.begin(); it!=vox.end();)
    {
        sum+=it->second*(it->first-vertebra[vi]->avgValue)*(it->first-vertebra[vi]->avgValue);
        it++;
    }
    vertebra[vi]->sigmaValue=sqrt(sum/count);

    vertebra[vi]->lowThreshold=vox.begin()->first;
    vertebra[vi]->highThreshold=vox.rbegin()->first+noiseUnit;
    //InternalImageType::PixelType newHigh=vox.rbegin()->first+noiseUnit, newLow=vox.begin()->first;
    //if (newLow<vertebra[vi]->lowThreshold)
    //    vertebra[vi]->lowThreshold=newLow;
    //if (newHigh>vertebra[vi]->highThreshold)
    //    vertebra[vi]->highThreshold=newHigh;
    //if (newLow==vertebra[vi]->lowThreshold && withinEpsilonEnvironment(vertebra[vi]->highThreshold, newHigh, 0.2))
    //    vertebra[vi]->highThreshold=0.9*vertebra[vi]->highThreshold+0.1*newHigh;
    //if (newHigh==vertebra[vi]->highThreshold && withinEpsilonEnvironment(vertebra[vi]->lowThreshold, newLow, 0.2))
    //    vertebra[vi]->lowThreshold=0.9*vertebra[vi]->lowThreshold+0.1*newLow;
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

//create a matrix of wieghts which contains influence of free vertices on subdivided vertices
void updateWeightMatrix(Eigen::MatrixXd& weights, myLU<Eigen::MatrixXd>& LU,
    Cell *qe, unsigned maxLevel, unsigned freeLevels, unsigned& freeCount)
{
    QString levelInfo=QString::number(maxLevel)+"_"+QString::number(freeLevels);
    QDir cache("./cache_LU");
    if (!cache.exists())
        cache.mkpath(".");
    QFile fw("./cache_LU/W_"+levelInfo+vertebraMeshFilename);
    QFile flu("./cache_LU/LU_"+levelInfo+vertebraMeshFilename);
    QFile fp("./cache_LU/P_"+levelInfo+vertebraMeshFilename);
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
        LU=myLU<Eigen::MatrixXd>(freeCount);//LU.resize(freeCount, freeCount);
#ifdef _DEBUG
        LU.markInitialized();
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

//returns average distance from center to surface
double MainLogic::segmentOneButterfly(double previousDist, bool useMRM, unsigned call)
{
    InternalImageType::SpacingType sp = current->GetSpacing();
    double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);
    double minSpacing=min(min(sp[0], sp[1]), sp[2]);
    vertebra[vi]->qe=objReadCell(vertebraMeshFilename);
    calcVertexValence(vertebra[vi]->qe, false);
    CellVertexIterator it3(vertebra[vi]->qe);
    Vertex *v;
    if (useMRM) //put into right position and orientation
        while ((v = it3.next()) != 0)
            v->pos=vertebra[vi]->mrm*v->pos;
    else
    {
        rotate(vertebra[vi]->qe, vec3(0,0,1), vertebra[vi]->axis);
        translate(vertebra[vi]->qe, vertebra[vi]->center);
    }
    copyPos2SubdivPos(vertebra[vi]->qe);
    vertebra[vi]->actor=0;
    if (mainForm.actionInteractive_inflation->isChecked())
    {
        vertebra[vi]->actor=poly2actor(vertebra[vi]->getPoly(), mainForm.actionShow_wireframe->isChecked());
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[vi]->actor);
        QApplication::processEvents();
    }

    const unsigned dependentLevels=2;
    unsigned maxLevel=dependentLevels, freeCount;
    for (unsigned i=1; i<=maxLevel; i++)
        topologicallySubdivide(vertebra[vi]->qe, i);

    Eigen::Matrix<double,Eigen::Dynamic,3> P;
    Eigen::MatrixXd W;
    myLU<Eigen::MatrixXd> LU;
    P.resize((vertebra[vi]->qe)->countVertices(), 3);
    updateWeightMatrix(W, LU, vertebra[vi]->qe, maxLevel, 0, freeCount);

    calculateSubdivisionPositions(vertebra[vi]->qe, maxLevel,0);
    copySubdivPos2Pos(vertebra[vi]->qe);

    vertebra[vi]->calcAvgDistance();
    vertebra[vi]->calcNewCenter();
    double oldVolume, oldDist2, oldDist=vertebra[vi]->radius;
    vec3 origCenter=vertebra[vi]->center;
    float edgeLen, stdDev, minEdge, maxEdge;
    edgeLen=getAverageEdgeLength(vertebra[vi]->qe, minEdge, maxEdge, stdDev);
    int iter=0, extraIter=5; //extraIter=minimum number of iterations
    //extraIter=35.5/avgSpacing; //average VB size divided by spacing
    int maxIter=250;
    do
    {
        if (useMRM)
        {
            if (vertebra[vi]->mrm)
                vertebra[vi]->mrm->Delete();
            vertebra[vi]->mrm=vertebra[vi]->calcRigidTransform(maxLevel-dependentLevels);
        }

        inflate(min(1.5*minSpacing, minSpacing*0.5 + vertebra[vi]->radius - oldDist));
        //inflate(minSpacing); //constant step size, about twice faster
        {
            CellVertexIterator cvi0(vertebra[vi]->qe);
            while ((v = cvi0.next())!=0)
            {
                P(v->realIndex(),0)=v->pos[0];
                P(v->realIndex(),1)=v->pos[1];
                P(v->realIndex(),2)=v->pos[2];
            }
        }

        oldDist2=oldDist;
        oldDist=vertebra[vi]->radius;
        oldVolume=vertebra[vi]->volume;
        vertebra[vi]->calcAvgDistance();
        vertebra[vi]->calcNewCenter();

        normalizeHeuristic(vertebra[vi]->qe, maxLevel, heuristicLevels); //implicit smoothing
        //heuristic normalization for a few finest levels, then optimal for all levels
        {
            Eigen::Matrix<double,Eigen::Dynamic,3> N;
            N=LU.solve(W*P); //N=WWTinv*(W*P);
            CellVertexIterator cvi(vertebra[vi]->qe);
            while ((v = cvi.next())!=0)
                if (v->level==0)
                {
                    v->pos[0]=N(v->realIndex(),0);
                    v->pos[1]=N(v->realIndex(),1);
                    v->pos[2]=N(v->realIndex(),2);
                }
            calculateSubdivisionPositions(vertebra[vi]->qe, maxLevel, 0);
        }

        edgeLen=getAverageEdgeLength(vertebra[vi]->qe, minEdge, maxEdge, stdDev);
        if (stdDev>1.0*minSpacing && (vertebra[vi]->center-origCenter).length()>2*call)
        //large difference between shortest and longest edge
        //and center has moved at least 2mm per recursive call (prevents infinite loops)
        {
            if (mainForm.actionInteractive_inflation->isChecked())
            {
                mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
                vertebra[vi]->actor->Delete();
            }
            Cell::kill(vertebra[vi]->qe); //maybe recalculate thresholds first?
            QApplication::processEvents();
            return segmentOneButterfly(previousDist, useMRM, call+1);
            //restart segmentation of this vertebra with much better center estimate
        }

        if (edgeLen>avgSpacing*1.5)
        {
            topologicallySubdivide(vertebra[vi]->qe, ++maxLevel);
            calculateSubdivisionPositions(vertebra[vi]->qe, maxLevel, maxLevel-dependentLevels);
            P.resize((vertebra[vi]->qe)->countVertices(), 3);
            updateWeightMatrix(W, LU, vertebra[vi]->qe, maxLevel, maxLevel-dependentLevels, freeCount);
        }

        if (mainForm.actionInteractive_inflation->isChecked())
        {
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
            vertebra[vi]->actor->Delete();
            vertebra[vi]->actor=poly2actor(vertebra[vi]->getPoly(), mainForm.actionShow_wireframe->isChecked());
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[vi]->actor);
            mainForm.vis->GetRenderWindow()->Render();
            QApplication::processEvents();
        }

        if (extraIter>0)
            extraIter--;
        iter++;
    } while (extraIter>0 || //if extra iterations given OR
        iter<maxIter && vertebra[vi]->radius<previousDist*1.5 //precautionary conditions
        && (vertebra[vi]->radius - oldDist>percentMoreInflated*minSpacing
            || vertebra[vi]->radius - oldDist2>percentMoreInflated*minSpacing)
        );

    vertebra[vi]->mrm=vertebra[vi]->calcRigidTransform(maxLevel-dependentLevels);

    if (mainForm.actionInteractive_inflation->isChecked())
    {
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
        vertebra[vi]->actor->Delete();
    }
    
    vertebra[vi]->actor=poly2actor(vertebra[vi]->getPoly(), mainForm.actionShow_wireframe->isChecked());
    if (mainForm.actionShow_mask_overlays->isChecked())
    {
        vertebra[vi]->getMask();
        vertebra[vi]->addMask();
    }
    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[vi]->actor);
    mainForm.vis->GetRenderWindow()->Render();
    QApplication::processEvents();
    return vertebra[vi]->radius;
}

void MainLogic::segmentationLoop(int direction, const vec3& lastEndplateVector, double sizeGrowth)
{
    if (vi > 29 || vertebra.size()>29)
		return; //precaution
    Vertebra *last1, *last2, *last3=0, *thisVB=new Vertebra(*this);
    if (direction>0)
    {
        last1=vertebra[vi];
        last2=vertebra[vi-1];
        if (vi>1)
            last3=vertebra[vi-2];
        vi++;
        vertebra.push_back(thisVB);
    }
    else
    {
        last1=vertebra[0];
        last2=vertebra[1];
        if (vertebra.size()>2)
            last3=vertebra[2];
        vi=0;
        vertebra.insert(vertebra.begin(), thisVB);
    }

    vec3 ccVec, ccn, axis;

    ccVec=last1->center - last2->center;
    if (last3==0) //linear fit + sizeGrowth
        vertebra[vi]->center=last1->center+lastEndplateVector*ccVec.length()*(1.0+sizeGrowth);
    else //TODO: superlinear polynomial fit using all available vertebrae (+ sizeGrowth?)
    {        
        vertebra[vi]->center=last1->center+lastEndplateVector*ccVec.length()*(1.0+sizeGrowth);
    }

    time=QTime::currentTime();

    thisVB->labelIndex=last1->labelIndex+direction;
    thisVB->centerIndex=vec3ToIndex(thisVB->center, current);
    thisVB->highThreshold=last1->highThreshold;
    thisVB->lowThreshold=last1->lowThreshold;
    thisVB->avgValue=last1->avgValue;
    thisVB->sigmaValue=last1->sigmaValue;
    calcDistanceField();
    mainForm.statusbar->showMessage(QString("Segmenting vertebra ")+Vertebra::labels[vertebra[vi]->labelIndex]);
    QApplication::processEvents();
    
    segmentOneButterfly(last1->radius, false, 1);

    f<<time.toString("hh:mm:ss.zzz\t").toStdString()<<setw(5)<<time.msecsTo(QTime::currentTime())
        <<"\tSegmenting vertebra "<<Vertebra::labels[vertebra[vi]->labelIndex]<<endl;

    if (withinEpsilonEnvironment(last1->radius, thisVB->radius, 0.35) &&
        thisVB->volume/last1->volume>0.5 && thisVB->volume/last1->volume<2 &&
        withinEpsilonEnvironment(vertebra1->radius, thisVB->radius, 0.8) &&
        thisVB->volume/vertebra1->volume>0.08 && thisVB->volume/vertebra1->volume<12.5 &&
        withinEpsilonEnvironment(vertebra1->volume/vertebra1->surface, thisVB->volume/thisVB->surface, 0.50) &&
        (thisVB->center-last1->center).length()>thisVB->radius /*crude overlap check*/)
    {
        clipWrite(vertebra[vi]->qe, current, (filename_base+"_"+Vertebra::labels[vertebra[vi]->labelIndex]+".obj").c_str());
        mainForm.statusbar->showMessage("Recalculating thresholds");
        
        recalcThresholds(thisVB->qe, thisVB->center, thisVB->radius);
        mainForm.statusbar->showMessage("Estimating axis");
        axis=thisVB->calcAxis(lastEndplateVector);
        thisVB->axis=direction*axis*thisVB->radius;
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(thisVB->axis2vtk(false));
        ccn=thisVB->center-last1->center;
        sizeGrowth=(ccn.length()-ccVec.length())/ccn.length(); //expect the same growth in this direction
        segmentationLoop(direction, axis, sizeGrowth); //recursion
    }
    else
    {
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
        vertebra[vi]->actor->Delete();
        if (direction>0)
            vertebra.pop_back();
        else
            vertebra.erase(vertebra.begin());
    }
    mainForm.vis->GetRenderWindow()->Render();
    QApplication::processEvents();
}

//segment first and second vertebra, then call segmentation loop
void MainLogic::segmentVertebrae(vec3& originalCenter, bool segmentOtherVertebrae)
{
    InternalImageType::SpacingType sp = current->GetSpacing();
    vec3 lastAdjecantVertebraVector, firstAdjecantVertebraVector;
    vec3 center, lastCenter, ccVec, axis;

    mainForm.statusbar->showMessage(QString("Segmenting initial vertebra ")+Vertebra::labels[vertebra1->labelIndex]);
    QTime time=QTime::currentTime();
    vertebra[vi]->centerIndex=vec3ToIndex(originalCenter, current);
    calcDistanceField();
    vertebra[vi]->center=originalCenter;
    //first vertebra, axis is just a guess - and it does not influence result that much
    vertebra[vi]->axis=vec3(0,0,1);
    segmentOneButterfly(dist2D, false, 1);
    
    mainForm.statusbar->showMessage("Recalculating thresholds");
    recalcThresholds(vertebra[vi]->qe, vertebra[vi]->center, vertebra[vi]->radius);

    f<<time.toString("hh:mm:ss.zzz\t").toStdString()
        <<setw(5)<<time.msecsTo(QTime::currentTime())<<"\tInitial vertebra "
        <<Vertebra::labels[vertebra1->labelIndex]<<endl;

    clipWrite(vertebra[vi]->qe, current, (filename_base+"_"+Vertebra::labels[vertebra1->labelIndex]+".obj").c_str());
    
    if (!segmentOtherVertebrae)
        return;

    mainForm.statusbar->showMessage("Estimating axis");
    axis=vertebra[vi]->calcAxis();
    
    //segment first adjacent vertebra
    vi++;
    vertebra.push_back(new Vertebra(*this));
    vertebra[vi]->labelIndex=vertebra[vi-1]->labelIndex/**/+1;
    vertebra[vi]->highThreshold=vertebra[vi-1]->highThreshold;
    vertebra[vi]->lowThreshold=vertebra[vi-1]->lowThreshold;
    vertebra[vi]->avgValue=vertebra[vi-1]->avgValue;
    vertebra[vi]->sigmaValue=vertebra[vi-1]->sigmaValue;
    const int numberOfFactors=4;
    double distanceFactors[numberOfFactors]={2.0,1.8,2.2,1.6};
    
    bool found=false;
    for (int i=0; i<numberOfFactors; i++)
    {
        vec3 p=axis;
        p.normalize();
        p*=dist2D*distanceFactors[i];
        p+=originalCenter; //candidate position of adjecant vertebra
        time=QTime::currentTime();
        float useless;
        if (!pos2val(current, p, useless))
            continue;
        vertebra[vi]->centerIndex=vec3ToIndex(p, current);
        calcDistanceField();
        vertebra[vi]->center=p;
        vertebra[vi]->axis=axis;
        mainForm.statusbar->showMessage(QString("Segmenting first adjacent vertebra ")
            +Vertebra::labels[vertebra[vi]->labelIndex]);
        segmentOneButterfly(vertebra[vi-1]->radius, false, 1);
        f<<time.toString("hh:mm:ss.zzz\t").toStdString()
            <<setw(5)<<time.msecsTo(QTime::currentTime())
            <<"\tAttempt to locate adjacent vertebra "
            <<Vertebra::labels[vertebra[vi]->labelIndex]<<endl;
        if (withinEpsilonEnvironment(vertebra[vi-1]->radius, vertebra[vi]->radius, 0.35) &&
            vertebra[vi]->volume/vertebra1->volume>0.5 &&
            vertebra[vi]->volume/vertebra1->volume<2 &&
            //within 20% of the first one by distance (cubed for volume)
            withinEpsilonEnvironment(vertebra1->volume/vertebra1->surface,
                                        vertebra[vi]->volume/vertebra[vi]->surface, 0.3) )
        {
            firstAdjecantVertebraVector=(vertebra[vi]->center-originalCenter);
            firstAdjecantVertebraVector.normalize();
            lastAdjecantVertebraVector=firstAdjecantVertebraVector;
            center=vertebra[vi]->center; lastCenter=originalCenter;
            found=true;
            break;
        }
        else
        {
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
            vertebra[vi]->actor->Delete();
            vertebra[vi]->actor=0;
        }
    }

    if (!found)
    {
        if (vertebra[vi]->actor)
        {
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(vertebra[vi]->actor);
            vertebra[vi]->actor->Delete();
        }
        vertebra.pop_back();
        return; //we only managed to segment the first vertebra
    }

    vertebra1->axis=axis*vertebra[vi]->radius;
    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra1->axis2vtk(false));
    clipWrite(vertebra[vi]->qe, current, (filename_base+"_"+Vertebra::labels[vertebra[vi]->labelIndex]+".obj").c_str());

    mainForm.statusbar->showMessage("Recalculating thresholds");
    recalcThresholds(vertebra[vi]->qe, vertebra[vi]->center, vertebra[vi]->radius);
    
    mainForm.statusbar->showMessage("Estimating axis");
    axis=vertebra[vi]->calcAxis(lastAdjecantVertebraVector);

    vertebra[vi]->axis=axis*vertebra[vi]->radius;    
    mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[vi]->axis2vtk(false));

    segmentationLoop(+1, axis, 0/*initial guess*/);

    double sizeGrowth=pow((double)(vertebra[vertebra.size()-1]->center-vertebra[vertebra.size()-2]->center).length()
                     /(vertebra[1]->center-vertebra[0]->center).length(), 1.0/(vertebra.size()-1))-1;
    //find step of geometric progression

    //now go in the opposite direction (from the first vertebra)
    segmentationLoop(-1, -firstAdjecantVertebraVector, -sizeGrowth);
}

void MainLogic::clearMasks()
{
    masks=VisualizingImageType::New();
    masks->CopyInformation(visualizing);
    masks->SetRegions(masks->GetLargestPossibleRegion());
    masks->Allocate();
    masks->FillBuffer(0);
    mainForm.updateOverlay();
}

void MainLogic::postSegmentation()
{
    QApplication::processEvents();  
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Starting ICP"<<endl;
    for (int i=0; i<vertebra.size(); i++)
    {
        mainForm.statusbar->showMessage("Determining 1st pass vertebral body orientations ("+
            QString::number(i+1)+"/"+QString::number(vertebra.size())+")");
        vtkSmartPointer<vtkMatrix4x4> mrs=vertebra[i]->calcICP(true); //needs scale removal
        float scale=sqrt(mrs->Element[0][0]*mrs->Element[0][0]+
                         mrs->Element[0][1]*mrs->Element[0][1]+
                         mrs->Element[0][2]*mrs->Element[0][2]);
        vtkSmartPointer<vtkMatrix4x4> is=vtkSmartPointer<vtkMatrix4x4>::New();
        for (int e=0; e<3; e++)
            is->Element[e][e]/=scale;
        //vtkMatrix4x4::Multiply4x4(mrs,is,vertebra[i]->mrm);
        vertebra[i]->mrm->DeepCopy(mrs);
        
        //determine axis
        Real p[4]={0,0,1,0};
        vertebra[i]->mrm->MultiplyPoint(p, p);
        vertebra[i]->axis=vec3(p);
        vertebra[i]->axis.normalize();

        if (mainForm.actionSave_binary_masks->isChecked())
        {
            mainForm.statusbar->showMessage(QString("Saving ")+Vertebra::labels[vertebra[i]->labelIndex]+" binary mask");
            //if (!mainForm.actionShow_mask_overlays->isChecked())
                vertebra[i]->getMask();
            writeImage(vertebra[i]->mask, filename_base+"_"+Vertebra::labels[vertebra[i]->labelIndex]+".mha");
        }
    }

    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Finished ICP"<<endl;

    if (vertebra.size()>1) //=1 => error or dummy test
    {
        //fit positions and orientations to polyline. initial guess: a0=x0, b0=y0; other ai,bi =0
        Eigen::VectorXd x,y;
        polynomialFit(4, vertebra,x,y, false);

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

//calculate distance field in region of interest around centerIndex
void MainLogic::calcDistanceField()
{
    InternalImageType::SpacingType sp = current->GetSpacing();
    typedef itk::RegionOfInterestImageFilter<InternalImageType, InternalImageType> RoiType;
    RoiType::Pointer roi=RoiType::New();
    InternalImageType::RegionType region;
    region.SetIndex(0, max<double>(0.0, vertebra[vi]->centerIndex[0]-dist2D*2/sp[0]));
    region.SetIndex(1, max<double>(0.0, vertebra[vi]->centerIndex[1]-dist2D*2/sp[1]));
    region.SetIndex(2, max<double>(0.0, vertebra[vi]->centerIndex[2]-dist2D*2.5/sp[2]));
    region.SetSize(0, min<unsigned long long>(dist2D*4/sp[0], current->GetLargestPossibleRegion().GetSize(0)-region.GetIndex(0)));
    region.SetSize(1, min<unsigned long long>(dist2D*4/sp[1], current->GetLargestPossibleRegion().GetSize(1)-region.GetIndex(1)));
    region.SetSize(2, min<unsigned long long>(dist2D*5/sp[2], current->GetLargestPossibleRegion().GetSize(2)-region.GetIndex(2)));
    roi->SetRegionOfInterest(region);
    roi->SetInput(current);
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(roi->GetOutput(), filename_base+"_roi"+QString::number(vertebra.size()).toStdString()+".mha");

    typedef itk::BinaryThresholdImageFilter<InternalImageType, InternalImageType> BinarizorType;
    BinarizorType::Pointer binTh=BinarizorType::New();
    binTh->SetInput(roi->GetOutput());
    binTh->SetUpperThreshold(vertebra[vi]->highThreshold);
    binTh->SetLowerThreshold(vertebra[vi]->lowThreshold);
    binTh->SetInsideValue(1.0f);
    binTh->Update();

    mainForm.statusbar->showMessage("Filling holes in the regional binary image");
    typedef itk::BinaryBallStructuringElement<float, 3> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(2);
    structuringElement.CreateStructuringElement();
    typedef itk::BinaryClosingByReconstructionImageFilter<InternalImageType, StructuringElementType> BinFillType;
    BinFillType::Pointer binFill=BinFillType::New();
    binFill->SetKernel(structuringElement);
    binFill->SetForegroundValue(1.0f);
    binFill->SetInput(binTh->GetOutput());
    binFill->Update();

    typedef itk::SignedMaurerDistanceMapImageFilter<InternalImageType, InternalImageType> DistanceMapType;
    mainForm.statusbar->showMessage("Calculating distance field from thresholded image...");
    DistanceMapType::Pointer dm2=DistanceMapType::New();
    dm2->SetInput(binFill->GetOutput());
    dm2->SetUseImageSpacing(true);
    dm2->SquaredDistanceOff();
    dm2->SetInsideIsPositive(true);
    dm2->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(dm2->GetOutput(), filename_base+"_dfTh"+QString::number(vertebra.size()).toStdString()+".mha");

    typedef itk::UnaryFunctorImageFilter<InternalImageType,InternalImageType,dfEdgeClassifier> dfClassifier;
    dfClassifier::Pointer dfCl=dfClassifier::New();
    dfCl->SetInput(dm2->GetOutput());
    dfCl->Update();
    vertebra[vi]->cl_dfTh=dfCl->GetOutput();

    typedef itk::UnaryFunctorImageFilter<InternalImageType,InternalImageType, valueClassifier> valClassifier;
    valClassifier::Pointer valCl=valClassifier::New();
    valCl->SetInput(roi->GetOutput());
    valCl->GetFunctor().avgValue=vertebra[vi]->avgValue;
    valCl->GetFunctor().stdDev=vertebra[vi]->sigmaValue;
    valCl->Update();
    vertebra[vi]->cl_value=valCl->GetOutput();

    std::vector<InternalImageType::Pointer> CombinedCl;
    CombinedCl.push_back(vertebra[vi]->cl_dfTh);
    CombinedCl.push_back(cl_dfCanny);
    CombinedCl.push_back(cl_dfMult);
    CombinedCl.push_back(cl_LH2);
    CombinedCl.push_back(binFill->GetOutput());
    vertebra[vi]->clCombined=combineClassifiers(CombinedCl, mean);

    if (mainForm.actionSave_debug_images->isChecked())
    {
        writeImage(vertebra[vi]->cl_dfTh, filename_base+"_dfThCl"+QString::number(vertebra.size()).toStdString()+".mha");
        writeImage(vertebra[vi]->cl_value, filename_base+"_valCl"+QString::number(vertebra.size()).toStdString()+".mha");
        writeImage(vertebra[vi]->clCombined, filename_base+"_clCombined"+QString::number(vertebra.size()).toStdString()+".mha");
        writeImage(binFill->GetOutput(), filename_base+"_filled"+QString::number(vertebra.size()).toStdString()+".mha");
        writeImage(binTh->GetOutput(), filename_base+"_thresh"+QString::number(vertebra.size()).toStdString()+".mha");
    }
}

QPointF centerOfRegion(const QPainterPath &path)
{
    double xc=0, yc=0;
	unsigned long long count=0;
	QRectF lh(path.boundingRect());
	for (int i=lh.y(); i<lh.y()+lh.height(); i++)
		for (int k=lh.x(); k<lh.x()+lh.width() ; k++)
			if (path.contains(QPointF(k,i)))
			{
				xc+=k;
				yc+=i;
				count++;
			}

	if (count>0)
        return QPointF(xc/count, yc/count);
    else
        return QPointF(0,0);
}

void MainLogic::vertebraCentersPicked(int sliceNumber, int labelIndex)
{
	const std::vector<QPointF> & vertebraCenters = mainForm.painter->getVertebraCenters();
    vi=0;
    vertebra.clear();
    
    mainForm.on_actionSave_initialization_triggered(volume_filename.c_str());
    f.open((filename_base+"_times.csv").c_str());
    f<<"StartTime	Duration|ms	Description\n";

	for (std::vector<QPointF>::const_iterator it = vertebraCenters.begin(), end = vertebraCenters.end(); it != end; ++it)
	{
		QPointF neighbourVec;
				
		if (vertebraCenters.size() > 1)
		{
			//look for the nearest neighbour vertebra
			
			double dist = DBL_MAX, currentDist;
					
			for (std::vector<QPointF>::const_iterator compIt = vertebraCenters.begin(), compEnd = vertebraCenters.end(); compIt != compEnd; ++compIt)
			{
				if (*it != *compIt)
				{
					QPointF distVec = *compIt - *it;
					currentDist = sqrt( distVec.x()*distVec.x() + distVec.y()*distVec.y() );

					if (currentDist < dist)
					{
						dist = currentDist;
						neighbourVec = distVec;
					}					
				}
			}
		}
		else
		{
			neighbourVec = QPointF(0.0, -30.0);			
		}
		
		double distRatio	= 1.0 / 3.0, //ratio of half vertebral height to the distance between two neighboured vertebral centers
			   widthRatio	= 1.3,		 //estimated ration between width and height of a vertebra
			   safetyFactor = 0.66;		 //safety scale factor to ensure that we stay inside the vertebral body

		//estimate vertebra rectangle by taking the distance into account
		neighbourVec *= distRatio;
		QPointF orthoVec = QPointF(neighbourVec.y(), -neighbourVec.x()) * widthRatio;

		neighbourVec *= safetyFactor;
		orthoVec	 *= safetyFactor;
		
		QPointF p1 = *it + neighbourVec + orthoVec;
		QPointF p2 = *it - neighbourVec + orthoVec;
		QPointF p3 = *it - neighbourVec - orthoVec;
		QPointF p4 = *it + neighbourVec - orthoVec;

		//clamp the rectangle points to the image borders
		int w = mainForm.painter->get_slice(sliceNumber).width(),
			h = mainForm.painter->get_slice(sliceNumber).height();
		
		p1.setX(p1.x() < 0  ? 0		  : p1.x()); p1.setY(p1.y() < 0  ? 0	   : p1.y());
		p1.setX(p1.x() >= w ? w - 0.1 : p1.x()); p1.setY(p1.y() >= h ? h - 0.1 : p1.y());
		p2.setX(p2.x() < 0  ? 0		  : p2.x()); p2.setY(p2.y() < 0  ? 0	   : p2.y());
		p2.setX(p2.x() >= w ? w - 0.1 : p2.x()); p2.setY(p2.y() >= h ? h - 0.1 : p2.y());
		p3.setX(p3.x() < 0  ? 0		  : p3.x()); p3.setY(p3.y() < 0  ? 0	   : p3.y());
		p3.setX(p3.x() >= w ? w - 0.1 : p3.x()); p3.setY(p3.y() >= h ? h - 0.1 : p3.y());
		p4.setX(p4.x() < 0  ? 0		  : p4.x()); p4.setY(p4.y() < 0  ? 0	   : p4.y());
		p4.setX(p4.x() >= w ? w - 0.1 : p4.x()); p4.setY(p4.y() >= h ? h - 0.1 : p4.y());

		QPainterPath path(p1);
		path.lineTo(p2);
		path.lineTo(p3);
		path.lineTo(p4);
		path.lineTo(p1);
		
		startSegmentation(*it, path, sliceNumber, false, true, labelIndex++);
        vi++;
	}

    postSegmentation();
    f.close();
}

void MainLogic::threePointsPicked(int sliceNumber, int labelIndex)
{
	QPointF vertebraCenter, disk, cord;
	qreal circleRadius;
	mainForm.painter->getThreePointInitialization(vertebraCenter, disk, cord, circleRadius);

	QPainterPath path(vertebraCenter);
	path.addEllipse(vertebraCenter, circleRadius, circleRadius);

	InternalImageType::IndexType ind;	
	ind[2] = sliceNumber;

	ind[0] = disk.x();
	ind[1] = disk.y();
    InternalImageType::PixelType diskValue = current->GetPixel(ind);
	
	ind[0] = cord.x();
	ind[1] = cord.y();
	InternalImageType::PixelType cordValue = current->GetPixel(ind);

    vi=0;
    vertebra.clear();
	startSegmentation(vertebraCenter, path, sliceNumber, true, true, labelIndex, diskValue, cordValue);
}

void MainLogic::singlePointPicked(int sliceNumber, int labelIndex)
{
    QPointF center=mainForm.painter->getSinglePoint();
    QPainterPath path;
	path.addEllipse(center, PointRadius, PointRadius);

    vi=0;
    vertebra.clear();
    startSegmentation(center, path, sliceNumber, true, true, labelIndex);
}

void MainLogic::approximateBoundaryDrawn(int sliceNumber, int labelIndex)
{
    vi=0;
    vertebra.clear();
    startSegmentation(centerOfRegion(*mainForm.painter->getPath()), *mainForm.painter->getPath(), sliceNumber, true, true, labelIndex);
}

void MainLogic::startSegmentation(const QPointF & center2D, QPainterPath & path, int sliceNumber, bool segmentOtherVertebrae, 
                bool calculateThresholds, int labelIndex, InternalImageType::PixelType diskValue, InternalImageType::PixelType cordValue)
{
    useShape=false;

    vertebra1=new Vertebra(*this);
    vertebra.push_back(vertebra1);
    vertebra1->labelIndex=labelIndex;

	InternalImageType::PointType p;
    vertebra[vi]->centerIndex[0]=center2D.x();
    vertebra[vi]->centerIndex[1]=center2D.y();
    vertebra[vi]->centerIndex[2]=sliceNumber;
	current->TransformIndexToPhysicalPoint(vertebra[vi]->centerIndex, p);
	vec3 center(p[0], p[1], p[2]);

    filename_base=tempDir+itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);
    if (mainForm.painter->getInitMethod()!=mainForm.painter->PickVertebraCenters)
        mainForm.on_actionSave_initialization_triggered(volume_filename.c_str());
	InternalImageType::SpacingType sp = current->GetSpacing();
	dist2D=calcAvgDistance(path, center2D); //imprecise
    if (dist2D<=PointRadius+1 && dist2D>=PointRadius-1) //single point picked
        dist2D=35.5/2; //average human V.B. height=27mm width=44mm (~diameter)
    else
        dist2D*=1.3*pow(sp[0]*sp[1], 0.5); //geometric mean
        //1.3 factor compensates for VB geometry (wider than taller)
	if (calculateThresholds)
	{
		mainForm.statusbar->showMessage("Calculating threshold...");
		calcThresholds(path);
	}
    clearMasks(); //also clears user initialization
	//limit thresholds, if desired
	InternalImageType::PixelType middleValue = (vertebra[vi]->highThreshold + vertebra[vi]->lowThreshold) / 2.0;
	InternalImageType::PixelType spacing = 5.0;	
	if (diskValue  != -1)
	{
		if (diskValue	   < middleValue  && diskValue + spacing > vertebra[vi]->lowThreshold)
			vertebra[vi]->lowThreshold  = diskValue + spacing;
		else if (diskValue >= middleValue && diskValue - spacing < vertebra[vi]->highThreshold)
			vertebra[vi]->highThreshold = diskValue - spacing;
	}
	if (cordValue != -1)
	{
		if (cordValue	   < middleValue  && cordValue + spacing > vertebra[vi]->lowThreshold)
			vertebra[vi]->lowThreshold  = cordValue + spacing;
		else if (cordValue >= middleValue && cordValue - spacing < vertebra[vi]->highThreshold)
			vertebra[vi]->highThreshold = cordValue - spacing;
	}

	mainForm.statusbar->showMessage("Segmenting vertebrae...");
    if (mainForm.painter->getInitMethod()!=mainForm.painter->PickVertebraCenters)
    {
        f.open((filename_base+"_times.csv").c_str());
        f<<"StartTime	Duration|ms	Description\n";
    }
    //run the algorithm
    segmentVertebrae(center, segmentOtherVertebrae);
    if (mainForm.painter->getInitMethod()!=mainForm.painter->PickVertebraCenters)
        postSegmentation();

    if (mainForm.painter->getInitMethod()!=mainForm.painter->PickVertebraCenters)
        f.close(); //save execution times
    mainForm.vis->GetRenderWindow()->Render();
	mainForm.statusbar->showMessage("Ready");
}