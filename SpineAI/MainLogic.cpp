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
#include "itkBinaryDilateImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkMedianImageFilter.h"

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
#include "vtkSphereSource.h"

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
#include <QSettings>
#include <QtConcurrentRun>
#include <QtConcurrentMap>
#include <QFuture>

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/operations.hpp"
//#define OPENCV_GPU_SKIP_INCLUDE 1 //resolves Cannot open 'opencv2/opencv_modules.hpp'
//#include "opencv2/gpu/gpu.hpp"

#ifdef _WINDOWS
#include <windows.h>
#else
#include <unistd.h>
#define Sleep(x) usleep((x)*1000)
#endif

using namespace std;
using namespace cv;

MainLogic::MainLogic(MainWindow& mainWindow)
:mainForm(mainWindow)
{
    connect(&mainForm, SIGNAL(volumeOpened(std::string)), SLOT(volumeOpen(std::string)));
    connect(&mainForm, SIGNAL(dicomOpened(std::string)), SLOT(dicomOpen(std::string)));
	connect(mainForm.painter, SIGNAL(vertebraCentersPicked(int)), SLOT(vertebraCentersPicked(int)));
	connect(mainForm.painter, SIGNAL(clipWithBoundingGeometry(bool)), SLOT(clipWithBoundingGeometry(bool)));
	connect(mainForm.painter, SIGNAL(changeBGRadius(double)), SLOT(changeBGRadius(double)));
    connect(mainForm.painter, SIGNAL(centersChanged()), SLOT(recalcCenterline()));

    current=0;
    visualizing=0;
    masks=0;
    fullAuto=false;
    timer=new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(renderNeeded()));
    updateOverlay=false;
    timer->start(20);

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

//parameters
std::string MainLogic::tempDir;
//detection filtering ratios:
double outlierRadiusFactor=0.5;
double firstLastThreshold=2.5;
double smallLargeFactor=1.5;
double underDistance=0.75; //maybe decrease to 0.5?
double overDistance=1.5;
//image filtering radii:
double bilateralDomainSigma=10;
double bilateralRangeSigma=0.03;
double structureTensorSigma=1.0;
double cannyVariance=2;
//image filtering thresholds:
double lhEps=0.01;
double multBinaryThreshold=0.1;
double cannyLowerThreshold=0.05;
double cannyUpperThreshold=0.15;
//feature weights:
double clf[5]={0.2, 0.2, 0.2, 0.2, 0.2 }; //last one =1-sum(4)
//constraint importance factors:
double percentMoreInflated=0.05; //0.02=>2%
double smoothFactor=0.3; //0.1=10%
double cifEdges=3; //constraintImportanceFactor
double sizeGoalForceKappa=0.1;
//disease thresholds:
double crushedVertebraThreshold=0.2;
double significantSpondylolisthesis=0.25;

void readParameters()
{
    QSettings params((MainLogic::tempDir+"params.ini").c_str(), QSettings::IniFormat);
    #define readParam(p) p=params.value(#p, p).toDouble()

    //detection filtering ratios:
    readParam(outlierRadiusFactor);
    readParam(firstLastThreshold);
    readParam(smallLargeFactor);
    readParam(underDistance);
    readParam(overDistance);
    //image filtering radii:
    readParam(bilateralDomainSigma);
    readParam(bilateralRangeSigma);
    readParam(structureTensorSigma);
    readParam(cannyVariance);
    //image filtering thresholds:
    readParam(lhEps);
    readParam(multBinaryThreshold);
    readParam(cannyLowerThreshold);
    readParam(cannyUpperThreshold);
    //feature weights:
    double clfSum=0;
    for (int i=0; i<4; i++)
    {
        clf[i]=params.value(QString("classifier")+char('0'+i)+"weight", clf[i]).toDouble();
        clfSum+=clf[i];
    }
    clf[4]=1-clfSum;
    //constraint importance factors:
    readParam(percentMoreInflated);
    readParam(smoothFactor);
    readParam(cifEdges);
    readParam(sizeGoalForceKappa);
    //disease thresholds:
    readParam(crushedVertebraThreshold);
    readParam(significantSpondylolisthesis);
}

float weightedMean(std::vector<float> p)
{
    float sum=0;
    for (unsigned i=0; i<p.size(); i++)
        sum+=p[i]*clf[i];
    assert(sum<=1);
    return sum; //already ensured that it is <1
}

extern double round(double r);

void MainLogic::renderNeeded()
{
    if (render)
    {
        threadLock.lock();
        mainForm.vis->GetRenderWindow()->Render();
	    QApplication::processEvents();
        render=false;
        wc.wakeAll();
        threadLock.unlock();
    }
    if (updateOverlay)
    {
        updateOverlay=false;
        mainForm.updateOverlay();
    }
}

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

        if (directory[directory.length()-1]=='/' || directory[directory.length()-1]=='\\')
            directory=directory.substr(0, directory.length()-1);
        volume_filename=directory;
        filename_base=tempDir+itksys::SystemTools::GetFilenameName(volume_filename);
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
        filename_base=tempDir+itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);
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

template<class iT>
typename iT::Pointer readImage(std::string fileName)
{
    typedef  itk::ImageFileReader< iT >  ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( fileName.c_str() );
    reader->Update();
    return reader->GetOutput();
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

//apply median denoising, calculate maxValue
void MainLogic::afterOpen(const itk::MetaDataDictionary &metaData)
{
    mainForm.setWindowTitle(QString("Spine Analyzer - ")+QString::fromStdString(volume_filename));
    QApplication::processEvents();

    masks=0;
    mainForm.painter->centerline=QPainterPath();

    mainForm.painter->clearInitialization();

    typedef itk::MinimumMaximumImageCalculator< InternalImageType >  CalculatorType;
    CalculatorType::Pointer calculatorI = CalculatorType::New();
    calculatorI->SetImage( current );
    calculatorI->Compute();
    maxValue=maxValueOriginal = calculatorI->GetMaximum();
    float minValue=calculatorI->GetMinimum();
    if (minValue<0 || minValue>1)
    {
        mainForm.statusbar->showMessage("Intensity shifting"); QApplication::processEvents();
        //shift image intensities, so minimum ==0
        typedef itk::ShiftScaleImageFilter<InternalImageType,InternalImageType> shiftType;
        shiftType::Pointer shifter=shiftType::New();
        shifter->SetInput(current);
        shifter->SetShift(-minValue);
        shifter->Update();
        current=shifter->GetOutput();
        maxValue-=minValue;
        maxValueOriginal-=minValue;
    }
    //mainForm.updateVisualization();
    //mainForm.vis->GetRenderWindow()->Render();
    //detectCenters(); //detect centers using original image (in sagittal orientation)

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
        mainForm.statusbar->showMessage("Window intensity scaling"); QApplication::processEvents();
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
    
    mainForm.statusbar->showMessage("Median filtering"); QApplication::processEvents();
    typedef itk::MedianImageFilter < InternalImageType, InternalImageType> FilterType;
	FilterType::Pointer median = FilterType::New();
    median->SetInput(current);
    median->Update();
    current=median->GetOutput();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(current, tempDir+"median.mha");

    calculatorI->SetImage( current );
    calculatorI->ComputeMaximum();
    maxValue = calculatorI->GetMaximum();

    mainForm.updateVisualization();
    mainForm.vis->GetRenderWindow()->Render();
    detectCenters();
    autoInitSegmentation();
}

void MainLogic::calculateFeatures()
{
    readParameters();
    typedef itk::MinimumMaximumImageCalculator< InternalImageType >  CalculatorType;

    {//unnamed block to reduce peak memory usage
        CalculatorType::Pointer calculatorI = CalculatorType::New();


        mainForm.statusbar->showMessage("Bilateral filtering"); QApplication::processEvents();
        typedef itk::FastBilateralImageFilter<InternalImageType, InternalImageType> BilateralType;
        BilateralType::Pointer bilateral=BilateralType::New();
        bilateral->SetInput(current);
        bilateral->SetDomainSigma(bilateralDomainSigma);
        bilateral->SetRangeSigma(bilateralRangeSigma*maxValue);
        //FilterProgress::Pointer progress=FilterProgress::New(); //progress bar not needed for fast bilateral
        //bilateral->AddObserver(itk::ProgressEvent(),  progress);
        bilateral->Update();
        current=bilateral->GetOutput();
        current->DisconnectPipeline();

        calculatorI->SetImage( current );
        calculatorI->ComputeMaximum();
        maxValue = calculatorI->GetMaximum();
        if (mainForm.actionSave_debug_images->isChecked())
            writeImage(current, tempDir+"bilateral.mha");

        mainForm.updateVisualization();
    
        int i=1;
        while (pow(2.0f,i)-1<maxValue)
            i++;
        maxValue2=pow(2.0f,i)-1;

        mainForm.statusbar->showMessage("Calculating LH values"); QApplication::processEvents();
        lImage=InternalImageType::New();
        hImage=InternalImageType::New();
        calcLHvalues(current, lhEps*maxValue, lImage, hImage);
        if (mainForm.actionSave_LH_images_and_histogram->isChecked())
            saveLHfiles(i);
    } //release intermediate memory

    InternalImageType::Pointer surfel;

    {//unnamed block to reduce peak memory usage
        mainForm.statusbar->showMessage("Calculating structure tensors"); QApplication::processEvents();
        typedef itk::itkStructureTensorImageFilter<InternalImageType> TensorFilterType;
        TensorFilterType::Pointer tensorFilter=TensorFilterType::New();
        tensorFilter->SetInput(current);
        tensorFilter->SetSigma(structureTensorSigma);
        tensorFilter->Update();
	    
        mainForm.statusbar->showMessage("Conducting eigen-analysis"); QApplication::processEvents();
        typedef itk::Image<TensorFilterType::CVector, 3> CovariantImageType;
        typedef itk::SymmetricEigenAnalysisImageFilter<TensorFilterType::OutputImageType,
            CovariantImageType> EigenFilterType;
        EigenFilterType::Pointer ef=EigenFilterType::New();
        ef->SetInput(tensorFilter->GetOutput());
        ef->SetDimension(3);
        ef->OrderEigenValuesBy(EigenFilterType::FunctorType::OrderByMagnitude);
        ef->Update();

        typedef itk::UnaryFunctorImageFilter<CovariantImageType, InternalImageType,
            EigenvaluesToSurfel> Eigen2TypeFilter;
        Eigen2TypeFilter::Pointer e2t=Eigen2TypeFilter::New();
        e2t->SetInput(ef->GetOutput());
        e2t->Update();
        surfel=e2t->GetOutput();
        surfel->DisconnectPipeline();
        if (mainForm.actionSave_debug_images->isChecked())
            writeImage(surfel, tempDir+"eigen2surfel.mha");
    }//release intermediate memory

    mainForm.statusbar->showMessage("Calculating gradient"); QApplication::processEvents();
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
    mult->SetInput1(surfel);
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
    binMult->SetLowerThreshold(multBinaryThreshold*calculatorMult->GetMaximum());
    binMult->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(binMult->GetOutput(), tempDir+"multBin.mha");

    typedef itk::ShiftScaleImageFilter < InternalImageType, InternalImageType > RescaleImageFilterType;
    RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
    rescale->SetInput( mult->GetOutput() );
    rescale->SetScale(1/calculatorMult->GetMaximum());
    rescale->Update();

    mainForm.statusbar->showMessage("Canny edge detection"); QApplication::processEvents();
    typedef itk::CustomCanny<InternalImageType, InternalImageType> cannyType;
    cannyType::Pointer cannyFilter=cannyType::New();
    cannyFilter->SetVariance(cannyVariance);
    cannyFilter->SetUpperThreshold(cannyUpperThreshold);
    cannyFilter->SetLowerThreshold(cannyLowerThreshold);
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

    mainForm.statusbar->showMessage("Calculating distance field from edges"); QApplication::processEvents();
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

    mainForm.statusbar->showMessage("Executing classifiers"); QApplication::processEvents();
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

    //typedef itk::BinaryFunctorImageFilter<InternalImageType, InternalImageType,
    //    InternalImageType, lhClassifierFunctor> lhClassifier;
    //lhClassifier::Pointer lhCl=lhClassifier::New();
    //lhCl->GetFunctor().maxVal=maxValue;
    //lhCl->SetInput1(lImage);
    //lhCl->SetInput2(hImage);
    //lhCl->Update();
    //cl_LH1=lhCl->GetOutput();
    //if (mainForm.actionSave_debug_images->isChecked())
    //    writeImage(cl_LH1, tempDir+"clLH1.mha");

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

    if (mainForm.actionSave_debug_images->isChecked())
    {
        mainForm.statusbar->showMessage("Combining per-dataset classifiers");
        QApplication::processEvents();
        std::vector<InternalImageType::Pointer> classifiers(4);
        classifiers[0]=cl_LH1;
        classifiers[1]=cl_LH2;
        classifiers[2]=cl_dfCanny;
        classifiers[3]=cl_dfMult;
        InternalImageType::Pointer clCombinedDS=combineClassifiers(classifiers, mean);
        writeImage(clCombinedDS, tempDir+"clCombinedDS.mha");
    }
    
    mainForm.statusbar->showMessage("Calculating per vertebra features (in parallel)");
    QApplication::processEvents();
    vector<QFuture<void> > r;
	for (int i=0; i<vertebra.size(); i++)
	{       
        QFuture<void> future = QtConcurrent::run<void>(this, &MainLogic::calculatePerVertebraFeatures, i);
        r.push_back(future);
	}
    for (int i=0; i<vertebra.size(); i++)
        r[i].waitForFinished();

    mainForm.statusbar->showMessage("Ready");
    QApplication::processEvents();
}

void MainLogic::calculatePerVertebraFeatures(int i) //vertebra index
{
    VisualizingImageType::SpacingType sp = current->GetSpacing();
    QPainterPath path;
    path.addEllipse(vertebra[i]->centerIndex[0], vertebra[i]->centerIndex[1],
        circleRadius/current->GetSpacing()[0], circleRadius/current->GetSpacing()[1]);
    calcThresholds(path, vertebra[i]); //uses current
    //results in highThreshold, lowThreshold, avgValue, sigmaValue
    //calculate thresholds for each vertebra individually,
    //then smooth them out using medians to a degree (paramter)?
        
    InternalImageType::RegionType vRoI;
    vRoI.SetIndex(0, max<double>(0.0, vertebra[i]->centerIndex[0]-vertebra[i]->radius*1.5/sp[0]));
    vRoI.SetIndex(1, max<double>(0.0, vertebra[i]->centerIndex[1]-vertebra[i]->radius*1.5/sp[1]));
    vRoI.SetIndex(2, max<double>(0.0, vertebra[i]->centerIndex[2]-vertebra[i]->radius*2/sp[2]));
    vRoI.SetSize(0, min<unsigned long long>(vertebra[i]->radius*3/sp[0], current->GetLargestPossibleRegion().GetSize(0)-vRoI.GetIndex(0)));
    vRoI.SetSize(1, min<unsigned long long>(vertebra[i]->radius*3/sp[1], current->GetLargestPossibleRegion().GetSize(1)-vRoI.GetIndex(1)));
    vRoI.SetSize(2, min<unsigned long long>(vertebra[i]->radius*4/sp[2], current->GetLargestPossibleRegion().GetSize(2)-vRoI.GetIndex(2)));
        
    typedef itk::RegionOfInterestImageFilter<InternalImageType, InternalImageType> RoiType;
    RoiType::Pointer roi=RoiType::New();
    roi->SetRegionOfInterest(vRoI);
    roi->SetInput(current);
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(roi->GetOutput(), filename_base+"_roi"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");

    typedef itk::BinaryThresholdImageFilter<InternalImageType, InternalImageType> BinarizorType;
    BinarizorType::Pointer binTh=BinarizorType::New();
    binTh->SetInput(roi->GetOutput());
    binTh->SetUpperThreshold(vertebra[i]->highThreshold);
    binTh->SetLowerThreshold(vertebra[i]->lowThreshold);
    binTh->SetInsideValue(1.0f);
    binTh->Update();

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

    typedef itk::SignedMaurerDistanceMapImageFilter<InternalImageType, InternalImageType> InternalDistanceMapType;
    InternalDistanceMapType::Pointer dm2=InternalDistanceMapType::New();
    dm2->SetInput(binFill->GetOutput());
    dm2->SetUseImageSpacing(true);
    dm2->SquaredDistanceOff();
    dm2->SetInsideIsPositive(true);
    dm2->Update();
    if (mainForm.actionSave_debug_images->isChecked())
        writeImage(dm2->GetOutput(), filename_base+"_dfTh"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");

    typedef itk::UnaryFunctorImageFilter<InternalImageType,InternalImageType,dfEdgeClassifier> dfClassifier;
    dfClassifier::Pointer dfCl=dfClassifier::New();
    dfCl->SetInput(dm2->GetOutput());
    dfCl->Update();
    vertebra[i]->cl_dfTh=dfCl->GetOutput();

    //typedef itk::UnaryFunctorImageFilter<InternalImageType,InternalImageType, valueClassifier> valClassifier;
    //valClassifier::Pointer valCl=valClassifier::New();
    //valCl->SetInput(roi->GetOutput());
    //valCl->GetFunctor().avgValue=vertebra[i]->avgValue;
    //valCl->GetFunctor().stdDev=vertebra[i]->sigmaValue;
    //valCl->Update();
    //vertebra[i]->cl_value=valCl->GetOutput();
    //if (mainForm.actionSave_debug_images->isChecked())
    //    writeImage(vertebra[i]->cl_value, filename_base+"_valCl"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");

    std::vector<InternalImageType::Pointer> CombinedCl;
    CombinedCl.push_back(vertebra[i]->cl_dfTh);
    CombinedCl.push_back(cl_dfCanny);
    CombinedCl.push_back(cl_dfMult);
    CombinedCl.push_back(cl_LH2);
    CombinedCl.push_back(binFill->GetOutput());
    vertebra[i]->clCombined=combineClassifiers(CombinedCl, weightedMean); //remove later

    if (mainForm.actionSave_debug_images->isChecked())
    {
        writeImage(vertebra[i]->cl_dfTh, filename_base+"_dfThCl"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");        
        writeImage(vertebra[i]->clCombined, filename_base+"_clCombined"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");
        writeImage(binFill->GetOutput(), filename_base+"_filled"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");
        writeImage(binTh->GetOutput(), filename_base+"_thresh"+string(Vertebra::labels[vertebra[i]->labelIndex])+".mha");
    }
}

inline double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void MainLogic::inflate(Vertebra *vert, float stepSize, bool useShape, float force)
{
    typedef itk::LinearInterpolateImageFunction<InternalImageType,float> InterpolatorType;
    InterpolatorType::Pointer interp=InterpolatorType::New();
    interp->SetInputImage(vert->clCombined);
    InterpolatorType::ContinuousIndexType ind;

    InternalImageType::PixelType pixelFw, pixelBw;
    vec3 posFw, posBw;
    CellVertexIterator it3(vert->qe);
    Vertex *v;
    while ((v = it3.next()) != 0)
    {
        vec3 n=v->pos-vert->center; //center-vertex vector
        n.normalize();

        vec3 mp;
        if (useShape)//use model to steer inflation direction
        {
            mp=vert->mrm*Vertebra::model->vertices[v->realIndex()]->pos;
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

        itk::Point<float, 3> fw(posFw.ptr());
        if (vert->clCombined->TransformPhysicalPointToContinuousIndex<float>(fw, ind)) //in image?
        {
            pixelFw=interp->EvaluateAtContinuousIndex( ind );
            v->pos=posBw; //deflate
            v->wasDeflated=true;
        }
        else
            continue; //next vertex
        if (pixelFw+force>0.5)
        {
            v->pos=posFw; //inflate
            v->wasDeflated=false;
        }
        else
        {
            itk::Point<float, 3> bw(posBw.ptr());
            if (!vert->clCombined->TransformPhysicalPointToContinuousIndex<float>(bw, ind)) //outside?
                continue; //do not change wasDeflated property
            else
                pixelBw=interp->EvaluateAtContinuousIndex( ind );
            if (pixelFw<=pixelBw+force)
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

void MainLogic::calcThresholds(const QPainterPath &path, Vertebra *v)
{
    InternalImageType::IndexType ind;
    double xc=v->centerIndex[0]*0.2,
           yc=v->centerIndex[1]*0.2; //move 20% towards center
    int zc=v->centerIndex[2];
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
    v->avgValue=sum/count;
    
    sum=0; //now calculate standard deviation
    for (pCounter::iterator it=vox.begin(); it!=vox.end();)
    {
        sum+=it->second*(it->first-v->avgValue)*(it->first-v->avgValue);
        it++;
    }
    v->sigmaValue=sqrt(sum/count);

    v->highThreshold=vox.rbegin()->first+noiseUnit; //add 1 noiseUnit to round up
    v->lowThreshold=vox.begin()->first; //already rounded down
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

//inherit so I can add save/load without modifying library code
template <class MatrixType, int QRPreconditioner=2>
class MySVD: public Eigen::JacobiSVD<MatrixType, QRPreconditioner>
{
public:
    typedef typename MatrixType::Index Index;
    /** Saves the decomposition into a file */
    void save(const char *filename) const
    {
        ofstream f(filename, ios::binary);
        f.write((char *)&this->m_rows, sizeof(this->m_rows));
        f.write((char *)&this->m_cols, sizeof(this->m_cols));
        f.write((char *)&this->m_nonzeroSingularValues, sizeof(this->m_nonzeroSingularValues));
        f.write((char *)&this->m_singularValues[0], sizeof(Scalar)*this->m_diagSize);
        f.write((char *)&this->m_matrixU(0,0), sizeof(Scalar)*this->m_rows*this->m_cols);
        f.write((char *)&this->m_matrixV(0,0), sizeof(Scalar)*this->m_cols*this->m_cols);
        f.close();
    }

    /** Loads the decomposition from a file */
    void load(const char *filename)
    {
        this->m_isAllocated=false;
        ifstream f(filename, ios::binary);
        f.read((char *)&this->m_rows, sizeof(this->m_rows));
        f.read((char *)&this->m_cols, sizeof(this->m_cols));
        allocate(this->m_rows, this->m_cols, Eigen::ComputeThinU | Eigen::ComputeThinV);
        f.read((char *)&this->m_nonzeroSingularValues, sizeof(this->m_nonzeroSingularValues));
        f.read((char *)&this->m_singularValues[0], sizeof(Scalar)*this->m_diagSize);
        f.read((char *)&this->m_matrixU(0,0), sizeof(Scalar)*this->m_rows*this->m_cols);
        f.read((char *)&this->m_matrixV(0,0), sizeof(Scalar)*this->m_cols*this->m_cols);
        if (f.bad())
            throw std::exception();
        f.close();
        this->m_isInitialized=true;
    }
protected: //cannot access ancestor's private method, so copy it here
    void allocate(Index rows, Index cols, unsigned int computationOptions)
    {
        eigen_assert(rows >= 0 && cols >= 0);

        if (this->m_isAllocated &&
            rows == this->m_rows &&
            cols == this->m_cols &&
            computationOptions == this->m_computationOptions)
        {
            return;
        }

        this->m_rows = rows;
        this->m_cols = cols;
        this->m_isInitialized = false;
        this->m_isAllocated = true;
        this->m_computationOptions = computationOptions;
        this->m_computeFullU = (computationOptions & Eigen::ComputeFullU) != 0;
        this->m_computeThinU = (computationOptions & Eigen::ComputeThinU) != 0;
        this->m_computeFullV = (computationOptions & Eigen::ComputeFullV) != 0;
        this->m_computeThinV = (computationOptions & Eigen::ComputeThinV) != 0;
        eigen_assert(!(this->m_computeFullU && this->m_computeThinU) && "JacobiSVD: you can't ask for both full and thin U");
        eigen_assert(!(this->m_computeFullV && this->m_computeThinV) && "JacobiSVD: you can't ask for both full and thin V");
        eigen_assert(EIGEN_IMPLIES(this->m_computeThinU || this->m_computeThinV, MatrixType::ColsAtCompileTime==Eigen::Dynamic) &&
                    "JacobiSVD: thin U and V are only available when your matrix has a dynamic number of columns.");
        if (QRPreconditioner == Eigen::FullPivHouseholderQRPreconditioner)
        {
            eigen_assert(!(this->m_computeThinU || this->m_computeThinV) &&
                    "JacobiSVD: can't compute thin U or thin V with the FullPivHouseholderQR preconditioner. "
                    "Use the ColPivHouseholderQR preconditioner instead.");
        }
        this->m_diagSize = (std::min)(this->m_rows, this->m_cols);
        this->m_singularValues.resize(this->m_diagSize);
        this->m_matrixU.resize(this->m_rows, this->m_computeFullU ? this->m_rows
                                : this->m_computeThinU ? this->m_diagSize
                                : 0);
        this->m_matrixV.resize(this->m_cols, this->m_computeFullV ? this->m_cols
                                : this->m_computeThinV ? this->m_diagSize
                                : 0);
        this->m_workMatrix.resize(this->m_diagSize, this->m_diagSize);
    }
};

QMutex mutex; //do not calculate SVD in parallel
//create a matrix of wieghts which contains influence of free vertices on subdivided vertices
void calcSVD(MySVD<Eigen::MatrixXd>& svd,
    Cell *qe, unsigned maxLevel, unsigned controlLevel, unsigned& freeCount)
{
    QString levelInfo="_"+QString::number(maxLevel)+"_"+QString::number(controlLevel);
    levelInfo+="_"+QString::number(cifEdges,'f',16);//+"_"+QString::number(cifSmoothness,'f',16);
    QDir cache("./cache_SVD");
    if (!cache.exists())
        cache.mkpath(".");
    QString filename="./cache_SVD/"+QString(vertebraMeshFilename)+levelInfo+".svd";
    QFile fsvd(filename);

    mutex.lock();
    if (fsvd.exists())
    {
        mutex.unlock();
        svd.load(filename.toStdString().c_str());
        freeCount=svd.cols();
    }
    else //keep other threads blocked until this one has calculated SVD
    {
        Eigen::MatrixXd weights=createWeightMatrix(qe, maxLevel, controlLevel, freeCount, cifEdges); //, cifSmoothness);
        svd.compute(weights, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.save(filename.toStdString().c_str());
        mutex.unlock();
    }    
}

//returns average distance from center to surface
double MainLogic::segmentOneButterfly(Vertebra *v, unsigned call)
{
    //mainForm.statusbar->showMessage(QString("Segmenting vertebra ")+Vertebra::labels[v->labelIndex]);
    InternalImageType::SpacingType sp = current->GetSpacing();
    double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);
    double minSpacing=min(min(sp[0], sp[1]), sp[2]);
    calcVertexValence(v->qe, false);
    CellVertexIterator it3(v->qe);
    Vertex *vtx;
    copyPos2SubdivPos(v->qe);

    unsigned maxLevel=2, freeCount, controlLevel=0;
    for (unsigned i=1; i<=maxLevel; i++)
        topologicallySubdivide(v->qe, i);

    Eigen::Matrix<double,Eigen::Dynamic,3> P;
    MySVD<Eigen::MatrixXd> svd;
    calcSVD(svd, v->qe, maxLevel, controlLevel, freeCount);
    P.resize(svd.rows(), 3);

    calculateSubdivisionPositions(v->qe, maxLevel, controlLevel);
    copySubdivPos2Pos(v->qe);

    double origRadius=v->radius; //try to reach this size
    v->calcCurrentRadius();
    v->calcNewCenter();
    double oldVolume, oldDist2, oldDist=v->radius;
    vec3 origCenter=v->center;
    float edgeLen, stdDev, minEdge, maxEdge, force;
    edgeLen=getAverageEdgeLength(v->qe, minEdge, maxEdge, stdDev);
    int iter=0, extraIter=5; //extraIter=minimum number of iterations
    int maxIter=250;
    do
    {
        force=sizeGoalForceKappa*tanh((origRadius-v->radius)/origRadius)*(oldDist+minSpacing-v->radius)/minSpacing;
        inflate(v, minSpacing, false, force); //constant step size

        {
            CellVertexIterator cvi0(v->qe);
            while ((vtx = cvi0.next())!=0)
			{
                P.row(vtx->realIndex())[0]=vtx->pos[0];
				P.row(vtx->realIndex())[1]=vtx->pos[1];
				P.row(vtx->realIndex())[2]=vtx->pos[2];
			}


            unsigned i=v->qe->countVertices();
            CellVertexIterator cvi(v->qe);
            while ((vtx=cvi.next())!=0)
            {
                //assign3(P.row(i), cifSmoothness*v->valence*v->pos);
                for (int k=0; k<vtx->valence; k++)
				{
					P.row(i+k)[0]=cifEdges*2*vtx->pos[0];
					P.row(i+k)[1]=cifEdges*2*vtx->pos[1];
					P.row(i+k)[2]=cifEdges*2*vtx->pos[2];
				}
                i+=vtx->valence; //i+=v->valence+1;
            }
        }

        oldDist2=oldDist;
        oldDist=v->radius;
        oldVolume=v->volume;
        v->calcCurrentRadius();
        v->calcNewCenter();

        if (smoothFactor!=0.0)
            normalizeHeuristic(v->qe, maxLevel, controlLevel); //implicit smoothing

        {
            //optimal least sqaures normalization
            Eigen::Matrix<double,Eigen::Dynamic,3> N;
            N=svd.solve(P);
            CellVertexIterator cvi(v->qe);
            while ((vtx = cvi.next())!=0)
                if (vtx->level<=controlLevel)
				{
                    vtx->pos[0]=N.row(vtx->realIndex())[0];
                    vtx->pos[1]=N.row(vtx->realIndex())[1];
                    vtx->pos[2]=N.row(vtx->realIndex())[2];
				}
            swapSubdivPos(v->qe);
            calculateSubdivisionPositions(v->qe, maxLevel, controlLevel);
        }

        if (smoothFactor!=0.0)
        {
            //now combine optimal (subdivPos) and heuristic (pos)
            CellVertexIterator cvi0(v->qe);
            while ((vtx = cvi0.next())!=0)
            {
                //part from smooth heuristic, part from jagged optimal
                vtx->pos=vtx->pos*smoothFactor+vtx->subdivPos*(1-smoothFactor);
            }
            copyPos2SubdivPos(v->qe);
        }
        else
            copySubdivPos2Pos(v->qe);
        

        edgeLen=getAverageEdgeLength(v->qe, minEdge, maxEdge, stdDev);
        if (stdDev>1.0*minSpacing && (v->center-origCenter).length()>2*call)
        //large difference between shortest and longest edge
        //and center has moved at least 2mm per recursive call (prevents infinite loops)
        {
            Cell::kill(v->qe); //maybe recalculate thresholds first?
            v->qe=objReadCell(vertebraMeshFilename);
            v->radius=origRadius;
            rotate(v->qe, vec3(0,0,1), v->axis);
            translate(v->qe, v->center);
            return segmentOneButterfly(v, call+1);
            //restart segmentation of this vertebra with much better center estimate
        }

        if (maxLevel<maxRefinementLevel && (edgeLen>avgSpacing*1.5 ||
            v->radius - oldDist<percentMoreInflated*minSpacing))
        {
            topologicallySubdivide(v->qe, ++maxLevel);
            controlLevel=maxLevel/2-1;
            calculateSubdivisionPositions(v->qe, maxLevel, controlLevel);
            calcSVD(svd, v->qe, maxLevel, controlLevel, freeCount);
            P.resize(svd.rows(), 3);
            extraIter+=3;
        }

        if (mainForm.actionInteractive_inflation->isChecked())
        {
            vtkSmartPointer<vtkPolyData> pd=v->getPoly();
            threadLock.lock();
            ((vtkPolyDataMapper *)v->actor->GetMapper())->SetInputData(pd);
            v->actor->GetMapper()->Modified();
            render=true;
            threadLock.unlock();
            wc.wait(&v->renderLock);
        }

        if (extraIter>0)
            extraIter--;
        iter++;
    } while (extraIter>0 || //if extra iterations given OR
        iter<maxIter && v->radius<origRadius*1.5 //precautionary conditions
        && (v->radius - oldDist>percentMoreInflated*minSpacing
            || v->radius - oldDist2>percentMoreInflated*minSpacing)
        );

    v->mrm=v->calcRigidTransform(controlLevel);
   
    vtkSmartPointer<vtkPolyData> pd=v->getPoly();
    threadLock.lock();
    ((vtkPolyDataMapper *)v->actor->GetMapper())->SetInputData(pd);
    v->actor->GetMapper()->Modified();
    render=true;
    threadLock.unlock();

    clipWrite(v->qe, current, (filename_base+"_"+Vertebra::labels[v->labelIndex]+".obj").c_str());
    if (mainForm.actionShow_mask_overlays->isChecked())
        v->addMask();

    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Vertebra "<<Vertebra::labels[v->labelIndex]<<" segmented"<<endl;
    return v->radius;
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
    if (mainForm.actionDiagnosis_view->isChecked())
        mainForm.on_actionClear_polygonal_data_triggered();
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Starting ICP and binary mask saving"<<endl;
    //mask calculation takes the most time in this for loop (ICP/mask)
    for (int i=0; i<vertebra.size(); i++)
    {
        mainForm.statusbar->showMessage("Determining vertebral body orientations ("+
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
                vertebra[i]->getOriginalMask();
            writeImage(vertebra[i]->mask, filename_base+"_"+Vertebra::labels[vertebra[i]->labelIndex]+".mha", true);
        }
        //debug
        //Cell *t=Vertebra::model->deepCopy();
        //applyMatrix(t, mrs);
        //objWriteCell(t, (QString::fromStdString(filename_base)+"_"+QString::number(i)+"_tm.obj").toStdString().c_str());
        //mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(poly2actor(qe2vtk(t),true,vec3(0,1,0)));
        //vtkMatrix4x4::Multiply4x4(mrs,is,vertebra[i]->mrm);
        //Cell::kill(t);

        if (mainForm.actionDiagnosis_view->isChecked())
        {
            vtkSphereSource *sphere = vtkSphereSource::New(); 
            sphere->SetRadius(3.0); 
            sphere->SetThetaResolution(18); 
            sphere->SetPhiResolution(18);

            vtkPolyDataMapper *map = vtkPolyDataMapper::New(); 
            map->SetInputConnection(sphere->GetOutputPort());

            vtkActor *aSphere = vtkActor::New(); 
            aSphere->SetMapper(map); 
            aSphere->GetProperty()->SetColor(1,0,0);
            aSphere->SetUserMatrix(vertebra[i]->mrm);
            mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(aSphere);
        }
    }

    //ofstream cvr("D:/centerVolumeRadiusIndex.txt"); //for bugTester
    //for (int i=0; i<vertebra.size(); i++)
    //{
    //    //mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[i]->axis2vtk(true));
    //    cvr<<vertebra[i]->center[0]<<' '<<vertebra[i]->center[1]<<' '<<vertebra[i]->center[2];
    //    cvr<<' '<<vertebra[i]->volume<<' '<<vertebra[i]->radius<<' '<<vertebra[i]->labelIndex<<endl;
    //}
    //cvr.close();
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Finished ICP"<<endl;

    if (vertebra.size()>1) //=1 => error or dummy test
    {
        //fit positions to polyline
        polynomialFit(4, vertebra,x,y, false);
        drawCenterline();

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

//estimate linear fit to points (i,y[i])
void TheilSenEstimator(vector<double> y, double &slope, double &yIntercept)
{
    if (y.empty())
        return;
    if (y.size()==1)
    {
        slope=0;
        yIntercept=y[0];
        return;
    }

	vector<double> m;
	for (int i=0; i<y.size()-1; i++)
		for (int k=i+1; k<y.size(); k++)
			m.push_back((y[k]-y[i])/(k-i)); //slopes
	std::sort(m.begin(), m.end());
	slope=m[m.size()/2]; //median of slopes
		
	m.clear();
	for (int i=0; i<y.size(); i++)
		m.push_back(y[i]-slope*i); //y-intercepts
	std::sort(m.begin(), m.end());
	yIntercept=m[m.size()/2]; //median of y-intercepts
}

double distance3D(vec3 pointA, vec3 pointB)
{
    return std::sqrt((pointA[0]-pointB[0])*(pointA[0]-pointB[0])
               +(pointA[1]-pointB[1])*(pointA[1]-pointB[1])
               +(pointA[2]-pointB[2])*(pointA[2]-pointB[2]));
}

double curvature(Eigen::VectorXd P, double x) //evaluate curvature at given point x
{
    Eigen::VectorXd dP=derive(P);
    double dPe=polyEval(dP, x);
    return dEval(dP, x)/std::pow((1+dPe*dPe), 1.5);
}

void MainLogic::vertebraCentersPicked(int labelIndex)
{
    std::vector<CenterRadiusType> centers = mainForm.painter->centers;
    if (centers.size()<2)
    {
        mainForm.statusbar->showMessage("At least 2 vertebra required");
        return;
    }
    vertebra.clear();
    clearMasks();

    time=QTime::currentTime();
    f.open((filename_base+"_times.csv").c_str());
    f<<"Time    Event\n";
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Started segmentation"<<endl;

    //re-order centers bottom to top
    for (int i=0; i<centers.size()-1; i++)
        for (int k=i+1; k<centers.size(); k++)
            if (centers[i].y<centers[k].y)
                swap(centers[i], centers[k]);

    //fill in missing size
    vector<double> sizes;
    for (int i=0; i<centers.size(); i++)
        if (centers[i].r!=0)
            sizes.push_back(centers[i].r);
    double slope, fS;
    if (sizes.size()<2) //at least 2 points required for linear fitting
    {
        slope=0;
        fS=12; //24mm average lumbar VB height, radius is half that
    } else
        TheilSenEstimator(sizes, slope, fS);
    for (int i=0; i<centers.size(); i++)
        if (centers[i].r==0)
            centers[i].r=slope*i+fS;

    InternalImageType::SpacingType sp = current->GetSpacing();
    double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);
    maxRefinementLevel=round(log(centers[0].r/(1.3*avgSpacing))/log(2.0));
    //1.3 is ratio or radius to average edge lentgh in base mesh (level 0)
    //edge length of the largest VB will be between 1s and 2s, s=avgSpacing
    
    //initialize vertebra array
    for (int i=0; i<centers.size(); i++)
    {
        Vertebra *vertebra1=new Vertebra(*this);
        vertebra.push_back(vertebra1);
        vertebra1->labelIndex=labelIndex++;

	    InternalImageType::PointType p;
        vertebra1->centerIndex[0]=centers[i].x;
        vertebra1->centerIndex[1]=centers[i].y;
        vertebra1->centerIndex[2]=centers[i].z;
	    visualizing->TransformIndexToPhysicalPoint(vertebra1->centerIndex, p);	    
        vertebra1->center=vec3(p[0], p[1], p[2]);
        current->TransformPhysicalPointToIndex(p, vertebra1->centerIndex);
        vertebra1->radius=centers[i].r;
    }

    polynomialFit(3, vertebra,x,y, true); //for initial axis estimates
    drawCenterline();
    for (int i=0; i<centers.size(); i++)
        vertebra[i]->axis=vec3(dEval(x, vertebra[i]->center[2]),
                               dEval(y, vertebra[i]->center[2]), 1);

    //initialize other stuff
    filename_base=tempDir+itksys::SystemTools::GetFilenameWithoutLastExtension(volume_filename);
    mainForm.saveInitialization(Vertebra::labels[vertebra[0]->labelIndex], volume_filename.c_str());
    calculateFeatures(); //needed for segmentation
    f<<QTime::currentTime().toString("hh:mm:ss.zzz\t").toStdString()<<"Features calculated"<<endl;

    for (int i=0; i<vertebra.size(); i++)
    {
        vertebra[i]->qe=objReadCell(vertebraMeshFilename);
        rotate(vertebra[i]->qe, vec3(0,0,1), vertebra[i]->axis);
        translate(vertebra[i]->qe, vertebra[i]->center);
        vertebra[i]->actor=poly2actor(vertebra[i]->getPoly(), mainForm.actionShow_wireframe->isChecked());
        mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(vertebra[i]->actor);
    }
    mainForm.vis->GetRenderWindow()->Render(); //debug
    QApplication::processEvents(); //debug

    vector<QFuture<double> > r;
	for (int i=0; i<vertebra.size(); i++)
	{       
        QFuture<double> future = QtConcurrent::run<double>(this, &MainLogic::segmentOneButterfly, vertebra[i], 1);
        r.push_back(future);
        vertebra[i]->renderLock.lock(); //needs to be locked for QWaitCondition
	}

    mainForm.statusbar->showMessage("Segmenting vertebrae in parallel");    
    vector<bool> ToDo(vertebra.size(), true);
    bool anyToDo=true;
    while (anyToDo)
    {
        anyToDo=false;
        bool sleepNeeded=true;
        for (int i=0; i<vertebra.size(); i++)
            if (ToDo[i])
                if (r[i].isFinished())
                {          
                    vertebra[i]->actor->GetProperty()->SetColor(1,0.8,0); //orange
                    render=true; //request render
                    updateOverlay=true; //request overlay to be updated
                    QString message="Segmenting";
                    for (int i2=0; i2<vertebra.size(); i2++)
                    {
                        if (r[i2].isFinished())
                            continue; //finished
                        else
                            message+=" "+QString(Vertebra::labels[vertebra[i2]->labelIndex]);
                        if (r[i2].isPaused())
                            message+="p";
                    }
                    mainForm.statusbar->showMessage(message);
                    ToDo[i]=false;
                    sleepNeeded=false;
                }
                else
                    anyToDo=true;
        if (sleepNeeded)
            Sleep(10); //10 milliseconds
        QApplication::processEvents();
    }

    QApplication::processEvents();
    postSegmentation(); //diagnosis
    f.close();
    mainForm.statusbar->showMessage("Ready");
    if (fullAuto)
        exit(0); //normal exit
}

double distanceToPoly(vec3 p, Eigen::VectorXd x, Eigen::VectorXd y, vec3& closestPoint)
{
    closestPoint[0]=polyEval(x, p[2]);
    closestPoint[1]=polyEval(y, p[2]);
    closestPoint[2]=p[2];
    vec3 pa=closestPoint-p;
    double dot, dAlen, paLen;
    do
    {
        vec3 dA(dEval(x, closestPoint[2]), dEval(y, closestPoint[2]), 1);
        dot=dA*pa;
        dAlen=dA.length();
        paLen=pa.length();
        closestPoint-=dA*(dot/dAlen);
        pa=closestPoint-p;
    } while (pa.length()<paLen && abs(dot)/(dAlen*paLen)>0.1);
    //process converges and cos>0.1
    return pa.length();
}

double distanceToPoly2D(vec3 p, Eigen::VectorXd x, Eigen::VectorXd y, vec3& closestPoint)
{
    p[0]=0;
    Eigen::VectorXd x0(x.size());
    for (int i=0; i<x0.size(); i++)
        x0[i]=0;
    double dist=distanceToPoly(p,x0,y,closestPoint);
    closestPoint[0]=polyEval(x, closestPoint[2]);
    return dist;
}

double findClosest(vec3 p, vector<vec3> points, int &closest)
{
    closest=0;
    double minDist=distance3D(p, points[0]);
    for (int i=1; i<points.size(); i++)
        if (distance3D(p, points[i])<minDist)
        {
            minDist=distance3D(p, points[i]);
            closest=i;
        }
    return minDist;
}

void MainLogic::removeVertebra(unsigned index)
{
    assert(index<vertebra.size());
    if (index<vertebra.size()-1)
        for (int i=index; i<vertebra.size()-1; i++)
        {
            vertebra[i]=vertebra[i+1];
            mainForm.painter->centers[i]=mainForm.painter->centers[i+1];
        }
    vertebra.pop_back();
    mainForm.painter->centers.pop_back();
}

// class for grouping object candidates, detected by Cascade Classifier, HOG etc.
// instance of the class is to be passed to cv::partition (see cxoperations.hpp)
//class SimilarRects
//{
//public:
//    SimilarRects(double _eps) : eps(_eps) {}
//    inline bool operator()(const Rect& r1, const Rect& r2) const
//    {
//        double delta = eps*(std::min(r1.width, r2.width) + std::min(r1.height, r2.height))*0.5;
//        return std::abs(r1.x - r2.x) <= delta &&
//        std::abs(r1.y - r2.y) <= delta &&
//        std::abs(r1.x + r1.width - r2.x - r2.width) <= delta &&
//        std::abs(r1.y + r1.height - r2.y - r2.height) <= delta;
//    }
//    double eps;
//};

//// This function splits the input sequence or set into one or more equivalence classes and
//// returns the vector of labels - 0-based class indexes for each element.
//// predicate(a,b) returns true if the two sequence elements certainly belong to the same class.
////
//// The algorithm is described in "Introduction to Algorithms"
//// by Cormen, Leiserson and Rivest, the chapter "Data structures for disjoint sets"
//template<typename _Tp, class _EqPredicate> int
//partition( const std::vector<_Tp>& _vec, std::vector<int>& labels,
//           _EqPredicate predicate=_EqPredicate())
//{
//    int i, j, N = (int)_vec.size();
//    const _Tp* vec = &_vec[0];
//
//    const int PARENT=0;
//    const int RANK=1;
//
//    std::vector<int> _nodes(N*2);
//    int (*nodes)[2] = (int(*)[2])&_nodes[0];
//
//    // The first O(N) pass: create N single-vertex trees
//    for(i = 0; i < N; i++)
//    {
//        nodes[i][PARENT]=-1;
//        nodes[i][RANK] = 0;
//    }
//
//    // The main O(N^2) pass: merge connected components
//    for( i = 0; i < N; i++ )
//    {
//        int root = i;
//
//        // find root
//        while( nodes[root][PARENT] >= 0 )
//            root = nodes[root][PARENT];
//
//        for( j = 0; j < N; j++ )
//        {
//            if( i == j || !predicate(vec[i], vec[j]))
//                continue;
//            int root2 = j;
//
//            while( nodes[root2][PARENT] >= 0 )
//                root2 = nodes[root2][PARENT];
//
//            if( root2 != root )
//            {
//                // unite both trees
//                int rank = nodes[root][RANK], rank2 = nodes[root2][RANK];
//                if( rank > rank2 )
//                    nodes[root2][PARENT] = root;
//                else
//                {
//                    nodes[root][PARENT] = root2;
//                    nodes[root2][RANK] += rank == rank2;
//                    root = root2;
//                }
//                CV_Assert( nodes[root][PARENT] < 0 );
//
//                int k = j, parent;
//
//                // compress the path from node2 to root
//                while( (parent = nodes[k][PARENT]) >= 0 )
//                {
//                    nodes[k][PARENT] = root;
//                    k = parent;
//                }
//
//                // compress the path from node to root
//                k = i;
//                while( (parent = nodes[k][PARENT]) >= 0 )
//                {
//                    nodes[k][PARENT] = root;
//                    k = parent;
//                }
//            }
//        }
//    }
//
//    // Final O(N) pass: enumerate classes
//    labels.resize(N);
//    int nclasses = 0;
//
//    for( i = 0; i < N; i++ )
//    {
//        int root = i;
//        while( nodes[root][PARENT] >= 0 )
//            root = nodes[root][PARENT];
//        // re-use the rank as the class label
//        if( nodes[root][RANK] >= 0 )
//            nodes[root][RANK] = ~nclasses++;
//        labels[i] = ~nodes[root][RANK];
//    }
//
//    return nclasses;
//}

void filterRectangles(vector<Rect>& rectList, vector<int>& slices)
{
    const int groupThreshold=4;
    const double eps=0.2;
    vector<int> labels;
    int nclasses = partition(rectList, labels, SimilarRects(0.2));

    vector<Rect> rrects(nclasses);
    vector<int> rweights(nclasses, 0);
    vector<int> sliceMid(nclasses, 0);

    int i, j, nlabels = (int)labels.size();
    for( i = 0; i < nlabels; i++ )
    {
        int cls = labels[i];
        rrects[cls].x += rectList[i].x;
        rrects[cls].y += rectList[i].y;
        rrects[cls].width += rectList[i].width;
        rrects[cls].height += rectList[i].height;
        rweights[cls]++;
        sliceMid[cls]+=slices[i];
    }

    for( i = 0; i < nclasses; i++ )
    {
        Rect r = rrects[i];
        float s = 1.f/rweights[i];
        rrects[i] = Rect(saturate_cast<int>(r.x*s),
             saturate_cast<int>(r.y*s),
             saturate_cast<int>(r.width*s),
             saturate_cast<int>(r.height*s));
        sliceMid[i]=round(sliceMid[i]/rweights[i]);
    }

    rectList.clear();
    slices.clear();

    for( i = 0; i < nclasses; i++ )
    {
        Rect r1 = rrects[i];
        int n1 = rweights[i];
        if( n1 <= groupThreshold )
            continue;
        // filter out small face rectangles inside large rectangles
        for( j = 0; j < nclasses; j++ )
        {
            Rect r2 = rrects[j];

            int dx = saturate_cast<int>( r2.width * eps );
            int dy = saturate_cast<int>( r2.height * eps );

            if( i != j &&
                r1.x >= r2.x - dx &&
                r1.y >= r2.y - dy &&
                r1.x + r1.width <= r2.x + r2.width + dx &&
                r1.y + r1.height <= r2.y + r2.height + dy )
                break;
        }

        if( j == nclasses )
        {
            rectList.push_back(r1);
            slices.push_back(sliceMid[i]);
        }
    }
}

void MainLogic::detectCenters()
{
    mainForm.statusbar->showMessage("Autodetecting vertebra centers");
    CascadeClassifier cascade;
    cascade.load("cascade.xml");
    //gpu::CascadeClassifier_GPU cascadeGPU;
    //cascadeGPU.load("cascade.xml"); //only the old cascade format supported for now

    VisualizingImageType::SizeType size=visualizing->GetLargestPossibleRegion().GetSize();
    InternalImageType::SpacingType sp=visualizing->GetSpacing();
    double xres=sp[0], yres=sp[1];
    unsigned int sbf=std::max(1.0, round(2*sqrt(xres*yres)/sp[2])); //slice binning factor

    vector<int> slices;
    vector<Rect> rects;
    Mat img;

    for (int slice=0; slice<size[2]-sbf; slice+=sbf)
    {       
        if (sbf==1)
        {
            unsigned char *p=visualizing->GetBufferPointer()+slice*visualizing->GetOffsetTable()[2];
            img=Mat(size[1], size[0], CV_8UC1, p); //just take current slice
        }
        else //slice binning increases SNR and reduces number of slices to be processed
        {
            img=Mat::zeros(size[1], size[0], CV_8UC1);
            for (int i=0; i<sbf; i++)
            {
                unsigned char *p=visualizing->GetBufferPointer()+(slice+i)*visualizing->GetOffsetTable()[2];
                img+=(1.0/sbf)*Mat(size[1], size[0], CV_8UC1, p);
            }
        }

        Size scaledImageSize( size[0]*mainForm.painter->mulX, size[1]*mainForm.painter->mulY );
        Mat scaledImage( scaledImageSize, CV_8U);
        resize( img, scaledImage, scaledImageSize);
		//equalizeHist(scaledImage, scaledImage);
        //if (slice>=size[2]/2-sbf && slice<=size[2]/2+sbf) //only middle
        //{
        //    cv::imwrite("D:/img.png",img);
        //    cv::imwrite("D:/imgScaled.png",scaledImage);
        //}

        vector<Rect> candidates; //cpu version
        cascade.detectMultiScale( scaledImage, candidates, 1.2, 1, 0,
            Size(mainForm.painter->mulX*20/xres, mainForm.painter->mulY*20/yres),
            Size(mainForm.painter->mulX*50/xres, mainForm.painter->mulY*50/yres) );
        //size 20-50 mm
        rects.insert(rects.end(), candidates.begin(), candidates.end());
        slices.insert(slices.end(), candidates.size(), slice+sbf/2);

        //gpu::GpuMat objbuf; //gpu version
        //gpu::GpuMat scaledImageGPU(scaledImage);
        //cv::imwrite("D:/imgGPU.png",scaledImageGPU);
        //int numRect=cascadeGPU.detectMultiScale( scaledImageGPU, objbuf, 1.2, 1,
        //    Size(mainForm.painter->mulX*20/xres, mainForm.painter->mulY*20/yres) );
        //Mat obj_host;
        //// download only detected number of rectangles
        //objbuf.colRange(0, numRect).download(obj_host);
        //Rect* gpuRects = obj_host.ptr<Rect>();
        //rects.reserve(rects.size()+numRect);
        //for (int i=0; i<numRect; i++)
        //    rects.insert(rects.end(), gpuRects[i]);
        //slices.insert(slices.end(), numRect, slice);
    }

    mainForm.painter->centers.clear();
    if (rects.size()==0)
        return; //no vertebral bodies detected
    filterRectangles(rects, slices);
    for (int i=0; i<rects.size(); i++)
    {
        rects[i].x/=mainForm.painter->mulX;
        rects[i].y/=mainForm.painter->mulY;
        rects[i].width/=mainForm.painter->mulX;
        rects[i].height/=mainForm.painter->mulY;
    }
    int centerSlice=0;
    for (int i=0; i<rects.size(); i++)
        centerSlice+=slices[i];
    centerSlice=round(centerSlice/rects.size());
    mainForm.painter->sliceSlider->setValue(centerSlice);

    for (int i=0; i<rects.size(); i++)
    {
        CenterRadiusType cr;
        cr.x=rects[i].x+rects[i].width/2;
        cr.y=rects[i].y+rects[i].height/2;
        cr.z=slices[i];
        cr.r=sqrt((double)rects[i].width*xres*rects[i].height*yres)/2;
        mainForm.painter->centers.push_back(cr);
    }

    //re-order centers bottom to top
    for (int i=0; i<mainForm.painter->centers.size()-1; i++)
        for (int k=i+1; k<mainForm.painter->centers.size(); k++)
            if (mainForm.painter->centers[i].y<mainForm.painter->centers[k].y)
                swap(mainForm.painter->centers[i], mainForm.painter->centers[k]);
    mainForm.painter->setVertebraCenters(mainForm.painter->centers);
    mainForm.statusbar->showMessage("Ready");
    QApplication::processEvents();
}

void MainLogic::autoInitSegmentation()
{
    if (mainForm.painter->centers.size()<4)
        return; //3rd order polynomial will be junk
    mainForm.statusbar->showMessage("Auto-initializing segmentation");

    centersToVertebrae();

    //filter obvious misdetections - far away candidates
    while (true) //exit with break
    {
        polynomialFit(3, vertebra,x,y, true);
        drawCenterline();
        int index;
        double maxDist=-1;
        for (int i=0; i<vertebra.size(); i++)
        {
            vec3 unused;
            double dist=distanceToPoly2D(vertebra[i]->center, x, y, unused);
            if ( dist>maxDist)
            {
                maxDist=dist;
                index=i;
            }
        }        
        if (maxDist<=vertebra[index]->radius*outlierRadiusFactor)
            break; //no erroneus candidates remain
        if (vertebra.size()<4) //4 points minimum for poly3 fitting
            break; //prevent crashing
        removeVertebra(index); //else remove this vertebra candidate
    }

    //remove small/large detections
    vector<double> radius;
    for (int i=0; i<vertebra.size(); i++)
        radius.push_back(vertebra[i]->radius);

    double radiusK, radiusN;
	//distance formula estimate: radius(Vi)=radiusK*i+radiusN
	TheilSenEstimator(radius, radiusK, radiusN);
    
    for (int i=0; i<vertebra.size(); i++)
        if (vertebra[i]->radius<(radiusK*i+radiusN)/smallLargeFactor ||
            vertebra[i]->radius>smallLargeFactor*(radiusK*i+radiusN))
        {
            removeVertebra(i);
            i--;
            drawCenterline();
        }

    radius.clear();
    for (int i=0; i<vertebra.size(); i++)
        radius.push_back(vertebra[i]->radius);
	TheilSenEstimator(radius, radiusK, radiusN); //recalculate it with more correct data

    double lowerZ=vertebra[0]->center[2];
    double upperZ=vertebra[vertebra.size()-1]->center[2];
    //remove first and/or last if needed
	double avgDist=distance3D(vertebra[0]->center, vertebra[vertebra.size()-1]->center)/vertebra.size();
    if (distance3D(vertebra[0]->center, vertebra[1]->center)>firstLastThreshold*avgDist)
        removeVertebra(0);
    if (distance3D(vertebra[vertebra.size()-1]->center,
                   vertebra[vertebra.size()-2]->center) > firstLastThreshold*avgDist)
        vertebra.pop_back(); //remove last vertebra candidate

    polynomialFit(3, vertebra,x,y, true); //fresh poly fitting
    drawCenterline();

    vec3 tp; //temporary point
    distanceToPoly(vertebra[0]->center, x, y, tp);
    lowerZ=tp[2];//lowerZ=vertebra[0]->center[2];
    distanceToPoly(vertebra[vertebra.size()-1]->center, x, y, tp);
    upperZ=tp[2];//upperZ=vertebra[vertebra.size()-1]->center[2];

    ////show polyline for debugging purposes
    //vtkActor *polyLine=centerline(x,y,lowerZ,upperZ, vec3(0,0,1)); //blue
    //mainForm.vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(polyLine);
    //mainForm.vis->GetRenderWindow()->Render();
    //QApplication::processEvents();

	vector<double> d; //distances between centers
	for (int i=1; i<vertebra.size(); i++)
		d.push_back(distance3D(vertebra[i-1]->center, vertebra[i]->center));

	double slope, yIntercept;
	//distance formula estimate: d(Vi,Vi+1)=slope*i+yIntercept
    TheilSenEstimator(d, slope, yIntercept);
	
	//check for extraneous detections (small distances)
	//first check for paired small distances
	for (int i=0; i<d.size()-1; i++)
		if (d[i]<underDistance*(slope*i+yIntercept) &&
			d[i+1]<underDistance*(slope*(i+1)+yIntercept) )
		{
			d[i]+=d[i+1];
			d.erase(d.begin()+i+1);
			removeVertebra(i+1);
			TheilSenEstimator(d, slope, yIntercept); //recalculate fitting
            drawCenterline();
		}
	//check first and last vertebrae for small distance
	if (d[0]<yIntercept*underDistance)
	{
		d.erase(d.begin());
		removeVertebra(0);
		TheilSenEstimator(d, slope, yIntercept); //recalculate fitting
        drawCenterline();
	}
	if (d[d.size()-1]<underDistance*(slope*(d.size()-1)+yIntercept))
	{
		d.pop_back();
		vertebra.pop_back();
		TheilSenEstimator(d, slope, yIntercept); //recalculate fitting
        drawCenterline();
	}
	polynomialFit(3, vertebra,x,y, true); //fresh poly fitting
    drawCenterline();

	//fill gaps
	for (int i=d.size()-1; i>=0; i--)
		if (d[i]>overDistance*(slope*i+yIntercept))
		{
			double estDist=slope*i+yIntercept;
			int n=round(d[i]/estDist);
            vec3 A=vertebra[i]->center;
			vec3 B=vertebra[i+1]->center;

			for (int k=n-1; k>0; k--)
			{
				Vertebra *v=new Vertebra(*this);
				vec3 c=(A*(n-k)+B*k)/n;
				distanceToPoly(c, x, y, c);
				v->center=vec3(polyEval(x, c[2]), polyEval(y, c[2]), c[2]);
				v->radius=radiusK*i+radiusN;
				vertebra.insert(vertebra.begin()+i+1, v);
                d.insert(d.begin()+i+1,d[i]/n);
			}
            d[i]/=n;

			TheilSenEstimator(d, slope, yIntercept); //recalculate fitting
            vertebraeToCenters();

            radius.clear();
            for (int i=0; i<vertebra.size(); i++)
                radius.push_back(vertebra[i]->radius);
	        TheilSenEstimator(radius, radiusK, radiusN);
		}

    for (int i=0; i<vertebra.size(); i++)
        vertebra[i]->labelIndex=i+4; //start from S1
    mainForm.painter->bottomLabel->setCurrentIndex(4); //S1

	mainForm.statusbar->showMessage("Ready");

    if (fullAuto)
        vertebraCentersPicked(vertebra[0]->labelIndex);
}

void MainLogic::vertebraeToCenters()
{
    mainForm.painter->centers.clear();
    for (int i=0; i<vertebra.size(); i++)
    {
        InternalImageType::PointType p;
        p[0]=vertebra[i]->center[0];
        p[1]=vertebra[i]->center[1];
        p[2]=vertebra[i]->center[2];
        visualizing->TransformPhysicalPointToIndex(p, vertebra[i]->centerIndex);
        mainForm.painter->centers.push_back(CenterRadiusType(vertebra[i]->centerIndex[0], vertebra[i]->centerIndex[1], vertebra[i]->centerIndex[2], vertebra[i]->radius));
    }
    mainForm.painter->setVertebraCenters(mainForm.painter->centers);
    drawCenterline();
}

void MainLogic::drawCenterline()
{
    if (vertebra.size()<2 || x.size()<3 || y.size()<3)
        return; //nothing to do

    double fromZ=vertebra[0]->center[2];
    double toZ=vertebra[vertebra.size()-1]->center[2];
    if (fromZ>toZ)
        std::swap(fromZ,toZ);
    double step=1.0, z=fromZ;

    itk::ContinuousIndex<double, 3> ind;        
    VisualizingImageType::PointType p;
    mainForm.painter->centerline=QPainterPath(); //clear path
    
    while (z<toZ)
    {
        p[0]=polyEval(x,z);
        p[1]=polyEval(y,z);
        p[2]=z;
        visualizing->TransformPhysicalPointToContinuousIndex(p, ind);

        if (mainForm.painter->centerline.elementCount()==0)
            mainForm.painter->centerline.moveTo(
                ind[0]*mainForm.painter->mulX,
                ind[1]*mainForm.painter->mulY);
        else
            mainForm.painter->centerline.lineTo(
                ind[0]*mainForm.painter->mulX,
                ind[1]*mainForm.painter->mulY);

        z+=step;
    }
    mainForm.painter->paintSlice();
    QApplication::processEvents();
}

void MainLogic::centersToVertebrae()
{
    vertebra.clear();
    for (int i=0; i<mainForm.painter->centers.size(); i++)
    {        
        Vertebra *vert=new Vertebra(*this);
        vert->radius=mainForm.painter->centers[i].r;
        vert->centerIndex[0]=mainForm.painter->centers[i].x;
        vert->centerIndex[1]=mainForm.painter->centers[i].y;
        vert->centerIndex[2]=mainForm.painter->centers[i].z;
        InternalImageType::PointType p;
        visualizing->TransformIndexToPhysicalPoint(vert->centerIndex, p);
        vert->center=vec3(p[0], p[1], p[2]);
        vertebra.push_back(vert);        
    }
}

void MainLogic::recalcCenterline()
{
    centersToVertebrae();
    polynomialFit(3, vertebra,x,y, true);
    drawCenterline();
}

//Clip using the bounding geometry
void MainLogic::clipWithBoundingGeometry(bool triggered)
{	
    if(vertebra.size()!=0)
		mainForm.updateVisualization();
	else
		mainForm.statusbar->showMessage("No vertebra data available",3000);
}

//Change the current radius of the bounding geometry.
void MainLogic::changeBGRadius(double radius)
{
    if (mainForm.painter->clip_checkBox->isChecked())
        clipWithBoundingGeometry(true);
}
