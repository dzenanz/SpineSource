#include "mainWindow.h"
#include "MainLogic.h"
#include <QFileDialog>
#include <QSettings>
#include <QGLWidget>
#include <qmessagebox.h>

#include "itkDatImageIOFactory.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkImageImport.h"

#include <vtkRenderWindow.h>
#include <vtkCamera.h>
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
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkOpenGLRenderer.h>

#include <fstream>

MainWindow::MainWindow(QWidget *parent)
:QMainWindow(parent)
{
	QGLFormat f;
    f.setAlpha( true );
    QGLFormat::setDefaultFormat( f );
	setupUi(this); //instantiates QGLWidget

	QActionGroup *subdiv = new QActionGroup(menuOptions);
	subdiv->addAction(actionSubdivisionDisabled);
	subdiv->addAction(actionSubdivisionHeuristic);
	subdiv->addAction(actionSubdivisionOptimal);

    vtkOpenGLRenderer *renderer = vtkOpenGLRenderer::New();
    vis->GetRenderWindow()->AddRenderer(renderer);
    vis->GetRenderWindow()->SetAlphaBitPlanes(1);
    renderer->SetBackground(0.5f,0.5f,1.0f);
    renderer->Delete();
    volume=0;
	on_actionReset_TF_triggered();

	QSettings settings;
    settings.beginGroup("MainWindow");
    resize(settings.value("size", QSize(800, 600)).toSize());
    move(settings.value("pos", QPoint(100, 100)).toPoint());
	actionInteractive_inflation->setChecked(settings.value("interactiveInflation", false).toBool());
	actionSave_binary_masks->setChecked(settings.value("saveMasks", false).toBool());
	actionShow_wireframe->setChecked(settings.value("wireframe", true).toBool());
	actionSubdivisionDisabled->setChecked(settings.value("subdivisionDisabled", false).toBool());
	actionSubdivisionHeuristic->setChecked(settings.value("subdivisionHeuristic", false).toBool());
	actionSubdivisionOptimal->setChecked(settings.value("subdivisionOptimal", false).toBool());
    actionSave_LH_images_and_histogram->setChecked(settings.value("saveLHimages", true).toBool());
    actionMaintain_vertex_surface_density->setChecked(settings.value("maintainVertexDensity", false).toBool());
    settings.endGroup();

	itk::ObjectFactoryBase::RegisterFactory( itk::DatImageIOFactory::New() );
	logic=new MainLogic(*this);
}
MainWindow::~MainWindow()
{
	delete logic;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
	QSettings settings;
    settings.beginGroup("MainWindow");
    settings.setValue("size", size());
    settings.setValue("pos", pos());
	settings.setValue("interactiveInflation", actionInteractive_inflation->isChecked());
	settings.setValue("saveMasks", actionSave_binary_masks->isChecked());
	settings.setValue("wireframe", actionShow_wireframe->isChecked());
	settings.setValue("subdivisionDisabled", actionSubdivisionDisabled->isChecked());
	settings.setValue("subdivisionHeuristic", actionSubdivisionHeuristic->isChecked());
	settings.setValue("subdivisionOptimal", actionSubdivisionOptimal->isChecked());
    settings.setValue("saveLHimages", actionSave_LH_images_and_histogram->isChecked());
    settings.setValue("maintainVertexDensity", actionMaintain_vertex_surface_density->isChecked());
    settings.endGroup();
}

void MainWindow::on_actionReset_TF_triggered()
{
	opacity=0.05;
	vThLow=90;
	vThHigh=240;
	constructTF();
}

void MainWindow::on_actionTF_shift_Up_triggered()
{
	if (vThHigh<240)
		vThHigh+=15;
	else if (vThLow<225)
		vThLow+=15;
	constructTF();
}

void MainWindow::on_actionTF_shift_Down_triggered()
{
	if (vThLow>15)
		vThLow-=15;
	else if (vThHigh>30)
		vThHigh-=15;
	constructTF();
}

void MainWindow::on_actionOpacityPlus_triggered()
{
	if (opacity<0.8)
		opacity*=1.25;
	constructTF();
}

void MainWindow::on_actionOpacityMinus_triggered()
{
	if (opacity>0.01)
		opacity*=0.8;
	constructTF();
}

void MainWindow::constructTF()
{
	transFunc=vtkSmartPointer<vtkVolumeProperty>::New();
	transFunc->SetInterpolationTypeToLinear();
    vtkColorTransferFunction *colorFun = vtkColorTransferFunction::New();
    vtkPiecewiseFunction *opacityFun = vtkPiecewiseFunction::New();
    transFunc->SetColor( colorFun );
    transFunc->SetScalarOpacity( opacityFun );
    
    colorFun->AddRGBPoint( vThLow, .3, .9, .3 );
    colorFun->AddRGBPoint( vThHigh, .9, .3, .3 );
    opacityFun->AddPoint(0, 0);
    opacityFun->AddPoint(vThLow, .01 );
	opacityFun->AddPoint(vThHigh, opacity);
	opacityFun->AddPoint(255, opacity );
    transFunc->ShadeOn();
    
    if (volume) //volume is loaded and set up
    {
        volume->SetProperty( transFunc );
        volume->GetProperty()->GetRGBTransferFunction()->Modified();
        volume->GetProperty()->GetScalarOpacity()->Modified();
        vis->GetRenderWindow()->Render();
    }
}

void MainWindow::on_actionQuick_help_triggered()
{
    QMessageBox::information(this, "Keyboard shortcuts","While 3D view widget has the focus:\n\
        s: switch polygonal models to surface view\n\
        w: switch polygonal models to wireframe view\n\
		PageUp/PageDown: adjust transfer function\n\
		+/-: adjust transparency");
}

void MainWindow::on_actionLoad_polygonal_mesh_triggered()
{
    std::string s = QFileDialog::getOpenFileName(this, tr("Open polygonal surface mesh"), "", 
        //tr("Polygonal surfaces (*.obj *.osg *.stl *.sta *.iv *.wrl);;All files (*.*)")).toStdString();
        tr("Wavefront surface (*.obj);;All files (*.*)")).toStdString();
	if (s!="")
    {
        statusbar->showMessage("Loading polygonal model "+QString::fromStdString(s));
        vtkSmartPointer<vtkOBJReader> reader = //switch depending on the extension
        vtkSmartPointer<vtkOBJReader>::New();
        reader->SetFileName(s.c_str());
        reader->Update();

        vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());

        vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
        vis->GetRenderWindow()->Render();
        vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->ResetCamera();
        QApplication::processEvents();
        vis->GetRenderWindow()->Render();
        QApplication::processEvents();
        statusbar->showMessage("Ready");
    }
}

void MainWindow::on_actionOpen_triggered() //"Open" clicked on menu
{
	std::string s = QFileDialog::getOpenFileName(this, tr("Open Volume"), "",
        tr("3D images (*.dat *.mha *.mhd *.nrrd *.nii *.hdr *.dcm);;All files (*.*)")).toStdString();
    if (s!="")
        emit volumeOpened(s);
}

QPainterPath loadPath(std::ifstream& f)
{
	float x,y;
	f>>x>>y;
	QPainterPath p(QPointF(x,y));
	while (f>>x>>y)
		p.lineTo(QPointF(x,y));
	f.close();
	return p;
}

void MainWindow::openInitialization(const char * initFilename)
{
    char buffer[2048];
    std::ifstream f(initFilename);
    f.getline(buffer, 2048); //"file"
    f.getline(buffer, 2048);
    if (strcmp(logic->volume_filename.c_str(),buffer))
        emit volumeOpened(buffer);
    QApplication::processEvents();
    
    f.getline(buffer, 2048);
    QString slice(buffer);
    this->painter->sliceSlider->setValue(slice.toInt());

    f.getline(buffer, 2048);
    std::string label(buffer);
    int vertebraIndex=0;
    while (Vertebra::labels[vertebraIndex]!=label)
        vertebraIndex++;
    painter->bottomLabel->setCurrentIndex(vertebraIndex);

    f.getline(buffer, 2048); //"vertebra outline:"
	painter->path=loadPath(f);
	painter->update();
    f.close();
	logic->majorSliceUpdate(painter->sliceSlider->value(), vertebraIndex);
}

void MainWindow::on_actionOpen_initialization_triggered()
{
	std::string s = QFileDialog::getOpenFileName(this, tr("Open Initialization"), "", tr("Spine analyzer initializations (*.init)")).toStdString();
	if (s!="")
	    openInitialization(s.c_str());        
}

void MainWindow::on_actionSave_slices_triggered()
{
	std::string s = QFileDialog::getSaveFileName(this, tr("Save slices"), "", tr("2D images (*.png *.bmp *.ppm *.tiff *.xpm *.xbm *.jpg)")).toStdString();
	if (s!="")
	{
		statusbar->showMessage("Saving all the slices of the dataset...");
		this->painter->saveSlices(s);
		statusbar->showMessage("Ready");
	}
}

void MainWindow::on_actionScreenshot_triggered()
{
	std::string s = QFileDialog::getSaveFileName(this, tr("Save screenshot"), "", tr("2D images (*.png *.bmp *.jpg)")).toStdString();
	if (s!="")
	{
        statusbar->showMessage("Saving hi resolution screenshot");        
        vtkWindowToImageFilter *windowGrabber=vtkWindowToImageFilter::New();
	    windowGrabber->SetInput(vis->GetRenderWindow());
        windowGrabber->SetInputBufferTypeToRGBA();
        windowGrabber->SetMagnification(4); //hi-res screenshot
        
        vtkPNGWriter *writer=vtkPNGWriter::New();
	    writer->SetInputConnection(windowGrabber->GetOutputPort());
        writer->SetFileName(s.c_str());
	    writer->Write();
        
        windowGrabber->Delete();
        writer->Delete();
        vis->GetRenderWindow()->Render();
        statusbar->showMessage("Ready");
	}
}

void MainWindow::defaultCameraPos()
{
	VisualizingImageType::PointType p;
	unsigned xsize=logic->visualizing->GetLargestPossibleRegion().GetSize(0),
			 ysize=logic->visualizing->GetLargestPossibleRegion().GetSize(1),
			 zsize=logic->visualizing->GetLargestPossibleRegion().GetSize(2);
	VisualizingImageType::IndexType ind;

	ind[0]=xsize/2; ind[2]=zsize/2;
	ind[1]=0;
	logic->visualizing->TransformIndexToPhysicalPoint(ind, p);
    vec3 top(p[0], p[1], p[2]);

	ind[0]=xsize-1; ind[1]=ysize/2; ind[2]=zsize/2;
	logic->visualizing->TransformIndexToPhysicalPoint(ind, p);
    vec3 right(p[0], p[1], p[2]);

	ind[0]=xsize/2; ind[1]=ysize/2; ind[2]=zsize/2;
	logic->visualizing->TransformIndexToPhysicalPoint(ind, p);
    vec3 center(p[0], p[1], p[2]);

	vec3 up(top-center);
    vec3 eye((right-center)^up); //view from left
    eye=center+4*eye/up.length();
    vtkCamera *cam=vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
    cam->SetPosition(eye[0], eye[1], eye[2]);
    cam->SetViewUp(up[0], up[1], up[2]);
    cam->SetFocalPoint(center[0], center[1], center[2]);
}

QVector<QPixmap> extractSlices(VisualizingImageType::Pointer itkImage, QColor color)
{
    itk::Vector<double, 3> spacing = itkImage->GetSpacing();
	VisualizingImageType::SizeType size = itkImage->GetLargestPossibleRegion().GetSize();
	QVector<QPixmap> result(size[2]);
	VisualizingImageType::IndexType ind;
	unsigned char val;
	QRgb rgb;

	for (int z=0; z<size[2]; z++)
	{
		ind[2]=z;
		QImage img(size[0], size[1], QImage::Format_ARGB32);
		for (int y=0; y<size[1]; y++)
		{
			ind[1]=y;
			for (int x=0; x<size[0]; x++)
			{
				ind[0]=x;
				val=itkImage->GetPixel(ind);
				if (color.alpha()==255)
					rgb=QColor(color.red()*val/255.0, color.green()*val/255.0, color.blue()*val/255.0).rgba();
				else
					rgb=QColor(color.red(), color.green(), color.blue(), color.alpha()*val).rgba();
				img.setPixel(x,y,rgb);
			}
		}
        img.setDotsPerMeterX(1000/spacing[0]);
        img.setDotsPerMeterY(1000/spacing[1]);
		result[z]=QPixmap::fromImage(img);
	}
	return result;
}

void MainWindow::updateVisualization()
{
	statusbar->showMessage("Rescaling intensities to 8 bits");
    typedef itk::ShiftScaleImageFilter < InternalImageType, VisualizingImageType > RescaleImageFilterType;
    RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
    rescale->SetInput( logic->current );
    rescale->SetScale(255/logic->maxValue);
    rescale->Update();
    logic->visualizing=rescale->GetOutput();
    
    statusbar->showMessage("Converting itkImage to vtkVolume");
    typedef itk::ImageToVTKImageFilter<VisualizingImageType> itkVtkConverter;
    itkVtkConverter::Pointer conv=itkVtkConverter::New();
    conv->SetInput(logic->visualizing);
    conv->Update();
    vtkSmartPointer<vtkImageData> image=vtkSmartPointer<vtkImageData>::New();
    image->ShallowCopy(conv->GetOutput());
    //shallow copy is vtk's equivalent of disconnect pipeline

    vtkSmartPointer<vtkGPUVolumeRayCastMapper> mapper = 
        vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
    mapper->SetInputData(image);

    vtkRenderer *renderer = vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

    if (volume)
    {
        renderer->RemoveViewProp(volume);
        volume->Delete();
    }
    volume=vtkVolume::New();
	constructTF();
    volume->SetMapper( mapper );

    VisualizingImageType::DirectionType d=logic->visualizing->GetDirection();
    VisualizingImageType::PointType origin=logic->visualizing->GetOrigin();
    vtkSmartPointer<vtkMatrix4x4> mat=vtkSmartPointer<vtkMatrix4x4>::New();
    for (int i=0; i<3; i++)
        for (int k=0; k<3; k++)
            mat->SetElement(i,k, d(i,k));
    for (int i=0; i<3; i++)
        mat->SetElement(i,3, origin[i]);
    vtkSmartPointer<vtkMatrix4x4> t=vtkSmartPointer<vtkMatrix4x4>::New();
    for (int i=0; i<3; i++)
        t->SetElement(i,3, -origin[i]); //invert(Mtrans)
    vtkMatrix4x4::Multiply4x4(mat,t,mat); //keep it all in UserMatrix
    volume->SetUserMatrix(mat);

    //renderer->RemoveAllViewProps();
    renderer->AddVolume( volume );
    defaultCameraPos();    
    vis->GetRenderWindow()->Render();
    QApplication::processEvents();

    updateOverlay();
	statusbar->showMessage("Ready");
	QApplication::processEvents();
}

void MainWindow::updateOverlay()
{
	QVector<QPixmap> base=extractSlices(logic->visualizing, Qt::white);
    //if (logic->masks)
    //{
    //    Vertebra::maskLock.lock();
    //    QVector<QPixmap> overlay=extractSlices(logic->masks, QColor(255,255,0,128));
    //    Vertebra::maskLock.unlock();
    //    for (int i=0; i<base.size(); i++)
	   // {
		  //  QPainter p(&base[i]);
    //        p.drawPixmap(0,0,overlay[i]);
    //        p.end();
	   // }
    //    overlay.clear();
    //}
    painter->set_slices(base);
}

bool MainWindow::checkImageLoaded()
{
	if (!logic->current)
	{
		QMessageBox::warning(this, "No image loaded", "You have to open an image before using this operation.");
		return false;
	}
	return true;
}

void MainWindow::applyFilter(InternalImageType::Pointer &image, itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer filter)
{
	if (!checkImageLoaded())
		return;
	statusbar->showMessage("Applying filter: "+QString(filter->GetNameOfClass()));
	filter->SetInput( image );
	filter->Update();
	image=filter->GetOutput();
	statusbar->showMessage("Ready");
}

void MainWindow::on_actionAnisotropic_diffusion_triggered()
{
	typedef itk::CurvatureAnisotropicDiffusionImageFilter < InternalImageType, InternalImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetNumberOfIterations(5);
	filter->SetTimeStep( 0.0625 );
	filter->SetConductanceParameter( 2.5 );
	applyFilter(logic->current, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
	updateVisualization();
	updateOverlay();
}

void MainWindow::on_actionMedian_denoising_triggered()
{
	typedef itk::MedianImageFilter < InternalImageType, InternalImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	applyFilter(logic->current, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
	updateVisualization();
	updateOverlay();
}

void MainWindow::on_actionGaussian_smoothing_triggered()
{
	typedef itk::SmoothingRecursiveGaussianImageFilter < InternalImageType, InternalImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetSigmaArray(logic->current->GetSpacing());
	applyFilter(logic->current, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
	updateVisualization();
	updateOverlay();
}

void MainWindow::on_actionLogarithmic_rescaling_triggered()
{
	//statusbar->showMessage("Applying logarithmic rescaling (x=log(1+x) voxelwise)");
	Log1Type::Pointer filter = Log1Type::New();
	applyFilter(logic->current, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    applyFilter(logic->lImage, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    applyFilter(logic->hImage, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    logic->maxValue=log(1.0+logic->maxValue);
	updateVisualization();
	updateOverlay();
}
class SchlickFunctor
{
public:
   float operator()( float input )
   {
	   const float b=1.005; //brightness
	   return b*input/(1+(b-1)*input);
   }
};
void MainWindow::on_actionSchlick_URQ_rescaling_triggered()
{
	statusbar->showMessage("Applying Schlick rescaling");
	typedef itk::UnaryFunctorImageFilter<InternalImageType, InternalImageType, SchlickFunctor> FilterType;
	FilterType::Pointer filter = FilterType::New();
	applyFilter(logic->current, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    applyFilter(logic->lImage, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    applyFilter(logic->hImage, (itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer)filter);
    SchlickFunctor sf;
    logic->maxValue=sf(logic->maxValue);
	updateVisualization();
	updateOverlay();
}

void MainWindow::on_actionClear_polygonal_data_triggered()
{
    vtkProp *volume=vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetVolumes()->GetLastProp();
    vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetViewProps()->RemoveAllItems();
    vis->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddVolume(volume);
    vis->GetRenderWindow()->Render();
}