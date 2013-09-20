#define _SCL_SECURE_NO_WARNINGS
#include "mainWindow.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QPainter>
#include <QColor>

#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDatImageIOFactory.h"
#include "itkPGMImageIOFactory.h"
#include "itkOrientImageFilter.h"

typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Image<unsigned char, 2> SliceImageType;
typedef itk::Image<float, 3> InternalImageType;
double scaleFactorX[3], scaleFactorY[3];

using namespace std;

MainWindow::MainWindow(QWidget *parent)
:QMainWindow(parent)
{
	setupUi(this);

	connect(base, SIGNAL(clicked()),SLOT(baseOpen()));
	connect(seg1, SIGNAL(clicked()),SLOT(seg1Open()));
	connect(seg2, SIGNAL(clicked()),SLOT(seg2Open()));
	connect(slice, SIGNAL(valueChanged(int)),SLOT(sliceChanged(int)));
	connect(save, SIGNAL(clicked()),SLOT(saveClicked()));
	connect(dsc, SIGNAL(clicked()),SLOT(calcDSC()));
    connect(axialButton,    SIGNAL(toggled(bool)), SLOT(viewDirRadioButtonToggled(bool)));
    connect(coronalButton,  SIGNAL(toggled(bool)), SLOT(viewDirRadioButtonToggled(bool)));
    connect(saggitalButton, SIGNAL(toggled(bool)), SLOT(viewDirRadioButtonToggled(bool)));

	itk::ObjectFactoryBase::RegisterFactory( itk::DatImageIOFactory::New() );
    itk::ObjectFactoryBase::RegisterFactory( itk::PGMImageIOFactory::New() );
} 

QVector<QImage> extractSlices(VisualizingImageType::Pointer itkImage, QColor color, ViewDirection view)
{
    itk::OrientImageFilter<VisualizingImageType, VisualizingImageType>::Pointer orientationFilter =
        itk::OrientImageFilter<VisualizingImageType, VisualizingImageType>::New();

    //swap / flip axes, if desired

    //this call determines the image direction by checking the orientation matrix of the image
    orientationFilter->UseImageDirectionOn();

    if (view == AxialView)
       orientationFilter->SetDesiredCoordinateOrientationToAxial();
    else if (view == CoronalView)
       orientationFilter->SetDesiredCoordinateOrientationToCoronal();
    else if (view == SaggitalView)
       orientationFilter->SetDesiredCoordinateOrientationToSagittal();
  
    orientationFilter->SetInput(itkImage);
        
    orientationFilter->UpdateLargestPossibleRegion();

    VisualizingImageType::Pointer rotatedImage = orientationFilter->GetOutput();

    VisualizingImageType::SizeType size = rotatedImage->GetLargestPossibleRegion().GetSize();

	QVector<QImage> result(size[2]);
	VisualizingImageType::IndexType ind;
	unsigned char val;
	QRgb rgb;

    for (int z=0; z<size[2]; z++)
        result[z]=QImage(size[0], size[1], QImage::Format_ARGB32);

	for (int z=0; z<size[2]; z++)
	{		
        ind[2]=z;

		for (int y=0; y<size[1]; y++)
		{			
            ind[1]=y;

			for (int x=0; x<size[0]; x++)
			{
                ind[0]=x;

                val=rotatedImage->GetPixel(ind);

                if (color.alpha()==255)
					rgb=QColor(color.red()*val/255.0, color.green()*val/255.0, color.blue()*val/255.0).rgba();
				else
					rgb=QColor(color.red(), color.green(), color.blue(), color.alpha()*val/255.0).rgba();

                result[z].setPixel(x, y, rgb);
			}
		}
	}


    //maintain voxel spacing (this is not automatically rotated along with the orientation!)
    
    itk::Vector<double, 3> spacing = itkImage->GetSpacing();

    itk::FixedArray<unsigned int> permuteAxes = orientationFilter->GetPermuteOrder();

    double scaleX = spacing[permuteAxes[0]],
           scaleY = spacing[permuteAxes[1]];
   
    double minSize = std::min<double>(scaleX, scaleY);

    //the minimal size of both is mapped to 1 pixel
    scaleX /= minSize;
    scaleY /= minSize;
    
    scaleFactorX[view] = scaleX;
    scaleFactorY[view] = scaleY;

	return result;
}

Slices openHelper(string filename, QColor color, bool sagOnly=false)
{
	typedef  itk::ImageFileReader< InternalImageType >  ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( filename.c_str() );
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		QMessageBox mb(QMessageBox::Critical, QString("Error opening file"),
			QString("Error opening file ")+QString::fromStdString(filename)+"\nException content:\n"+e.GetDescription(),
            QMessageBox::Ok);
		mb.exec();
	}
	typedef itk::RescaleIntensityImageFilter < InternalImageType, VisualizingImageType > RescaleImageFilterType;
	RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
	rescale->SetInput( reader->GetOutput() );
	rescale->Update();

    Slices slcs;

    slcs.saggital = extractSlices(rescale->GetOutput(), color, SaggitalView);
    if (!sagOnly)
    {
        slcs.axial    = extractSlices(rescale->GetOutput(), color, AxialView);
        slcs.coronal  = extractSlices(rescale->GetOutput(), color, CoronalView);
    }

	return slcs;
}

void MainWindow::baseOpen()
{
	string filename=QFileDialog::getOpenFileName(this, tr("Open original image"), "", tr("3D images (*.dat *.mha *.mhd *.nrrd *.nii *.hdr *.dcm);;All files (*.*)")).toStdString();
	if (filename=="")
		return;
	slices=openHelper(filename, Qt::white);
	setWindowTitle(QString("Segmentation Comparer - ")+QString::fromStdString(filename));
	
    //setup viewDir and gui elements
    viewDirRadioButtonToggled(true);
}

void MainWindow::seg1Open()
{
	if (slices.axial.size()==0)
		QMessageBox::warning(0, "Error", "Open base (original) image first!");
	else
	{
		string filename = QFileDialog::getOpenFileName(this, tr("Open segmentation 1"), ".", tr("3D images (*.dat *.mha *.mhd *.nrrd *.nii *.hdr *.dcm);;All files (*.*)")).toStdString();
		if (filename=="")
			return;
		Slices seg=openHelper(filename, c1); 
		blend(seg);
	}
}

void MainWindow::seg2Open()
{
	if (slices.axial.size()==0)
		QMessageBox::warning(0, "Error", "Open base (original) image first!");
	else
	{
		string filename=QFileDialog::getOpenFileName(this, tr("Open segmentation 2"), ".", tr("3D images (*.dat *.mha *.mhd *.nrrd *.nii *.hdr *.dcm);;All files (*.*)")).toStdString();
		if (filename=="")
			return;
		Slices seg=openHelper(filename, c2);
		blend(seg);
	}
}

void MainWindow::sliceChanged(int index)
{
	lcd->display(index);

    QVector<QImage> * currentSlices;
    if (viewDir == AxialView)
        currentSlices = &slices.axial;
    else if (viewDir == CoronalView)
        currentSlices = &slices.coronal;
    else
        currentSlices = &slices.saggital;

	if (currentSlices->size()>0)
	{
        QImage img = currentSlices->at(index);

        int w = img.width();
        int h = img.height();

        //this step is necessary since otherwise we would get a (correctly) stretched and auto - smoothed image,
        //but we want to see the nearest neighbour upsampling instead, which is only available for scaling images manually
        img = img.scaled(w * scaleFactorX[viewDir], h * scaleFactorY[viewDir],
                         Qt::IgnoreAspectRatio, Qt::FastTransformation);

        image->setPixmap(QPixmap::fromImage(img));
	}
}

void MainWindow::blend(Slices seg)
{
    for (int i=0; i<seg.axial.size(); i++)    
    {
	    QPainter p(&slices.axial[i]);
        p.drawImage(0,0,seg.axial[i]);
    }

    for (int i=0; i<seg.coronal.size(); i++)    
    {
	    QPainter p(&slices.coronal[i]);
        p.drawImage(0,0,seg.coronal[i]);
    }

    for (int i=0; i<seg.saggital.size(); i++)    
    {
	    QPainter p(&slices.saggital[i]);
        p.drawImage(0,0,seg.saggital[i]);
    }

	sliceChanged(slice->value());
}

void MainWindow::saveClicked()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save comparison slices"), ".", tr("2D images (*.png *.bmp *.jpg)"));
	if (filename!="")
		image->pixmap()->save(filename);
	//{
 //       int pos=filename.toStdString().find_last_of('.');
	//	if (pos==filename.toStdString().npos)
	//		QMessageBox::warning(this, tr("Filename problem"), tr("Could not determine extension from the given filename."));
	//	else
	//	{
 //           QString fname=filename.left(pos)+"%1"+filename.right(filename.length()-pos);
 //           for (int i=0; i<slices.saggital.size(); i++)
 //           {                
 //               QImage img = slices.saggital[i];
 //               img=img.scaled(img.width() * scaleFactorX[SaggitalView],
 //                              img.height() * scaleFactorY[SaggitalView],
 //                   Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
 //               img.save(fname.arg(i));
 //           }
	//	}
	//}
}

float MainWindow::DiceSimilarityCoefficient()
{
    QVector<QImage> * currentSlices;
    if (viewDir == AxialView)
        currentSlices = &slices.axial;
    else if (viewDir == CoronalView)
        currentSlices = &slices.coronal;
    else
        currentSlices = &slices.saggital;
    
    unsigned s0=0, s1=0, s2=0;

	for (int z=0; z<currentSlices->size(); z++)
		for (int y=0; y<(*currentSlices)[0].height(); y++)
			for (int x=0; x<(*currentSlices)[0].width(); x++)
			{
				QColor c((*currentSlices)[z].pixel(x,y));
				if (c.saturation()!=0)
					if (c.hue()==c1.hue())
						s1++; //redish pixel
					else if (c.hue()==c2.hue())
						s2++; //bluish pixel
					else
						s0++; //purplish pixel = intersection
			}
    return (2.0*s0)/(2.0*s0+s1+s2);
}

void MainWindow::calcDSC()
{
	if (slices.axial.size()==0)
		QMessageBox::warning(0, "Error", "Open base (original), and 2 segmented images first!");
	else
	{
	    QMessageBox::information(0, "DSC","Dice Similarity Coefficient = "
            +QString::number(100*DiceSimilarityCoefficient(),'g',4)+"%");
	}
}

void MainWindow::viewDirRadioButtonToggled(bool)
{
    if (axialButton->isChecked())
        viewDir = AxialView;
    else if (coronalButton->isChecked())
        viewDir = CoronalView;
    else
        viewDir = SaggitalView;

    QVector<QImage> * currentSlices;
    if (viewDir == AxialView)
        currentSlices = &slices.axial;
    else if (viewDir == CoronalView)
        currentSlices = &slices.coronal;
    else
        currentSlices = &slices.saggital;

    slice->setMaximum(currentSlices->size()-1);
	slice->setValue(currentSlices->size()/2); //does not trigger valueChanged if old and new value are the same
	sliceChanged(currentSlices->size()/2);
    
    image->setFixedSize(currentSlices->at(0).width()  * scaleFactorX[viewDir],
                        currentSlices->at(0).height() * scaleFactorY[viewDir]);

    //@todo: proper window resizing
    //...
    this->updateGeometry();
}
