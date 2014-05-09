#pragma once
#include "declarations.h"
#include "vertebra.h"
#include <vtkVolume.h>
#include <vtkMutexLock.h>
#include <vtkActor.h>

#include <QtCore/QObject>
#include <QTime>
#include <QTimer>
#include <QMutex>
#include <QWaitCondition>
#include "mainWindow.h"


class MainLogic : public QObject
{
	Q_OBJECT

public:
	MainLogic(MainWindow& mainWindow);
	InternalImageType::Pointer current, lImage, hImage;
    std::vector<class Vertebra *> vertebra;
	VisualizingImageType::Pointer visualizing, masks;
    InternalImageType::PixelType maxValueOriginal, maxValue, maxValue2, maxGM;
	std::string volume_filename, filename_base, series_description;
    static std::string tempDir;
	MainWindow& mainForm;
    QTime time;
    ofstream f; //for outputting time information
    bool fullAuto;
    int maxRefinementLevel;

    //for multithreaded segmentation
    QTimer *timer;
    QMutex threadLock;
    QWaitCondition wc;
    bool render;
    bool updateOverlay;
    Eigen::VectorXd x, y;

    //classifier output images
    InternalImageType::Pointer cl_LH1, cl_LH2, cl_dfMult, cl_dfCanny; //per dataset

private:
	MainLogic();
																						
	double segmentOneButterfly(Vertebra *v, unsigned call);
    void inflate(Vertebra *v, float stepSize, bool useShape=false, float force=0);
    void calcThresholds(const QPainterPath &path, Vertebra *v);
	void afterOpen(const itk::MetaDataDictionary &metaData);
    void calculateFeatures();
    void calculatePerVertebraFeatures(int vertebraIndex);
	void saveLHfiles(int i);
    void postSegmentation();
    void removeVertebra(unsigned index);
    void vertebraeToCenters();
    void centersToVertebrae();
    void drawCenterline();

public slots:
    void recalcCenterline();
    void renderNeeded();
    void clearMasks();
    void volumeOpen(std::string filename);
	void dicomOpen(std::string filename);
    void detectCenters();
    void autoInitSegmentation();
	void vertebraCentersPicked(int labelIndex);
	void clipWithBoundingGeometry(bool triggered);
	void changeBGRadius(double radius);
};