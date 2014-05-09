#pragma once
#include "declarations.h"
#include "vertebra.h"
#include <vtkVolume.h>

#include <QtCore/QObject>
#include <QTime>
#include "mainWindow.h"

class MainLogic : public QObject
{
	Q_OBJECT

public:
	MainLogic(MainWindow& mainWindow);
	InternalImageType::Pointer current, lImage, hImage;
    std::vector<class Vertebra *> vertebra;
    Vertebra *vertebra1;
    unsigned vi; //vertebra index
	VisualizingImageType::Pointer visualizing, binVertebra, masks;
    InternalImageType::PixelType maxValueOriginal, maxValue, maxValue2, maxGM;
	std::string volume_filename, filename_base, series_description;
    static std::string tempDir;
	MainWindow& mainForm;
    QTime time;
    double dist2D;
    ofstream f; //for outputting time information
    bool useShape;

    //classifier output images
    InternalImageType::Pointer cl_LH1, cl_LH2, cl_dfMult, cl_dfCanny; //per dataset

private:
	MainLogic();
																						
	void startSegmentation(const QPointF & center2D, QPainterPath & path, int sliceNumber, bool segmentOtherVertebrae, bool calculateThresholds,
						   int labelIndex, InternalImageType::PixelType diskValue = -1, InternalImageType::PixelType cordValue = -1);	//-1 => values not set

    void segmentationLoop(int direction, const vec3& lastEndplateVector, double sizeGrowth);
	void segmentVertebrae(vec3& originalCenter, bool segmentOtherVertebrae = true);
	double segmentOneButterfly(double dist2d, bool useMRM, unsigned call);
    void inflate(float stepSize);
    void calcDistanceField();
    void calcThresholds(const QPainterPath &path);
	void calcStats(const QPainterPath &path,
	            InternalImageType::ValueType& minv, InternalImageType::ValueType& maxv,
	            InternalImageType::ValueType& avg, InternalImageType::ValueType& stddev);
    void recalcThresholds(Cell *qe, vec3 center, double lastDist);
	void afterOpen(const itk::MetaDataDictionary &metaData);
	void saveLHfiles(int i);
    void postSegmentation();

public slots:
    void clearMasks();
    void volumeOpen(std::string filename);
	void dicomOpen(std::string filename);
	void approximateBoundaryDrawn(int sliceNumber, int labelIndex);
	void vertebraCentersPicked(int sliceNumber, int labelIndex);
	void threePointsPicked(int sliceNumber, int labelIndex);
    void singlePointPicked(int sliceNumber, int labelIndex);
};