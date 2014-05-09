#include "declarations.h"
#include "vertebra.h"
#include <vtkVolume.h>
#include <vtkActor.h>

#include <QtCore/QObject>
#include <QTime>
#include "mainWindow.h"

class MainLogic : public QObject
{
	Q_OBJECT

public:
	MainLogic(MainWindow& mainWindow);
	InternalImageType::Pointer current;
    std::vector<Vertebra *> vertebra;
    Vertebra *vertebra0;
	VisualizingImageType::Pointer visualizing;
	InternalImageType::Pointer lImage, hImage;
    InternalImageType::PixelType maxValue, maxValue2;
	std::string volume_filename, filename_base;
    static std::string tempDir;

private:
	MainLogic();
	MainWindow& mainForm;
    void segmentationLoop(const vec3& last2Center, const Vertebra& last,
                          const vec3& lastAdjecantVertebraVector, double sizeGrowth, int labelIndexDirection);
	void segmentVertebrae(vec3& originalCenter, double dist2D, int labelIndex);
    void postSegmentation();
	double segmentOne(Vertebra *vertebra);
	double segmentOneButterfly(Vertebra *vertebra);
    void inflate(Cell *qe, vec3& center, float avgEdgeLen);

    InternalImageType::PixelType highThreshold, lowThreshold;
    QTime time;
    ofstream f; //for outputting time information

public slots:
	void volumeOpen(std::string filename);
	void majorSliceUpdate(int sliceNumber, int labelIndex);
};