#ifndef __main_window_h_
#define __main_window_h_
#include "declarations.h"
#include "ui_mainWindow.h"
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
#include <vtkSmartPointer.h>
#include <string>

class MainLogic;

class MainWindow : public QMainWindow, public Ui::mainWindow
{
	Q_OBJECT

public:
	using Ui::mainWindow::centralWidget;
	MainWindow(QWidget *parent = 0);
	~MainWindow();
	void defaultCameraPos();
	void updateVisualization();
    void updateOverlay();
    void openInitialization(const char * initFilename);
	MainLogic *logic;

protected:
	void closeEvent(QCloseEvent *event);
	bool checkImageLoaded();
	void applyFilter(InternalImageType::Pointer &image, itk::ImageToImageFilter<InternalImageType,InternalImageType>::Pointer filter);

    vtkVolume *volume;
	vtkSmartPointer<vtkVolumeProperty> transFunc;
	unsigned char vThLow, vThHigh; //visualization threshold (TF parameter)
	float opacity;
	void constructTF();

signals:
	void volumeOpened(std::string filename);
	void dicomOpened(std::string directory);

public slots:
	void on_actionAnisotropic_diffusion_triggered();
	void on_actionMedian_denoising_triggered();
    void on_actionClear_polygonal_data_triggered();
	void saveInitialization(std::string vertebraLabel, const char * volumeFilename);

private slots:
	void on_actionOpen_triggered();
	void on_actionOpen_DICOM_series_triggered();
	void on_actionOpen_initialization_triggered();	
	void on_actionSave_slices_triggered();
	void on_actionScreenshot_triggered();
	void on_actionGaussian_smoothing_triggered();
	void on_actionLogarithmic_rescaling_triggered();
	void on_actionSchlick_URQ_rescaling_triggered();
	void on_actionLoad_polygonal_mesh_triggered();
    void on_actionQuick_help_triggered();
	void on_actionTF_shift_Up_triggered();
	void on_actionTF_shift_Down_triggered();
	void on_actionOpacityPlus_triggered();
	void on_actionOpacityMinus_triggered();
	void on_actionReset_TF_triggered();
};
#endif