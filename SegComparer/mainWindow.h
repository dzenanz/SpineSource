#ifndef __main_window_h_
#define __main_window_h_
#include "ui_mainWindow.h" 


enum ViewDirection
{
    AxialView,
    CoronalView,
    SaggitalView
};


struct Slices
{
    QVector<QImage> axial;
    QVector<QImage> coronal;
    QVector<QImage> saggital;
};

const QColor c1(255,0,0,96); //translucent red
const QColor c2(0,255,0,64); //translucent green
//so overlap is yellow

class MainWindow : public QMainWindow, public Ui::mainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
    Slices slices;
    
public slots:
	void baseOpen();
	void seg1Open();
	void seg2Open();
	void sliceChanged(int index);
	void blend(Slices seg);
	void saveClicked();
	void calcDSC();
    void viewDirRadioButtonToggled(bool);
    float DiceSimilarityCoefficient();

private:
    ViewDirection viewDir;
};


#endif //__main_window_h_