#ifndef __slice_painter_h_
#define __slice_painter_h_
#include "ui_slicePainter.h"
#include "declarations.h"
#include <QPushButton>


class slicePainter : public QWidget, public Ui::slicePainter
{
	Q_OBJECT
		
public:
	slicePainter(QWidget *parent);

    std::vector<CenterRadiusType> centers;

    void setVertebraCenters(const std::vector<CenterRadiusType> & centers);
    void clearInitialization();

	int getSliceNumber() const { return sliceSlider->value(); }
    float spacingX, spacingY, mulX, mulY;
    QPainterPath centerline;

public slots:
	void set_slices(QVector<QPixmap> volume_slices, float spacingX, float spacingY);
	QPixmap get_slice(int index);
	void set_slice(int index, QPixmap slice);
	void saveSlices(std::string filename);
    void startSegmentationButtonClicked();
	void paintSlice();
	//The slots for bounding geometry and clipping.
	void clipIsChecked();
	void radius_Slider_moved(int i);
	//
signals:    
	void vertebraCentersPicked(int labelIndex);
	void showStatusMessage(const QString & msg, int timeout = 0);
	//The signals for communication with the MainLogic-Class.
	void clipWithBoundingGeometry(bool triggered);
	void changeBGRadius(double delta);
    void centersChanged();

private:
	QVector<QPixmap> slices;
    //spacings and aspect ratio corrective multiplicative factors

private slots:
	void slider_moved();

protected:	
	void mousePressEvent(QMouseEvent * event);	
	void mouseReleaseEvent(QMouseEvent * event);	
	void mouseMoveEvent(QMouseEvent * event);
    int movingPoint;	//-1 = no, n = centers[n]
	QPointF lastMousePosOnSlice;
};
#endif