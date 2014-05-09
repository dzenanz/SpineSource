#ifndef __slice_painter_h_
#define __slice_painter_h_
#include "ui_slicePainter.h" 


class slicePainter : public QWidget, public Ui::slicePainter
{
	Q_OBJECT
		
public:
	slicePainter(QWidget *parent);
	
	~slicePainter();

	QPainterPath * getPath() { return path; }

	void setPath(QPainterPath * p) { delete path; path = p; }

	void setPath(const QPainterPath & p) { delete path; path = new QPainterPath(p); }


	//init method getters

    const QPointF & getSinglePoint() const { return initMethod_PlaceSinglePoint.point; }

	const std::vector<QPointF> & getVertebraCenters() const { return initMethod_PickVertebraCenters.centers; }
    
	void getThreePointInitialization(QPointF & vertebraCenter, QPointF & disk, QPointF & cord, qreal & circleRadius) const
	{
		vertebraCenter = initMethod_PickThreePoints.vertebraCenter;
		disk = initMethod_PickThreePoints.disk;
		cord = initMethod_PickThreePoints.cord;
		circleRadius = initMethod_PickThreePoints.circleRadius;
	}

    //this getter is just used to save the initialization, the logic actually just needs the path
    void getCenterAndTwoBoundaryPoints(QPointF & c, QPointF & b1, QPointF & b2)
    {
        c = initMethod_PlaceCenterAndTwoBP.center;
        b1 = initMethod_PlaceCenterAndTwoBP.boundaryPoint1;
        b2 = initMethod_PlaceCenterAndTwoBP.boundaryPoint2;
    }


	// init method setters (useful when loading an initialization)

	void setVertebraOutline(const QPainterPath & outlinePath);

    void setSinglePoint(const QPointF & p);

	void setCenterAndTwoBoundaryPoints(const QPointF & c, const QPointF & b1, const QPointF & b2);
	
    void setVertebraCenters(const std::vector<QPointF> & centers);

	void setThreePointInitialization(const QPointF & vertebraCenter, const QPointF & disk, const QPointF & cord, qreal circleRadius);


	void clearInitialization();

	enum InitializationMethod
	{
		DrawVertebraOutline,
        PlaceSinglePoint,
		PlaceCenterAndTwoBP,	//place center + 2 boundary points
		PickVertebraCenters,
		PickThreePoints     	//pick points: vertebral center, disk, spinal cord
	};

	InitializationMethod getInitMethod(){ return initializationMethod; }

	int getSliceNumber() const { return sliceSlider->value(); }

public slots:
	void set_slices(QVector<QPixmap> volume_slices);
	QPixmap get_slice(int index);
	void set_slice(int index, QPixmap slice);
	void saveSlices(std::string filename);
    void startButtonClicked();

signals:    
	void approximateBoundaryDrawn(int sliceNumber, int labelIndex);	 
    void singlePointPicked(int sliceNumber, int labelIndex);
	void vertebraCentersPicked(int sliceNumber, int labelIndex);
	void threePointsPicked(int sliceNumber, int labelIndex);
	void showStatusMessage(const QString & msg, int timeout = 0);

private:
	QVector<QPixmap> slices;
	QPainterPath * path;

	//! The number of steps that was performed during the user initialization.
	//! 0 indicates that no initialization has been done.
	//! The number of steps that are needed to start the segmentation depend on
	//! the particular initialization method.
	int stepsPerformed;

	//! The Method that is used by the user to initialize the segmentation
	InitializationMethod initializationMethod;

private slots:
	void slider_moved();

    void on_placeSinglePointRadioButton_toggled(bool checked);
	void on_drawVertebraOutlineRadioButton_toggled(bool checked);
	void on_placeCenterAndTwoBPRadioButton_toggled(bool checked);
	void on_pickVertebraCentersRadioButton_toggled(bool checked);
	void on_pickThreePointsRadioButton_toggled(bool checked);

protected:	
	void mousePressEvent(QMouseEvent * event);	
	void mouseReleaseEvent(QMouseEvent * event);	
	void mouseMoveEvent(QMouseEvent * event);
	void paintSlice();

private:
	QPointF lastMousePosOnSlice;

    struct InitMethod_PlaceSinglePoint
	{
        QPointF point;
    } initMethod_PlaceSinglePoint;

	struct InitMethod_PlaceCenterAndTwoBP
	{
		QPointF center, boundaryPoint1, boundaryPoint2;

		int movingPoint;	//-1 = no, 0 = center, 1 = bp1, 2 = bp2
	
		//! guess a shape that approximates a default vertebra
		QPainterPath * createBoundaryPath();
	} initMethod_PlaceCenterAndTwoBP;	

	struct InitMethod_PickVertebraCenters
	{
		std::vector<QPointF> centers;

		int movingPoint;	//-1 = no, n = centers[n]
	} initMethod_PickVertebraCenters;

	struct InitMethod_PickThreePoints
	{
		QPointF vertebraCenter, disk, cord;

		int movingPoint;	//-1 = no, 0 = body, 1 = disk, 2 = cord, 3 = scaling circle
		
		qreal circleRadius;
	} initMethod_PickThreePoints;
};
#endif