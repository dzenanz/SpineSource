#ifndef __slice_painter_h_
#define __slice_painter_h_
#include "ui_slicePainter.h" 

class slicePainter : public QWidget, public Ui::slicePainter
{
	Q_OBJECT

public:
	slicePainter(QWidget *parent);
	QSize sizeHint () const;

public slots:
	void set_slices(QVector<QPixmap> volume_slices);
	QPixmap get_slice(int index);
	void set_slice(int index, QPixmap slice);
	void saveSlices(std::string filename);

signals:
	void majorUpdate(int sliceNumber, int labelIndex);

private slots:
	void slider_moved();

private:
	QSize prefSize;
	int VertebraBrushSize;
	bool painting;
public:
	QPainterPath path;
	QVector<QPixmap> slices;

protected:
	void resizeEvent(QResizeEvent *re);
	void mousePressEvent(QMouseEvent * event);
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent * event);
	void paintSlice();
	void paintEvent(QPaintEvent * event);
};
#endif