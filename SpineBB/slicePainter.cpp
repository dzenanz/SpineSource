#include "slicePainter.h"
#include <QPainter>
#include <QCursor>
#include <QMouseEvent>
#include <QString>
#include <QMessageBox>
#include "itkMacro.h"
slicePainter::slicePainter(QWidget *parent)
:QWidget(parent)
{
	setupUi(this);

	//connect signals and slots
	connect(sliceSlider, SIGNAL(valueChanged(int)), SLOT(slider_moved()));
	connect(sliceSlider, SIGNAL(sliderMoved(int)), SLOT(slider_moved()));
	
	prefSize=layoutWidget->size();
	VertebraBrushSize=1;
	painting=false;
}

void slicePainter::slider_moved()
{
	lcdSlice->display(sliceSlider->value());
	if (slices.size()>0)
		slice->setPixmap(slices.at(sliceSlider->value()));
}

QPixmap slicePainter::get_slice(int index)
{
	return slices[index];
};

void slicePainter::set_slice(int index, QPixmap slice)
{
	slices[index]=slice;
	slider_moved();
}

void slicePainter::saveSlices(std::string filename)
{
	if (slices.size()>0)
	{
		int pos=filename.find_last_of('.');
		if (pos==filename.npos)
			QMessageBox::warning(this, tr("Filename problem"), tr("Could not determine extension from the given filename."));
			
		else
		{
			QString fname=QString::fromStdString(filename.substr(0, pos)+"%1"+filename.substr(pos));

			int pad=1;
			if (slices.size()<=100)
				pad=2;
			else if (slices.size()<=1000)
				pad=3;
			else 
				pad=4;

			for (int i=0; i<slices.size(); i++)
				slices[i].save(fname.arg(i,pad,10, QChar('0')));
		}
	}
	else
		QMessageBox::warning(this, tr("No slices present"), tr("There are no slices currently present. You should load a dataset first."));
}

void slicePainter::set_slices(QVector<QPixmap> volume_slices)
{
	slices=volume_slices;

	int nw=slices[0].width(), nh=slices[0].height();
	int ow=slice->width(), oh=slice->height();
	if (nw!=ow || nh!=oh)
		prefSize=QSize(layoutWidget->width()-ow+nw, layoutWidget->height()-oh+nh);

	this->sliceSlider->setMaximum(slices.size()-1);
	sliceSlider->setValue(slices.size()/2);
	slider_moved();
	layoutWidget->updateGeometry();
	this->updateGeometry();
}

QSize slicePainter::sizeHint() const
{
	return prefSize;
}

void slicePainter::resizeEvent(QResizeEvent *re)
{
	if (slices.size()>0)
	{
		slice->setMinimumSize(slices[0].width(),slices[0].width());
		layoutWidget->adjustSize();
		layoutWidget->updateGeometry();
		this->adjustSize();
		this->updateGeometry();
	}
}

void slicePainter::paintSlice()
{
	QPixmap newSlice(slices[sliceSlider->value()]);
	QPainter p(&newSlice);
	QPen pen(QColor(255,0,0,128),0);
	p.setPen(pen);
	//p.setPen(Qt::NoPen);
	//p.setBrush(QBrush(Qt::red));
	//path.setFillRule(Qt::WindingFill);
	p.drawPath(path);
	slice->setPixmap(newSlice);
	//emit minorUpdate(sliceSlider->value());
}

void slicePainter::paintEvent(QPaintEvent *event)
{
	if (slices.size()>0)
		paintSlice();
}

void slicePainter::mousePressEvent(QMouseEvent *event)
{
	if (VertebraBrushSize>0 && slices.size()>0 && slice->geometry().contains(event->pos()))
	{
		path=QPainterPath(event->pos()-slice->pos());
		painting=true;
		paintSlice();
	}
	else
		QApplication::beep();
}

void slicePainter::mouseReleaseEvent(QMouseEvent *event)
{
	if (VertebraBrushSize>0 && slices.size()>0 && painting)
	{
		painting=false;
		//QMessageBox::warning(this, "path",QString::number(path.elementCount())); //debug
		try
		{
            emit majorUpdate(sliceSlider->value(), bottomLabel->currentIndex());
		}
		catch (itk::ExceptionObject & exc)
		{
			QMessageBox::warning(0, "Exception", exc.what());
		}
		catch (std::exception & exc)
		{
			QMessageBox::warning(0, "Exception", exc.what());
		}
	}
}

void slicePainter::mouseMoveEvent(QMouseEvent *event)
{
	if (painting && slices.size()>0 && slice->geometry().contains(event->pos()))
	{
		path.lineTo(event->pos()-slice->pos());
		paintSlice();
	}
}