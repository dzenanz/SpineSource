#include "slicePainter.h"
#include <QPainter>
#include <QCursor>
#include <QMouseEvent>
#include <QString>
#include <QMessageBox>
#include <QPushButton>
#include "itkMacro.h"
#include <math.h>
#include <assert.h>

#include <assert.h>
#include "declarations.h"

slicePainter::slicePainter(QWidget *parent)
:QWidget(parent)
{
	setupUi(this);

    connect(startButton, SIGNAL(clicked()), SLOT(startSegmentationButtonClicked()));
	connect(sliceSlider, SIGNAL(valueChanged(int)), SLOT(slider_moved()));
	connect(clip_checkBox, SIGNAL(stateChanged(int)), SLOT(clipIsChecked()));
	connect(radius_Slider, SIGNAL(valueChanged(int)), SLOT(radius_Slider_moved(int)));
}

void slicePainter::setVertebraCenters(const std::vector<CenterRadiusType> & centers)
{
    //clearInitialization();

    if (!slices.empty())
	{
        this->centers = centers;
		startButton->setEnabled(true);

		paintSlice();
    }
}

void slicePainter::clearInitialization()
{
	startButton->setEnabled(false);

	centers.clear();
    centerline=QPainterPath();
	movingPoint = -1;

	paintSlice();
}

void slicePainter::slider_moved()
{
	if (slices.size()>0)
	{
        slice->setPixmap(slices.at(sliceSlider->value())
            .scaled(slices[0].width()*mulX,slices[0].height()*mulY,
                    Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
        paintSlice();
	}
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

void slicePainter::set_slices(QVector<QPixmap> volume_slices, float spacingX, float spacingY)
{
	slices=volume_slices;
    this->spacingX=spacingX;
    this->spacingY=spacingY;
    float minSp=std::min(spacingX, spacingY);
    if (this->spacingY/minSp*slices[0].height()>=513)
    {
        minSp=sqrt(spacingX*spacingY);
    }
    this->mulX=spacingX/minSp;
    this->mulY=spacingY/minSp;

    sliceSlider->setMaximum(slices.size() - 1);
	slider_moved();
    this->layout()->update();
	this->updateGeometry();
}

void slicePainter::startSegmentationButtonClicked()
{
	if (slices.size()>0)
	{
		//QMessageBox::warning(this, "path",QString::number(path->elementCount())); //debug
		try
		{
            emit vertebraCentersPicked(bottomLabel->currentIndex());
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

void slicePainter::paintSlice()
{
    if (!slices.empty())
	{
        QPixmap newSlice(slices[sliceSlider->value()].scaled(slices[0].width()*mulX,slices[0].height()*mulY,
                         Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
        QPainter p(&newSlice);
        p.setRenderHint(QPainter::Antialiasing);

        QPen penC(Qt::blue);
        penC.setWidth(2);
        p.setPen(penC);
        p.drawPath(centerline);

        p.setOpacity(0.5);
        QPen pen(Qt::yellow);
        pen.setWidth(2);
        p.setPen(pen);

        for (std::vector<CenterRadiusType>::iterator it = centers.begin(),
		        end = centers.end(); it != end; ++it)
        {
            if (it->z==sliceSlider->value())
                p.setBrush(Qt::red);
            else
                p.setBrush(Qt::transparent);
            p.drawEllipse(QPointF(it->x*mulX, it->y*mulY), mulX*circleRadius/spacingX, mulY*circleRadius/spacingY);
        }
        p.end();

        slice->setPixmap(newSlice);
	}
}

void slicePainter::mousePressEvent(QMouseEvent *event)
{
    if (slices.size()>0 && centers.size() > 0)
	{
		movingPoint = -1;

		QPointF mousePos = slice->mapFromGlobal(event->globalPos());
		QPointF distVec;
		qreal   dist;

		unsigned int n = centers.size();
		for (unsigned int i = 0; i < n; ++i)
		{
            CenterRadiusType cr=centers.at(i);
            distVec = (mousePos - QPointF(cr.x*mulX, cr.y*mulY));
			dist	= sqrt(distVec.x()*distVec.x()+distVec.y()*distVec.y());

			if (dist <= circleRadius/sqrt(spacingX*spacingY/(mulX*mulY)))
            {
				movingPoint = i;
				break;
			}
		}
	}
}

void slicePainter::mouseReleaseEvent(QMouseEvent * event)
{
    if (movingPoint!=-1 && event->button()==Qt::RightButton)
    {
        centers.erase(centers.begin()+movingPoint);
        if (centers.size()==0)
        {
            startButton->setEnabled(false);
            movingPoint=-1;
        }
    }
    else if (movingPoint == -1 && slices.size()>0
        && slice->geometry().contains(event->pos()) && event->button()!=Qt::RightButton)
	{
        QPointF p=slice->mapFromGlobal(event->globalPos());
        p.setX(p.x()/mulX);
        p.setY(p.y()/mulY);
		centers.push_back(
            CenterRadiusType(p.x(), p.y(), sliceSlider->value(), 0));

		startButton->setEnabled(true);
	}
    
    //re-order centers bottom to top
    for (int i=0; i<centers.size()-1; i++)
        for (int k=i+1; k<centers.size(); k++)
            if (centers[i].y<centers[k].y)
                std::swap(centers[i], centers[k]);

    emit centersChanged();
    paintSlice();
}

void slicePainter::mouseMoveEvent(QMouseEvent *event)
{
    if (centers.size() > 0 && slice->geometry().contains(event->pos()) &&
        event->button()!=Qt::RightButton)
	{
		if (movingPoint >= 0)
		{
			QPointF mousePos = slice->mapFromGlobal(event->globalPos());
            mousePos.setX(mousePos.x()/mulX);
            mousePos.setY(mousePos.y()/mulY);
            centers.at(movingPoint).x = mousePos.x();
            centers.at(movingPoint).y = mousePos.y();

            emit centersChanged();
			paintSlice();
		}
	}

}

void slicePainter::clipIsChecked()
{
	if (clip_checkBox->checkState() == Qt::Checked)
		emit clipWithBoundingGeometry(true);
	else
		emit clipWithBoundingGeometry(false);
}

void slicePainter::radius_Slider_moved(int i)
{
	emit changeBGRadius(i);
}