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

const qreal DefaultVertebraCenterCircleRadius = 8.0;

slicePainter::slicePainter(QWidget *parent)
:QWidget(parent),
 path(NULL),
 stepsPerformed(0),
 initializationMethod(DrawVertebraOutline)
{
	setupUi(this);

    path = new QPainterPath();

	//connect signals and slots
	connect(sliceSlider, SIGNAL(valueChanged(int)), SLOT(slider_moved()));
    connect(startButton, SIGNAL(clicked()), SLOT(startButtonClicked()));

	//setup default initialization method
	if (drawVertebraOutlineRadioButton->isChecked())
		on_drawVertebraOutlineRadioButton_toggled(true);
	else if (placeCenterAndTwoBPRadioButton->isChecked())
		on_placeCenterAndTwoBPRadioButton_toggled(true);
	else if (pickVertebraCentersRadioButton->isChecked())
		on_pickVertebraCentersRadioButton_toggled(true);
	else if (pickThreePointsRadioButton->isChecked())
		on_pickThreePointsRadioButton_toggled(true);

	clearInitialization();
}

slicePainter::~slicePainter()
{
	delete path;
}


//initialization setters

void slicePainter::setVertebraOutline(const QPainterPath & outlinePath)
{
	clearInitialization();

	initializationMethod = DrawVertebraOutline;

    drawVertebraOutlineRadioButton->toggle();

	if (!slices.empty())
	{
		path = new QPainterPath(outlinePath);
		
		stepsPerformed = 1;
		slice->setCursor(Qt::CrossCursor);
		startButton->setEnabled(true);

		paintSlice();
	}
}

void slicePainter::setSinglePoint(const QPointF & p)
{
    clearInitialization();

    initializationMethod = PlaceSinglePoint;

    placeSinglePointRadioButton->toggle();

    if (!slices.empty())
	{
        initMethod_PlaceSinglePoint.point = p;
        
        stepsPerformed = 1;
		startButton->setEnabled(true);

		paintSlice();
    }
}

void slicePainter::setCenterAndTwoBoundaryPoints(const QPointF & c, const QPointF & b1, const QPointF & b2)
{
	clearInitialization();

	initializationMethod = PlaceCenterAndTwoBP;

    placeCenterAndTwoBPRadioButton->toggle();

	if (!slices.empty())
	{
		initMethod_PlaceCenterAndTwoBP.center = c;
		initMethod_PlaceCenterAndTwoBP.boundaryPoint1 = b1;
		initMethod_PlaceCenterAndTwoBP.boundaryPoint2 = b2;
		
		delete path;
		path = initMethod_PlaceCenterAndTwoBP.createBoundaryPath();
		
		stepsPerformed = 1;
		startButton->setEnabled(true);

		paintSlice();
	}
}

void slicePainter::setVertebraCenters(const std::vector<QPointF> & centers)
{
    clearInitialization();

    initializationMethod = PickVertebraCenters;

    pickVertebraCentersRadioButton->toggle();

    if (!slices.empty())
	{
        initMethod_PickVertebraCenters.centers = centers;
            
		stepsPerformed = 1;
		startButton->setEnabled(true);

		paintSlice();
    }
}

void slicePainter::setThreePointInitialization(const QPointF & vertebraCenter, const QPointF & disk, const QPointF & cord, qreal circleRadius)
{
    clearInitialization();

    initializationMethod = PickThreePoints;

    pickThreePointsRadioButton->toggle();

    if (!slices.empty())
	{
        initMethod_PickThreePoints.vertebraCenter = vertebraCenter;
        initMethod_PickThreePoints.disk           = disk;
        initMethod_PickThreePoints.cord           = cord;
	    initMethod_PickThreePoints.circleRadius   = circleRadius;

        stepsPerformed = 3;
        startButton->setEnabled(true);

		paintSlice();
    }
}

void slicePainter::clearInitialization()
{
	stepsPerformed = 0;
	
	startButton->setEnabled(false);

	delete path;
	path = new QPainterPath();

	initMethod_PickVertebraCenters.centers.clear();	
	initMethod_PickVertebraCenters.movingPoint = -1;

	initMethod_PlaceCenterAndTwoBP.movingPoint = -1;
	
	initMethod_PickThreePoints.movingPoint = -1;

	initMethod_PickThreePoints.circleRadius = DefaultVertebraCenterCircleRadius;

	paintSlice();
}

void slicePainter::slider_moved()
{
	if (slices.size()>0)
	{
		slice->setPixmap(slices.at(sliceSlider->value()));
		clearInitialization();
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

void slicePainter::set_slices(QVector<QPixmap> volume_slices)
{
	slices=volume_slices;
	sliceSlider->setMaximum(slices.size() - 1);
	slider_moved();
    this->layout()->update();
	this->updateGeometry();
}

void slicePainter::on_placeSinglePointRadioButton_toggled(bool checked)
{
    if (checked)
	{
		initializationMethod = PlaceSinglePoint;

		clearInitialization();

		emit showStatusMessage("Please place a single point inside one of the vertebrae.");
	}
}

void slicePainter::on_drawVertebraOutlineRadioButton_toggled(bool checked)
{
	if (checked)
	{
		initializationMethod = DrawVertebraOutline;

		clearInitialization();

		emit showStatusMessage("Please draw an approximate outline inside one of the middle vertebrae.");
	}
}

void slicePainter::on_placeCenterAndTwoBPRadioButton_toggled(bool checked)
{
	if (checked)
	{
		initializationMethod = PlaceCenterAndTwoBP;

		clearInitialization();
		
		slice->setCursor(Qt::CrossCursor);

		emit showStatusMessage("Please place the center point in one of the middle vertebrae (Step 1 / 2) ...");
	}
}

void slicePainter::on_pickVertebraCentersRadioButton_toggled(bool checked)
{
	if (checked)
	{
		initializationMethod = PickVertebraCenters;

		clearInitialization();
		
		slice->setCursor(Qt::CrossCursor);

		emit showStatusMessage("Please place the center points of all vertebrae.");
	}
}

void slicePainter::on_pickThreePointsRadioButton_toggled(bool checked)
{
	if (checked)
	{
		initializationMethod = PickThreePoints;

		clearInitialization();
		
	slice->setCursor(Qt::CrossCursor);

		emit showStatusMessage("Please place a point at a vertebral center (Step 1 / 4) ...");
	}
}

void slicePainter::startButtonClicked()
{
	if (slices.size()>0)
	{
		//QMessageBox::warning(this, "path",QString::number(path->elementCount())); //debug
		try
		{
            int labelIndex=bottomLabel->currentIndex();
            //int labelIndex=((QPushButton *)sender())->toolTip().toInt();
			switch(initializationMethod)
			{
                case PlaceSinglePoint:
                    emit singlePointPicked(sliceSlider->value(), labelIndex);
                    break;

				case DrawVertebraOutline:
				case PlaceCenterAndTwoBP:
					emit approximateBoundaryDrawn(sliceSlider->value(), labelIndex);
					break;

				case PickVertebraCenters:
					emit vertebraCentersPicked(sliceSlider->value(), labelIndex);
					break;
				
				case PickThreePoints:
					emit threePointsPicked(sliceSlider->value(), labelIndex);
					break;
			}	
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
	QPixmap newSlice(slices[sliceSlider->value()]);
	QPainter p(&newSlice);
		p.setRenderHint(QPainter::Antialiasing);
		p.setOpacity(0.5);		
		
		switch (initializationMethod)
		{
			case DrawVertebraOutline:
			{				
				p.setPen(Qt::red);
				p.drawPath(*path);

				break;
			}

            case PlaceSinglePoint:
            {
                if (stepsPerformed > 0)
				{
                    p.setPen(Qt::white);
				    p.setBrush(Qt::red);

                    p.drawEllipse(initMethod_PlaceSinglePoint.point, PointRadius, PointRadius);
                }

                break;
            }

			case PlaceCenterAndTwoBP:
			{
				if (stepsPerformed > 0)
				{	
					p.setPen(Qt::white);
															
					p.setBrush(Qt::green);
					p.drawLine(initMethod_PlaceCenterAndTwoBP.center, initMethod_PlaceCenterAndTwoBP.boundaryPoint1);
					p.drawEllipse(initMethod_PlaceCenterAndTwoBP.boundaryPoint1, PointRadius, PointRadius);
					p.setBrush(Qt::blue);				
					p.drawLine(initMethod_PlaceCenterAndTwoBP.center, initMethod_PlaceCenterAndTwoBP.boundaryPoint2);
					p.drawEllipse(initMethod_PlaceCenterAndTwoBP.boundaryPoint2, PointRadius, PointRadius);
					p.setBrush(Qt::lightGray);
					p.drawEllipse(initMethod_PlaceCenterAndTwoBP.center,		 PointRadius, PointRadius);
					
					p.setPen(Qt::red);
					p.setBrush(Qt::NoBrush);					
					p.drawPath(*path);
				}

				break;
			}

			case PickVertebraCenters:
			{
				p.setPen(Qt::white);
				p.setBrush(Qt::red);

				for (std::vector<QPointF>::iterator it = initMethod_PickVertebraCenters.centers.begin(),
					 end = initMethod_PickVertebraCenters.centers.end(); it != end; ++it)
				{	
					p.drawEllipse(*it, PointRadius, PointRadius);
				}

				break;
			}

			case PickThreePoints:
			{
				if (stepsPerformed > 0)
				{	
					p.setPen(Qt::white);
															
					p.setBrush(Qt::red);					
					p.drawEllipse(initMethod_PickThreePoints.vertebraCenter, PointRadius, PointRadius);
					p.setPen(Qt::red);
					p.setBrush(Qt::NoBrush);
					p.drawEllipse(initMethod_PickThreePoints.vertebraCenter,
								  initMethod_PickThreePoints.circleRadius, initMethod_PickThreePoints.circleRadius);
					p.drawText(initMethod_PickThreePoints.vertebraCenter, "vertebral center");
					p.setPen(Qt::white);

					if (stepsPerformed > 1)
					{
						p.setBrush(Qt::green);					
						p.drawEllipse(initMethod_PickThreePoints.disk, PointRadius, PointRadius);
						p.setPen(Qt::green);
						p.drawText(initMethod_PickThreePoints.disk, "disk");
						p.setPen(Qt::white);

						if (stepsPerformed > 2)
						{
							QColor lightBlue = QColor(Qt::blue).lighter();
							p.setBrush(lightBlue);
							p.drawEllipse(initMethod_PickThreePoints.cord, PointRadius, PointRadius);
							p.setPen(lightBlue);
							p.drawText(initMethod_PickThreePoints.cord, "spinal cord");
							p.setPen(Qt::white);
						}
					}
				}

				break;
			}

			default:;
		}
				
	slice->setPixmap(newSlice);
	}
}

void slicePainter::mousePressEvent(QMouseEvent *event)
{
	switch (initializationMethod)
	{
		case DrawVertebraOutline:
		{
			if (slices.size()>0 && slice->geometry().contains(event->pos()))
			{
				delete path;
				path = new QPainterPath(slice->mapFromGlobal(event->globalPos()));

				stepsPerformed = 1;
				slice->setCursor(Qt::CrossCursor);
				startButton->setEnabled(true);				
				paintSlice();
			}
			else
				QApplication::beep();

			break;
		}

        case PlaceSinglePoint:
        {
            //do nothing (-> point is placed on mouseRelease / moved on mouseMove)
			break;
        }

		case PlaceCenterAndTwoBP:
		{
			if (slices.size()>0 && stepsPerformed > 0)
			{
				initMethod_PlaceCenterAndTwoBP.movingPoint = -1;

				QPointF mousePos = slice->mapFromGlobal(event->globalPos());

				QPointF * pts[] = {&initMethod_PlaceCenterAndTwoBP.center,
								   &initMethod_PlaceCenterAndTwoBP.boundaryPoint1,
								   &initMethod_PlaceCenterAndTwoBP.boundaryPoint2};

				QPointF distVec;
				qreal   dist;

				for (unsigned int i = 0; i < 3; ++i)
				{
					distVec = (mousePos - *pts[i]);
					dist	= sqrt(distVec.x()*distVec.x() + distVec.y()*distVec.y());

					if (dist <= PointRadius)
					{
						initMethod_PlaceCenterAndTwoBP.movingPoint = i;
						break;
					}
				}
			}

			break;
		}

		case PickVertebraCenters:
		{
			if (slices.size()>0 && stepsPerformed > 0)
			{
				initMethod_PickVertebraCenters.movingPoint = -1;

				QPointF mousePos = slice->mapFromGlobal(event->globalPos());
				
				QPointF distVec;
				qreal   dist;

				unsigned int n = initMethod_PickVertebraCenters.centers.size();
				for (unsigned int i = 0; i < n; ++i)
				{
					distVec = (mousePos - initMethod_PickVertebraCenters.centers.at(i));
					dist	= sqrt(distVec.x()*distVec.x() + distVec.y()*distVec.y());

					if (dist <= PointRadius)
					{
						initMethod_PickVertebraCenters.movingPoint = i;
						break;
					}
				}
			}

			break;
		}

		case PickThreePoints:
		{
			if (slices.size()>0 && stepsPerformed == 3)
			{
				initMethod_PickThreePoints.movingPoint = -1;

				QPointF mousePos = slice->mapFromGlobal(event->globalPos());

				//scale circle
				if (event->button() == Qt::RightButton)
				{					
					initMethod_PickThreePoints.movingPoint = 3;
					lastMousePosOnSlice = mousePos;
				}
				//move points
				else
				{
					QPointF * pts[] = {&initMethod_PickThreePoints.vertebraCenter,
									   &initMethod_PickThreePoints.disk,
									   &initMethod_PickThreePoints.cord};

					QPointF distVec;
					qreal   dist;

					for (unsigned int i = 0; i < 3; ++i)
					{
						distVec = (mousePos - *pts[i]);
						dist	= sqrt(distVec.x()*distVec.x() + distVec.y()*distVec.y());

						if (dist <= PointRadius)
						{
							initMethod_PickThreePoints.movingPoint = i;
							break;
						}
					}
				}
			}

			break;
		}

		default:;
	}
}

void slicePainter::mouseReleaseEvent(QMouseEvent * event)
{
	switch (initializationMethod)
	{
		case DrawVertebraOutline:
		{	
			if (slices.size()>0)
			{
				path->closeSubpath();
				path->setFillRule(Qt::WindingFill);
				paintSlice();

				slice->setCursor(Qt::ArrowCursor);
			}
			
			break;
		}

        case PlaceSinglePoint:
        {
            if (slices.size()>0 && slice->geometry().contains(event->pos()))
			{
				initMethod_PlaceSinglePoint.point = slice->mapFromGlobal(event->globalPos());

				stepsPerformed = 1;
				startButton->setEnabled(true);
				paintSlice();
			}

            break;
        }

		case PlaceCenterAndTwoBP:
		{
			if (slices.size()>0 && stepsPerformed == 0 && slice->geometry().contains(event->pos()))
			{	
				initMethod_PlaceCenterAndTwoBP.center = slice->mapFromGlobal(event->globalPos());
				initMethod_PlaceCenterAndTwoBP.boundaryPoint1 = initMethod_PlaceCenterAndTwoBP.center + QPoint(32,  0);
				initMethod_PlaceCenterAndTwoBP.boundaryPoint2 = initMethod_PlaceCenterAndTwoBP.center + QPoint(0, -23);
				
				delete path;
				path = initMethod_PlaceCenterAndTwoBP.createBoundaryPath();
				
				++stepsPerformed;
				startButton->setEnabled(true);

				emit showStatusMessage("Adjust the control points and start the segmentation if the result is satisfying (Step 2 / 2).");

				paintSlice();
			}

			break;
		}

		case PickVertebraCenters:
		{	
			if (initMethod_PickVertebraCenters.movingPoint == -1 && slices.size()>0 && slice->geometry().contains(event->pos()))
			{
				initMethod_PickVertebraCenters.centers.push_back(slice->mapFromGlobal(event->globalPos()));

				stepsPerformed = 1;
				startButton->setEnabled(true);
				paintSlice();
			}

			break;
	}

		case PickThreePoints:
		{
			if (slices.size()>0 && slice->geometry().contains(event->pos()) && stepsPerformed < 3)
			{
				if (stepsPerformed == 0)
				{
					initMethod_PickThreePoints.vertebraCenter = slice->mapFromGlobal(event->globalPos());
					emit showStatusMessage("Please place a point inside a disk (Step 2 / 4) ...");					
				}
				else if (stepsPerformed == 1)
				{
					initMethod_PickThreePoints.disk = slice->mapFromGlobal(event->globalPos());
					emit showStatusMessage("Please place a point inside the spinal cord (Step 3 / 4) ...");
				}
				else
				{
					assert(stepsPerformed == 2);
					initMethod_PickThreePoints.cord = slice->mapFromGlobal(event->globalPos());
					emit showStatusMessage("Adjust the points (LMB) and the circle size (RMB) and start the segmentation if the result is satisfying (Step 4 / 4).");
					startButton->setEnabled(true);					
				}

				++stepsPerformed;
				paintSlice();
			}
			else
			{
				initMethod_PickThreePoints.movingPoint = -1;
			}

			break;
		}

		default:;
	}
}

void slicePainter::mouseMoveEvent(QMouseEvent *event)
{
	switch (initializationMethod)
	{
		case DrawVertebraOutline:
		{
			if (stepsPerformed > 0 && slice->geometry().contains(event->pos()))
			{
				path->lineTo(slice->mapFromGlobal(event->globalPos()));
		        paintSlice();
	        }

			break;
		}

        case PlaceSinglePoint:
        {
            if (slice->geometry().contains(event->pos()) && (event->buttons() != Qt::NoButton))
			{
                QPointF mousePos = slice->mapFromGlobal(event->globalPos());
			    initMethod_PlaceSinglePoint.point = mousePos;
				
			    stepsPerformed = 1;
			    startButton->setEnabled(true);
			    paintSlice();
            }

            break;
        }

		case PlaceCenterAndTwoBP:
		{
			if (stepsPerformed > 0 && slice->geometry().contains(event->pos()))
			{
				if (initMethod_PlaceCenterAndTwoBP.movingPoint >= 0)
				{
					switch (initMethod_PlaceCenterAndTwoBP.movingPoint)
					{
						case 0:
						{
							QPointF newCenter = slice->mapFromGlobal(event->globalPos());
							QPointF diffVec   = newCenter - initMethod_PlaceCenterAndTwoBP.center;

							//when moving the center, move the other two points along
							initMethod_PlaceCenterAndTwoBP.center		  += diffVec;
							initMethod_PlaceCenterAndTwoBP.boundaryPoint1 += diffVec;
							initMethod_PlaceCenterAndTwoBP.boundaryPoint2 += diffVec;

							break;
						}

						case 1:
						{
							initMethod_PlaceCenterAndTwoBP.boundaryPoint1 = slice->mapFromGlobal(event->globalPos());
							break;
						}

						case 2:
						{
							initMethod_PlaceCenterAndTwoBP.boundaryPoint2 = slice->mapFromGlobal(event->globalPos());							
							break;
						}

						default:;
					}

					delete path;
					path = initMethod_PlaceCenterAndTwoBP.createBoundaryPath();

					paintSlice();
				}
			}

			break;
		}

		case PickVertebraCenters:
		{
			if (stepsPerformed > 0 && slice->geometry().contains(event->pos()))
			{
				if (initMethod_PickVertebraCenters.movingPoint >= 0)
				{					
					QPointF mousePos = slice->mapFromGlobal(event->globalPos());
					initMethod_PickVertebraCenters.centers.at(initMethod_PickVertebraCenters.movingPoint) = mousePos;
					
					paintSlice();
				}				
			}

			break;
		}

		case PickThreePoints:
		{
			if (stepsPerformed == 3 && slice->geometry().contains(event->pos()))
			{
				if (initMethod_PickThreePoints.movingPoint >= 0)
				{					
					switch (initMethod_PickThreePoints.movingPoint)
					{
						case 0:
						{
							initMethod_PickThreePoints.vertebraCenter = slice->mapFromGlobal(event->globalPos());
							break;
						}

						case 1:
						{
							initMethod_PickThreePoints.disk = slice->mapFromGlobal(event->globalPos());
							break;
						}

						case 2:
						{
							initMethod_PickThreePoints.cord = slice->mapFromGlobal(event->globalPos());
							break;
						}
						//scaling circle
						case 3:
						{
							QPointF mousePos = slice->mapFromGlobal(event->globalPos());
							QPointF deltaVec = mousePos - lastMousePosOnSlice;

							initMethod_PickThreePoints.circleRadius += deltaVec.y();

							lastMousePosOnSlice = mousePos;
							break;
						}

						default:;
					}
					
					paintSlice();
				}
			}

			break;
		}

		default:;
	}	
}

QPainterPath * slicePainter::InitMethod_PlaceCenterAndTwoBP::createBoundaryPath()
{
	//if desired, this could be a property given by the user
	qreal pincushionFactor = 0.2;

	QPointF v1 = boundaryPoint1 - center;
	QPointF v2 = boundaryPoint2 - center;

	QPointF p1 = center - v1 - v2,
			p2 = center + v1 - v2,
			p3 = center + v1 + v2,
			p4 = center - v1 + v2;
				 
	qreal r = 1.0 - pincushionFactor;

	QPointF q1 = center - r * v2,
			q2 = center + r * v1,
			q3 = center + r * v2,
			q4 = center - r * v1;
		
	QPainterPath * painterPath = new QPainterPath(p1);
		
	painterPath->quadTo(q1, p2);
	painterPath->quadTo(q2, p3);
	painterPath->quadTo(q3, p4);
	painterPath->quadTo(q4, p1);
	
	return painterPath;
}