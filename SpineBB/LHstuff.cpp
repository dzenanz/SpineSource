#include "declarations.h"
#include <QProgressDialog>
#include <QDateTime>
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"

#ifndef round
double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}
#endif

typedef itk::LinearInterpolateImageFunction<InternalImageType, float> ValueInterpolatorType;

float calcLHvalue(InternalImageType::Pointer in, ValueInterpolatorType::Pointer valInterp,
			  InternalImageType::Pointer gm, GradientImageType::Pointer gv,
			  InternalImageType::IndexType voxel, short dir, const float eps)
//in is the image, gm is gradient magnitude, gv is gradient vector
//dir is either +1 or -1 (direction along or opposite to gradient)
{
	typedef itk::Point<float, 3> PointType;
	typedef itk::Vector<float, 3> VectorType;
	typedef itk::CovariantVector<float, 3> CovariantType;
	typedef itk::ContinuousIndex<float, 3> ContinuousIndexType;
	VectorType v;
	CovariantType cv;
	ContinuousIndexType cInd=voxel;
	unsigned iter=0;

	float oldVal, val=in->GetPixel(voxel);
	float gmVal=gm->GetPixel(voxel);
	
	if (gmVal<eps) //voxel not on boundary, do not waste time on any computation
		return val;

	while (iter<100)
	{
		oldVal=val;
		iter++;
		cv=gv->GetPixel(voxel);
		cv.Normalize();
		cv*=dir; //uphill or downhill
		v.SetVnlVector(cv.GetVnlVector());
        cInd+=v;
        if (cInd[0]<=-0.5 || cInd[1]<=-0.5 || cInd[2]<=-0.5 ||
            cInd[0]>=in->GetLargestPossibleRegion().GetSize()[0]-0.5 ||
            cInd[1]>=in->GetLargestPossibleRegion().GetSize()[1]-0.5 ||
            cInd[2]>=in->GetLargestPossibleRegion().GetSize()[2]-0.5) //we are outside of image
		{
			return oldVal;
		}
        
        voxel[0]=round(cInd[0]);
        voxel[1]=round(cInd[1]);
        voxel[2]=round(cInd[2]);
		gmVal=gm->GetPixel(voxel);
		val=valInterp->EvaluateAtContinuousIndex(cInd);
		if (dir>0&&val<=oldVal || dir<0&&val>=oldVal) //exiting local extremum
            return oldVal;
		else if (gmVal<eps) //entering homogeneous region
            return val;
	}
	return val;
}

void calcLHvalues(InternalImageType::Pointer in, const float eps, bool gpu,
				  InternalImageType::Pointer l, InternalImageType::Pointer h)
{
	l->CopyInformation(in);
	l->SetRegions(in->GetLargestPossibleRegion().GetSize());
	l->Allocate();
	h->CopyInformation(in);
	h->SetRegions(in->GetLargestPossibleRegion().GetSize());
	h->Allocate();
	InternalImageType::SizeType size=in->GetLargestPossibleRegion().GetSize();
	InternalImageType::SpacingType sp = in->GetSpacing();
	double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);

	typedef itk::GradientRecursiveGaussianImageFilter<InternalImageType, GradientImageType> GradientFilterType;
	GradientFilterType::Pointer gFilter=GradientFilterType::New();
	gFilter->SetInput(in);
	gFilter->SetSigma(avgSpacing);
    gFilter->SetUseImageDirection(false); //do not work in physical space (work in index space)
	gFilter->Update();
	GradientImageType::Pointer gv=gFilter->GetOutput();

    typedef itk::VectorMagnitudeImageFilter<GradientImageType, InternalImageType> GtMtype;
	GtMtype::Pointer gtmFilter=GtMtype::New();
	gtmFilter->SetInput(gv);
	gtmFilter->Update();
	InternalImageType::Pointer gm=gtmFilter->GetOutput();
    	
	ValueInterpolatorType::Pointer valInterp=ValueInterpolatorType::New();
	valInterp->SetInputImage(in);

	typedef itk::ImageRegionIterator< InternalImageType > IteratorType;
	IteratorType it( in, in->GetLargestPossibleRegion() );
	InternalImageType::IndexType pos;

	for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		pos=it.GetIndex();
		float lval=calcLHvalue(in, valInterp, gm, gv, pos, -1, eps);
		l->SetPixel(pos, lval);
		float hval=calcLHvalue(in, valInterp, gm, gv, pos, +1, eps);
		h->SetPixel(pos, hval);
	}
}

void calc2DJointHistogram(VisualizingImageType::Pointer x_aka_l, VisualizingImageType::Pointer y_aka_h, std::string savefilename)
{
	typedef itk::Image<unsigned char> iucType;
	iucType::Pointer hist8=iucType::New();
	typedef itk::Image<unsigned long long> ullType;
	ullType::Pointer hist=ullType::New();
	ullType::SizeType size;
	size[0]=256; size[1]=256;
	hist->SetRegions(size);
	hist->Allocate();
	hist->FillBuffer(0);

	typedef itk::ImageRegionConstIterator< VisualizingImageType > ConstIteratorType;
	ConstIteratorType l_it( x_aka_l, x_aka_l->GetLargestPossibleRegion() );
	ConstIteratorType h_it( y_aka_h, y_aka_h->GetLargestPossibleRegion() );
	ullType::IndexType index;
	for ( l_it.GoToBegin(), h_it.GoToBegin(); !l_it.IsAtEnd(); ++l_it, ++h_it )
	{
		index[0]=l_it.Get();
		index[1]=255-h_it.Get();
		(*hist)[index]++;
	}

    typedef itk::UnaryFunctorImageFilter<ullType, InternalSliceImageType, Log1Functor> Log1Type;
	Log1Type::Pointer log_filter = Log1Type::New();
	log_filter->SetInput( hist );

	typedef itk::RescaleIntensityImageFilter < InternalSliceImageType, iucType > RescaleImageFilterType;
	RescaleImageFilterType::Pointer rescale = RescaleImageFilterType::New();
	rescale->SetInput( log_filter->GetOutput() );

	typedef itk::ImageFileWriter<iucType> WriterType;
	WriterType::Pointer writer1=WriterType::New();
	writer1->SetFileName(savefilename);
	writer1->SetInput(rescale->GetOutput());
	writer1->Update();
}