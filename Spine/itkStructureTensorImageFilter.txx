#ifndef __itkStructureTensorImageFilter_txx
#define __itkStructureTensorImageFilter_txx

#include <itkGradientImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include "itkStructureTensorImageFilter.h"

namespace itk
{

template< class TInputImage >
void itkStructureTensorImageFilter< TInputImage >
::GenerateData()
{
	typedef Image<CVector, TInputImage::ImageDimension> CovariantImageType;

	typedef GradientImageFilter<TInputImage> GradientType; //for already smooth images
	//typedef GradientRecursiveGaussianImageFilter<TInputImage, CovariantImageType> GradientType;
	typename GradientType::Pointer grad=GradientType::New();
	//grad->SetSigma(m_Sigma);
	grad->SetInput(this->GetInput());
  
	typedef itk::UnaryFunctorImageFilter< CovariantImageType, OutputImageType, CovariantVectorToTensorFunctor > FilterType;
	typename FilterType::Pointer filter=FilterType::New();
	filter->SetInput(grad->GetOutput());

	typedef RecursiveGaussianImageFilter<OutputImageType, OutputImageType> SmoothingType;
	typename SmoothingType::Pointer smoothing=SmoothingType::New();
	smoothing->SetInput(filter->GetOutput());
	smoothing->SetSigma(m_Sigma);
	smoothing->Update();
	this->GraftOutput(smoothing->GetOutput());
}

}// end namespace

#endif //__itkStructureTensorImageFilter_txx