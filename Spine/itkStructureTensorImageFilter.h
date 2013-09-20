#ifndef __itkStructureTensorImageFilter_h
#define __itkStructureTensorImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>

namespace itk
{

/** \class StructureTensorImageFilter
 * \brief Constructs a structure tensor for each voxel.
 *
 * \ingroup ImageFilters
 */
template< class TInputImage >
class itkStructureTensorImageFilter:public ImageToImageFilter< TInputImage, Image<SymmetricSecondRankTensor<float, TInputImage::ImageDimension>, TInputImage::ImageDimension> >
{
public:
	/** Standard class typedefs. */
	typedef itkStructureTensorImageFilter             Self;
	typedef ImageToImageFilter< TInputImage, Image<SymmetricSecondRankTensor<float, TInputImage::ImageDimension>, TInputImage::ImageDimension> > Superclass;
	typedef SmartPointer< Self >        Pointer;
	typedef Image<SymmetricSecondRankTensor<float, TInputImage::ImageDimension>, TInputImage::ImageDimension>	OutputImageType;
	typedef CovariantVector<float, TInputImage::ImageDimension> CVector;
	typedef SymmetricSecondRankTensor<float, TInputImage::ImageDimension> TensorType;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(itkStructureTensorImageFilter, ImageToImageFilter);
  
	itkSetMacro(Sigma, unsigned int);
	itkGetConstMacro(Sigma, const unsigned int);

protected:
	itkStructureTensorImageFilter()
	{
		m_Sigma=2;
	}
	~itkStructureTensorImageFilter(){}

	virtual void GenerateData();

	class CovariantVectorToTensorFunctor
	{
	public:
		TensorType operator()( CVector in )
		{
			TensorType result;
			for (unsigned i=0; i<TInputImage::ImageDimension; i++)
				for (unsigned k=i; k<TInputImage::ImageDimension; k++)
					result(i,k)=in[i]*in[k];
			return result;
		}
	};

private:
	itkStructureTensorImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

	unsigned int m_Sigma;
};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStructureTensorImageFilter.txx"
#endif

#endif // __itkStructureTensorImageFilter_h