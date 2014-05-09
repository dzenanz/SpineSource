#include "declarations.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageFileWriter.h"
#if ITK_VERSION_MAJOR==4
#include "itkVectorMagnitudeImageFilter.h"
#else
#include "itkGradientToMagnitudeImageFilter.h"
#endif

typedef itk::LinearInterpolateImageFunction<InternalImageType, float> ValueInterpolatorType;

#ifndef round
double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}
#endif

//multithreaded calculation of LH values
class lhImageFilter:public itk::ImageToImageFilter< InternalImageType, InternalImageType >
{
public:
    /** Standard class typedefs. */
    typedef lhImageFilter             Self;
    typedef ImageToImageFilter< InternalImageType, InternalImageType > Superclass;
    typedef itk::SmartPointer< Self >        Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(lhImageFilter, ImageToImageFilter);

    itkSetMacro(Epsilon, float);
    itkGetConstMacro(Epsilon, const float);

protected:
    lhImageFilter()
    {
        this->SetNumberOfRequiredOutputs(2);
        this->SetNumberOfRequiredInputs(1);
 
        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );
        m_Epsilon=0.01; //should be set by caller
    }
    ~lhImageFilter(){}

    //dir is either +1 or -1 (direction along or opposite to gradient)
    float calcLHvalue(InternalImageType::IndexType voxel, short dir)
    {
        typedef itk::Point<float, 3> PointType;
        typedef itk::Vector<float, 3> VectorType;
        typedef itk::CovariantVector<float, 3> CovariantType;
        typedef itk::ContinuousIndex<float, 3> ContinuousIndexType;
        VectorType v;
        CovariantType cv;
        ContinuousIndexType cInd=voxel;
        unsigned iter=0;

        float oldVal, val=this->GetInput()->GetPixel(voxel);
        float gmVal=gm->GetPixel(voxel);
    
        if (gmVal<m_Epsilon) //voxel not on boundary, do not waste time on any computation
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
                cInd[0]>=GetInput()->GetLargestPossibleRegion().GetSize()[0]-0.5 ||
                cInd[1]>=GetInput()->GetLargestPossibleRegion().GetSize()[1]-0.5 ||
                cInd[2]>=GetInput()->GetLargestPossibleRegion().GetSize()[2]-0.5) //we are outside of image
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
            else if (gmVal<m_Epsilon) //entering homogeneous region
                return val;
        }
        return val;
    }

    virtual void BeforeThreadedGenerateData()
    {
        InternalImageType::ConstPointer in=this->GetInput();
        this->GetOutput(0)->CopyInformation(in);
        this->GetOutput(0)->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        this->GetOutput(1)->CopyInformation(in);
        this->GetOutput(1)->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        InternalImageType::SpacingType sp = in->GetSpacing();
        double avgSpacing=pow(sp[0]*sp[1]*sp[2], 1/3.0);

        typedef itk::GradientRecursiveGaussianImageFilter<InternalImageType, GradientImageType> GradientFilterType;
        GradientFilterType::Pointer gFilter=GradientFilterType::New();
        gFilter->SetInput(in);
        gFilter->SetSigma(avgSpacing);
        gFilter->SetUseImageDirection(false); //do not work in physical space (work in index space)
        gFilter->Update();
        gv=gFilter->GetOutput();

#if ITK_VERSION_MAJOR==4
        typedef itk::VectorMagnitudeImageFilter<GradientImageType, InternalImageType> GtMtype;
#else
        typedef itk::GradientToMagnitudeImageFilter<GradientImageType, InternalImageType> GtMtype;
#endif
        GtMtype::Pointer gtmFilter=GtMtype::New();
        gtmFilter->SetInput(gv);
        gtmFilter->Update();
        gm=gtmFilter->GetOutput();
        
        valInterp=ValueInterpolatorType::New();
        valInterp->SetInputImage(this->GetInput());
    }

#if ITK_VERSION_MAJOR==4
    virtual void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, itk::ThreadIdType thread)
#else
    virtual void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, int thread)
#endif
    {
        InternalImageType::IndexType pos;
        InternalImageType::ConstPointer input = this->GetInput();
        InternalImageType::Pointer l = this->GetOutput(0);
        InternalImageType::Pointer h = this->GetOutput(1);

        itk::ImageRegionIterator<InternalImageType> outL(l, outputRegionForThread);
        itk::ImageRegionIterator<InternalImageType> outH(h, outputRegionForThread);
        itk::ImageRegionConstIterator<InternalImageType> it(input, outputRegionForThread);
        itk::ProgressReporter progress(this, thread, outputRegionForThread.GetNumberOfPixels());

        while(!outL.IsAtEnd())
        {
            pos=it.GetIndex();
            outL.Set(calcLHvalue(pos, -1));
            outH.Set(calcLHvalue(pos, +1));
            ++it; ++outL; ++outH;
            progress.CompletedPixel();
        }
    }

    itk::DataObject::Pointer MakeOutput(unsigned int idx)
    {
      itk::DataObject::Pointer output;
 
      switch ( idx )
        {
        case 0:
          output = ( InternalImageType::New() ).GetPointer();
          break;
        case 1:
          output = ( InternalImageType::New() ).GetPointer();
          break;
        default:
          std::cerr << "No output " << idx << std::endl;
          output = NULL;
          break;
        }
      return output.GetPointer();
    }
private:
    lhImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

    float m_Epsilon;
    double avgSpacing;
    GradientImageType::Pointer gv;
    InternalImageType::Pointer gm;
    ValueInterpolatorType::Pointer valInterp;
};

void calcLHvalues(InternalImageType::Pointer in, const float eps,
                  InternalImageType::Pointer l, InternalImageType::Pointer h)
{
    lhImageFilter::Pointer lh=lhImageFilter::New();
    lh->SetInput(in);
    lh->SetEpsilon(eps);
    lh->Update();
    l->Graft(lh->GetOutput(0));
    h->Graft(lh->GetOutput(1));
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