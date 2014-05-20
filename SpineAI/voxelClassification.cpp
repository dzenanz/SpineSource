#include "voxelClassification.h"
#include "itkImageRegionIterator.h"
#include <algorithm>

float mean(std::vector<float> p)
{
    float sum=0;
    for (unsigned i=0; i<p.size(); i++)
        sum+=p[i];
    return sum/p.size();
}

float median(std::vector<float> p)
{
    std::sort(p.begin(), p.end());
    if (p.size()%2==1)
        return p[p.size()/2];
    else
        return (p[p.size()/2-1]+p[p.size()/2])/2;
}

float product(std::vector<float> p)
{
    float prod1=1, prod0=1;
    for (unsigned i=0; i<p.size(); i++)
    {
        prod1*=p[i];
        prod0*=1-p[i];
    }
    return prod1/(prod1+prod0);
}

InternalImageType::Pointer combineClassifiers
    (std::vector<InternalImageType::Pointer> classifiers,
    CombineFunction combinationMethod)
{
    InternalImageType::Pointer output=InternalImageType::New();
    output->CopyInformation(classifiers[0]);
    output->SetRegions(output->GetLargestPossibleRegion());
    output->Allocate();
    itk::ImageRegionIterator<InternalImageType> out(output, output->GetLargestPossibleRegion());

    unsigned i=0,nIn=classifiers.size();
    std::vector<float> probabilities(nIn); //individual classifier results
    InternalImageType::PointType p;
    InternalImageType::IndexType ind0, indN;

    while(!out.IsAtEnd())
    {
        ind0=out.GetIndex();
        probabilities[0]=classifiers[0]->GetPixel(ind0);
        output->TransformIndexToPhysicalPoint(ind0, p);

        for (i=1; i<nIn; i++) //collect from input images into std::vector
        {            
            classifiers[i]->TransformPhysicalPointToIndex(p, indN);
            probabilities[i]=classifiers[i]->GetPixel(indN);
        }
        
        out.Set(combinationMethod(probabilities)); //calculate output using chosen method
        ++out;
    }
    return output;
}