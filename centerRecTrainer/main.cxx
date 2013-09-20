#include <string>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryNotImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkOrImageFilter.h>
#include <itkOrientImageFilter.h>

const string labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::OrImageFilter<VisualizingImageType> orType;

string fnBase, refBase;
int ind1=0, ind2=0;

VisualizingImageType::Pointer openExact(string filename)
{
    typedef itk::ImageFileReader<VisualizingImageType> ReaderType;
    ReaderType::Pointer reader=ReaderType::New();
    reader->SetFileName(filename);
    reader->Update();
    return reader->GetOutput();
}

VisualizingImageType::Pointer openWithRescale(string filename)
{
    typedef itk::ImageFileReader<InternalImageType> ReaderType;
    ReaderType::Pointer reader=ReaderType::New();
    reader->SetFileName(filename);
    reader->Update();
    typedef itk::RescaleIntensityImageFilter<InternalImageType,VisualizingImageType> scalerType;
    scalerType::Pointer rescale=scalerType::New();
    rescale->SetInput(reader->GetOutput());
    rescale->Update();
    return rescale->GetOutput();
}

void writeImage(InternalImageType::Pointer image, std::string fileName)
{
    typedef itk::ImageFileWriter<InternalImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    writer->Update();
}

void writeImage(VisualizingImageType::Pointer image, std::string fileName)
{
    typedef itk::ImageFileWriter<VisualizingImageType> WriterType;
    WriterType::Pointer writer=WriterType::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(image);
    //writer->SetUseCompression(true);
    writer->Update();
}

void boundingBox(VisualizingImageType::Pointer image,
    VisualizingImageType::IndexType& indLow, VisualizingImageType::IndexType& indHigh)
{
    VisualizingImageType::IndexType ind;
    typedef itk::ImageRegionIteratorWithIndex<VisualizingImageType> itType;
    itType it(image, image->GetLargestPossibleRegion());
    it.GoToReverseBegin();
    indLow=it.GetIndex();
    it.GoToBegin();
    indHigh=it.GetIndex();

    while (!it.IsAtEnd())
    {
        if (it.Value()!=0)
        {
            ind=it.GetIndex();
            if (ind[0]<indLow[0])
                indLow[0]=ind[0];
            if (ind[1]<indLow[1])
                indLow[1]=ind[1];
            if (ind[2]<indLow[2])
                indLow[2]=ind[2];
            if (ind[0]>indHigh[0])
                indHigh[0]=ind[0];
            if (ind[1]>indHigh[1])
                indHigh[1]=ind[1];
            if (ind[2]>indHigh[2])
                indHigh[2]=ind[2];
        }
        ++it;
    }
}

double percentFilled(VisualizingImageType::Pointer image, int z,
    VisualizingImageType::RegionType bb)
{
    typedef itk::ImageRegionIteratorWithIndex<VisualizingImageType> itType;
    itType it(image, bb);
    unsigned long long c0=0, c1=0;

    while (!it.IsAtEnd())
    {
        if (it.Value()==0)
            c0++;
        else
            c1++;
        ++it;
    }
    return double(c1)/(c0+c1);
}

void expandDimension(VisualizingImageType::RegionType& bb, int dim,
    VisualizingImageType::SizeType size, int diff)
{
    //if image is smaller than bounding box, bb will have negative indices!
    int l,h;
    if (bb.GetIndex(dim)-diff/2<0)
        l=0;
    else
        l=bb.GetIndex(dim)-diff/2;

    h=l+bb.GetSize(dim)+diff;
    bb.SetSize(dim, bb.GetSize(dim)+diff);

    if (h<size[dim])
        bb.SetIndex(dim, l);
    else
        bb.SetIndex(dim, size[dim]-bb.GetSize(dim));    
}

void makeSquare(VisualizingImageType::RegionType& boundingBoxSlice,
    VisualizingImageType::SizeType imageSize)
{
    if (boundingBoxSlice.GetSize(0)<boundingBoxSlice.GetSize(1))
        expandDimension(boundingBoxSlice, 0, imageSize,
        boundingBoxSlice.GetSize(1)-boundingBoxSlice.GetSize(0)); //expand x
    else if (boundingBoxSlice.GetSize(0)>boundingBoxSlice.GetSize(1))
        expandDimension(boundingBoxSlice, 1, imageSize,
        boundingBoxSlice.GetSize(0)-boundingBoxSlice.GetSize(1)); //expand y
    //else already square
}

void ensureSagittal(VisualizingImageType::Pointer& image)
{
    if (itk::SpatialOrientationAdapter().FromDirectionCosines(image->GetDirection())
        !=itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL) //reorient to sagittal
    {
        itk::OrientImageFilter<VisualizingImageType, VisualizingImageType>::Pointer orientationFilter =
            itk::OrientImageFilter<VisualizingImageType, VisualizingImageType>::New();

        orientationFilter->UseImageDirectionOn();
        orientationFilter->SetDesiredCoordinateOrientationToSagittal();  
        orientationFilter->SetInput(image);
        orientationFilter->UpdateLargestPossibleRegion();
        image=orientationFilter->GetOutput();
    }

}

int main( int argc, char *argv[] )
{
    if (argc<5)
        return 1; //error

    fnBase=string(argv[1]);
    refBase=string(argv[2]);

    string v1(argv[3]); //start vertebra
    string v2(argv[4]); //end vertebra
    while (labels[ind1]!=v1)
        ind1++;
    while (labels[ind2]!=v2)
        ind2++;
    if (ind1>ind2)
        swap(ind1, ind2);

    VisualizingImageType::Pointer ref, refAll;
    VisualizingImageType::IndexType pos1, pos2;
    refAll=openWithRescale(refBase+labels[ind2]+".mha");
    ensureSagittal(refAll);
    VisualizingImageType::SizeType size=refAll->GetLargestPossibleRegion().GetSize();
    refAll->FillBuffer(0);

    VisualizingImageType::SpacingType sp=refAll->GetSpacing();
    float minSp=std::min(sp[0], sp[1]);
    float mulX=sp[0]/minSp;
    float mulY=sp[1]/minSp;
    size[0]*=mulX;
    size[1]*=mulY;

    ofstream f("info.txt", ios::app);
    for (int vertIndex=ind2; vertIndex>=ind1; vertIndex--)
    {
        ref=openWithRescale(refBase+labels[vertIndex]+".mha");
        ensureSagittal(ref);
        boundingBox(ref, pos1, pos2);
        cout<<"BoundingBox of "<<refBase+labels[vertIndex]+".dcm"<<" is:"<<endl;
        cout<<"\tLower bound: "<<pos1<<endl;
        cout<<"\tUpper bound: "<<pos2<<endl;
        pos1[0]*=mulX;
        pos1[1]*=mulY;
        pos2[0]*=mulX;
        pos2[1]*=mulY;
        for (int z=pos1[2]; z<pos2[2]; z++)
        {
            VisualizingImageType::RegionType bb; //slice of a bounding box
            bb.SetIndex(pos1);
            bb.SetIndex(2, z);
            bb.SetSize(0, pos2[0]-pos1[0]);
            bb.SetSize(1, pos2[1]-pos1[1]);
            bb.SetSize(2, 1);
            //double p=percentFilled(ref, z, bb);

            //if (p>0.40) //filled more than 40%
            {
                //cout<<"Slice "<<z<<" filled "<<p*100<<'%'<<endl;
                makeSquare(bb, size);
                f<<"Slices/"+fnBase<<z<<".png 1 "<<bb.GetIndex(0)<<' '<<bb.GetIndex(1)
                    <<' '<<bb.GetSize(0)<<' '<<bb.GetSize(1)<<endl;
            }
        }

        orType::Pointer orFilter=orType::New();
        orFilter->SetInPlace(true);
        orFilter->SetInput1(refAll);
        orFilter->SetInput2(ref);
        orFilter->Update();
        refAll=orFilter->GetOutput();
    }
    f.close();

    ofstream batch("cropBg.cmd", ios::app);
    f.open("bg.txt", ios::app);
    //generate background
    boundingBox(refAll, pos1, pos2);
    pos1[0]*=mulX;
    pos1[1]*=mulY;
    pos2[0]*=mulX;
    pos2[1]*=mulY;

    for (int z=0; z<pos1[2]; z++)
        f<<"Slices/"+fnBase<<z<<".png"<<endl;
    for (int z=pos2[2]+1; z<size[2]; z++)
        f<<"Slices/"+fnBase<<z<<".png"<<endl;

    for (int z=pos1[2]; z<=pos2[2]; z++)
    {
        if (pos1[0]>=16)
            f<<"Backgrounds/Left_"+fnBase<<z<<".png"<<endl;
        if (size[0]-pos2[0]>=16)
            f<<"Backgrounds/Right_"+fnBase<<z<<".png"<<endl;
        if (size[1]-pos2[1]>=16)
            f<<"Backgrounds/Bottom_"+fnBase<<z<<".png"<<endl;

        if (pos1[0]>=16)
            batch<<"nconvert -crop 0 0 "<<pos1[0]<<' '<<size[1]<<" -out png "
                <<"-o Backgrounds/Left_"+fnBase<<z<<".png Slices/"+fnBase<<z<<".png"<<endl;
        if (size[0]-pos2[0]>=16)
            batch<<"nconvert -crop "<<pos2[0]<<" 0 "<<(size[0]-pos2[0])<<' '<<size[1]<<" -out png "
                <<"-o Backgrounds/Right_"+fnBase<<z<<".png Slices/"+fnBase<<z<<".png"<<endl;
        if (size[1]-pos2[1]>=16)
            batch<<"nconvert -crop "<<pos1[0]<<" "<<pos2[1]<<" "<<(pos2[0]-pos1[0])<<' '<<size[1]-pos2[1]<<" -out png "
                <<"-o Backgrounds/Bottom_"+fnBase<<z<<".png Slices/"+fnBase<<z<<".png"<<endl;
    }
    f.close();
    batch.close();
    return 0;
}
