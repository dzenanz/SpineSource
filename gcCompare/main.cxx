#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#ifdef _WINDOWS
#include <io.h> //for mktemp
#include <direct.h> //for _mkdir
#endif
using namespace std;

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryNotImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkBinaryMask3DMeshSource.h>
#include <itkConstantPadImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

#include "itkRelabelComponentImageFilter.h" //local version
#include "itkPGMImageIOFactory.h"

const string labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::OrImageFilter<VisualizingImageType> orType;
typedef itk::BinaryBallStructuringElement
    <VisualizingImageType::PixelType, 3> StructuringElementType;
VisualizingImageType::Pointer centers=0, fg=0, bg=0, roi=0, maskPW=0, maskGC=0, seeds=0, tempBG=0;

typedef itk::Mesh<float,3> MeshType;
typedef MeshType::CellType CellType;
typedef CellType::CellAutoPointer CellAutoPointer;
typedef itk::TriangleCell< CellType > TriangleType;
typedef MeshType::PointsContainerIterator  PointIterator;
typedef MeshType::CellsContainerIterator  CellIterator;

string base, folderManual, autoBase, filename;
const int bufsize=99999, radiusBG=20, moreBGseeds=2;
char buf[bufsize+1], *tmp;
int seedID=2; //first V.B. seed id =2, backgroundID=1, 0 rest (to be segmented)
int ind1=0, ind2=0;

#define setImpreciseTolerances(filter)\
filter->SetCoordinateTolerance(1e-2);\
filter->SetDirectionTolerance(1e-2);

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

void writeMeshObj(MeshType::Pointer mesh, const char * fileName)
{
    if (mesh->GetNumberOfPoints()==0)
        return; //avoid errors
	std::ofstream f(fileName);

	f<<"# "<<mesh->GetNumberOfPoints()<<" vertices\n";
	PointIterator pointIterator = mesh->GetPoints()->Begin();
	PointIterator pointEnd      = mesh->GetPoints()->End();
	while( pointIterator != pointEnd )
	{
		f<<"v "<<pointIterator.Value()[0]<<' '<<pointIterator.Value()[1]<<' '<<pointIterator.Value()[2]<< std::endl;
		++pointIterator;
	}

	f<<"# "<<mesh->GetNumberOfCells()<<" faces\n";	
	CellIterator cellIterator = mesh->GetCells()->Begin();
	CellIterator cellEnd      = mesh->GetCells()->End();
	while( cellIterator != cellEnd )
	{
		CellType * cell = cellIterator.Value();
		if ( cell->GetType() == CellType::TRIANGLE_CELL)
		{
			TriangleType * t=(TriangleType *)cell;
			TriangleType::PointIdIterator pit = t->PointIdsBegin();
			f<<"f "<<(*pit+1);
			pit++;
			f<<' '<<(*pit+1);
			pit++;
			f<<' '<<(*pit+1)<<'\n';
		}
		++cellIterator;
	}

	f.close();
}

float dsc(VisualizingImageType::Pointer image, unsigned char val, string refFilename)
{
    VisualizingImageType::Pointer refImage=openExact(refFilename);
    long long ref=0, seg=0, inter=0;
    itk::ImageRegionConstIterator<VisualizingImageType> it(image, image->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<VisualizingImageType> itRef(refImage, refImage->GetLargestPossibleRegion());
    it.GoToBegin();
    itRef.GoToBegin();
    while (!it.IsAtEnd())
    {
        if (it.Value()==val && itRef.Value()!=0)
            inter++; //intersection
        else if (it.Value()!=val && itRef.Value()!=0)
            ref++; //reference - segmented
        else if (it.Value()==val && itRef.Value()==0)
            seg++; //segmented - reference
        ++it; ++itRef;
    }
    return 2.0*inter/(2.0*inter+seg+ref); //dice coefficent
}

string compare1GCPW(int index, VisualizingImageType::Pointer mask, string prefix)
{
    string gcPoly=autoBase+prefix+labels[index]+".obj";
    stringstream result;

    result<<','<<dsc(mask, 1+index-ind1, folderManual+labels[index]+".mha");

    cout<<"Executing marching cubes"<<endl;
    typedef itk::ConstantPadImageFilter<VisualizingImageType,VisualizingImageType> padType;
    padType::Pointer pf=padType::New();
    pf->SetInput(mask);
    pf->SetConstant(0);
    VisualizingImageType::SizeType pad;
    pad[0]=1; pad[1]=1; pad[2]=1;
    pf->SetPadBound(pad);
    pf->Update();

    typedef itk::BinaryMask3DMeshSource<VisualizingImageType,MeshType> mcType;
    mcType::Pointer mcf=mcType::New();
    mcf->SetInput(pf->GetOutput());
    mcf->SetObjectValue(1+index-ind1);
    mcf->Update();
    writeMeshObj(mcf->GetOutput(),gcPoly.c_str());

    string ref=folderManual+labels[index]+".obj";
    cout<<"Executing Metro"<<endl;
    system(("metro.exe "+ref+" "+gcPoly+" -O -e -f >"+tmp).c_str());
    
    float dmax,dmean,drms;
    ifstream f(tmp);
    int debug_counter=0;
    while (f.good())
    {
        debug_counter++;
        f.getline(buf,bufsize);
        if (strcmp(buf,"distances:"))
            continue; //next line
        f.ignore(9); f>>dmax; f.getline(buf,bufsize);
        f.ignore(9); f>>dmean; f.getline(buf,bufsize);
        f.ignore(9); f>>drms; f.getline(buf,bufsize);
        result<<','<<dmax<<','<<dmean<<','<<drms;
    }

    f.close();
    remove(tmp);
    return result.str();
}

string compare1(int index)
{
    string ref=folderManual+labels[index]+".obj";
	string seg=autoBase+filename+'_'+labels[index]+".obj";
    stringstream result;
    VisualizingImageType::Pointer img;

    try
    {
        img=openExact(autoBase+filename+'_'+labels[index]+".mha");
    }
    catch (...)
    {
        return "";
    }
    result<<','<<dsc(img, 1, folderManual+labels[index]+".mha");

    cout<<"Executing Metro"<<endl;
    system(("metro.exe "+ref+" "+seg+" -O -e -f >"+tmp).c_str());
    
    float dmax,dmean,drms;
    ifstream f(tmp);
    int debug_counter=0;
    while (f.good())
    {
        debug_counter++;
        f.getline(buf,bufsize);
        if (strcmp(buf,"distances:"))
            continue; //next line
        f.ignore(9); f>>dmax; f.getline(buf,bufsize);
        f.ignore(9); f>>dmean; f.getline(buf,bufsize);
        f.ignore(9); f>>drms; f.getline(buf,bufsize);
        result<<','<<dmax<<','<<dmean<<','<<drms;
    }

    f.close();
    remove(tmp);
    return result.str();
}

void center1(int index)
{
    VisualizingImageType::Pointer ref=openExact(folderManual+labels[index]+".mha");
    orType::Pointer orFilter=orType::New();
    orFilter->SetInPlace(true);
    orFilter->SetInput1(roi);
    orFilter->SetInput2(ref);
    setImpreciseTolerances(orFilter)
    orFilter->Update();
    roi=orFilter->GetOutput();

    cout<<"Executing SignedMaurerDistanceMapImageFilter "<<labels[index]<<endl;
    typedef itk::SignedMaurerDistanceMapImageFilter
        <VisualizingImageType,InternalImageType> dfType;
    dfType::Pointer df=dfType::New();
    df->SetInput(ref);
    df->InsideIsPositiveOn();
    df->Update();
    //writeImage(df->GetOutput(), "sdf.mha"); //debug

    typedef itk::MinimumMaximumImageCalculator<InternalImageType> minmaxType;
    minmaxType::Pointer minmax=minmaxType::New();
    minmax->SetImage(df->GetOutput());
    minmax->ComputeMaximum();
    VisualizingImageType::IndexType maxIndex=minmax->GetIndexOfMaximum();
    float max=minmax->GetMaximum();
    max=unsigned(max);

    typedef itk::BinaryThresholdImageFilter<InternalImageType,VisualizingImageType> thType;
    thType::Pointer th=thType::New();
    th->SetInput(df->GetOutput());

    th->SetLowerThreshold(max);
    th->SetInsideValue(seedID++);
    th->Update();
    //writeImage(th->GetOutput(), "th.mha"); //debug

    orType::Pointer orFilter1=orType::New();
    orFilter1->SetInPlace(true);
    orFilter1->SetInput1(centers);
    orFilter1->SetInput2(th->GetOutput());
    setImpreciseTolerances(orFilter1)
    orFilter1->Update();
    centers=orFilter1->GetOutput();
}

void constructBg()
{
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radiusBG);
    structuringElement.CreateStructuringElement();
    typedef itk::BinaryDilateImageFilter
        <VisualizingImageType,VisualizingImageType,StructuringElementType> enlargeType;
    enlargeType::Pointer enlarge=enlargeType::New();
    enlarge->SetInput(roi);
    enlarge->SetKernel(structuringElement);
    enlarge->Update();
    fg=enlarge->GetOutput();

    typedef itk::BinaryNotImageFilter<VisualizingImageType> notType;
    notType::Pointer notFilter=notType::New();
    notFilter->SetInput(fg);
    notFilter->Update();
    typedef itk::BinaryThresholdImageFilter<VisualizingImageType,VisualizingImageType> thType;
    thType::Pointer th=thType::New();
    th->SetInput(notFilter->GetOutput());
    th->SetLowerThreshold(1);
    th->SetInsideValue(1);
    th->Update();
    bg=th->GetOutput();
}

std::vector<VisualizingImageType::IndexType> indices(VisualizingImageType::Pointer image, unsigned char value)
{
    std::vector<VisualizingImageType::IndexType> indices;
    itk::ImageRegionConstIteratorWithIndex<VisualizingImageType> it(image, image->GetLargestPossibleRegion());
    it.GoToBegin();
    while (!it.IsAtEnd())
    {
        if (it.Value()==value)
            indices.push_back(it.GetIndex());
        ++it;
    }
    return indices;
}

VisualizingImageType::IndexType createInit(int index)
{
    VisualizingImageType::IndexType retVal, indC;
    std::vector<VisualizingImageType::IndexType> ind;

    VisualizingImageType::SizeType radiusFG;
    radiusFG.Fill(5); radiusFG.SetElement(2,2);
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radiusFG);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter<VisualizingImageType,VisualizingImageType,
        StructuringElementType> DilateType;
    DilateType::Pointer enlarge=DilateType::New();
    enlarge->SetKernel(structuringElement);
    enlarge->SetInput(seeds);

    seeds->FillBuffer(0);

    for (int i=ind1; i<=ind2; i++)
    {
        ind=indices(centers, 2+i-ind1); //random center within this vertebra
        indC=ind[((long long)rand()*(long long)rand())%ind.size()];
        if (i==index)
            retVal=indC;
        seeds->SetPixel(indC, 2+i-ind1);
        enlarge->SetInput(seeds);
        enlarge->SetDilateValue(2+i-ind1);
        enlarge->Update();
        seeds=enlarge->GetOutput();
        seeds->DisconnectPipeline();
    }

    typedef itk::AndImageFilter<VisualizingImageType> andType;
    andType::Pointer andFilter=andType::New();
    andFilter->SetInput1(seeds);
    andFilter->SetInput2(fg);
    andFilter->SetInPlace(true);
    andFilter->Update();
    seeds=andFilter->GetOutput();
    seeds->DisconnectPipeline();

    ind=indices(bg,1);
    tempBG->FillBuffer(0);
    for (int i=0; i<(ind2-ind1+1)*moreBGseeds; i++)
    {
        indC=ind[((long long)rand()*(long long)rand())%ind.size()];
        tempBG->SetPixel(indC, 1);
    }

    structuringElement.SetRadius(radiusBG);
    structuringElement.CreateStructuringElement();
    enlarge->SetKernel(structuringElement);
    enlarge->SetInput(tempBG);
    enlarge->SetDilateValue(1);
    enlarge->Update();
    tempBG=enlarge->GetOutput();

    andFilter->SetInput1(tempBG);
    andFilter->SetInput2(bg);
    andFilter->SetInPlace(true);
    andFilter->Update();
    tempBG=andFilter->GetOutput();
    tempBG->DisconnectPipeline();

    orType::Pointer orFilter=orType::New();
    orFilter->SetInPlace(true);
    orFilter->SetInput1(seeds);
    orFilter->SetInput2(tempBG);
    orFilter->Update();
    seeds=orFilter->GetOutput();

    return retVal;
}

void allocateImages()
{
    fg=VisualizingImageType::New();
    fg->CopyInformation(centers);
    fg->SetRegions(centers->GetLargestPossibleRegion());
    fg->Allocate();
    tempBG=VisualizingImageType::New();
    tempBG->CopyInformation(centers);
    tempBG->SetRegions(centers->GetLargestPossibleRegion());
    tempBG->Allocate();
    roi=VisualizingImageType::New();
    roi->CopyInformation(centers);
    roi->SetRegions(centers->GetLargestPossibleRegion());
    roi->Allocate();
    maskPW=VisualizingImageType::New();
    maskPW->CopyInformation(centers);
    maskPW->SetRegions(centers->GetLargestPossibleRegion());
    maskPW->Allocate();
    maskGC=VisualizingImageType::New();
    maskGC->CopyInformation(centers);
    maskGC->SetRegions(centers->GetLargestPossibleRegion());
    maskGC->Allocate();
    seeds=VisualizingImageType::New();
    seeds->CopyInformation(centers);
    seeds->SetRegions(centers->GetLargestPossibleRegion());
    seeds->Allocate();
}

VisualizingImageType::Pointer compressLabels(VisualizingImageType::Pointer mask)
{
    typedef itk::RelabelComponentImageFilter<VisualizingImageType,VisualizingImageType> chLabType;
    chLabType::Pointer relabel=chLabType::New();
    relabel->SetInput(mask);
    relabel->SetSortByObjectSize(false);
    relabel->Update();
    mask=relabel->GetOutput();
    return mask;
}

int main( int argc, char **argv )
{
    itk::ObjectFactoryBase::RegisterFactory( itk::PGMImageIOFactory::New() );

    if (argc<6)
        return 1; //error

    unsigned randSeed=(unsigned)time( NULL );
    srand( randSeed );
    ofstream rs("randSeed.txt");
    rs<<randSeed;
    rs.close();
    base=string(argv[1]); //original image

    centers=openWithRescale(base);
    writeImage(centers, "base.pgm"); //needed for PW
    allocateImages();
    centers->FillBuffer(0);
    fg->FillBuffer(0);
    roi->FillBuffer(0);
    maskPW->FillBuffer(0);

    string v1(argv[4]); //start vertebra
    string v2(argv[5]); //end vertebra
    while (labels[ind1]!=v1)
        ind1++;
    while (labels[ind2]!=v2)
        ind2++;
    if (ind1>ind2)
        swap(ind1, ind2);

    folderManual=string(argv[2]);
    string manualInit=folderManual.substr(0,folderManual.length()-1)+".init";
    autoBase=string(argv[3]);
    _mkdir(autoBase.c_str());
    int slash=base.find_last_of("/\\");
    int dot=base.find_last_of('.');
    filename=base.substr(slash+1,dot-slash-1);
 
    for (int i=ind1; i<=ind2; i++)
        center1(i);
    constructBg();

    VisualizingImageType::Pointer temp;
    //temp=openImage("D:\\Temp\\clCombinedDS.mha");
    //writeImage(temp, "clCombDS.pgm");


    char fntemplate[11]="tmp_XXXXXX";
    tmp=_mktemp(fntemplate);
    //VisualizingImageType::IndexType indSSC=createInit(vertIndex);
    createInit(4); //index irrelevant and ignored

    //system("../Spine/Release/Spine.exe SpineSeg.init");
    STARTUPINFO siStartupInfo; 
    PROCESS_INFORMATION piProcessInfo; 
    memset(&siStartupInfo, 0, sizeof(siStartupInfo)); 
    memset(&piProcessInfo, 0, sizeof(piProcessInfo)); 
    siStartupInfo.cb = sizeof(siStartupInfo);
    string cmdLine=" "+manualInit+" D:/Temp/";//+" fullAuto"; //fullAuto ignores initialization points
    CreateProcess("../SpineAI/Release/Spine.exe",&cmdLine[0],0,0,false,0,0,"../Spine",&siStartupInfo,&piProcessInfo);
    WaitForSingleObject(piProcessInfo.hProcess,INFINITE);
    //add timing before/after?
        
    system(("move D:\\Temp\\"+filename+"*.* "+autoBase).c_str());
    ofstream fss(autoBase+"SS.csv");
    fss<<"Vertebra,DSC,maxRA,meanRA,rmsRA,maxAR,meanAR,rmsAR"<<endl;       
    for (int i=ind2; i>=ind1; i--)
        fss<<labels[i]<<compare1(i)<<endl;
    fss.close();
    return 0; //we do not need to redo PW and GC segmentations

    writeImage(seeds, "seeds.mha");
    writeImage(seeds, "seeds.pgm");
    try
    {
        cout<<"Executing PowerWatersheds"<<endl;
        system("powerwatsegm.exe -a 2 -i base.pgm -m seeds.pgm -o pw.pgm");
        temp=openExact("pw.pgm");
        if (temp->GetPixelContainer()->Size()!=maskPW->GetPixelContainer()->Size())
        {
            cout<<"Error: base.pgm has different number of voxels from pw.pgm!";
            exit(1);
        }
        memcpy(maskPW->GetBufferPointer(),temp->GetBufferPointer(),maskPW->GetPixelContainer()->Size());
        maskPW=compressLabels(maskPW);
        writeImage(maskPW, "maskPW.mha");
    }
    catch(...)
    {
        cout<<"PowerWatersheds failed"<<endl;
    }

    try
    {
        cout<<"Executing GraphCut"<<endl;
        system("powerwatsegm.exe -a 1 -i base.pgm -m seeds.pgm -o gc.pgm");
        temp=openExact("gc.pgm");
        if (temp->GetPixelContainer()->Size()!=maskGC->GetPixelContainer()->Size())
        {
            cout<<"Error: base.pgm has different number of voxels from gc.pgm!";
            exit(1);
        }
        memcpy(maskGC->GetBufferPointer(),temp->GetBufferPointer(),maskGC->GetPixelContainer()->Size());
        maskGC=compressLabels(maskGC);
        writeImage(maskGC, "maskGC.mha");
    }
    catch(...)
    {
        cout<<"GraphCut failed"<<endl;
    }

    system(("copy base.pgm "+autoBase).c_str());
    //system(("copy clCombDS.pgm "+autoBase).c_str());
    system(("copy randSeed.txt "+autoBase).c_str());
    system(("move seeds.pgm "+autoBase).c_str());
    system(("move seeds.mha "+autoBase).c_str());
    system(("move proba.pgm "+autoBase).c_str());
    system(("move pw.pgm "+autoBase).c_str());
    system(("move maskPW.mha "+autoBase).c_str());
    system(("move gc.pgm "+autoBase).c_str());
    system(("move maskGC.mha "+autoBase).c_str());

    ofstream fpw(autoBase+"PW.csv");
    fpw<<"Vertebra,DSC,maxRA,meanRA,rmsRA,maxAR,meanAR,rmsAR"<<endl;
    for (int i=ind2; i>=ind1; i--)
        fpw<<labels[i]<<compare1GCPW(i,maskPW, "PW_")<<endl;
    fpw.close();

    ofstream fgc(autoBase+"GC.csv");
    fgc<<"Vertebra,DSC,maxRA,meanRA,rmsRA,maxAR,meanAR,rmsAR"<<endl;
    for (int i=ind2; i>=ind1; i--)
        fgc<<labels[i]<<compare1GCPW(i,maskGC, "GC_")<<endl;
    fgc.close();

    return 0;
}
