#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <QSettings>
#include <direct.h> //for _mkdir
#include <nlopt.hpp>
#include <Windows.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <QtConcurrentRun>
#include <QFuture>

using namespace std;

const string labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

const unsigned nParams=16;
typedef itk::Image<unsigned char, 3> VisualizingImageType;

struct dsInfo
{
    string dsName;
    unsigned indL, indU;
    string refPath;
};

vector<dsInfo> ds;

double test(const vector<double> &x, vector<double> &grad, void* f_data)
{
    double r=0;
    for (int i=0; i<x.size(); i++)
        r+=x[i]*x[i];
    return -r; //max f=0 for Xi=0
}

VisualizingImageType::Pointer openExact(string filename)
{
    typedef itk::ImageFileReader<VisualizingImageType> ReaderType;
    ReaderType::Pointer reader=ReaderType::New();
    reader->SetFileName(filename);
    reader->Update();
    return reader->GetOutput();
}

//double distanceError()
//{
//    string ref=folderManual+labels[index]+".obj";
//	string seg=folderAuto2+filename+'_'+labels[index]+".obj";
//    stringstream result;
//    VisualizingImageType::Pointer img;
//
//    try
//    {
//        img=openExact(folderAuto2+filename+'_'+labels[index]+".mha");
//    }
//    catch (...)
//    {
//        return "";
//    }
//    result<<','<<dsc(img, 1, folderManual+labels[index]+".mha");
//
//    cout<<"Executing Metro"<<endl;
//    system(("metro.exe "+ref+" "+seg+" -O -e -f >"+tmp).c_str());
//    
//    float dmax,dmean,drms;
//    ifstream f(tmp);
//    int debug_counter=0;
//    while (f.good())
//    {
//        debug_counter++;
//        f.getline(buf,bufsize);
//        if (strcmp(buf,"distances:"))
//            continue; //next line
//        f.ignore(9); f>>dmax; f.getline(buf,bufsize);
//        f.ignore(9); f>>dmean; f.getline(buf,bufsize);
//        f.ignore(9); f>>drms; f.getline(buf,bufsize);
//        result<<','<<dmax<<','<<dmean<<','<<drms;
//    }
//
//    f.close();
//    remove(tmp);
//    return result.str();
//}

double calcDSC(dsInfo ds1, int index, string tempDir)
{
    VisualizingImageType::Pointer image;
    VisualizingImageType::Pointer refImage=openExact(ds1.refPath+labels[index]+".mha");
    try
    {
        image=openExact(tempDir+ds1.dsName+'_'+labels[index]+".mha");
    }
    catch(...)
    {
        return 0;
    }
    
    long long ref=0, seg=0, inter=0;
    itk::ImageRegionConstIterator<VisualizingImageType> it(image, image->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<VisualizingImageType> itRef(refImage, refImage->GetLargestPossibleRegion());
    it.GoToBegin();
    itRef.GoToBegin();

    while (!it.IsAtEnd())
    {
        if (it.Value()!=0 && itRef.Value()!=0)
            inter++; //intersection
        else if (it.Value()==0 && itRef.Value()!=0)
            ref++; //reference - segmented
        else if (it.Value()!=0 && itRef.Value()==0)
            seg++; //segmented - reference
        ++it; ++itRef;
    }
    return 2.0*inter/(2.0*inter+seg+ref); //dice coefficent
}

double dsc1(dsInfo ds1, const vector<double> &x)
{
    char buf[65536];
    _getcwd(buf,65535);
    string cwd=buf;
    string init="D:/Repo/Segmented/"+ds1.dsName+".init";
    string tempDir=cwd+"/Temp/"+ds1.dsName+'/';
    string clParam="C:/Spine/Release/Spine.exe "+init+' '+tempDir;//+" fullAuto";
	//fullAuto mislabels many vertebrae leading to wrong results

    _mkdir(tempDir.c_str());
    QSettings params((tempDir+"params.ini").c_str(), QSettings::IniFormat);
    params.setValue("percentMoreInflated", x[0]);
    params.setValue("smoothFactor", x[1]);
    params.setValue("bilateralDomainSigma", x[2]);
    params.setValue("bilateralRangeSigma", x[3]);
    params.setValue("lhEps", x[4]);
    params.setValue("structureTensorSigma", x[5]);
    params.setValue("multBinaryThreshold", x[6]);
    params.setValue("cannyVariance", x[7]);
    params.setValue("cannyLowerThreshold", x[8]);
    params.setValue("cannyUpperThreshold", x[9]);
    for (int i=0; i<4; i++)
        params.setValue(QString("classifier")+char('0'+i)+"weight", x[10+i]);
    params.setValue("cifEdges", x[14]);
    params.setValue("sizeGoalForceKappa", x[15]);
    
    params.sync();

    unsigned long exitCode=0;
    STARTUPINFO siStartupInfo; 
    PROCESS_INFORMATION piProcessInfo; 
    memset(&siStartupInfo, 0, sizeof(siStartupInfo)); 
    memset(&piProcessInfo, 0, sizeof(piProcessInfo)); 
    siStartupInfo.cb = sizeof(siStartupInfo);
    CreateProcess(0,&clParam[0],0,0,false,0,0,"../Spine",&siStartupInfo,&piProcessInfo);
    if (WaitForSingleObject(piProcessInfo.hProcess,300000)==WAIT_TIMEOUT) //5 minutes
    {
        TerminateProcess(piProcessInfo.hProcess, 123);
        return 0; //maybe something has been segmented?
    }
    else
    {
        GetExitCodeProcess(piProcessInfo.hProcess, &exitCode);
        if (exitCode)
            return 0; //it finished with an error
    }

    double total=0;
    int count=0;
    for (int i=ds1.indU; i>=ds1.indL; i--)
    {
        double d1=calcDSC(ds1,i,tempDir);
        if (d1>0)
        {
            count++;
            total+=d1;
        }
    }

    if (count!=0)
        total/=count;
    return total;
}

double dsc(const vector<double> &x, vector<double> &grad, void* f_data)
{
    vector<QFuture<double>> r;
    for (int i=0; i<ds.size(); i++)
    {
        QFuture<double> future = QtConcurrent::run<double>(dsc1,ds[i],x);
        r.push_back(future);
        //let the first one run alone, so it computes the decompositions
        if (i==0)
            future.waitForFinished();
    }

    double total=0;
    int count=0;
    for (int i=0; i<r.size(); i++)
        if (r[i].result()>0)
        {
            count++;
            total+=r[i].result();
        }
    if (count!=0)
        total/=count;

	ofstream f("params.txt", ios::app);
    f<<"DSC: "<<total<<"  Params: ";
    for (int i=0; i<x.size(); i++)
        f<<' '<<x[i];
	f<<endl;
    f.close();
    return total; //average dsc
}

int main( int argc, char **argv )
{
    ifstream f("datasets.txt");
    char buf[1000];
    f.getline(buf,999);
    while (!f.eof())
    {
        dsInfo ds1;
        string line=buf;
        unsigned pos=line.find_first_of(' ');
        ds1.dsName=line.substr(0,pos);
        line=line.substr(pos+1);

        pos=line.find_first_of(' ');
        ds1.refPath=line.substr(0,pos);
        line=line.substr(pos+1);

        pos=line.find_first_of(' ');
        ds1.indU=0;
        while (labels[ds1.indU]!=line.substr(0,pos))
            ds1.indU++;
        line=line.substr(pos+1);

        ds1.indL=0;
        while (labels[ds1.indL]!=line.substr(0,pos))
            ds1.indL++;
        
        ds.push_back(ds1);
        f.getline(buf,999);
    }

    //nlopt::opt maximizer(nlopt::GN_DIRECT_L, nParams); //global optimizer
    nlopt::opt maximizer(nlopt::LN_COBYLA, nParams);

    maximizer.set_max_objective(&dsc,0);

    //maximizer.set_maxtime(5); //5 seconds
    //maximizer.set_maxeval(1000);
    maximizer.set_xtol_abs(0.005);
    maximizer.set_stopval(0.95);

    vector<double> lb(nParams);
    vector<double> x(nParams);
    vector<double> ub(nParams);

    lb[0]=0;
    x[0]=0.0833333333333333; //percentMoreInflated
    ub[0]=0.1;

    lb[1]=0;
    x[1]=0.5; //smoothFactor
    ub[1]=1;
    
    lb[2]=1;
    x[2]=5; //bilateralDomainSigma
    ub[2]=25;
    
    lb[3]=0.01;
    x[3]=0.05; //bilateralRangeSigma
    ub[3]=0.25;
    
    lb[4]=0.005;
    x[4]=0.0125; //lhEps
    ub[4]=0.05;
    
    lb[5]=0.5;
    x[5]=1.75; //structureTensorSigma
    ub[5]=3;
    
    lb[6]=0.05;
    x[6]=0.2; //multBinaryThreshold
    ub[6]=0.35;
    
    lb[7]=0.5;
    x[7]=1.25; //cannyVariance
    ub[7]=5;
    
    lb[8]=0.01;
    x[8]=0.05; //cannyLowerThreshold
    ub[8]=0.25;
    
    lb[9]=0.03;
    x[9]=0.265; //cannyUpperThreshold
    ub[9]=0.5;

    for (int i=0; i<4; i++)
    {
        lb[10+i]=0.05;
        x[10+i]=0.2;
        ub[10+i]=1;
    }

    lb[14]=0.1;
    x[14]=8.35; //cifEdges
    ub[14]=10;

    lb[15]=0;
    x[15]=0.1; //sizeGoalForceKappa
    ub[15]=0.5;
    

    maximizer.set_upper_bounds(ub);
    maximizer.set_lower_bounds(lb);

    double fval=0;
    try
    {
        maximizer.optimize(x, fval);
    }
    catch (std::exception exc)
    {
        cout<<"std::exception occured\n"<<exc.what()<<endl;
    }
    catch (...)
    {
        cout<<"Exception occured\n";
    }

    cout<<"Maximum is "<<fval<<", and it is attained for\nx=[";
    for (int i=0; i<x.size(); i++)
        cout<<' '<<x[i];
    cout<<" ]";
}