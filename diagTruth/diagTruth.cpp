#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>
#include "lp_lib.h" //http://lpsolve.sourceforge.net/5.5/

#include "obj.hh"
#include "vec3.h"
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

#include <itkImageFileReader.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryNotImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::OrImageFilter<VisualizingImageType> orType;
typedef itk::BinaryBallStructuringElement
    <VisualizingImageType::PixelType, 3> StructuringElementType;
//VisualizingImageType::Pointer centers=0, fg=0, bg=0, roi=0, maskPW=0, maskGC=0, seeds=0, tempBG=0;


string base, folderManual, autoBase, filename;
const int bufsize=99999;
char buf[bufsize+1], *tmp;
int seedID=2; //first V.B. seed id =2, backgroundID=1, 0 rest (to be segmented)
int ind1=0, ind2=0;

#define setImpreciseTolerances(filter)\
filter->SetCoordinateTolerance(1e-2);\
filter->SetDirectionTolerance(1e-2);

using namespace Eigen;

const unsigned degree=3;
struct Vertebra
{
    vec3 center;
    double volume, radius, surface;
    int labelIndex;
    static const char* labels[];
    Cell *qe;
    int crushedFlag, spondylFlag;
    double crushed, spondyl;

    double Vertebra::calcCurrentRadius() //average distance of vertices to center
    {
        CellVertexIterator cvi(qe);
        Vertex *v;
        radius=0;
        while ( (v=cvi.next())!=0 )
            radius+=(center-v->pos).length();
        radius/=qe->countVertices();
        return radius;
    }

    //calculates center, volume and surface area
    void Vertebra::calcNewCenter()
    {
        double tv;
        volume=0;
        surface=0;
        vec3 newC(0,0,0), tc, a,b,c;
        CellFaceIterator cfi(qe);
        Face *f;
        while ( (f=cfi.next()) != 0 )
        {
            FaceEdgeIterator fei(f);
            Edge *e=fei.next();
            a=e->Org()->pos;
            e=fei.next();
            b=e->Org()->pos;
            e=fei.next();
            c=e->Org()->pos;
            tv=(a-center) * ((b-center)^(c-center)) / 6.0;
            //volume of tetrahedron consisting of this face and the old center
            tc=(a+b+c+center)/4.0; //center of mass of this tetrahedron
            newC+=(tc*tv);
            volume+=tv;
            surface+=((b-a)^(c-a)).length()/2; //vector product is double surface area
        }
        center=newC/volume;
    }
};

const char* Vertebra::labels[]=
    {"S5", "S4", "S3", "S2", "S1", "L5", "L4", "L3", "L2", "L1",
    "T12", "T11", "T10", "T9", "T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1",
    "C7", "C6", "C5", "C4", "C3", "C2", "C1"};

vector<Vertebra *> vertebra;

//disease thresholds:
double crushedVertebraThreshold=0.2;
double significantSpondylolisthesis=0.25;

using namespace std;

double polyEval(const Eigen::VectorXd P, const double x) //evaluate using Horner's method
{
    double result=0;
    for (int i=P.size()-1; i>=0; i--)
        result=result*x+P[i];
    return result;
}

double dEval(const Eigen::VectorXd P, const double x) //evaluate polynomial's derivative
{
    double result=0;
    for (int i=P.size()-1; i>0; i--)
        result=result*x+i*P[i];
    return result;
}

Eigen::VectorXd derive(const Eigen::VectorXd P) //symbolic derivative of a polynomial
{
    Eigen::VectorXd result=Eigen::VectorXd::Zero(P.size()-1);
    for (int i=0; i<result.size(); i++)
        result[i]=P[i+1]*(i+1);
    return result;
}

void polynomialFitL1Norm(unsigned degree, std::vector<double> x, std::vector<double> y, VectorXd &p)
{
//2D L-1 norm minimization formulated as a linear program:
//
//minimize    e0+e1+...+en
//
//subject to  p0+x0*p1+x0^2*p2+x0^3*p3+e0>=y0
//            p0+x0*p1+x0^2*p2+x0^3*p3-e0<=y0
//            ...
//            p0+xn*p1+xn^2*p2+xn^3*p3+en>=yn
//            p0+xn*p1+xn^2*p2+xn^3*p3-en<=yn

    lprec *lp=make_lp(0, degree+1+y.size());
    
    char buffer[3]="pX";
    for (int k=0; k<=degree; k++)
    {
        buffer[1]='0'+k;
        set_col_name(lp, 1+k, buffer);
        set_unbounded(lp, 1+k);
    }
    buffer[0]='e';
    for (int k=0; k<y.size(); k++)
    {
        if (k<10)
            buffer[1]='0'+k;
        else
            buffer[1]='A'+k-10; //hexadecimal digits
        set_col_name(lp, 1+degree+1+k, buffer);
    }
    
    set_add_rowmode(lp, TRUE);
    std::vector<REAL> row(degree+1+y.size()+1);
    for (int i=0; i<y.size(); i++)
    {
        std::fill(row.begin(), row.end(), 0);
        double zk=1;
        for (int k=0; k<=degree; k++)
        {
            row[1+k]=zk;
            zk*=x[i];
        }

        row[1+degree+1+i]=1;
        add_constraint(lp,&row[0],GE,y[i]);

        row[1+degree+1+i]=-1;
        add_constraint(lp,&row[0],LE,y[i]);
    }
    set_add_rowmode(lp, FALSE);

    std::fill(row.begin(), row.end(), 0);
    for (int k=0; k<y.size(); k++)
        row[1+degree+1+k]=1;
    set_obj_fn(lp,&row[0]);
    
    set_minim(lp);
    write_lp(lp, "p.lp"); //debug
    set_verbose(lp, 1); //only critical errors
    solve(lp);

    p.resize(degree+1);
    get_variables(lp, &row[0]);
    for (int k=0; k<=degree; k++)
        p[k]=row[k];
    delete_lp(lp);
}

void polynomialFitL1Norm(unsigned degree, std::vector<double> volumes, VectorXd &p)
{
    std::vector<double> x(volumes.size());
    for (int i=0; i<x.size(); i++)
        x[i]=i;
    polynomialFitL1Norm(degree, x, volumes, p);
}

void polynomialFitL1Norm(unsigned degree, std::vector<Vertebra *> vertebra, VectorXd &xpoly, VectorXd &ypoly)
{
    std::vector<double> x(vertebra.size());
    for (int i=0; i<x.size(); i++)
        x[i]=vertebra[i]->center[2];
    
    std::vector<double> y(vertebra.size());
    for (int i=0; i<y.size(); i++)
        y[i]=vertebra[i]->center[0];
    polynomialFitL1Norm(degree, x, y, xpoly); //x-axis polynomial
    
    for (int i=0; i<y.size(); i++)
        y[i]=vertebra[i]->center[1];
    polynomialFitL1Norm(degree, x, y, ypoly); //y-axis polynomial
}

unsigned binomial(unsigned n, unsigned k)
{
    if (k==0 || k==n)
        return 1;
    else
        return binomial(n-1,k-1)+binomial(n-1,k);
}

VectorXd translatePolynomial(const VectorXd P, double minusTz)
{
    VectorXd r(P.size()); //allocate result
    for (int i=0; i<P.size(); i++)
    {
        r[i]=0;
        for (int k=P.size()-1; k>=i; k--)
            r[i]=r[i]*minusTz+binomial(k,i)*P[k];
    }
    return r;
}

void polynomialFit(unsigned degree, const std::vector<Vertebra *> vertebra, Eigen::VectorXd &x, Eigen::VectorXd &y, bool ignoreOrientations)
//x,y - output polynomial coefficients for A0+A1*z+A2*z^2+...+An*Z^n
{
    vec3 cTrans;
    for (int i=0; i<vertebra.size(); i++)
        cTrans+=vertebra[i]->center;
    cTrans/=vertebra.size();
    for (int i=0; i<vertebra.size(); i++)
        vertebra[i]->center-=cTrans; //centroid translation, numerically stable

    polynomialFitL1Norm(degree, vertebra, x, y);
    //if (ignoreOrientations)
    //    polynomialFitLinear(degree, vertebra, x, y);
    //else
    //    polynomialFitNonLinear(degree, vertebra, x, y);

    for (int i=0; i<vertebra.size(); i++)
        vertebra[i]->center+=cTrans; //counter centroid translation
    x=translatePolynomial(x,-cTrans[2]);
    x[0]+=cTrans[0];
    y=translatePolynomial(y,-cTrans[2]);
    y[0]+=cTrans[1];
}

void TheilSenEstimator(vector<double> y, double &slope, double &yIntercept)
{
    if (y.empty())
        return;
    if (y.size()==1)
    {
        slope=0;
        yIntercept=y[0];
        return;
    }

	vector<double> m;
	for (int i=0; i<y.size()-1; i++)
		for (int k=i+1; k<y.size(); k++)
			m.push_back((y[k]-y[i])/(k-i)); //slopes
	std::sort(m.begin(), m.end());
	slope=m[m.size()/2]; //median of slopes
		
	m.clear();
	for (int i=0; i<y.size(); i++)
		m.push_back(y[i]-slope*i); //y-intercepts
	std::sort(m.begin(), m.end());
	yIntercept=m[m.size()/2]; //median of y-intercepts
}

double distanceToPoly(vec3 p, Eigen::VectorXd x, Eigen::VectorXd y, vec3& closestPoint)
{
    closestPoint[0]=polyEval(x, p[2]);
    closestPoint[1]=polyEval(y, p[2]);
    closestPoint[2]=p[2];
    vec3 pa=closestPoint-p;
    double dot, dAlen, paLen;
    do
    {
        vec3 dA(dEval(x, closestPoint[2]), dEval(y, closestPoint[2]), 1);
        dot=dA*pa;
        dAlen=dA.length();
        paLen=pa.length();
        closestPoint-=dA*(dot/dAlen);
        pa=closestPoint-p;
    } while (pa.length()<paLen && abs(dot)/(dAlen*paLen)>0.1);
    //process converges and cos>0.1
    return pa.length();
}

void diagnose(const std::vector<Vertebra *> vertebra, const Eigen::VectorXd x, const Eigen::VectorXd y, std::string filename_base)
{
    ofstream diag((filename_base+".txt").c_str());
    double lowerZ=vertebra[0]->center[2];
    double upperZ=vertebra[vertebra.size()-1]->center[2];
    
    //crushed vertebral bodies, v1
    std::vector<double> vols;
    for (int i=0; i<vertebra.size(); i++)
        vols.push_back(vertebra[i]->volume);
    Eigen::VectorXd volF; //expected volume function
    polynomialFitL1Norm(3, vols, volF);
    for (int i=vertebra.size()-1; i>=0; i--)
    {
        double expectedVolume=polyEval(volF,i);
        vertebra[i]->crushed=100*abs(expectedVolume-vertebra[i]->volume)/(expectedVolume);      

        if (vertebra[i]->volume-(expectedVolume)<-crushedVertebraThreshold*vertebra[i]->volume)
            vertebra[i]->crushedFlag=1; //crushed
        else if (vertebra[i]->volume-(expectedVolume)>crushedVertebraThreshold*vertebra[i]->volume)
            vertebra[i]->crushedFlag=0; //oversegmented
        else
            vertebra[i]->crushedFlag=0; //nothing special
    }

    //spondylolisthesis
    vertebra[vertebra.size()-1]->spondyl=0;
    vertebra[vertebra.size()-1]->spondylFlag=0;
    for (int i=vertebra.size()-2; i>=0; i--)
    {
        vec3 pp1,pp2; //poly point
        double dm1=distanceToPoly(vertebra[i]->center,x,y,pp1);
        distanceToPoly(vertebra[i+1]->center,x,y,pp2);
        vec3 d1=pp1-vertebra[i]->center, d2=vertebra[i+1]->center-pp2;
        vec3 cd=d1+d2;

        if (cd.length()<2*vertebra[i]->radius*significantSpondylolisthesis)
            vertebra[i]->spondylFlag=0;
        else
            vertebra[i]->spondylFlag=1;

        vertebra[i]->spondyl=50*cd.length()/vertebra[i]->radius;
    }

    //scoliosis, cobb angle >20° warrents further tracking, >30° warrents action
    int maxi=0,maxk=0;
    double maxCobb=0;
    for (int i=1; i<vertebra.size()-2; i++) //skip first and last
    {
        vec3 vi=vec3(dEval(x,vertebra[i]->center[2]),0,1); //projection to XZ plane
        for (int k=i+1; k<vertebra.size()-1; k++)
        {
            vec3 vk=vec3(dEval(x,vertebra[k]->center[2]),0,1); //projection to XZ plane
            double cobb=(180/M_PI)*acos(vi*vk/(vi.length()*vk.length()));
            if (cobb>=maxCobb)
            {
                maxCobb=cobb;
                maxi=i;
                maxk=k;
            }
        }
    }
    
    diag.precision(1);
    diag<<std::fixed;
    diag<<maxCobb<<" "<<vertebra[maxk]->labelIndex<<' ';
    diag<<vertebra[maxi]->labelIndex<<"  CobbAngle° startIndex endIndex";

    //output per-vertebra info
    for (int i=vertebra.size()-1; i>=0; i--)
    {
        diag<<endl<<vertebra[i]->volume<<" ";
        diag.precision(4);
        diag<<vertebra[i]->center[0]<<" "<<vertebra[i]->center[1]<<" "<<vertebra[i]->center[2]<<" ";
        diag.precision(1);
        diag<<vertebra[i]->crushed<<" "<<vertebra[i]->spondyl<<" ";
        diag<<vertebra[i]->crushedFlag<<" "<<vertebra[i]->spondylFlag;
        diag<<"  ("<<Vertebra::labels[vertebra[i]->labelIndex]<<") volume center crushed% spondyl% crushedFlag spondylFlag";
    }

    diag.close();
}

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

vec3 VBcenter(VisualizingImageType::Pointer ref)
{
    itk::ImageRegionConstIteratorWithIndex<VisualizingImageType> it(ref, ref->GetLargestPossibleRegion());
    it.GoToBegin();
    itk::ContinuousIndex<double,3> ind;
    unsigned long long i=0,j=0,k=0, sumCount=0;
    while (!it.IsAtEnd())
    {
        if (it.Value()!=0)
        {
            VisualizingImageType::PixelType val=it.Value();
            ind=it.GetIndex();
            i+=ind[0]*val;
            j+=ind[1]*val;
            k+=ind[2]*val;
            sumCount+=val;
        }
        ++it;
    }
    ind[0]=(long double)i/sumCount;
    ind[1]=(long double)j/sumCount;
    ind[2]=(long double)k/sumCount;
    VisualizingImageType::PointType p;
    ref->TransformContinuousIndexToPhysicalPoint(ind,p);
    return vec3(p[0],p[1],p[2]);
}

int main( int argc, char **argv )
{
    if (argc<4)
        return 1; //error

    string v1(argv[2]); //start vertebra
    string v2(argv[3]); //end vertebra
    while (Vertebra::labels[ind1]!=v1)
        ind1++;
    while (Vertebra::labels[ind2]!=v2)
        ind2++;
    if (ind1>ind2)
        swap(ind1, ind2);

    folderManual=string(argv[1]);

    for (int ind=ind1; ind<=ind2; ind++)
    {
        VisualizingImageType::Pointer ref=openExact(folderManual+Vertebra::labels[ind]+".mha");
        vertebra.push_back(new Vertebra);
        vertebra[ind-ind1]->labelIndex=ind;
        vertebra[ind-ind1]->center=VBcenter(ref); //from mask
        vertebra[ind-ind1]->qe=objReadCell((folderManual+Vertebra::labels[ind]+".obj").c_str());
        //objWriteCell(vertebra[ind-ind1]->qe, (folderManual+Vertebra::labels[ind]+".obj").c_str()); //debug
        vertebra[ind-ind1]->calcNewCenter(); //from poly
        vertebra[ind-ind1]->calcCurrentRadius();
    }

    Eigen::VectorXd x,y;
    polynomialFit(4,vertebra, x,y, true);
    diagnose(vertebra, x,y, folderManual+"diag");
    cout<<"Finished "+folderManual;
}