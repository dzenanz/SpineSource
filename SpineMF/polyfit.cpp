#include <declarations.h>
#include "vertebra.h"
#include <vtkCylinderSource.h>
#include <vtkSmartPointer.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include "lp_lib.h" //http://lpsolve.sourceforge.net/5.5/
#include <string>
#include <vector>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkSphere.h>
#include <vtkClipPolyData.h>
#include <vtkVolume.h>
#include <vtkClipVolume.h>
#include <vtkImageData.h>
#include <vtkTubeFilter.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>

using namespace Eigen;
extern double crushedVertebraThreshold;
extern double significantSpondylolisthesis;

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

struct fit_functor
{
    enum
    {
        InputsAtCompileTime = Dynamic,//2*(degree+1),
        ValuesAtCompileTime = Dynamic
    };
    typedef double Scalar;
    typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;

    unsigned degree; //instance variable
    std::vector<Vertebra *> vert; //instance variable

    fit_functor(unsigned polynomialDegree, std::vector<Vertebra *> vertebra)
        : degree(polynomialDegree), vert(vertebra) {}//constructor
    
    int inputs() const { return 2*(degree+1); }
    int values() const { return 3*vert.size(); }

    int operator()(const VectorXd &x, VectorXd &fvec) const //function evaluation
    {
        vec3 ti,&vi=vert[0]->axis;
        double zi;
        for (unsigned i=0; i<vert.size(); i++)
        {
            zi=vert[i]->center[2];
            vi=vert[i]->axis;
            ti=vec3(dEval(x.segment(0,degree+1),zi),dEval(x.segment(degree+1,degree+1),zi),1);
            fvec[3*i]=vi*ti-vi.length()*ti.length();

            fvec[3*i+1]=polyEval(x.segment(0,degree+1),zi)-vert[i]->center[0];
            fvec[3*i+2]=polyEval(x.segment(degree+1,degree+1),zi)-vert[i]->center[1];
        }

        return 0;
    }

    //LM with this symbolic Jacobian does not work (result=non-sense)
    int df(const VectorXd &x, MatrixXd &fjac) const //Jacobian evaluation
    {
        vec3 &vi=vert[0]->axis;
        double zi,da,db,viRoot,zik;
        for (unsigned i=0; i<vert.size(); i++)
        {
            zi=vert[i]->center[2];
            da=dEval(x.segment(0,degree+1),zi);
            db=dEval(x.segment(degree+1,degree+1),zi);
            vi=vert[i]->axis;
            viRoot=vi.length()/sqrt(da*da+db*db+1);
            zik=1; //values up to now are cached to avoid recomputation
            for (unsigned k=0; k<degree+1; k++)
            {
                zik*=zi; //zi^k
                fjac(3*i,k)=k*zik*(vi[0]-viRoot*da);
                fjac(3*i,degree+1+k)=k*zik*(vi[1]-viRoot*db);

                fjac(3*i+1,k)=zik;
                fjac(3*i+1,degree+1+k)=0; //dA/db=0
                fjac(3*i+2,k)=0; //dB/da=0
                fjac(3*i+2,degree+1+k)=zik;
            }
        }        
        return 0;
    }   
};

void polynomialFitNonLinear(unsigned degree, std::vector<Vertebra *> vertebra, VectorXd &x, VectorXd &y)
//fit positions and orientations to polyline, numerically unstable
{  
    VectorXd xvec;
    //initial guess: a0=x0, b0=y0; other ai,bi =0/1
    xvec.setConstant(2*(degree+1), 0.1);
    xvec[0]=vertebra[0]->center[0];
    xvec[degree+1]=vertebra[0]->center[1];

    // do the computation
    fit_functor functor(degree, vertebra);
    LevenbergMarquardt<fit_functor> lm(functor);
    //int infoDer = lm.lmder1(xvec);
    LevenbergMarquardt<fit_functor>::Index nfev;
    int infoDif = LevenbergMarquardt<fit_functor>::lmdif1(functor,xvec,&nfev);
    
    x=xvec.segment(0,degree+1);
    y=xvec.segment(degree+1,degree+1);
}

void polynomialFitLinear(unsigned degree, std::vector<Vertebra *> vertebra, VectorXd &x, VectorXd &y)
{
    Eigen::MatrixXd A(vertebra.size(), degree+1);
    for (int i=0; i<vertebra.size(); i++)
    {
        double zik=1;
        for (int k=0; k<=degree; k++)
        {
            A(i,k)=zik;
            zik*=vertebra[i]->center[2]; //z
        }
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::VectorXd b(vertebra.size());
    for (int i=0; i<vertebra.size(); i++)
        b[i]=vertebra[i]->center[0]; //x
    x=svd.solve(b);
    for (int i=0; i<vertebra.size(); i++)
        b[i]=vertebra[i]->center[1]; //y
    y=svd.solve(b);
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
    //write_lp(lp, "p.lp"); //debug
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

vtkSmartPointer<vtkPolyData> centerlineTube(Eigen::VectorXd &x, Eigen::VectorXd &y, double fromZ, double toZ, double radius)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    if (fromZ>toZ)
        std::swap(fromZ,toZ);
    double z=fromZ-radius/2;
    double step=std::min(1.0, radius/10);
    while (z<toZ+radius/2)
    {
        points->InsertNextPoint(polyEval(x,z), polyEval(y,z), z);
        z+=step;
    }

    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
	spline->SetPoints(points);

	vtkSmartPointer<vtkParametricFunctionSource> function = vtkSmartPointer<vtkParametricFunctionSource>::New();
	function->SetParametricFunction(spline);
	function->Update();

	vtkTubeFilter *tube = vtkTubeFilter::New();
	tube->SetInputConnection(function->GetOutputPort());
	tube->SetRadius(radius);
	tube->SetNumberOfSides(12);
	tube->CappingOn();
	tube->Update();

	return tube->GetOutput();
}

vtkActor * centerline(Eigen::VectorXd &x, Eigen::VectorXd &y, double fromZ, double toZ, const vec3 color)
{
    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(centerlineTube(x,y,fromZ,toZ,1));
    vtkActor *actor =vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color[0],color[1],color[2]);
    return actor;
}

double centerlineLength(Eigen::VectorXd x, Eigen::VectorXd y, double lowerZ, double upperZ)
{
    double len=0, step=(upperZ-lowerZ)/10000;
    double xp=polyEval(x, lowerZ), yp=polyEval(y, lowerZ);
    for (double d=lowerZ+step; d<=upperZ; d+=step)
    {
        double xd=polyEval(x, d);
        double yd=polyEval(y, d);
        len+=sqrt((xd-xp)*(xd-xp)+(yd-yp)*(yd-yp)+step*step);
        xp=xd; yp=yd; //save current points into old points
    }
    return len;
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
    ofstream diag((filename_base+".csv").c_str());
    ofstream diagHuman((filename_base+".txt").c_str());
    double lowerZ=vertebra[0]->center[2];
    double upperZ=vertebra[vertebra.size()-1]->center[2];
    diagHuman<<"Centerline length: "<<std::fixed<<centerlineLength(x, y, lowerZ, upperZ)<<" mm"<<endl;
    for (int i=0; i<vertebra.size(); i++)
        vertebra[i]->actor->GetProperty()->SetColor(0.35,0.35,0.35); //dark grey
    double r,g,b; //used later

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
        vertebra[i]->actor->GetProperty()->GetColor(r,g,b);
        diagHuman<<Vertebra::labels[vertebra[i]->labelIndex]<<", V: ";
        diagHuman<<vertebra[i]->volume<<" V-v(i): "<<vertebra[i]->volume-(expectedVolume);
        if (vertebra[i]->volume-(expectedVolume)<-crushedVertebraThreshold*vertebra[i]->volume)
        {
            vertebra[i]->crushedFlag=1; //crushed
            diagHuman<<" Crushed vertebra, "<<100*abs(expectedVolume-vertebra[i]->volume)/(expectedVolume)<<"%";           
            vertebra[i]->actor->GetProperty()->SetColor(r,g,1); //add blue
        }
        else if (vertebra[i]->volume-(expectedVolume)>crushedVertebraThreshold*vertebra[i]->volume)
        {
            vertebra[i]->crushedFlag=0; //oversegmented
            diagHuman<<" Oversegmented vertebra, "<<100*abs(expectedVolume-vertebra[i]->volume)/(expectedVolume)<<"%";           
            vertebra[i]->actor->GetProperty()->SetColor(r,g,0.5); //add blue
        }
        else
            vertebra[i]->crushedFlag=0; //nothing special
        diagHuman<<endl;
    }
    diagHuman<<endl;

    if (vertebra.size()>2)
    {
        //spondylolisthesis
        vertebra[vertebra.size()-1]->spondyl=0;
        vertebra[vertebra.size()-1]->spondylFlag=0;
        for (int i=vertebra.size()-2; i>=0; i--)
        {
            diagHuman<<Vertebra::labels[vertebra[i]->labelIndex]<<":\n";

            vec3 pp1,pp2; //poly point
            double dm1=distanceToPoly(vertebra[i]->center,x,y,pp1);
            distanceToPoly(vertebra[i+1]->center,x,y,pp2);
            diagHuman<<"\tCenter distance from smooth polyline: "<<dm1<<"mm.\n";
            vec3 d1=pp1-vertebra[i]->center, d2=vertebra[i+1]->center-pp2;
            vec3 cd=d1+d2;
            diagHuman<<"\tCompound displacement magnitude (";
            diagHuman<<Vertebra::labels[vertebra[i+1]->labelIndex];
            diagHuman<<"->"<<Vertebra::labels[vertebra[i]->labelIndex];
            diagHuman<<"): "<<cd.length()<<endl;

            //calculate angle
            vec3 axis=vertebra[i]->axis;
            vec3 c=vertebra[i]->center;
            vec3 v=vec3(dEval(x,c[2]),dEval(y,c[2]),1);
            double angle=acos(v*axis/(v.length()*axis.length()));
            diagHuman<<"\tAngle between vertebra axis and polyline at center:"<<180*angle/M_PI<<"°\n";
        
            diagHuman<<"\tDiagnosis:";
            if (cd.length()<2*vertebra[i]->radius*significantSpondylolisthesis)
            {
                vertebra[i]->spondylFlag=0;
                diagHuman<<" Insignificant ";
            }
            else
            {
                vertebra[i]->spondylFlag=1;
                diagHuman<<" Significant ";
                vertebra[i]->actor->GetProperty()->GetColor(r,g,b);
                vertebra[i]->actor->GetProperty()->SetColor(r,1,b); //add green
            }
            vertebra[i]->spondyl=50*cd.length()/vertebra[i]->radius;
            diagHuman<<"spondylolisthesis "<<float(vertebra[i]->spondyl)<<"%";
            diagHuman<<endl;
        }

        //scoliosis, cobb angle >20° warrents further tracking, >30° warrents action
        diagHuman<<"\nScoliosis: ";
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
        if (maxCobb<20)
            diagHuman<<"Insignificant, ";
        else if (maxCobb>=45)
            diagHuman<<"Severe, ";
        else //20-45
            diagHuman<<"Present, ";

        diagHuman<<"Cobb angle="<<maxCobb<<"°";
        diagHuman<<" ("<<Vertebra::labels[vertebra[maxk]->labelIndex]<<',';
        diagHuman<<Vertebra::labels[vertebra[maxi]->labelIndex]<<')';
        //also: kyphosis, lordosis?

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
    }
    diag.close();

  //  //crushed vertebra diagnosis, v2
  //  diagHuman.precision(2);
  //  diagHuman.setf(ios::fixed, ios::floatfield);
  //  diagHuman<<'\n'<<endl;
  //  for (int i=vertebra.size()-1; i>=0; i--)
  //  {
  //      diagHuman<<Vertebra::labels[vertebra[i]->labelIndex]<<", V: ";
  //      diagHuman<<vertebra[i]->volume;
  //      
  //      assert(strcmp(vertebraMeshFilename,"subdivBody.obj")==0); //only works with this model
  //      Vertex *topC=vertebra[i]->qe->vertices[31];
  //      Vertex *bottomC=vertebra[i]->qe->vertices[30];

  //      //maxLevel not preserved, we need to recalculate it
  //      CellVertexIterator cvi(vertebra[i]->qe);
  //      Vertex *v;
  //      int maxLevel=0;
  //      while ((v=cvi.next())!=0)
  //      {
  //          if (v->level>maxLevel)
  //              maxLevel=v->level;
  //      }

  //      double dt=0;
  //      VertexEdgeIterator vei(topC);
		//Edge *e;
		//while ((e = vei.next())!=0)
  //      {
  //          v=stepAlong(e, exp2i(maxLevel - topC->level + 1))->Dest();
  //          dt+=distance3D(v->pos, vertebra[i]->center);
  //      }
  //      dt/=topC->valence;

  //      VertexEdgeIterator vei2(bottomC);
		//double db=0;
		//while ((e = vei2.next())!=0)
  //      {
  //          v=stepAlong(e, exp2i(maxLevel - bottomC->level + 1))->Dest();
  //          db+=distance3D(e->Dest()->pos, vertebra[i]->center);
  //      }
  //      db/=bottomC->valence;

  //      double distTop=distance3D(topC->pos, vertebra[i]->center);
  //      double distBottom=distance3D(bottomC->pos, vertebra[i]->center);
  //      vertebra[i]->actor->GetProperty()->GetColor(r,g,b);
  //      
  //      if (db+dt<(1+crushedVertebraThreshold)*(distTop+distBottom))
  //      {
  //          diagHuman<<" Crushed vertebra,";
  //          vertebra[i]->actor->GetProperty()->SetColor(1,g,b); //add red
  //      }
  //      else if (db+dt>(1+2*crushedVertebraThreshold)*(distTop+distBottom))
  //      {
  //          diagHuman<<" Oversegmented vertebra, ";
  //          vertebra[i]->actor->GetProperty()->SetColor(0.5,g,b); //add red
  //      }
  //      diagHuman<<" outD/inD="<<(db+dt)/(distTop+distBottom);
  //      diagHuman<<" ("<<((db+dt)/(distTop+distBottom)-(1+1.5*crushedVertebraThreshold))*100<<"%)";
  //      diagHuman<<endl;
  //  }

    diagHuman.close();
}
//http://www.umm.edu/patiented/articles/how_scoliosis_diagnosed_000068_6.htm
//The severity of scoliosis and need for treatment is usually determined by two factors:
//
//The extent of the spinal curvature (scoliosis is diagnosed when the curve measures 11 degrees or more)
//The angle of the trunk rotation (ATR)
//
//Both are measured in degrees. These two factors are usually related. For example, a person with a spinal curve of 20 degrees will usually have a trunk rotation (ATR) of 5 degrees. These two measurements, in fact, used to be the cutoff for recommending treatment. However, the great majority of 20-degree curves do not get worse. Patients do not usually need medical attention until the curve reaches 30 degrees, and the ATR is 7 degrees.

//http://www.iscoliosis.com/symptoms.html
//Once suspected, scoliosis is usually confirmed with an x-ray, spinal radiograph, CT scan, MRI or bone scan of the spine. The curve is then measured by the Cobb Method and is discussed in terms of degrees. Generally speaking, a curve is considered significant if it is greater than 25 to 30 degrees. Curves exceeding 45 to 50 degrees are considered severe and often require more aggressive treatment.