#pragma once
#include "declarations.h"
#include <vector>

typedef float (*CombineFunction)(std::vector<float>); //define type of "pointer to function"

float mean(std::vector<float> p);
float median(std::vector<float> p);
float product(std::vector<float> p);

//return image is of the same size as first classifier
InternalImageType::Pointer combineClassifiers
    (std::vector<InternalImageType::Pointer> classifiers,
    CombineFunction combinationMethod);

class EigenvaluesToSurfel //similar to Fernandez 2003
{
public:
    float operator()( itk::CovariantVector<float,3> ev )
    {
        float p1,p2,p3;
        p1=(ev[2]-ev[1])/ev[2];
        p2=(ev[1]-ev[0])/ev[2];
        p3=ev[0]/ev[2];

        if (p1>p2 && p1>p3)
            return p1-p2 + p1-p3;
        else
            return 0;
    }
};

class valueClassifier
{
public:
    float avgValue, stdDev;
    float operator()(float val)
    {
        //normal distribution with peak fixed at 100% and width dependent on sigma
        return exp(-(val-avgValue)*(val-avgValue)/(2*stdDev*stdDev));
    }
};

class lhClassifierFunctor
{
public:
    float maxVal;
    float operator()( float l, float h )
    {
        return (maxVal-(h-l))/maxVal;
    }
};

class lhClassifierFunctor2
{
public:
    float maxVal;
    float operator()( float l, float h, float val )
    {
        if (h-l<maxVal/20) //the difference is too small
            return 1;
        float dl=val-l, dh=h-val; //difference current-low and high-current pixel
        return abs(dh-dl)/maxVal;
        //return 1-(1-abs(dh-dl)/maxVal)*(h-l)/maxVal; //similar to classifier 1
    }
};


class dfEdgeClassifier
{
public:
    float operator()(float dfVal)
    {
        if (dfVal>=10)
            return 0.95; //95%
        if (dfVal<=0)
            return -0.1/(dfVal-2); //5% and lower inside the edge itself
        return dfVal/10; //10%-90%
    }
};