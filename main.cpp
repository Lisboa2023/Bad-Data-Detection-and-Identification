#include<iostream>
#include"ChiSquareTest.h"

int main()
{   
    const int size = 5, df = 4;
    const float cl = 0.99;

    const float measurement[] = {1.0, 1.2, 1.6, 1.0, 0.9};
    const float estimatedMeasurements[] = {0.5, 2.4, 0.8, 0.8, 4.0};
    const float covarianceMatrix[][size] = {{1, 0,0,0,0},
                                             {0,1,0,0,0},
                                             {0,0,1,0,0},
                                             {0,0,0,1,0},
                                             {0,0,0,0,1}};

    const float *cmPtr = covarianceMatrix[0];
    
/* 
    const int size = 3, df = 1;
    const float cl = 0.95;

    const float measurement[] = {-4.70, 2.04, -1.90};
    const float estimatedMeasurements[] = {-4.70, 1.97, -1.97};
    const double covarianceMatrix[][size] = {{0.004, 0,0},
                                             {0,0.002,0},
                                             {0,0,0.002}};

    const double *cmPtr = covarianceMatrix[0];  */                                      

    ChiSquareDistribution X2(size,df,cl);
    X2.ChiSquareTest(measurement, estimatedMeasurements, cmPtr);

    return 0;
}



