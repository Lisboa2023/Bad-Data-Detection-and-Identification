#include<iostream>
#include<cmath>
#include<iomanip>

#include"LargestNormalizedResidualTest.h"

int main(){
    //Largest Normalized Residual Test ========================================= 
   /*  const int size = 5;
    const float threshold = 3.0;

    float measurementArray[] = {1.0285, 1.0121, 0.9893, 0.4796, 0.3891};
    float estimatedMeasurement[] = {0.9999, 0.9886, 0.9833, 0.4856, 0.3821};
    float covarianceMatrix[][size] = {{0.00004637, 0,0,0,0},
                                       {0,0.00003285,0,0,0},
                                       {0,0, 0.000006805,0,0},
                                       {0,0,0, 0.000006805,0},
                                       {0,0,0,0, 0.0000405}};

    float *mPtr = measurementArray;
    float *emPtr = estimatedMeasurement;
    float *cmPtr = covarianceMatrix[0];

    NormalizedResidual LNRTest(size,threshold);
    LNRTest.LargestNormalizedResidualTest(mPtr,emPtr,cmPtr); */
    // =============================================================================

    // Teste 2 ========================================================================
    const int size = 4;
    const float threshold = 3.0;
    const int number_of_bus = 2;

    float measurementArray[] = {3.90, -4.05, -0.48, 2.04};
    float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};
    float covarianceMatrix[][size] = {{0.001, 0,0,0},
                                      {0,0.004,0,0},
                                      {0,0, 0.001,0},
                                      {0,0,0, 0.002}};

    float jacobianMatrix[][number_of_bus] = {{-50, -100},
                                             {150,-100},
                                             {-100, 200},
                                             {0, 100}};

    float gainMatrix[][number_of_bus] = {{18125000, -18750000},
                                         {-18750000, 57500000}};

    float *mPtr = measurementArray;
    float *emPtr = estimatedArray;
    float *cmPtr = covarianceMatrix[0];
    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];

    NormalizedResidual LNRTest(size,threshold);
    LNRTest.LargestNormalizedResidualTest(mPtr,emPtr,jPtr,gPtr,cmPtr,number_of_bus);
    // =================================================================================

    return 0;
}

