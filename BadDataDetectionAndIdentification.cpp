#include<iostream>

#include"ChiSquareTest.h"
#include"LargestNormalizedResidualTest.h"
#include"HypothesisTest.h"

int main()
{   
 /*    const int degrees_of_freedom = 4;
    const float condidence_level = 0.99;
    const int number_of_measurements = 5;
    const int number_of_bus = 2;
    const float threshold = 3.0;
    const float n_beta = -2.32, n_maximus = 3.0;
        
    float measurementArray[] = {1.0, 1.2, 1.6, 1.0, 0.9};
    float estimatedArray[] = {0.5, 2.4, 0.8, 0.8, 4.0};
    float covarianceMatrix[][number_of_measurements] = {{1,0,0,0,0},
                                             {0,1,0,0,0},
                                             {0,0,1,0,0},
                                             {0,0,0,1,0},
                                             {0,0,0,0,1}}; */

    
    const int degrees_of_freedom = 4;
    const float condidence_level = 0.99;
    const int number_of_measurements = 4;
    const int number_of_bus = 2;
    const float threshold = 3.0;
    const float n_beta = -2.32, n_maximus = 3.0;
     
    
    float measurementArray[] = {3.90, -4.05, -0.48, 2.04};

    float estimatedArray[] = {3.992, -3.61, -0.374, 2.09};

    float jacobianMatrix[][number_of_bus] = {{-50, -100},
                                             {150,-100},
                                             {-100, 200},
                                             {0, 100}};

    float gainMatrix[][number_of_bus] = {{18125000, -18750000},
                                         {-18750000, 57500000}};

    float covarianceMatrix[][number_of_measurements] = {{0.001,0,0,0},
                                                        {0,0.004,0,0},
                                                        {0,0, 0.001,0},
                                                        {0,0,0, 0.002}};

    float *mPtr = measurementArray;
    float *emPtr = estimatedArray;
    float *jPtr = jacobianMatrix[0];
    float *gPtr = gainMatrix[0];
    float *cmPtr = covarianceMatrix[0];

    ChiSquareDistribution X2(number_of_measurements,degrees_of_freedom,condidence_level);
    X2.ChiSquareTest(mPtr, emPtr, cmPtr);

    HypothesisTest HT(n_beta,n_maximus,number_of_measurements,threshold);
    HT.BadDataIdentification(mPtr,emPtr,jPtr,gPtr,cmPtr,number_of_bus);

    return 0;
}



