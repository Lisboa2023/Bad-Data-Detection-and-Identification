#include <iostream>
#include <cmath>
#include <iomanip>

#include "LargestNormalizedResidualTest.h"
#include "HypothesisTest.h"

HypothesisTest::HypothesisTest(const float nBeta, const float nMaximus, const int SIZE, const float THRESHOLD):NormalizedResidual(SIZE, THRESHOLD){
    setNBeta(nBeta);
    setNMaximus(nMaximus);
    setNumberSelectedMeasurements(0);
    new_number_suspect_measurements = 0;
}

HypothesisTest::~HypothesisTest(){

    delete [] suspect_selected_measurements;
    delete [] suspect_residual_covariance_matrix;
    delete [] sensitivity_matrix_SS;
    delete [] suspect_residual_measurements;
    delete [] inverse_sensitivity_matrix_ss;
    delete [] estimated_error_measurements;
    delete [] N_measurements;
    delete [] threshold_measurements;

}

//Funcoes SET ============================================================
void HypothesisTest::setNBeta(float nBeta){
    N_beta = nBeta;
}

void HypothesisTest::setNMaximus(float nMaximus){
    N_maximus = nMaximus;
}

void HypothesisTest::setNumberSelectedMeasurements(int numberSelectedMeasurements){
    number_selected_measurements = numberSelectedMeasurements;
}

//Funcoes GET ===========================================================

float HypothesisTest::setNBeta() const{
    return N_beta;
}

float HypothesisTest::setNMaximus() const{
    return N_maximus;
}

int HypothesisTest::getNumberSelectedMeasurements() const{
    return number_selected_measurements;
}

float *HypothesisTest::getSuspectSelectedMeasurements()const{
    return suspect_selected_measurements;
}

float *HypothesisTest::getSuspectResidualCovarianceMatrix()const{
    return suspect_residual_covariance_matrix;
}

float *HypothesisTest::getSensitivityMatrixSS()const{
    return sensitivity_matrix_SS;
}

float *HypothesisTest::getInverseSensitivityMatrixSS()const{
    return inverse_sensitivity_matrix_ss;
}

float *HypothesisTest::getSuspectResidualMeasurements()const{
    return suspect_residual_measurements;
}

//=========================================================================  


//=======================================================================

void HypothesisTest::SelectSuspectMeasurements(){
    
    float num = getNumberOfMeasurements();
    float *temp = getNormalizedMeasurements();
    float normalizedThreshold = getThreshold();

    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){
            number_selected_measurements++;
        }
    }

    suspect_selected_measurements = new float[getNumberSelectedMeasurements()];
    
    int count = 0;
    for(int i = 0; i < num; i++){
        if(temp[i] > normalizedThreshold){   
            suspect_selected_measurements[count] = i;
            count++;
        }
    }
}

void HypothesisTest::SelectSuspectResidualCovarianceMatrix(){
    int num = getNumberOfMeasurements();
    float *temp = getResidualCovarianceMatrix();
    int k;

    suspect_residual_covariance_matrix = new float[getNumberSelectedMeasurements()*getNumberSelectedMeasurements()];

    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        for(int j = 0; j < getNumberSelectedMeasurements(); j++){
            if(i == j){
                k = suspect_selected_measurements[i];
                suspect_residual_covariance_matrix[i*getNumberSelectedMeasurements() + j] = temp[k*num + k];
            }
            else{
                suspect_residual_covariance_matrix[i*getNumberSelectedMeasurements() + j] = 0;
            }
        }
    }
}

void HypothesisTest::SelectSensitivityMatrixSS(){
    int num = getNumberOfMeasurements();
    float *temp = getSensitivityMatrix();
    int k;

    sensitivity_matrix_SS = new float[getNumberSelectedMeasurements()*getNumberSelectedMeasurements()];
    
    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        for(int j = 0; j < getNumberSelectedMeasurements(); j++){
            if(i==j){
                k = suspect_selected_measurements[i];
                sensitivity_matrix_SS[i*getNumberSelectedMeasurements() + j] = temp[k*num + k];
            }
            else{
                sensitivity_matrix_SS[i*getNumberSelectedMeasurements() + j] = 0;
            }

        }
    }
}

void HypothesisTest::CalculateInverseSensitivityMatrixSS(){

    inverse_sensitivity_matrix_ss = new float[getNumberSelectedMeasurements()*getNumberSelectedMeasurements()];

    inverse_sensitivity_matrix_ss = CalculateInverseMatrix(sensitivity_matrix_SS, getNumberSelectedMeasurements());
    
}

void HypothesisTest::SelectSuspectResidualMeasurements(){
    
    float *temp = getResidualMeasurements();
    int k;
    
    suspect_residual_measurements = new float[getNumberSelectedMeasurements()];

    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        k = suspect_selected_measurements[i];
        suspect_residual_measurements[i] = temp[k];
    }
}

void HypothesisTest::CalculateEstimatedErrorMeasurements(){

    estimated_error_measurements = new float[getNumberSelectedMeasurements()];
    
    estimated_error_measurements = MultiplyArray(getInverseSensitivityMatrixSS(),getSuspectResidualMeasurements(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements(),1);

}

void HypothesisTest::CalculateNMeasurements(){

    N_measurements = new float[getNumberSelectedMeasurements()];

    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        N_measurements[i] = fabs(estimated_error_measurements[i]); 
        N_measurements[i] += suspect_residual_covariance_matrix[i*getNumberSelectedMeasurements() + i]*sqrt(inverse_sensitivity_matrix_ss[i*getNumberSelectedMeasurements()+i]-1)*N_beta;
        N_measurements[i] /= suspect_residual_covariance_matrix[i*getNumberSelectedMeasurements()+i]*sqrt(inverse_sensitivity_matrix_ss[i*getNumberSelectedMeasurements()+i]);
    }
}

void HypothesisTest::CalculateThersholdMeasurements(){

    threshold_measurements = new float[getNumberSelectedMeasurements()];

    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        threshold_measurements[i] = suspect_residual_covariance_matrix[i*getNumberSelectedMeasurements()+i]*sqrt(inverse_sensitivity_matrix_ss[i*getNumberSelectedMeasurements()+i]);
        threshold_measurements[i] *= N_measurements[i];
    }
}

void HypothesisTest::SelectNewSuspectMeasurements(){
   
    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            new_number_suspect_measurements++;
        }
    } 

    float *new_suspect_selected_measurements = new float[new_number_suspect_measurements];

    int j = 0; 
    for(int i = 0; i < getNumberSelectedMeasurements(); i++){
        if(fabs(estimated_error_measurements[i]) > threshold_measurements[i]){
            new_suspect_selected_measurements[j] = i;
            j++;
        }
    }

    print(new_suspect_selected_measurements,1,new_number_suspect_measurements);

    float *new_suspect_residual_measurements = new float[new_number_suspect_measurements];
    float *new_suspect_residual_covariance_matrix = new float[new_number_suspect_measurements*new_number_suspect_measurements];
    float *new_sensitivity_matrix_ss = new float[new_number_suspect_measurements*new_number_suspect_measurements];

    int k;
    if(getNumberSelectedMeasurements() > new_number_suspect_measurements){
        
        //New Residual Measurements
        for(int i = 0; i < new_number_suspect_measurements; i++){
            k = new_suspect_selected_measurements[i];
            new_suspect_residual_measurements[i] = suspect_residual_measurements[k];
        }

        delete [] suspect_residual_measurements;

        suspect_residual_measurements = new float[new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            suspect_residual_measurements[i] = new_suspect_residual_measurements[i];
        }

        //New Residual Covariance Matrix
        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                if(i == j){
                    k = new_suspect_selected_measurements[i];
                    new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = suspect_residual_covariance_matrix[k*getNumberSelectedMeasurements() + k];
                }

                else{
                    new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = 0;
                }
            }
        }

        delete [] suspect_residual_covariance_matrix;

        suspect_residual_covariance_matrix = new float [new_number_suspect_measurements*new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j] = new_suspect_residual_covariance_matrix[i*new_number_suspect_measurements + j];
            }
        }

        //New Residual Sensitivity Matrix SS
        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                if(i==j){
                    k = new_suspect_selected_measurements[i];
                    new_sensitivity_matrix_ss[i*new_number_suspect_measurements + j] = sensitivity_matrix_SS[k*getNumberSelectedMeasurements() + k];
                }

                else{
                    new_sensitivity_matrix_ss[i*new_number_suspect_measurements + j] = 0;
                }
            }
        }

        delete [] sensitivity_matrix_SS;

        sensitivity_matrix_SS = new float[new_number_suspect_measurements];

        for(int i = 0; i < new_number_suspect_measurements; i++){
            for(int j = 0; j < new_number_suspect_measurements; j++){
                sensitivity_matrix_SS[i*new_number_suspect_measurements + j] = new_sensitivity_matrix_ss[i*getNumberSelectedMeasurements() + j];
            }
        }

        //New Number of Measurements
        number_selected_measurements = new_number_suspect_measurements;

        delete [] new_suspect_selected_measurements;
        delete [] new_suspect_residual_measurements;
        delete [] new_suspect_residual_covariance_matrix;
        delete [] new_sensitivity_matrix_ss;

        delete [] estimated_error_measurements;
        delete [] N_measurements;
        delete [] threshold_measurements;

    }

}

void HypothesisTest::HypothesisTestIdentification(float *residual_array,float *normalized_residual_array, float *sensitivity_matrix, float *residual_covariance_matrix){
    
    setSensitivityMatrix(sensitivity_matrix);
    setResidualCovarianceMatrix(residual_covariance_matrix);
    setResidualMeasurements(residual_array);
    setNormalizedMeasurements(normalized_residual_array);

    SelectSuspectMeasurements();
    SelectSuspectResidualCovarianceMatrix();
    SelectSuspectResidualMeasurements();
    std::cout << "Suspect Residual Measurements" << std::endl;
    print(getSuspectResidualMeasurements(),1,getNumberSelectedMeasurements());

    SelectSensitivityMatrixSS();
    std::cout << "Sensitivity Matrix SS" << std::endl;
    print(getSensitivityMatrixSS(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements());

    while (getNumberSelectedMeasurements() > new_number_suspect_measurements)
    {
    
        CalculateInverseSensitivityMatrixSS();
        std::cout << "Inverse Sensitivity Matrix SS" << std::endl;
        print(getInverseSensitivityMatrixSS(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements());
        CalculateEstimatedErrorMeasurements();
        std::cout << "Estimated Error" << std::endl;
        print(estimated_error_measurements,1,getNumberSelectedMeasurements());
        CalculateNMeasurements();
        std::cout << "N measurements" << std::endl;
        print(N_measurements,1,getNumberSelectedMeasurements());
        CalculateThersholdMeasurements();
        std::cout << "Threshold" << std::endl;    
        print(threshold_measurements,1,getNumberSelectedMeasurements());
        std::cout << "Measurements Identified with Error" << std::endl;
        SelectNewSuspectMeasurements();
        
    }

    std::cout << "Message: All the errors in the measurement set were identified" << std::endl;
    
}

void HypothesisTest::HypothesisTestIdentification(float *measurementArray, float *estimatedArray, float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix,const int length){
    
    CalculateHatMatrix(jacobianMatrix, gainMatrix, covarianceMatrix,length);
    CalculateSensitivityMatrix();
    CalculateResidualCovarianceMatrix(covarianceMatrix);
    CalculateResidualMeasurements(measurementArray,estimatedArray);
    CalculateNormalizedResidualMeasurements();

    SelectSuspectMeasurements();
    SelectSuspectResidualCovarianceMatrix();
    SelectSuspectResidualMeasurements();
    std::cout << "Suspect Residual Measurements" << std::endl;
    print(getSuspectResidualMeasurements(),1,getNumberSelectedMeasurements());

    SelectSensitivityMatrixSS();
    std::cout << "Sensitivity Matrix SS" << std::endl;
    print(getSensitivityMatrixSS(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements());

    while (getNumberSelectedMeasurements() > new_number_suspect_measurements)
    {
    
        CalculateInverseSensitivityMatrixSS();
        std::cout << "Inverse Sensitivity Matrix SS" << std::endl;
        print(getInverseSensitivityMatrixSS(),getNumberSelectedMeasurements(),getNumberSelectedMeasurements());
        CalculateEstimatedErrorMeasurements();
        std::cout << "Estimated Error" << std::endl;
        print(estimated_error_measurements,1,getNumberSelectedMeasurements());
        CalculateNMeasurements();
        std::cout << "N measurements" << std::endl;
        print(N_measurements,1,getNumberSelectedMeasurements());
        CalculateThersholdMeasurements();
        std::cout << "Threshold" << std::endl;    
        print(threshold_measurements,1,getNumberSelectedMeasurements());
        std::cout << "Measurements Identified with Error" << std::endl;
        SelectNewSuspectMeasurements();
        
    }

    std::cout << "Message: All the errors in the measurement set were identified" << std::endl << std::endl;
    
}