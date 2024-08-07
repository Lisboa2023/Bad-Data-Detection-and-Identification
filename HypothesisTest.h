#ifndef HTI_H
#define HTI_H

#include"LargestNormalizedResidualTest.h"

class HypothesisTest:public NormalizedResidual{
    public:
        HypothesisTest(const float, const float, const int, const float);
        ~HypothesisTest();

        //Funcoes SET ============================================================
        void setNBeta(const float);
        void setNMaximus(const float);
        void setNumberSelectedMeasurements(int);

        //=======================================================================     

        //Funcoes GET ===========================================================
        int getNumberSelectedMeasurements() const;
        float getNBeta() const;
        float getNMaximus() const;

        float *getSuspectSelectedMeasurements()const;
        float *getSuspectResidualCovarianceMatrix()const;
        float *getSensitivityMatrixSS()const;
        float *getInverseSensitivityMatrixSS()const;
        float *getSuspectResidualMeasurements()const;

        //=========================================================================  

        void SelectSuspectMeasurements();

        void SelectSuspectResidualCovarianceMatrix();
        void SelectSensitivityMatrixSS();
        void SelectSuspectResidualMeasurements();
        
        void CalculateInverseSensitivityMatrixSS();
        void CalculateEstimatedErrorMeasurements();
        void CalculateNMeasurements();
        void CalculateThersholdMeasurements();

        void SelectNewSuspectMeasurements();

        void HypothesisTestIdentification(float*,float*,float*,float*);
        void HypothesisTestIdentification(float*,float*,float*,float*,float*,const int);
        void BadDataIdentification(float*,float*,float*,float*,float*,const int);

    private:
        int number_selected_measurements;
        float N_beta;
        float N_maximus;
        float *suspect_selected_measurements;
        float *suspect_residual_covariance_matrix;
        float *sensitivity_matrix_SS;
        float *inverse_sensitivity_matrix_ss;
        float *suspect_residual_measurements;
        float *estimated_error_measurements;
        float *N_measurements;
        float *threshold_measurements;
        int new_number_suspect_measurements;

};

#endif