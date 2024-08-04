#ifndef LNRTEST_H
#define LNRTEST_H

class NormalizedResidual{
    public:
        NormalizedResidual(const int, const float);
        ~NormalizedResidual();

        //Funcoes SET ============================================================
        void setNumberOfMeasurements(const int);
        void setThreshold(const float);

        void setHatMatrix(float *);
        void setSensitivityMatrix(float *);
        void setResidualCovarianceMatrix(float *);

        void setResidualMeasurements(float *);
        void setNormalizedMeasurements(float *);
        //=======================================================================

        //Funcoes GET============================================================
        int getNumberOfMeasurements() const;
        float getThreshold() const;

        float *getHatMatrix() const;
        float *getSensitivityMatrix() const;
        float *getResidualCovarianceMatrix() const;
        
        float *getResidualMeasurements() const;
        float *getNormalizedMeasurements() const;

        //=======================================================================
        
        float *CalculateInverseMatrix(float *, const int);
        float *CalculateTransposedMatrix(const float *, const int, const int);
        float *MultiplyArray(const float *, const float *, const int, const int, const int, const int);

        void CalculateHatMatrix(float *, float *,float *, const int);
        void CalculateSensitivityMatrix();
        void CalculateResidualCovarianceMatrix(float *);

        void CalculateResidualMeasurements(const float *,const float *);
        void CalculateNormalizedResidualMeasurements();

        void FindLargestResidual(float &, int &);
        void DeleteError(const int, const float, const int,float *,float *);
        void print(const float*, const int, const int);

        void LargestNormalizedResidualTest(float *, float *, float *);
        void LargestNormalizedResidualTest(float *, float *, float *, float *, float*, const int);
        
    private:
        int number_of_measurements;
        float threshold;

        float *hat_matrix;
        float *sensitivity_matrix;
        float *residual_covariance_matrix;

        float *residual_measurements;
        float *normalized_measurements;

};

#endif