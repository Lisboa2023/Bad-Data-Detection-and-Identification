#ifndef CHISQUARETEST_H
#define CHISQUARETEST_H

class ChiSquareDistribution{
    public:
        ChiSquareDistribution(const int, const int, const float);
        ~ChiSquareDistribution();
        
        //Funções SET======================================================
        void setNumberOfMeasurements(const int);
        void setDegreesOfFreedom(const int);
        void setConfidenceLevel(const float);
        void setChiSquaredVariable(const float);
                
        //Funções GET========================================================
        int getNumberOfMeasurements() const;
        int getDegreesOfFreedom() const;
        float getConfidenceLevel() const;
        float getChiSquaredVariable() const;
        float getProbability() const;
        //===================================================================
        
        void CalculateChiSquaredVariable(const float *, const float *, const float *);
        float ChiSquareProbabilityDensityFunction(const float);
        void ChiSquareCumulativeDistributionFunction(const float, const float);
        
        void ChiSquareTest(const float *, const float *, const float *);
        void ChiSquareTest(float);
        void printResult() const;

    private:
        int number_of_measurements;
        int degrees_of_freedom;
        float confidence_level;
        float chi_squared_variable;
        float probability;
};

#endif
