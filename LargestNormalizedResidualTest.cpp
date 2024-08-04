#include<iostream>
#include<cmath>
#include<iomanip>
#include"LargestNormalizedResidualTest.h"

NormalizedResidual::NormalizedResidual(const int SIZE, const float THRESHOLD){
    
    setNumberOfMeasurements(SIZE);
    setThreshold(THRESHOLD);

    hat_matrix = new float[SIZE*SIZE];
    sensitivity_matrix = new float[SIZE*SIZE];
    residual_covariance_matrix = new float[SIZE*SIZE];

    residual_measurements = new float[SIZE];
    normalized_measurements = new float[SIZE];
}

NormalizedResidual::~NormalizedResidual(){

    delete [] residual_measurements;
    delete [] normalized_measurements;
    delete [] hat_matrix;
    delete [] sensitivity_matrix;
    delete [] residual_covariance_matrix;
}

//Funcoes SET ==========================================================================================
void NormalizedResidual::setNumberOfMeasurements(const int SIZE){
    number_of_measurements = SIZE;
}

void NormalizedResidual::setThreshold(const float THRESHOLD){
    threshold = THRESHOLD;
}

void NormalizedResidual::setHatMatrix(float *hMatrix){
    hat_matrix = hMatrix;
}

void NormalizedResidual::setSensitivityMatrix(float *sMatrix){
    sensitivity_matrix = sMatrix;
}

void NormalizedResidual::setResidualCovarianceMatrix(float *rcMatrix){
    residual_covariance_matrix = rcMatrix;
}

void NormalizedResidual::setResidualMeasurements(float *rArray){
    residual_measurements = rArray;
}

void NormalizedResidual::setNormalizedMeasurements(float *nArray){
    normalized_measurements = nArray;
}


//=============================================================================================


//Funcoes GET============================================================
int NormalizedResidual::getNumberOfMeasurements() const{
    return number_of_measurements;
}

float NormalizedResidual::getThreshold() const{
    return threshold;
}

float *NormalizedResidual::getHatMatrix() const{
    return hat_matrix;
}

float *NormalizedResidual::getSensitivityMatrix() const{
    return sensitivity_matrix;
}

float *NormalizedResidual::getResidualCovarianceMatrix() const{
    return residual_covariance_matrix;
}

float *NormalizedResidual::getResidualMeasurements() const{
    return residual_measurements;
}

float *NormalizedResidual::getNormalizedMeasurements() const{
    return normalized_measurements;
}

//=======================================================================

float *NormalizedResidual::CalculateInverseMatrix(float *matrix, const int length){
    float *inverseMatrix = new float[length*length];
    inverseMatrix = matrix;

    //permutacao de linhas
    float temp;
    for(int i = 0; i < length; i++)
    { 
        for(int j = 0; j < length; j++)
        {
            if(i==j && inverseMatrix[i*length + j] == 0)
            {   
                for(int p = i+1; p < length; p++)
                {   
                    if(inverseMatrix[p*length + j] != 0)
                    {
                        for(int n = 0; n < length; n++)
                        {
                            temp = inverseMatrix[i*length + n];
                            inverseMatrix[i*length + n] = inverseMatrix[p*length + n];
                            inverseMatrix[p*length + n] = temp;
                        }
                        break;
                    }
                }
            }
        }
    }

    //algoritmo de eliminacao de gauss-jordan
    float identityMatrix[length][length];

    
    //Criando matriz identidade
    for(int i = 0; i < length; i++){
        for(int j = 0;j < length; j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    
    for (int i = 0; i < length; i++)
	{
		// Dividindo a linha atual pelo elemento diagonal correspondente
		float pivot = inverseMatrix[i*length + i];
		for ( int j = 0; j < length; j++)
		{
			inverseMatrix[i*length + j] /= pivot;
			identityMatrix[i][j] /= pivot; //As mesmas operações são feitas na matriz identidade para se obter a inversa
		}

		// Reduzindo as outras linhas
		for (int j = 0; j < length; j++)
		{
			if (j != i)
			{
				float a = inverseMatrix[j*length + i];
				for (int k = 0; k < length; k++)
				{
					inverseMatrix[j*length + k] -= a * inverseMatrix[i*length + k]; //Colocando 0 abaixo e acima dos pivos
					identityMatrix[j][k] -= a * identityMatrix[i][k];//Repetindo operações na matriz identidade
				}
			}
		}
	}

    for(int i = 0; i < length; i++){
        for(int j = 0; j < length; j++){
            inverseMatrix[i*length + j] = identityMatrix[i][j];
        }
    }

    return inverseMatrix;

}

float *NormalizedResidual::CalculateTransposedMatrix(const float *matrix, const int ROWS, const int COLUMNS){

    float *transposedMatrix = new float[ROWS*COLUMNS];

    for(int i = 0; i < ROWS; i++){
        for(int j = 0; j < COLUMNS; j++){
            transposedMatrix[j*ROWS + i] = matrix[i*COLUMNS + j];
        }
    }

    return transposedMatrix;
}

float *NormalizedResidual::MultiplyArray(const float *array_a,const float *array_b, const int rows_a, const int columns_a, const int rows_b, const int columns_b){

    float *temp = new float[rows_a*columns_b];

    for(int i = 0; i < rows_a; i++){
        for(int j = 0; j < columns_b; j++){
            for(int k = 0; k < rows_b; k++){
                if(k == 0){
                    temp[i*columns_b + j] = array_a[i*columns_a + k]*array_b[k*columns_b + j];
                }

                else{   
                    temp[i*columns_b + j] += array_a[i*columns_a + k]*array_b[k*columns_b + j];
                }
            }
        }
    }

    return temp;
}

void NormalizedResidual::CalculateHatMatrix(float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix, const int length){
    
    //hatMatrix = matriz jacobiana * matriz de ganho invertida * matriz jacobiana transposta * matriz de covariancia invertida
    float *tempInverse = CalculateInverseMatrix(gainMatrix,length);
    float *temp = MultiplyArray(jacobianMatrix,tempInverse,getNumberOfMeasurements(),length,length,length);

    float *tempTransposed = CalculateTransposedMatrix(jacobianMatrix,getNumberOfMeasurements(),length);
    temp = MultiplyArray(temp,tempTransposed,getNumberOfMeasurements(),length,length,getNumberOfMeasurements());

    tempInverse = CalculateInverseMatrix(covarianceMatrix,getNumberOfMeasurements());
    hat_matrix = MultiplyArray(temp,tempInverse,getNumberOfMeasurements(),getNumberOfMeasurements(),getNumberOfMeasurements(),getNumberOfMeasurements());

    //Retorna as matrizes de Ganhio e Covariancia para seus valores originais
    CalculateInverseMatrix(covarianceMatrix,getNumberOfMeasurements());
}

void NormalizedResidual::CalculateSensitivityMatrix(){

    //Criando Matrix Identidade =============================================
    bool identityMatrix[getNumberOfMeasurements()][getNumberOfMeasurements()];
    for(int i = 0; i < getNumberOfMeasurements(); i++){
        for(int j = 0; j < getNumberOfMeasurements(); j++){
            if(i==j){
                identityMatrix[i][j] = 1;
            }
            else{
                identityMatrix[i][j] = 0;
            }
        }
    }
    //======================================================================

    for(int i = 0; i < getNumberOfMeasurements(); i++){
        for(int j = 0; j < getNumberOfMeasurements(); j++){
            sensitivity_matrix[i*getNumberOfMeasurements() + j] = identityMatrix[i][j] - hat_matrix[i*getNumberOfMeasurements() + j];
        }
    } 

}

void NormalizedResidual::CalculateResidualCovarianceMatrix(float *covarianceMatrix){

    residual_covariance_matrix = MultiplyArray(sensitivity_matrix,covarianceMatrix,getNumberOfMeasurements(),getNumberOfMeasurements(),getNumberOfMeasurements(),getNumberOfMeasurements());

}

void NormalizedResidual::CalculateResidualMeasurements(const float *measurements,const float *estimated_measurements){

    for(int i = 0; i < getNumberOfMeasurements(); i++){
        residual_measurements[i] = measurements[i] - estimated_measurements[i];
    }

}

void NormalizedResidual::CalculateNormalizedResidualMeasurements(){

    for(int i = 0; i < getNumberOfMeasurements(); i++){
        normalized_measurements[i] = fabs(residual_measurements[i])/sqrt(residual_covariance_matrix[i*getNumberOfMeasurements() + i]);
    }

}

void NormalizedResidual::FindLargestResidual(float &temp, int &pos){
    temp = *normalized_measurements;
    for(int i = 0; i < getNumberOfMeasurements(); i++){
        if(normalized_measurements[i] >= temp){
            temp = normalized_measurements[i];
            pos = i;
        }
    }
}

void NormalizedResidual::DeleteError(const int threshold, const float lg, const int p,float* measurements,float *estimated_measurements){
    if(lg > threshold){
        residual_measurements[p] = 0;
        normalized_measurements[p] = 0;
        measurements[p] = 0;
        estimated_measurements[p] = 0;
    }  
}

void NormalizedResidual::print(const float *array, const int rows, const int columns){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            std::cout << std::setw(20) << array[i*columns + j];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, float *residualCovarianceMatrix){

    setResidualCovarianceMatrix(residualCovarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i < getNumberOfMeasurements(); i++){
        CalculateResidualMeasurements(measurementArray,estimatedArray);
        CalculateNormalizedResidualMeasurements();
        if(i==0){
            std::cout << "Residual Measurements Set: " << std::endl; 
            print(getResidualMeasurements(), 1,getNumberOfMeasurements());
            std::cout << std::endl << "Normalized Residual Measurements Set: " << std::endl; 
            print(getNormalizedMeasurements(), 1,getNumberOfMeasurements());
            std::cout << std::endl << "Threshold: " << std::setprecision(3) <<getThreshold() << std::endl;
        }
        FindLargestResidual(largestResidual, position);
        if (largestResidual > getThreshold()){ 
            DeleteError(getThreshold(), largestResidual, position,measurementArray,estimatedArray);
            print(getNormalizedMeasurements(), 1, getNumberOfMeasurements());
        }
        else{
            std::cout << std::endl << "Message: Measurement set free of errors!" << std::endl <<
            std::endl;
            break;
        }
    }
}

void NormalizedResidual::LargestNormalizedResidualTest(float *measurementArray, float *estimatedArray, float *jacobianMatrix, float *gainMatrix, float *covarianceMatrix,const int length){

    CalculateHatMatrix(jacobianMatrix, gainMatrix, covarianceMatrix,length);
    CalculateSensitivityMatrix();
    CalculateResidualCovarianceMatrix(covarianceMatrix);

    float largestResidual;
    int position;

    for(int i=0; i < getNumberOfMeasurements(); i++){
        CalculateResidualMeasurements(measurementArray,estimatedArray);
        CalculateNormalizedResidualMeasurements();
        if(i==0){
            std::cout << "Residual Measurement Set: " << std::endl; 
            print(getResidualMeasurements(), 1,getNumberOfMeasurements());
            std::cout << std::endl << "Normalized Residual Measurement Set: " << std::endl; 
            print(getNormalizedMeasurements(), 1, getNumberOfMeasurements());
            std::cout << std::endl << "Threshold: " << std::setprecision(3) << getThreshold() << std::endl;
        }
        FindLargestResidual(largestResidual, position);
        if (largestResidual > getThreshold()){ 
            DeleteError(getThreshold(), largestResidual, position,measurementArray,estimatedArray);
            print(getNormalizedMeasurements(), 1, getNumberOfMeasurements());
        }
        else{
            std::cout << std::endl << "Message: Measurement set free of errors!" << std::endl <<
            std::endl;
            break;
        }
    }
}
