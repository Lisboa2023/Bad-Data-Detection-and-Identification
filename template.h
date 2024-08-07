//MATRIZ INVERSA
template <class T>
T *CalculateInverseMatrix(T *matrix, const int length){
    T *inverseMatrix = new T[length*length];
    inverseMatrix = matrix;

    //permutacao de linhas
    T temp;
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
    T identityMatrix[length][length];

    
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
		T pivot = inverseMatrix[i*length + i];
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
				T a = inverseMatrix[j*length + i];
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

//MATRIZ TRANSPOSTA ======================================================================================================
template<class T>
T *CalculateTransposedMatrix(const T *matrix, const int ROWS, const int COLUMNS){

    float *transposedMatrix = new float[ROWS*COLUMNS];

    for(int i = 0; i < ROWS; i++){
        for(int j = 0; j < COLUMNS; j++){
            transposedMatrix[j*ROWS + i] = matrix[i*COLUMNS + j];
        }
    }

    return transposedMatrix;
}

//MULTIPLICAÇÃO DE MATRIZES================================================================================
template<class T>
T *MultiplyArray(const T *array_a,const T *array_b, const int rows_a, const int columns_a, const int rows_b, const int columns_b){

    T *temp = new T[rows_a*columns_b];

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
