#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512
float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];
float V4[N];

void InitData(){
	int i,j;
	srand(334411);
	for( i = 0; i < N; i++ )
		for( j = 0; j < N; j++ ){
 			Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
 			if ( (abs(i - j) <= 3) && (i != j))
 				MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 			else if ( i == j )
 				MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 			else MatDD[i][j] = 0.0;
 		}
	for( i = 0; i < N; i++ ){
 		V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 		V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 		V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
	}
}


//Activitat 1

void PrintVect( float vect[N], int from, int numel ){
	for ( int i = from; i < numel; i++){
		printf("%f ",vect[i]);
	}
	printf("\n");
}

//Activitat 2


void PrintRow( float mat[N][N], int row, int from, int numel ){
	for (int i = from; i < numel; i++){
		printf("%f",mat[row][i]);
	}
printf("\n");
}

//Activitat 3

void MultEscalar( float vect[N], float vectres[N], float alfa ){
	for (int i = 0; i < N; i++){
		vectres[i] = alfa * vect[i];
		printf("%f",vectres[i]);
	}
printf("\n");
}

//Activitat 4

float Scalar( float vect1[N], float vect2[N]) {
	float res = 0;
	for ( int i = 0; i < N; i++) {
		res += vect1[i] * vect2[i];
	}
	return res;
}

//Activitat 5

float Magnitude( float vect[N] ) {
	float res = Scalar( vect, vect );
	float sol = sqrt(res);
	return sol;
}

//Activitat 6

int Ortogonal( float vect1[N], float vect2[N] ) {
	float res = Scalar( vect1, vect2 );
	if ( res == 0 )
		return 1;
	else
		return 0;
}
//Activitat 7

int Projection( float vect1[N], float vect2[N], float vectres[N] ) {
	float res = Scalar( vect1, vect2);
	float magn =  Magnitude(vect1);
	float div = (res / magn) * vect1[N];
	return div;
}

//Activitat 8

float Infininorm( float M[N][N] ) {
        float sum = 0;
        float max = 0;
        for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                        sum += fabs(M[i][j]);
                }
        }
        if (sum > max)
                max = sum;
        return max;
}

//Activitat 9

float Onenorm( float M[N][N] ) {
        float sum = 0;
        float max = 0;
        for (int j = 0; j < N; j++){
                for (int i = 0; i < N; i++){
                        sum += fabs(M[i][j]);
                }
        	if (sum > max){
                	max = sum;
		}
	}
        return max;
}

//Activitat 10

float NormFrobenius( float M[N][N] ) {
        float sum = 0;
        for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                                sum += fabs(M[i][j]*M[i][j]);
                }
        }
        float res = sqrt(sum);
        return res;
}

//Activitat 11

int DiagonalDom( float M[N][N] ) {
        for (int i = 0; i < N; i++){
                float diagonal = fabs(M[i][i]);
                float sum = 0;
                for (int j = 0; j < N; j++){
                        if (i!=j){
                                sum += fabs(M[i][j]);
                        }
                }
        	if ( diagonal < sum ){
                	return 0;
		}
        }
	return 1;
}

//Activitat 12

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ) {
        for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++) {
			vectres[i] += M[i][j]*vect[j];
                }
        }
}

//Activitat 13

int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned iter ){
    float nou[N];
    int i, j;
    float suma, diferencia;

    for (i = 0; i < N; i++) {
        vectres[i] = 0.0;
    }

    for (unsigned k = 0; k < iter; k++) {
        for (i = 0; i < N; i++) {
            //suma `vect[i] - sum(A[i][j] * x[j])` per a totes les j ≠ i
            suma = vect[i];
            for (j = 0; j < N; j++) {
                if (j != i) {
                    suma -= M[i][j] * vectres[j];
                }
            }
            nou[i] = suma / M[i][i];
        }
        //comprovar convergència(comparar nous valors i antics)
        diferencia = 0.0;
        for (i = 0; i < N; i++) {
            diferencia += fabs(nou[i] - vectres[i]);
            vectres[i] = nou[i];
        }

        if (diferencia < 1e-6) { //si diferencia menor que tolerancia, sistema convergit
            return 1;
        }
    }

    return 0;
}
//Main
int main(){
	InitData();
	
	//a
        PrintVect(V1,0,10);
	PrintVect(V1,256,10);
	printf("\n");

	//b
	PrintRow(Mat,0,0,10);
	PrintRow(Mat,100,0,10);
	printf("\n");

	//c
	PrintRow(MatDD,0,0,10);
	PrintRow(MatDD,100,95,10);
	printf("\n");

	//d
	printf("Infininorm(Mat): %f\n", Infininorm(Mat));
	printf("Onenorm(Mat): %f\n", Onenorm(Mat));
	printf("NormFrobenius(Mat): %f\n\n", NormFrobenius(Mat));
	
	printf("Infininorm(MatDD): %f\n", Infininorm(MatDD));
        printf("Onenorm(MatDD): %f\n", Onenorm(MatDD));
        printf("NormFrobenius(MatDD): %f\n\n", NormFrobenius(MatDD));
}
