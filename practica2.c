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
	for (int i = from; mat[i = row][N] < numel; i++){
		printf("%f",mat[row][i]);
	}
printf("\n");
}

//Activitat 3

void MultEscalar( float vect[N], float vectres[N], float alfa ){
	for (int i = 0; i < N; i++){
		float res = alfa * vect[i];
		printf("%f",res);
	}
printf("\n");
}

//Activitat 4

float Scalar( float vect1[N], float vect2[N]) {
	float res = 0;
	for ( int i = 0; i < N; i++) {
		res += vect1[i] * vect2[i];
	return res;
	}
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

void Projection( float vect1[N], float vect2[N], float vectres[N] ) {
	float res = Scalar( vect1, vect2);
	float magn =  Magnitude(vect1);
	float div = (res / magn) * vect1[N];
	return div;
}

//Activitat 8

float Infininorm( float M[N][N] );
	sum = 0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			sum += fabs(M[i][j]);
		}
	}
	return sum;
}

//Activitat 9

float Onenorm( float M[N][N] );
	sum = 0;
	for (int j = 0; j < N; j++){
        	for (int i = 0; i < N; i++){
        		sum += fabs(M[i][j]);
                }
        }
        return sum;
}

//Activitat 10

float NormFrobenius( float M[N][N] );
	sum = 0;
        for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                        sum += M[i][j] ** 2;
                }
        }
	res = sqrt(sum);
        return res;
}

//Activitat 11

int DiagonalDom( float M[N][N] );



//Activitat 12

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] )

//Activitat 13

int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned iter )

//Main
void main(){
	InitData();
	Escalar(V3, 9, 2);
}
