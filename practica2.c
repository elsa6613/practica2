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
	float scalar = res / magn;
	for (int i = 0; i < N; i++) {
        	vectres[i] = scalar * vect2[i];
	}
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
        vectres[i] = vect[i];
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
	printf("Activitat a\n");
        PrintVect(V1,0,10);
	PrintVect(V1,256,10);
	printf("\n");

	//b
	printf("Activitat b\n");
	PrintRow(Mat,0,0,10);
	PrintRow(Mat,100,0,10);
	printf("\n");

	//c
	printf("Activitat c\n");
	PrintRow(MatDD,0,0,10);
	PrintRow(MatDD,100,95,10);
	printf("\n");

	//d
	printf("Activitat d\n");
	printf("Infininorm(Mat): %f\n", Infininorm(Mat));
	printf("Onenorm(Mat): %f\n", Onenorm(Mat));
	printf("NormFrobenius(Mat): %f\n\n", NormFrobenius(Mat));
	if (DiagonalDom(Mat) == 1){
		printf("La matriu Mat és diagonal dominant\n");
	}else{
		printf("La matriu Mat no és diagonal dominant\n");
	}

	printf("Infininorm(MatDD): %f\n", Infininorm(MatDD));
        printf("Onenorm(MatDD): %f\n", Onenorm(MatDD));
        printf("NormFrobenius(MatDD): %f\n\n", NormFrobenius(MatDD));
	if (DiagonalDom(MatDD) == 1){
                printf("La matriu MatDD és diagonal dominant\n");
        }else{
                printf("La matriu MatDD no és diagonal dominant\n");
	}

	//e
	printf("Activitat e\n");
	printf("%f\n",Scalar(V1,V2));
	printf("%f\n",Scalar(V1,V3));
	printf("%f\n",Scalar(V2,V3));

	//f
	printf("Activitat f\n");
	printf("%f\n",Magnitude(V1));
	printf("%f\n",Magnitude(V2));
	printf("%f\n",Magnitude(V3));

	//g
	printf("Activitat g\n");
	if (Ortogonal(V1,V2) == 1){
                printf("V1 i V2 són ortogonals\n");
        }else{
                printf("V1 i V2 no són ortogonals\n");
        }

	if (Ortogonal(V1,V3) == 1){
                printf("V1 i V3 són ortogonals\n");
        }else{
                printf("V1 i V3 no són ortogonals\n");
        }

	 if (Ortogonal(V2,V3) == 1){
                printf("V2 i V3 són ortogonals\n");
        }else{
                printf("V2 i V3 no són ortogonals\n");
        }

	//h
	printf("Activitat h\n");
	float Escalar[N];
	MultEscalar(V3,Escalar,2.0);
	PrintVect(Escalar,0,10);
	PrintVect(Escalar,256,266);

	//i
	printf("Activitat i\n");
	float projeccio[N];
        Projection(V2,V3,projeccio);
        PrintVect(projeccio,0,11);
        float projeccio2[N];
        Projection(V1,V2,projeccio2);
        PrintVect(projeccio2,0,11);

	//j
	printf("Activitat j\n");
	float vectres[N];
	Matriu_x_Vector(Mat,V2,vectres);
	PrintVect(vectres,0,11);

	//k
	printf("Activitat k\n");
	float vectorres[N];
	Jacobi(MatDD,V3,vectorres,1);
	PrintVect(vectorres,0,11);

	float vectorres2[N];
        Jacobi(MatDD,V3,vectorres2,1000);
        PrintVect(vectorres,0,11);
}

