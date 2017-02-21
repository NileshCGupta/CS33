#include "func.h"
#include <omp.h>

void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
    int i, j, k, r, num;

    #pragma omp parallel for private(i, j, k, r) collapse(1)
    for (i=0; i<TSIZE; i++) {
	for (k=0; k<TSIZE; k++) {
	    r = a[i][k];
	    for (j=0; j<TSIZE; j++) {
		c[i][j]+=r*b[k][j];
	    }
	}
    }
}


/*void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
    int i, j, k, i1, j1, k1, r, num;
#define block 256
#pragma omp parallel for private(i, j, k) collapse(3)
    for (i=0; i<TSIZE; i+=block) {
	for (j=0; j<TSIZE; j+=block) {
	    for (k=0; k<TSIZE; k+=block) {
#pragma omp parallel for private(i1, j1, k1, r) collapse(2)
		for (i1=i; i1<i+block; i1++){
		    for (k1=k; k1<k+block; k1++){
			r = a[i1][k1];
			for (j1=j; j1<j+block; j1++)
			    c[i1][j1]+=r*b[k1][j1];
		    }
		}
	    }
	}
    }
    }*/


void Func2(int d[][MSIZE], int c[][MSIZE])
{
    /*int i,j;
    #pragma omp parallel for private(i, j) collapse(2)
    for(i=0; i<MSIZE; i++) 
	for(j=0; j<MSIZE; j++)
	d[j][i]=c[i][j];*/


    int i, j, k, l, kn;
    #define blocksize 256
    int n = MSIZE;

#pragma omp parallel for private(i, j) collapse(2)
    for (i = 0; i < n; i += blocksize) {
	for (j = 0; j < n; j += blocksize) {
	    // transpose the block beginning at [i,j]
#pragma omp parallel for private(k, l) collapse(2)
	    for (k = i; k < i + blocksize; ++k) {
//		kn = k*n;
//#pragma omp parallel for private(l)
		for (l = j; l < j + blocksize; ++l) {
		    d[l][k] = c[k][l];
		}
	    }
	}
    }
}

/*void Func3(int z[][MSIZE], int d[][MSIZE])
{
	int y, x;
	int near = 2;  		// The weight of neighbor
	int itself = 84; 	// The weight of center
	//int current, xup, xdown, yup, ydown, bothup, bothdown, it;

#pragma omp parallel for collapse(2)
	for (y=0; y<MSIZE; y++) {
		for (x=0; x<MSIZE; x++) {
		    //it = (d[y][x] << 6) + (d[y][x] << 2) + (d[y][x] << 4);
		    //current = d[y][x] << 1;
			if (y==0) {
				if (x==0) {
				    //xup = d[y][x+1] << 1;
				    //yup = d[y+1][x] << 1;
				    //bothup = d[y+1][x+1] << 1;
				    z[y][x] = 	near * d[y][x] +
						near * d[y][x+1] +
						near * d[y+1][x] +
						near * d[y+1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				    // z[y][x] = current + xup + yup + bothup + current + xup + current + yup + it;
				} else if (x==MSIZE-1) {
				    
					z[y][x] = 	near * d[y][x-1] +
						near * d[y][x] +
						near * d[y+1][x-1] +
						near * d[y+1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else {
					z[y][x] = 	near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y+1][x-1] +
						near * d[y+1][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				}
			} else if (y==MSIZE-1) {
				if (x==0) {
					z[y][x] = 	near * d[y-1][x] +
						near * d[y-1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				} else {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y][x] +
						itself * d[y][x];
				}
			} else {
				if (x==0) {
					z[y][x] = 	near * d[y-1][x] +
						near * d[y-1][x+1] +
						near * d[y+1][x] +
						near * d[y+1][x+1] +
						near * d[y][x] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x] +
						near * d[y+1][x-1] +
						near * d[y+1][x] +
						near * d[y][x-1] +
						near * d[y][x] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				} else {
					z[y][x] = 	near * d[y-1][x-1] +
						near * d[y-1][x+1] +
						near * d[y+1][x-1] +
						near * d[y+1][x+1] +
						near * d[y][x-1] +
						near * d[y][x+1] +
						near * d[y-1][x] +
						near * d[y+1][x] +
						itself * d[y][x];
				}
			}
			z[y][x]/=100;
		}
	}
}*/

void Func3(int z[][MSIZE], int d[][MSIZE])
{
	int y, x;
	int near = 2;  		// The weight of neighbor
	int itself = 84; 	// The weight of center
	//int current, xup, xdown, yup, ydown, bothup, bothdown, it;

#pragma omp parallel for collapse(2)
	for (y=0; y<MSIZE; y++) {
		for (x=0; x<MSIZE; x++) {
		    //it = (d[y][x] << 6) + (d[y][x] << 2) + (d[y][x] << 4);
		    //current = d[y][x] << 1;
			if (y==0) {
				if (x==0) {
				    //xup = d[y][x+1] << 1;
				    //yup = d[y+1][x] << 1;
				    //bothup = d[y+1][x+1] << 1;
				    z[y][x] = 	(d[y][x] +
						 d[y][x+1] +
						 d[y+1][x] +
						 d[y+1][x+1] +
						 d[y][x] +
						 d[y][x+1] +
						 d[y][x] +
						 d[y+1][x]) * near +
					itself * d[y][x];
				    // z[y][x] = current + xup + yup + bothup + current + xup + current + yup + it;
				} else if (x==MSIZE-1) {
				    
				    z[y][x] = 	(d[y][x-1] +
						 d[y][x] +
						 d[y+1][x-1] +
						 d[y+1][x] +
						d[y][x-1] +
					        d[y][x] +
					        d[y][x] +
						 d[y+1][x]) * near +
						itself * d[y][x];
				} else {
				    z[y][x] = 	(d[y][x-1] +
						d[y][x+1] +
						d[y+1][x-1] +
						d[y+1][x+1] +
					       d[y][x-1] +
						d[y][x+1] +
						d[y][x] +
						 d[y+1][x]) * near +
						itself * d[y][x];
				}
			} else if (y==MSIZE-1) {
				if (x==0) {
				    z[y][x] = 	(d[y-1][x] +
					        d[y-1][x+1] +
						d[y][x] +
						d[y][x+1] +
					       d[y][x] +
						d[y][x+1] +
						d[y-1][x] +
						 d[y][x]) * near +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
				    z[y][x] = 	(d[y-1][x-1] +
						d[y-1][x] +
					        d[y][x-1] +
					       d[y][x] +
					        d[y][x-1] +
						d[y][x] +
					       d[y-1][x] +
						 d[y][x]) * near +
						itself * d[y][x];
				} else {
				    z[y][x] = 	(d[y-1][x-1] +
						d[y-1][x+1] +
					       d[y][x-1] +
						d[y][x+1] +
						d[y][x-1] +
						d[y][x+1] +
					        d[y-1][x] +
						 d[y][x]) * near +
						itself * d[y][x];
				}
			} else {
				if (x==0) {
				    z[y][x] = 	(d[y-1][x] +
						d[y-1][x+1] +
					       d[y+1][x] +
						d[y+1][x+1] +
						d[y][x] +
						d[y][x+1] +
						d[y-1][x] +
						 d[y+1][x]) * near +
						itself * d[y][x];
				} else if (x==MSIZE-1) {
				    z[y][x] = 	(d[y-1][x-1] +
						d[y-1][x] +
						d[y+1][x-1] +
					        d[y+1][x] +
					        d[y][x-1] +
						d[y][x] +
					        d[y-1][x] +
						 d[y+1][x]) * near +
						itself * d[y][x];
				} else {
				    z[y][x] = 	(d[y-1][x-1] +
						d[y-1][x+1] +
						d[y+1][x-1] +
						d[y+1][x+1] +
						d[y][x-1] +
						d[y][x+1] +
						d[y-1][x] +
						 d[y+1][x]) * near +
						itself * d[y][x];
				}
			}
			z[y][x]/=100;
		}
	}
}


/*void Func3(int z[][MSIZE], int d[][MSIZE])
{
	int y, x;
	int near = 2;  		// The weight of neighbor
	int itself = 84; 	// The weight of center
	//int current, xup, xdown, yup, ydown, bothup, bothdown, it;

	//first row
	z[0][0] = near * d[0][0] +
	    near * d[0][1] +
	    near * d[1][0] +
	    near * d[1][1] +
	    near * d[0][0] +
	    near * d[0][1] +
	    near * d[0][0] +
	    near * d[1][0] +
	    itself * d[0][0];
	z[0][0]/=100;

	z[0][MSIZE-1] = near * d[0][MSIZE-2] +
	    near * d[0][MSIZE-1] +
	    near * d[1][MSIZE-2] +
	    near * d[1][MSIZE-1] +
	    near * d[0][MSIZE-2] +
	    near * d[0][MSIZE-1] +
	    near * d[0][MSIZE-1] +
	    near * d[1][MSIZE-1] +
	    itself * d[0][MSIZE-1];
	z[0][MSIZE-1]/=100;

#pragma omp parallel for private(x)
	for (x=1; x<MSIZE-1; x++){	  
	    z[0][x] = 	near * d[0][x-1] +
		near * d[0][x+1] +
		near * d[1][x-1] +
		near * d[1][x+1] +
		near * d[0][x-1] +
		near * d[0][x+1] +
		near * d[0][x] +
		near * d[1][x] +
		itself * d[0][x];
	    z[0][x]/=100;
	}

	//last row
	z[MSIZE-1][0] = near * d[MSIZE-2][0] +
	    near * d[MSIZE-2][1] +
	    near * d[MSIZE-1][0] +
	    near * d[MSIZE-1][1] +
	    near * d[MSIZE-1][0] +
	    near * d[MSIZE-1][1] +
	    near * d[MSIZE-2][0] +
	    near * d[MSIZE-1][0] +
	    itself * d[MSIZE-1][0];
	z[MSIZE-1][0]/=100;

	z[MSIZE-1][MSIZE-1] = near * d[MSIZE-2][MSIZE-2] +
	    near * d[MSIZE-2][MSIZE-1] +
	    near * d[MSIZE-1][MSIZE-2] +
	    near * d[MSIZE-1][MSIZE-1] +
	    near * d[MSIZE-1][MSIZE-2] +
	    near * d[MSIZE-1][MSIZE-1] +
	    near * d[MSIZE-2][MSIZE-1] +
	    near * d[MSIZE-1][MSIZE-1] +
	    itself * d[MSIZE-1][MSIZE-1];
	z[MSIZE-1][MSIZE-1]/=100;

#pragma omp parallel for private(x)
	for (x=0; x<MSIZE; x++){
	    z[MSIZE-1][x] =   near * d[MSIZE-2][x-1] +
		near * d[MSIZE-2][x+1] +
		near * d[MSIZE-1][x-1] +
		near * d[MSIZE-1][x+1] +
		near * d[MSIZE-1][x-1] +
		near * d[MSIZE-1][x+1] +
		near * d[MSIZE-2][x] +
		near * d[MSIZE-1][x] +
		itself * d[MSIZE-1][x];
	    z[MSIZE-1][x]/=100;
	}

	//sides

#pragma omp parallel for private(y)
	for (y=1; y<MSIZE-1; y++){
	     z[y][0] =       near * d[y-1][0] +
		 near * d[y-1][1] +
		 near * d[y+1][0] +
		 near * d[y+1][1] +
		 near * d[y][0] +
		 near * d[y][1] +
		 near * d[y-1][0] +
		 near * d[y+1][0] +
		 itself * d[y][0];
	     z[y][0]/=100;
	}

#pragma omp parallel for private(y)
	for (y=1; y<MSIZE-1; y++){
	     z[y][MSIZE-1] =       near * d[y-1][MSIZE-2] +
		 near * d[y-1][MSIZE-1] +
		 near * d[y+1][MSIZE-2] +
		 near * d[y+1][MSIZE-1] +
		 near * d[y][MSIZE-2] +
		 near * d[y][MSIZE-1] +
		 near * d[y-1][MSIZE-1] +
		 near * d[y+1][MSIZE-1] +
		 itself * d[y][MSIZE-1];
	     z[y][MSIZE-1]/=100;
	}

	//body
#pragma omp parallel for collapse(2)
	for (y=1; y<MSIZE-1; y++) {
		for (x=1; x<MSIZE-1; x++) {
		    //it = (d[y][x] << 6) + (d[y][x] << 2) + (d[y][x] << 4);
		    //current = d[y][x] << 1;
		    z[y][x] = 	near * d[y-1][x-1] +
			near * d[y-1][x+1] +
			near * d[y+1][x-1] +
			near * d[y+1][x+1] +
			near * d[y][x-1] +
			near * d[y][x+1] +
			near * d[y-1][x] +
			near * d[y+1][x] +
			itself * d[y][x];
		    z[y][x]/=100;
		}
	}
	}*/

void Func4(int b[], int a[])
{
	int chuck_size=MSIZE;	 
	int array_size=VSIZE/chuck_size;
	int chuck[chuck_size];
	int i, j, num, jarray=0;
#pragma omp parallel for private(i, jarray)
	for(j=0; j<chuck_size; j++) {
	    //a[j*array_size]=b[j*array_size];
	    jarray=j*array_size;
	    b[jarray]=a[jarray];
	    //#pragma omp parallel for
		for (i=1; i<array_size; i++) {
		    //b[j*array_size+i]=b[j*array_size+i-1]+a[j*array_size+i];
		    b[jarray+i]=b[jarray+i-1]+a[jarray+i];
		}
		//jarray += array_size;
		chuck[j]=b[(jarray+array_size)-1];
	}
// #pragma omp parallel for private(j)
	/*for(j=1; j<chuck_size; j++) {
		chuck[j]=chuck[j-1]+chuck[j];
		}*/
	for(j=0; j<chuck_size-1; j+=8)
	{
	    chuck[j+1]=chuck[j]+chuck[j+1];
	    chuck[j+2]=chuck[j+1]+chuck[j+2];
	    chuck[j+3]=chuck[j+2]+chuck[j+3];
	    chuck[j+4]=chuck[j+3]+chuck[j+4];
	    chuck[j+5]=chuck[j+4]+chuck[j+5];
	    chuck[j+6]=chuck[j+5]+chuck[j+6];
	    chuck[j+7]=chuck[j+6]+chuck[j+7];
	    chuck[j+8]=chuck[j+7]+chuck[j+8];
	}


       
/*//#pragma omp parallel for private(i, num)
	for(j=1; j<chuck_size; j++) {
	    num=chuck[j-1];
//#pragma omp parallel for private(i, num)
		for (i=0; i<VSIZE / chuck_size; i++) {
			b[j*array_size+i]+=num/(j+1);
		}
		}*/
#pragma omp parallel for private(i, jarray, num)
	for(j = 1; j<chuck_size; j++){
	    jarray=j*array_size;
	    num=chuck[j-1]/(j+1);
	    for (i = 0; i<array_size; i++){
		b[jarray + i] += num;
	    }
	}
}
#undef block
void Func5(int b[], int a[])
{
    int i=0, j, j1,  stride, stride2, step;
    int temp;
    long log_N=log2(VSIZE);

#define block 64
#pragma omp parallel private(j,j1)
#pragma omp for schedule (guided)
	for(j=0; j<VSIZE; j+=block) {
	    for(j1=j; j1 < j+block; j1+=2){
		b[j1]=a[j1];
		b[j1+1] = a[j1] + a[j1+1];
	    }
	}
	for(i=4; i<VSIZE; i*=2) {
	    for(j=0; j<VSIZE; j+=i) {
		b[j+i-1] = b[j+i/2-1] + b[j+i-1];
		}
	}
	
	b[VSIZE-1]=0;
	    
//#pragma omp parallel for private(i)
	for(i=(log_N-1); i>=0; i--) {
		stride2=(2<<i)-1;
		stride=(1<<i)-1;
		step=stride2+1;
#pragma omp parallel for private(temp)
		for(j=0; j<VSIZE; j+=step) {
		    //temp=b[j+(int)(pow(2, i))-1];
		    temp=b[j+stride];
		    //b[j+(int)(pow(2, i))-1] = b[j+(int)(pow(2, i+1))-1];
		    b[j+stride] = b[j+stride2];
		    
		    //b[j+(int)(pow(2, i+1))-1] = temp+b[j+(int)(pow(2, i+1))-1];
		    b[j+stride2] = temp + b[j+stride2];
		}
	}
}

