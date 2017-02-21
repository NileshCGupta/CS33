#include "func.h"
#include "omp.h"

void Func1(int c[][TSIZE], int a[][TSIZE], int b[][TSIZE])
{
    int i, j, k, r;
#pragma omp parallel for private(k, r, j)
    for (i=0; i<TSIZE; i++) {
	for (k=0; k<TSIZE; k++) {
	    r = a[i][k];
	    for (j=0; j<TSIZE; j++) {

		c[i][j]+=r*b[k][j];
	    }
	}
    }
}
void Func2(int d[][MSIZE], int c[][MSIZE])
{
    int i,j,i1,j1;
    int B = MSIZE/500;
#pragma omp parallel for private(j,i1,j1)
    for (i = 0; i < MSIZE; i+=B)
	for (j = 0; j < MSIZE; j+=B)
	    for (i1 = i; i1 < i+B; i1++)
		for (j1 = j; j1 < j+B; j1++)
		    d[i1][j1]=c[j1][i1];
}

void Func3(int z[][MSIZE], int d[][MSIZE])
{

    int y, x;
    int near = 2;  // The weight of neighbor
    int itself = 84; // The weight of center
#pragma omp parallel for collapse(2)
    for (y=0; y<MSIZE; y++) {
	for (x=0; x<MSIZE; x++) {
	    if (y==0) {
		if (x==0) {
		    z[y][x] = near * d[y][x] +
			  near * d[y][x+1] +
			  near * d[y+1][x] +
			  near * d[y+1][x+1] +
			  near * d[y][x] +
			  near * d[y][x+1] +
			  near * d[y][x] +
			  near * d[y+1][x] +
			itself * d[y][x];
		} else if (x==MSIZE-1) {
		    z[y][x] = near * d[y][x-1] +
			  near * d[y][x] +
			  near * d[y+1][x-1] +
			  near * d[y+1][x] +
			  near * d[y][x-1] +
			  near * d[y][x] +
			  near * d[y][x] +
			  near * d[y+1][x] +
			itself * d[y][x];
		} else {
		    z[y][x] = near * d[y][x-1] +
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
		    z[y][x] = near * d[y-1][x] +
			  near * d[y-1][x+1] +
			  near * d[y][x] +
			  near * d[y][x+1] +
			  near * d[y][x] +
			  near * d[y][x+1] +
			  near * d[y-1][x] +
			  near * d[y][x] +
			itself * d[y][x];
		} else if (x==MSIZE-1) {
		    z[y][x] = near * d[y-1][x-1] +
			  near * d[y-1][x] +
			  near * d[y][x-1] +
			  near * d[y][x] +
			  near * d[y][x-1] +
			  near * d[y][x] +
			  near * d[y-1][x] +
			  near * d[y][x] +
			itself * d[y][x];
		} else {
		    z[y][x] = near * d[y-1][x-1] +
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
		    z[y][x] = near * d[y-1][x] +
			  near * d[y-1][x+1] +
			  near * d[y+1][x] +
			  near * d[y+1][x+1] +
			  near * d[y][x] +
			  near * d[y][x+1] +
			  near * d[y-1][x] +
			  near * d[y+1][x] +
			itself * d[y][x];
		} else if (x==MSIZE-1) {
		    z[y][x] = near * d[y-1][x-1] +
			  near * d[y-1][x] +
			  near * d[y+1][x-1] +
			  near * d[y+1][x] +
			  near * d[y][x-1] +
			  near * d[y][x] +
			  near * d[y-1][x] +
			  near * d[y+1][x] +
			itself * d[y][x];
		} else {
		    z[y][x] = near * d[y-1][x-1] +
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
}
/*
void Func3(int z[][MSIZE], int d[][MSIZE])
{
  int y, x;
  int near = 2;  // The weight of neighbor
  int itself = 84; // The weight of center
  //4 corner cases
  //y == 0 && x == 0
  #pragma omp parallel
  z[0][0] = near * d[0][0] +
        near * d[0][1] +
        near * d[1][0] +
        near * d[1][1] +
        near * d[0][0] +
        near * d[0][1] +
        near * d[0][0] +
        near * d[1][0] +
    itself * d[0][0];
 

   #pragma omp parallel
  // y == 0 && x == MSIZE -1

      z[0][MSIZE-1] = near * d[0][MSIZE-1-1] +
          near * d[0][MSIZE-1] +
	      near * d[0+1][MSIZE-1-1] +
	          near * d[0+1][MSIZE-1] +
		      near * d[0][MSIZE-1-1] +
		          near * d[0][MSIZE-1] +
			      near * d[0][MSIZE-1] +
			          near * d[0+1][MSIZE-1] +
				  itself * d[0][MSIZE-1];

      // y == MSIZE - 1 X == 0
       #pragma omp parallel
      z[MSIZE - 1][0] = near * d[MSIZE - 1-1][0] +
          near * d[MSIZE - 1-1][0+1] +
	      near * d[MSIZE - 1][0] +
	          near * d[MSIZE - 1][0+1] +
		      near * d[MSIZE - 1][0] +
		          near * d[MSIZE - 1][0+1] +
			      near * d[MSIZE - 1-1][0] +
			          near * d[MSIZE - 1][0] +
				  itself * d[MSIZE - 1][0];

      //  y == MSIZE-1 && X == MSIZE-1
       #pragma omp parallel
      z[MSIZE-1][MSIZE-1] = near * d[MSIZE-1-1][MSIZE-1-1] +
          near * d[MSIZE-1-1][MSIZE-1] +
	      near * d[MSIZE-1][MSIZE-1-1] +
	          near * d[MSIZE-1][MSIZE-1] +
		      near * d[MSIZE-1][MSIZE-1-1] +
		          near * d[MSIZE-1][MSIZE-1] +
			      near * d[MSIZE-1-1][MSIZE-1] +
			          near * d[MSIZE-1][MSIZE-1] +
				  itself * d[MSIZE-1][MSIZE-1];
#pragma omp parallel for
      for (x = 1 ; x < MSIZE-1 ; x ++)
      {
        z[0][x] = near * d[0][x-1] +
	            near * d[0][x+1] +
		                near * d[0+1][x-1] +
				            near * d[0+1][x+1] +
					                near * d[0][x-1] +
							            near * d[0][x+1] +
								                near * d[0][x] +
										            near * d[0+1][x] +
											        itself * d[0][x];
												}
#pragma omp parallel for
      for (x = 1; x < MSIZE-1; x++)
      {
        z[MSIZE-1][x] = near * d[MSIZE-1-1][x-1] +
	            near * d[MSIZE-1-1][x+1] +
		                near * d[MSIZE-1][x-1] +
				            near * d[MSIZE-1][x+1] +
					                near * d[MSIZE-1][x-1] +
							            near * d[MSIZE-1][x+1] +
								                near * d[MSIZE-1-1][x] +
										            near * d[MSIZE-1][x] +
											        itself * d[MSIZE-1][x];
												}
#pragma omp parallel for
      for ( y = 1; y < MSIZE-1; y++)
      {
        z[y][0] = near * d[y-1][0] +
	            near * d[y-1][0+1] +
		                near * d[y+1][0] +
				            near * d[y+1][0+1] +
					                near * d[y][0] +
							            near * d[y][0+1] +
								                near * d[y-1][0] +
	            near * d[y+1][0] +
		    itself * d[y][0];
		    }
#pragma omp parallel for
for ( y = 1 ; y < MSIZE-1 ; y++)
{
  z[y][MSIZE-1] = near * d[y-1][MSIZE-1-1] +
              near * d[y-1][MSIZE-1] +
              near * d[y+1][MSIZE-1-1] +
              near * d[y+1][MSIZE-1] +
              near * d[y][MSIZE-1-1] +
              near * d[y][MSIZE-1] +
              near * d[y-1][MSIZE-1] +
              near * d[y+1][MSIZE-1] +
	      itself * d[y][MSIZE-1];
	      }

	      #pragma omp parallel for collapse(2)
	      for (y=1; y<MSIZE-1; y++)
	      {
	      for (x=1; x<MSIZE-1; x++)
	      {
      z[y][x] = near * d[y-1][x-1] +
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
		  #pragma omp parallel for collapse(2)
		  for ( x = 0 ; x < MSIZE; x ++)
		  {
		  for ( y = 0; y < MSIZE; y++)
		  {
		  z[y][x]/=100;
		  }
		  }
		  }
*/


void Func4(int b[], int a[])
{
    int chuck_size=MSIZE;
    int array_size=VSIZE/chuck_size;
    //int array_size = chuck_size;
    int chuck[chuck_size];
    int i, j;
    int prev;
    int curr;
#pragma omp parallel for private(i)
    for(j=0; j<chuck_size; j++) {
	b[j*array_size]=a[j*array_size];
	//#pragma omp parallel for
	for (i=1; i<array_size; i++)
	{
	    // int prev = b[j*array_size+i-1];
	    // int curr = a[j*array_size+i];
	    //maybe use flush since we're changing value at memory level?


	    b[j*array_size+i]=b[j*array_size+i-1]+a[j*array_size+i];
	}
	chuck[j]=b[(j+1)*array_size-1];
    }
    //#pragma omp parallel for
    for(j=1; j<chuck_size; j++)
    {
	// prev = chuck[j-1];
	//#pragma omp atomic
	//chuck[j] += prev;

	chuck[j]=chuck[j-1]+chuck[j];
    }
    /*
      for(j=0; j<chuck_size-1; j+=8) {
        chuck[j+1]=chuck[j]+chuck[j+1];
	  chuck[j+2]=chuck[j+1]+chuck[j+2];
	    chuck[j+3]=chuck[j+2]+chuck[j+3];
	      chuck[j+4]=chuck[j+3]+chuck[j+4];
	        chuck[j+5]=chuck[j+4]+chuck[j+5];
		  chuck[j+6]=chuck[j+5]+chuck[j+6];
		    chuck[j+7]=chuck[j+6]+chuck[j+7];
		      chuck[j+8]=chuck[j+7]+chuck[j+8];
		      }
    */
#pragma omp parallel for private(i)
    for(j=1; j<chuck_size; j++) {
	//#pragma omp parallel for
	for (i=0; i<VSIZE/chuck_size; i++) {
	    b[j*array_size+i]+=chuck[j-1]/(j+1);
	}
    }
}


/*
void Func4(int b[], int a[])
{
  int chuck_size=MSIZE;
  int array_size=chuck_size;
  int chuck[chuck_size];
  int i, j,temp;
  #pragma omp parallel for private(i)
  for(j=0; j<chuck_size; j++) {
    b[j*array_size]=a[j*array_size];
         #pragma omp parallel for
    for (i=1; i<chuck_size; i++) {
      b[j * array_size+i]=b[j * array_size+i-1]+a[j * array_size+i];
    }
    chuck[j]=b[(j+1)*array_size-1];
  }

  for(j=0; j<chuck_size-1; j+=8) {
    chuck[j+1]=chuck[j]+chuck[j+1];
    chuck[j+2]=chuck[j+1]+chuck[j+2];
    chuck[j+3]=chuck[j+2]+chuck[j+3];
    chuck[j+4]=chuck[j+3]+chuck[j+4];
    chuck[j+5]=chuck[j+4]+chuck[j+5];
    chuck[j+6]=chuck[j+5]+chuck[j+6];
    chuck[j+7]=chuck[j+6]+chuck[j+7];
    chuck[j+8]=chuck[j+7]+chuck[j+8];
  }

#pragma omp parallel for private(i)
  for(j=1; j<chuck_size; j++) {
       #pragma omp parallel for
    for (i=0; i<chuck_size; i++) {
      b[j * array_size+i]+=chuck[j-1]/(j+1);
    }
  }
}
*/
#undef block
void Func5(int b[], int a[])
{
    int i=0, j, j1, stride, stride2, step;
    int temp;
    long log_N=log2(VSIZE);

    # define block 64
#pragma omp parallel shared(a,b) private(j,j1)
#pragma omp for schedule(guided)
    for(j=0; j<VSIZE; j+=block)
    {
	for (j1 = j; j1 < j + block; j1 += 2)
	{
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
    for(i=(log_N-1); i>=0; i--) {
	stride2=pow(2,i+1)-1;
	stride=pow(2,i)-1;
	step=stride2+1;
#pragma omp parallel for private(temp)
	for(j=0; j<VSIZE; j+= step) {
	    int js1 = j + stride;
	    int js2 = j + stride2;
	    temp=b[js1];
	    b[js1] = b[js2];
	    b[js2] = temp+b[js2];
	}
    }
}
