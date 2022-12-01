#include <matrix.h>
#include <math.h>   
#include <mex.h>   
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{    
    const int  *dim_array;
    int number_of_dims;
    int i,t,x,y,change,k;
    int index;
    
    dim_array = mxGetDimensions(prhs[0]);
    number_of_dims = mxGetNumberOfDimensions(prhs[0]);
    int numARows = *dim_array; // h 
    int numACols = *(dim_array+1);// w 
    int numALays;
    if (number_of_dims ==2)
        numALays = 1;
    else
        numALays = *(dim_array+2);
    
    double *AllFprob = mxGetPr(prhs[0]);
    double *Alllabel = mxGetPr(prhs[1]);
	unsigned int superpixels = ( unsigned int )( ( double * )mxGetData( prhs[ 3 ] ) )[ 0 ];
    
    mxArray *massNumMxArray = mxCreateNumericMatrix( superpixels, 1, mxSINGLE_CLASS, mxREAL );
	mxArray *meanFprobMxArray = mxCreateNumericMatrix( superpixels, 1, mxSINGLE_CLASS, mxREAL );
//     float *labels = ( float * )mxGetData( massNumMxArray );
    float *snum = ( float * )mxGetData( massNumMxArray );
    float *slabels = ( float * )mxGetData( meanFprobMxArray );
//     for (x=0;x<3;x++){
        for(y=0;y<superpixels;y++){
//             labels[x*superpixels+y] = 0;
//             snum[x*superpixels+y] = 0;
            slabels[y] = 0;
            snum[y] = 0;
        }
//     }
    unsigned int *superpixelMap;
    for( t=0; t<numALays;t++ ){ 
        superpixelMap = ( unsigned int * )( mxGetData( mxGetCell( prhs[ 2 ], t ) ) );
        for( x=0; x<numACols;x++ ){
            for( y=0; y<numARows;y++ ){             
//                 index = AllFprob[numACols*numARows*t + numARows*x + y]*superpixels+ superpixelMap[numARows*x + y]-1;
                index = superpixelMap[numARows*x + y]-1;
                slabels[index] += AllFprob[numACols*numARows*t + numARows*x + y];
                snum[index] += 1;
            }
        }
    }
    
    for(y=0;y<superpixels;y++){
        slabels[y] = slabels[y]/(snum[y]+0.0001);
    }

    plhs[ 0 ] = meanFprobMxArray;
}

