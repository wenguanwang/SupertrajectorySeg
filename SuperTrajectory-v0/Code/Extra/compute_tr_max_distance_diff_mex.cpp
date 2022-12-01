// function A=compute_tr_max_measurement_diff_mex(X,T,tr_start, tr_end, ntr, sigma, pi, pj);
// Input:
// X (double)      : space time point x
// Y (double)      : space time point y
// T (double)         : space time point t
// npoint          : number of points
// Ux (double)     : x speed
// Uy (double)     : y speed
// tr_start(double)   : the space time point start index of all tracks
// tr_end(double)     : the space time point end index of all tracks
// ntr             : number of tracks
// aff_var         : affinity variance
// [pi, pj]        : index pair representation for MATLAB sparse matrix
// verbose         : if 2, display info

// Out:
//  A = affinity with max velocity difference and max distance difference at [pi,pj]

//  Weiyu Zhang, Jan 12 2012
//  Based on Katerina's code


# include "mex.h"
# include "math.h"

double max( double a, double b ) 
{
    double big;
    if( a > b )
        big = a;
    else
        big = b;
    return big;
}

double min( double a, double b ) 
{
    double small;
    if( a < b )
        small = a;
    else
        small = b;
    return small;
}

double compute_affinities( double *X,  double *T, double *tr_start, double *tr_end, 
	double sigma, int np, int nd, int tr1, int tr2 )
{
    
    double mmax,  mc;
    int st1c, st2c;
    double aff;
    
    int common_start = ( int )max( T[ ( int )tr_start[ tr1 ] ], T[ ( int )tr_start[ tr2 ] ] );
    int common_end   = ( int )min( T[ ( int )tr_end[ tr1 ] ], T[ ( int )tr_end[ tr2 ] ] );
    int common_len    = common_end - common_start + 1;
    
    mmax = 0;
    for( int frameindex = common_start; frameindex <= common_end; frameindex ++ )
	{
        st1c = ( int )tr_start[ tr1 ] + frameindex - ( int )T[ ( int )tr_start[ tr1 ] ];
        st2c = ( int )tr_start[ tr2 ] + frameindex - ( int )T[ ( int )tr_start[ tr2 ] ];
        mc = 0;
        for( int k = 0; k < nd; k ++ )
		{
            mc = mc + pow( X[ st1c + k * np ] - X[ st2c + k * np ], 2 );
		}
        mc = sqrt( mc );
        if( mc > mmax )
            mmax = mc;
    }
    
    aff = exp( - sigma * mmax );

    if ( aff < 0.01 ) 
		aff = 0; 

    return aff;
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[])
{
    /* declare variables */

    double *X, *T, *tr_start, *tr_end, sigma;
    int    ntr, np, nd, npairs;
    int    i, j, k, total;
    unsigned int *pi, *pj;
    double *w;
    mwIndex *ir, *jc;
    
    X = mxGetPr( prhs[ 0 ] );
    np = ( int )mxGetM( prhs[ 0 ] );
    nd = ( int )mxGetN( prhs[ 0 ] );
    mexPrintf( "np = %d , \t nd=%d \t", np, nd );
    T = mxGetPr( prhs[ 1 ] );
  
    tr_start = mxGetPr( prhs[ 2 ] );
    tr_end = mxGetPr( prhs[ 3 ] );
    ntr = ( int )mxGetScalar( prhs[ 4 ] );
    sigma = mxGetScalar( prhs[ 5 ] );
    
    
    pi = ( unsigned int * )mxGetData( prhs[ 6 ] );
    pj = ( unsigned int * )mxGetData( prhs[ 7 ] );
    npairs      = ( int ) mxGetScalar( prhs[ 8 ] );
    mexPrintf( "\t\t\t %d comparisons \n", npairs );
    
    if ( !mxIsUint32( prhs[ 6 ] ) | !mxIsUint32( prhs[ 7 ] ) )
	{
        mexErrMsgTxt( "Index pair shall be of type UINT32" );
    }
    
    if ( nlhs > 0 )
	{
        plhs[ 0 ] = mxCreateDoubleMatrix( npairs, 1, mxREAL );
    }
    if ( plhs[ 0 ] == NULL )
	{
        mexErrMsgTxt( "Not enough memory for the plhsput matrix" );
    }

    w = mxGetPr( plhs[ 0 ] );
    total = 0;
    for( k = 0; k < npairs; k ++ )
	{
		w[ total ] = compute_affinities( X, T, tr_start, tr_end, sigma, np, nd, pi[ k ], pj[ k ] );
        total = total + 1;
    }

	return;
}
