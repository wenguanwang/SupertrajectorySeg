#include <matrix.h>
#include <mex.h>

//#define DEBUG_MODE


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	unsigned int frame, frames, height, width, elements, i, j, point, superpixel, superpixels, size, bins;
	double *image;
	unsigned int *superpixelMap, *sizes;
	float *sHist1,*sHist2,*sHist3;
// 		/** Assert that the input cell arrays are of the same length */
	frames = mxGetNumberOfElements( prhs[ 0 ] );
// 		if( frames != mxGetNumberOfElements( prhs[ 1 ] ) )
// 			mexErrMsgTxt( USAGE_NOTE );
// 
// 		/** Assert that input cell arrays are not empty */
	if( frames == 0 )
		mexErrMsgTxt( "THe number of frames is 0!\n" );
// 		
// 		/** Assert that corresponding cell contents have the same number of elements */
    for( frame = 0; frame < 1; frame++ )
    {
        height = mxGetM( mxGetCell( prhs[ 0 ], frame ) );
		width = mxGetN( mxGetCell( prhs[ 0 ], frame ) ) / 3;
		elements = height * width;
	}

	superpixels = ( unsigned int )( ( double * )mxGetData( prhs[ 2 ] ) )[ 0 ];
    bins = ( unsigned int )( ( double * )mxGetData( prhs[ 3 ] ) )[ 0 ];
	#ifdef DEBUG_MODE
	mexPrintf( "getSuperpixelStats: Number of frames: %i\n", frames );
	mexPrintf( "getSuperpixelStats: Frame size: %ix%i\n", height, width );
	mexPrintf( "getSuperpixelStats: Total number of superpixels: %i\n", superpixels );
	mexEvalString( "pause(0.001)" );
	#endif
	
	#ifdef DEBUG_MODE
	mexPrintf( "getSuperpixelStats: Allocating memory for outputs...\n" );
	mexEvalString( "pause(0.001)" );
	#endif
	mxArray *sizesMxArray = mxCreateNumericMatrix( superpixels, 1, mxUINT32_CLASS, mxREAL );
	mxArray *histMxArrayR = mxCreateNumericMatrix( superpixels, bins, mxSINGLE_CLASS, mxREAL );
	mxArray *histMxArrayG = mxCreateNumericMatrix( superpixels, bins, mxSINGLE_CLASS, mxREAL );
    mxArray *histMxArrayB = mxCreateNumericMatrix( superpixels, bins, mxSINGLE_CLASS, mxREAL );
    
	sizes = ( unsigned int * )mxGetData( sizesMxArray );
	sHist1 = ( float * )mxGetData( histMxArrayR );
	sHist2 = ( float * )mxGetData( histMxArrayG );
    sHist3 = ( float * )mxGetData( histMxArrayB );
    
	#ifdef DEBUG_MODE
	mexPrintf( "getSuperpixelStats: Computing superpixel stats...\n" );
	mexEvalString( "pause(0.001)" );
	#endif
	for( frame = 0; frame < frames; frame++ )
	{
		image = ( double * )( mxGetData ( mxGetCell( prhs[ 0 ], frame ) ) );
		superpixelMap = ( unsigned int * )( mxGetData( mxGetCell( prhs[ 1 ], frame ) ) );
		
		for( i = 0; i < height; i++ )
		{
			for( j = 0; j < width; j++ )
			{
				point = j * height + i;
				superpixel = superpixelMap[ point ] - 1;
// 				
// 				if( superpixel < 0 && superpixel >= superpixels )
// 					mexErrMsgTxt( "getSuperpixelStats: Superpixel labels outside given range" );
// 				
				sHist1[ ( int )(superpixel + superpixels*image[ point ])] += 1;
				sHist2[ ( int )(superpixel + superpixels*image[ point + elements])] += 1;
				sHist3[ ( int )(superpixel + superpixels*image[ point + elements * 2])] += 1;
				sizes[ superpixel ]++;
			}
		}
	}
	
	for( superpixel = 0; superpixel < superpixels; superpixel++ )
	{
		size = ( float )sizes[ superpixel ];
        for( i = 0; i<bins; i++)
        {
            if( size > 0 )
            {
                sHist1[ superpixel + i * superpixels] /= size;
                sHist2[ superpixel + i * superpixels] /= size;
                sHist3[ superpixel + i * superpixels] /= size;
            }
        }
        
	}
	
	plhs[ 0 ] = histMxArrayR;
    plhs[ 1 ] = histMxArrayG;
    plhs[ 2 ] = histMxArrayB;
	
}

