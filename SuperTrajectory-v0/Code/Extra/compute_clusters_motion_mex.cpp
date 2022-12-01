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
# include "string.h"

double max(double a, double b) {
    double big;
    if(a>b)
        big = a;
    else
        big = b;
    return big;
}

double min(double a, double b) {
    double small;
    if(a<b)
        small = a;
    else
        small = b;
    return small;
}

int* findLabels(double* trLabels, int labelLen, int labelVal)
{
	// find the labels whose value equal to the labelVal
	int* labelFound;
	int count = 0,j = 0;
	// count the number
	for(int i = 0; i < labelLen; i++)
		if((int)trLabels[i] == labelVal)
			count++;
	// record the index
	labelFound = (int*)malloc(count*sizeof(int));
	for(int i = 0; i < labelLen; i++)
		if((int)trLabels[i] == labelVal)
			labelFound[j++] = i;
	labelFound[j] = -1;
	return labelFound;
}

// input: labelNum trLabels Atr
// output: motion

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]) {
    /* declare variables */
	double *Atr, *motion, *trLabels;
	unsigned int labelNum, *countLabel;
	int nAtr, nL;
	mxArray * motionArray;

    int    i, j, k, total;
    
    
     if (nrhs != 3) {
         mexErrMsgTxt("Not enough input arguments!! 13 Input Argument expected");
     }
    
	trLabels = mxGetPr(prhs[0]);
	nL = mxGetN(prhs[0]);
	labelNum = (int)mxGetScalar(prhs[1]);
	Atr = mxGetPr(prhs[2]);
	nAtr = (int)mxGetM(prhs[2]);	

    mexPrintf("\t\t %d labels \n",labelNum);
    
	int dims[2] = {labelNum, labelNum};
	if(nlhs>0){
		plhs[0] = mxCreateDoubleMatrix(labelNum*labelNum, 1, mxREAL);
		//plhs[1] = mxCreateNumericArray(1,dims,mxUINT32_CLASS,mxREAL);
	}else{
		 mexErrMsgTxt("Not enough output arguments!\n");
	}

	motion = (double*)mxGetPr(plhs[0]);
	//countLabel = (unsigned int*)mxGetPr(plhs[1]);

	for(i = 0; i<labelNum; i++)
	{
		for(j = i+1; j<labelNum; j++)
		{
			int* L = findLabels(trLabels, nL, i+1);
			int* R = findLabels(trLabels, nL, j+1);
			double tmpMotion = 0.0;
			int tmpCount = 0;

			for(int l = 0; L[l] != -1; l++)
			{
				int idxL = L[l];
				for(int r = 0; R[r] != -1; r++)
				{
					int idxR = R[r];
					tmpMotion += Atr[ idxR*nAtr+idxL ];
					/*if(Atr[ idxR*nAtr+idxL ])
						tmpCount++;*/
				}
			}
			motion[j*labelNum+i] = tmpMotion;
			/*if(tmpCount)
				countLabel[j*labelNum+i] = tmpCount;*/
		}
	}
	mexPrintf("\t\t Motion Computation Done! \n");
	
}
