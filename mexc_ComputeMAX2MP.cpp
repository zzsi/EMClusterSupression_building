/*
 * matching pursuit on SUM2 maps
 *
 *    [activations] = mexc_ComputeMAX2MP( S2Map, templateRadius );
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <vector>
using namespace std;

#define ABS(x) ((x)>0? (x):(-(x)))
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))
#define PI 3.1415926
#define NEGMAX -1e10

 
/* variable declaration */
const float** S2Map;            /* SUM2 maps for multiple templates. All SUM2 maps are of the same size. */
int height, width;              /* size of S2Map maps */
int supressionRadius;        /* radius of surround supression */
int nTemplate;                      /* number of S2 templates */
vector<float> activations;			/* row, col, and iTemplate that corresponds to activated templates */
int S2Thres; 					/* cut-off value for S2 score */

void Compute()
{
    int x, y;
    float maxResponse, r;
    int bestX, bestY;
    int jT;
    int bestLocationshift, bestTemplate;
	int i, ind;
	bool* explained;
	bool found;
	float* maxOverTemplates;
	int* templateTrace;
	
	/* compute a single S2 map by maximizing over all S2 maps at corresponding pixels */
	maxOverTemplates = (float*)mxCalloc( height*width, sizeof(float) );
	templateTrace = (int*)mxCalloc( height*width, sizeof(int) );
	jT = 0;
	for (ind=0; ind<height*width; ++ind)
	{
		maxOverTemplates[ind] = S2Map[jT][ind];
		templateTrace[ind] = jT;
	}
	for( jT = 1; jT < nTemplate; ++jT )
	{
		for (ind=0; ind<height*width; ++ind)
		{
			r = S2Map[jT][ind];
			if( r > maxOverTemplates[ind] )
			{
				maxOverTemplates[ind] = r;
				templateTrace[ind] = jT;
			}
		}
	}
	
	activations.clear();
	explained = (bool*)mxCalloc( height*width, sizeof(bool) );
	for( i = 0; i < height * width; ++i )
	{
		explained[i] = false;
	}
	for( y = 0; y < width; ++y )
	{
		for( x = 0; x < supressionRadius; ++x )
		{
			explained[x+y*height] = true;
		}
		for( x = height-supressionRadius; x < height; ++x )
		{
			explained[x+y*height] = true;
		}
	}
	for( y = 0; y < supressionRadius; ++y )
	{
		for( x = 0; x < height; ++x )
		{
			explained[x+y*height] = true;
		}
	}
	for( y = width-supressionRadius; y < width; ++y )
	{
		for( x = 0; x < height; ++x )
		{
			explained[x+y*height] = true;
		}
	}
	
	found = true;
	while ( found )
	{
		found = false;
		/* find the global maximum of the SUM2 maps */
		maxResponse = NEGMAX;
		for (y=0; y<width; y++)
		{
			for (x=0; x<height; x++)
			{
				ind = x+y*height;
				r = maxOverTemplates[ind];
				if (!(explained[ind]) && r>maxResponse )
				{
					maxResponse = r;
					bestX = x; bestY = y;
					bestTemplate = templateTrace[ind];
					found = true;
				}				
			}
		}
		if(!found || maxResponse < S2Thres)
		{
			break;
		}
		/* record it */
		activations.push_back( bestX );
		activations.push_back( bestY );
		activations.push_back( bestTemplate );
		activations.push_back( maxResponse );
		/* inhibition */
		for( x = bestX-supressionRadius; x < bestX+supressionRadius; ++x )
		{
			for( y = bestY-supressionRadius; y <  bestY+supressionRadius; ++y )
			{
				if ((x>=0)&&(x<height)&&(y>=0)&&(y<width))
				{
					explained[x+y*height] = true;
				}
			}
		}
    }
}


/* entry point */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int ind, i, dataDim, bytes_to_copy;
    const mxArray *f;
    const mxArray *pAS2Map;
    mxArray *outPA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;

    /*
	 * input variable 0: S2 maps
	 */
    pAS2Map = prhs[0];
    nTemplate = mxGetM(pAS2Map) * mxGetN(pAS2Map);
    S2Map = (const float**)mxCalloc( nTemplate, sizeof(*S2Map) );   /* SUM1 maps */
    for (i=0; i<nTemplate; ++i)
    {
        f = mxGetCell(pAS2Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        S2Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        height = mxGetM(f);    /* overwriting is ok, since it is constant, or we would only input one image */
        width = mxGetN(f);
    }
    
    /*
     * input variable 1: location shift radius
     */
    supressionRadius = (int)mxGetScalar(prhs[1]);
    
    /*
     * input variable 2: cut-off activation S2 score
     */
    S2Thres = (int)mxGetScalar(prhs[2]);

    Compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
	/*
	 * output variable 0: activations
	 */
	dimsOutput[1] = (int)(activations.size()/4); dimsOutput[0] = 4;
    outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    float* a = (float*)mxGetData(outPA);
	for( int i = 0; i < (int)activations.size(); ++i )
	{
		a[i] = activations[i];
	}
    plhs[0] = outPA;
} 
