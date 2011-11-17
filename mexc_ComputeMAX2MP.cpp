/*
 * Matching pursuit on SUM2 maps of a single image (at multiple resolutions).
 *
 *    [activations] = mexc_ComputeMAX2MP( S2Map, templateRadius, S2Thres );
 *
 *  In matlab, S2Map is defined as:  SUM2map = cell(numTemplate,numResolution);
 *  The S2 templates are assumed to have the same size. 
 *	If the image has multiple resolutions, make sure it is arranged from low resolution to high resolution.
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
int *height, *width;              /* sizes of S2Map maps */
int supressionRadius;        /* radius of surround supression */
int nTemplate;                      /* number of S2 templates */
int nResolution;					/* number of image resolutions */
vector<float> activations;			/* row, col, and iTemplate that corresponds to activated templates */
int S2Thres; 					/* cut-off value for S2 score */

void Compute()
{
    int x, y;
	float x2, y2;
    float maxResponse, r, scaling;
    int bestX, bestY;
    int iT, iR; /* index for template and reolution */
    int bestLocationshift, bestTemplate, bestResolution;
	int i, ind;
	bool* explained;
	bool found;
	float* maxOverTemplates;
	int* templateTrace, *resolutionTrace;
	
	/* compute a single MAX2 map and ARGMAX map by maximizing over all S2 maps at corresponding pixels */
	maxOverTemplates = (float*)mxCalloc( height[nResolution-1]*width[nResolution-1], sizeof(float) ); /* single the image is at multiple resolutions, use the finest (largest) resolution */
	templateTrace = (int*)mxCalloc( height[nResolution-1]*width[nResolution-1], sizeof(int) );
	resolutionTrace = (int*)mxCalloc( height[nResolution-1]*width[nResolution-1], sizeof(int) );
	
	ind = 0;
	for (y=0; y<width[nResolution-1]; ++y)
	{
		for (x=0; x<height[nResolution-1]; ++x)
		{
			maxOverTemplates[ind] = NEGMAX;
			templateTrace[ind] = -1;
			resolutionTrace[ind] = -1;
			++ind;
		}
	}
	
	i = 0; /* index for SUM2 map */
	for (iR=0; iR<nResolution; ++iR)
	{
		scaling = (float)height[iR] / (float)height[nResolution-1]; /* scaling < 1 */
		for (iT=0; iT<nTemplate; ++iT)
		{
			ind = 0; /* index for pixel inside the largest SUM2 map */
			x2 = 0; y2 = 0;
			for (y=0; y<width[nResolution-1]; ++y )
			{
				for (x=0; x<height[nResolution-1]; ++x )
				{
					r = S2Map[i][ (int)floor(x2)+(int)floor(y2)*height[iR] ];
					if( r > maxOverTemplates[ind] )
					{
						maxOverTemplates[ind] = r;
						templateTrace[ind] = iT;
						resolutionTrace[ind] = iR;
					}
					ind++;
					x2 += scaling;
				}
				y2 += scaling; /* (x2,y2): pixel location in the SUM2 map at resolution iR */
			}
			++i;
		}
	}
	
	activations.clear();
	explained = (bool*)mxCalloc( height[nResolution-1]*width[nResolution-1], sizeof(bool) );
	for( i = 0; i < height[nResolution-1] * width[nResolution-1]; ++i )
	{
		explained[i] = false;
	}
	/* deal with boundary condition */
	/*
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
	*/
	
	found = true;
	while ( found )
	{
		found = false;
		/* find the global maximum of the SUM2 maps */
		maxResponse = NEGMAX;
		ind = 0;
		for (y=0; y<width[nResolution-1]; y++)
		{
			for (x=0; x<height[nResolution-1]; x++)
			{
				r = maxOverTemplates[ind];
				if (!(explained[ind]) && r>maxResponse )
				{
					bestResolution = resolutionTrace[ind];
					scaling = (float)height[bestResolution] / (float)height[nResolution-1];
					maxResponse = r;
					bestX = (int)floor(x*scaling+.5); bestY = (int)floor(x*scaling+.5); /* (bestX,bestY) is the location in the corresponding resolution */
					bestTemplate = templateTrace[ind];
					found = true;
				}
				ind++;
			}
		}
		if(!found || maxResponse < S2Thres)
		{
			break;
		}
		/* record it */
		activations.push_back( (float)bestX );
		activations.push_back( (float)bestY );
		activations.push_back( (float)bestResolution );
		activations.push_back( (float)bestTemplate );
		activations.push_back( maxResponse );
		/* inhibition */
		scaling = (float)height[nResolution-1] / (float)height[bestResolution];
		for( y = (int)floor(bestY-supressionRadius*scaling); y <  bestY+supressionRadius*scaling; ++y )
		{
			for( x = (int)floor(bestX-supressionRadius*scaling); x < bestX+supressionRadius*scaling; ++x )
			{
				if ((x>=0)&&(x<height[nResolution-1])&&(y>=0)&&(y<width[nResolution-1]))
				{
					explained[x+y*height[nResolution-1]] = true;
				}
			}
		}
    }
}


/* entry point */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int ind, i, j, dataDim, bytes_to_copy;
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
    nTemplate = (int)mxGetM(pAS2Map);
	nResolution = (int)mxGetN(pAS2Map);
    S2Map = (const float**)mxCalloc( nTemplate*nResolution, sizeof(*S2Map) );   /* SUM2 maps */
	ind = 0;
	for (i=0; i < nResolution; ++i)
    {
		for (j=0; j<nTemplate; ++j)
		{
			f = mxGetCell(pAS2Map, ind);
			datatype = mxGetClassID(f);
			if (datatype != mxSINGLE_CLASS)
				mexErrMsgTxt("warning !! single precision required.");
			S2Map[ind] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
			++ind;
		}
		height[i] = mxGetM(f);
        width[i] = mxGetN(f);
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
	dimsOutput[1] = (int)(activations.size()/5); dimsOutput[0] = 5;
    outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    float* a = (float*)mxGetData(outPA);
	for( int i = 0; i < (int)activations.size(); ++i )
	{
		a[i] = activations[i];
	}
    plhs[0] = outPA;
} 
