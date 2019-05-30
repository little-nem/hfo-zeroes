#include <cmath>
#include <vector>

#include "mex.h"

#define min(x,y) (x < y ? x : y)
#define max(x,y) (x < y ? y : x)

/* Size of the region considered to check for local minimality, should be odd */
#define SQUARE_SIZE 15

/* 
 * The function find_zeros_mex operates on a surface and computes the
 * coordinates of local minima it encounters on the surface. It outputs the
 * index of such minima in two arrays, one for the row indexes and the other
 * the column indexes. To check for local minimality, the function uses a
 * square frame of size SQUARE_SIZE * SQUARE_SIZE centered on the candidate
 * and decide that the candidate is indeed a minimum when it is minimum amongst 
 * the frame.
 * 
 * Input:
 *  A matrix of real numbers representing the surface we want to study.
 * 
 * Outputs:
 *  Two arrays X and Y of integers of the same length, such as for every i,
 *  surface[X[i]][Y[i]] is a minimum of the surface.
 */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 1) {
        mexPrintf("Wrong argument: find_zeros expects a tfrq map\n");
        for(int i = 0; i < nlhs; i++) {
            plhs[i] = mxCreateDoubleMatrix (0, 0, mxREAL);
        }
    } else if(nlhs != 2) {
        mexPrintf("Two return values expected...\n");
    } else {
        int square_size = SQUARE_SIZE;
        int half_square = floor(square_size/2);

        int line_nb = mxGetM(prhs[0]);
        int row_nb = mxGetN(prhs[0]);
        
        std::vector<int> Is = std::vector<int>(), Js = std::vector<int>();

        double *tfrq_map = mxGetPr(prhs[0]);

        /* Search for local minima, using a square grid of size SQUARE_SIZE.
         * Could be sped up by taking advantage of the information aquired in
         * previous runs (ie. by moving the grid with bigger steps when 
         * possible), but the speed is already fine. */

        for(int i = half_square; i < line_nb - half_square; i++) {
            for(int j = half_square; j < row_nb - half_square; j++) {
                double a = tfrq_map[j*line_nb + i];
                double mini = a;
                for(int u = 0; u < square_size; u++) {
                    for(int v = 0; v < square_size; v++) {
                        int index = (j-half_square+v)*line_nb + i-half_square+u;
                        mini = min(mini, tfrq_map[index]);
                        if(mini < a) break;
                    }
                    if(mini < a) break;
                }

                if(mini == a) {
                    Is.push_back(i);
                    Js.push_back(j);
                }
            }
        }

        mwSize *dims = (mwSize*)mxMalloc(2 * sizeof(mwSize));
        dims[0] = 1; dims[1] = Is.size();

        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); 

        double *Is_output = (double*)mxGetPr(plhs[0]);
        double *Js_output = (double*)mxGetPr(plhs[1]);

        for(int k; k < Is.size(); k++) {
           Is_output[k] = Is[k]+1;
           Js_output[k] = Js[k]+1;
        }

        mxFree(dims);
    }
}
