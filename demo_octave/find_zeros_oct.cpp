#include <octave/oct.h>
#include <vector>

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

DEFUN_DLD (find_zeros_oct, args, nargout,
           "Returns the indexes of zeros (seen as local minima)")
{
    octave_value_list retval(2);

    if(args.length() != 1) {
        octave_stdout <<"Wrong argument number: " \
                       << "find_zeros_oct expects a tfrq map" << std::endl;
        
        retval(0) = octave_value(Matrix());
        retval(1) = octave_value(Matrix());
    } else if(nargout != 2) {
        octave_stdout << "find_zeros_oct returns two values." << std::endl;
    } else {
        int square_size = SQUARE_SIZE;
        int half_square = floor(square_size/2);
        
        NDArray tfrq_map = args(0).array_value ();
        int line_nb = tfrq_map.dims()(0);
        int row_nb = tfrq_map.dims()(1);

        std::vector<int> Is = std::vector<int>(), Js = std::vector<int>();

        for(int i = half_square; i < line_nb - half_square; i++) {
            for(int j = half_square; j < row_nb - half_square; j++) {
                double a = tfrq_map(i,j);
                double mini = a;
                for(int u = 0; u < square_size; u++) {
                    for(int v = 0; v < square_size; v++) {
                        mini = min(mini, 
                                   tfrq_map(i-half_square+u, j-half_square+v));
                        if(mini < a) break;
                    }
                    if(mini < a) break;
                }

                /* We have found a minimum, so store its position */
                if(mini == a) {
                    Is.push_back(i);
                    Js.push_back(j);
                }
            }
        }

        /* Output : the two arrays storing the positions of the minima */

        Matrix output_Is(1,Is.size());
        Matrix output_Js(1,Js.size());
        
        for(int k = 0; k < Is.size(); k++) {
            output_Is(0,k) = Is[k]+1;
            output_Js(0,k) = Js[k]+1;
        }

        retval(0) = octave_value(output_Is);
        retval(1) = octave_value(output_Js);
    }

    return retval;
}
