#include <cmath>
#include <vector>
#include <stack>
#include <utility>
#include <map>
#include <limits>

#include "mex.h"

#define min(x,y) (x < y ? x : y)
#define max(x,y) (x < y ? y : x)

#define point(x,y) std::make_pair(x,y)
#define pair_of_points(x,y,u,v) std::make_pair(point(x,y),point(u,v))
#define leftmost_index(t, i, j) (t[i] < t[j] ? i : j)
#define rightmost_index(t, i, j) (t[i] < t[j] ? j : i)

/* 
 * The function find_components_mex operates on a set of triangles and computes
 * the adjacency components of this set, ie. it labels triangles according to
 * their belonging to one of the component of the plane. It also computes some
 * information about each component: its highest, lowest, leftmost and rightmost
 * points.
 * 
 * Inputs:
 *  Two matrices of size 3×nb_triangles, the first one for the x coordinates 
 *  of the vertices of each triangle, columnwise, and the second one for
 *  their y coordinates.
 * 
 * Outputs:
 *  - first output: an array of integers, of size nb_triangles, where each
 *  triangle is associated to the index of the component it belongs to.
 *  - second output: an array of size 4×nb_components, each column carrying
 *  informations for one component; the data is columnwise organised as follows:
 *                [ min_time,
 *                  max_time, 
 *                  min_freq, 
 *                  max_freq ]
 * 
 * WARNING : There may be some issues with types (platform dependant?). The
 * following code assume that the input data can be loaded as an array of
 * floats. 
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs != 2) {
        mexPrintf("find_components_mex expects two coords lists as input\n");
        for(int i = 0; i < nlhs; i++) {
            plhs[i] = mxCreateDoubleMatrix (0, 0, mxREAL);
        }
    } else if(nlhs < 2) {
        mexPrintf("find_components_mex returns two values that should be kept");
    } else {
        float *X = (float*)mxGetData(prhs[0]);
        float *Y = (float*)mxGetData(prhs[1]);

        int triangle_number = mxGetN(prhs[0]);

        std::vector<int> components(triangle_number);
        std::vector<bool> seen(triangle_number);

        std::map<std::pair<std::pair<float,float>,std::pair<float,float> >, 
                                             std::vector<int> > edge_to_triangle;

        /* Build the graph of adjacency of the triangles 
         * We use a "natural" ordering to store the coordinates of the 
         * extremal points of a segment to ensure that there is a single 
         * way to store each segment (from left to right) */

        for(int t = 0; t < triangle_number; t++) {
            int left, right;
            float x[3], y[3];

            /* Read the 3 points of the triangle */
            for(int k = 0; k < 3; k++) {
                x[k] = X[3*t+k];
                y[k] = Y[3*t+k];
            }

            /* Associate the edges with triangles using our ordering */
            for(int k = 0; k < 3; k++) {
                left = leftmost_index(x, k, (k+1)%3);
                right = rightmost_index(x, k, (k+1)%3);

                edge_to_triangle[pair_of_points(x[left],y[left],
                                             x[right], y[right])].push_back(t);
                                             
            }
        }

        int component_nb = 0;
        std::pair<std::pair<float,float>, std::pair<float, float> > edge;

        /* Now we simply do a graph traversal to label each triangle according 
         * to the connex component it belongs to */
        for(int triangle = 0; triangle < triangle_number; triangle++) {
            if(!seen[triangle]) {
                component_nb++;

                std::stack<int> current_component;
                current_component.push(triangle);

                /* While we can find unexplored triangles in the component */
                while(!current_component.empty()) {
                    int current_triangle = current_component.top();
                    current_component.pop();
                    
                    if(!seen[current_triangle]) {
                        seen[current_triangle] = true;
                        components[current_triangle] = component_nb; 

                        int left, right;
                        float x[3], y[3];

                        for(int k = 0; k < 3; k++) {
                            x[k] = X[3*current_triangle+k];
                            y[k] = Y[3*current_triangle+k];
                        }

                        /* For each edge of the current triangle */
                        for(int k = 0; k < 3; k++) {
                            left = leftmost_index(x, k, (k+1)%3);
                            right = rightmost_index(x, k, (k+1)%3);

                            edge = pair_of_points(x[left], y[left],
                                                  x[right], y[right]);

                            /* We find the neighbours of the current edge */
                            int neighbour_nb = edge_to_triangle[edge].size();
                            for(int neighbour : edge_to_triangle[edge]) {
                            /*for(int neigh = 0; neigh = neighbour_nb; neigh++) {*/
                                current_component.push(neighbour);
                                //current_component.push(edge_to_triangle[edge][neigh]);
                            }
                        }
                    }
                }
            }
        }

        /* Output */
      
        /* We compute some additional information for each component : 
         * ...[0] : min_time, ...[1] : max_time,
         * ...[2] : min_freq, ...[3] : max_freq */

        std::vector<std::vector<double> > caract;

        /* Initialize it with some infinity so we can handle max/min easily */
        caract = std::vector<std::vector<double> >(component_nb, 
               std::vector<double>(4, std::numeric_limits<double>::infinity()));
        for(int c = 0; c < component_nb; c++) {
            caract[c][1] *= -1;
            caract[c][3] *= -1;
        }


        for(int t = 0; t < triangle_number; t++) {
            float x[3], y[3];
            float min_freq, max_freq, min_time, max_time;

            for(int k = 0; k < 3; k++) {
                x[k] = X[3*t+k];
                y[k] = Y[3*t+k];
            }
            
            int t_comp = components[t];
            min_freq = min(min(y[0], y[1]), y[2]);
            max_freq = max(max(y[0], y[1]), y[2]);
            min_time = min(min(x[0], x[1]), x[2]);
            max_time = max(max(x[0], x[1]), x[2]);

            caract[t_comp-1][0] = min(caract[t_comp-1][0], min_time);
            caract[t_comp-1][1] = max(caract[t_comp-1][1], max_time);
            caract[t_comp-1][2] = min(caract[t_comp-1][2], min_freq);
            caract[t_comp-1][3] = max(caract[t_comp-1][3], max_freq);
        }

        /* First output: each triangle is labelled with its component number */
        mwSize *dims = (mwSize *)mxMalloc(2 * sizeof(mwSize));
        dims[0] = 1; dims[1] = triangle_number;

        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        double *triangle_labels = (double*)mxGetPr(plhs[0]);
        
        for(int t = 0; t < triangle_number; t++) {
            triangle_labels[t] = components[t];
        }
      
        /* Second output: 4 attributes for each components */ 
        dims[0] = 4; dims[1] = component_nb;
        plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); 

        double *attributes = (double*)mxGetPr(plhs[1]);

        for(int c = 0; c < component_nb; c++) {
            for(int k = 0; k < 4; k++) {
                attributes[4*c + k] = caract[c][k];
            }
        }

        mxFree(dims);
    }
}
