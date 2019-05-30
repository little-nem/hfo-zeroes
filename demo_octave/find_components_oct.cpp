#include <octave/oct.h>
#include <vector>
#include <stack>
#include <utility>
#include <map>
#include <limits>

#define min(x,y) (x < y ? x : y)
#define max(x,y) (x < y ? y : x)

#define point(x,y) std::make_pair(x,y)
#define pair_of_points(x,y,u,v) std::make_pair(point(x,y),point(u,v))
#define leftmost_index(t, i, j) (t[i] < t[j] ? i : j)
#define rightmost_index(t, i, j) (t[i] < t[j] ? j : i)

/* 
 * The function find_components_oct operates on a set of triangles and computes
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

DEFUN_DLD (find_components_oct, args, nargout,
           "Return the connected components of the selected triangulation")
{
    octave_value_list retval(2);

    if(args.length() != 2) {
        octave_stdout <<"find_components_oct expects a two coords lists input."\
                                                               << std::endl;
        
        retval(0) = octave_value(Matrix());
        retval(1) = octave_value(Matrix());
    } else if(nargout != 2) {
        octave_stdout << "find_components_oct returns two values." << std::endl;
    } else {
        NDArray X = args(0).array_value ();
        NDArray Y = args(1).array_value ();

        int triangle_number = X.dims()(1);

        std::vector<int> components(triangle_number);
        std::vector<bool> seen(triangle_number);

        std::map<std::pair<std::pair<double,double>,std::pair<double,double>>, 
                                            std::vector<int>> edge_to_triangle;

        /* Build the graph of adjacency of the triangles 
         * We use a "natural" ordering to store the coordinates of the 
         * extremal points of a segment to ensure that there is a single 
         * way to store each segment (from left to right) */

        for(int t = 0; t < triangle_number; t++) {
            int left, right;
            float x[3], y[3];

            for(int k = 0; k < 3; k++) {
                x[k] = X(k, t);
                y[k] = Y(k, t);
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
                            x[k] = X(k, current_triangle);
                            y[k] = Y(k, current_triangle);
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
                                current_component.push(neighbour);
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
        
        std::vector<std::vector<double>> caract;
        caract = std::vector<std::vector<double>>(component_nb, 
               std::vector<double>(4, std::numeric_limits<double>::infinity()));

        for(int c = 0; c < component_nb; c++) {
            caract[c][1] *= -1;
            caract[c][3] *= -1;
        }


        for(int t = 0; t < triangle_number; t++) {
            float x[3], y[3];
            float min_freq, max_freq, min_time, max_time;

            for(int k = 0; k < 3; k++) {
                x[k] = X(k, t);
                y[k] = Y(k, t);
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

        Matrix output_comp(1,triangle_number);
        Matrix output_caract(4,component_nb);

        /* First ouput: we label each triangle with the id of the component 
         * it belongs to */
        for(int t = 0; t < triangle_number; t++) {
            output_comp(0,t) = components[t];
        }

        /* Second output: we store 4 attributes for each component */
        for(int c = 0; c < component_nb; c++) {
            for(int k = 0; k < 4; k++) {
                output_caract(k,c) = caract[c][k];
                output_caract(k,c) = caract[c][k];
                output_caract(k,c) = caract[c][k];
                output_caract(k,c) = caract[c][k];
            }
        }

        retval(0) = octave_value(output_comp);
        retval(1) = octave_value(output_caract);
    }

    return retval;
}
