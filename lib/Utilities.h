//
// v1.0.0
//

#ifndef EFGE_UTILITIES_H
#define EFGE_UTILITIES_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "Constants.h"

using namespace std;

int parse_arguments(int argc, char** argv, string &corpus_file, string &output_file, string &method_name,
                     int &window_size, int &negative_sample_size, int &dimension,
                     double &starting_alpha, double &min_alpha, double &decay_rate, int &num_iters,
                     vector <double> &optionalParams);


#endif //EFGE_UTILITIES_H
