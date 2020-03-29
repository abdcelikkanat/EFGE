//
// v1.0.0
//
#include "Utilities.h"



int parse_arguments(int argc, char** argv, string &corpus_file, string &output_file, string &method_name,
                     int &window_size, int &negative_sample_size, int &dimension,
                     double &starting_alpha, double &min_alpha, double &decay_rate, int &num_iters,
                     vector <double> &optionalParams) {

    vector <string> parameter_names{"--help",
                              "--corpus", "--output", "--method",
                              "--window_size", "--neg_sample", "--dim",
                              "--starting_alpha", "--min_alpha", "--decay_rate", "--num_iters",
                              "--sigma"
    };

    string arg_name;
    stringstream help_msg, help_msg_required, help_msg_opt;
    vector <string> constraintValues = {"bern", "pois", "norm"};

    // Set the help message
    help_msg_required << "\nUsage: ./" << EFGEConstants::ProgramName;
    help_msg_required << " " << parameter_names[1] << " CORPUS_FILE "
                             << parameter_names[2] << " OUTPUT_FILE "
                             << parameter_names[3] << " METHOD_NAME{" << constraintValues[0] << ","
                             << constraintValues[1] << "," << constraintValues[2] << "}\n";
    help_msg_opt << "\nOptional parameters:\n";
    help_msg_opt << "\t[ " << parameter_names[4] << " (Default: " << window_size << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[5] << " (Default: " << negative_sample_size << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[6] << " (Default: " << dimension <<") ]\n";
    help_msg_opt << "\t[ " << parameter_names[7] << " (Default: " << starting_alpha << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[8] << " (Default: " << min_alpha << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[9] << " (Default: " << decay_rate << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[10] << " (Default: " << num_iters << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[11] << " (Default: " << optionalParams[0] << ") ]\n";
    help_msg_opt << "\t[ " << parameter_names[0] << ", -h ] Shows this message";

    help_msg << "" << help_msg_required.str() << help_msg_opt.str();

    // Read the argument values
    for(int i=1; i<argc; i=i+2) {

        arg_name.assign(argv[i]);

        if (arg_name.compare(parameter_names[1]) == 0) {
            corpus_file = argv[i + 1];
        } else if (arg_name.compare(parameter_names[2]) == 0) {
            output_file = argv[i + 1];
        } else if (arg_name.compare(parameter_names[3]) == 0) {
            method_name = argv[i + 1];
        } else if (arg_name.compare(parameter_names[4]) == 0) {
            window_size = stoi(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[5]) == 0) {
            negative_sample_size = stoi(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[6]) == 0) {
            dimension = stoi(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[7]) == 0) {
            starting_alpha = stod(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[8]) == 0) {
            min_alpha = stod(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[9]) == 0) {
            decay_rate = stod(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[10]) == 0) {
            num_iters = stoi(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[11]) == 0) {
            optionalParams[0] = stod(argv[i + 1]);
        } else if (arg_name.compare(parameter_names[0]) == 0 or arg_name.compare("-h") == 0) {
            cout << help_msg.str() << endl;
            return 1;
        } else {
            cout << "Invalid argument name: " << arg_name << endl;
            return -2;
        }
        arg_name.clear();

    }

    // Check if the required parameters were set or not
    if(corpus_file.empty() or output_file.empty() or method_name.empty()) {
        cout << "Please enter the required parameters: ";
        cout << help_msg_required.str() << endl;

        return -4;
    }

    /* Check if the constraints are satisfied */
    if( parameter_names[3].compare(constraintValues[0]) == 0 or
        parameter_names[3].compare(constraintValues[1]) == 0 or
        parameter_names[3].compare(constraintValues[2]) == 0 ) {
        cout << "Wrong method name: \'" << method_name << "\'. The possible values are " << "("
             << constraintValues[0] << ", " << constraintValues[1] << ", and " << constraintValues[2] << ")" << endl;
        return -5;
    }

    return 0;

}
