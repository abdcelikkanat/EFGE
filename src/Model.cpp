#include "Model.h"

Model::Model(string f_path, string m_name, double s_alpha, double m_alpha, double d_rate, int dim, int neg, int w_size,  int num_iters, vector <double> &optionalParams) {

    corpus_path = f_path;
    method_name = m_name;

    window_size = w_size;
    dim_size = dim;
    negative_sample_size = neg;

    starting_alpha = s_alpha;
    decay_rate = d_rate;
    min_alpha = m_alpha;
    num_of_iters = num_iters;

    std_dev = optionalParams[0];

    Vocabulary vocab(corpus_path);
    node2Id = vocab.getNode2Id();
    total_nodes = vocab.getTotalNodes();
    vocab_size = vocab.getVocabSize();
    vocab_items = vocab.getVocabItems();

    // Set up sampling class
    vector <int> counts = vocab.getNodeCounts();
    uni = Unigram(vocab_size, counts, 0.75);

    emb0 = new double*[vocab_size];
    emb1 = new double*[vocab_size];
    for(int i=0; i<vocab_size; i++) {
        emb0[i] = new double[dim_size];
        emb1[i] = new double[dim_size];
    }

}

Model::~Model() {

    for(int i=0; i<vocab_size; i++) {
        delete [] emb0[i];
        delete [] emb1[i];
    }
    delete emb0;
    delete emb1;

}

double Model::sigmoid(double z) {

    if(z > 10)
        return 1.0;
    else if(z < -10)
        return 0.0;
    else
        return 1.0 / ( 1.0 +  exp(-z));

}


void Model::bernoulli_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < (int)contextIds.size(); i++) {
        z = 0.0;
        for (int d = 0; d < dim_size; d++)
            z += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = sigmoid(z);

        g = alpha * (labels[i] - z);

        for (int d = 0; d < dim_size; d++) {
            neule[d] += g * emb1[contextIds[i]][d];
        }

        for (int d = 0; d < dim_size; d++)
            emb1[contextIds[i]][d] += g * emb0[centerId][d];
    }
    for (int d = 0; d < dim_size; d++)
        emb0[centerId][d] += neule[d];


    delete[] neule;
}



void Model::poisson_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double eta, g, z;
    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < (int)contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = -exp( eta ) + labels[i];

        g = alpha * z;

        for (int d = 0; d < dim_size; d++) {
            neule[d] += g * emb1[contextIds[i]][d];
        }

        for (int d = 0; d < dim_size; d++)
            emb1[contextIds[i]][d] += g * emb0[centerId][d];
    }
    for (int d = 0; d < dim_size; d++)
        emb0[centerId][d] += neule[d];


    delete[] neule;
}





void Model::gaussian_known_var_exp(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;
    //double std_dev = 1.0;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < (int)contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = -(labels[i]/std_dev)*exp(-eta) + exp(-2.0*eta);
        //z = ( labels[i] - pow(eta, 2.0) / 2.0 ) * eta;

        g = alpha * z;

        for (int d = 0; d < dim_size; d++) {
            neule[d] += g * emb1[contextIds[i]][d];
        }

        for (int d = 0; d < dim_size; d++)
            emb1[contextIds[i]][d] += g * emb0[centerId][d];
    }
    for (int d = 0; d < dim_size; d++)
        emb0[centerId][d] += neule[d];


    delete[] neule;
}




void Model::run() {

    //default_random_engine generator;
    normal_distribution<double> normal_distr(0.0, 1.0);

    /* */
    // Initialize parameters
    uniform_real_distribution<double> real_distr(-0.5 /dim_size , 0.5/dim_size);

    for(int node=0; node<vocab_size; node++) {
        for(int d=0; d<dim_size; d++) {
            emb0[node][d] = real_distr(generator);
            emb1[node][d] = 0.0;
        }
    }


    fstream fs(corpus_path, fstream::in);
    if(fs.is_open()) {

        string line, token, center_node, context_node;
        vector <string> nodesInLine;
        int context_start_pos, context_end_pos;
        vector <double> x;
        vector <int> contextIds;
        int centerId;
        int *neg_sample_ids;
        double alpha;
        int processed_node_count = 0;

        alpha = starting_alpha;

        cout << "--> The update of the model parameters has started." << endl;

        for(int iter=0; iter<num_of_iters; iter++) {

            fs.clear();
            fs.seekg(0, ios::beg);
            cout << "    + Iteration: " << iter+1 << endl;

            while (getline(fs, line)) {
                stringstream ss(line);
                while (getline(ss, token, ' '))
                    nodesInLine.push_back(token);

                for (int center_pos = 0; center_pos < (int)nodesInLine.size(); center_pos++) {

                    // Update alpha
                    if (processed_node_count % 10000 == 0) {
                        alpha = starting_alpha * (1.0 - decay_rate * ((float) processed_node_count / (total_nodes*num_of_iters)));

                        if (alpha < min_alpha)
                            alpha = min_alpha;

                        cout << "\r    + Current alpha: " << setprecision(4) << alpha;
                        cout << " and " << processed_node_count-(total_nodes*iter) << "" << setprecision(3) << "("
                             << 100 * (float) ( processed_node_count-(total_nodes*iter) ) / total_nodes << "%) "
                             << "nodes in the file have been processed.";
                        cout << flush;
                    }

                    context_start_pos = max(0, center_pos - window_size);
                    context_end_pos = min(center_pos + window_size, (int)nodesInLine.size() - 1);

                    center_node = nodesInLine[center_pos];
                    centerId = node2Id[center_node];

                    // Resize
                    contextIds.resize((int) negative_sample_size + 1);
                    x.resize((int) negative_sample_size + 1);
                    neg_sample_ids = new int[negative_sample_size];
                    for (int context_pos = context_start_pos; context_pos <= context_end_pos; context_pos++) {

                        if (center_pos != context_pos) {
                            context_node = nodesInLine[context_pos];

                            contextIds[0] = node2Id[context_node];
                            uni.sample(negative_sample_size, neg_sample_ids);
                            for (int i = 0; i < negative_sample_size; i++)
                                contextIds[i + 1] = neg_sample_ids[i];
                            x[0] = 1.0;
                            fill(x.begin() + 1, x.end(), 0.0);

                            if(method_name.compare("bern") == 0) {
                                bernoulli_update(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("pois") == 0) {

                                x[0] = 1.0;
                                poisson_update(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("norm") == 0) {

                                x[0] = 1.0;
                                for(int n=1; n<=negative_sample_size; n++)
                                    x[n] = 0.0;
                                gaussian_known_var_exp(alpha, x, centerId, contextIds);

                            }else {
                                cout << "Not a valid method name" << endl;
                            }

                        }

                    }

                    // Increase the node count
                    processed_node_count++;

                    // Clear the vectors
                    contextIds.clear();
                    x.clear();
                    delete[] neg_sample_ids;
                }


                nodesInLine.clear();

            }
            cout << endl;

        }
        fs.close();

        cout << endl << "Done" << endl;

    } else {
        cout << "An error occurred during reading file!" << endl;
    }


}


void Model::save_embeddings(string file_path) {

    fstream fs(file_path, fstream::out);
    if(fs.is_open()) {
        fs << vocab_size << " " << dim_size << endl;
        for(int node=0; node<vocab_size; node++) {
            fs << vocab_items[node].node << " ";
            for(int d=0; d<dim_size; d++) {
                fs << emb0[node][d] << " ";
            }
            fs << endl;
        }

        fs.close();

    } else {
        cout << "An error occurred during opening the file!" << endl;
    }

}
