//
// Created by abdulkadir on 23/04/19.
//

#include "Model.h"

Model::Model(string f_path, int w_size, int neg, double s_alpha, double d_rate, double m_alpha, int num_iters, int dim, string m_name) {

    method_name = m_name;

    file_path = f_path;
    window_size = w_size;
    dim_size = dim;
    negative_sample_size = neg;

    Vocabulary vocab(file_path);
    node2Id = vocab.getNode2Id();
    total_nodes = vocab.getTotalNodes();
    vocab_size = (int)vocab.getVocabSize();
    vocab_items = vocab.getVocabItems();

    starting_alpha = s_alpha;
    decay_rate = d_rate;
    min_alpha = m_alpha;
    num_of_iters = num_iters;

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


vector <unordered_map <int, int>> Model::getCoOccurenceCount() {

    string line, token, center_node, context_node;
    vector <string> nodesInLine;
    int context_start_pos, context_end_pos;
    vector <int> contextIds;
    int centerId, contextId;
    vector <unordered_map <int, int>> cooccurences;


    fstream fs(file_path, fstream::in);
    if(fs.is_open()) {

        cooccurences.resize(vocab_size);

        while (getline(fs, line)) {
            stringstream ss(line);
            while (getline(ss, token, ' '))
                nodesInLine.push_back(token);

            for (int center_pos = 0; center_pos < nodesInLine.size(); center_pos++) {


                context_start_pos = max(0, center_pos - window_size);
                context_end_pos = min(center_pos + window_size, (int) nodesInLine.size() - 1);

                center_node = nodesInLine[center_pos];
                centerId = node2Id[center_node];

                for (int context_pos = context_start_pos; context_pos <= context_end_pos; context_pos++) {

                    if (center_pos != context_pos) {

                        context_node = nodesInLine[context_pos];
                        contextId = node2Id[context_node];

                        auto search = cooccurences[centerId].find(contextId);
                        if(search == cooccurences[centerId].end()) { // if not exist
                            cooccurences[centerId][contextId] = 1;

                        } else {
                            cooccurences[centerId][contextId]++;
                        }

                    }

                }

                // Clear the vectors
                contextIds.clear();
            }

            nodesInLine.clear();

        }
        fs.close();


    } else {
        cout << "An error occurred during reading file!" << endl;
    }

    return cooccurences;

}

vector <unordered_map <int, double>> Model::getRelativeFreq(vector <unordered_map <int, int>> cooccurences) {

    double row_sum;
    vector <unordered_map <int, double>> relative_freq;
    relative_freq.resize(vocab_size);

    for(int i=0; i<vocab_size; i++) {
        row_sum = 0.0;
        for(auto it=cooccurences[i].begin(); it != cooccurences[i].end(); ++it)
            row_sum += (float) cooccurences[i][it->first];

        for(auto it=cooccurences[i].begin(); it != cooccurences[i].end(); ++it)
            relative_freq[i][it->first] = (float) cooccurences[i][it->first] / row_sum;

    }

    return relative_freq;

}


vector <unordered_map <int, double>> Model::getNormalizedFreq(vector <unordered_map <int, int>> cooccurences) {

    double row_sum, mean, bias_std_dev;
    vector <unordered_map <int, double>> relative_freq;
    relative_freq.resize(vocab_size);

    for(int i=0; i<vocab_size; i++) {
        row_sum = 0.0;
        for(auto it=cooccurences[i].begin(); it != cooccurences[i].end(); ++it)
            row_sum += (float) cooccurences[i][it->first];
        mean = row_sum / cooccurences[i].size();

        bias_std_dev = 0.0;
        for(auto it=cooccurences[i].begin(); it != cooccurences[i].end(); ++it)
            bias_std_dev += pow( (float) cooccurences[i][it->first] - mean, 2 );
        bias_std_dev = sqrt( bias_std_dev / cooccurences[i].size() );

        if(cooccurences[i].size() == 1)
            bias_std_dev = 1.0;
        //cout << cooccurences[i].size() << " " << bias_std_dev << endl;

        for(auto it=cooccurences[i].begin(); it != cooccurences[i].end(); ++it)
            relative_freq[i][it->first] = (float) (cooccurences[i][it->first] - mean) / bias_std_dev;

    }

    return relative_freq;

}


void Model::bernoulli_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
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


void Model::bernoulli_update_v2(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double eta, g, delta;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        delta = exp( labels[i]*eta ) * sigmoid(-eta) * ( labels[i] - sigmoid(eta) );

        g = alpha * delta;

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


void Model::poisson_update_v1(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double eta, g, z;
    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
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

void Model::poisson_update_v2(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double eta, g, z;
    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        if(labels[i] > 0.0)
            z = -exp( eta );
        else
            z = - log(sigmoid( exp( eta ) ));

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

void Model::poisson_update_v3(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    // link function, log
    double *neule;
    double eta, g, z;
    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        if(labels[i] > 0.0)
            z = -eta;
        else
            z = -log(sigmoid( eta ));

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


void Model::poisson_update_v4(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {




    // link function, log
    double *neule;
    double eta, g, z;
    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++) {
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];
        }

        eta = exp(eta);
        if(labels[i] > 0.0)
            z = -(1.0 / (1.0 - exp(eta)) )*eta;
            //z = exp(-eta) / (exp(-eta) - 1.0);
        else
            z = -eta;


        g = alpha * z;

        for (int d = 0; d < dim_size; d++) {
            neule[d] += g * emb1[contextIds[i]][d];
        }

        for (int d = 0; d < dim_size; d++) {
            emb1[contextIds[i]][d] += g * emb0[centerId][d];
            //emb1[contextIds[i]][d] = max(0.0, emb1[contextIds[i]][d]);
        }
    }
    for (int d = 0; d < dim_size; d++) {
        emb0[centerId][d] += neule[d];
        //emb0[centerId][d] = max(0.0, emb0[centerId][d]);
    }




    delete[] neule;
}


void Model::exponential_update_v1(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = labels[i] + 1 / eta;

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

void Model::gaussian_known_var(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;
    double std_dev = 1.0;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = labels[i]/std_dev - eta;
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

void Model::gaussian_known_var_exp(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;
    double std_dev = 1.0;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = -labels[i]*exp(-eta) + exp(-2.0*eta);
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


void Model::pareto_update(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;
    double x_min = 1;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = log(labels[i]) + 1.0 / (1.0 + eta) + log(x_min);

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

void Model::pareto_update_v2(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;
    double x_min = 1;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = exp(eta)*log(labels[i]) + sigmoid(eta);

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

void Model::gaussian_my_prior(double alpha, vector <double> labels, int centerId, vector <int> contextIds) {
    double *neule;
    double z, g, eta;

    neule = new double[dim_size];
    for (int d = 0; d < dim_size; d++)
        neule[d] = 0.0;

    for (int i = 0; i < contextIds.size(); i++) {
        eta = 0.0;
        for (int d = 0; d < dim_size; d++)
            eta += emb0[centerId][d] * emb1[contextIds[i]][d];

        z = labels[i] - 2.0*eta;

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


    /*
    // Initialize parameters
    uniform_real_distribution<double> pois_uni(-2000.0 , 2000.0);

    for(int node=0; node<vocab_size; node++) {
        for(int d=0; d<dim_size; d++) {
            //emb0[node][d] = max(0.0, pois_uni(generator));
            //emb1[node][d] =  max(0.0, pois_uni(generator));
            emb0[node][d] = pois_uni(generator);
            emb1[node][d] =  pois_uni(generator);
        }
    }
    */



    vector <unordered_map <int, int>> cooccurences = getCoOccurenceCount();
    vector <unordered_map <int, double>> relative_freq = getRelativeFreq(cooccurences);
    vector <unordered_map <int, double>> normalized_freq = getNormalizedFreq(cooccurences);

    fstream fs(file_path, fstream::in);
    if(fs.is_open()) {

        string line, token, center_node, context_node;
        vector <string> nodesInLine;
        int context_start_pos, context_end_pos;
        vector <double> x;
        vector <int> contextIds;
        int centerId;
        double z, g, *neule;
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

                for (int center_pos = 0; center_pos < nodesInLine.size(); center_pos++) {

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
                    context_end_pos = min(center_pos + window_size, (int) nodesInLine.size() - 1);

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

                            //////////////////////////
                            //x[0] = cooccurences[centerId][contextIds[0]] / vocab_items[centerId].count;
                            //x[0] = 15.0*relative_freq[centerId][contextIds[0]];

                            /*
                            x[0]= 0;
                            for(int count_pos=context_start_pos; count_pos<=context_end_pos; count_pos++) {
                                if(center_pos != count_pos && node2Id[context_node] == node2Id[nodesInLine[count_pos]]) {
                                    x[0] += 1;
                                }
                            }
                            //cout << x[0] << endl;
                            poisson_update_v1(alpha, x, centerId, contextIds);
                            */

                            //bernoulli_update(alpha, x, centerId, contextIds);
                            //////////////////////////
                            //gaussian_known_var(alpha, x, centerId, contextIds);
                            //gaussian_my_prior(alpha, x, centerId, contextIds);

                            //x[0] = pow(1.0*relative_freq[centerId][contextIds[0]], 0.15);

                            /*
                            x[0] =  relative_freq[centerId][contextIds[0]]+ 0.0;
                            x[1] = normal_distr(generator)-2.0; //-2 best
                            x[2] = normal_distr(generator)-2.0;
                            x[3] = normal_distr(generator)-2.0;
                            x[4] = normal_distr(generator)-2.0;
                            x[5] = normal_distr(generator)-2.0;

                            gaussian_known_var(alpha, x, centerId, contextIds);
                            */

                            //x[0] = 2.0*cooccurences[centerId][contextIds[0]] / vocab_items[centerId].count;
                            //x[1] = normal_distr(generator)-0.0;
                            //x[2] = normal_distr(generator)-0.0;
                            //x[3] = normal_distr(generator)-0.0;
                            //x[4] = normal_distr(generator)-0.0;
                            //x[5] = normal_distr(generator)-0.0;

                            //poisson_update_v4(alpha, x, centerId, contextIds);
                            /////////////////////////////
                            //exponential_update_v1(alpha, x, centerId, contextIds);
                            ////////////////////////////////////////
                            //poisson_update_v2(alpha, x, centerId, contextIds);
                            //poisson_update_v3(alpha, x, centerId, contextIds);
                            ///////////////////////////////////////////////

                            if(method_name.compare("method1") == 0) {
                                //cout << "method1" << endl;
                                bernoulli_update(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("method2") == 0) {
                                //cout << "method2" << endl;
                                /* */
                                x[0] = 1;
                                poisson_update_v1(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("method3") == 0) {
                                //cout << "method3" << endl;

                            } else if(method_name.compare("method4_eski") == 0) {
                                //cout << "method4" << endl;

                                x[0] =  relative_freq[centerId][contextIds[0]]+ 0.0;
                                x[1] = normal_distr(generator)-5.0;
                                x[2] = normal_distr(generator)-5.0;
                                x[3] = normal_distr(generator)-5.0;
                                x[4] = normal_distr(generator)-5.0;
                                x[5] = normal_distr(generator)-5.0;

                                gaussian_known_var(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("method4") == 0) {
                                x[0] = 10.0;
                                for(int n=1; n<=negative_sample_size; n++)
                                    x[n] = 5.0;
                                gaussian_known_var(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("method4_exp") == 0) {
                                x[0] = 1.0;
                                for(int n=1; n<=negative_sample_size; n++)
                                    x[n] = 0.0;
                                gaussian_known_var_exp(alpha, x, centerId, contextIds);


                            } else if(method_name.compare("pareto") == 0) {


                               for(int n=0; n < negative_sample_size + 1; n++) {
                                   if(cooccurences[centerId][contextIds[0]] < 1)
                                       x[n] = 1;
                                   else
                                       x[n] = cooccurences[centerId][contextIds[n]]; // / (double)vocab_items[centerId].count;
                               }

                                pareto_update(alpha, x, centerId, contextIds);

                            } else if(method_name.compare("pareto2") == 0) {


                                for(int n=0; n < negative_sample_size + 1; n++) {
                                    if(cooccurences[centerId][contextIds[0]] < 1)
                                        x[n] = 1;
                                    else
                                        x[n] = cooccurences[centerId][contextIds[n]]; // / (double)vocab_items[centerId].count;
                                }
                                pareto_update_v2(alpha, x, centerId, contextIds);

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
