//
// Created by abdulkadir on 13/11/18.
//

#include "Model_son.h"


template <typename A, typename B>
void zip(
        const std::vector<A> &a,
        const std::vector<B> &b,
        std::vector<std::pair<A,B>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

// Write the first and second element of the pairs in
// the given zipped vector into a and b. (This assumes
// that the vectors have equal length)
template <typename A, typename B>
void unzip(
        const std::vector<std::pair<A, B>> &zipped,
        std::vector<A> &a,
        std::vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}





Model::Model(Graph graph, int dim) {
    g = graph;
    dim_size = dim;
    num_of_nodes = graph.getNumOfNodes();
    adj_list = graph.getAdjList();
    getNeighbors();

    emb0 = new double*[num_of_nodes];
    emb1 = new double*[num_of_nodes];
    for(int i=0; i<num_of_nodes; i++) {
        emb0[i] = new double[dim_size];
        emb1[i] = new double[dim_size];
    }



}



Model::~Model() {

    for(int i=0; i<num_of_nodes; i++) {
        delete [] emb0[i];
        delete [] emb1[i];
    }
    delete emb0;
    delete emb1;


}

void Model::readGraph(string file_path, string filetype, bool directed) {

    g.readGraph(file_path, filetype, directed);
    num_of_nodes = g.getNumOfNodes();
    adj_list = g.getAdjList();

}


void Model::getNeighbors() {

    //int sample_size = 64;
    int sample_size = 128;

    nb_list.resize(num_of_nodes);

    int nb, nb_nb, nb_nb_nb, nb_nb_nb_nb;

    vector<int> nb_candidates, nb_nb_candidates, nb_nb_nb_candidates;
    vector<int> counts;

    // Get average node degree, max, min
    int max_degree = 0, min_degree = num_of_nodes;
    float avg_degree = 0;

    for (int i = 0; i < num_of_nodes; i++) {

        avg_degree += (float) adj_list[i].size();

        if (adj_list[i].size() > max_degree) {
            max_degree = adj_list[i].size();
        }

        if (adj_list[i].size() < min_degree) {
            min_degree = adj_list[i].size();
        }
    }

    cout << "Max degree: " << max_degree << endl;
    cout << "Min degree: " << min_degree << endl;
    cout << "Avg degree: " << avg_degree / (float) num_of_nodes << endl;

    vector<pair<int, int>> nb_temp_pair, nb_nb_temp_pair, nb_nb_nb_temp_pair;
    for (int node = 0; node < num_of_nodes; node++) {

        nb_candidates.clear();
        for (int nb_inx = 0; nb_inx < adj_list[node].size(); nb_inx++) {
            nb = adj_list[node][nb_inx]; // Get 1-neighbor
            nb_candidates.push_back(nb);
        }
        counts.clear();
        for (int nb_inx = 0; nb_inx < nb_candidates.size(); nb_inx++) {
            nb = nb_candidates[nb_inx];
            //counts.push_back((int)adj_list[nb].size());
            //counts.push_back((int) g.getClusteringCoefficient(node, nb));
            counts.push_back((int) g.getCommonNeighbours(node, nb).size());
        }
        nb_temp_pair.clear();

        zip(nb_candidates, counts, nb_temp_pair);
        sort(begin(nb_temp_pair), end(nb_temp_pair), [&](const auto &a, const auto &b) { return a.second < b.second; });
        unzip(nb_temp_pair, nb_candidates, counts);

        // Now add the nb list
        for (int i = 0; i < sample_size/1 && i < nb_candidates.size(); i++) {
            nb = nb_candidates[i];
            nb_list[node].push_back(nb); // Set nb
            nb_list[node].push_back(nb);
            nb_list[node].push_back(nb);

            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb); // Set nb
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);
            //nb_list[node].push_back(nb);

            nb_nb_candidates.clear();
            // For each nb_nb, apply the same procedure
            for (int nb_nb_inx = 0; nb_nb_inx < adj_list[nb].size(); nb_nb_inx++) {
                nb_nb = adj_list[nb][nb_nb_inx]; // Get 2-neighbor
                nb_nb_candidates.push_back(nb_nb);
            }
            counts.clear();
            for (int nb_nb_inx = 0; nb_nb_inx < nb_nb_candidates.size(); nb_nb_inx++) {
                nb_nb = nb_nb_candidates[nb_nb_inx];
                //counts.push_back((int)adj_list[nb_nb].size());
                //counts.push_back((int) g.getClusteringCoefficient(node, nb_nb));
                counts.push_back((int) g.getCommonNeighbours(node, nb_nb).size());
            }
            nb_nb_temp_pair.clear();

            zip(nb_nb_candidates, counts, nb_nb_temp_pair);
            sort(begin(nb_nb_temp_pair), end(nb_nb_temp_pair),
                 [&](const auto &a, const auto &b) { return a.second < b.second; });
            unzip(nb_nb_temp_pair, nb_nb_candidates, counts);


            // Now add the nb list
            for (int k = 0; k < sample_size/1 && k < nb_nb_candidates.size(); k++) {
                //nb_list[node].push_back(nb);
                nb_nb = nb_nb_candidates[k];
                nb_list[node].push_back(nb_nb); // Set nb
                //nb_list[node].push_back(nb_nb); // Set nb
                //nb_list[node].push_back(nb_nb); // Set nb
                //nb_list[node].push_back(nb_nb); // Set nb
                //nb_list[node].push_back(nb_nb); // Set nb
                //nb_list[node].push_back(nb_nb); // Set nb
                /*
                nb_nb_nb_candidates.clear();
                // For each nb_nb, apply the same procedure
                for(int nb_nb_nb_inx = 0; nb_nb_nb_inx < adj_list[nb_nb].size(); nb_nb_nb_inx++) {
                    nb_nb_nb = adj_list[nb_nb][nb_nb_nb_inx]; // Get 2-neighbor
                    nb_nb_nb_candidates.push_back(nb_nb_nb);
                }
                counts.clear();
                for(int nb_nb_nb_inx = 0; nb_nb_nb_inx < nb_nb_nb_candidates.size(); nb_nb_nb_inx++) {
                    nb_nb_nb = nb_nb_nb_candidates[nb_nb_nb_inx];
                    counts.push_back((int) g.getCommonNeighbours(node, nb_nb_nb).size());
                }
                nb_nb_nb_temp_pair.clear();
                zip(nb_nb_nb_candidates, counts, nb_nb_nb_temp_pair);
                sort(begin(nb_nb_nb_temp_pair), end(nb_nb_nb_temp_pair),
                     [&](const auto &a, const auto &b) { return a.second < b.second; });
                unzip(nb_nb_nb_temp_pair, nb_nb_nb_candidates, counts);

                for(int kk = 0; kk < sample_size/1 && kk < nb_nb_nb_candidates.size(); kk++) {
                    nb_nb_nb = nb_nb_nb_candidates[kk];
                    nb_list[node].push_back(nb_nb_nb); // Set nb_nb
                    //nb_list[node].push_back(nb_nb_nb); // Set nb_nb
                    //nb_list[node].push_back(nb_nb_nb); // Set nb_nb
                }
                */
            }

        }

    }
}

/*
void Model::getNeighbors() {

    int sample_size = 32;

    nb_list.resize(num_of_nodes);

    int nb, nb_nb, nb_nb_nb;

    vector <int> nb_candidates, nb_nb_candidates;
    vector <int> counts;

    // Get average node degree, max, min
    int max_degree=0, min_degree=num_of_nodes;
    float avg_degree=0;

    for(int i=0; i<num_of_nodes; i++) {

        avg_degree += (float)adj_list[i].size();

        if(adj_list[i].size() > max_degree) {
            max_degree = adj_list[i].size();
        }

        if(adj_list[i].size() < min_degree) {
            min_degree = adj_list[i].size();
        }
    }

    cout << "Max degree: " << max_degree << endl;
    cout << "Min degree: " << min_degree << endl;
    cout << "Avg degree: " << avg_degree/(float)num_of_nodes << endl;

    vector <pair <int, int>> nb_temp_pair, nb_nb_temp_pair;
    for(int node=0; node<num_of_nodes; node++) {

        nb_candidates.clear();
        for(int nb_inx=0; nb_inx<adj_list[node].size(); nb_inx++) {
            nb = adj_list[node][nb_inx]; // Get 1-neighbor
            nb_candidates.push_back(nb);
        }
        counts.clear();
        for(int nb_inx=0; nb_inx<nb_candidates.size(); nb_inx++) {
            nb = nb_candidates[nb_inx];
            //counts.push_back((int)adj_list[nb].size());
            counts.push_back((int)g.getClusteringCoefficient(node, nb));
        }
        nb_temp_pair.clear();

        zip(nb_candidates, counts, nb_temp_pair);
        sort(begin(nb_temp_pair), end(nb_temp_pair), [&](const auto& a, const auto& b) { return a.second < b.second; });
        unzip(nb_temp_pair, nb_candidates, counts);

        // Now add the nb list
        for(int i=0; i<sample_size && i<nb_candidates.size(); i++) {
            nb = nb_candidates[i];
            nb_list[node].push_back(nb); // Set nb


            nb_nb_candidates.clear();
            // For each nb_nb, apply the same procedure
            for(int nb_nb_inx=0; nb_nb_inx<adj_list[nb].size(); nb_nb_inx++) {
                nb_nb = adj_list[nb][nb_nb_inx]; // Get 2-neighbor
                nb_nb_candidates.push_back(nb_nb);
            }

            counts.clear();
            for(int nb_nb_inx=0; nb_nb_inx<nb_nb_candidates.size(); nb_nb_inx++) {
                nb_nb = nb_nb_candidates[nb_nb_inx];
                //counts.push_back((int)adj_list[nb_nb].size());
                counts.push_back((int)g.getClusteringCoefficient(node, nb_nb));
            }

            nb_nb_temp_pair.clear();

            zip(nb_nb_candidates, counts, nb_nb_temp_pair);
            sort(begin(nb_nb_temp_pair), end(nb_nb_temp_pair), [&](const auto& a, const auto& b) { return a.second < b.second; });
            unzip(nb_nb_temp_pair, nb_nb_candidates, counts);


            for(int j=0; j<sample_size && j<nb_nb_candidates.size(); j++) {
                nb_nb_nb = nb_nb_candidates[j];
                nb_list[node].push_back(nb_nb_nb); // Set nb_nb
            }


        }

//
//        unsigned long int size = adj_list[node].size();

//        for(int i=0; i<size; i++) {
//            for(int j=0; j<adj_list[nb_list[node][i]].size(); j++)
 //               nb_list[node].push_back(nb_list[node][i]);

 //           unsigned seed = chrono::system_clock::now().time_since_epoch().count();
 //           shuffle(nb_list[node].begin(), nb_list[node].end(), default_random_engine(seed) );

//        }
 //

 //
 //       int count = 0;
  //      vector <pair <int, int>> temp_pair;
   //     for(int i=0; i<candidates.size(); i++) {
//            count = (int)g.getCommonNeighbours(node, candidates[i]).size();
 //           counts.push_back(count);
 //       }
 //       zip(candidates, counts, temp_pair);
 //       sort(begin(temp_pair), end(temp_pair), [&](const auto& a, const auto& b) { return a.second > b.second; });
 //       unzip(temp_pair, candidates, counts);

 //       int adding_count=0;
 //       default_random_engine generator;
 //       discrete_distribution<int> node_distr(counts.begin(), counts.end());

 //       while(adding_count < sample_size) {
 //
 //
 //           if(adding_count < candidates.size()) {
 //               nb_list[node].push_back(candidates[adding_count]);
  //          } else {
 //               nb_list[node].push_back(candidates[node_distr(generator)]);
 //           }
  //
  //
  //          nb_list[node].push_back(candidates[node_distr(generator)]);

   //         adding_count++;

//        }

//        candidates.clear();
//        counts.clear();
//

//    }

//}

*/

/* YEDEK
void Model::getNeighbors() {

    nb_list.resize(num_of_nodes);

    int nb, nb_nb;

    for(int node=0; node<num_of_nodes; node++) {

        for(int nb_inx=0; nb_inx<adj_list[node].size(); nb_inx++) {

            nb = adj_list[node][nb_inx]; // Get 1-neighbor
            nb_list[node].push_back(nb); // Set nb

            for(int nb_nb_inx=0; nb_nb_inx<adj_list[nb].size(); nb_nb_inx++) {

                nb_nb = adj_list[nb][nb_nb_inx]; // Get 2-neighbor
                nb_list[node].push_back(nb_nb); // Set nb_nb
                // ############
                // Add one more time for each 2-neighbor
                nb_list[node].push_back(nb);
            }

        }

    }

}

*/




double Model::sigmoid(double z) {

    if(z > 6)
        return 1.0;
    else if(z < -6)
        return 0.0;
    else
        return 1.0 / ( 1.0 +  exp(-z));

}

void Model::run(double starting_alpha, int num_of_iters, int negative, int save_step, string save_file) {

    // Initialize parameters
    uniform_real_distribution<double> real_distr(-5.0 /dim_size , 5.0/dim_size);

    for(int node=0; node<num_of_nodes; node++) {
        for(int d=0; d<dim_size; d++) {
            emb0[node][d] = real_distr(generator);
            emb1[node][d] = 0.0;
        }
    }

    //save_embeddings("./citeseer0.embedding");

    // Set up sampling class
    vector <int> freq = g.getDegreeSequence();
    Unigram uni(num_of_nodes, freq, 0.75);

    stringstream embed_file_path;
    int *samples, *labels;
    int target, neg_sample_size, pos_sample_size, target_seq_size, label;
    double alpha, z, g, *neule;

    alpha = starting_alpha;

    for(int iter=0; iter<num_of_iters; iter++) {

        cout << "Iteration: " << iter << endl;

        for(int node=0; node<num_of_nodes; node++) {

            neg_sample_size = (int)nb_list[node].size()*negative;
            pos_sample_size = (int)nb_list[node].size();
            target_seq_size =  neg_sample_size + pos_sample_size;

            samples = new int[target_seq_size];
            labels = new int[target_seq_size];

            // set negative samples
            uni.sample(neg_sample_size, samples);
            //
            for(int l=0; l<neg_sample_size; l++)
                labels[l] = 0;
            //memset(labels, 0, (size_t)neg_sample_size);

            // set positive samples
            for(int j=neg_sample_size; j<target_seq_size; j++)
                samples[j] = nb_list[node][j-neg_sample_size];
            //
            for(int l=neg_sample_size; l<target_seq_size; l++)
                labels[l] = 1;
            //memset(labels+neg_sample_size, 1, (size_t)pos_sample_size);

            neule = new double[dim_size];
            //
            for(int d=0; d<dim_size; d++)
                neule[d] = 0.0;
            //memset(neule, 0, (size_t)dim_size);

            for(int j=0; j<target_seq_size; j++) {

                target = samples[j];
                label = labels[j];

                z = 0.0;
                for(int d=0; d<dim_size; d++)
                    z += emb0[node][d] * emb1[target][d];

                z = sigmoid(z);

                g = alpha * (label - z);

                for(int d=0; d<dim_size; d++) {
                    neule[d] += g * emb1[target][d];
                }

                for(int d=0; d<dim_size; d++)
                    emb1[target][d] += g*emb0[node][d];

            }

            for(int d=0; d<dim_size; d++)
                emb0[node][d] += neule[d];


            if


            delete [] samples;
            delete [] labels;
            delete [] neule;
        }


        if((iter+1) % save_step == 0) {
            embed_file_path.str("");
            embed_file_path << save_file + "_" << iter+1 << ".embedding";
            save_embeddings(embed_file_path.str());
        }

    }

}


void Model::save_embeddings(string file_path) {

    fstream fs(file_path, fstream::out);
    if(fs.is_open()) {
        fs << num_of_nodes << " " << dim_size << endl;
        for(int node=0; node<num_of_nodes; node++) {
            fs << node << " ";
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