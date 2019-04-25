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
    cout << "BURADAss" << endl;
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
    /* */
    nb_list.resize(num_of_nodes);
    cout << num_of_nodes;
    string file_path = "/home/abdulkadir/Desktop/citeseer_undirected.nblist";
    string line;
    string token;
    int node = 0;

    fstream fs(file_path, fstream::in);
    if(fs.is_open()) {

        while( getline(fs, line) ) {
            stringstream ss(line);
            while( getline(ss, token, ' ') ) {
                //cout << token << " ";
                nb_list[node].push_back(stoi(token));
            }
            //cout << endl;
            //cout << nb_list[node].size() << endl;
            node++;
        }
        fs.close();

    } else {
        cout << "An error occurred during reading file!" << endl;
    }

}







double Model::sigmoid(double z) {

    if(z > 10)
        return 1.0;
    else if(z < -10)
        return 0.0;
    else
        return 1.0 / ( 1.0 +  exp(-z));

}

void Model::run(double starting_alpha, double min_alpha, double decay_rate, int num_of_iters, int negative, int save_step, string save_file) {

    // Initialize parameters
    uniform_real_distribution<double> real_distr(-0.5 /dim_size , 0.5/dim_size);

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


            delete [] samples;
            delete [] labels;
            delete [] neule;

        }

        /* Update alpha */
        alpha = max(min_alpha, alpha*(1.0 - decay_rate*iter));
        /* ----------- */

        if((iter+1) % save_step == 0) {
            embed_file_path.str("");
            embed_file_path << save_file + "deneme_" << iter+1 << ".embedding";
            save_embeddings(embed_file_path.str());
            cout << "Alpha: " << alpha << endl;
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

