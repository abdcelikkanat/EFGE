//
// Created by abdulkadir on 13/11/18.
//

#include "Graph.h"

Graph::Graph() {
    cout << "Constructor" << endl;
}

Graph::~Graph() {

}


vector <int> Graph::getCommonNeighbours(int u, int v) {

    vector <int> common_node_list;
    int u_nb_inx= 0, v_nb_inx=0;

    while(u_nb_inx<adjlist[u].size() && v_nb_inx<adjlist[v].size()) {
        if (adjlist[u][u_nb_inx] < adjlist[v][v_nb_inx]) {

            u_nb_inx++;

        } else {

            if(adjlist[u][u_nb_inx] == adjlist[v][v_nb_inx]) {
                common_node_list.push_back(adjlist[u][u_nb_inx]);
            }

            v_nb_inx++;
        }

    }

    return common_node_list;
}

double Graph::getClusteringCoefficient(int v, int u) {

    int u_nb_inx= 0, v_nb_inx=0;
    double common_count = 0;

    while(u_nb_inx<adjlist[u].size() && v_nb_inx<adjlist[v].size()) {
        if (adjlist[u][u_nb_inx] < adjlist[v][v_nb_inx]) {

            u_nb_inx++;

        } else {

            if(adjlist[u][u_nb_inx] == adjlist[v][v_nb_inx]) {
                common_count += 1.0;
            }

            v_nb_inx++;
        }

    }

    if(adjlist[v] > adjlist[u]) {

        return common_count / (double)adjlist[u].size();

    } else {

        return common_count / (double)adjlist[v].size();
    }

}


void Graph::vector2Adjlist(bool directed) {

    adjlist.resize(num_of_nodes);

    for(int j=0; j<num_of_edges; j++) {
        adjlist[edges[j][0]].push_back(edges[j][1]);
        if( !directed ) {
            adjlist[edges[j][1]].push_back(edges[j][0]);
        }
    }

    // Sort the nodes in the adjacency list
    for(int i=0; i<num_of_nodes; i++)
        sort(adjlist[i].begin(), adjlist[i].end());

}

void Graph::readGraph(string file_path, string filetype, bool directed) {

    if( filetype == "edgelist" ) {

        readEdgeList(file_path);
        // Convert edges into adjacent list
        vector2Adjlist(directed);

    } else {
        cout << "Unknown file type!" << endl;
    }



}

void Graph::readEdgeList(string file_path) {

    fstream fs(file_path, fstream::in);
    if(fs.is_open()) {
        int u, v;
        int maxNodeId=0, minNodeId = 0;


        while(fs >> u >> v) {
            vector <int> edge{u, v};
            edges.push_back(edge);
            num_of_edges++;

            if (u > maxNodeId) { maxNodeId = u; }
            if (v > maxNodeId) { maxNodeId = v; }
        }
        fs.close();

        num_of_nodes = (unsigned int ) maxNodeId - minNodeId + 1; // Fix the minNodeId

    } else {
        cout << "An error occurred during reading file!" << endl;
    }

}

unsigned int Graph::getNumOfNodes() {

    return num_of_nodes;
}

unsigned int Graph::getNumOfEdges() {

    return num_of_edges;
}

void Graph::getEdges() {

    for(int i=0; i<edges.size(); i++) {
        cout << edges[i][0] << " " << edges[i][1] << endl;
    }

}
void Graph::printAdjList() {

    for(int i=0; i<num_of_nodes; i++) {
        cout << i << ": ";
        for(int j=0; j<adjlist[i].size(); j++) {
            cout << j << " ";
        }
        cout << endl;
    }

}

vector <vector <int>> Graph::getAdjList() {

    return adjlist;

}


vector <int> Graph::getDegreeSequence() {

    vector <int> degree_seq(num_of_nodes);
    //degree_seq.resize(num_of_nodes);
    for(unsigned int i=0; i<num_of_nodes; i++) {
        degree_seq[i] = (unsigned int) adjlist[i].size();
    }



    return degree_seq;

}