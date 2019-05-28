//
//
//

#include "Vocabulary.h"



Vocabulary::Vocabulary(string file_path) {

    fstream fs(file_path, fstream::in);
    if(fs.is_open()) {

        string line, token;

        total_nodes = 0;
        while( getline(fs, line) ) {
            stringstream ss(line);
            while( getline(ss, token, ' ') ) {

                // if the node was not added before
                auto token_search = node2Id.find(token);
                if (token_search == node2Id.end()) {

                    node2Id[token] = (int)vocab_items.size();
                    vocab_items.push_back(Node(token));
                }

                vocab_items[node2Id[token]].count++;
                total_nodes++;

            }

        }
        fs.close();

    } else {
        cout << "An error occurred during reading file." << endl;
    }

    cout << "--> The input file has been successfully read!" << endl;
    cout << "    + The total number of nodes: " << total_nodes << endl;
    cout << "    + The total number of distinct nodes: " << vocab_items.size() << endl;
}

unsigned long Vocabulary::getVocabSize() {

    return node2Id.size();

}

int Vocabulary::getTotalNodes() {

    return total_nodes;

}

vector <int> Vocabulary::getNodeCounts() {

    vector <int> nodeCounts;
    nodeCounts.resize(vocab_items.size());

    for(int i=0; i<vocab_items.size(); i++)
        nodeCounts[i] = vocab_items[i].count;

    return nodeCounts;

}

vector <Node> Vocabulary::getVocabItems() {

    return vocab_items;

}

unordered_map <string, int> Vocabulary::getNode2Id() {

    return node2Id;

}