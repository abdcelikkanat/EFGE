//
//
//

#ifndef EFGE_VOCABULARY_H
#define EFGE_VOCABULARY_H

#include "Graph.h"
#include <string>
#include <sstream>
#include <unordered_map>
#include "Node.h"
#include <vector>

using namespace std;

class Vocabulary {
private:
    unordered_map <string, int> node2Id;
    vector <Node> vocab_items;
    int total_nodes;

public:
    Vocabulary(string file_path);
    int getTotalNodes();
    int getVocabSize();
    vector <int> getNodeCounts();
    vector <Node> getVocabItems();
    unordered_map <string, int> getNode2Id();




};


#endif //EFGE_VOCABULARY_H
