//
// Created by abdulkadir on 23/04/19.
//

#ifndef FAST_BERN_NODE_H
#define FAST_BERN_NODE_H

#include <string>

using namespace std;

class Node {

public:
    string node;
    int count;
private:

    int path, code;

public:
    Node(string &n);

};


#endif //FAST_BERN_NODE_H
