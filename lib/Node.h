

#ifndef EFGE_NODE_H
#define EFGE_NODE_H

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


#endif // EFGE_NODE_H
