#ifndef PATHWAY_H
#define PATHWAY_H
#include "standardincludes.h"
#include "reaction.h"

using namespace std;

class Pathway
{
    public:
    Pathway(string id, string name)
    {
        this->id = id;
        this->name = name;
    }
    string id;
    string name;
    vector<Reaction *> reactions;
};

#endif // PATHWAY_H
