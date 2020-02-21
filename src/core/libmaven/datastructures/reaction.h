#ifndef REACTION_H
#define REACTION_H
#include "Compound.h"

class Reaction
{
    public:
    Reaction(string db, string id, string name)
    {
        this->db = db;
        this->id = id;
        this->name = name;
        reversable = false;
    }
    void addReactant(Compound *r, int s)
    {
        reactants.push_back(r);
        stoichiometry[r] = s;
    }
    void addProduct(Compound *p, int s)
    {
        products.push_back(p);
        stoichiometry[p] = s;
    }
    void setReversable(bool r) { reversable = r; }

    string db;
    string id;
    string name;
    deque<Compound *> reactants;
    deque<Compound *> products;
    map<Compound *, int> stoichiometry;
    bool reversable;
};

#endif // REACTION_H
