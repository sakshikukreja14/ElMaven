#ifndef ISOTOPE_H
#define ISOTOPE_H
#include "standardincludes.h"

using namespace std;

class Isotope
{
    public:
    string name;
    double mass;
    double abundance;
    int N15;
    int C13;
    int S34;
    int H2;

    Isotope(string name, double mass, int c = 0, int n = 0, int s = 0, int h = 0)
    {
        this->mass = mass;
        this->name = name;
        C13 = c;
        N15 = n;
        S34 = s;
        H2 = h;
        abundance = 0; //naman: was unintialized
    }

    Isotope()
    {
        mass = 0;
        abundance = 0;
        N15 = 0;
        C13 = 0;
        S34 = 0;
        H2 = 0;
    }

    Isotope(const Isotope &b)
    {
        name = b.name; //naman: Consider performing initialization in initialization list.
        mass = b.mass;
        abundance = b.abundance;
        N15 = b.N15;
        S34 = b.S34;
        C13 = b.C13;
        H2 = b.H2;
    }

    Isotope &operator=(const Isotope &b)
    {
        name = b.name;
        mass = b.mass;
        abundance = b.abundance;
        N15 = b.N15;
        S34 = b.S34;
        C13 = b.C13;
        H2 = b.H2;
        return *this;
    }
};

#endif // ISOTOPE_H
