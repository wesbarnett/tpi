
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <string>
using namespace std;

class Atomtype {
    private:
        string name;
        // These are from the interactions with the test particle
        double c6;
        double c12;
    public:
        Atomtype();
        Atomtype(string name, double c6, double c12);
        string GetName();
        double GetC6();
        double GetC12();
};

#endif
