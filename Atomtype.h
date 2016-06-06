
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include "gmxcpp/coordinates.h"
#include "gmxcpp/triclinicbox.h"
#include "gmxcpp/Utils.h"
#include <string>
using namespace std;

class Atomtype {
    private:
        string name;
        // These are from the interactions with the test particle
        double c6;
        double c12;
        double rcut2;
        double tail_factor;
    public:
        Atomtype();
        Atomtype(string name, double c6, double c12, double rcut2);
        string GetName();
        double GetC6();
        double GetC12();
        double CalcLJ(coordinates a, coordinates b, triclinicbox box);
        double CalcTail(double rho);
};

#endif
