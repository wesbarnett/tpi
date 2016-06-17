
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include "gmxcpp/Trajectory.h"
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
        int n;
        double CalcLJ(coordinates &a, coordinates &b, triclinicbox &box);
        __m256 rcut2_8;
        __m256 c12_8;
        __m256 c6_8;
    public:
        Atomtype();
        Atomtype(Trajectory &trj, string name, double c6, double c12, double rc2);
        double CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, cubicbox_m256 &box, double vol);
};

#endif
