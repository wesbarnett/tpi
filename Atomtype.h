
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <string>
#include <vector>
#include "gmxcpp/Trajectory.h"

class Atomtype {

    private:

        double c6;
        double c12;
        double rcut2;
        double tail_factor;
        int n;
        string name;
        __m256 rcut2_8;
        __m256 c12_8;
        __m256 c6_8;

    public:

        Atomtype();
        Atomtype(Trajectory &trj, string name, double sig1, double eps1, double sig2, double eps2, double rc2, double epsfact);
        double CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, cubicbox_m256 &box, double vol) const;
        double GetC6() const;
        double GetC12() const;

};

#endif
