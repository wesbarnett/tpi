
#ifndef ATOMTYPE_H
#define ATOMTYPE_H

#include <string>
#include <vector>
#include "gmxcpp/Trajectory.h"
#include <immintrin.h>
using namespace std;

class Atomtype {

    private:

        float c6;
        float c12;
        float rcut2;
        double tail_factor;
        int n;
        string name;
/* TODO
        __m256 rcut2_8;
        __m256 c12_8;
        __m256 c6_8;
*/

    public:

        Atomtype();
        Atomtype(const Trajectory &trj, string name, float sig1, float eps1, float sig2, float eps2, float rc2, float epsfact);
        double CalcPE(int frame_i, const Trajectory &trj, const coordinates &rand_xyz, const cubicbox_m256 &box, double vol) const;
        double GetC6() const;
        double GetC12() const;

};

#endif
