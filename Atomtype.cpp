
#include "Atomtype.h"

Atomtype::Atomtype()
{

}

Atomtype::Atomtype(Trajectory &trj, string name, double c6, double c12, double rc2)
{
    this->name = name;
    this->c6 = c6;
    this->c12 = c12;
    this->rcut2 = rc2;
    this->n = trj.GetNAtoms(this->name);

    double ri6 = 1.0/(pow(rc2,3));
    this->tail_factor = 2.0/3.0 * M_PI * ri6*(1.0/3.0*this->c12*ri6 - this->c6);

}

// Get the total PE contribution for a test particle for all atom of this atom type
double Atomtype::CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, triclinicbox &box)
{
    double pe = 0.0;

    for (int atom_i = 0; atom_i < this->n; atom_i++)
    {
        coordinates atom_xyz = trj.GetXYZ(frame_i, this->name, atom_i);
        pe += CalcLJ(rand_xyz, atom_xyz, box);
    }

    return pe;
}

// Lennard-Jones interaction between an atom of this atomtype and the test particle
double Atomtype::CalcLJ(coordinates &a, coordinates &b, triclinicbox &box)
{
    double r2 = distance2(a, b, box);
    if (r2 < this->rcut2)
    {
        double ri6 = pow(1.0/r2,3);
        return ri6*(this->c12*ri6 - this->c6);
    }

    return 0.0;
}

// LJ tail correction for this atom type (depends on density)
double Atomtype::CalcTail(double vol)
{
    return this->n/vol * tail_factor;
}
