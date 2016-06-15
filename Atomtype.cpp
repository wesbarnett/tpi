
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
double Atomtype::CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, triclinicbox &box, double vol)
{
    double pe = 0.0;

    for (int atom_i = 0; atom_i < this->n; atom_i++)
    {
        coordinates atom_xyz = trj.GetXYZ(frame_i, this->name, atom_i);
        double r2 = distance2(atom_xyz, rand_xyz, box);
        if (r2 < this->rcut2)
        {
            double ri6 = 1.0/(pow(r2,3));
            pe += ri6*(this->c12*ri6 - this->c6);
        }
    }

    pe += this->n/vol * this->tail_factor;;

    return pe;
}
