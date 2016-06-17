
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

    this->rcut2_8 = _mm256_set1_ps(rcut2);
    this->c12_8 = _mm256_set1_ps(c12);
    this->c6_8 = _mm256_set1_ps(c6);
}

// Get the total PE contribution for a test particle for all atom of this atom type
double Atomtype::CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, cubicbox_m256 &box, double vol)
{
    coordinates8 rand_xyz8(rand_xyz), atom_xyz;
    __m256 r2_8, mask, r6, ri6, pe_tmp;
    __m256 pe_sum = _mm256_setzero_ps();
    float result[n] __attribute__((aligned (16)));
    double pe = 0.0;
    int atom_i = 0;
    for (; atom_i < this->n-8; atom_i+=8)
    {
        atom_xyz = trj.GetXYZ8(frame_i, this->name, atom_i);
        r2_8 = distance2(atom_xyz, rand_xyz8, box);
        mask = _mm256_cmp_ps(r2_8, rcut2_8, _CMP_LT_OS);
        r6 = _mm256_and_ps(mask, _mm256_mul_ps(_mm256_mul_ps(r2_8, r2_8), r2_8));
        ri6 = _mm256_and_ps(mask, _mm256_rcp_ps(r6));
        pe_tmp = _mm256_and_ps(mask, _mm256_mul_ps(ri6, _mm256_sub_ps(_mm256_mul_ps(c12_8, ri6), c6_8)));
        pe_sum = _mm256_add_ps(pe_tmp, pe_sum);
    }
    _mm256_store_ps(result, pe_sum);
    for (; atom_i < this->n; atom_i++)
    {
        coordinates atom_xyz = trj.GetXYZ(frame_i, this->name, atom_i);
        double r2 = distance2(atom_xyz, rand_xyz, cubicbox(box));
        if (r2 < this->rcut2)
        {
            double ri6 = 1.0/(pow(r2,3));
            pe += ri6*(this->c12*ri6 - this->c6);
        }
    }
    for (int i = 0; i < 8; i++)
    {
        pe += result[i];
    }

    pe += this->n/vol * this->tail_factor;;

    return pe;
}
