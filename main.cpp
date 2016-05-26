/* 
 * Test particle insertion
 * James W. Barnett
 * May 26, 2016
 * 
 * Limitations:
 *  - Must use a cubic box (does not have to be equal on all sides)
 *  - Particle being inserted has no charge, so electrostatics excluded
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <iostream>
#include <random>

#include "gmxcpp/Trajectory.h"

using namespace std;

const double R = 8.3144598e-3; // kJ/(mol*K) gas constant

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

Atomtype::Atomtype()
{

}

Atomtype::Atomtype(string name, double c6, double c12)
{
    this->name = name;
    this->c6 = c6;
    this->c12 = c12;
}

string Atomtype::GetName()
{
    return this->name;
}

double Atomtype::GetC6()
{
    return this->c6;
}

double Atomtype::GetC12()
{
    return this->c12;
}

// Units end up being kJ/mol if using GROMACS epsilons and sigmas
double lj(coordinates a, coordinates b, Atomtype at, triclinicbox box, double rcut2);
double tail(Atomtype at, double rc2, double vol);

int main(int argc, char* argv[])
{

    const int block_n = 5;

    /* BEGIN CONFIGURATION FILE PARSING */

    if (argc != 2)
    {
        cout << "ERROR: Configuration file should be first command line argument." << endl;
        return -1;
    }

    cout << "Using '" << argv[1] << "' for configuration file..." << endl;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(argv[1], pt);
    char *endptr;

    const string xtcfile = pt.get<std::string>("xtcfile","prd.xtc");
    cout << "xtcfile = " << xtcfile << endl;

    const string ndxfile = pt.get<std::string>("ndxfile","index.ndx");
    cout << "ndxfile = " << ndxfile << endl;

    const int rand_n = strtol(pt.get<std::string>("rand_n","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rand_n' must be an integer." << endl;
        return -1;
    }
    cout << "rand_n = " << rand_n << endl;

    const double rcut = strtod(pt.get<std::string>("rcut","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rcut' must be a double." << endl;
        return -1;
    }
    cout << "rcut = " << rcut << endl;
    const double rcut2 = rcut*rcut;

    const double T = strtod(pt.get<std::string>("T","298.15").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'T' must be a double." << endl;
        return -1;
    }
    cout << "T = " << T << endl;
    const double beta = 1.0/(R * T); // kJ/mol

    const int atomtypes = strtol(pt.get<std::string>("atomtypes","1").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'atomtypes' must be an integer." << endl;
        return -1;
    }
    cout << "atomtypes = " << atomtypes << endl;
    vector <Atomtype> at(atomtypes);
    istringstream iss; 
    string atomtype_name;
    double atomtype_sigma;
    double atomtype_epsilon;
    double randtype_sigma;
    double randtype_epsilon;

    string randtype_str = pt.get<std::string>("randtype");
    iss.clear();
    iss.str(randtype_str);
    iss >> randtype_sigma;
    iss >> randtype_epsilon;

    for (int i = 0; i < atomtypes; i++)
    {
        string atomtype_str = pt.get<std::string>("atomtype"+to_string(i+1));
        iss.clear();
        iss.str(atomtype_str);
        iss >> atomtype_name;
        iss >> atomtype_sigma;
        iss >> atomtype_epsilon;
        cout << atomtype_name << " " << atomtype_sigma << " " << atomtype_epsilon << endl;

        double eps = sqrt(randtype_epsilon * atomtype_epsilon);
        double sig = 0.5 * (randtype_sigma + atomtype_sigma);
        double c6 = 4.0 * eps * pow(sig, 6);
        double c12 = 4.0 * eps * pow(sig, 12);
        Atomtype at_tmp(atomtype_name, c6, c12);
        at.at(i) = at_tmp;
    }

    /* END CONFIGURATION FILE PARSING */

    Trajectory trj(xtcfile, ndxfile);
    const int frame_n = trj.GetNFrames();

    random_device rd;
    mt19937 gen(rd());
    vector <double> V_exp_pe(frame_n);
    vector <double> V(frame_n);

    #pragma omp parallel for schedule(guided, 15)
    for (int frame_i = 0; frame_i < frame_n; frame_i++)
    {

        cout << frame_i << endl;

        triclinicbox box = trj.GetBox(frame_i);
        double box_x = box.at(X).at(X);
        double box_y = box.at(Y).at(Y);
        double box_z = box.at(Z).at(Z);

        uniform_real_distribution<double> distrib_x(0.0, box_x);
        uniform_real_distribution<double> distrib_y(0.0, box_y);
        uniform_real_distribution<double> distrib_z(0.0, box_z);
        double vol = volume(box);

        double V_exp_pe_tmp = 0.0;

        for (int rand_i = 0; rand_i < rand_n; rand_i++)
        {
            coordinates rand_xyz(distrib_x(gen), distrib_y(gen), distrib_z(gen));
            double pe = 0.0;

            for (int atomtype_i = 0; atomtype_i < atomtypes; atomtype_i++)
            {

                string grp = at.at(atomtype_i).GetName();
                int atom_n = trj.GetNAtoms(grp);
                for (int atom_i = 0; atom_i < atom_n; atom_i++)
                {
                    coordinates atom_xyz = trj.GetXYZ(frame_i, grp, atom_i);
                    pe += lj(rand_xyz, atom_xyz, at.at(atomtype_i), box, rcut2);
                }

                pe += tail(at.at(atomtype_i), rcut2, (double)atom_n/vol);

            }

            V_exp_pe_tmp += vol * exp(-pe * beta);

        }

        V_exp_pe.at(frame_i) = V_exp_pe_tmp/(double)rand_n;
        V.at(frame_i) = vol;
    }

    vector <double> chem_pot_block(block_n);
    double chem_pot = 0.0;

    double V_avg = 0.0;
    double V_exp_pe_avg = 0.0;
    double chem_pot_block_avg = 0.0;

    for (int block_i = 0; block_i < block_n; block_i++)
    {

        int block_end;
        int block_start = block_i*frame_n/block_n;
        double V_block = 0.0;
        double V_exp_pe_block = 0.0;

        if (block_i == block_n-1)
        {
            block_end = frame_n;
        }
        else
        {
            block_end = (block_i+1)*frame_n/block_n;
        }

        for (int frame_i = block_start; frame_i < block_end; frame_i++)
        {
            V_block += V.at(frame_i);
            V_exp_pe_block += V_exp_pe.at(frame_i);
            V_avg += V.at(frame_i);
            V_exp_pe_avg += V_exp_pe.at(frame_i);
        }

        V_block /= (block_end-block_start);
        V_exp_pe_block /= (block_end-block_start);

        chem_pot_block.at(block_i) = -log(V_exp_pe_block/V_block) / beta;
        chem_pot_block_avg += -log(V_exp_pe_block/V_block) / beta;

    }
    chem_pot_block_avg /= block_n;

    double chem_pot_var = 0.0;
    for (int block_i = 0; block_i < block_n; block_i++)
    {
        chem_pot_var += pow(chem_pot_block_avg - chem_pot_block.at(block_i), 2);
    }
    chem_pot_var /= (block_n-1);

    V_avg /= frame_n;
    V_exp_pe_avg /= frame_n;
//TODO: divide by volume? Is that correct?
    chem_pot = -log(V_exp_pe_avg/V_avg) / beta;

    cout << chem_pot << " +/- " << sqrt(chem_pot_var) << " kJ / mol" << endl;

    return 0;

}

double lj(coordinates a, coordinates b, Atomtype at, triclinicbox box, double rcut2)
{
    double r2 = distance2(a, b, box);
    if (r2 < rcut2)
    {
        double ri6 = 1.0/(pow(r2,3));
        double ri12 = pow(ri6,2);
        return at.GetC12()*ri12 - at.GetC6()*ri6;
    }

    return 0.0;
}

double tail(Atomtype at, double rc2, double rho)
{
    double ri6 = 1.0/(pow(rc2,3));
    double ri12 = pow(ri6,2);
    return (2.0/3.0)*M_PI*rho*((at.GetC12()*ri12)/3.0 - at.GetC6()*ri6);

}
