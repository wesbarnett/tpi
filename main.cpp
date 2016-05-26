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

const double N_A = 6.0221412927e23; // Avogadro
const double k_b = 1.380648813e-23; // Boltzmann

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

double lj(coordinates a, coordinates b, Atomtype at, triclinicbox box, double rcut2);

int main(int argc, char* argv[])
{

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
    const double beta = 1.0/(k_b * T);

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
    double exp_pe = 0.0;

    #pragma omp parallel for schedule(guided, 15) reduction(exp_pe:+)
    for (int frame_i = 0; frame_i < frame_n; frame_i++)
    {

        triclinicbox box = trj.GetBox(frame_i);
        double box_x = box.at(X).at(X);
        double box_y = box.at(Y).at(Y);
        double box_z = box.at(Z).at(Z);

        uniform_real_distribution<double> distrib_x(0.0, box_x);
        uniform_real_distribution<double> distrib_y(0.0, box_y);
        uniform_real_distribution<double> distrib_z(0.0, box_z);

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

            }

            exp_pe += exp(-pe * beta);

        }

    }


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
