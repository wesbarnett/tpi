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
const double k_B = 1.380648813e-23; // Boltzmann

class Atomtype {
    private:
        string name;
        double sigma;
        double epsilon;
    public:
        Atomtype();
        Atomtype(string name, double sigma, double epsilon);
        string GetName();
        double GetSigma();
        double GetEpsilon();
};

Atomtype::Atomtype()
{

}

Atomtype::Atomtype(string name, double sigma, double epsilon)
{
    this->name = name;
    this->sigma = sigma;
    this->epsilon = epsilon;
}

string Atomtype::GetName()
{
    return this->name;
}

double Atomtype::GetSigma()
{
    return this->sigma;
}

double Atomtype::GetEpsilon()
{
    return this->epsilon;
}

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

    const int atomtypes = strtol(pt.get<std::string>("atomtypes","1").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'atomtypes' must be an integer." << endl;
        return -1;
    }
    cout << "atomtypes = " << atomtypes << endl;
    vector <Atomtype> at(atomtypes);

    for (int i = 0; i < atomtypes; i++)
    {
        string atomtype_str = pt.get<std::string>("atomtype"+to_string(i+1));
        string atomtype_name;
        double atomtype_sigma;
        double atomtype_epsilon;
        istringstream iss; 
        iss.str(atomtype_str);
        iss >> atomtype_name;
        iss >> atomtype_sigma;
        iss >> atomtype_epsilon;
        Atomtype at_tmp(atomtype_name, atomtype_sigma, atomtype_epsilon);
        at.at(i) = at_tmp;
    }

    /* END CONFIGURATION FILE PARSING */

    Trajectory trj(xtcfile, ndxfile);
    const int frame_n = trj.GetNFrames();

    random_device rd;
    mt19937 gen(rd());

    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 15)
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
                coordinates rand_point(distrib_x(gen), distrib_y(gen), distrib_z(gen));
            }

        }

        #pragma omp critical
        {

        }

    }

    return 0;

}
