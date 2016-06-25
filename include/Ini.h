
#ifndef INI_H
#define INI_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

const double R = 8.3144598e-3; // kJ/(mol*K) gas constant

class Ini
{

public:

    string outfile;
    string xtcfile;
    string ndxfile;
    int frame_freq;
    int chunk_size;
    int rand_n;
    float rand_ni;
    int block_n;
    int boot_n;
    double epsfact;
    double rcut;
    double rcut2;
    double T;
    double beta;
    double betai;
    double testtype_sigma;
    double testtype_epsilon;
    int atomtypes;
    vector <string> atomtype_name;
    vector <double> atomtype_sigma;
    vector <double> atomtype_epsilon;

    Ini(char file[]);

};

#endif
