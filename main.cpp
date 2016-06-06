/* 
 * Test particle insertion
 * James W. Barnett
 * May 26, 2016
 * 
 * Limitations:
 *  - Must use a cubic box (does not have to be equal on all sides)
 *  - Particle being inserted has no charge, so electrostatics excluded
 *  - Every atom type needs to have its own index group!
 *  - Isothermal-isobaric ensemble
 */

#define BLOCK 15

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "omp.h"
#include <string>
#include <random>

#include "Atomtype.h"
#include "gmxcpp/Trajectory.h"

using namespace std;

const double R = 8.3144598e-3; // kJ/(mol*K) gas constant
const double oneThird = 1.0/3.0;
const double twoThirdPi = 2.0*oneThird*M_PI;
const double kcal = 4.184; // kJ

int main(int argc, char* argv[])
{
    chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now(); 

    /* BEGIN CONFIGURATION FILE PARSING */

    if (argc != 2)
    {
        cout << "ERROR: Configuration file should be first command line argument." << endl;
        return -1;
    }

    cout << "Using '" << argv[1] << "' for configuration file..." << endl;

    int nthreads;

    #pragma omp parallel
    #pragma omp master
    nthreads = omp_get_num_threads();

    cout << "Using " << nthreads << " OpenMP threads." << endl;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(argv[1], pt);
    char *endptr;

    const string outfile = pt.get<std::string>("outfile","tpi.dat");
    cout << "outfile = " << outfile << endl;

    ofstream ofs(outfile.c_str());
    ofs << scientific << setprecision(6) << left;

    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << "Test particle insertion program -- Wes Barnett" << endl;
    ofs << "-----------------------------------------------------------------------" << endl;
    time_t start_time = chrono::system_clock::to_time_t(start);
    ofs << setw(40) << "Started computation at:" << setw(20) << ctime(&start_time);
    ofs << setw(40) << "Number of OMP threads:" << setw(20) << nthreads << endl;

    ofs << setw(40) << "Configuration file:" << setw(20) << argv[1] << endl;
    ofs << setw(40) << "Output file:" << setw(20) << outfile << endl;

    const string xtcfile = pt.get<std::string>("xtcfile","prd.xtc");
    cout << "xtcfile = " << xtcfile << endl;
    ofs << setw(40) << "Compressed trajectory file:" << setw(20) << xtcfile << endl;

    const string ndxfile = pt.get<std::string>("ndxfile","index.ndx");
    cout << "ndxfile = " << ndxfile << endl;
    ofs << setw(40) << "Index file:" << setw(20) << ndxfile << endl;

    const int frame_freq = strtol(pt.get<std::string>("frame_freq","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'frame_freq' must be an integer." << endl;
        return -1;
    }
    cout << "frame_freq = " << frame_freq << endl;

    const int rand_n = strtol(pt.get<std::string>("rand_n","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rand_n' must be an integer." << endl;
        return -1;
    }
    cout << "rand_n = " << rand_n << endl;
    ofs << setw(40) << "Insertions of particle per frame:" << setw(20) << rand_n  << endl;

    const int block_n = strtol(pt.get<std::string>("block_n","5").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'block_n' must be an integer." << endl;
        return -1;
    }
    cout << "block_n = " << block_n << endl;
    ofs << setw(40) << "Blocks used in uncertainty analysis:" << setw(20) << block_n << endl;

    const double rcut = strtod(pt.get<std::string>("rcut","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rcut' must be a double." << endl;
        return -1;
    }
    cout << "rcut = " << rcut << endl;
    ofs << setw(40) << "Cutoff distance (nm):" << setw(20) << rcut << endl;
    const double rcut2 = rcut*rcut;

    const double T = strtod(pt.get<std::string>("T","298.15").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'T' must be a double." << endl;
        return -1;
    }
    cout << "T = " << T << endl;
    ofs << setw(40) << "System temperature (K):" << setw(20) << T << endl;
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
    double testtype_sigma;
    double testtype_epsilon;

    string testtype_str = pt.get<std::string>("testtype");
    iss.clear();
    iss.str(testtype_str);
    iss >> testtype_sigma;
    iss >> testtype_epsilon;

    ofs << setw(40) << "Test particle sigma (nm):" << setw(20) << testtype_sigma << endl;
    ofs << setw(40) << "Test particle epsilon (kJ/mol):" << setw(20) << testtype_epsilon << endl;
    ofs << setw(40) << "Number of atom types:" << setw(20) << atomtypes << endl;
    ofs << setw(20) << "INDEX NAME" << setw(20) << "SIGMA (nm)" << setw(20) << "EPSILON (kJ/mol)" << setw(20) << "C6" << setw(20) << "C12" << endl;

    for (int i = 0; i < atomtypes; i++)
    {
        string atomtype_str = pt.get<std::string>("atomtype"+to_string(i+1));
        iss.clear();
        iss.str(atomtype_str);
        iss >> atomtype_name;
        iss >> atomtype_sigma;
        iss >> atomtype_epsilon;
        cout << atomtype_name << " " << atomtype_sigma << " " << atomtype_epsilon << endl;

        double eps = sqrt(testtype_epsilon * atomtype_epsilon);
        double sig = 0.5 * (testtype_sigma + atomtype_sigma);
        double c6 = 4.0 * eps * pow(sig, 6);
        double c12 = 4.0 * eps * pow(sig, 12);
        ofs << setw(20) << atomtype_name << setw(20) << atomtype_sigma << setw(20) << atomtype_epsilon << setw(20) << c6 << setw(20) << c12 << endl;
        Atomtype at_tmp(atomtype_name, c6, c12, rcut2);
        at.at(i) = at_tmp;
    }

    /* END CONFIGURATION FILE PARSING */

    Trajectory trj(xtcfile, ndxfile);
    const int frame_n = trj.GetNFrames();
    ofs << setw(40) << "Number of frames used:" << setw(20) << frame_n << endl;

    /* BEGIN MAIN ANALYSIS */

    random_device rd;
    mt19937 gen(rd());
    vector <double> V_exp_pe(frame_n);
    vector <double> V(frame_n);

    #pragma omp parallel for schedule(guided, BLOCK)
    for (int frame_i = 0; frame_i < frame_n; frame_i++)
    {

        int thread_id = omp_get_thread_num();

        if (frame_i % frame_freq == 0)
        {
            cout << "Thread: " << thread_id << " Frame: " << frame_i << endl;
        }

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
                    pe += at.at(atomtype_i).CalcLJ(rand_xyz, atom_xyz, box);
                }

                pe += at.at(atomtype_i).CalcTail((double)atom_n/vol);

            }

            V_exp_pe_tmp += vol * exp(-pe * beta);

        }

        V_exp_pe.at(frame_i) = V_exp_pe_tmp/(double)rand_n;
        V.at(frame_i) = vol;
    }

    /* END MAIN ANALYSIS */

    /* BEGIN ERROR ANALYSIS */

    vector <double> chem_pot_block(block_n);
    double chem_pot = 0.0;

    double V_avg = 0.0;
    double V_exp_pe_avg = 0.0;
    double chem_pot_block_avg = 0.0;

    //TODO: bootstrap?
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
    chem_pot = -log(V_exp_pe_avg/V_avg) / beta;

    /* END ERROR ANALYSIS */

    cout << chem_pot << " ± " << sqrt(chem_pot_var) << " kJ / mol" << endl;
    cout << chem_pot/kcal << " ± " << sqrt(chem_pot_var)/kcal << " kcal / mol" << endl;
    cout << chem_pot*beta << " ± " << sqrt(chem_pot_var)*beta << " kT" << endl;

    end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);
    ofs << setw(40) << "Finished computation at:" << setw(20) << ctime(&end_time);
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << "FINAL RESULT - Excess chemical potential of adding test particle" << endl;
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << fixed;
    ofs << chem_pot << " ± " << sqrt(chem_pot_var) << " kJ / mol" << endl;
    ofs << chem_pot/kcal << " ± " << sqrt(chem_pot_var)/kcal << " kcal / mol" << endl;
    ofs << chem_pot*beta << " ± " << sqrt(chem_pot_var)*beta << " kT" << endl;
    ofs.close();

    return 0;

}
