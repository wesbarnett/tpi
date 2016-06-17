/* 
 * Test particle insertion
 * Used SIMD intrinsics (at least AVX necessary)
 * James W. Barnett
 * May 26, 2016
 * 
 * Limitations:
 *  - Must use a cubic box (does not have to be equal on all sides)
 *  - Particle being inserted has no charge, so electrostatics excluded
 *  - Every atom type needs to have its own index group!
 *  - Isothermal-isobaric ensemble or NVT
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "omp.h"
#include <string>
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
        double rcut2;
        double tail_factor;
        int n;

        __m256 rcut2_8;
        __m256 c12_8;
        __m256 c6_8;

    public:

        Atomtype() { }
    
        Atomtype(Trajectory &trj, string name, double c6, double c12, double rc2)
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

        double CalcPE(int frame_i, Trajectory &trj, coordinates &rand_xyz, cubicbox_m256 &box, double vol)
        {
            double pe = 0.0;
            int atom_i = 0;

            /* BEGIN SIMD SECTION */

            coordinates8 rand_xyz8(rand_xyz), atom_xyz;
            __m256 r2_8, mask, r6, ri6, pe_tmp;
            __m256 pe_sum = _mm256_setzero_ps();
            float result[n] __attribute__((aligned (16)));

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
            for (int i = 0; i < 8; i++)
            {
                pe += result[i];
            }

            /* END SIMD SECTION */

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

            pe += this->n/vol * this->tail_factor;;

            return pe;
        }
};

int main(int argc, char* argv[])
{

    chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now(); 


    /***** BEGIN CONFIGURATION FILE PARSING *****/

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
    ofs << "         Test particle insertion program -- Wes Barnett" << endl;
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
    const float rand_ni = 1.0/rand_n;
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

    const int boot_n = strtol(pt.get<std::string>("boot_n","200").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'boot_n' must be an integer." << endl;
        return -1;
    }
    cout << "boot_n = " << boot_n << endl;
    ofs << setw(40) << "Bootstrapped iterations in uncertainty analysis:" << setw(20) << boot_n << endl;

    const double epsfact = strtod(pt.get<std::string>("epsfact","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'epsfact' must be a double." << endl;
        return -1;
    }
    cout << "epsfact = " << epsfact << endl;
    ofs << setw(40) << "Well depth factor:" << setw(20) << epsfact << endl;

    const double rcut = strtod(pt.get<std::string>("rcut","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rcut' must be a double." << endl;
        return -1;
    }
    cout << "rcut = " << rcut << endl;
    ofs << setw(40) << "Cutoff distance (nm):" << setw(20) << rcut << endl;
    const double rcut2 = rcut*rcut;

    const double T = strtod(pt.get<std::string>("T","298.15").c_str(), &endptr); // Kelvin
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'T' must be a double." << endl;
        return -1;
    }
    cout << "T = " << T << endl;
    ofs << setw(40) << "System temperature (K):" << setw(20) << T << endl;
    const double beta = 1.0/(R * T); // kJ/mol
    const double betai = (R * T);

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

    Trajectory trj(xtcfile, ndxfile);

    for (int i = 0; i < atomtypes; i++)
    {
        string atomtype_str = pt.get<std::string>("atomtype"+to_string(i+1));
        iss.clear();
        iss.str(atomtype_str);
        iss >> atomtype_name;
        iss >> atomtype_sigma;
        iss >> atomtype_epsilon;
        cout << atomtype_name << " " << atomtype_sigma << " " << atomtype_epsilon << endl;

        double eps = epsfact * sqrt(testtype_epsilon * atomtype_epsilon);
        double sig = 0.5 * (testtype_sigma + atomtype_sigma);
        double c6 = 4.0 * eps * pow(sig, 6);
        double c12 = 4.0 * eps * pow(sig, 12);
        ofs << setw(20) << atomtype_name << setw(20) << atomtype_sigma << setw(20) << atomtype_epsilon << setw(20) << c6 << setw(20) << c12 << endl;
        Atomtype at_tmp(trj, atomtype_name, c6, c12, rcut2);
        at[i] = at_tmp;
    }

    /***** END CONFIGURATION FILE PARSING *****/


    /***** BEGIN MAIN ANALYSIS *****/

    const int frame_n = trj.GetNFrames();
    ofs << setw(40) << "Number of frames used:" << setw(20) << frame_n << endl;

    random_device rd;
    mt19937 gen(rd());
    double V_avg = 0.0;
    double V_exp_pe_avg = 0.0;

    vector <double> V_exp_pe(frame_n,0.0);
    vector <double> V(frame_n);

    #pragma omp parallel
    {
        cubicbox_m256 box;
        int thread_id;
        double pe;
        vector <coordinates> rand_xyz;

        #pragma omp for schedule(static) reduction(+:V_avg,V_exp_pe_avg)
        for (int frame_i = 0; frame_i < frame_n; frame_i++)
        {

            box = trj.GetCubicBoxM256(frame_i);
            gen_rand_box_points(rand_xyz, box, rand_n);
            V[frame_i] = volume(box);

            for (int rand_i = 0; rand_i < rand_n; rand_i++)
            {
                pe = 0.0;
                for (int atomtype_i = 0; atomtype_i < atomtypes; atomtype_i++)
                {
                    pe += at[atomtype_i].CalcPE(frame_i, trj, rand_xyz[rand_i], box, V[frame_i]);
                }

                V_exp_pe[frame_i] += exp(-pe * beta);
            }

            V_exp_pe[frame_i] *= V[frame_i] * rand_ni;
            V_avg += V[frame_i];
            V_exp_pe_avg += V_exp_pe[frame_i];

            if (frame_i % frame_freq == 0)
            {
                thread_id = omp_get_thread_num();
                printf("Thread: %-1d ", thread_id); 
                printf("Frame: %-8d ", frame_i);
                printf("μ = %-12.6f ", -log(V_exp_pe[frame_i]/V[frame_i]) * betai);
                printf("<μ> = %-12.6f\n", -log(V_exp_pe_avg/V_avg) * betai);
            }

        }

    }

    /***** END MAIN ANALYSIS *****/


    /***** BEGIN ERROR ANALYSIS *****/

    vector <double> chem_pot_boot(boot_n);
    double chem_pot_boot_avg = 0.0;
	uniform_int_distribution<int> dist(0,block_n-1);

    #pragma omp parallel for reduction(+:chem_pot_boot_avg)
    for (int boot_i = 0; boot_i < boot_n; boot_i++)
    {

        double V_boot = 0.0;
        double V_exp_pe_boot = 0.0;

        for (int block_i = 0; block_i < block_n; block_i++)
        {

            int block = dist(gen);
            int block_end;
            int block_start = (double)block/block_n * frame_n;

            if (block == block_n-1)
            {
                block_end = frame_n;
            }
            else
            {
                block_end = (double)(block+1)/block_n * frame_n;
            }

            for (int frame_i = block_start; frame_i < block_end; frame_i++)
            {
                V_boot += V[frame_i];
                V_exp_pe_boot += V_exp_pe[frame_i];
            }

        }

        chem_pot_boot[boot_i] = -log(V_exp_pe_boot/V_boot) * betai;
        chem_pot_boot_avg += chem_pot_boot[boot_i];

    }

    chem_pot_boot_avg /= boot_n;

    double chem_pot_boot_var = 0.0;
    for (int boot_i = 0; boot_i < boot_n; boot_i++)
    {
        chem_pot_boot_var += pow(chem_pot_boot_avg - chem_pot_boot[boot_i], 2);
    }
    chem_pot_boot_var /= (boot_n-1);
    double chem_pot = -log(V_exp_pe_avg/V_avg) * betai;

    /***** END ERROR ANALYSIS *****/


    cout << "---------------------------------------------------------" << endl;
    cout << "    μ (kJ / mol) = " << chem_pot << " ± " << sqrt(chem_pot_boot_var) << endl;
    cout << "---------------------------------------------------------" << endl;

    end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);
    ofs << setw(40) << "Finished computation at:" << setw(20) << ctime(&end_time);
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << "     FINAL RESULT - Excess chemical potential of test particle" << endl;
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << fixed;
    ofs << "μ (kJ / mol) = " << chem_pot << " ± " << sqrt(chem_pot_boot_var) << endl;
    ofs.close();

    return 0;

}
