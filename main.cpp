/* 
 * Test particle insertion
 * Used SIMD intrinsics (at least AVX necessary)
 * James W. Barnett
 * May 26, 2016
 * 
 * Limitations:
 *  - Must use a cubic box (does not have to be equal on all sides)
 *  - Particle being inserted has no charge, so electrostatics excluded
 *  - Every atom type with vdw interactions needs to have its own index group!
 *  - Isothermal-isobaric ensemble or NVT
 */

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include "gmxcpp/Trajectory.h"
#include "gmxcpp/Utils.h"

#include "omp.h"

#include "Atomtype.h"
#include "Ini.h"
#include "immintrin.h"

using namespace std;

double do_uncertainty(int boot_n, int block_n, int frame_total, double betai, vector <double> &V, vector <double> V_exp_pe);
double do_chempot(Trajectory &trj, Atomtype at[], vector <double> &V_exp_pe, vector <double> &V, int &frame_total, int frame_freq, int rand_n, double rand_ni, int chunk, double beta, double betai, int atomtypes);
void do_output(Ini &ini, Atomtype at[], double chem_pot, double chem_pot_uncertainty, time_t start_time, time_t end_time, int nthreds);

int main(int argc, char* argv[])
{

    if (argc != 2)
    {
        cout << "ERROR: Configuration file should be first command line argument." << endl;
        return -1;
    }

    chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now(); 
    time_t start_time = chrono::system_clock::to_time_t(start);

    cout << "Using '" << argv[1] << "' for configuration file..." << endl;
    Ini ini(argv[1]);

    int nthreads;
    #pragma omp parallel
    #pragma omp master
    nthreads = omp_get_num_threads();

    Trajectory trj(ini.xtcfile, ini.ndxfile);
    Atomtype at[ini.atomtypes];

    for (int i = 0; i < ini.atomtypes; i++)
    {
        at[i] = (Atomtype(trj, ini.atomtype_name[i], ini.testtype_sigma, ini.testtype_epsilon, ini.atomtype_sigma[i], ini.atomtype_epsilon[i], ini.rcut2, ini.epsfact));
    }

    vector <double> V_exp_pe;
    vector <double> V;
    int frame_total = 0;
    double chem_pot = do_chempot(trj, at, V_exp_pe, V, frame_total, ini.frame_freq, ini.rand_n, ini.rand_ni, nthreads*ini.chunk_size, ini.beta, ini.betai, ini.atomtypes);
    double chem_pot_uncertainty = do_uncertainty(ini.boot_n, ini.block_n, frame_total, ini.betai, V, V_exp_pe);

    cout << "---------------------------------------------------------" << endl;
    cout << "    μ (kJ / mol) = " << chem_pot << " ± " << chem_pot_uncertainty << endl;
    cout << "---------------------------------------------------------" << endl;

    end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);

    do_output(ini, at, chem_pot, chem_pot_uncertainty, start_time, end_time, nthreads);

    return 0;

}

double do_chempot(Trajectory &trj, Atomtype at[], vector <double> &V_exp_pe, vector <double> &V, int &frame_total, int frame_freq, int rand_n, double rand_ni, int chunk, double beta, double betai, int atomtypes)
{

    random_device rd;
    mt19937 gen(rd());
    double V_avg = 0.0;
    double V_exp_pe_avg = 0.0;

    int frames_read = -1;

    cout << left;

    // Read in frames chunk by chunk (saves on RAM).
    while (frames_read != 0)
    {
        cout << "Reading in " << chunk << " frames..." << endl;
        frames_read = trj.read_next(chunk);
        cout << "Finished reading in " << frames_read << " frames." << endl;

        // Expand for additional frames read in
        V_exp_pe.resize(V_exp_pe.size()+frames_read, 0.0);
        V.resize(V_exp_pe.size()+frames_read);

        // Now split up the frames that we read in to each processor
        #pragma omp parallel
        {

            cubicbox_m256 box;
            int thread_id;
            double pe;
            vector <coordinates> rand_xyz;

            #pragma omp for schedule(static) reduction(+:V_avg,V_exp_pe_avg)
            for (int frame_i = 0; frame_i < frames_read; frame_i++)
            {

                // frame_i is in relation to the current set of frames read in
                // so i is where we are in the total count of frames
                int i = frame_i + frame_total;

                // Gets the current box and saves it 8 times in memory for use
                // with SIMD instructions
                box = trj.GetCubicBoxM256(frame_i);

                gen_rand_box_points(rand_xyz, box, rand_n);
                V[i] = volume(box);

                for (int rand_i = 0; rand_i < rand_n; rand_i++)
                {
                    pe = 0.0;
                    for (int atomtype_i = 0; atomtype_i < atomtypes; atomtype_i++)
                    {
                        pe += at[atomtype_i].CalcPE(frame_i, trj, rand_xyz[rand_i], box, V[i]);
                    }

                    V_exp_pe[i] += exp(-pe * beta);
                }

                V_exp_pe[i] *= V[i] * rand_ni;

                // These averages are for this current thread, but will be
                // updated outside of the parallel region using reduction
                V_avg += V[i];
                V_exp_pe_avg += V_exp_pe[i];

                if (frame_i % frame_freq == 0)
                {
                    thread_id = omp_get_thread_num();
                    cout << "Thread: " << setw(3) << thread_id;
                    cout << " Frame: " << setw(6) << i;
                    cout << " μ = " << setw(12) << -log(V_exp_pe[i]/V[i]) * betai;
                    cout <<  "Thread <μ> = " << setw(12) << -log(V_exp_pe_avg/V_avg) * betai;
                    cout << endl;
                }

            }
        }

        frame_total += frames_read;
        cout << "Total <μ> = " << -log(V_exp_pe_avg/V_avg) * betai << endl;
    }

    return -log(V_exp_pe_avg/V_avg) * betai;
}

double do_uncertainty(int boot_n, int block_n, int frame_total, double betai, vector <double> &V, vector <double> V_exp_pe)
{
    random_device rd;
    mt19937 gen(rd());
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
            int block_start = (double)block/block_n * frame_total;

            if (block != block_n-1)
            {
                block_end = (double)(block+1)/block_n * frame_total;
            }
            else
            {
                block_end = frame_total;
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

    return sqrt(chem_pot_boot_var);
}

void do_output(Ini &ini, Atomtype at[], double chem_pot, double chem_pot_uncertainty, time_t start_time, time_t end_time, int nthreads)
{

    ofstream ofs(ini.outfile.c_str());
    ofs << scientific << setprecision(6) << left;
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << "         Test particle insertion program -- Wes Barnett" << endl;
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << setw(40) << "Started computation at:" << setw(20) << ctime(&start_time);
    ofs << setw(40) << "Finished computation at:" << setw(20) << ctime(&end_time);
    ofs << setw(40) << "Number of OMP threads:" << setw(20) << nthreads << endl;
    ofs << setw(40) << "Output file:" << setw(20) << ini.outfile << endl;
    ofs << setw(40) << "Compressed trajectory file:" << setw(20) << ini.xtcfile << endl;
    ofs << setw(40) << "Index file:" << setw(20) << ini.ndxfile << endl;
    ofs << setw(40) << "Insertions of particle per frame:" << setw(20) << ini.rand_n  << endl;
    ofs << setw(40) << "Blocks used in uncertainty analysis:" << setw(20) << ini.block_n << endl;
    ofs << setw(40) << "Bootstrapped iterations in uncertainty analysis:" << setw(20) << ini.boot_n << endl;
    ofs << setw(40) << "Well depth factor:" << setw(20) << ini.epsfact << endl;
    ofs << setw(40) << "Cutoff distance (nm):" << setw(20) << ini.rcut << endl;
    ofs << setw(40) << "System temperature (K):" << setw(20) << ini.T << endl;
    ofs << setw(40) << "Test particle sigma (nm):" << setw(20) << ini.testtype_sigma << endl;
    ofs << setw(40) << "Test particle epsilon (kJ/mol):" << setw(20) << ini.testtype_epsilon << endl;
    ofs << setw(40) << "Number of atom types:" << setw(20) << ini.atomtypes << endl;
    ofs << setw(20) << "INDEX NAME" << setw(20) << "SIGMA (nm)" << setw(20) << "EPSILON (kJ/mol)" << setw(20) << "C6" << setw(20) << "C12" << endl;
    for (int i = 0; i < ini.atomtypes; i++)
    {
        ofs << setw(20) << ini.atomtype_name[i] << setw(20) << ini.atomtype_sigma[i] << setw(20) << ini.atomtype_epsilon[i] << setw(20) << at[i].GetC6() << setw(20) << at[i].GetC12() << endl;
    }
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << "     FINAL RESULT - Excess chemical potential of test particle" << endl;
    ofs << "-----------------------------------------------------------------------" << endl;
    ofs << fixed;
    ofs << "μ (kJ / mol) = " << chem_pot << " ± " << chem_pot_uncertainty << endl;
    ofs.close();

    return;
}
