#include "Ini.h"

Ini::Ini(char file[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(file, pt);
    char *endptr;

    outfile = pt.get<std::string>("outfile","tpi.dat");
    xtcfile = pt.get<std::string>("xtcfile","prd.xtc");
    ndxfile = pt.get<std::string>("ndxfile","index.ndx");
    frame_freq = strtol(pt.get<std::string>("frame_freq","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'frame_freq' must be an integer." << endl;
    }
    chunk_size = strtol(pt.get<std::string>("chunk_size","1000").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'chunk_size' must be an integer." << endl;
    }
    rand_n = strtol(pt.get<std::string>("rand_n","1000").c_str(), &endptr, 10);
    rand_ni = 1.0/rand_n;
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rand_n' must be an integer." << endl;
    }
    block_n = strtol(pt.get<std::string>("block_n","5").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'block_n' must be an integer." << endl;
    }
    boot_n = strtol(pt.get<std::string>("boot_n","200").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'boot_n' must be an integer." << endl;
    }
    epsfact = strtod(pt.get<std::string>("epsfact","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'epsfact' must be a double." << endl;
    }
    rcut = strtod(pt.get<std::string>("rcut","1.0").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'rcut' must be a double." << endl;
    }
    rcut2 = rcut*rcut;
    T = strtod(pt.get<std::string>("T","298.15").c_str(), &endptr); // Kelvin
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'T' must be a double." << endl;
    }
    beta = 1.0/(R * T); // kJ/mol
    betai = (R * T);
    atomtypes = strtol(pt.get<std::string>("atomtypes","1").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0)
    {
        cout << "ERROR: 'atomtypes' must be an integer." << endl;
    }

    istringstream iss; 

    string testtype_str = pt.get<std::string>("testtype");
    iss.clear();
    iss.str(testtype_str);
    iss >> testtype_sigma;
    iss >> testtype_epsilon;

    atomtype_name.resize(atomtypes);
    atomtype_sigma.resize(atomtypes);
    atomtype_epsilon.resize(atomtypes);

    for (int i = 0; i < atomtypes; i++)
    {
        string atomtype_str = pt.get<std::string>("atomtype"+to_string(i+1));
        iss.clear();
        iss.str(atomtype_str);
        iss >> atomtype_name[i];
        iss >> atomtype_sigma[i];
        iss >> atomtype_epsilon[i];
    }

    cout << "outfile = "    << outfile << endl;
    cout << "xtcfile = "    << xtcfile << endl;
    cout << "ndxfile = "    << ndxfile << endl;
    cout << "frame_freq = " << frame_freq << endl;
    cout << "chunk_size = " << chunk_size << endl;
    cout << "rand_n = "     << rand_n << endl;
    cout << "block_n = "    << block_n << endl;
    cout << "boot_n = "     << boot_n << endl;
    cout << "epsfact = "    << epsfact << endl;
    cout << "rcut = "       << rcut << endl;
    cout << "T = "          << T << endl;
    cout << "atomtypes = "  << atomtypes << endl;
    for (int i = 0; i < atomtypes; i++)
    {
        cout << atomtype_name[i] << " " << atomtype_sigma[i] << " " << atomtype_epsilon[i] << endl;
    }

}
