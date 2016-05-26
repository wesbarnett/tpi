/* 
 * Test particle insertion
 * James W. Barnett
 * May 26, 2016
 * 
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "gmxcpp/Trajectory.h"
using namespace std;

int main(int argc, char* argv[])
{
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

    return 0;

}
