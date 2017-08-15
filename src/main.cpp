#include <iostream>
#include <sstream>
#include "../include/Lattice.h"
#include "../include/Export.h"
#include "../include/SizingProcess.h"
#include <boost\thread.hpp>

const int NUM_PROTEIN_INIT = 0;
const int MIN_DISSOCIATION_SIZE = 60;
const int INSERTION_MULTIPLIER = 30;
const int DISSOCIATION_INTERVAL = 1250;
const int EXPORT_INTERVAL = 5000;
const int CLUMP_START_SIZE = 10;
const int CLUMP_MIN_SIZE = 3;
const int MAX_PROTEINS = 20000;
const int LATTICE_SIZE = 800;

const double R_I = 0.2;
const double R_D = 0.5;
const double DELTA_T = 0.0002328;

const bool EXPORT_LATTICE = 1;
const bool EXPORT_HISTOGRAM = 0;
const bool EXPORT_CLUMP = 0;
const bool EXPORT_AMAX = 0;

double gK;
int gIterations;
int gCodeItr;
std::string gOutputPath("");
boost::thread* gpExpThread;


bool ParseArguments(int argc, const char* argv[])
{
    int paramNum = 2; //Min number of parameters

    if (argc < paramNum * 2)
    {
        std::cout << "Not enough cmd line inputs" << std::endl;
    }
    else
    {
        for (int i = 1; i < argc;)
        {
            if (strcmp(argv[i], "-k") == 0)
            {
                gK = std::stod(argv[i + 1]);
                i += 2;
            }
            else if (strcmp(argv[i], "-i") == 0)
            {
                gIterations = std::atoi(argv[i + 1]);
                i += 2;
            }
            else if (strcmp(argv[i], "-n") == 0)
            {
                gCodeItr = std::atoi(argv[i + 1]);
                i += 2;
            }
            else if (strcmp(argv[i], "-o") == 0)
            {
                gOutputPath = argv[i + 1];
                i += 2;
            }
            else
            {
                std::cout << "Invalid input parameter" << std::endl;
                return false;
            }
        }
    }

    return true;
}

int main(int argc, const char* argv[])
{
    if (!ParseArguments(argc, argv))
    {
        std::cout << "Arguments could not be parsed";
        return 0;
    }

    Lattice lattice = Lattice(LATTICE_SIZE, NUM_PROTEIN_INIT, INSERTION_MULTIPLIER, DISSOCIATION_INTERVAL, R_D, R_I, DELTA_T, gK, MAX_PROTEINS, MIN_DISSOCIATION_SIZE);
    SizingProcess sizing_process = SizingProcess(&lattice, CLUMP_START_SIZE, CLUMP_MIN_SIZE);
    Export exp = Export(&lattice, gOutputPath, gCodeItr, gIterations, gK, EXPORT_LATTICE, EXPORT_HISTOGRAM, EXPORT_CLUMP, EXPORT_AMAX);
    exp.WriteParameters(MIN_DISSOCIATION_SIZE, INSERTION_MULTIPLIER, CLUMP_START_SIZE, CLUMP_MIN_SIZE, R_I, R_D, DELTA_T, gK, gIterations);

    for (int itr = 0; itr < gIterations; itr++)
    {
        std::cout << "Itr: " << itr << std::endl;

        if (itr % EXPORT_INTERVAL == 0)
        {
            if (gpExpThread != nullptr)
            {
                gpExpThread->join();
                delete gpExpThread;
            }

            //Need to copy lattice before spinning a thread so the lattice is unchanged
            std::vector<std::vector<int>> l_lattice;
            lattice.GetLatticeCopy(l_lattice);
            gpExpThread = new boost::thread(&Export::Run, &exp, itr, l_lattice);
        }

        if (itr % DISSOCIATION_INTERVAL == 0)
        {
            sizing_process.Run(itr);
            lattice.CheckDissociation();
        }

        lattice.CheckInsertion();

        lattice.RunIteration();
    }

    exp.WriteClumps();

    {
        std::vector<std::vector<int>> l_lattice;
        lattice.GetLatticeCopy(l_lattice);
        exp.Run(gIterations, l_lattice);
    }
}
