#include <iostream>
#include <sstream>
#include <queue>
#include <utility>
#include <memory>
#include "../include/Export.h"
#include <boost/thread.hpp>
#include <boost/chrono.hpp>

const int NUM_PROTEIN_INIT = 0;
const int MIN_DISSOCIATION_SIZE = 60;
const int INSERTION_MULTIPLIER = 30;
const int DISSOCIATION_INTERVAL = 2500;
const int EXPORT_INTERVAL = 500;
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

const int THREAD_SLEEP_TIMEOUT = 10;

double gK;
int gIterations;
int gCodeItr;
bool gExportInit;
bool gRunning;
std::queue<std::pair<int, std::vector<std::vector<int>>>> gExportQueue;
std::shared_ptr<Lattice> gpLattice;
std::string gOutputPath("");
boost::thread gExpThread;


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

void ExportWorkerThread()
{
    Export exp = Export(gpLattice.get(), gOutputPath, gCodeItr, gK, EXPORT_LATTICE, EXPORT_HISTOGRAM, EXPORT_CLUMP, EXPORT_AMAX);
    exp.WriteParameters(MIN_DISSOCIATION_SIZE, INSERTION_MULTIPLIER, CLUMP_START_SIZE, CLUMP_MIN_SIZE, R_I, R_D, DELTA_T, gK, gIterations);

    while (gRunning || !gExportQueue.empty())
    {
        if (gExportQueue.empty())
        {
            boost::this_thread::sleep_for(boost::chrono::milliseconds(THREAD_SLEEP_TIMEOUT));
        }
        else
        {
            exp.Run(gExportQueue.front().first, gExportQueue.front().second);
            gExportQueue.pop();
        }
    }
}

int main(int argc, const char* argv[])
{
    if (!ParseArguments(argc, argv))
    {
        std::cout << "Arguments could not be parsed";
        return 0;
    }
    gRunning = true;

    gpLattice.reset(new Lattice(LATTICE_SIZE, NUM_PROTEIN_INIT, INSERTION_MULTIPLIER, R_I, DELTA_T, gK, MAX_PROTEINS));

    gExpThread = boost::thread(ExportWorkerThread);

    for (int itr = 0; itr < gIterations; itr++)
    {
        std::cout << "Itr: " << itr << std::endl;

        if (itr % EXPORT_INTERVAL == 0)
        {
            std::vector<std::vector<int>> l_lattice;
            gpLattice->GetLattice(l_lattice);
            gExportQueue.push(std::make_pair(itr, l_lattice));
        }

        gpLattice->CheckInsertion();

        gpLattice->RunIteration();
    }

    {
        std::vector<std::vector<int>> l_lattice;
        gpLattice->GetLattice(l_lattice);
        gExportQueue.push(std::make_pair(gIterations, l_lattice));
    }

    gRunning = false;

    if (gExpThread.joinable())
    {
        gExpThread.join();
    }
}
