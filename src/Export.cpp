#include "../include/Export.h"
#include <iostream>
#include <sstream>
#include <fstream>

Export::Export(Lattice* rLatticeObj, std::string outputPath, int codeItr, int maxIteration, double K, bool exportLattice, bool exportHistogram, bool exportClump, bool exportAmax)
{
    auto facet = new boost::posix_time::time_facet("%Y-%m-%d-T%H%M%S");

    std::stringstream out_path_ss;
    out_path_ss.imbue(std::locale(std::locale::classic(), facet));
    out_path_ss << outputPath << "K_" << K << "_";
    out_path_ss << boost::posix_time::second_clock::local_time();
    out_path_ss << "_" << codeItr;
    mOutputPath = out_path_ss.str();

    boost::filesystem::path dir(mOutputPath);
    if (boost::filesystem::create_directory(dir))
    {
        std::cout << "Created directory";
    }

    mpLatticeObj = rLatticeObj;

    mExportLattice = exportLattice;
    mExportHistogram = exportHistogram;
    mExportClump = exportClump;
    mExportAmax = exportAmax;

    mMaxIteration = maxIteration;
}

void Export::Run(int iteration, std::vector<std::vector<int>>& rLattice)
{
    if (mExportLattice)
    {
        ExportLattice(iteration, rLattice);
    }
}

void Export::ExportLattice(int iteration, std::vector<std::vector<int>>& rLattice)
{
    std::ofstream lattice_file;
    lattice_file.open(mOutputPath + "\\Lattice_" + std::to_string(iteration) + ".txt");

    for (size_t row = 0; row < rLattice.size(); row++)
    {
        for (size_t col = 0; col < rLattice[0].size(); col++)
        {
            lattice_file << rLattice[row][col] << ",";
        }
        lattice_file << "\n";
    }

    lattice_file.close();
}

void Export::WriteParameters(int minDissSize, int instMult, int clumpStrtSize, int clumpMinSize,
                             double rI, double rD, double deltaT, double K, int iterations)
{
    std::ofstream param_file;
    param_file.open(mOutputPath + "\\Parameters.txt");

    param_file << "Minimum dissociation size \t" << minDissSize << "\n";
    param_file << "Insertion multiplier \t" << instMult << "\n";
    param_file << "Clump start size \t" << clumpStrtSize << "\n";
    param_file << "Clump minimum size \t" << clumpMinSize << "\n";
    param_file << "Insertion rate \t" << rI << "\n";
    param_file << "Dissociation rate \t" << rD << "\n";
    param_file << "Time interval \t" << deltaT << "\n";
    param_file << "K \t" << K << "\n";
    param_file << "Iterations \t" << iterations << "\n";

    param_file.close();
}

void Export::WriteClumps()
{
    std::vector<Clump*> clump_vector = mpLatticeObj->GetClumpReferences();

    std::ofstream temp;
    temp.open(mOutputPath + "\\Clump.txt");

    for (auto clump : clump_vector)
    {
        std::vector<std::pair<int, int>> size_history = clump->GetSizeHistoryReference();

        std::ofstream clump_file;
        clump_file.open(mOutputPath + "\\Clump_" + std::to_string(clump->GetID()) + ".txt");

        switch (clump->GetClumpState())
        {
            case ClumpState::PRODUCTIVE:
            {
                clump_file << "productive\n";
                temp << "productive\n";
                break;
            }
            case ClumpState::ACTIVE:
            {
                if (size_history[size_history.size() - 1].first == mMaxIteration)
                {
                    clump_file << "active\n";
                    temp << "active\n";
                }
                else
                {
                    clump_file << "aborted\n";
                    temp << "aborted\n";
                }
                break;
            }
            default:
                break;
        }

        for (auto size : size_history)
        {
            clump_file << size.first << "," << size.second << "\n";
        }

        clump_file.close();
    }
}