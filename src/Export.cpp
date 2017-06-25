#include "../include/Export.h"
#include <iostream>
#include <sstream>
#include <fstream>

Export::Export(Lattice* rLatticeObj, std::string outputPath, int codeItr, double K, bool exportLattice, bool exportHistogram, bool exportClump, bool exportAmax)
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

    mLatticeObj = rLatticeObj;

    mExportLattice = exportLattice;
    mExportHistogram = exportHistogram;
    mExportClump = exportClump;
    mExportAmax = exportAmax;
}

void Export::Run(int iteration)
{
    if (mExportLattice)
    {
        ExportLattice(iteration);
    }
}

void Export::ExportLattice(int iteration)
{
    std::vector<std::vector<int>> lattice = mLatticeObj->GetLattice();

    std::ofstream lattice_file;
    lattice_file.open(mOutputPath + "\\Lattice_" + std::to_string(iteration) + ".txt");

    for (size_t row = 0; row < lattice.size(); row++)
    {
        for (size_t col = 0; col < lattice[0].size(); col++)
        {
            lattice_file << lattice[row][col] << ",";
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
