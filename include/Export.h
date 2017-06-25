#include "Lattice.h"
#include <boost\filesystem.hpp>
#include "boost\date_time\posix_time\posix_time.hpp"

class Export
{
public:
    Export() {};

    Export(Lattice* rLatticeObj, std::string outputPath, int codeItr, double K, bool exportLattice, bool exportHistogram, bool exportClump, bool exportAmax);

    void Run(int itr);

    void WriteParameters(int minDissSize, int instMult, int clumpStrtSize, int clumpMinSize,
                         double rI, double rD, double deltaT, double K, int iterations);

private:

    void ExportLattice(int iteration);

    Lattice* mLatticeObj;

    bool mExportLattice;
    bool mExportHistogram;
    bool mExportClump;
    bool mExportAmax;

    std::string mOutputPath;
};
