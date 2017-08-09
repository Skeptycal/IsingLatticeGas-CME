#ifndef SIZING_PROCESS_H
#define SIZING_PROCESS_H

#include "Lattice.h"
#include "Protein.h"

class SizingProcess
{
public:
    SizingProcess(Lattice* rLatticeObj, int clumptStartSize, int clumpMinSize);
    ~SizingProcess() {};

    void Run(int iteration);
private:

    void GenerateCheckedSizeVecs(int numProteins);

    void SizeClumpRecursive(int proteinIndex, std::vector<std::vector<int>>& rLattice, std::vector<Protein*>& rInclusiveProteins);

    void GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY);

    void ClassifyClump(int iteration, std::vector<Protein*>& inclusive_proteins);

    Lattice* mpLatticeObj;

    std::vector<std::vector<bool>> mCheckedLatSite;
    std::vector<bool> mCheckedProt;

    int mLatticeSize;
    int mClumpStartSize;
    int mClumpMinSize;
};

#endif //SIZING_PROCESS_H