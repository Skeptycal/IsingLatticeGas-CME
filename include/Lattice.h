#ifndef LATTICE_H
#define LATTICE_H

#include <algorithm>
#include <vector>
#include "Protein.h"
#include "Clump.h"


/* Coordinate System

    y-->
   x[][][][]....
   |[][][][]....
   v[][][][]....
    [][][][]....
    : : : :
    : : : :

*/

enum Direction
{
    NORTH = 0,
    EAST,
    SOUTH,
    WEST,
    NONE,
};

class Lattice
{
public:
    Lattice() {};

    Lattice(int latticeSize, int numProteinInit, int insertionMultiplier, int dissociationInterval, double dissociationRate, double insertionRate, double deltaTime, double k, int maxProteins, int minDissociationSize);
    ~Lattice();

    int GetNumProteins();

    bool CheckInsertion();

    void CheckDissociation();

    void RunIteration();

    void GetLatticeCopy(std::vector<std::vector<int>>& rLattice)
    {
        rLattice = mLattice;
    };

    std::vector<std::vector<int>>& GetLatticeReference()
    {
        return mLattice;
    };

    Protein* GetProteinReference(int index)
    {
        return mProteins[index];
    };

    std::vector<Clump*>& GetClumpReferences()
    {
        return mClumps;
    };

    void GetProteinAddress(int index, int& x, int& y)
    {
        mProteins[index]->GetAddress(x, y);
    };

    int GetLatticeSize()
    {
        return mLatticeSize;
    };

    int GetNextClumpId()
    {
        return mNextClumpID;
    };

    int GetProteinIndexAtLatSite(int x, int y);

    void CreateNewClump(int iteration, std::vector<Protein*>& rInclusiveProteins);
    void UpdateClump(int iteration, int id, std::vector<Protein*>& rInclusiveProteins);

    double RandomSample();

private:
    void GenerateLattice(int numProteinInit);

    void GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY);

    Direction GetNextAddress(int x, int y, int& newX, int& newY);

    int CalculateEnergyChange(int currX, int currY, int newX, int newY);

    std::vector<std::vector<int>> mLattice;
    std::vector<std::vector<Protein*>> mProteinPtrLattice;
    std::vector<Protein*> mProteins;
    std::vector<Clump*> mClumps;

    int mLatticeSize;
    int mMaxProteins;
    int mNextProtID;
    int mNextClumpID;
    int mDissociationInterval;
    int mMinDissociationSize;

    double mInsertionProb;
    double mMvmtProb[3];
    double mDissociationProb;
};

#endif //LATTICE_H