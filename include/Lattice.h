#ifndef LATTICE_H
#define LATTICE_H

#include <algorithm>
#include <vector>
#include "Protein.h"


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

    Lattice(int latticeSize, int numProteinInit, int insertionMultiplier, double insertionRate, double deltaTime, double k, int maxProteins);

    int GetNumProteins();

    bool CheckInsertion();

    void RunIteration();

    void GetLatticeCopy(std::vector<std::vector<int>>& rLattice)
    {
        rLattice = mLattice;
    };

    std::vector<std::vector<int>>& GetLatticeReference()
    {
        return mLattice;
    };

    Protein& GetProteinReference(int index)
    {
        return mProteins[index];
    };

    void GetProteinAddress(int index, int& x, int& y)
    {
        mProteins[index].GetAddress(x, y);
    };

    int GetLatticeSize()
    {
        return mLatticeSize;
    };

    int GetProteinIndexAtLatSite(int x, int y);

private:
    void GenerateLattice(int numProteinInit);

    double RandomSample();

    void GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY);

    Direction GetNextAddress(int x, int y, int& newX, int& newY);

    int CalculateEnergyChange(int currX, int currY, int newX, int newY);

    std::vector<std::vector<int>> mLattice;
    std::vector<std::vector<Protein*>> mProteinPtrLattice;
    std::vector<Protein> mProteins;

    int mLatticeSize;
    int mMaxProteins;
    int mNextProtID;

    double mInsertionProb;
    double mMvmtProb[3];
};

#endif //LATTICE_H