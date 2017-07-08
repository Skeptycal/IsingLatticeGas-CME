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


class Lattice
{
public:
    enum Direction
    {
        NORTH = 0,
        EAST,
        SOUTH,
        WEST,
        NONE,
    };

    Lattice() {};

    Lattice(int latticeSize, int numProteinInit, int insertionMultiplier, double insertionRate, double deltaTime, double k, int maxProteins);

    int GetNumProteins();

    bool CheckInsertion();

    void RunIteration();

    std::vector<std::vector<int>>& GetLattice()
    {
        return mLattice;
    };

private:
    void GenerateLattice(int numProteinInit);

    double RandomSample();

    void GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY);

    Direction GetNextAddress(int x, int y, int& newX, int& newY);

    int CalculateEnergyChange(int currX, int currY, int newX, int newY);

    std::vector<std::vector<int>> mLattice;
    std::vector<Protein> mProteins;

    int mLatticeSize;
    int mMaxProteins;

    double mInsertionProb;
    double mMvmtProb[3];
};
