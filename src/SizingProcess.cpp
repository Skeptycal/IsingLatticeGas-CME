#include "../include/SizingProcess.h"

SizingProcess::SizingProcess(Lattice* rLatticeObj, int clumptStartSize, int clumpMinSize)
{
    mpLatticeObj = rLatticeObj;

    mLatticeSize = mpLatticeObj->GetLatticeSize();

    mClumpStartSize = clumptStartSize;

    mClumpMinSize = clumpMinSize;
}

void SizingProcess::Run()
{
    std::vector<std::vector<int>>& rLattice = mpLatticeObj->GetLatticeReference();

    GenerateCheckedSizeVecs(mpLatticeObj->GetNumProteins());

    for (size_t i = 0; i < mCheckedProt.size(); i++)
    {
        std::vector<Protein*> inclusive_proteins;

        if (mCheckedProt[i])
        {
            continue;
        }

        SizeClumpRecursive(i, rLattice, inclusive_proteins);

        ClassifyClump(inclusive_proteins);
    }

}

void SizingProcess::SizeClumpRecursive(int proteinIndex, std::vector<std::vector<int>>& rLattice, std::vector<Protein*>& rInclusiveProteins)
{
    //Add reference to vector of clump inclusive proteins
    rInclusiveProteins.emplace_back(&(mpLatticeObj->GetProteinReference(proteinIndex)));

    int x, y;
    mpLatticeObj->GetProteinAddress(proteinIndex, x, y);
    mCheckedLatSite[x][y] = true;
    mCheckedProt[proteinIndex] = true;

    for (int i = 0; i < 4; i++)
    {
        int new_x, new_y;
        GetNeighborAddress(x, y, static_cast<Direction>(i), new_x, new_y);

        if (!mCheckedLatSite[new_x][new_y] &&
            rLattice[new_x][new_y])
        {
            SizeClumpRecursive(mpLatticeObj->GetProteinIndexAtLatSite(new_x, new_y), rLattice, rInclusiveProteins);
        }
    }
}

void SizingProcess::GenerateCheckedSizeVecs(int numProteins)
{
    mCheckedProt.clear();
    mCheckedProt = std::vector<bool>(numProteins, false);

    mCheckedLatSite.clear();
    for (int i = 0; i < mLatticeSize; i++)
    {
        std::vector<bool> l_chckd_lat_row(mLatticeSize, false);
        mCheckedLatSite.emplace_back(l_chckd_lat_row);
    }
}

void SizingProcess::GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY)
{
    //Implement Periodic Boundary Conditions
    if (dir == Direction::NORTH)
    {
        if (currX == 0)
        {
            newX = mLatticeSize - 1;
            newY = currY;
        }
        else
        {
            newX = currX - 1;
            newY = currY;
        }
    }
    else if (dir == Direction::EAST)
    {
        if (currY == (mLatticeSize - 1))
        {
            newX = currX;
            newY = 0;
        }
        else
        {
            newX = currX;
            newY = currY + 1;
        }
    }
    else if (dir == Direction::SOUTH)
    {
        if (currX == (mLatticeSize - 1))
        {
            newX = 0;
            newY = currY;
        }
        else
        {
            newX = currX + 1;
            newY = currY;
        }
    }
    else if (dir == Direction::WEST)
    {
        if (currY == 0)
        {
            newX = currX;
            newY = mLatticeSize - 1;
        }
        else
        {
            newX = currX;
            newY = currY - 1;
        }
    }
}

void SizingProcess::ClassifyClump(std::vector<Protein*>& inclusive_proteins)
{

}