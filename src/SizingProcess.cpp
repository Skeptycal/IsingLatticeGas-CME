#include "../include/SizingProcess.h"

SizingProcess::SizingProcess(Lattice* rLatticeObj, int clumptStartSize, int clumpMinSize)
{
    mpLatticeObj = rLatticeObj;

    mLatticeSize = mpLatticeObj->GetLatticeSize();

    mClumpStartSize = clumptStartSize;

    mClumpMinSize = clumpMinSize;
}

void SizingProcess::Run(int iteration)
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

        ClassifyClump(iteration, inclusive_proteins);
    }

}

void SizingProcess::SizeClumpRecursive(int proteinIndex, std::vector<std::vector<int>>& rLattice, std::vector<Protein*>& rInclusiveProteins)
{
    //Add reference to vector of clump inclusive proteins
    rInclusiveProteins.emplace_back(mpLatticeObj->GetProteinReference(proteinIndex));

    int x, y;
    mpLatticeObj->GetProteinAddress(proteinIndex, x, y);
    mCheckedLatSite[x][y] = true;
    mCheckedProt[proteinIndex] = true;

    for (int i = 0; i < 4; i++)
    {
        int new_x, new_y;
        GetNeighborAddress(x, y, static_cast<Direction>(i), new_x, new_y);

        if (!mCheckedLatSite[new_x][new_y] &&
            rLattice[new_x][new_y] == 1)
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

void SizingProcess::ClassifyClump(int iteration, std::vector<Protein*>& rInclusiveProteins)
{
    // Count the clump ids from every protein
    std::vector<int> clump_id_count(mpLatticeObj->GetNextClumpId() + 1, 0);
    for (auto itr : rInclusiveProteins)
    {
        clump_id_count[itr->GetClumpId() + 1]++;
    }

    // Determine the most frequent value and use that as the clump id
    int max_index = 0;
    for (size_t i = 1; i < clump_id_count.size(); i++)
    {
        if (clump_id_count[i] > clump_id_count[max_index])
        {
            max_index = i;
        }
    }

    int clump_id = max_index - 1;
    if (clump_id == -1) //Create a new clump if the id is 0 and above start size
    {
        if (rInclusiveProteins.size() >= mClumpStartSize)
        {
            mpLatticeObj->CreateNewClump(iteration, rInclusiveProteins);
        }
    }
    else // Otherwise abort or update the clump depending on size
    {
        if (rInclusiveProteins.size() < mClumpMinSize)
        {
            /*
                This could either happen when the clump dissolves or if a protein that was apart of a clump moves
                away from it. There is no way to identify the case so we reset the protein clump id and continue.
            */
            for (auto itr : rInclusiveProteins)
            {
                itr->SetClumpId(-1);
            }
        }
        else
        {
            mpLatticeObj->UpdateClump(iteration, clump_id, rInclusiveProteins);
        }
    }
}