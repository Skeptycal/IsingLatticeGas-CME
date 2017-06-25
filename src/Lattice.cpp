#include "../include/Lattice.h"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <math.h>

Lattice::Lattice(int latticeSize, int numProteinInit, int insertionMultiplier, double insertionRate, double deltaTime, double k, int maxProteins)
{
    mLatticeSize = latticeSize;
    mMaxProteins = maxProteins;

    //Initialize lattice
    GenerateLattice(numProteinInit);

    //Calculate insertion probability
    mInsertionProb = static_cast<double>(insertionMultiplier) * insertionRate * deltaTime;

    //Calculate movement probability
    mMvmtProb[0] = std::exp((-1.0 * k) * 1.0);
    mMvmtProb[1] = std::exp((-1.0 * k) * 2.0);
    mMvmtProb[2] = std::exp((-1.0 * k) * 3.0);
}

void Lattice::GenerateLattice(int numProteinInit)
{
    int prot_to_insert = numProteinInit;

    //Initialize lattice to zero
    for (int i = 0; i < mLatticeSize; i++)
    {
        std::vector<int> l_lat_row(mLatticeSize, 0);
        mLattice.emplace_back(l_lat_row);
    }

    //Seed rand with system time
    std::srand(static_cast<int>(time(0)));

    //Randomly insert protein to the lattice
    while (prot_to_insert > 0)
    {
        int x = std::rand() % mLatticeSize ;
        int y = std::rand() % mLatticeSize ;
        if (mLattice[x][y] == 0)
        {
            Protein l_protein(x, y);
            mProteins.emplace_back(l_protein);
            mLattice[x][y] = 1;
            prot_to_insert--;
        }
    }
}

int Lattice::GetNumProteins()
{
    int total = 0;

    for (size_t i = 0; i < mLattice.size(); i++)
    {
        for (size_t j = 0; j < mLattice[0].size(); j++)
        {
            if (mLattice[i][j] > 0)
            {
                total++;
            }
        }
    }

    if (total != mProteins.size())
    {
        std::cout << "Protein vector disagrees with proteins in lattice." << std::endl;
    }

    return total;
}

bool Lattice::CheckInsertion()
{
    //Test insertion by sampling a random number between 0..1
    if (mProteins.size() < mMaxProteins &&
        mInsertionProb > RandomSample())
    {
        while (true)
        {
            int x = std::rand() % mLatticeSize;
            int y = std::rand() % mLatticeSize;
            if (mLattice[x][y] == 0)
            {
                Protein l_protein(x, y);
                mProteins.emplace_back(l_protein);
                mLattice[x][y] = 1;
                return true;
            }
        }
    }

    return false;
}

void Lattice::RunIteration()
{
    for (size_t i = 0; i < mProteins.size(); i++)
    {
        Direction l_direction;
        int curr_x, curr_y;
        int new_x, new_y;
        mProteins[i].GetAddress(curr_x, curr_y);

        // Find new direction for protein to move
        l_direction = GetNextAddress(curr_x, curr_y, new_x, new_y);

        if (l_direction == Direction::NONE)
        {
            continue;
        }
        else
        {
            // Calculate the energy change for the protein movement
            int d_U = CalculateEnergyChange(curr_x, curr_y, new_x, new_y);
            if (d_U <= 0 ||
                mMvmtProb[d_U - 1] > RandomSample())
            {
                mLattice[curr_x][curr_y] = 0;
                mLattice[new_x][new_y] = 1;
                mProteins[i].SetAddress(new_x, new_y);
            }
        }
    }
}

double Lattice::RandomSample()
{
    return static_cast<double>(rand()) / RAND_MAX;
}

auto Lattice::GetNextAddress(int x, int y, int& newX, int& newY) -> Direction
{
    Direction dir;
    bool chckd_direction[4] = { false, false, false, false };
    int new_x, new_y;

    for (int i = 0; i < 4; i)
    {
        dir = static_cast<Direction>(std::rand() % 4);

        if (chckd_direction[dir] == true)
        {
            continue;
        }

        GetNeighborAddress(x, y, dir, new_x, new_y);

        if (mLattice[new_x][new_y] == 1)
        {
            chckd_direction[dir] = true;
            i++;
        }
        else
        {
            newX = new_x;
            newY = new_y;
            return dir;
        }
    }

    newX = x;
    newY = y;
    return Direction::NONE;
}

void Lattice::GetNeighborAddress(int currX, int currY, Direction dir, int& newX, int& newY)
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

int Lattice::CalculateEnergyChange(int currX, int currY, int newX, int newY)
{
    int x[2] = { currX, newX };
    int y[2] = { currY, newY };
    int u[2] = { 0, 1 }; //u[1] starts with 1 to account for the protein no longer being at its previous location
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            int l_x, l_y;
            GetNeighborAddress(x[i], y[i], static_cast<Direction>(j), l_x, l_y);
            u[i] -= mLattice[l_x][l_y];
        }
    }
    return u[1] - u[0];
}
