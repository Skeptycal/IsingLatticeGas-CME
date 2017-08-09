#include "../include/Clump.h"
#include <iostream>

Clump::Clump(int id, int iteration, int startSize)
{
    mID = id;
    mClumpState = ClumpState::ACTIVE;
    mLifetimeSize.emplace_back(std::make_pair(iteration, startSize));
}

bool Clump::UpdateSize(int iteration, int size)
{
    // Checks to see if this clump is already classified as a completed clump
    if (mClumpState != ClumpState::ACTIVE)
    {
        switch (mClumpState)
        {
            case ClumpState::ABORTED:
                std::cout << "Clump has already aborted. New clump will be created." << std::endl;
                break;
            case ClumpState::PRODUCTIVE:
                std::cout << "Clump has already dissociated. New clump will be created." << std::endl;
                break;
            default:
                break;
        }
        return false;
    }

    mLifetimeSize.emplace_back(std::make_pair(iteration, size));
    return true;
}

void Clump::SetClumpState(ClumpState state)
{
    mClumpState = state;
}

int Clump::GetClumpSize()
{
    return mLifetimeSize[mLifetimeSize.size() - 1].second;
}