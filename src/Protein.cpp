#include "../include/Protein.h"

Protein::Protein(int x, int y, int id)
{
    mID = id;
    mClumpID = -1;
    mX = x;
    mY = y;
}

void Protein::GetAddress(int& x, int& y)
{
    x = mX;
    y = mY;
}

void Protein::SetAddress(int x, int y)
{
    mX = x;
    mY = y;
}

int Protein::GetId()
{
    return mID;
}

int Protein::GetClumpId()
{
    return mClumpID;
}

void Protein::SetClumpId(int id)
{
    mClumpID = id;
}