#include "../include/Protein.h"

Protein::Protein(int x, int y)
{
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
