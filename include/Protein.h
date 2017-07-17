#ifndef PROTEIN_H
#define PROTEIN_H

class Protein
{
public:
    Protein() {};

    Protein(int x, int y, int id);

    void GetAddress(int& x, int& y);

    void SetAddress(int x, int y);

    int GetId()
    {
        return mID;
    }

private:

    int mID;
    int mX;
    int mY;

};

#endif //PROTEIN_H