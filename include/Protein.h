#ifndef PROTEIN_H
#define PROTEIN_H

class Protein
{
public:
    Protein() {};

    Protein(int x, int y, int id);

    void GetAddress(int& x, int& y);
    void SetAddress(int x, int y);

    int GetId();

    int GetClumpId();
    void SetClumpId(int id);

private:

    int mID;
    int mX;
    int mY;
    int mClumpID;
};

#endif //PROTEIN_H