
class Protein
{
public:
    Protein() {};

    Protein(int x, int y);

    void GetAddress(int& x, int& y);

    void SetAddress(int x, int y);

private:

    int mX;
    int mY;

};
