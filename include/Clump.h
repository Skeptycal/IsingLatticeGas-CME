#ifndef CLUMP_H
#define CLUMP_H

#include <vector>

enum ClumpState
{
    ACTIVE = 0,
    ABORTED,
    PRODUCTIVE,
};

class Clump
{
public:
    Clump(int id, int iteration, int startSize);

    ~Clump() {};

    bool UpdateSize(int iteration, int size);

    void SetClumpState(ClumpState state);

    int GetClumpSize();

    int GetID()
    {
        return mID;
    };

    ClumpState GetClumpState()
    {
        return mClumpState;
    };

    std::vector<std::pair<int, int>>& GetSizeHistoryReference()
    {
        return mLifetimeSize;
    };

private:

    int mID;
    ClumpState mClumpState;
    std::vector<std::pair<int, int>> mLifetimeSize;
};

#endif //CLUMP_H