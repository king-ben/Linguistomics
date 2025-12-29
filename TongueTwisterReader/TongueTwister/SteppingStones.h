#ifndef SteppingStones_H
#define SteppingStones_H

#include <vector>
class Stone;



class SteppingStones {

    public:
                            SteppingStones(std::vector<double> pows);
                           ~SteppingStones(void);
        void                addSample(int powIdx, double lnL);
        void                clearSamples(int powIdx);
        double              marginalLikelihood(void);
        int                 numSamplesPerStone(void);
        int                 numStones(void) { return (int)stones.size(); }

    protected:
        std::vector<Stone*> stones;
};

#endif
