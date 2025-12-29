#ifndef TreeSamples_hpp
#define TreeSamples_hpp

#include <vector>
class McmcSummary;
class Tree;



class TreeSamples {

    public:
                            TreeSamples(void) = delete;
                            TreeSamples(McmcSummary& samples, double bf);
                           ~TreeSamples(void);
        void                print(void);
   
    private:
        std::vector<Tree*>  trees;
};

#endif
