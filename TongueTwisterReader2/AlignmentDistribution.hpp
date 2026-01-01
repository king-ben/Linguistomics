#ifndef AlignmentDistribution_hpp
#define AlignmentDistribution_hpp

#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include "Alignment.hpp"
#include "json.hpp"
class Partition;
class RandomVariable;


class AlignmentDistribution {

    public:
                                                    AlignmentDistribution(void) = delete;
                                                    AlignmentDistribution(RandomVariable* r, std::vector<Alignment*>& alns, std::string n);
                                                   ~AlignmentDistribution(void);
        AlignmentDistribution&                      operator+=(const AlignmentDistribution& rhs);
        void                                        addAlignment(Alignment* aln);
        void                                        addAlignmentCopy(Alignment* aln);
        const std::map<Alignment*,int>&             getAlignmentSamples(void) const { return samples; }
        std::map<int,std::pair<Alignment*,float>>   getSortedAlignments(void);
        int                                         ciSize(void);
        std::string                                 getName(void) { return name; }
        Alignment*                                  getMapAlignment(void);
        int                                         numSamples(void);
        void                                        print(void);
        void                                        print(bool showAlignments);
        void                                        print(bool showAlignments, Partition* part);
        void                                        setName(std::string s) { name = s; }
        size_t                                      size(void) { return samples.size(); }
        nlohmann::json                              toJson(double credibleSetSize, int maxAlignment, std::ostream& file);

    private:
        void                                        print(Alignment* aln, Partition* part);
        RandomVariable*                             rv;
        std::string                                 name;
        std::map<Alignment*,int>                    samples;
};

#endif
