#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>
#include "json.hpp"


enum SubstitutionModel { jc69, f81, gtr, custom };
enum MatrixExpBackend { autoBackend, cpuThreaded, cpuBatched, gpuBatched };

class UserSettings {

    public:
        static UserSettings&    userSettings(void) {
                                    static UserSettings us;
                                    return us;
                                }
        bool                    getCalculateMarginalLikelihood(void) { return calculateMarginalLikelihood; }
        int                     getCheckPointFrequency(void) { return checkPointFrequency; }
        std::string             getDataFile(void) { return dataFile; }
        double                  getBranchLengthLambda(void) { return branchLengthLambda; }
        int                     getBurnLength(void) { return burnLength; }
        int                     getFirstBurnLength(void) { return firstBurnLength; }
        MatrixExpBackend        getMatrixExpBackend(void) { return matrixExpBackend; }
        std::string             getOutFile(void) { return outFile; }
        int                     getNumChains(void) { return numChains; }
        int                     getNumMcmcCycles(void) { return numMcmcCycles; }
        int                     getNumIndelCategories(void) { return numIndelCategories; }
        int                     getNumRateCategories(void) { return numRateCategories; }
        int                     getPreburninLength(void) { return preburninLength; }
        int                     getPrintFrequency(void) { return printFrequency; }
        int                     getSampleFrequency(void) { return sampleFrequency; }
        int                     getSampleLength(void) { return sampleLength; }
        uint32_t                getSeed(void) { return seed; }
        SubstitutionModel       getSubstitutionModel(void) { return substitutionModel; }
        double                  getTemperature(void) { return temperature; }
        int                     getTuneLength(void) { return tuneLength; }
        void                    print(void);
        void                    readCommandLineArguments(int argc, char* argv[]);
        bool                    getUseOnlyCompleteWords(void) { return useOnlyCompleteWords; }
        bool                    getUseClockConstraint(void) { return useClockConstraint; }
        bool                    getSampleAlignments(void) { return sampleAlignments; }
    
    private:
                                UserSettings(void);
                                UserSettings(const UserSettings& s) = delete;
        void                    readJsonSettings(void);
        std::string             getVariable(nlohmann::json& settings, const char *name);
        std::string             outFile;
        std::string             dataFile;
        std::string             executablePath;
        double                  branchLengthLambda;
        int                     checkPointFrequency;
        int                     numIndelCategories;
        int                     numMcmcCycles;
        int                     numRateCategories;
        int                     printFrequency;
        SubstitutionModel       substitutionModel;
        MatrixExpBackend        matrixExpBackend;
        bool                    useOnlyCompleteWords;
        bool                    calculateMarginalLikelihood;
        int                     preburninLength;
        uint32_t                seed;
        int                     numChains;
        double                  temperature;
        bool                    useClockConstraint;
        int                     firstBurnLength;
        int                     tuneLength;
        int                     burnLength;
        int                     sampleLength;
        int                     sampleFrequency;
        bool                    sampleAlignments;
};

#endif
