#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>
#include <vector>
enum SubstitutionModel { jc69, gtr, custom };



class UserSettings {

    public:
        static UserSettings&        userSettings(void) {
                                        static UserSettings us;
                                        return us;
                                    }
        bool                        getFullOutput(void) { return fullOutput; }
        std::vector<std::string>    getInputDirectoryName(void) { return inFile; }
        std::string                 getOutFile(void) { return outFile; }
        std::string                 getNytrilOutputFileName(void) { return nytrilOutFile; }
        std::string                 getROutFile(void) { return rOutFile; }
        double                      getBurnFraction(void) { return burnFraction; }
        void                        print(void);
        void                        readCommandLineArguments(int argc, char* argv[]);
    
    private:
                                    UserSettings(void);
                                    UserSettings(const UserSettings& s) = delete;
        void                        usage(void);
        std::vector<std::string>    inFile;
        std::string                 outFile;
        std::string                 nytrilOutFile;
        std::string                 rOutFile;
        double                      burnFraction;
        bool                        fullOutput;
};

#endif
