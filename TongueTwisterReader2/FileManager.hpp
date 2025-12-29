#ifndef FileManager_hpp
#define FileManager_hpp

#include <string>
#include <vector>

typedef std::vector<std::pair<std::string,std::string>> DirectoryList;



class FileManager {

    public:
                                    FileManager(void) = delete;
                                    FileManager(std::string dn);
        std::vector<std::string>    filesWithExtension(std::string extension);
        void                        print(void);
    
    private:
        std::string                 directoryName;
        DirectoryList               directoryFileContents;
};

#endif
