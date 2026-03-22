#ifndef FileManager_hpp
#define FileManager_hpp

#include <string>
#include <vector>

typedef std::vector<std::pair<std::string,std::string>> DirectoryList;



class FileManager {

    public:
                                    FileManager(void) = delete;
                                    FileManager(std::string dn);
        std::string                 getDirectoryBaseName(void) const;
        std::string                 getDirectoryPath(void) const { return directoryName; }
        std::vector<std::string>    filesWithExtension(const std::string& extension) const;
        size_t                      numFiles(void) const { return directoryFileContents.size(); }
        void                        print(void) const;
        static bool                 directoryExists(const std::string& path);
        static bool                 ensureDirectoryExists(const std::string& path);
    
    private:
        std::string                 directoryName;
        DirectoryList               directoryFileContents;
};

#endif
