#include <filesystem>
#include <iomanip>
#include <iostream>
#include "FileManager.hpp"



FileManager::FileManager(std::string dn) : directoryName(dn) {

    for (const auto& entry : std::filesystem::directory_iterator(directoryName))
        {
        std::filesystem::path fp = entry.path();
        #ifdef _CONSOLE
        auto filePath      = fp.string();
        auto fileExtension = fp.extension();
        #else
        std::string filePath      = fp;
        std::string fileExtension = fp.extension();
        #endif
        
        std::erase(fileExtension, '.');
        std::pair<std::string,std::string> nameExt;
        nameExt.first = filePath;
        nameExt.second = fileExtension;
        directoryFileContents.push_back(nameExt);
        }
}

std::vector<std::string> FileManager::filesWithExtension(std::string extension) {

    std::vector<std::string> fileList;
    for (auto f : directoryFileContents)
        {
        if (f.second == extension)
            fileList.push_back(f.first);
        }
    return fileList;
}

void FileManager::print(void) {

    size_t len = 0;
    for (auto f : directoryFileContents)
        {
        if (f.first.length() > len)
            len = f.first.length();
        }

    int n = 0;
    for (auto f : directoryFileContents)
        {
        std::cout << std::setw(4) << ++n << " -- " << f.first << " ";
        for (size_t i=0; i<len-f.first.length(); i++)
            std::cout << " ";
        std::cout << f.second << std::endl;
        }
}
