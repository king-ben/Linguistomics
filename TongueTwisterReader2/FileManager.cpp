#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include "FileManager.hpp"



FileManager::FileManager(std::string dn) : directoryName(std::move(dn)) {

    for (const auto& entry : std::filesystem::directory_iterator(directoryName))
        {
        std::filesystem::path fp = entry.path();
        std::string filePath = fp.string();
        std::string fileExtension = fp.extension().string();
        
        // Remove leading dot from extension
        if (!fileExtension.empty() && fileExtension[0] == '.')
            fileExtension.erase(0, 1);
        
        directoryFileContents.emplace_back(filePath, fileExtension);
        }
    std::sort(directoryFileContents.begin(), directoryFileContents.end());
}

std::string FileManager::getDirectoryBaseName(void) const {

    std::filesystem::path p(directoryName);
    if (p.filename().empty())
        p = p.parent_path();
    return p.filename().string();
}

std::vector<std::string> FileManager::filesWithExtension(const std::string& extension) const {

    std::vector<std::string> fileList;
    for (const auto& f : directoryFileContents)
        {
        if (f.second == extension)
            fileList.push_back(f.first);
        }
    return fileList;
}

void FileManager::print(void) const {

    size_t len = 0;
    for (const auto& f : directoryFileContents)
        {
        if (f.first.length() > len)
            len = f.first.length();
        }

    int n = 0;
    for (const auto& f : directoryFileContents)
        {
        std::cout << std::setw(4) << ++n << " -- " << f.first << " ";
        for (size_t i = 0; i < len - f.first.length(); i++)
            std::cout << " ";
        std::cout << f.second << std::endl;
        }
}
