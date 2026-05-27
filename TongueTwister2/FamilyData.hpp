#pragma once
#include <string>
#include <vector>
#include "json.hpp"

class FamilyData {

    public:
                        FamilyData(void) = default;
        void            loadFromFile(const std::string& path);
        void            loadFromJson(const nlohmann::json& j);
        
        template<typename T>
        T               getValue(const std::string& key) const
                            { return data.at(key).get<T>(); }
        bool            hasKey(const std::string& key) const
                            { return data.find(key) != data.end(); }
        const nlohmann::json* getJsonPointer(const std::string& key) const;
        std::string     dump(void) const { return data.dump(4); }
        const std::string& getFilePath(void) const { return filePath; }

    private:
        nlohmann::json  data;
        std::string     filePath;
};