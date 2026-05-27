#include <fstream>
#include "FamilyData.hpp"
#include "Msg.hpp"

void FamilyData::loadFromFile(const std::string& path) {

    filePath = path;
    std::ifstream ifs(path);
    if (!ifs.is_open())
        Msg::error("Could not open family data file: " + path);
    try
        {
        data = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::error("Error parsing JSON file \"" + path + "\" at byte "
                   + std::to_string(ex.byte));
        }
}

void FamilyData::loadFromJson(const nlohmann::json& j) {

    data = j;
}

const nlohmann::json* FamilyData::getJsonPointer(const std::string& key) const {

    auto it = data.find(key);
    if (it == data.end())
        return nullptr;
    return &(*it);
}