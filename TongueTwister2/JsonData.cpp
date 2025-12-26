#include <fstream>
#include <iostream>
#include "JsonData.hpp"
#include "Msg.hpp"



JsonData::JsonData(void) {

}

bool JsonData::hasKey(const std::string& key) {

    return findValueForKey(j, key) != nullptr;
}

nlohmann::json* JsonData::getJsonPointer(const std::string& key) {

    return findValueForKey(j, key);
}

nlohmann::json* JsonData::findValueForKey(nlohmann::json& obj, const std::string& key) {

    // If this is an object, check its keys
    if (obj.is_object())
        {
        // Does this object directly contain the key?
        if (obj.contains(key))
            {
            return &obj[key];
            }
        // Otherwise recurse into child values
        for (auto& [k, v] : obj.items())
            {
            nlohmann::json* result = findValueForKey(v, key);
            if (result != nullptr)
                return result;
            }
        }
    // If this is an array, search each element
    else if (obj.is_array())
        {
        for (auto& element : obj)
            {
            nlohmann::json* result = findValueForKey(element, key);
            if (result != nullptr)
                return result;
            }
        }
    
    return nullptr;
}

size_t JsonData::numKeys(void) {

    return j.size();
}

void JsonData::print(void) {

    std::string jDump = j.dump(3);
    std::cout << jDump << std::endl;
}

void JsonData::readJsonFile(std::string fn) {

    std::cout << "   Input File Information" << std::endl;
    std::cout << "   * File name: " << fn << std::endl << std::endl;
    
    std::ifstream ifs(fn);
    if (!ifs)
        Msg::error("Could not find file \"" + fn + "\"");
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::error("Error parsing JSON file at byte " + std::to_string(ex.byte));
        }
}
