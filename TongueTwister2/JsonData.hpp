#ifndef JsonData_hpp
#define JsonData_hpp

#include <string>
#include <optional>
#include "json.hpp"



class JsonData {
    
    public:
        static JsonData&    jsonInstance(void)
                                {
                                static JsonData jsonData;
                                return jsonData;
                                }
        bool                hasKey(const std::string& key);                          // check if key exists anywhere in hierarchy (recursive)
        template <typename T>
        T                   getValue(const std::string& key);                        // get value with automatic type conversion (throws if not found)        
        template <typename T> 
        T                   getValue(const std::string& key, const T& defaultValue); // get value with default if not found
        template <typename T>
        std::optional<T>    getOptionalValue(const std::string& key);                // get optional value (returns std::nullopt if not found)
        nlohmann::json*     getJsonPointer(const std::string& key);                  // get pointer to json object at key (nullptr if not found)
        
                            // utility functions
        std::string         dump(void) { return j.dump(); }
        size_t              numKeys(void);
        void                readJsonFile(std::string fn);
        void                print(void);
    
    private:
                            JsonData(void);
                            JsonData(const JsonData&);
        JsonData&           operator=(const JsonData&);
        nlohmann::json*     findValueForKey(nlohmann::json& obj, const std::string& key);
        nlohmann::json      j;
};

template <typename T>
T JsonData::getValue(const std::string& key) {

    nlohmann::json* p = findValueForKey(j, key);
    if (p == nullptr)
        throw std::runtime_error("Key not found: " + key);
    return p->get<T>();
}

template <typename T>
T JsonData::getValue(const std::string& key, const T& defaultValue) {

    nlohmann::json* p = findValueForKey(j, key);
    if (p == nullptr)
        return defaultValue;
    return p->get<T>();
}

template <typename T>
std::optional<T> JsonData::getOptionalValue(const std::string& key){

    nlohmann::json* p = findValueForKey(j, key);
    if (p == nullptr)
        return std::nullopt;
    return p->get<T>();
}

#endif
