#ifndef ParameterStatistics_hpp
#define ParameterStatistics_hpp

#include "json.hpp"
#include "Statistics.hpp"
#include <string>
#include <vector>



class ParameterStatistics {

    public:
                                    ParameterStatistics(void);
                                    ParameterStatistics(const ParameterStatistics& ps);
        double                      operator[](size_t idx);
        double                      operator[](size_t idx) const;
        ParameterStatistics&        operator+=(const ParameterStatistics& rhs);
        void                        addValue(float x) { values.push_back(x); }
        std::string                 getName(void) { return name; }
        const std::vector<float>&   getValues(void) const { return values; }
        int                         size(void) { return (int)values.size(); }
        int                         size(void) const { return (int)values.size(); }
        void                        setName(std::string s) { name = s; }
        nlohmann::json              toJson(void);
        nlohmann::json              toJson(float burnFraction);
        void                        toFile(std::ostream& findex);

    private:
        std::string                 name;
        std::vector<float>          values;
};

#endif
