#include <algorithm>
#include <iostream>
#include "Msg.hpp"
#include "ParameterStatistics.hpp"


ParameterStatistics::ParameterStatistics(void) {

    name = "";
    values.clear();
}

ParameterStatistics::ParameterStatistics(const ParameterStatistics& ps) {

    name = ps.name;
    values = ps.values;
}

double ParameterStatistics::operator[](size_t idx) {

    if (idx >= values.size())
        Msg::error("Subscript index out-of-range");
    return values[idx];
}

double ParameterStatistics::operator[](size_t idx) const {

    if (idx >= values.size())
        Msg::error("Subscript index out-of-range");
    return values[idx];
}

ParameterStatistics& ParameterStatistics::operator+=(const ParameterStatistics& rhs) {

    for (int i=0; i<rhs.size(); i++)
        {
        double x = rhs[i];
        this->addValue(x);
        }
    return *this;
}

nlohmann::json ParameterStatistics::toJson(void) {

    return toJson(0.0);
}

nlohmann::json ParameterStatistics::toJson(float burnFraction) {

    if (burnFraction < 0.0 || burnFraction > 0.95)
        Msg::error("Burn fraction must be between 0 and 0.95");
    if (values.size() * (1.0 - burnFraction) < 3)
        Msg::error("Can only calclate statistics for at least three samples");

    double m = Statistics::getMeanAndVariance(values, burnFraction).first;
    CredibleInterval i = Statistics::getCredibleInterval(values, burnFraction);
    nlohmann::json j = nlohmann::json::object();
    j["cognate"] = name;
    j["mean"] = m;
    j["lower"] = i.lower;
    j["upper"] = i.upper;
    return j;
}

void ParameterStatistics::toFile(std::ostream& file) {

    double m = Statistics::getMeanAndVariance(values).first;
    CredibleInterval i = Statistics::getCredibleInterval(values);
    file << "new(" << i.lower << "," << m << "," << i.upper << ")";
}
