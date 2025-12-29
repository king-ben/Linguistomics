#include "Parameter.hpp"



Parameter::Parameter(RandomVariable* r, Model* m, std::string n) {

    rv = r;
    modelPtr = m;
    parmName = n;
}
