#include "Msg.hpp"
#include "Parameter.hpp"


Parameter::Parameter(Model* m, RandomVariable* r, std::string n) : modelPtr(m), rng(r), name(n) {

}
