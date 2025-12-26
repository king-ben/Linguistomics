#ifndef Parameter_hpp
#define Parameter_hpp

#include <string>
class Model;
class RandomVariable;



class Parameter {

    public:
                                    Parameter(void) = delete;
                                    Parameter(Model* m, RandomVariable* r, std::string n);
        virtual                    ~Parameter(void) { }
        std::string&                getName(void) { return name; }
        virtual void                keep(void) = 0;
        virtual double              lnPriorProbability(void) = 0;
        virtual void                print(void) = 0;
        virtual void                restore(void) = 0;
        
    protected:
        Model*                      modelPtr;
        RandomVariable*             rng;
        std::string                 name;
};

#endif
