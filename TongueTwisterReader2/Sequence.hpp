#ifndef Sequence_hpp
#define Sequence_hpp

#include <cstdint>
#include <string>
#include <vector>



class Sequence {

    public:
                                        Sequence(void) = delete;
                                        Sequence(std::string n, std::vector<int8_t> d);
        int8_t&                         operator[](size_t i) { return data[i]; }
        const int8_t&                   operator[](size_t i) const { return data[i]; }
        const std::string&              getName(void) const { return name; }
        const std::vector<int8_t>&      getSequence(void) const { return data; }
        size_t                          size(void) const { return data.size(); }
        
    private:
        std::string                     name;
        std::vector<int8_t>             data;
};

#endif
