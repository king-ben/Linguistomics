#include "JsonData.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "States.hpp"



States::States(void) {

    initialzieStates();
}

States::~States(void) {

    if (statePartitions != nullptr)
        delete statePartitions;
}

void States::initialzieStates(void) {

    JsonData& jd = JsonData::jsonInstance();
 
    // read the number of states
    if (jd.hasKey("NumberOfStates") == false)
        Msg::error("Could not find the number of states in JSON object");
    numStates = jd.getValue<int>("NumberOfStates");

    // add the partitioning of the state set into natural classes
    statePartitions = nullptr;
    if (jd.hasKey("PartitionSets") == true)
        statePartitions = new Partition(*jd.getJsonPointer("PartitionSets"));
}

void States::print(void) {

    std::cout << "   State Information" << std::endl;
    std::cout << "   * Number of states: " << numStates << std::endl;
    if (statePartitions != nullptr)
        statePartitions->print();
}
