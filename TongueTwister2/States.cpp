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
    
    // read the information on the state names, if available
    if (jd.hasKey("SegmentIPA") == true)
        segmentIPA = jd.getValue<std::vector<std::string>>("SegmentIPA");
    if (jd.hasKey("SegmentName") == true)
        stateNames = jd.getValue<std::vector<std::string>>("SegmentName");

    // add the partitioning of the state set into natural classes
    statePartitions = nullptr;
    if (jd.hasKey("PartitionSets") == true)
        statePartitions = new Partition(*jd.getJsonPointer("PartitionSets"));
}

void States::print(void) {

    std::cout << "   State Information" << std::endl;
    std::cout << "   * Number of states: " << numStates << std::endl;
    if (stateNames.size() == numStates || segmentIPA.size() == numStates)
        {
        std::cout << "   * State labels:" << std::endl;
        for (size_t i=0; i<numStates; i++)
            {
            std::cout << "     " << i << " : ";
            if (stateNames.size() == numStates)
                std::cout << stateNames[i] << " ";
            if (segmentIPA.size() == numStates == numStates)
                std::cout << segmentIPA[i] << " ";
            std::cout << std::endl;
            }
        }
    if (statePartitions != nullptr)
        statePartitions->print();
    std::cout << std::endl;
}
