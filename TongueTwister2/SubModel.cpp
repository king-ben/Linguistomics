#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include "Ctmc.hpp"
#include "CtmcAccelerated.hpp"
#include "FamilyData.hpp"
#include "LikelihoodCalculator.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterExchangeabilities.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "States.hpp"
#include "SubModel.hpp"
#include "Threads.hpp"
#include "TransitionProbabilitiesCpu.hpp"
#include "TransitionProbabilitiesGpu.hpp"
#include "UserSettings.hpp"

SubModel::SubModel(Model* parentModel,
                   RandomVariable* rng,
                   FamilyData* familyData,
                   int familyIndex,
                   ThreadPool* pool)
    : parentModel(parentModel), rng(rng), familyData(familyData),
      familyIndex(familyIndex), pool(pool), numStates(0),
      tree(nullptr), indelRates(nullptr), frequencies(nullptr),
      exchangeabilities(nullptr), rateMatrix(nullptr), tiProbs(nullptr),
      ownsIndelRates(false), ownsFrequencies(false),
      ownsExchangeabilities(false), ownsRateMatrix(false) {
}

SubModel::~SubModel(void) {

    // Do NOT delete alignments, tree, indelRates, frequencies,
    // exchangeabilities, or calculators here.
    //
    // They are registered with Model and deleted by Model.

    if (tiProbs != nullptr)
        {
        delete tiProbs;
        tiProbs = nullptr;
        }

    // Only delete rateMatrix here if this SubModel created and owns it.
    // In the linked model, the shared rate matrix is owned by Model.
    if (ownsRateMatrix && rateMatrix != nullptr)
        {
        delete rateMatrix;
        rateMatrix = nullptr;
        }
}

bool SubModel::initializeAlignments(ParameterIndelRates* sharedIndelRates) {

    UserSettings& settings = UserSettings::userSettings();
    bool onlyComplete = settings.getUseOnlyCompleteWords();
    
    if (familyData->hasKey("Taxa") == false)
        Msg::error("Family " + std::to_string(familyIndex) + ": could not find \"Taxa\"");
    if (familyData->hasKey("Words") == false)
        Msg::error("Family " + std::to_string(familyIndex) + ": could not find \"Words\"");
    if (familyData->hasKey("NumberOfStates") == false)
        Msg::error("Family " + std::to_string(familyIndex) + ": could not find \"NumberOfStates\"");

    numStates = familyData->getValue<size_t>("NumberOfStates");
    
    std::vector<std::string> canonicalTaxonList =
        familyData->getValue<std::vector<std::string>>("Taxa");
    
    if (canonicalTaxonList.size() > sizeof(unsigned) * 8)
        Msg::error("Too many taxa in family " + std::to_string(familyIndex));

    // Set up indel rates for this family.
    // Linked model: use the shared object passed in.
    // Unlinked model: create a family-specific object.
    if (sharedIndelRates != nullptr)
        {
        indelRates = sharedIndelRates;
        }
    else
        {
        indelRates = new ParameterIndelRates(
            parentModel,
            rng,
            "Insertion/Deletion rates family " + std::to_string(familyIndex + 1),
            2.0,
            2.0
        );

        parentModel->registerParameter(indelRates);
        }

    const nlohmann::json* jWords = familyData->getJsonPointer("Words");
    size_t numWords = jWords->size();
    
    for (size_t i=0; i<numWords; i++)
        {
        nlohmann::json jw = (*jWords)[i];
        std::string alnName = jw["Name"];
        ParameterAlignment* aln = new ParameterAlignment(parentModel, rng, alnName);
        aln->initializeFromJson(jw, numStates, canonicalTaxonList);
        
        if (aln->getNumTaxa() != canonicalTaxonList.size() && onlyComplete)
            {
            delete aln;
            }
        else
            {
            // register with parent model so UpdateManager sees it
            parentModel->registerParameter(aln);
            alignments.push_back(aln);
            }
        }
    
    // find unique taxon combinations
    for (ParameterAlignment* p : alignments)
        {
        const unsigned& mask = p->getTaxonMask();
        uniqueTaxonCombinations.insert(mask);
        }
    
    std::cout << "   Family " << familyIndex + 1 << " alignments: "
              << alignments.size() << " cognates" << std::endl;
    
    return true;
}

bool SubModel::initializeTree(void) {

    if (familyData->hasKey("Tree") == false)
        Msg::error("Could not find \"Tree\" key in family "
                   + std::to_string(familyIndex + 1) + " JSON object");

    if (familyData->hasKey("Taxa") == false)
        Msg::error("Could not find \"Taxa\" key in family "
                   + std::to_string(familyIndex + 1) + " JSON object");

    std::vector<std::string> taxa =
        familyData->getValue<std::vector<std::string>>("Taxa");

    std::string newickString =
        familyData->getValue<std::string>("Tree");

    // default, matching ParameterTree constructor default
    double lambda = 10.0;

    // accept several possible names, since your JSONs currently use BrlenPriorVal
    if (familyData->hasKey("BranchLengthLambda") == true)
        lambda = familyData->getValue<double>("BranchLengthLambda");
    else if (familyData->hasKey("BrlenPriorVal") == true)
        lambda = familyData->getValue<double>("BrlenPriorVal");
    else if (familyData->hasKey("TreeLengthPriorVal") == true)
        lambda = familyData->getValue<double>("TreeLengthPriorVal");

    tree = new ParameterTree(parentModel,
                             rng,
                             "Tree_" + std::to_string(familyIndex + 1));

    tree->initialize(uniqueTaxonCombinations,
                     taxa,
                     newickString,
                     lambda);

    parentModel->registerParameter(tree);

    std::cout << "   Family " << familyIndex + 1 << " tree: "
              << tree->getNumTaxa() << " taxa" << std::endl;

    return true;
}

bool SubModel::initializeSubstitutionModel(
        ParameterFrequencies*       sharedFreqs,
        ParameterExchangeabilities* sharedExch,
        RateMatrix*                 sharedRateMatrix,
        SubstitutionModel           modelType) {

    UserSettings& settings = UserSettings::userSettings();
    MatrixExpBackend matExpBackend = settings.getMatrixExpBackend();

    Ctmc* subModel = nullptr;
    rateMatrix = nullptr;
    
    if (sharedFreqs != nullptr)
        {
        // linked model: use shared parameters, don't own them
        frequencies       = sharedFreqs;
        exchangeabilities = sharedExch;
        rateMatrix        = sharedRateMatrix;
        ownsFrequencies       = false;
        ownsExchangeabilities = false;
        ownsRateMatrix        = false;
        
        if (modelType == jc69)
            subModel = new CtmcJC69(numStates);
        else if (modelType == f81)
            subModel = new CtmcF81(numStates, frequencies);
        else if (modelType == gtr || modelType == custom)
            subModel = new CtmcAccelerated(numStates, rateMatrix);
        else
            Msg::error("Unknown model type");
        }
    else
        {
        // unlinked model: create own parameters
        ownsFrequencies       = true;
        ownsExchangeabilities = true;
        ownsRateMatrix        = true;
        
        if (modelType == jc69)
            {
            subModel = new CtmcJC69(numStates);
            ownsFrequencies       = false;
            ownsExchangeabilities = false;
            ownsRateMatrix        = false;
            }
        else if (modelType == f81)
            {
            frequencies = new ParameterFrequencies(parentModel, rng,
                "Frequencies_" + std::to_string(familyIndex), numStates);
            parentModel->registerParameter(frequencies);
            subModel = new CtmcF81(numStates, frequencies);
            ownsExchangeabilities = false;
            ownsRateMatrix        = false;
            }
        else if (modelType == gtr)
            {
            frequencies = new ParameterFrequencies(parentModel, rng,
                "Frequencies_" + std::to_string(familyIndex), numStates);
            exchangeabilities = new ParameterExchangeabilities(parentModel, rng,
                "Exchangeabilities_" + std::to_string(familyIndex), numStates);
            parentModel->registerParameter(frequencies);
            parentModel->registerParameter(exchangeabilities);
            rateMatrix = new RateMatrixGTR(numStates, exchangeabilities, frequencies);
            rateMatrix->updateRateMatrix();
            subModel = new CtmcAccelerated(numStates, rateMatrix);
            }
        else if (modelType == custom)
            {
            States* states = parentModel->getStates();
            Partition* part = states->getPartition();
            frequencies = new ParameterFrequencies(parentModel, rng,
                "Frequencies_" + std::to_string(familyIndex), numStates);
            exchangeabilities = new ParameterExchangeabilities(parentModel, rng,
                "Exchangeabilities_" + std::to_string(familyIndex),
                numStates, part->numSubsets());
            parentModel->registerParameter(frequencies);
            parentModel->registerParameter(exchangeabilities);
            rateMatrix = new RateMatrixCustom(numStates, exchangeabilities,
                                              frequencies, part);
            rateMatrix->updateRateMatrix();
            subModel = new CtmcAccelerated(numStates, rateMatrix);
            }
        else
            Msg::error("Unknown model type");
        }
    
    if (subModel == nullptr)
        Msg::error("Failed to create substitution model for family "
                   + std::to_string(familyIndex));

    // create transition probability manager for this family's tree
    switch (matExpBackend)
        {
        case cpuThreaded:
            tiProbs = new TransitionProbabilitiesCpu(pool, tree, subModel);
            break;
        case cpuBatched:
            {
            TransitionProbabilitiesGpu* g = new TransitionProbabilitiesGpu(pool, tree, subModel);
            g->setComputeBackend(ComputeBackend::CPUBatched);
            tiProbs = g;
            }
            break;
        case gpuBatched:
            {
            TransitionProbabilitiesGpu* g = new TransitionProbabilitiesGpu(pool, tree, subModel);
            g->setComputeBackend(ComputeBackend::GPUBatched);
            tiProbs = g;
            }
            break;
        case autoBackend:
        default:
            {
            TransitionProbabilitiesGpu* g = new TransitionProbabilitiesGpu(pool, tree, subModel);
            g->setComputeBackend(ComputeBackend::Auto);
            tiProbs = g;
            }
            break;
        }
    
    tiProbs->updateAllBranches();
    return true;
}

bool SubModel::initializeCalculators(ParameterIndelRates* indelRates,
                                     ParameterFrequencies* freqs) {

    size_t cognateBase = parentModel->getNumCalculators();
    
    for (size_t i=0; i<alignments.size(); i++)
        {
        ParameterAlignment* aln = alignments[i];
        LikelihoodCalculator* calc = new LikelihoodCalculator(
            tiProbs, aln, tree, indelRates, freqs);
        calculators.push_back(calc);
        parentModel->registerCalculator(calc);
        
        aln->setCognateIndex(cognateBase + i);
        }
    
    return true;
}