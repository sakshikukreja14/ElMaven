#ifndef CHARGESSPECIES_H
#define CHARGESSPECIES_H
#include "scan.h"

/**
 * @class ChargedSpecies
 * @ingroup libmaven
 * @brief Wrapper class for a charged species.
 * @author Elucidata
 */
class ChargedSpecies
{
    public:
    ChargedSpecies()
    {
        deconvolutedMass = 0;
        minZ = 0;
        maxZ = 0;
        countMatches = 0;
        error = 0;
        upCount = 0;
        downCount = 0;
        scan = NULL;
        totalIntensity = 0;
        isotopicClusterId = 0;
        minRt = maxRt = meanRt = 0;
        istotopicParentFlag = false;
        minIsotopeMass = 0;
        maxIsotopeMass = 0;
        massDiffectCluster = -100;
        filterOutFlag = false;
        qscore = 0;
    }

    /**
     * parent mass guess
     */
    float deconvolutedMass;
    float error;
    int minZ;
    int maxZ;
    int countMatches;
    float totalIntensity;
    float qscore;
    int upCount;
    int downCount;
    Scan *scan;

    vector<float> observedMzs;
    vector<float> observedIntensities;
    vector<float> observedCharges;

    int isotopicClusterId;
    int massDiffectCluster;
    bool istotopicParentFlag;
    bool filterOutFlag;

    float meanRt;
    float minRt;
    float maxRt;

    vector<ChargedSpecies *> isotopes;
    int minIsotopeMass;
    int maxIsotopeMass;

    /**
         * [compIntensity ]
         * @method compIntensity
         * @param  a             []
         * @param  b             []
         * @return []
         */
    static bool compIntensity(ChargedSpecies *a, ChargedSpecies *b) { return a->totalIntensity > b->totalIntensity; }

    /**
         * [compMass ]
         * @method compMass
         * @param  a        []
         * @param  b        []
         * @return []
         */
    static bool compMass(ChargedSpecies *a, ChargedSpecies *b) { return a->deconvolutedMass < b->deconvolutedMass; }

    /**
         * [compMatches ]
         * @method compMatches
         * @param  a           []
         * @param  b           []
         * @return []
         */
    static bool compMatches(ChargedSpecies *a, ChargedSpecies *b) { return a->countMatches > b->countMatches; }

    /**
         * [compRt ]
         * @method compRt
         * @param  a      []
         * @param  b      []
         * @return []
         */
    static bool compRt(ChargedSpecies *a, ChargedSpecies *b) { return a->meanRt < b->meanRt; }

    /**
         * [compMetaGroup ]
         * @method compMetaGroup
         * @param  a             []
         * @param  b             []
         * @return []
         */
    static bool compMetaGroup(ChargedSpecies *a, ChargedSpecies *b)
    {
        return (a->isotopicClusterId * 100000 + a->deconvolutedMass < b->isotopicClusterId * 100000 + b->deconvolutedMass);
    }

    /**
         * [updateIsotopeStatistics ]
         * @method updateIsotopeStatistics
         */
    void updateIsotopeStatistics()
    {
        minIsotopeMass = maxIsotopeMass = deconvolutedMass;
        for (unsigned int i = 0; i < isotopes.size(); i++)
        {
            if (isotopes[i]->deconvolutedMass < minIsotopeMass)
                minIsotopeMass = isotopes[i]->deconvolutedMass;
            if (isotopes[i]->deconvolutedMass > maxIsotopeMass)
                maxIsotopeMass = isotopes[i]->deconvolutedMass;
        }
    }
};
#endif // CHARGESSPECIES_H
