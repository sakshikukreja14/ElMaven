#ifndef PEAKGROUP_H
#define PEAKGROUP_H

#include "Fragment.h"
#include "Peak.h"
#include "datastructures/mzSlice.h"
#include "standardincludes.h"

class mzSample;
class Isotope;
class MassCalculator;
class Compound;
class Peak;
class Scan;
class EIC;
class MassCutoff;
class Adduct;

using namespace std;

class PeakGroup
{
    private:
        Adduct* _adduct;
        mzSlice _slice;
        bool _sliceSet;
        string _tableName;

    public:
        enum class GroupType {
            None = 0,
            C13 = 1,
            Adduct = 2,
            Covariant = 4,
            Isotope = 5
        };
        enum QType {
            AreaTop = 0,
            Area = 1,
            Height = 2,
            AreaNotCorrected = 3,
            RetentionTime = 4,
            Quality = 5,
            SNRatio = 6,
            AreaTopNotCorrected = 7
        };
        /**
         * @brief PeakGroup Default Constructor.
         */
        PeakGroup();
        /**
         * @brief PeakGroup Copy Constructor.
         * @param o
         */
        PeakGroup(const PeakGroup& o);
        /**
         * @brief operator =    Assignment operator overloading
         * @param o     PeakGroup object.
         * @return
         */
        PeakGroup& operator = (const PeakGroup& o);
        /**
         * @brief operator ==   Equals to operator overloading.
         * @param o     PeakGroup object.
         * @return
         */
        bool operator==(const PeakGroup* o);
        /**
         * @brief copyObj   Copies of object one PeakGroup
         * to another.
         * @param o
         */
        void copyObj(const PeakGroup& o);
        /**
         * @brief Destructor of class PeakGroup.
         */
        ~PeakGroup();

        PeakGroup* parent;
        /**
         * @brief This parent group represents a group for
         * the parent (primary) ion, from which this adduct
         * group would have formed.
         */
        PeakGroup* parentIon;
        vector<Peak> peaks;
        vector<PeakGroup> children;
        vector<PeakGroup> childAdducts;
        vector<PeakGroup> childrenBarPlot;
        vector<PeakGroup> childrenIsoWidget;
        /**
         * @brief samples Sample which has been used
         * for peak detection.
         */
        vector<mzSample*> samples;
        string srmId;
        string tagString;
        /**
         * classification label
         */
        char label;

        bool isFocused;

        QType quantitationType;
        GroupType _type;

        int groupId;
        int metaGroupId;
        int clusterId;

        bool deletedFlag;

        float maxIntensity;
        float maxAreaTopIntensity;
        float maxAreaIntensity;
        float maxHeightIntensity;
        float maxAreaNotCorrectedIntensity;
        float maxAreaTopNotCorrectedIntensity;
        float currentIntensity;
        float meanRt;
        float meanMz;
        float expectedMz;
        int totalSampleCount;

        int ms2EventCount;

        FragmentationMatchScore fragMatchScore;
        Fragment fragmentationPattern;

        // isotopic information
        float expectedAbundance;
        int isotopeC13count;

        double minIntensity;

        float minRt;
        float maxRt;
        float minMz;
        float maxMz;

        float blankMax;
        float blankMean;
        unsigned int blankSampleCount;

        int sampleCount;
        float sampleMean;
        float sampleMax;

        unsigned int maxNoNoiseObs;
        unsigned int maxPeakOverlap;
        float maxQuality;
        float avgPeakQuality;
        //@Kailash: group quality computed using neural network
        int markedGoodByCloudModel = 0;
        int markedBadByCloudModel = 0;
        float groupQuality;
        float weightedAvgPeakQuality;
        int predictedLabel;
        double minQuality;
        float maxPeakFracionalArea;
        float maxSignalBaseRatio;
        float maxSignalBaselineRatio;
        int goodPeakCount;
        float groupRank;

        // for sample contrasts  ratio and pvalue
        float changeFoldRatio;
        float changePValue;
        /**
         * @brief returns group name
         * @method getName.
         * @detail returns compound name, tagString, srmID,
         * meanMz@meanRt or groupId in this order of preference.
         * @return string
         */
        string getName();
        /**
         * @brief isMS1 Checks if scan level is 1.
         * @return  bool value.
         */
        bool isMS1();
        /**
         * @brief hasSrmId  Checks if peakgroup has its SrmId set.
         * @return
         */
        bool hasSrmId() const
        {
            return srmId.empty();
        }
        /**
         * @brief hasCompoundLink   Checks if compound is linked
         * with the slice.
         * @return  bool value.
         */
        bool hasCompoundLink() const;

        /**
         * @brief isEmpty   Checks if 'peaks' vector is empty.
         * @return
         */
        inline bool isEmpty() const
        {
            if (peaks.size() == 0)
                return true;
            return false;
        }
        /**
         * @brief peakCount Returns the peaks vector size.
         * @return
         */
        inline unsigned int peakCount() const
        {
            return peaks.size();
        }
        /**
         * @brief childCount    Returns the children count of peakgroup.
         * @return
         */
        inline unsigned int childCount() const
        {
            return children.size();
        }
        /**
         * @brief childCountBarPlot Returns childBarPlot vector size.
         * @return
         */
        inline unsigned int childCountBarPlot() const
        {
            return childrenBarPlot.size();
        }
        /**
         * @brief childCountIsoWidget   Returns childIsoWidget vector size.
         * @return
         */
        inline unsigned int childCountIsoWidget() const
        {
            return childrenIsoWidget.size();
        }
        /**
         * @brief Check whether a slice has previosuly been set for this group.
         * @return true if a slice has been set, false otherwise.
         */
        bool hasSlice() const;
        /**
         * @brief Check whether both bounds of the group's slice are close to
         * zero in either m/z or rt dimensions.
         * @return `true` if either slice's `mzmin` and `mzmax` are both close
         * to zero or slice's `rtmin` and `rtmax` are both close to zero, false
         * otherwise.
         */
        bool sliceIsZero() const;
        /**
         * @brief getCompound   Return associated slice's compound.
         * @return
         */
        Compound* getCompound() const;
        /**
         * @brief setCompound   Sets the compound for the slice.
         * @param compound
         */
        void setCompound(Compound* compound);

        /**
         * @brief setSlice  Sets the value for the slice of peakgroup.
         * @param slice
         */
        void setSlice(const mzSlice& slice);
        /**
         * @brief getSlice Returns slice.
         * @return
         */
        const mzSlice& getSlice() const;

        /**
         * @brief getParent Returns parent for the PeakGroup.
         * @return
         */
        inline PeakGroup* getParent()
        {
            return parent;
        }
        /**
         * @brief getPeaks Returns peak vector.
         * @return
         */
        inline vector<Peak>& getPeaks()
        {
            return peaks;
        }
        /**
         * @brief getChildren   Returns children vecotr of the peakgroup.
         * @return
         */
        inline vector<PeakGroup>& getChildren()
        {
            return children;
        }

        vector<Scan*> getRepresentativeFullScans();
        /**
         * @brief find all MS2 scans for this group
         * @return vector of all MS2 scans for this group
         */
        vector<Scan*> getFragmentationEvents();
        /**
         * @brief build a consensus fragment spectra for this group
         * @param productPpmTolr ppm tolerance for fragment m/z
         */
        void computeFragPattern(float productPpmTolr);

        Scan* getAverageFragmentationScan(float productPpmTolr);

        void matchFragmentation(float ppmTolerance, string scoringAlgo);

        double getExpectedMz(int charge);
        /**
         * @brief setParent Sets the parent of the peakgroup and
         * assigns its group type to 'Isotope'.
         * @param p
         */
        inline void setParent(PeakGroup* p)
        {
            parent = p;
            if (parent != nullptr)
                _type = GroupType::Isotope;
        }
        /**
         * @brief setLabel  Sets label value for the peakgroup.
         * @param label
         */
        void setLabel(char label);
        /**
         * @brief Set the adduct form for this `PeakGroup`.
         * @details If the adduct is a type of parent ion, then this group's
         * `_type` attribute is set to `GroupType::Adduct`.
         * @param adduct Pointer to an `Adduct` object to be assigned.
         */
        void setAdduct(Adduct* adduct);
        /**
         * @brief Get the adduct form for this `PeakGroup`.
         * @return Pointer to the `Adduct` object set for this group.
         */
        Adduct* getAdduct() const;
        /**
         * [ppmDist ]
         * @method ppmDist
         * @param  cmass   []
         * @return []
         */
        float massCutoffDist(float cmass, MassCutoff* massCutoff);
        /**
         * @brief addPeak   Adds the peak in peaks vector.
         * @param peak
         */
        void addPeak(const Peak& peak);
        /**
         * @brief addChild  Adds child to the children vector.
         * @param child
         */
        inline void addChild(const PeakGroup& child)
        {
            children.push_back(child);
            children.back().parent = this;
        }
        /**
         * @brief addChildBarPlot   Adds childBarPlot to its vector.
         * @param child
         */
        inline void addChildBarPlot(const PeakGroup& child)
        {
            childrenBarPlot.push_back(child);
            childrenBarPlot.back().parent = this;
        }
        /**
         * @brief addChildIsoWidget Adds childIsoWidget to its vector.
         * @param child
         */
        inline void addChildIsoWidget(const PeakGroup& child)
        {
            childrenIsoWidget.push_back(child);
            childrenIsoWidget.back().parent = this;
        }
        /**
         * @brief setGroupIdForChildren Sets groupId for children.
         */
        inline void setGroupIdForChildren()
        {
            for (auto& child : children)
                child.groupId = groupId;
        }

        /**
         * @brief getPeak   Return the peak of particular sample
         * from PeakGroup.
         * @param sample
         * @return
         */
        Peak* getPeak(mzSample* sample);
        /**
         * @brief type  Returns the group type of the PeakGroup.
         * @return
         */
        inline GroupType type() const
        {
            return _type;
        }
        /**
         * @brief setType   Sets the group type for the peakGroup.
         * @param t
         */
        inline void setType(GroupType t)
        {
            _type = t;
        }
        /**
         * @brief setQuantitationType   Sets the Quantitation type
         * for the peakgroup.
         * @param type
         */
        void setQuantitationType(QType type)
        {
            quantitationType = type;
        }
        /**
         * @brief isIsotope Checks if the peakgroup has its group id
         * as Isotope.
         * @return
         */
        inline bool isIsotope() const
        {
            return _type == GroupType::Isotope;
        }
        /**
         * @brief isAdduct  Checks if the peakgroup type is Adduct.
         * @return
         */
        inline bool isAdduct() const
        {
            return _type == GroupType::Adduct;
        }
        /**
         * @brief summary   Prints the summary of the PeakGroup.
         */
        void summary();
        /**
         * @brief groupStatistics Calculates various statistical value for
         * data members of a peakgroup.
         */
        void groupStatistics();
        /**
         * @brief updateQuality Updates the quality of the peakgroup.
         */
        void updateQuality();
        /**
         * @brief medianRt  Calculates the median of Rts.
         * @return
         */
        float medianRt();
        /**
         * @brief Obtain the deviation of this peak-group from its expected
         * retention time.
         * @details The RT difference between expected and observed is returned
         * as an absolute value. In case there was no expected value for this
         * group, then the value returned is -1.0.
         * @return A floating point value denoting absolute RT deviation in
         * minutes.
         */
        float expectedRtDiff();
        /**
         * @brief reduce Make sure there is only one peak per sample.
         */
        void reduce();
        /**
         * @brief fillInPeaks   Adds the missing peaks in the peak vector according
         * to eics vector.
         * @param eics
         */
        void fillInPeaks(const vector<EIC*>& eics);
        /**
         * [computeAvgBlankArea ]
         * @method computeAvgBlankArea
         * @param  eics                []
         */
        void computeAvgBlankArea(const vector<EIC*>& eics);
        /**
         * @brief groupOverlapMatrix    Checks overlap of the peaks in peakgroup
         * and updates the value for groupOverlapFrac
         */
        void groupOverlapMatrix();
        /**
         * @brief deletePeaks   Clears the value of peaks in peakgroup.
         */
        void deletePeaks();
        /**
         * @brief clear Sets the values of pointers to Null.
         */
        void clear();
        /**
         * [deleteChildren ]
         * @method deleteChildren
         */
        void deleteChildren();
        /**
         * [deleteChild ]
         * @method deleteChild
         * @param  child       []
         * @return []
         */
        bool deleteChild(PeakGroup* child);
        /**
         * [copyChildren ]
         * @method copyChildren
         * @param  other        []
         */
        void copyChildren(const PeakGroup& other);
        /**
         * @brief getOrderedIntensityVector Return intensity vectory
         * ordered by samples.
         * @param samples   mzSample object.
         * @param type      Group type.
         * @return
         */
        vector<float> getOrderedIntensityVector(vector<mzSample*>& samples,
                                                QType type);
        /**
         * @brief reorderSamples    Sorts the samples according to peak intensity.
         */
        void reorderSamples();
        /**
         * @brief compRt    Compares mean Rt of the two peakgroups.
         * @param a         PeakGroup object.
         * @param b         PeakGroup object.
         * @return          Boolean value.
         */
        static bool compRt(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.meanRt < b.meanRt);
        }
        /**
         * @brief compMz    Compares mean Mz of the two peakgroups.
         * @param a         PeakGroup object.
         * @param b         PeakGroup object.
         * @return          Boolean value.
         */
        static bool compMz(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.meanMz > b.meanMz);
        }
        /**
         * @brief compIntensity Compares max intensity of the two
         * peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compIntensity(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.maxIntensity > b.maxIntensity);
        }
        /**
         * @brief compArea  Compares max peak fractional area of two
         * peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compArea(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.maxPeakFracionalArea > b.maxPeakFracionalArea);
        }
        /**
         * @brief compQuality   Compares max quality of the two peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compQuality(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.maxQuality > b.maxQuality);
        }
        /**
         * @brief compRank  Compares group rank of the two peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compRank(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.groupRank < b.groupRank);
        }
        /**
         * @brief compRatio     Compares change fold ration of the two peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compRatio(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.changeFoldRatio < b.changeFoldRatio);
        }
        /**
         * @brief compPvalue    Compares p value of the two peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compPvalue(const PeakGroup* a, const PeakGroup* b)
        {
            return (a->changePValue < b->changePValue);
        }
        /**
         * @brief compC13   Compares isotopeC13count of the two peakgroups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compC13(const PeakGroup* a, const PeakGroup* b)
        {
            return (a->isotopeC13count < b->isotopeC13count);
        }
        /**
         * @brief compMetaGroup Compares metaGroupId for the groups.
         * @param a     PeakGroup object.
         * @param b     PeakGroup object.
         * @return      Boolean value.
         */
        static bool compMetaGroup(const PeakGroup& a, const PeakGroup& b)
        {
            return (a.metaGroupId < b.metaGroupId);
        }
        /**
         * @brief operator <    '<' operating overloading.
         * @param b
         * @return
         */
        bool operator<(const PeakGroup* b) const
        {
            return this->maxIntensity < b->maxIntensity;
        }

        void calGroupRank(bool deltaRtCheckFlag,
                          int qualityWeight,
                          int intensityWeight,
                          int deltaRTWeight);
        /**
         * @brief take list of sample and filter those which are marked as selected
         * @details this method used for assigning samples to this group based on
         * whether that samples are marked as selected.
         */
        void setSelectedSamples(vector<mzSample*> vsamples);

        /**
         * @brief Obtain the name of peak-table this group belongs to.
         * @return A string containing peak-table name.
         */
        string tableName() const;

        /**
         * @brief Set the peak-table (name) in which this group is listed.
         * @param tableName A string containing peak-table name.
         */
        void setTableName(string tableName);
};
#endif
