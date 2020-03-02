#ifndef PEAK_H
#define PEAK_H

#include "standardincludes.h"

class mzSample;
class EIC;
class Scan;
class mzLink;

using namespace std;

class Peak
{
    public:
    /**
     * @brief Peak Non-parameterised constructor.
     */
    Peak();

    /**
     * @brief Peak  Parameterised constructor.
     * @param e     Object of class EIC.
     * @param p
     */
    Peak(EIC* e, int p);

    /**
     * @brief Peak  Copy constructor.
     * @param p     Object of class peak.
     */
    Peak(const Peak& p);

    /**
     * @brief operator =  Overloading assignment
     * operator.
     * @param o  Object of class peak.
     * @return
     */
    Peak& operator = (const Peak& o);

    /**
     * @brief findCovariants Returns the Covariants.
     * @return
     */
    vector<mzLink> findCovariants();

    /**
     * @brief hasEIC Checks whether peak has value for
     * eic set.
     * @return bool value.
     */
    inline bool hasEIC() const
    {
        return (_eic != NULL);
    }

    /**
     * @brief hasSample Checks where the peak has value
     * for sample set.
     * @return bool value.
     */
    inline bool hasSample() const
    {
        return (_sample != NULL);
    }

    /**
     * @brief Getters and Setters for private dataMembers.
     */
    void setEIC(EIC* e)
    {
        _eic = e;
    }

    inline EIC* eic()
    {
        return _eic;
    }

    Scan* getScan();

    void setSample(mzSample* s)
    {
        _sample = s;
    }

    void setScan(int s)
    {
        this->_scan = s;
    }

    inline int scan() const
    {
        return this->_scan;
    }

    inline mzSample* sample()
    {
        return _sample;
    }

    void setLabel(char label)
    {
        this->_label = label;
    }

    inline char label() const
    {
        return _label;
    }

    void setPos(int pos)
    {
        this->_pos = pos;
    }

    inline int pos() const
    {
        return this->_pos;
    }

    void setMinpos(int minpos)
    {
        this->_minpos = minpos;
    }

    inline int minpos() const
    {
        return this->_minpos;
    }

    void setMaxpos(int maxpos)
    {
        this->_maxpos = maxpos;
    }

    inline int maxpos() const
    {
        return this->_maxpos;
    }

    void setSplineminpos(int min)
    {
        this->_splineminpos = min;
    }

    inline int splineminpos() const
    {
        return this->_splineminpos;
    }

    void setSplinemaxpos(int max)
    {
        this->_splinemaxpos = max;
    }

    inline int splinemaxpos() const
    {
        return this->_splinemaxpos;
    }

    void setRt(float rt)
    {
        this->_rt = rt;
    }

    inline float rt() const
    {
        return this->_rt;
    }

    void setRtmin(float min)
    {
        this->_rtmin = min;
    }

    inline float rtmin() const
    {
        return this->_rtmin;
    }

    void setRtmax(float max)
    {
        this->_rtmax = max;
    }

    inline float rtmax() const
    {
        return this->_rtmax;
    }

    void setMzmin(float min)
    {
        this->_mzmin = min;
    }

    inline float mzmin() const
    {
        return this->_mzmin;
    }

    void setMzmax(float max)
    {
        this->_mzmax = max;
    }

    inline float mzmax() const
    {
        return this->_mzmax;
    }

    void setMinscan(int min)
    {
        this->_minscan = min;
    }

    inline int minscan() const
    {
        return this->_minscan;
    }

    void setMaxscan(int max)
    {
        this->_maxscan = max;
    }

    inline int maxscan() const
    {
        return this->_maxscan;
    }

    void setPeakArea(float area)
    {
        this->_peakArea = area;
    }

    inline float peakArea() const
    {
        return this->_peakArea;
    }

    void setPeakSplineArea(float area)
    {
        this->_peakSplineArea = area;
    }

    inline float peakSplineArea() const
    {
        return this->_peakSplineArea;
    }

    void setPeakAreaCorrected(float corrected)
    {
        this->_peakAreaCorrected = corrected;
    }

    inline float peakAreaCorrected() const
    {
        return this->_peakAreaCorrected;
    }

    void setPeakAreaTop(float top)
    {
        this->_peakAreaTop = top;
    }

    inline float peakAreaTop() const
    {
        return this->_peakAreaTop;
    }

    void setPeakAreaTopCorrected(float corrected)
    {
        this->_peakAreaTopCorrected = corrected;
    }

    inline float peakAreaTopCorrected() const
    {
        return this->_peakAreaTopCorrected;
    }

    void setPeakAreaFractional(float fractional)
    {
        this->_peakAreaFractional = fractional;
    }

    inline float peakAreaFractional() const
    {
        return this->_peakAreaFractional;
    }

    void setPeakRank(float rank)
    {
        this->_peakRank = rank;
    }

    inline float peakRank() const
    {
        return this->_peakRank;
    }

    void setPeakIntensity(float intensity)
    {
        this->_peakIntensity = intensity;
    }

    inline float peakIntensity() const
    {
        return this->_peakIntensity;
    }

    void setPeakBaseLineLevel(float level)
    {
        this->_peakBaseLineLevel = level;
    }

    inline float peakBaseLineLevel() const
    {
        return this->_peakBaseLineLevel;
    }

    void setPeakMz(float mz)
    {
        this->_peakMz = mz;
    }

    inline float peakMz() const
    {
        return this->_peakMz;
    }

    void setMedianMz(float mz)
    {
        this->_medianMz = mz;
    }

    inline float medianMz() const
    {
        return this->_medianMz;
    }

    void setBaseMz(float mz)
    {
        this->_baseMz = mz;
    }

    inline float baseMz() const
    {
        return this->_baseMz;
    }

    void setQuality(float quality)
    {
        this->_quality = quality;
    }

    inline float quality() const
    {
        return this->_quality;
    }

    void setGaussFitSigma(float sigma)
    {
        this->_gaussFitSigma = sigma;
    }

    inline float gaussFitSigma() const
    {
        return this->_gaussFitSigma;
    }

    void setGaussFitR2(float r2)
    {
        this->_gaussFitR2 = r2;
    }

    inline float gaussFitR2() const
    {
        return this->_gaussFitR2;
    }

    void setWidth(unsigned int width)
    {
        this->_width = width;
    }

    inline int width() const
    {
        return this->_width;
    }

    void setGroupNum(int num)
    {
        this->_groupNum = num;
    }

    inline int groupNum() const
    {
        return this->_groupNum;
    }

    void setNoNoiseObs(unsigned int num)
    {
        this->_noNoiseObs = num;
    }

    inline unsigned int noNoiseObs() const
    {
        return this->_noNoiseObs;
    }

    void setNoNoiseFraction(float noise)
    {
        this->_noNoiseFraction = noise;
    }

    inline float noNoiseFraction() const
    {
        return this->_noNoiseFraction;
    }

    void setSymmetry(float symmetry)
    {
        this->_symmetry = symmetry;
    }

    inline float symmetry() const
    {
        return this->_symmetry;
    }

    void setSignalBaselineRatio(float ratio)
    {
        this->_signalBaselineRatio = ratio;
    }

    inline float signalBaselineRatio() const
    {
        return this->_signalBaselineRatio;
    }

    void setSignalBaselineDifference(float baselineDif)
    {
        this->_signalBaselineDifference = baselineDif;
    }

    inline float signalBaselineDifference() const
    {
        return this->_signalBaselineDifference;
    }

    void setGroupOverlap(float overlap)
    {
        this->_groupOverlap = overlap;
    }

    inline float groupOverlap() const
    {
        return this->_groupOverlap;
    }

    void setGroupOverlapFrac(float fraction)
    {
        this->_groupOverlapFrac = fraction;
    }

    inline float groupOverlapFrac() const
    {
        return this->_groupOverlapFrac;
    }

    void setLocalMaxFlag(bool maxFlag)
    {
        this->_localMaxFlag = maxFlag;
    }

    inline bool localMaxFlag() const
    {
        return this->_localMaxFlag;
    }

    void setFromBlankSample(bool sample)
    {
        this->_fromBlankSample = sample;
    }

    inline bool fromBlankSample() const
    {
        return this->_fromBlankSample;
    }

    static bool compRtMin(const Peak& a, const Peak& b)
    {
        return a._rtmin < b._rtmin;
    }
    /**
     * [compRt ]
     * @method compRt
     * @param  a      []
     * @param  b      []
     * @return []
     */
    static bool compRt(const Peak& a, const Peak& b)
    {
        return a._rt < b._rt;
    }

    /**
     * [compIntensity ]
     * @method compIntensity
     * @param  a             []
     * @param  b             []
     * @return []
     */
    static bool compIntensity(const Peak& a, const Peak& b)
    {
        return b._peakIntensity < a._peakIntensity;
    }

    /**
     * [compArea ]
     * @method compArea
     * @param  a        []
     * @param  b        []
     * @return []
     */
    static bool compArea(const Peak& a, const Peak& b);

    /**
     * [compMz ]
     * @method compMz
     * @param  a      []
     * @param  b      []
     * @return []
     */
    static bool compMz(const Peak& a, const Peak& b);

    /**
     * [compSampleName ]
     * @method compSampleName
     * @param  a              []
     * @param  b              []
     * @return []
     */
    static bool compSampleName(const Peak& a, const Peak& b);

    /**
     * [compSampleOrder ]
     * @method compSampleOrder
     * @param  a               []
     * @param  b               []
     * @return []
     */
    static bool compSampleOrder(const Peak& a, const Peak& b);

    /**
     * [overlap ]
     * @method overlap
     * @param  a       []
     * @param  b       []
     * @return []
     */
    static float overlap(const Peak& a, const Peak& b);

    private:
    /**
     * pointer to eic
     */
    EIC* _eic;

    /**
     * pointer to sample
     */
    mzSample* _sample;
    unsigned int _pos;
    unsigned int _minpos;        // left bound of peak
    unsigned int _maxpos;        // right bound of peak
    unsigned int _splineminpos;  // left bound of spline peak
    unsigned int _splinemaxpos;  // right bound of spline peak

    float _rt;
    float _rtmin;
    float _rtmax;
    float _mzmin;
    float _mzmax;

    unsigned int _scan;
    unsigned int _minscan;
    unsigned int _maxscan;

    /**
     * @brief Non corrected sum of all intensities.
     */
    float _peakArea;
    /**
     * @brief Area of spline
     */
    float _peakSplineArea;
    /**
     * @brief Baseline substracted area.
     */
    float _peakAreaCorrected;
    /**
     * @brief Top 3 points of the peak
     */
    float _peakAreaTop;
    /**
     * @brief Top 3 points of the peak
     * baseline subtracted.
     */
    float _peakAreaTopCorrected;
    /**
     * @brief Area of the peak divided by total
     * area in the EIC.
     */
    float _peakAreaFractional;
    /**
     * @brief Peak rank (sorted by peakAreaCorrected).
     */
    float _peakRank;
    /**
     * @brief Not corrected intensity at top of the peak.
     */
    float _peakIntensity;
    /**
     * @brief Baseline level below the highest point.
     */
    float _peakBaseLineLevel;
    /**
     * @brief mz value at the top of the peak.
     */
    float _peakMz;
    /**
     * @brief Averaged mz value across all points in the peak.
     */
    float _medianMz;
    /**
     * @brief mz value across base of the peak.
     */
    float _baseMz;
    /**
     * @brief From 0 to 1. indicator of peak goodness.
     */
    float _quality;
    /**
     * @brief Width of the peak at the baseline.
     */
    unsigned int _width;
    /**
     * @brief Fit to gaussian curve.
     */
    float _gaussFitSigma;
    /**
     * @brief Fit to gaussian curve.
     */
    float _gaussFitR2;
    int _groupNum;

    unsigned int _noNoiseObs;
    float _noNoiseFraction;
    float _symmetry;
    float _signalBaselineRatio;
    float _signalBaselineDifference;
    /**
     * @brief 0 no overlap, 1 perfect overlap.
     */
    float _groupOverlap;
    float _groupOverlapFrac;

    bool _localMaxFlag;

    /**
     * @brief True if peak is from blank sample.
     */
    bool _fromBlankSample;

    /**
     * @brief Classification label.
     */
    char _label;

    /**
     * @brief copyObj Copies one peak to another.
     * @param o  Object of class peak.
     */
    void copyObj(const Peak& o);
};
#endif
