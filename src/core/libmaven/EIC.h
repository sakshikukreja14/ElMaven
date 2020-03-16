/**
 * @class EIC
 * @ingroup libmaven
 * @brief Wrapper class for a eic.
 * @author Sabu George
 */
#ifndef MZEIC_H
#define MZEIC_H

#include <Eigen>

#include "standardincludes.h"

class Peak;
class PeakGroup;
class mzSample;
class mzPoint;
class Scan;
class Compound;
class mzSlice;

using namespace std;

class EIC
{
    public:
    /**
     *  Default constructor.
     */
    EIC();

    /**
     *  Destructor
     */
    ~EIC();

    /**
     * @brief The SmootherType enum Enumeration to select
     * the smoothing algorithm.
     */
    enum SmootherType
    {
        SAVGOL = 0,
        GAUSSIAN = 1,
        AVG = 2
    };

    /**
     * @brief The EicType enum Enumeration to select how
     * intensity and/or mass is calculated at a particular
     * retention time.
     */
    enum EicType
    { MAX = 0,
      SUM = 1 };

    /**
     * @brief scannum  Store all scan numbers in an EIC.
     */
    vector<int> scannum;

    /**
     * @brief rt  Store all retention times in an EIC.
     */
    vector<float> rt;

    /**
     * @brief mz  Store all mass/charge ratios in an EIC.
     */
    vector<float> mz;

    /**
     * @brief intensity  Store all intensities in an EIC.
     */
    vector<float> intensity;

    /**
     * @brief peaks  Store all peak objects in an EIC.
     */
    vector<Peak> peaks;

    /**
     * @brief color Color of the eic line, [r,g,b, alpha].
     */
    float color[4];

    /**
     * @brief spline Pointer to smoothed intensity array.
     */
    float* spline;

    /**
     * @brief baseline Pointer to baseline array.
     */
    float* baseline;

    enum class BaselineMode { Threshold, AsLSSmoothing };

    /**
     * @brief find peak positions after smoothing, baseline calculation and peak
     * filtering
     * @param smoothWindow number of scans used for smoothing in each iteration
     */
    void getPeakPositions(int smoothWindow);

    /**
     * @brief set values for all members of a peak object
     * @param  peak peak object
     */
    void getPeakDetails(Peak& peak);



    void setBaselineMode(BaselineMode b)
    {
        _baselineMode = b;
    }

    /**
     * @brief Calculate baseline for the current baseline mode.
     */
    void computeBaseline();

    /**
     * @brief calculate spline of the EIC
     * @details smoothen intensity data according to selected algorithm. stores
     * it as spline
     * @param  smoothWindow  number of scans used for smoothing in each
     * iteration
     */
    void computeSpline(int smoothWindow);

    /**
     * @brief find the first and last position of a peak
     * @param  peak peak object
     */
    void findPeakBounds(Peak& peak);



    /**
     * @brief get vector of all intensity points in a peak
     * @param peak peak object
     * @return mzPoint vector of intensity points in the peak
     */
    vector<mzPoint> getIntensityVector(Peak& peak);

    /**
     * @brief print parameter values of an EIC in log window
     */
    void summary();

    /**
     * @brief set smoothing algorithm
     * @param  x SmootherType
     */
    void setSmootherType(EIC::SmootherType x)
    {
        smootherType = x;
    }

    EIC::SmootherType getSmootherType(){
        return smootherType;
    }
    /**
     * @brief set smoothing window for baseline
     * @param  x number of scans used for smoothing in one iteration
     */
    void setBaselineSmoothingWindow(int x)
    {
        baselineSmoothingWindow = x;
    }

    /**
     * @brief set percentage of top intensity points to remove for setting
     * baseline
     * @param  x percentage of top intensity points to remove
     */
    void setBaselineDropTopX(int x)
    {
        baselineDropTopX = x;
    }

    /**
     * @brief Set smoothness (λ) to be used for default AsLS baseline
     * estimation.
     * @param s smoothness (will be mutated to 10^s when actually used)
     */
    void setAsLSSmoothness(int s)
    {
        _aslsSmoothness = s;
    }

    /**
     * @brief Set asymmetry (p) to be used for default AsLS baseline estimation.
     * @param a asymmetry value (will be divided by 100 when actually used).
     */
    void setAsLSAsymmetry(int a)
    {
        _aslsAsymmetry = a;
    }

    /**
     * @brief set minimum signal baseline difference for every peak
     * @param x signal baseline difference threshold for every peak
     */
    void setFilterSignalBaselineDiff(double x)
    {
        _filterSignalBaselineDiff = x;
    }

    /**
     * @brief get EIC of a sample using given mass/charge and retention time
     * range
     * @details
     * @param
     * @return bool true if EIC is pulled. false otherwise
     */
    bool makeEICSlice(mzSample* sample,
                      float mzmin,
                      float mzmax,
                      float rtmin,
                      float rtmax,
                      int mslevel,
                      int eicType,
                      string filterline);

    void getRTMinMaxPerScan();

    void normalizeIntensityPerScan(float scale);


//    void clearEICContents();
 //   void interpolate();
    /**
     * [size ]
     * @method size
     * @return []
     */
    inline unsigned int size()
    {
        return intensity.size();
    }

    /**
     * @return sample associated with the EIC
     */
    inline mzSample* getSample()
    {
        return _sample;
    }

    /**
     * @brief return list of groups given a set of EICs
     * @details assigns every peak to a group based on the best matching merged
     *EIC
     * @param smoothingWindow number of scans read at a time for smoothing
     * @param maxRtDiff maximum retention time difference between peaks in a
     *group
     * @param minQuality minimum peak quality for every group. used for
     *calculation of good peaks
     * @param distXWeight weight of rt difference between a peak and merged EIC
     * @param distYWeight weight of intensity difference between a peak and
     *merged EIC
     * @param overlapWeight weight of peak overlap between a peak and merged EIC
     * @param userOverlap flag to determine which group score formula is used
     * @param minSignalBaselineDifference minimum difference between peak and
     *baseline for peak to be marked
     * @return vector of peak groups found
     **/
    static vector<PeakGroup> groupPeaks(vector<EIC*>& eics,
                                        mzSlice* slice,
                                        int smoothingWindow,
                                        float maxRtDiff,
                                        double minQuality,
                                        double distXWeight,
                                        double distYWeight,
                                        double overlapWeight,
                                        bool useOverlap,
                                        double minSignalBaselineDifference,
                                        float fragmentPpmTolerance,
                                        string scoringAlgo);
    /**
     * [eicMerge ]
     * @method eicMerge
     * @param  eics     []
     * @return []
     */
    static EIC* eicMerge(const vector<EIC*>& eics);

    /**
     * [remove Low Rank Groups ]
     * @method removeLowRankGroups
     * @param  groups              [vector of peak groups]
     * @param  rankLimit           [group rank limit ]
     */
    static void removeLowRankGroups(vector<PeakGroup>& groups,
                                    unsigned int rankLimit);

    /**
     * [compare Max Intensity]
     * @method compMaxIntensity
     * @param  a                [EIC a]
     * @param  b                [EIC b]
     * @return [true or false]
     */
    static bool compMaxIntensity(EIC* a, EIC* b)
    {
        return a->_maxIntensity > b->_maxIntensity;
    }

    /**
     * @brief Getters and Setters for private data members.
     */


    BaselineMode baselineMode()
    {
        return _baselineMode;
    }

    void setSampleName(string name)
    {
        _sampleName = name;
    }

    string sampleName()
    {
        return _sampleName;
    }

    void setSample(mzSample * sample)
    {
        _sample = sample;
    }

    void setMaxIntensity(float intensity)
    {
        _maxIntensity = intensity;
    }

    float maxIntensity()
    {
        return _maxIntensity;
    }

    void setRtAtMaxIntensity(float rt)
    {
        _rtAtMaxIntensity = rt;
    }

    float rtAtMaxIntensity()
    {
        return _rtAtMaxIntensity;
    }

    void setMzAtMaxIntensity(float mz)
    {
        _mzAtMaxIntensity = mz;
    }

    float mzAtMaxIntensity()
    {
        return _mzAtMaxIntensity;
    }

    void setMaxAreaTopIntensity(float intensity)
    {
        _maxAreaTopIntensity = intensity;
    }

    float maxAreaTopIntensity()
    {
        return _maxAreaTopIntensity;
    }

    void setMaxAreaIntensity(float intensity)
    {
        _maxAreaIntensity = intensity;
    }

    float maxAreaIntensity()
    {
        return _maxAreaIntensity;
    }

    void setMaxAreaNotCorrectedIntensity(float intensity)
    {
        _maxAreaNotCorrectedIntensity = intensity;
    }

    float maxAreaNotCorrectedIntensity()
    {
        return _maxAreaNotCorrectedIntensity;
    }

    void setMaxAreaTopNotCorrectedIntensity(float intensity)
    {
        _maxAreaTopNotCorrectedIntensity = intensity;
    }

    float maxAreaTopNotCorrectedIntensity()
    {
        return _maxAreaTopNotCorrectedIntensity;
    }

    float filterSignalBaselineDiff()
    {
        return _filterSignalBaselineDiff;
    }

    void setTotalIntensity(float intensity)
    {
        _totalIntensity = intensity;
    }

    float totalIntensity()
    {
        return _totalIntensity;
    }

    void setEic_noNoiseObs(int num)
    {
        _eic_noNoiseObs = num;
    }

    int eic_noNoiseObs()
    {
        return _eic_noNoiseObs;
    }

    void setMzmin(float min)
    {
        _mzmin = min;
    }

    float mzmin()
    {
        return _mzmin;
    }

    void setMzmax(float max)
    {
        _mzmax = max;
    }

    float mzmax()
    {
        return _mzmax;
    }

    void setRtmin(float min)
    {
        _rtmin = min;
    }

    float rtmin()
    {
        return _rtmin;
    }

    void setRtmax(float max)
    {
        _rtmax = max;
    }

    float rtmax()
    {
        return _rtmax;
    }

    private:

    /**
     * Name of selected smoothing algorithm
     */
    SmootherType smootherType;

    /**
     * @brief _baselineMode decides which algorithm to use for computing
     * baseline.
     */
    BaselineMode _baselineMode;

    /**
     * Sets the number of scans used for smoothing in one iteration.
     */
    int baselineSmoothingWindow;

    /*
     * Percentage of top intensity points to remove before computing baseline.
     */
    int baselineDropTopX;

    /**
     * @brief Smoothness parameter for AsLS Smoothing algorithm
     */
    int _aslsSmoothness;

    /**
     * @brief Asymmetry parameter for AsLS Smoothing algorithm
     */
    int _aslsAsymmetry;

    /**
     * @brief sampleName  Store name of the sample associated
     * with the EIC.
     */
    string _sampleName;

    /**
     * @brief sample  Pointer to originating sample.
     */
    mzSample* _sample;

    /**
     * @brief maxIntensity  Maximum intensity of all scans.
     */
    float _maxIntensity;

    /**
     * @brief rtAtMaxIntensity  Rt value of maximum intensity.
     */
    float _rtAtMaxIntensity;

    /**
     * @brief mzAtMaxIntensity  Mz value of maximum intensity.
     */
    float _mzAtMaxIntensity;

    /**
     * @brief maxAreaTopIntensity  Maximum areaTop intensity
     * (after baseline correction) out of all peaks.
     */
    float _maxAreaTopIntensity;

    /**
     * @brief maxAreaIntensity  Maximum area intensity
     * (after baseline correction) out of all peaks.
     */
    float _maxAreaIntensity;

    /**
     * @brief maxAreaNotCorrectedIntensity  Maximum area intensity
     * (without baseline correction) out of all peaks.
     */
    float _maxAreaNotCorrectedIntensity;

    /**
     * @brief maxAreaTopNotCorrectedIntensity  Maximum areaTop intensity
     * (without baseline correction) out of all peaks.
     */
    float _maxAreaTopNotCorrectedIntensity;

    /**
     * @brief filterSignalBaselineDiff  Minimum threshold for peak
     * signal-baseline difference.
     */
    double _filterSignalBaselineDiff;

    /**
     * @brief totalIntensity  Sum of all intensities in an EIC.
     */
    float _totalIntensity;

    /**
     * @brief eic_noNoiseObs  Number of observations above baseline.
     */
    int _eic_noNoiseObs;

    /**
     * @brief mzmin  Minimum mass/charge for pulling an EIC.
     */
    float _mzmin;

    /**
     * @brief mzmax  Maximum mass/charge for pulling an EIC.
     */
    float _mzmax;

    /**
     * @brief rtmin  Minimum retention time for pulling an EIC.
     */
    float _rtmin;

    /**
     * @brief rtmax  Maximum retention time for pulling an EIC.
     */
    float _rtmax;

    /**
     * @brief add peak object to vector
     * @details create a peak object for given peak position and append to
     * vector
     * @param  peakPos position of peak in spline
     * @return pointer to newly added peak object in the vector
     */
    Peak* addPeak(int peakPos);

    /**
     * @brief width of given peak in terms of number of scans
     * @details find the number of scans where peak intensity is above baseline
     * @param  peak peak object
     */
    void getPeakWidth(Peak& peak);

    /**
     * @brief find all parameter values for every peak in an EIC
     * @method getPeakStatistics
     */
    void getPeakStatistics();

    /**
     * @brief find all peaks in an EIC
     * @details find all local maxima in an EIC and save them as objects
     */
    void findPeaks();

    /**
     * @brief remove peaks with parameter values below user-set thresholds
     */
    void filterPeaks();

    /**
     * brief
     * @param  peak             [peak]
     */
    void checkGaussianFit(Peak& peak);

 //   void subtractBaseLine();
    /**
     * @brief Clear the baseline if exists and reallocate memory for a new one.
     * @return Whether the baseline should be processed further or not.
     */
    bool _clearBaseline();

    /**
     * @brief Computes a baseline using naive thresholding method.
     * @param smoothingWindow is the size of window used for 1D guassian
     * smoothing.
     * @param dropTopX percent of the highest intensities will be truncated.
     */
    void _computeThresholdBaseline(const int smoothingWindow,
                                   const int dropTopX);

    /**
     * @brief Computes a baseline using Asymmetric Least Squares Smoothing
     * techinique.
     * @details A (Whittaker) smoother is used to get a slowly varying estimate
     * of the baseline. In contrast to ordinary least squares smoothing,
     * however, positive deviations with respect to baseline estimate are
     * weighted (much) less than negative ones.
     *
     * Ref: Baseline Correction with Asymmetric Least Squares Smoothing,
     * P. Eilers, H. Boelens, 2005
     *
     * @param lambda for smoothness. Typical values of lambda for MS data range
     * from 10^2 to 10^9, depending on dataset. But since we resample the
     * intensity signal, we can limit this range to [10^0, 10^3]. The exponent
     * value should be passed here as integer, i.e. lambda should be in range
     * [0, 3].
     * @param p for asymmetry. Values between 0.01 to 0.10 work reasonable well
     * for MS data.
     * @param numIterations for the number of iterations that should be
     * performed (since this is an iterative optimization algorithm).
     */
    void _computeAsLSBaseline(const float lambda,
                              const float p,
                              const int numIterations = 10);
};
#endif  // MZEIC_H
