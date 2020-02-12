#ifndef SCAN_H
#define SCAN_H

#include <QRegExp>
#include <QString>
#include <QStringList>

#include "standardincludes.h"

class mzSample;
class mzPoint;
class ChargedSpecies;
class MassCutoff;

using namespace std;

/**
* @class Scan
* @ingroup libmaven
* @brief Scan contains information about retention time, m/z, intensity
* @details Each scan stores  m/z and it's corresponding intensity for a
* particular Retention time.
* There are multiple m/z and intenstities in one scan
*/

class Scan
{
  public:

      /**
       * @brief Scan Default Constructor
       */
      Scan(){

      }
    /**
     * @brief Scan  Parameterised constructor.
     * @param sample
     * @param scannum
     * @param mslevel
     * @param rt
     * @param precursorMz
     * @param polarity
     */
    Scan(mzSample *sample, int scannum, int mslevel,
         float rt, float precursorMz, int polarity);

    /**
     * @brief deepcopy Copies one object to another
     * along with pointer values.
     * @param b object to be copied
     */
    void deepcopy(Scan *b);

    /**
     * @brief return number of m/z's(number of
     * observatiosn) recorded in a scan.
     */
    inline unsigned int nobs() const{
        return mz.size();
    }

    /**
     * @brief Obtain the smallest m/z value stored.
     * @return Fractional m/z value.
     */
    inline float minMz() {
        if(nobs() > 0)
            return *(std::min_element(begin(mz),
                                      end(mz)));
        return 0.0f;
    }

    /**
     * @brief Obtain the largest m/z value stored.
     * @return Fractional m/z value.
     */
    inline float maxMz() {
        if(nobs() > 0)
            return *(std::max_element(begin(mz),
                                      end(mz)));
        return 0.0f; }

    /**
     *@brief return the corresponding sample
     */
    inline mzSample* sample(){
        return _sample;
    }

    //TODO why int vector?
    /**
     * @brief findMatchingMzs   Finds the vector of mz in given range.
     * @param mzmin             Defines lower bound of the range.
     * @param mzmax             Defines upper bound of the range.
     * @return
     */
    vector<int> findMatchingMzs(float mzmin, float mzmax);

    /**
     *@brief In a given m/z range (mz + ppm) return the position
     * of Highest intensity present in the scan
     *@param m/z and ppm(parts per million). Together they define
     * the m/z range
     */
    int findHighestIntensityPos(float mz, MassCutoff *massCutoff);

    //TODO: Sahil, Added while merging point
    /**
     * @brief findClosestHighestIntensityPos  Finds the highst closest mz to
     * the given mz.
     * @param mz
     * @param massCutoff
     * @return
     */
    int findClosestHighestIntensityPos(float mz, MassCutoff *massCutoff);

    /**
    * @brief checks if an input m/z value is present in the scan.
    * @param  mz(input m/z value) and ppm(prats per million).
    * @return returns true if input m/z exists, otherwise false.
    */
    bool hasMz(float mz, MassCutoff *massCutoff);

    /**
    * @brief check if the data is centroided.
    * @return return true if data is centroided else false.
    */
    bool isCentroided() const{
        return _centroided;
    }

    /**
     * @brief setCentroided set _centroided variable
     * @param centre
     */
    void setCentroided(bool centre){
        _centroided = centre;
    }

    /**
     * @brief polarity of the scan defines whether the metabolites
     * were positively or negatively charged
     * @return returns the polarity of the scan (+1 for positive, -1 for
     * negative and 0 for neutral)
     */
    inline int getPolarity() const {
        return _polarity;
    }

    /**
    * @brief set the polarity of scan
    * @param  x is the polarity to be set (+1 for positive, -1
    * for negative and 0 for neutral)
    */
    void setPolarity(int x) {
        _polarity = x;
    }

    /**
    * @brief Calculate the sum of all the intensities for a scan
    * @return return total intensity
    */
    int totalIntensity() const
    {
        int sum = 0;
        for (unsigned int i = 0; i < intensity.size(); i++)
            sum += intensity[i];
        return sum;
    }


    //TODO: basePeakIntensity value in every scan of mzXml file represents
    //maxIntensity. Use it rather than looping over all the  intensities
    /**
    * @brief return the maxIntensity in scan
    */

    float maxIntensity()
    {
        float max = 0;
        for (unsigned int i = 0; i < intensity.size(); i++)
            if (intensity[i] > max)
                max = intensity[i];
        return max;
    }

    /**
    * @brief return pairs of m/z, intensity values for top intensities.
    * @details intensities are normalized to a maximum intensity in a
    * scan * 100]
    * @param minFracCutoff specfies mininum relative intensity; for example
    * 0.05, filters out all intensites below 5% of maxium scan intensity]
    */
    vector<pair<float, float> > getTopPeaks(float minFracCutoff,
                                          float minSigNoiseRatio, int dropTopX);

    vector<int> assignCharges(MassCutoff *massCutoffTolr);

    ChargedSpecies *deconvolute(float mzfocus,
                                float noiseLevel,
                                MassCutoff *massCutoffMerge,
                                float minSigNoiseRatio,
                                int minDeconvolutionCharge,
                                int maxDeconvolutionCharge,
                                int minDeconvolutionMass,
                                int maxDeconvolutionMass,
                                int minChargedStates);

    string toMGF();

    /**
    *@brief  return position of intensities in descending
    * order(highest to lowest)
    */
    vector<int> intensityOrderDesc();

    /**
         * intensities found in one scan
         */
    vector<float> intensity;
    /**
         * m/z's found in one scan
         */
    vector<float> mz;

    /**
    *@brief centroid the data
    *@see also refer to findLocalMaximaInIntensitySpace
    */
    void simpleCentroid();

    /**
    * @brief removes intensities from scan that are  lower than minIntensity
    */
    void intensityFilter(int minIntensity);

    /**
    * @brief removes intensities from scan that are  lower than minQuantile
    */
    void quantileFilter(int minQuantile);

    /**
    * @brief adjusts precursor m/z for MS2 scans
    * @details precursor m/z available in the MS2 scans have
    * lower precision than m/z values recorded in the MS1 fullscans
    * To ensure correct mapping of MS2 scans to parent scans, precursor
    * m/z needs to be adjusted
    * @param ppm lowest ppm range for finding the correct
    * precursor m/z
    */
    void recalculatePrecursorMz(float ppm);

    /**
     * @brief calculates purity of the spectra
     * @details if the parent full scan has multiple readings
     * within a precursor m/z window
     * the fragmentation scan would be a mixture of fragments
     * from all those species
     */ 
    double getPrecursorPurity(float ppm = 10.0);

    /**
    *@brief print the info present in a scan
    */
    void summary();



    /**
     * @brief compare total intensity of two scans
     * @return true if Scan a has a higher totalIntensity
     * than b, else false
     */
    static bool compIntensity(Scan* a, Scan* b){
        return a->totalIntensity() > b->totalIntensity();
    }

    /**
    * @brief compare retention times of two scans
    * @return return True if Scan a has lower retention time
    * than Scan b, else false
    */
    static bool compRt(Scan *a, Scan *b){
        return a->_rt < b->_rt;
    }

    /**
    * @brief compare precursor m/z of two samples
    * @return return True if Scan a has lower precursor
    * m/z than Scan b, else false
    */
    static bool compPrecursor(Scan *a, Scan *b){
        return a->_precursorMz < b->_precursorMz;
    }

    bool operator<(const Scan &b) const { return _rt < b._rt; }

    /**
     * @brief Getters and Setters for the
     * private datamembers.
     * @return
     */
    int mslevel() const
    {
        return _mslevel;
    }

    void setMslevel(int level)
    {
        _mslevel = level;
    }

    void setScannum(int scan)
    {
        _scannum = scan;
    }

    int scannum()   const
    {
        return _scannum;
    }

    void setProductMz(float mz){
        _productMz = mz;
    }

    float productMz()   const
    {
        return _productMz;
    }

    void setFilterLine(string filter){
        _filterLine = filter;
    }

    string filterLine() const
    {
        return _filterLine;
    }

    string scanType()   const
    {
        return _scanType;
    }
    void setScanType(string type)
    {
        _scanType = type;
    }

    float rt()  const
    {
        return _rt;
    }

    void setRt(float rt)
    {
        _rt = rt;
    }

    float precursorMz() const
    {
        return _precursorMz;
    }

    void setPrecursorMz(float mz)
    {
        _precursorMz = mz;
    }

    float precursorIntensity()  const
    {
        return _precursorIntensity;
    }

    void setPrecursorIntensity(float it)
    {
        _precursorIntensity = it;
    }

    int precursorCharge()   const
    {
        return _precursorCharge;
    }

    void setPrecursorCharge(int charge)
    {
        _precursorCharge = charge;
    }

    int precursorScanNum()  const
    {
        return _precursorScanNum;
    }

    void setPrecursorScanNum(int num)
    {
        _precursorScanNum = num;
    }

    void setIsolationWindow(float win)
    {
        _isolationWindow = win;
    }

    float isolationWindow() const
    {
        return _isolationWindow;
    }

    void setCollisionEnergy(float energy)
    {
        _collisionEnergy = energy;
    }

    float collisionEnergy() const
    {
        return _collisionEnergy;
    }

    float originalRt()  const
    {
        return _originalRt;
    }

    void setOriginalRt(float rt)
    {
        _originalRt = rt;
    }

  private:
        int _mslevel;
        bool _centroided;
        /**
         * retention time at which the scan was recorded
         */
        float _rt;
        /**
         * originalRt will hold original retention time when rt
         * is modified in case of alignment
         */
        float _originalRt;
        int _scannum;

        float _precursorMz;
        float _precursorIntensity;
        int _precursorCharge;
        int _precursorScanNum;
        /**
         * @brief precursor mass resolution for fragmentation event
         */
        float _isolationWindow;
        float _productMz;
        float _collisionEnergy;

        string _scanType;
        string _filterLine;
        /**
         * sample corresponding to the scan
         */
        mzSample* _sample;
        /**
         * +1 for positively charged, -1 for negatively
         * charged, 0 for neutral
         */
        int _polarity;



        float _parentPeakIntensity;

        struct _BrotherData
        {
            float expectedMass;
            int countMatches;
            float totalIntensity;
            int upCount;
            int downCount;
            int minZ;
            int maxZ;
        };

        ofstream file;
        _BrotherData *brotherdata, b;
        void initialiseBrotherData(int z, float mzfocus);
        void updateBrotherDataIfPeakFound(int loopdirection,
                                          int ii, bool *flag,
                                          bool *lastMatched,
                                          float *lastIntensity,
                                          float noiseLevel,
                                          MassCutoff *massCutoffMerge);

        void updateChargedSpeciesDataAndFindQScore(ChargedSpecies *x, int z,
                                                   float mzfocus, float noiseLevel,
                                                   MassCutoff *massCutoffMerge,
                                                   int minChargedStates);

        void findBrotherPeaks(ChargedSpecies *x, float mzfocus,
                              float noiseLevel,
                              MassCutoff *massCutoffMerge,
                              int minDeconvolutionCharge,
                              int maxDeconvolutionCharge,
                              int minDeconvolutionMass,
                              int maxDeconvolutionMass,
                              int minChargedStates);

        bool setParentPeakData(float mzfocus, float noiseLevel,
                               MassCutoff *massCutoffMerge,
                               float minSigNoiseRatio);

        void findError(ChargedSpecies *x);

        /**
         * @brief gets the previous MS1 scan till historySize
         */
        Scan* getLastFullScan(int historySize = 50);

        /**
         * @brief gets the full-scan m/z-int readings that fall within
         * the isolation window of the precursor
         */
        vector<mzPoint> getIsolatedRegion(float isolationWindowAmu = 1.0);

        vector<float> smoothenIntensitites();
        /**
        *@brief compare the current intensity with the left and the right
        * intensity and store the largest amongest them
        */
        void findLocalMaximaInIntensitySpace(int vsize, vector<float> *cMz,
                                             vector<float> *cIntensity,
                                             vector<float> *spline);

        void updateIntensityWithTheLocalMaximas(vector<float> *cMz,
                                                vector<float> *cIntensity);
};
#endif
