#ifndef MZSAMPLE_H
#define MZSAMPLE_H

#include <chrono_io.h>
#include <date.h>

#include "assert.h"
#include "mzUtils.h"
#include "pugixml.hpp"
#include "standardincludes.h"

#ifdef ZLIB
#include <zlib.h>
#endif

#ifdef CDFPARSER
#include "ms10.h"
#endif

#if defined(WIN32) || defined(WIN64)
#ifndef strncasecmp
#define strncasecmp strnicmp
#endif
#define isnanwin(x) ((x) = (x))
#endif /* Def WIN32 or Def WIN64 */

class Scan;
class Peak;
class PeakGroup;
class EIC;
class Compound;
class Adduct;
class mzSlice;
class MassCalculator;
class MassCutoff;
class ChargedSpecies;

using namespace pugi;
using namespace mzUtils;
using namespace std;

/**
 * @brief Parses input sample files and stores related metadata
 *
 * @details Loads and parses sample files from mass spectrometry.
 * Supported file formats: .mzxml, .mzml, .mzdata, .mzcsv and .cdf
 *
 */
class mzSample
{
    public:
    /**
     * @brief Constructor for class mzSample
     */
    mzSample();

    /**
     * @brief Destructor for class mzSample
     */
    ~mzSample();

    /**
     * @brief Load sample (supported formats: .mzxml, .mzml,
     *  .mzdata, .mzcsv and .cdf)
     * @param filename Sample file name
     */
    void loadSample(string filename);

    /**
     * @brief Average time difference between scans
     * @return Average scan time
     */
    float getAverageFullScanTime();

    /**
     * @brief Find correlation between two EICs
     * @param mz1 m/z for first EIC
     * @param mz2 m/z for second EIC
     * @param ppm ppm window
     * @param rt1 Retention time for first EIC
     * @param rt2 Retention time for second EIC
     * @param eicType Type of EIC (max or sum)
     * @param filterline selected filterline
     * @return correlation
     */
    float correlation(float mz1,
                      float mz2,
                      MassCutoff* massCutoff,
                      float rt1,
                      float rt2,
                      int eicType,
                      string filterline);

    /**
     * @brief Get normalization constant
     * @return Normalization constant
     */
    float normalizationConstant() const
    {
        return _normalizationConstant;
    }

    /**
     * @brief Set normalization constant
     * @param x Normalization constant
     */
    void setNormalizationConstant(float x)
    {
        _normalizationConstant = x;
    }

    /**
     * @brief Get scan for a given scan number
     * @param scanNum Scan number
     * @return Scan class object
     * @see Scan
     */
    Scan* getScan(unsigned int scanNum);

    /**
     * @brief Get Average Scan
     * @param rtmin Minimum retention time
     * @param rtmax Maximum retention time
     * @param mslevel MS Level
     * @param polarity Ionization mode/polarity
     * @param resolution Mass resolution of the MS machine
     * @return Scan class object
     */
    Scan* getAverageScan(float rtmin,
                         float rtmax,
                         int mslevel,
                         int polarity,
                         float resolution);

    /**
     * @brief Get EIC based on minMz, maxMz, minRt, maxRt, mslevel
     * @param mzmin Minimum m/z
     * @param mzmax Maximum m/z
     * @param rtmin Minimum retention time
     * @param rtmax Maximum retention time
     * @param mslevel MS Level. MS Level is 1 for MS data and 2 for
     * MS/MS data
     * @param eicType Type of EIC (max or sum)
     * @param filterline selected filterline
     * @return EIC class object
     * @see EIC
     */
    EIC* getEIC(float mzmin,
                float mzmax,
                float rtmin,
                float rtmax,
                int mslevel,
                int eicType,
                string filterline);

    /**
     * @brief Get EIC based on srmId
     * @param srmId Filterline
     * @param eicType Type of EIC (max or sum)
     * @return EIC class object
     * @see EIC
     */
    EIC* getEIC(string srmId, int eicType);

    /**
     * @brief Get EIC for MS-MS dataset
     * @param precursorMz m/z of precursor Ion
     * @param collisionEnergy collision Energy
     * @param productMz m/z of product Ion]
     * @param eicType Type of EIC (max or sum)
     * @param filterline selected filterline
     * @param amuQ1 delta difference in Q1
     * @param amuQ3 delta difference in Q3
     * @return EIC class object
     */
    EIC* getEIC(float precursorMz,
                float productMz,
                int eicType,
                string filterline,
                float amuQ1,
                float amuQ3);

    /**
     * @brief Get Total Ion Chromatogram
     * @param rtmin Minimum retention time
     * @param rtmax Maximum retention time
     * @param mslevel MS level of the MS machine
     * @return EIC class object
     */
    EIC* getTIC(float rtmin, float rtmax, int mslevel);

    /**
     * @brief Get Base Peak Chromatogram
     * @param rtmin Minimum retention time
     * @param rtmax Maximum retention time
     * @param mslevel MS level of the MS machine
     * @return EIC class object
     */
    EIC* getBIC(float rtmin, float rtmax, int mslevel);

    /**
     * @brief save alignment state before next alignment is performed
     **/
    void saveCurrentRetentionTimes();

    /**
     * @brief restore last state of alignment
     * @details last saved state is restored on cancelling alignment
     **/
    void restorePreviousRetentionTimes();

    // TODO: Sahil, Added while merging projectdockwidget
    void applyPolynomialTransform();

    /**
     * @brief addScan Adds the parsed scan in scans vector.
     * @param s Scan to be added.
     */
    void addScan(Scan* s);

    /**
     * @brief getPolarity returns polarity
     * @return
     */
    int getPolarity();

    /**
     * @brief scanCount Returns the size of "scans" vector.
     * @return
     */
    inline unsigned int scanCount() const
    {
        return (scans.size());
    }

    /**
     * @brief Get the number of MS1 level scans detected in this sample.
     * @return MS1 scan count as unsigned int.
     */
    inline unsigned int ms1ScanCount()
    {
        return _numMS1Scans;
    }

    /**
     * @brief Get the number of MS2 level scans detected in this sample.
     * @return MS2 scan count as unsigned int.
     */
    inline unsigned int ms2ScanCount()
    {
        return _numMS2Scans;
    }

    /**
     * @brief Obtain the unique sample ID for the sample.
     * @return Sample ID as an integer.
     */
    inline int getSampleId()
    {
        return _id;
    }

    /**
     * @brief Set the sample ID for this sample.
     * @param id An integer which is not an ID for another sample.
     */
    inline void setSampleId(const int id)
    {
        _id = id;
    }

    /**
     * @brief getSampleName returns sampleName.
     * @return
     */
    inline string getSampleName() const
    {
        return sampleName;
    }

    /**
     * @brief setSampleOrder Sets sample order.
     * @param x
     */
    void setSampleOrder(int x)
    {
        _sampleOrder = x;
    }

    /**
     * @brief getSampleOrder Returns sampleOrder.
     * @return
     */
    inline int getSampleOrder() const
    {
        return _sampleOrder;
    }

    /**
     * @brief the order in  which the samples are inserted into the LC column.
     * @method setInjectionOrder
     * @param time
     */
    void setInjectionOrder(int order)
    {
        injectionOrder = order;
    }

    /**
     * @brief return the injection order.
     * @method getInjectionOrder.
     * @return injectionOrder.
     */

    int getInjectionOrder()
    {
        return injectionOrder;
    }

    /**
     * @brief getSetName Returns setName.
     * @return
     */
    inline string getSetName() const
    {
        return _setName;
    }

    /**
     * @brief getSetName sets setName.
     * @return
     */
    void setSetName(string x)
    {
        _setName = x;
    }

    /**
     * @brief setSampleName sets sampleName.
     * @param x
     */
    void setSampleName(string x)
    {
        sampleName = x;
    }

    /**
     * @brief getMaxRt Returns maxRt from vector of samples.
     * @param samples
     * @return
     */
    static float getMaxRt(const vector<mzSample*>& samples);

    /**
     * @brief find all MS2 scans within the slice.
     * @return vector of all matching MS2 scans.
     */
    vector<Scan*> getFragmentationEvents(mzSlice* slice);

    /**
     * [C13Labeled?]
     * @method C13Labeled
     * @return [true or false]
     */
    bool C13Labeled() const
    {
        return _C13Labeled;
    }

    /**
     * [N15Labeled]
     * @method N15Labeled
     * @return [true or false]
     */
    bool N15Labeled() const
    {
        return _N15Labeled;
    }

    /**
     * @brief compSampleOrder True if sample order for sample a is
     * less than b else false.
     * @param a
     * @param b
     * @return
     */
    static bool compSampleOrder(const mzSample* a, const mzSample* b)
    {
        return a->_sampleOrder < b->_sampleOrder;
    }

    /**
     * @brief Compare sample order of two samples
     * @param a object of class mzSample
     * @param b object of class mzSample
     * @return True if sample order of 'a' is higher than 'b' else false
     **/
    static bool compRevSampleOrder(const mzSample* a, const mzSample* b)
    {
        return a->_sampleOrder > b->_sampleOrder;
    }

    /**
     * @brief compSampleSort Compares AlphaNumeric SampleNames.
     * @param a object of class mzSample.
     * @param b object of class mzSample.
     * @return
     */
    static bool compSampleSort(const mzSample* a, const mzSample* b)
    {
        return mzUtils::strcasecmp_withNumbers(a->sampleName, b->sampleName);
    }

    /**
     * @brief Compare injection time (in epoch seconds) of two samples.
     * @param a object of class mzSample.
     * @param b object of class mzSample.
     * @return True if mzSample a has lower injection time than mzSample b.
     */
    static bool compInjectionTime(const mzSample* a, const mzSample* b)
    {
        return a->injectionTime < b->injectionTime;
    }

    /**
     * @brief setFilter_minIntensity sets filter_minIntensity.
     * @param x
     */
    static void setFilter_minIntensity(int x)
    {
        filter_minIntensity = x;
    }

    /**
     * @brief setFilter_centroidScans   sets filter_centroidScans.
     * @param x
     */
    static void setFilter_centroidScans(bool x)
    {
        filter_centroidScans = x;
    }

    /**
     * @brief setFilter_intensityQuantile sets filter_intensity_Quantile.
     * @param x
     */
    static void setFilter_intensityQuantile(int x)
    {
        filter_intensityQuantile = x;
    }

    /**
     * @brief setFilter_mslevel sets filter_mslevel.
     * @param x
     */
    static void setFilter_mslevel(int x)
    {
        filter_mslevel = x;
    }

    /**
     * @brief setFilter_polarity sets filter_polarity.
     * @param x
     */
    static void setFilter_polarity(int x)
    {
        filter_polarity = x;
    }

    /**
     * @brief getFilter_minIntensity returns filter_minIntensity.
     * @return
     */
    static int getFilter_minIntensity()
    {
        return filter_minIntensity;
    }

    /**
     * @brief getFilter_intensityQuantile returns filter_intensityQuantile.
     * @return
     */
    static int getFilter_intensityQuantile()
    {
        return filter_intensityQuantile;
    }

    /**
     * @brief getFilter_centroidScans  returns filter_centroidScans.
     * @return
     */
    static int getFilter_centroidScans()
    {
        return filter_centroidScans;
    }

    /**
     * @brief getFilter_mslevel returns filter_mslevel.
     * @return
     */
    static int getFilter_mslevel()
    {
        return filter_mslevel;
    }

    /**
     * @brief getFilter_polarity returns filter_polarity.
     * @return
     */
    static int getFilter_polarity()
    {
        return filter_polarity;
    }

    /**
     * @brief getIntensityDistribution
     * @param mslevel
     * @return
     */
    vector<float> getIntensityDistribution(int mslevel);

    deque<Scan*> scans;
    string sampleName;
    string fileName;
    bool isSelected;
    bool isBlank;

    /**
     * sample display color, [r,g,b, alpha]
     */
    float color[4];
    float minMz;
    float maxMz;
    float minRt;
    float maxRt;
    float maxIntensity;
    float minIntensity;
    float totalIntensity;
    // Sample display order
    int _sampleOrder;

    int sampleNumber;

    bool _C13Labeled;
    bool _N15Labeled;
    // Feng note: added to track S34 labeling state
    bool _S34Labeled;
    // Feng note: added to track D2 labeling state
    bool _D2Labeled;

    float _normalizationConstant;
    string _setName;

    unsigned long int injectionTime;
    int injectionOrder;

    map<string, vector<int>> srmScans;

    /**
     * tags associated with this sample
     */
    map<string, string> instrumentInfo;

    vector<float> lastSavedRTs;
    vector<double> polynomialAlignmentTransformation;

    private:
    int _id;
    unsigned int _numMS1Scans;
    unsigned int _numMS2Scans;

    /**
     * @brief Parse mzData file format
     * @param char* mzData file name
     */
    void parseMzData(string);

    /**
     * @brief Parse mzCSV file format
     * @param char* mzCSV file name
     */
    void parseMzCSV(string);

    /**
     * @brief Parse mzXML file format
     * @param char* mzXML file name
     */
    void parseMzXML(string);

    /**
     * @brief Parse mzML file format
     * @param char* mzML file name
     */
    void parseMzML(string);

    /**
     * @brief Parse netcdf file format
     * @param filename netcdf file name
     * @param is_verbose verbose or not
     * @return Returns 0 if error in loading else 1
     */
    int parseCDF(string filename, int is_verbose);

    void sampleNaming(string filename);

    void checkSampleBlank(string filename);

    void setInstrumentSettigs(xml_node spectrumstore);

    void parseMzXMLData(const xml_node& spectrumstore);

    /**
     * @brief Parse scan in mzXml file format
     * @param scan xml_node object of pugixml library
     * @param scannum scan number
     */
    void parseMzXMLScan(const xml_node& scan, const int& scannum);

    xml_node getmzXMLSpectrumData(xml_document& doc, string filename);

    float parseRTFromMzXML(xml_attribute& attr);

    static int parsePolarityFromMzXML(xml_attribute& attr);

    static int getPolarityFromfilterLine(string filterLine);

    vector<float> parsePeaksFromMzXML(const xml_node& scan);

    void populateMzAndIntensity(const vector<float>& mzint, Scan* _scan);

    void populateFilterline(const string& filterLine, Scan* _scan);

    /**
     * @brief Update injection time stamp
     * @param xml_node xml_node object of pugixml library
     */
    void parseMzMLInjectionTimeStamp(const xml_attribute&);

    /**
     * @brief Parse mzML chromatogrom list
     * @param xml_node xml_node object of pugixml library
     */
    void parseMzMLChromatogramList(const xml_node&);

    /**
     * @brief Get cv parameters of an xml_node from mzML format
     * @param node xml_node object of pugixml lirary
     * @return Returns a map of name and corresponding value of an
     * xml_node object
     */
    static map<string, string> mzML_cvParams(xml_node node);

    int getSampleNoChromatogram(const string& chromatogramId);

    void cleanFilterLine(string& filterline);

    /**
     * @brief Parse mzML spectrum list
     * @param xml_node xml_node object of pugixml library
     */
    void parseMzMLSpectrumList(const xml_node&);

    /**
     * @brief Map scan numbers to filterline
     * @details Update map srmScans where key is the filterline and value
     * is int vector.
     * int vector contains scan numbers
     * @see mzSample:srmScans
     */
    void enumerateSRMScans();

    /**
     * @brief Compute min and max values for mz and rt
     * @details Compute min and max values for mz and rt by iterating over
     * the scans. Also while iterating over scans, calculate min and max
     * Intensity
     */
    void calculateMzRtRange();

    // TODO: This should be moved
    static string getFileName(const string& filename);
    static int filter_minIntensity;
    static bool filter_centroidScans;
    static int filter_intensityQuantile;
    static int filter_mslevel;
    static int filter_polarity;

    vector<string> filterChromatogram{"sample",
                                      "start",
                                      "end",
                                      "index",
                                      "scan"};
};

#endif
