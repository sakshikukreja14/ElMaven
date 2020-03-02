#include "Peak.h"
#include "EIC.h"
#include "Scan.h"
#include "doctest.h"
#include "mzSample.h"
#include "mzUtils.h"
#include "mavenparameters.h"
#include "mzSample.h"
#include "classifierNeuralNet.h"
#include "databases.h"
#include "PeakDetector.h"
#include "masscutofftype.h"
#include "PeakGroup.h"

using namespace mzUtils;

Peak::Peak()
{
    _pos = 0;
    _eic = NULL;
    _sample = NULL;
    _baseMz = 0;
    _fromBlankSample = false;
    _groupOverlap = 0;
    _groupNum = 0;
    _groupOverlapFrac = 0;
    _label = 0;
    _localMaxFlag = 0;
    _maxscan = 0;
    _medianMz = 0;
    _minpos = _pos;
    _maxpos = _pos;
    _splineminpos = _pos;
    _splinemaxpos = _pos;
    _minscan = 0;
    _mzmax = 0;
    _mzmin = 0;
    _noNoiseObs = 0;
    _noNoiseFraction = 0;
    _peakArea = 0;
    _peakSplineArea = 0;
    _peakAreaCorrected = 0;
    _peakAreaFractional = 0;
    _peakAreaTop = 0;
    _peakAreaTopCorrected = 0;
    _peakIntensity = 0;
    _peakMz = 0;
    _peakBaseLineLevel = 0;
    _quality = 0;
    _rt = 0;
    _rtmax = 0;
    _rtmin = 0;
    _scan = 0;
    _signalBaselineRatio = 0;
    _signalBaselineDifference = 0;
    _symmetry = 0;
    _width = 0;
    _gaussFitSigma = 10;
    _gaussFitR2 = FLT_MAX;
    _peakRank = INT_MAX;
}

Peak::Peak(EIC* e, int p)
{
    _pos = p;
    _eic = e;
    _sample = NULL;
    _baseMz = 0;
    _fromBlankSample = false;
    _groupOverlap = 0;
    _groupNum = 0;
    _groupOverlapFrac = 0;
    _label = 0;
    _localMaxFlag = 0;
    _maxscan = 0;
    _medianMz = 0;
    _minpos = _pos;
    _maxpos = _pos;
    _splineminpos = _pos;
    _splinemaxpos = _pos;
    _minscan = 0;
    _mzmax = 0;
    _mzmin = 0;
    _noNoiseObs = 0;
    _noNoiseFraction = 0;
    _peakArea = 0;
    _peakSplineArea = 0;
    _peakAreaCorrected = 0;
    _peakAreaFractional = 0;
    _peakAreaTop = 0;
    _peakAreaTopCorrected = 0;
    _peakIntensity = 0;
    _peakMz = 0;
    _peakBaseLineLevel = 0;
    _quality = 0;
    _rt = 0;
    _rtmax = 0;
    _rtmin = 0;
    _scan = 0;
    _signalBaselineRatio = 0;
    _signalBaselineDifference = 0;
    _symmetry = 0;
    _width = 0;
    _gaussFitSigma = 10;
    _gaussFitR2 = FLT_MAX;
    _peakRank = INT_MAX;
    if (_sample == NULL && _eic != NULL)
        _sample = _eic->sample;
}

void Peak::copyObj(const Peak& o)
{
    _pos = o._pos;
    _eic = o._eic;
    _sample = o._sample;
    _baseMz = o._baseMz;
    _fromBlankSample = o._fromBlankSample;
    _groupOverlap = o._groupOverlap;
    _groupOverlapFrac = o._groupOverlapFrac;
    _groupNum = o._groupNum;
    _label = o._label;
    _localMaxFlag = o._localMaxFlag;
    _maxscan = o._maxscan;
    _medianMz = o._medianMz;
    _minpos = o._minpos;
    _maxpos = o._maxpos;
    _splineminpos = o._splineminpos;
    _splinemaxpos = o._splinemaxpos;
    _minscan = o._minscan;
    _mzmax = o._mzmax;
    _mzmin = o._mzmin;
    _noNoiseObs = o._noNoiseObs;
    _noNoiseFraction = o._noNoiseFraction;
    _peakArea = o._peakArea;
    _peakSplineArea = o._peakSplineArea;
    _peakAreaCorrected = o._peakAreaCorrected;
    _peakAreaFractional = o._peakAreaFractional;
    _peakAreaTop = o._peakAreaTop;
    _peakAreaTopCorrected = o._peakAreaTopCorrected;
    _peakIntensity = o._peakIntensity;
    _peakMz = o._peakMz;
    _quality = o._quality;
    _peakBaseLineLevel = o._peakBaseLineLevel;
    _rt = o._rt;
    _rtmax = o._rtmax;
    _rtmin = o._rtmin;
    _scan = o._scan;
    _signalBaselineRatio = o._signalBaselineRatio;
    _signalBaselineDifference = o._signalBaselineDifference;
    _symmetry = o._symmetry;
    _width = o._width;
    _gaussFitSigma = o._gaussFitSigma;
    _gaussFitR2 = o._gaussFitR2;
    _peakRank = o._peakRank;
}

Peak& Peak::operator=(const Peak& o)
{
    copyObj(o);
    return *this;
}

Peak::Peak(const Peak& o)
{
    copyObj(o);
}

vector<mzLink> Peak::findCovariants()
{
    vector<mzLink> covariants;
    if (_sample == NULL)
        return covariants;

    // find scan range in which we will be checking for covariance
    vector<Scan*> scans;
    for (unsigned int i = 0; i < _sample->scans.size(); i++) {
        Scan* scan = _sample->scans[i];
        if (scan == NULL)
            continue;
        if (scan->mslevel != 1)
            continue;
        if (scan->rt < _rt - 0.1)
            continue;
        if (scan->rt > _rt + 0.1)
            break;
        scans.push_back(scan);
    }
    if (scans.size() == 0)
        return covariants;
    int scanCount = scans.size();

    // constuct Map
    map<int, vector<float>> M;
    map<int, vector<float>>::iterator itr;
    for (unsigned int i = 0; i < scans.size(); i++) {
        Scan* _scan = scans[i];
        for (unsigned int j = 0; j < _scan->nobs(); j++) {
            int rmz = int(_scan->mz[j] * 1000);
            if (M[rmz].size() == 0)
                M[rmz].resize(scanCount);
            M[rmz][i] = _scan->intensity[j];
        }
    }

    // merge ajacent slices
    int lastMz = 0;

    for (itr = M.begin(); itr != M.end(); ++itr) {
        // naman: Prefer prefix ++/-- operators for non-primitive types.
        float rmz = (*itr).first;
        if (lastMz != 0 && rmz - lastMz < 2) {
            // merge
            for (int j = 0; j < scanCount; j++) {
                M[rmz][j] = max(M[lastMz][j], M[rmz][j]);
                M[lastMz][j] = 0;
            }
        }
        lastMz = (*itr).first;
    }
    // refvector
    int refbinguess = int(_peakMz * 1000);
    int maxobs = 0;
    vector<float> yref;

    for (int i = refbinguess - 5; i < refbinguess + 5; i++) {
        if (M.count(i)) {
            vector<float> y = M[i];
            int nonzerCount = 0;
            for (int j = 0; j < scanCount; j++)
                if (y[j] != 0)
                    nonzerCount++;
            if (nonzerCount > maxobs) {
                yref = y;
                maxobs = nonzerCount;
            }
        }
    }

    if (yref.size() == 0)
        return covariants;

    for (itr = M.begin(); itr != M.end(); ++itr) {
        // naman: Prefer prefix ++/-- operators for
        //non-primitive types.
        int rmz = (*itr).first;
        vector<float> y = (*itr).second;
        // float score=matchScore(yref, y );
        float score = mzUtils::correlation(yref, y);
        if ((float)score < 0.5)
            continue;
        mzLink link;
        link.mz1 = _peakMz;
        link.mz2 = rmz / 1000.0;
        link.value1 = _pos;
        link.correlation = score;
        link.note = "Covariant";
        covariants.push_back(link);
    }
    return covariants;
}

Scan* Peak::getScan()
{
    if (_sample) {
        return _sample->getScan(_scan);
    } else {
        return NULL;
    }
}

bool Peak::compSampleName(const Peak& a, const Peak& b)
{
    return a._sample->getSampleName() < b._sample->getSampleName();
}

bool Peak::compSampleOrder(const Peak& a, const Peak& b)
{
    return a._sample->getSampleOrder() < b._sample->getSampleOrder();
}

float Peak::overlap(const Peak& a, const Peak& b)
{
    return (checkOverlap(a._rtmin, a._rtmax, b._rtmin, b._rtmax));
}

bool Peak::compMz(const Peak& a, const Peak& b)
{
    return a._peakMz < b._peakMz;
}

bool Peak::compArea(const Peak& a, const Peak& b)
{
    return b._peakAreaFractional < a._peakAreaFractional;
}

///////////////////////////Test case////////////////////////

class PeaksFixture
{
    private:
        vector<mzSample*> _samples;
        MavenParameters* _mavenparameters;
        Peak _peak;
        Databases _database;

        void _loadSamplesAndParameters(vector<mzSample*>& samplesToLoad,
                                         MavenParameters* mavenparameters)
        {

            ClassifierNeuralNet* clsf = new ClassifierNeuralNet();
            string loadmodel = "bin/default.model";
            clsf->loadModel(loadmodel);
            mavenparameters->compoundMassCutoffWindow->setMassCutoffAndType(10, "ppm");
            mavenparameters->clsf = clsf;
            mavenparameters->ionizationMode = -1;
            mavenparameters->matchRtFlag = true;
            mavenparameters->compoundRTWindow = 1;
            mavenparameters->samples = samplesToLoad;
            mavenparameters->eic_smoothingWindow = 10;
            mavenparameters->eic_smoothingAlgorithm = 1;
            mavenparameters->amuQ1 = 0.25;
            mavenparameters->amuQ3 = 0.30;
            mavenparameters->baseline_smoothingWindow = 5;
            mavenparameters->baseline_dropTopX = 80;
        }

        void _getGroupsFromProcessCompounds()
        {
            const char* loadCompoundDB = "bin/methods/KNOWNS.csv";
            _database.loadCompoundCSVFile(loadCompoundDB);
            vector<Compound*> compounds =
                _database.getCompoundsSubset("KNOWNS");

            _loadSamplesAndParameters(_samples, _mavenparameters);

            PeakDetector peakDetector;
            peakDetector.setMavenParameters(_mavenparameters);

            vector<mzSlice*> slices =
                peakDetector.processCompounds(compounds, "compounds");
            peakDetector.processSlices(slices, "compounds");
        }
    
    public:

        PeaksFixture()
        {
            mzSample* sample1 = new mzSample();
            sample1->loadSample("bin/methods/091215_120i.mzXML");
            _samples.push_back(sample1);
            _mavenparameters = new MavenParameters();

            _getGroupsFromProcessCompounds();

            auto allgroups = _mavenparameters->allgroups;
            _peak = allgroups[0].peaks[0];
        }

        Peak peaks() const
        {
            return _peak;
        }      
};

TEST_CASE_FIXTURE(PeaksFixture, "Testing peak class")
{
    Peak peak = peaks();

    SUBCASE("Testing copyConstructor"){
        Peak copyObj(peak);
        REQUIRE(copyObj.pos() == peak.pos());
        REQUIRE(copyObj.minpos() == peak.minpos());
        REQUIRE(copyObj.maxpos() == peak.maxpos());
        REQUIRE(copyObj.rt() == peak.rt());
        REQUIRE(copyObj.rtmin() == peak.rtmin());
        REQUIRE(copyObj.rtmax() == peak.rtmax());
        REQUIRE(copyObj.sample() == peak.sample());
        REQUIRE(copyObj.peakMz() == peak.peakMz());
        REQUIRE(copyObj.eic() == peak.eic());
        REQUIRE(copyObj.baseMz() == peak.baseMz());
        REQUIRE(copyObj.fromBlankSample() == peak.fromBlankSample());
        REQUIRE(copyObj.groupOverlap() == peak.groupOverlap());
        REQUIRE(copyObj.groupNum() == peak.groupNum());
        REQUIRE(copyObj.groupOverlapFrac() == peak.groupOverlapFrac());
        REQUIRE(copyObj.label() == peak.label());
        REQUIRE(copyObj.localMaxFlag() == peak.localMaxFlag());
        REQUIRE(copyObj.maxscan() == peak.maxscan());
        REQUIRE(copyObj.medianMz() == peak.medianMz());
        REQUIRE(copyObj.splineminpos() == peak.splineminpos());
        REQUIRE(copyObj.splinemaxpos() == peak.splinemaxpos());
        REQUIRE(copyObj.minscan() == peak.minscan());
        REQUIRE(copyObj.noNoiseObs() == peak.noNoiseObs());
        REQUIRE(copyObj.noNoiseFraction() == peak.noNoiseFraction());
        REQUIRE(copyObj.peakArea() == peak.peakArea());
        REQUIRE(copyObj.peakSplineArea() == peak.peakSplineArea());
        REQUIRE(copyObj.peakAreaCorrected() == peak.peakAreaCorrected());
        REQUIRE(copyObj.peakAreaFractional() == peak.peakAreaFractional());
        REQUIRE(copyObj.peakAreaTop() == peak.peakAreaTop());
        REQUIRE(copyObj.peakAreaTopCorrected() == peak.peakAreaTopCorrected());
        REQUIRE(copyObj.peakIntensity() == peak.peakIntensity());
        REQUIRE(copyObj.peakMz() == peak.peakMz());
        REQUIRE(copyObj.peakBaseLineLevel() == peak.peakBaseLineLevel());
        REQUIRE(copyObj.quality() == peak.quality());
        REQUIRE(copyObj.scan() == peak.scan());
        REQUIRE(copyObj.signalBaselineRatio() == peak.signalBaselineRatio());
        REQUIRE(copyObj.signalBaselineDifference() == peak.signalBaselineDifference());
        REQUIRE(copyObj.symmetry() == peak.symmetry());
        REQUIRE(copyObj.width() == peak.width());
        REQUIRE(copyObj.gaussFitSigma() == peak.gaussFitSigma());
        REQUIRE(copyObj.gaussFitR2() == peak.gaussFitR2());
        REQUIRE(copyObj.peakRank() == peak.peakRank());
    }

    SUBCASE("Testing Operator Overloading")
    {
        Peak copyObj;
        copyObj = peak;
        REQUIRE(copyObj.pos() == peak.pos());
        REQUIRE(copyObj.minpos() == peak.minpos());
        REQUIRE(copyObj.maxpos() == peak.maxpos());
        REQUIRE(copyObj.rt() == peak.rt());
        REQUIRE(copyObj.rtmin() == peak.rtmin());
        REQUIRE(copyObj.rtmax() == peak.rtmax());
        REQUIRE(copyObj.sample() == peak.sample());
        REQUIRE(copyObj.peakMz() == peak.peakMz());
        REQUIRE(copyObj.eic() == peak.eic());
        REQUIRE(copyObj.baseMz() == peak.baseMz());
        REQUIRE(copyObj.fromBlankSample() == peak.fromBlankSample());
        REQUIRE(copyObj.groupOverlap() == peak.groupOverlap());
        REQUIRE(copyObj.groupNum() == peak.groupNum());
        REQUIRE(copyObj.groupOverlapFrac() == peak.groupOverlapFrac());
        REQUIRE(copyObj.label() == peak.label());
        REQUIRE(copyObj.localMaxFlag() == peak.localMaxFlag());
        REQUIRE(copyObj.maxscan() == peak.maxscan());
        REQUIRE(copyObj.medianMz() == peak.medianMz());
        REQUIRE(copyObj.splineminpos() == peak.splineminpos());
        REQUIRE(copyObj.splinemaxpos() == peak.splinemaxpos());
        REQUIRE(copyObj.minscan() == peak.minscan());
        REQUIRE(copyObj.noNoiseObs() == peak.noNoiseObs());
        REQUIRE(copyObj.noNoiseFraction() == peak.noNoiseFraction());
        REQUIRE(copyObj.peakArea() == peak.peakArea());
        REQUIRE(copyObj.peakSplineArea() == peak.peakSplineArea());
        REQUIRE(copyObj.peakAreaCorrected() == peak.peakAreaCorrected());
        REQUIRE(copyObj.peakAreaFractional() == peak.peakAreaFractional());
        REQUIRE(copyObj.peakAreaTop() == peak.peakAreaTop());
        REQUIRE(copyObj.peakAreaTopCorrected() == peak.peakAreaTopCorrected());
        REQUIRE(copyObj.peakIntensity() == peak.peakIntensity());
        REQUIRE(copyObj.peakMz() == peak.peakMz());
        REQUIRE(copyObj.peakBaseLineLevel() == peak.peakBaseLineLevel());
        REQUIRE(copyObj.quality() == peak.quality());
        REQUIRE(copyObj.scan() == peak.scan());
        REQUIRE(copyObj.signalBaselineRatio() == peak.signalBaselineRatio());
        REQUIRE(copyObj.signalBaselineDifference() == peak.signalBaselineDifference());
        REQUIRE(copyObj.symmetry() == peak.symmetry());
        REQUIRE(copyObj.width() == peak.width());
        REQUIRE(copyObj.gaussFitSigma() == peak.gaussFitSigma());
        REQUIRE(copyObj.gaussFitR2() == peak.gaussFitR2());
        REQUIRE(copyObj.peakRank() == peak.peakRank());
    }

    SUBCASE("Testing getScan method")
    {
        Scan* scan = peak.getScan();
        REQUIRE(scan->polarity == -1);
        REQUIRE(scan->mslevel == 1);
        REQUIRE(doctest::Approx(scan->rt) == 8.57253);
        REQUIRE(scan->mz.size() == 843);
        REQUIRE(scan->intensity.size() == 843);
        REQUIRE(scan->precursorMz == 0);
        REQUIRE(doctest::Approx(scan->productMz) == 141.017);
        REQUIRE(scan->totalIntensity() == 3903634);
        REQUIRE(doctest::Approx(scan->mz[0]) == 85.0294);
        REQUIRE(doctest::Approx(scan->mz[1]) == 86.0611);
        REQUIRE(doctest::Approx(scan->intensity[0]) == 1546.34);
        REQUIRE(doctest::Approx(scan->intensity[1]) == 302.96);
    }

    SUBCASE("Testing compare Methods")
    {
        Peak p1;
        mzSample* sample = new mzSample();
        sample->loadSample("bin/methods/091215_120M.mzXML");
        p1.setPos(30);
        p1.setMinpos(20);
        p1.setMaxpos(40);
        p1.setRt(2.2312);
        p1.setRtmin(1.2233);
        p1.setRtmax(3.234);
        p1.setSample(sample);
        p1.setPeakMz(42.33);
        p1.setPeakIntensity(4.234);
        p1.setPeakAreaFractional(20.453);

        REQUIRE(peak.compSampleName(p1, peak) == true);
        REQUIRE(peak.compSampleOrder(peak, p1) == false);
        REQUIRE(peak.overlap(peak, p1) == false);
        REQUIRE(peak.compMz(p1, peak) == true);
        REQUIRE(peak.compArea(p1, peak) == true);
        REQUIRE(peak.compRtMin(p1, peak) == true);
        REQUIRE(peak.compRt(p1, peak) == true);
        REQUIRE(peak.compIntensity(p1, peak) == false);
    }

    SUBCASE("Testing Covariants")
    {
        peak.getScan();
        vector<mzLink> mzlink;
        mzlink = peak.findCovariants();
        REQUIRE(mzlink.size() == 218);
    }
}
