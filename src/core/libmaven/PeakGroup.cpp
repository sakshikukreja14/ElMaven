#include "doctest.h"
#include "PeakGroup.h"
#include "Compound.h"
#include "EIC.h"
#include "Scan.h"
#include "datastructures/adduct.h"
#include "datastructures/mzSlice.h"
#include "mzMassCalculator.h"
#include "mzSample.h"
#include "mavenparameters.h"
#include "classifierNeuralNet.h"
#include "databases.h"
#include "PeakDetector.h"
#include "masscutofftype.h"
#include "PeakGroup.h"

PeakGroup::PeakGroup()
{
    groupId = 0;
    metaGroupId = 0;
    clusterId = 0;
    groupRank = INT_MAX;

    maxIntensity = 0;
    maxAreaTopIntensity = 0;
    maxAreaIntensity = 0;
    maxHeightIntensity = 0;
    maxAreaNotCorrectedIntensity = 0;
    maxAreaTopNotCorrectedIntensity = 0;

    currentIntensity = 0;
    meanRt = 0;
    meanMz = 0;
    expectedMz = 0;

    ms2EventCount = 0;

    blankMax = 0;
    blankSampleCount = 0;
    blankMean = 0;

    sampleMax = 0;
    sampleCount = 0;
    sampleMean = 0;

    deletedFlag = false;

    totalSampleCount = 0;
    maxNoNoiseObs = 0;
    maxPeakFracionalArea = 0;
    maxSignalBaseRatio = 0;
    maxSignalBaselineRatio = 0;
    maxPeakOverlap = 0;
    maxQuality = 0;
    avgPeakQuality = 0;
    groupQuality = 0;
    weightedAvgPeakQuality = 0;
    predictedLabel = 0;
    minQuality = 0.2;
    minIntensity = 0;

    expectedAbundance = 0;
    isotopeC13count = 0;

    minRt = 0;
    maxRt = 0;

    minMz = 0;
    maxMz = 0;

    parent = NULL;
    parentIon = nullptr;

    _adduct = NULL;

    isFocused = false;
    label = 0;

    goodPeakCount = 0;
    _type = GroupType::None;
    _sliceSet = false;

    changePValue = 0;
    changeFoldRatio = 0;
    peaks.resize(0);
}

void PeakGroup::copyObj(const PeakGroup& o)
{
    groupId = o.groupId;
    metaGroupId = o.metaGroupId;
    clusterId = o.clusterId;
    groupRank = o.groupRank;

    minQuality = o.minQuality;
    minIntensity = o.minIntensity;
    maxIntensity = o.maxIntensity;
    maxAreaTopIntensity = o.maxAreaTopIntensity;
    maxAreaIntensity = o.maxAreaIntensity;
    maxHeightIntensity = o.maxHeightIntensity;
    maxAreaNotCorrectedIntensity = o.maxAreaNotCorrectedIntensity;
    maxAreaTopNotCorrectedIntensity = o.maxAreaTopNotCorrectedIntensity;
    currentIntensity = o.currentIntensity;
    meanRt = o.meanRt;
    meanMz = o.meanMz;
    expectedMz = o.expectedMz;

    ms2EventCount = o.ms2EventCount;
    fragMatchScore = o.fragMatchScore;
    fragmentationPattern = o.fragmentationPattern;
    _adduct = o.getAdduct();

    blankMax = o.blankMax;
    blankSampleCount = o.blankSampleCount;
    blankMean = o.blankMean;

    sampleMax = o.sampleMax;
    sampleCount = o.sampleCount;
    sampleMean = o.sampleMean;

    totalSampleCount = o.totalSampleCount;
    maxNoNoiseObs = o.maxNoNoiseObs;
    maxPeakFracionalArea = o.maxPeakFracionalArea;
    maxSignalBaseRatio = o.maxSignalBaseRatio;
    maxSignalBaselineRatio = o.maxSignalBaselineRatio;
    maxPeakOverlap = o.maxPeakOverlap;
    maxQuality = o.maxQuality;
    avgPeakQuality = o.avgPeakQuality;
    groupQuality = o.groupQuality;
    weightedAvgPeakQuality = o.weightedAvgPeakQuality;
    predictedLabel = o.predictedLabel;
    expectedAbundance = o.expectedAbundance;
    isotopeC13count = o.isotopeC13count;

    deletedFlag = o.deletedFlag;

    minRt = o.minRt;
    maxRt = o.maxRt;

    minMz = o.minMz;
    maxMz = o.maxMz;

    parent = o.parent;
    parentIon = o.parentIon;
    setSlice(o.getSlice());

    srmId = o.srmId;
    isFocused = o.isFocused;
    label = o.label;

    goodPeakCount = o.goodPeakCount;
    _type = o._type;
    _sliceSet = o.hasSlice();
    tagString = o.tagString;

    changeFoldRatio = o.changeFoldRatio;
    changePValue = o.changePValue;
    peaks = o.peaks;
    samples = o.samples;

    markedBadByCloudModel = o.markedBadByCloudModel;
    markedGoodByCloudModel = o.markedGoodByCloudModel;

    copyChildren(o);
}

PeakGroup::~PeakGroup()
{
    clear();
}

void PeakGroup::copyChildren(const PeakGroup& o)
{
    children = o.children;
    childrenBarPlot = o.childrenBarPlot;
    for (unsigned int i = 0; i < children.size(); i++)
        children[i].parent = this;
    for (unsigned int i = 0; i < childrenBarPlot.size(); i++)
        childrenBarPlot[i].parent = this;

    childAdducts = o.childAdducts;
    for (auto& adductGroup : childAdducts)
        adductGroup.parent = this;
}

void PeakGroup::clear()
{
    deletePeaks();
    deleteChildren();
    meanMz = 0;
    expectedMz = 0;
    groupRank = INT_MAX;
}

bool PeakGroup::isMS1()
{
    if (peaks.size() == 0)
        return false;
    Peak peak = peaks[0];
    if (peak.getSample()) {
        Scan* scan = peak.getSample()->getScan(peak.scan);
        if (scan && scan->mslevel == 1)
            return true;
    }
    return false;
}

bool PeakGroup::hasCompoundLink() const
{
    if (hasSlice() && _slice.compound != NULL)
        return true;
    return false;
}

Compound* PeakGroup::getCompound() const
{
    if (hasSlice() && _slice.compound != NULL) {
        return _slice.compound;
    }
    return NULL;
}

void PeakGroup::setCompound(Compound* compound)
{
    _slice.compound = compound;
    _sliceSet = true;
}

void PeakGroup::addPeak(const Peak& peak)
{
    peaks.push_back(peak);
    peaks.back().groupNum = groupId;
}

void PeakGroup::setSlice(const mzSlice& slice)
{
    _slice = slice;
    _sliceSet = true;
}

const mzSlice& PeakGroup::getSlice() const
{
    return _slice;
}

bool PeakGroup::hasSlice() const
{
    return _sliceSet;
}

bool PeakGroup::sliceIsZero() const
{
    if (((mzUtils::almostEqual(_slice.mzmin, 0.0f)
          && mzUtils::almostEqual(_slice.mzmax, 0.0f)))
        || (mzUtils::almostEqual(_slice.rtmin, 0.0f)
            && mzUtils::almostEqual(_slice.mzmin, 0.0f))) {
        return true;
    }
    return false;
}

void PeakGroup::deletePeaks()
{
    peaks.clear();
}

float PeakGroup::medianRt()
{
    float* rts = new float[peaks.size()];
    for (unsigned int i = 0; i < peakCount(); i++)
        rts[i] = peaks[i].rt;
    float medianValue = mzUtils::median(rts, peakCount());
    delete[] rts;
    return medianValue;
}

float PeakGroup::expectedRtDiff()
{
    auto associatedCompound = getCompound();
    if (associatedCompound != nullptr && associatedCompound->expectedRt > 0) {
        return abs(associatedCompound->expectedRt - meanRt);
    }
    return -1.0f;
}

void PeakGroup::deleteChildren()
{
    children.clear();
}

bool PeakGroup::deleteChild(PeakGroup* child)
{
    if (!child)
        return false;

    auto preDeletionChildCount = children.size();
    children.erase(remove_if(begin(children),
                             end(children),
                             [&](PeakGroup& group) { return child == &group; }),
                   children.end());

    // child was found and removed
    if (children.size() != preDeletionChildCount)
        return true;

    return false;
}

vector<float> PeakGroup::getOrderedIntensityVector(vector<mzSample*>& samples,
                                                   QType type)
{
    if (samples.size() == 0) {
        vector<float> x;
        return x;
    }  // empty vector;

    map<mzSample*, float> sampleOrder;
    vector<float> maxIntensity(samples.size(), 0);

    for (unsigned int j = 0; j < samples.size(); j++) {
        sampleOrder[samples[j]] = j;
        maxIntensity[j] = 0;
    }

    for (unsigned int j = 0; j < peaks.size(); j++) {
        Peak& peak = peaks.at(j);
        mzSample* sample = peak.getSample();

        if (sampleOrder.count(sample) > 0) {
            int s = sampleOrder[sample];
            float y = 0;
            switch (type) {
            case AreaTop:
                y = peak.peakAreaTopCorrected;
                break;
            case Area:
                y = peak.peakAreaCorrected;
                break;
            case Height:
                y = peak.peakIntensity;
                break;
            case AreaNotCorrected:
                y = peak.peakArea;
                break;
            case AreaTopNotCorrected:
                y = peak.peakAreaTop;
                break;
            case RetentionTime:
                y = peak.rt;
                break;
            case Quality:
                y = peak.quality;
                break;
            case SNRatio:
                y = peak.signalBaselineRatio;
                break;
            default:
                y = peak.peakIntensity;
                break;
            }

            // normalize
            if (sample)
                y *= sample->getNormalizationConstant();
            if (maxIntensity[s] < y) {
                maxIntensity[s] = y;
            }
        }
    }
    return maxIntensity;
}

void PeakGroup::computeAvgBlankArea(const vector<EIC*>& eics)
{
    if (peaks.size() == 0)
        return;

    // find range to fill in
    float rtmin = peaks[0].rtmin;
    float rtmax = peaks[0].rtmax;

    for (unsigned int i = 1; i < peaks.size(); i++) {
        if (peaks[i].rtmin < rtmin)
            rtmin = peaks[i].rtmin;
        if (peaks[i].rtmax > rtmax)
            rtmax = peaks[i].rtmax;
    }
    rtmin = rtmin - 0.25;
    rtmax = rtmax + 0.25;
    float sum = 0;
    int len = 0;
    for (unsigned int i = 0; i < eics.size(); i++) {
        EIC* eic = eics[i];
        if (eic->sample != NULL && eic->sample->isBlank == false)
            continue;
        for (unsigned int pos = 0; pos < eic->intensity.size(); pos++) {
            if (eic->rt[pos] >= rtmin && eic->rt[pos] <= rtmax
                && eic->intensity[pos] > 0) {
                sum += eic->intensity[pos];
                len++;
            }
        }
    }
    this->blankMean = 0;  // default zero
    if (len > 0)
        this->blankMean = (float)sum / len;
}

void PeakGroup::fillInPeaks(const vector<EIC*>& eics)
{
    if (peaks.size() == eics.size())
        return;
    if (peaks.size() == 0)
        return;

    // find range to fill in
    float rtmin = peaks[0].rtmin;
    float rtmax = peaks[0].rtmax;

    for (unsigned int i = 1; i < peaks.size(); i++) {
        if (peaks[i].rtmin < rtmin)
            rtmin = peaks[i].rtmin;
        if (peaks[i].rtmax > rtmax)
            rtmax = peaks[i].rtmax;
    }

    int filledInCount = 0;

    for (unsigned int i = 0; i < eics.size(); i++) {
        EIC* eic = eics[i];
        if (eic == NULL)
            continue;
        if (eic->spline == NULL)
            continue;
        if (eic->intensity.size() == 0)
            continue;

        bool missing = true;

        for (unsigned int j = 0; j < peaks.size(); j++) {
            if (peaks[j].getEIC() == eic) {
                missing = false;
                break;
            }
        }

        if (missing) {
            int maxpos = 0;
            for (unsigned int pos = 1; pos < eic->intensity.size() - 1; pos++) {
                if (eic != NULL && eic->intensity[pos] != 0 && eic->mz[pos] != 0
                    && eic->rt[pos] >= rtmin && eic->rt[pos] <= rtmax
                    && eic->spline[pos] > eic->spline[pos - 1]
                    && eic->spline[pos] > eic->spline[pos + 1]) {
                    if (maxpos != 0
                        && eic->intensity[pos] > eic->intensity[maxpos]) {
                        maxpos = pos;
                    } else {
                        maxpos = pos;
                    }
                }
            }

            if (maxpos != 0 && eic->intensity[maxpos] != 0) {
                Peak peak(eic, maxpos);
                eic->findPeakBounds(peak);
                eic->getPeakDetails(peak);
                this->addPeak(peak);
                filledInCount++;
            }
        }
    }
}

void PeakGroup::reduce()
{
    map<mzSample*, Peak> maxPeaks;
    map<mzSample*, Peak>::iterator itr;
    if (peaks.size() < 2)
        return;

    float groupMeanRt = 0;
    float totalWeight = 1;

    for (unsigned int i = 0; i < peaks.size(); i++) {
        totalWeight += peaks[i].peakIntensity;
    }
    for (unsigned int i = 0; i < peaks.size(); i++) {
        groupMeanRt += peaks[i].rt * peaks[i].peakIntensity / totalWeight;
    }

    // In each group, take peak that closest to the mean retention time of a
    // group
    for (unsigned int i = 0; i < peaks.size(); i++) {
        mzSample* c = peaks[i].getSample();
        // In each group, take the heighest peak
        if (maxPeaks.count(c) == 0
            || maxPeaks[c].peakIntensity < peaks[i].peakIntensity) {
            maxPeaks[c].copyObj(peaks[i]);
        }
    }
    peaks.clear();
    for (itr = maxPeaks.begin(); itr != maxPeaks.end(); ++itr) {
        const Peak& peak = (*itr).second;
        addPeak(peak);
    }
}

void PeakGroup::setLabel(char label)
{
    this->label = label;

    if (parent != nullptr && tagString == "C12 PARENT") {
        parent->setLabel(label);
        return;
    }

    for (auto& child : children) {
        if (child.tagString == "C12 PARENT" && child.label != label)
            child.setLabel(label);
    }
}

float PeakGroup::massCutoffDist(float cmass, MassCutoff* massCutoff)
{
    return mzUtils::massCutoffDist(cmass, meanMz, massCutoff);
}

void PeakGroup::updateQuality()
{
    maxQuality = 0;
    goodPeakCount = 0;

    float peakQualitySum = 0;
    float weightedSum = 0;
    float sumWeights = 0;
    for (const auto peak : peaks) {
        if (peak.quality > maxQuality)
            maxQuality = peak.quality;
        if (peak.quality > minQuality)
            goodPeakCount++;  // Sabu
        peakQualitySum += peak.quality;
        weightedSum += peak.quality * peak.peakIntensity;
        sumWeights += peak.peakIntensity;
    }
    avgPeakQuality = peakQualitySum / peaks.size();
    weightedAvgPeakQuality = weightedSum / sumWeights;
}

// TODO: Remove this function as expected mz should be calculated while creating
// the group - Sahil
double PeakGroup::getExpectedMz(int charge)
{
    float mz = 0;

    if (isIsotope() && childCount() == 0 && hasSlice()
        && _slice.compound != NULL && !_slice.compound->formula().empty()
        && _slice.compound->mass > 0) {
        return expectedMz;
    } else if (!isIsotope() && hasSlice() && _slice.compound != NULL
               && _slice.compound->mass > 0) {
        if (!_slice.compound->formula().empty() && _adduct != nullptr) {
            auto mass =
                MassCalculator::computeNeutralMass(_slice.compound->formula());
            mz = _adduct->computeAdductMz(mass);
        } else if (!_slice.compound->formula().empty()
                   || _slice.compound->neutralMass != 0.0f) {
            mz = _slice.compound->adjustedMass(charge);
        } else {
            mz = _slice.compound->mass;
        }
        return mz;
    } else if (hasSlice() && _slice.compound != NULL
               && _slice.compound->mass == 0
               && _slice.compound->productMz > 0) {
        mz = _slice.compound->productMz;
        return mz;
    }

    return -1;
}

void PeakGroup::groupStatistics()
{
    float rtSum = 0;
    float mzSum = 0;
    maxIntensity = 0;
    maxAreaTopIntensity = 0;
    maxAreaIntensity = 0;
    maxHeightIntensity = 0;
    maxAreaNotCorrectedIntensity = 0;
    maxAreaTopNotCorrectedIntensity = 0;
    currentIntensity = 0;
    totalSampleCount = 0;

    blankMax = 0;
    blankSampleCount = 0;

    sampleMax = 0;
    sampleCount = 0;
    sampleMean = 0;

    maxNoNoiseObs = 0;
    minRt = 0;
    maxRt = 0;
    minMz = 0;
    maxMz = 0;

    maxPeakFracionalArea = 0;
    maxQuality = 0;
    avgPeakQuality = 0;
    groupQuality = 0;
    weightedAvgPeakQuality = 0;
    predictedLabel = 0;
    goodPeakCount = 0;
    maxSignalBaselineRatio = 0;
    int nonZeroCount = 0;
    //@Kailash: Added for Avg Peak Quality and Intensity Weighted Peak Quality
    float peakQualitySum = 0;
    float highestIntensity = 0;
    float weightedSum = 0;
    float sumWeights = 0;

    for (unsigned int i = 0; i < peaks.size(); i++) {
        if (peaks[i].pos != 0 && peaks[i].baseMz != 0) {
            rtSum += peaks[i].rt;
            mzSum += peaks[i].baseMz;
            nonZeroCount++;
        }
        if (peaks[i].peakIntensity > 0)
            totalSampleCount++;

        float max;
        switch (quantitationType) {
        case AreaTop:
            max = peaks[i].peakAreaTopCorrected;
            break;
        case Area:
            max = peaks[i].peakAreaCorrected;
            break;
        case Height:
            max = peaks[i].peakIntensity;
            break;
        case AreaTopNotCorrected:
            max = peaks[i].peakAreaTop;
            break;
        case AreaNotCorrected:
            max = peaks[i].peakArea;
            break;
        default:
            max = peaks[i].peakIntensity;
            break;
        }

        if (peaks[i].peakAreaTopCorrected > maxAreaTopIntensity)
            maxAreaTopIntensity = peaks[i].peakAreaTopCorrected;
        if (peaks[i].peakAreaCorrected > maxAreaIntensity)
            maxAreaIntensity = peaks[i].peakAreaCorrected;
        if (peaks[i].peakIntensity > maxHeightIntensity)
            maxHeightIntensity = peaks[i].peakIntensity;
        if (peaks[i].peakArea > maxAreaNotCorrectedIntensity)
            maxAreaNotCorrectedIntensity = peaks[i].peakArea;
        if (peaks[i].peakAreaTop > maxAreaTopNotCorrectedIntensity)
            maxAreaTopNotCorrectedIntensity = peaks[i].peakAreaTop;

        if (max > maxIntensity) {
            maxIntensity = max;
            currentIntensity = max;
            meanMz = peaks[i].baseMz;
            meanRt = peaks[i].rt;
        }

        if (peaks[i].noNoiseObs > maxNoNoiseObs)
            maxNoNoiseObs = peaks[i].noNoiseObs;
        if (minRt == 0 || peaks[i].rtmin < minRt)
            minRt = peaks[i].rtmin;
        if (maxRt == 0 || peaks[i].rtmax > maxRt)
            maxRt = peaks[i].rtmax;
        if (minMz == 0 || peaks[i].mzmin < minMz)
            minMz = peaks[i].mzmin;
        if (maxMz == 0 || peaks[i].mzmax > maxMz)
            maxMz = peaks[i].mzmax;
        if (peaks[i].peakAreaFractional > maxPeakFracionalArea)
            maxPeakFracionalArea = peaks[i].peakAreaFractional;
        if (peaks[i].quality > maxQuality)
            maxQuality = peaks[i].quality;
        if (peaks[i].quality > minQuality)
            goodPeakCount++;  // Sabu
        if (peaks[i].signalBaselineRatio > maxSignalBaselineRatio)
            maxSignalBaselineRatio = peaks[i].signalBaselineRatio;

        if (peaks[i].fromBlankSample) {
            blankSampleCount++;
            if (peaks[i].peakIntensity > blankMax)
                blankMax = peaks[i].peakIntensity;
        } else {
            sampleMean += peaks[i].peakIntensity;
            sampleCount++;
            if (peaks[i].peakIntensity > sampleMax)
                sampleMax = peaks[i].peakIntensity;
        }

        weightedSum += peaks[i].quality * peaks[i].peakIntensity;
        sumWeights += peaks[i].peakIntensity;
        peakQualitySum += peaks[i].quality;
        if (peaks[i].peakIntensity > highestIntensity)
            highestIntensity = peaks[i].peakIntensity;
    }
    avgPeakQuality = peakQualitySum / peaks.size();
    weightedAvgPeakQuality = weightedSum / sumWeights;

    if (sampleCount > 0)
        sampleMean = sampleMean / sampleCount;
    if (nonZeroCount) {
        meanRt = rtSum / nonZeroCount;
        meanMz = mzSum / nonZeroCount;
    }
    groupOverlapMatrix();
}

void PeakGroup::groupOverlapMatrix()
{
    for (unsigned int i = 0; i < peaks.size(); i++)
        peaks[i].groupOverlapFrac = 0;

    for (unsigned int i = 0; i < peaks.size(); i++) {
        Peak& a = peaks[i];
        for (unsigned int j = i; j < peaks.size(); j++) {
            Peak& b = peaks[j];
            float overlap = checkOverlap(
                a.rtmin, a.rtmax, b.rtmin, b.rtmax);  // check for overlap
            if (overlap > 0) {
                b.groupOverlapFrac += log(overlap);
                a.groupOverlapFrac += log(overlap);
            }
        }
    }
    // normalize
    for (unsigned int i = 0; i < peaks.size(); i++)
        peaks[i].groupOverlapFrac /= peaks.size();
}

void PeakGroup::summary()
{
    cerr << tagString << endl;
    cerr << "\t"
         << "meanRt=" << meanRt << endl
         << "\t"
         << "meanMz=" << meanMz << endl
         << "\t"
         << "expectedMz=" << expectedMz << endl
         << "\t"
         << "goodPeakCount=" << goodPeakCount << endl
         << "\t"
         << "maxQuality=" << maxQuality << endl
         << "\t"
         << "maxNoNoiseObs=" << maxNoNoiseObs << endl
         << "\t"
         << "sampleCount=" << sampleCount << endl
         << "\t"
         << "maxSignalBaselineRatio=" << maxSignalBaselineRatio << endl
         << "\t"
         << "maxPeakFracionalArea=" << maxPeakFracionalArea << endl
         << "\t"
         << "blankMean=" << blankMean << endl
         << "\t"
         << "sampleMean=" << sampleMean << endl
         << "\t"
         << "maxIntensity=" << maxIntensity << endl
         << endl;

    for (unsigned int i = 0; i < peaks.size(); i++) {
        cerr << "\t\t"
             << "Q:" << peaks[i].quality << " "
             << "pAf:" << peaks[i].peakAreaFractional << " "
             << "noNf" << peaks[i].noNoiseFraction << " "
             << "noObs:" << peaks[i].noNoiseObs << " "
             << "w:" << peaks[i].width << " "
             << "sn:" << peaks[i].signalBaselineRatio << " "
             << "ovp:" << peaks[i].groupOverlapFrac << endl;
    }

    for (unsigned int i = 0; i < children.size(); i++)
        children[i].summary();
}

PeakGroup::PeakGroup(const PeakGroup& o)
{
    copyObj(o);
}

PeakGroup& PeakGroup::operator=(const PeakGroup& o)
{
    copyObj(o);
    return *this;
}

bool PeakGroup::operator==(const PeakGroup* o)
{
    if (this == o) {
        cerr << o << " " << this << endl;
        return true;
    }
    return false;
}

Peak* PeakGroup::getPeak(mzSample* s)
{
    if (s == NULL)
        return NULL;
    for (unsigned int i = 0; i < peaks.size(); i++) {
        if (peaks[i].getSample() == s) {
            return &peaks[i];
        }
    }
    return NULL;
}

void PeakGroup::reorderSamples()
{
    std::sort(peaks.begin(), peaks.end(), Peak::compIntensity);
    for (unsigned int i = 0; i < peaks.size(); i++) {
        mzSample* s = peaks[i].getSample();
        if (s != NULL)
            s->setSampleOrder(i);
    }
}

string PeakGroup::getName()
{
    string tag;
    // compound is assigned in case of targeted search
    if (hasSlice() && _slice.compound != NULL)
        tag = _slice.compound->name;
    // add name of external charged species fused with adduct
    if (_adduct != nullptr)
        tag += " | " + _adduct->getName();
    // add isotopic label
    if (!tagString.empty())
        tag += " | " + tagString;
    // add SRM ID for MS/MS data
    if (!srmId.empty())
        tag += " | " + srmId;
    // no compound in case of untargeted peak detection
    // group is referenced as MeanMz@MeanRT
    if (tag.empty() && meanMz && meanRt) {
        stringstream stream;
        stream << fixed << setprecision(6) << meanMz << "@" << setprecision(2)
               << meanRt;
        tag = stream.str();
    }
    // if all else fails, use group ID
    if (tag.empty())
        tag = integer2string(groupId);
    return tag;
}

vector<Scan*> PeakGroup::getRepresentativeFullScans()
{
    vector<Scan*> matchedscans;
    for (unsigned int i = 0; i < peaks.size(); i++) {
        mzSample* sample = peaks[i].getSample();
        if (sample == NULL)
            continue;
        Scan* scan = sample->getScan(peaks[i].scan);
        if (scan and scan->mslevel == 1)
            matchedscans.push_back(scan);
    }
    return matchedscans;
}

vector<Scan*> PeakGroup::getFragmentationEvents()
{
    vector<Scan*> matchedScans;
    if (!this->isMS1())
        return matchedScans;
    for (auto peak : peaks) {
        mzSample* sample = peak.getSample();
        if (sample == NULL)
            continue;
        mzSlice slice(minMz, maxMz, peak.rtmin, peak.rtmax);
        vector<Scan*> scans = sample->getFragmentationEvents(&slice);
        matchedScans.insert(matchedScans.end(), scans.begin(), scans.end());
    }
    return matchedScans;
}

void PeakGroup::computeFragPattern(float productPpmTolr)
{
    vector<Scan*> ms2Events = getFragmentationEvents();
    if (ms2Events.size() == 0)
        return;
    sort(ms2Events.begin(), ms2Events.end(), Scan::compIntensity);

    float minFractionalIntensity = 0.01;
    float minSignalNoiseRatio = 1;
    int maxFragmentSize = 1024;
    Fragment fragment(ms2Events[0],
                      minFractionalIntensity,
                      minSignalNoiseRatio,
                      maxFragmentSize);

    for (Scan* scan : ms2Events) {
        fragment.addBrotherFragment(new Fragment(scan,
                                                 minFractionalIntensity,
                                                 minSignalNoiseRatio,
                                                 maxFragmentSize));
    }

    fragment.buildConsensus(productPpmTolr);
    fragment.consensus->sortByMz();
    fragmentationPattern = fragment.consensus;
    ms2EventCount = ms2Events.size();
}

Scan* PeakGroup::getAverageFragmentationScan(float productPpmTolr)
{
    // build consensus ms2 specta
    computeFragPattern(productPpmTolr);
    Scan* avgScan = new Scan(NULL, 0, 2, 0, 0, 0);

    for (unsigned int i = 0; i < fragmentationPattern.mzValues.size(); i++) {
        avgScan->mz.push_back(fragmentationPattern.mzValues[i]);
        avgScan->intensity.push_back(fragmentationPattern.intensityValues[i]);
    }

    avgScan->precursorMz = meanMz;
    avgScan->rt = meanRt;

    return avgScan;
}

void PeakGroup::matchFragmentation(float ppmTolerance, string scoringAlgo)
{
    if (this->getCompound() == nullptr || this->isAdduct()
        || ms2EventCount == 0)
        return;

    fragMatchScore =
        getCompound()->scoreCompoundHit(&fragmentationPattern, ppmTolerance);
    fragMatchScore.mergedScore = fragMatchScore.getScoreByName(scoringAlgo);
}

void PeakGroup::calGroupRank(bool deltaRtCheckFlag,
                             int qualityWeight,
                             int intensityWeight,
                             int deltaRTWeight)
{
    float rtDiff = expectedRtDiff();

    // Peak Group Rank accoording to given weightage
    double A = (double)qualityWeight / 10;
    double B = (double)intensityWeight / 10;
    double C = (double)deltaRTWeight / 10;

    if (deltaRtCheckFlag && hasSlice() && _slice.compound != NULL
        && _slice.compound->expectedRt > 0) {
        groupRank = pow(rtDiff, 2 * C) * pow((1.1 - maxQuality), A)
                    * (1
                       / (pow(log(maxIntensity + 1),
                              B)));  // TODO Formula to rank groups
    } else {
        groupRank =
            pow((1.1 - maxQuality), A) * (1 / (pow(log(maxIntensity + 1), B)));
    }
}

void PeakGroup::setSelectedSamples(vector<mzSample*> vsamples)
{
    samples.clear();

    for (size_t i = 0; i < vsamples.size(); ++i) {
        if (vsamples[i]->isSelected) {
            samples.push_back(vsamples[i]);
        }
    }
}

void PeakGroup::setAdduct(Adduct* adduct)
{
    _adduct = adduct;
    if (_adduct != nullptr
        && _adduct->getName() != MassCalculator::MinusHAdduct->getName()
        && _adduct->getName() != MassCalculator::PlusHAdduct->getName()) {
        _type = GroupType::Adduct;
    }
}

Adduct* PeakGroup::getAdduct() const
{
    if (isIsotope() && parent != nullptr)
        return parent->getAdduct();
    return _adduct;
}

string PeakGroup::tableName() const
{
    return _tableName;
}

void PeakGroup::setTableName(string tableName)
{
    _tableName = tableName;
    for (auto& child : children)
        child.setTableName(tableName);
}

//////////////////////Test cases////////////////////////////
class PeakGroupFixture
{
    private:
    vector<mzSample*> _samples;
    MavenParameters* _mavenparameters;
    PeakGroup _peakgroup;
    Databases _database;
    PeakGroup _secondPeakGroup;

    void _loadSamplesAndParameters(vector<mzSample*>& samplesToLoad,
                                   MavenParameters* mavenparameters)
    {
        ClassifierNeuralNet* clsf = new ClassifierNeuralNet();
        string loadmodel = "bin/default.model";
        clsf->loadModel(loadmodel);
        mavenparameters->compoundMassCutoffWindow->setMassCutoffAndType(10,
                                                                        "ppm");
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

    void _getGroupsFromProcessCompounds(string loadCompoundDB)
    {
        _database.loadCompoundCSVFile(loadCompoundDB.c_str());
        vector<Compound*> compounds = _database.getCompoundsSubset("KNOWNS");

        _loadSamplesAndParameters(_samples, _mavenparameters);

        PeakDetector peakDetector;
        peakDetector.setMavenParameters(_mavenparameters);

        vector<mzSlice*> slices =
            peakDetector.processCompounds(compounds, "compounds");
        peakDetector.processSlices(slices, "compounds");
    }

    public:
    PeakGroupFixture()
    {
        mzSample* sample1 = new mzSample();
        sample1->loadSample("bin/methods/091215_120i.mzXML");
        _samples.push_back(sample1);
        mzSample * sample2 = new mzSample();
        sample2->loadSample("bin/methods/091215_120M.mzXML");
        _samples.push_back(sample2);
        mzSample * sample3 = new mzSample();
        sample3->loadSample("bin/methods/091215_240i.mzXML");
        _samples.push_back(sample3);
        mzSample * sample4 = new mzSample();
        sample4->loadSample("bin/methods/091215_240M.mzXML");
        _samples.push_back(sample4);
        _mavenparameters = new MavenParameters();

        _getGroupsFromProcessCompounds("bin/methods/KNOWNS.csv");

        auto allgroup = _mavenparameters->allgroups;
        _peakgroup = allgroup[0];
        _secondPeakGroup = allgroup[1];
    }

    PeakGroup peakgroup() const
    {
        return _peakgroup;
    }

    PeakGroup secondPeakGroup() const
    {
        return _secondPeakGroup;
    }

    vector<mzSample*> samples() const
    {
        return _samples;
    }

    ~PeakGroupFixture()
    {
        delete _samples[0];
        delete _samples[1];
        delete _samples[2];
        delete _samples[3];
        delete _mavenparameters;
    }
};

TEST_CASE_FIXTURE(PeakGroupFixture, "Testing PeakGroup class")
{
    PeakGroup peakGroup = peakgroup();

    SUBCASE("Testing CopyConsrtuctor")
    {
        PeakGroup o(peakGroup);

        REQUIRE(peakGroup.groupId == o.groupId);
        REQUIRE(peakGroup.metaGroupId == o.metaGroupId);
        REQUIRE(peakGroup.clusterId == o.clusterId);
        REQUIRE(peakGroup.groupRank == o.groupRank);
        REQUIRE(peakGroup.minQuality == o.minQuality);
        REQUIRE(peakGroup.minIntensity == o.minIntensity);
        REQUIRE(peakGroup.maxIntensity == o.maxIntensity);
        REQUIRE(peakGroup.maxAreaTopIntensity ==
                o.maxAreaTopIntensity);
        REQUIRE(peakGroup.maxAreaIntensity == o.maxAreaIntensity);
        REQUIRE(peakGroup.maxHeightIntensity == o.maxHeightIntensity);
        REQUIRE(peakGroup.maxAreaNotCorrectedIntensity ==
                o.maxAreaNotCorrectedIntensity);
        REQUIRE(peakGroup.maxAreaTopNotCorrectedIntensity ==
                o.maxAreaTopNotCorrectedIntensity);
        REQUIRE(peakGroup.currentIntensity == o.currentIntensity);
        REQUIRE(peakGroup.meanRt == o.meanRt);
        REQUIRE(peakGroup.meanMz == o.meanMz);
        REQUIRE(peakGroup.expectedMz == o.expectedMz);
        REQUIRE(peakGroup.ms2EventCount == o.ms2EventCount);
        REQUIRE(peakGroup.blankMax == o.blankMax);
        REQUIRE(peakGroup.blankSampleCount ==
                o.blankSampleCount);
        REQUIRE(peakGroup.blankMean == o.blankMean);
        REQUIRE(peakGroup.sampleMax == o.sampleMax);
        REQUIRE(peakGroup.sampleCount == o.sampleCount);
        REQUIRE(peakGroup.sampleMean == o.sampleMean);
        REQUIRE(peakGroup.totalSampleCount == o.totalSampleCount);
        REQUIRE(peakGroup.maxNoNoiseObs == o.maxNoNoiseObs);
        REQUIRE(peakGroup.maxPeakFracionalArea == o.maxPeakFracionalArea);
        REQUIRE(peakGroup.maxSignalBaseRatio == o.maxSignalBaseRatio);
        REQUIRE(peakGroup.maxSignalBaselineRatio == o.maxSignalBaselineRatio);
        REQUIRE(peakGroup.maxPeakOverlap == o.maxPeakOverlap);
        REQUIRE(peakGroup.maxQuality == o.maxQuality);
        REQUIRE(peakGroup.avgPeakQuality == o.avgPeakQuality);
        REQUIRE(peakGroup.groupQuality == o.groupQuality);
        REQUIRE(peakGroup.weightedAvgPeakQuality == o.weightedAvgPeakQuality);
        REQUIRE(peakGroup.predictedLabel == o.predictedLabel);
        REQUIRE(peakGroup.expectedAbundance == o.expectedAbundance);
        REQUIRE(peakGroup.isotopeC13count == o.isotopeC13count);
        REQUIRE(peakGroup.deletedFlag == o.deletedFlag);
        REQUIRE(peakGroup.minRt == o.minRt);
        REQUIRE(peakGroup.maxRt == o.maxRt);
        REQUIRE(peakGroup.minMz == o.minMz);
        REQUIRE(peakGroup.maxMz == o.maxMz);
        REQUIRE(peakGroup.parent == o.parent);
        REQUIRE(peakGroup.parentIon == o.parentIon);
        REQUIRE(peakGroup.srmId == o.srmId);
        REQUIRE(peakGroup.isFocused == o.isFocused);
        REQUIRE(peakGroup.label == o.label);
        REQUIRE(peakGroup.goodPeakCount ==
                o.goodPeakCount);
        REQUIRE(peakGroup._type == o._type);
        REQUIRE(peakGroup.hasSlice() == o.hasSlice());
        REQUIRE(peakGroup.tagString == o.tagString);
        REQUIRE(peakGroup.changeFoldRatio ==
                o.changeFoldRatio);
        REQUIRE(peakGroup.changePValue ==
                o.changePValue);
        REQUIRE(peakGroup.peaks.size() == o.peaks.size());
        REQUIRE(peakGroup.samples.size() == o.samples.size());
        REQUIRE(peakGroup.markedBadByCloudModel ==
                o.markedBadByCloudModel);
        REQUIRE(peakGroup.markedGoodByCloudModel ==
                o.markedGoodByCloudModel);
    }

    SUBCASE("Testing '=' operator overloading")
    {
        PeakGroup o;
        o = peakGroup;
        REQUIRE(peakGroup.groupId == o.groupId);
        REQUIRE(peakGroup.metaGroupId == o.metaGroupId);
        REQUIRE(peakGroup.clusterId == o.clusterId);
        REQUIRE(peakGroup.groupRank == o.groupRank);
        REQUIRE(peakGroup.minQuality == o.minQuality);
        REQUIRE(peakGroup.minIntensity == o.minIntensity);
        REQUIRE(peakGroup.maxIntensity == o.maxIntensity);
        REQUIRE(peakGroup.maxAreaTopIntensity ==
                o.maxAreaTopIntensity);
        REQUIRE(peakGroup.maxAreaIntensity == o.maxAreaIntensity);
        REQUIRE(peakGroup.maxHeightIntensity == o.maxHeightIntensity);
        REQUIRE(peakGroup.maxAreaNotCorrectedIntensity ==
                o.maxAreaNotCorrectedIntensity);
        REQUIRE(peakGroup.maxAreaTopNotCorrectedIntensity ==
                o.maxAreaTopNotCorrectedIntensity);
        REQUIRE(peakGroup.currentIntensity == o.currentIntensity);
        REQUIRE(peakGroup.meanRt == o.meanRt);
        REQUIRE(peakGroup.meanMz == o.meanMz);
        REQUIRE(peakGroup.expectedMz == o.expectedMz);
        REQUIRE(peakGroup.ms2EventCount == o.ms2EventCount);
        REQUIRE(peakGroup.blankMax == o.blankMax);
        REQUIRE(peakGroup.blankSampleCount ==
                o.blankSampleCount);
        REQUIRE(peakGroup.blankMean == o.blankMean);
        REQUIRE(peakGroup.sampleMax == o.sampleMax);
        REQUIRE(peakGroup.sampleCount == o.sampleCount);
        REQUIRE(peakGroup.sampleMean == o.sampleMean);
        REQUIRE(peakGroup.totalSampleCount == o.totalSampleCount);
        REQUIRE(peakGroup.maxNoNoiseObs == o.maxNoNoiseObs);
        REQUIRE(peakGroup.maxPeakFracionalArea == o.maxPeakFracionalArea);
        REQUIRE(peakGroup.maxSignalBaseRatio == o.maxSignalBaseRatio);
        REQUIRE(peakGroup.maxSignalBaselineRatio == o.maxSignalBaselineRatio);
        REQUIRE(peakGroup.maxPeakOverlap == o.maxPeakOverlap);
        REQUIRE(peakGroup.maxQuality == o.maxQuality);
        REQUIRE(peakGroup.avgPeakQuality == o.avgPeakQuality);
        REQUIRE(peakGroup.groupQuality == o.groupQuality);
        REQUIRE(peakGroup.weightedAvgPeakQuality == o.weightedAvgPeakQuality);
        REQUIRE(peakGroup.predictedLabel == o.predictedLabel);
        REQUIRE(peakGroup.expectedAbundance == o.expectedAbundance);
        REQUIRE(peakGroup.isotopeC13count == o.isotopeC13count);
        REQUIRE(peakGroup.deletedFlag == o.deletedFlag);
        REQUIRE(peakGroup.minRt == o.minRt);
        REQUIRE(peakGroup.maxRt == o.maxRt);
        REQUIRE(peakGroup.minMz == o.minMz);
        REQUIRE(peakGroup.maxMz == o.maxMz);
        REQUIRE(peakGroup.parent == o.parent);
        REQUIRE(peakGroup.parentIon == o.parentIon);
        REQUIRE(peakGroup.srmId == o.srmId);
        REQUIRE(peakGroup.isFocused == o.isFocused);
        REQUIRE(peakGroup.label == o.label);
        REQUIRE(peakGroup.goodPeakCount ==
                o.goodPeakCount);
        REQUIRE(peakGroup._type == o._type);
        REQUIRE(peakGroup.hasSlice() == o.hasSlice());
        REQUIRE(peakGroup.tagString == o.tagString);
        REQUIRE(peakGroup.changeFoldRatio ==
                o.changeFoldRatio);
        REQUIRE(peakGroup.changePValue ==
                o.changePValue);
        REQUIRE(peakGroup.peaks.size() == o.peaks.size());
        REQUIRE(peakGroup.samples.size() == o.samples.size());
        REQUIRE(peakGroup.markedBadByCloudModel ==
                o.markedBadByCloudModel);
        REQUIRE(peakGroup.markedGoodByCloudModel ==
                o.markedGoodByCloudModel);
    }

    SUBCASE("Testing copy childrens")
    {
        PeakGroup o;
        o.copyChildren(peakGroup);
        REQUIRE(o.children.size() == peakGroup.children.size());
        REQUIRE(o.childrenBarPlot.size() ==
                peakGroup.childrenBarPlot.size());
        for (unsigned int i = 0; i < o.children.size(); i++)
            REQUIRE(o.children[i].parent == &peakGroup);
        for (unsigned int i = 0; i < o.childrenBarPlot.size(); i++)
            REQUIRE(o.childrenBarPlot[i].parent == &peakGroup);

        for (auto& adductGroup : peakGroup.childAdducts)
            REQUIRE(adductGroup.parent == &peakGroup);
    }

    SUBCASE("Testing clear method")
    {
        PeakGroup o;
        o = peakGroup;
        o.clear();
        REQUIRE(o.peaks.size() == 0);
        REQUIRE(o.children.size() == 0);
        REQUIRE(o.meanMz == 0);
        REQUIRE(o.expectedMz == 0);
        REQUIRE(o.groupRank == INT_MAX);
    }

    SUBCASE("Testing Boolean functions")
    {
        REQUIRE(peakGroup.isMS1() == true);
        REQUIRE(peakGroup.hasCompoundLink() == true);
        REQUIRE(peakGroup.hasSlice() == true);
        REQUIRE(peakGroup.sliceIsZero() == false);
    }

    SUBCASE("Testing reduce functions")
    {
        peakGroup.reduce();
        REQUIRE(peakGroup.peaks.size() == 4);
    }

    SUBCASE("Testing Update Quality function")
    {
        peakGroup.updateQuality();
        REQUIRE(doctest::Approx(peakGroup.maxQuality) == 0.796583);
        REQUIRE(peakGroup.goodPeakCount == 4);
        REQUIRE(doctest::Approx(peakGroup.avgPeakQuality) == 0.752966);
        REQUIRE(doctest::Approx(peakGroup.weightedAvgPeakQuality) == 0.755977);
    }

    SUBCASE("Testing Getting Expected Mz")
    {
        double mz = peakGroup.getExpectedMz(1);
        REQUIRE(doctest::Approx(mz) == 87.0088);
    }

    SUBCASE("Testing Statistics function")
    {
        peakGroup.groupStatistics();

        REQUIRE(doctest::Approx(peakGroup.maxAreaTopIntensity) == 24856.1);
        REQUIRE(doctest::Approx(peakGroup.maxAreaIntensity) == 456148);
        REQUIRE(doctest::Approx(peakGroup.maxHeightIntensity) == 28994.2);
        REQUIRE(doctest::Approx(peakGroup.maxAreaNotCorrectedIntensity) == 559483);
        REQUIRE(doctest::Approx(peakGroup.maxAreaTopNotCorrectedIntensity) == 26187.1);
        REQUIRE(doctest::Approx(peakGroup.meanMz) == 87.0088);
        REQUIRE(doctest::Approx(peakGroup.meanRt) == 8.6007);
        REQUIRE(doctest::Approx(peakGroup.maxNoNoiseObs) == 75);
        REQUIRE(doctest::Approx(peakGroup.minRt) == 7.67238);
        REQUIRE(doctest::Approx(peakGroup.maxRt) == 8.98297);
        REQUIRE(doctest::Approx(peakGroup.minMz) == 87.0086);
        REQUIRE(doctest::Approx(peakGroup.maxMz) == 87.0089);
        REQUIRE(doctest::Approx(peakGroup.maxPeakFracionalArea) == 0.600747);
        REQUIRE(doctest::Approx(peakGroup.maxQuality) == 0.796583);
        REQUIRE(peakGroup.goodPeakCount == 4);
        REQUIRE(doctest::Approx(peakGroup.maxSignalBaselineRatio) == 9.45886);
        REQUIRE(peakGroup.blankMax == 0);
        REQUIRE(doctest::Approx(peakGroup.sampleMax) == 28994.2);
        REQUIRE(doctest::Approx(peakGroup.avgPeakQuality) == 0.752966);
        REQUIRE(doctest::Approx(peakGroup.weightedAvgPeakQuality) == 0.755977);
    }

    SUBCASE("Testing Get Name")
    {
        string str = peakGroup.getName();
        REQUIRE(str == "pyruvate | [M-H]-");

    }

    SUBCASE("Testing Comapre functions")
    {
        PeakGroup secondPeakgroup = secondPeakGroup();
        REQUIRE(peakGroup.compRt(peakGroup, secondPeakgroup) == true);
        REQUIRE(peakGroup.compIntensity(peakGroup, secondPeakgroup) == true);
        REQUIRE(peakGroup.compArea(peakGroup, secondPeakgroup) == true);
        REQUIRE(peakGroup.compQuality(peakGroup, secondPeakgroup) == true);
        REQUIRE(peakGroup.compRank(peakGroup, secondPeakgroup) == true);
        REQUIRE(peakGroup.compRatio(peakGroup, secondPeakgroup) == false);
        REQUIRE(peakGroup.compPvalue(&peakGroup, &secondPeakgroup) == false);
        REQUIRE(peakGroup.compC13(&peakGroup, &secondPeakgroup) == false);
        REQUIRE(peakGroup.compMetaGroup(peakGroup, secondPeakgroup) == false);
        bool res = peakGroup < &secondPeakgroup;
        REQUIRE(res == false);
    }

    SUBCASE("Testing getters and setters")
    {
        REQUIRE(peakGroup.isMS1() == true);
        REQUIRE(peakGroup.hasSrmId() == true);
        REQUIRE(peakGroup.hasCompoundLink() == true);
        REQUIRE(peakGroup.isEmpty() == false);
        REQUIRE(peakGroup.peakCount() == 4);
        REQUIRE(peakGroup.childCount() == 0);
        REQUIRE(peakGroup.childCountBarPlot() == 0);
        REQUIRE(peakGroup.childCountIsoWidget() == 0);
        REQUIRE(peakGroup.hasSlice() == true);
        REQUIRE(peakGroup.sliceIsZero() == false);

        Compound* a = new Compound("C00166", "UTP" ,
                                   "C9H15N2O14P3", 1);
        PeakGroup o (peakGroup);
        o.setCompound(a);

        Compound *b = o.getCompound();
        REQUIRE(b->name == "UTP");
        REQUIRE(b->id == "C00166");
        REQUIRE(b->formula() == "C9H15N2O14P3");
        REQUIRE(b->charge == 1);

        vector<Peak> peaks = o.getPeaks();
        REQUIRE(peaks.size() == 4);

        PeakGroup secondPeakgroup = secondPeakGroup();
        o.setParent(&secondPeakgroup);
        PeakGroup *parent = o.getParent();
        REQUIRE(parent != nullptr);

        mzSlice slice= o.getSlice();
        REQUIRE(doctest::Approx(slice.rt) == 0);
        REQUIRE(doctest::Approx(slice.mz) == 87.0088);
        REQUIRE(doctest::Approx(slice.mzmin) == 87.0079);
        REQUIRE(doctest::Approx(slice.mzmax) == 87.0096);
        REQUIRE(doctest::Approx(slice.rtmax) == 9.39);
        REQUIRE(doctest::Approx(slice.rtmin) == 7.39);
        REQUIRE(slice.srmId == "");

        Adduct *adductIn = new Adduct("[M-H]-", 1, -1, -1.00728);
        o.setAdduct(adductIn);
        Adduct *adduct = o.getAdduct();
        REQUIRE(adduct != nullptr);
        REQUIRE(adduct->getName() == "[M-H]-");
        REQUIRE(adduct->getCharge() == -1);
        REQUIRE(adduct->getNmol() == 1);
        REQUIRE(doctest::Approx(adduct->getMass()) == -1.00728);

        o.addChild(secondPeakgroup);
        REQUIRE(o.getChildren().size() == 1);

        o.setTableName("Table");
        REQUIRE(o.tableName() == "Table");

        o.setSelectedSamples(samples());
        REQUIRE(o.samples[0]->sampleName == "091215_120i");
        REQUIRE(o.samples[1]->sampleName == "091215_120M");
        REQUIRE(o.samples[2]->sampleName == "091215_240i");
        REQUIRE(o.samples[3]->sampleName == "091215_240M");

        vector<mzSample*> sample = samples();
        vector<float> res = o.getOrderedIntensityVector(sample, PeakGroup::QType::AreaTop);
        REQUIRE(res.size() == 4);
        REQUIRE(doctest::Approx(res[0]) == 24856.1);
        REQUIRE(doctest::Approx(res[1]) == 19315.3);
        REQUIRE(doctest::Approx(res[2]) == 18611.3);
        REQUIRE(doctest::Approx(res[3]) == 15058.9);
        REQUIRE(o.deleteChild(&o.children[0]) == true);
        o.deleteChildren();
        REQUIRE(o.getChildren().size() == 0);
    }

    SUBCASE("Testing ExpectedRtDiff")
    {
        float res = peakGroup.expectedRtDiff();
        REQUIRE(doctest::Approx(res) == 0.210696);
    }

    SUBCASE("Test computing Avg blank area")
    {
        vector<EIC*> eics;
        EIC *eic = new EIC();
        vector<float> intensities;
        intensities.push_back(7.68450);
        intensities.push_back(7.77906);
        intensities.push_back(7.87906);
        intensities.push_back(7.97906);
        intensities.push_back(8.57906);
        intensities.push_back(8.67906);
        eic->intensity = intensities;
        vector<float> rts;
        rts.push_back(7.68450);
        rts.push_back(7.77906);
        rts.push_back(7.87906);
        rts.push_back(7.97906);
        rts.push_back(8.57906);
        rts.push_back(8.67906);
        eic->rt = rts;
        eics.push_back(eic);
        peakGroup.computeAvgBlankArea(eics);
        REQUIRE(doctest::Approx(peakGroup.blankMean) == 8.09663);
    }

    SUBCASE("Test Filling missing peaks")
    {
        vector<EIC*> eics;
        EIC *eic = new EIC();
        eic = peakGroup.samples[0]->getEIC(peakGroup.minMz,
                                           peakGroup.maxMz,
                                           peakGroup.minRt,
                                           peakGroup.maxRt,
                                           1, 1, "");
        eics.push_back(eic);
        eic = peakGroup.samples[1]->getEIC(peakGroup.minMz,
                                           peakGroup.maxMz,
                                           peakGroup.minRt,
                                           peakGroup.maxRt,
                                           1, 1, "");
        eics.push_back(eic);
        eic = peakGroup.samples[2]->getEIC(peakGroup.minMz,
                                           peakGroup.maxMz,
                                           peakGroup.minRt,
                                           peakGroup.maxRt,
                                           1, 1, "");
        eics.push_back(eic);
        peakGroup.fillInPeaks(eics);
        REQUIRE(peakGroup.peaks.size() == 4);
    }

    SUBCASE("Testing mass cut off distance")
    {
        MassCutoff * massCutoff = new MassCutoff();
        massCutoff->setMassCutoff(10);
        massCutoff->setMassCutoffType("ppm");
        float res = peakGroup.massCutoffDist(120.0f, massCutoff);
        REQUIRE(doctest::Approx(res) == 274927);
    }

    SUBCASE("Testing group overlap matrix")
    {
        peakGroup.groupOverlapMatrix();
        std::sort(peakGroup.peaks.begin(), peakGroup.peaks.end(), Peak::compIntensity);
        REQUIRE(doctest::Approx(peakGroup.peaks[0].groupOverlapFrac) == -0.044771);
        REQUIRE(doctest::Approx(peakGroup.peaks[1].groupOverlapFrac) == -0.045815);
        REQUIRE(doctest::Approx(peakGroup.peaks[2].groupOverlapFrac) ==  -0.088632);
        REQUIRE(doctest::Approx(peakGroup.peaks[3].groupOverlapFrac) == 0);
    }

    SUBCASE("Test getting representative scans")
    {
        vector<Scan*> scan = peakGroup.getRepresentativeFullScans();
        REQUIRE(scan.size() == 4);
        sort(scan.begin(), scan.end(), Scan::compRt);
        REQUIRE(doctest::Approx(scan[0]->rt) == 8.57253);
        REQUIRE(doctest::Approx(scan[1]->rt) == 8.58375);
        REQUIRE(doctest::Approx(scan[2]->rt) == 8.60168);
        REQUIRE(doctest::Approx(scan[3]->rt) == 8.64482);
    }

    SUBCASE("Test Group Rank")
    {
        peakGroup.calGroupRank(false, 1, 1, 1);
        REQUIRE(doctest::Approx(peakGroup.groupRank) == 0.704177);
    }

    SUBCASE("Test getting fragmentation Events")
    {
        vector<Scan*> scan = peakGroup.getFragmentationEvents();
        REQUIRE(scan.size() == 0);

        Scan* fragmentScan = peakGroup.getAverageFragmentationScan(10);
        REQUIRE(doctest::Approx(fragmentScan->rt) == 8.6007);

        peakGroup.computeFragPattern(10);
    }
}
