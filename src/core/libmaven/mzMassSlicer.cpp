#include "EIC.h"
#include "mavenparameters.h"
#include "mzMassSlicer.h"
#include "mzSample.h"
#include "datastructures/mzSlice.h"
#include "masscutofftype.h"
#include "mzUtils.h"
#include "Matrix.h"
#include "PeakDetector.h"
#include "Scan.h"

using namespace mzUtils;

MassSlices::MassSlices()
{
    _maxSlices=INT_MAX;
    _minRt=FLT_MIN; _minMz=FLT_MIN; _minIntensity=FLT_MIN;
    _maxRt=FLT_MAX; _maxMz=FLT_MAX; _maxIntensity=FLT_MAX;
    _minCharge=0; _maxCharge=INT_MAX;
    massCutoff=NULL;
}

MassSlices::~MassSlices() { delete_all(slices); cache.clear(); }

void MassSlices::sendSignal(const string& progressText,
                unsigned int completed_samples,
                int total_samples)
{
    mavenParameters->sig(progressText, completed_samples, total_samples);
}

/**
 * MassSlices::algorithmA This is function is called when mass Slicing using 
 * AlgorithmB returns no slices. The slices here are created using the filterLine
 * in Mzml and Mzxml files.
 */
void MassSlices::algorithmA() {
    // clear cache
    delete_all(slices);
    slices.clear();
    cache.clear();
    map< string, int> seen;

    //#pragma omp parallel for ordered
    // Iterate over every sample
    for(unsigned int i=0; i < samples.size(); i++) {
        // Iterate over every scan 
        for(unsigned int j=0; j < samples[i]->scans.size(); j++ ) {
            // Make temprorary element scan with present scan
            Scan* scan = samples[i]->scans[j];

            // Check if filterLine(SRM transition for MS-MS, segment for LC-MS) for the scan is empty
            if ( scan->filterLine.empty() ) continue;

            if ( seen.count( scan->filterLine ) ) continue;

            //Create new slice for every filterLine
            mzSlice* s = new mzSlice(scan->filterLine);
            slices.push_back(s);
            seen[ scan->filterLine ]=1;
        }
    }
    cerr << "#algorithmA" << slices.size() << endl;
}

void MassSlices::stopSlicing() {
    if (slices.size() > 0) {
        delete_all(slices);
        slices.clear();
        cache.clear();
    }
}

void MassSlices::algorithmB(MassCutoff* massCutoff, int rtStep )
{
    // clear all previous data
    delete_all(slices);
    slices.clear();
    cache.clear();

    float rtWindow = 2.0f;
    this->massCutoff = massCutoff;

    int totalScans = 0;
    int currentScans = 0;

    // Calculate the total number of scans
    for (auto s : samples)
        totalScans += s->scans.size();

    // Calculating the rt window using average distance between RTs and
    // mutiplying it with rtStep (default 2.0)
    if (samples.size() > 0 and rtStep > 0) {
        rtWindow = accumulate(begin(samples),
                              end(samples),
                              0.0f,
                              [rtStep](float sum, mzSample* sample) {
                                  return sum + (sample->getAverageFullScanTime()
                                                * rtStep);
                              }) / static_cast<float>(samples.size());
    }
    cerr << "RT window used: " << rtWindow << endl;

    sendSignal("Status", 0 , 1);

    // #pragma omp parallel for ordered
    // Looping over every sample
    for (unsigned int i = 0; i < samples.size(); i++) {
        if (slices.size() > _maxSlices) break;

        // Check if Peak detection has been cancelled by the user
        if (mavenParameters->stop) {
            stopSlicing();
            break;
        }

        string num;
        if (i==0) {
            num = "st";
        } else if (i==1) {
            num = "nd";
        } else if (i==2) {
            num = "rd";
        } else {
            num = "th";
        }

        // updating progress on samples
        if (mavenParameters->showProgressFlag ) {
            string progressText = to_string(i + 1)
                                  + num
                                  + " out of "
                                  + to_string(mavenParameters->samples.size())
                                  + " sample(s) processing…";
            sendSignal(progressText, currentScans, totalScans);
        }

        // #pragma omp cancel for
        // for loop for iterating over every scan of a sample
        for (auto scan : samples[i]->scans) {
            // Check if Peak detection has been cancelled by the user
            if (mavenParameters->stop) {
                stopSlicing();
                break;
            }

            currentScans++;

            if (scan->mslevel != 1)
                continue;

            // Checking if RT is in the given min to max RT range
            if (_maxRt && !isBetweenInclusive(scan->rt, _minRt, _maxRt))
                continue;

            float rt = scan->rt;

            for (unsigned int k = 0; k < scan->nobs(); k++) {
                float mz = scan->mz[k];
                float intensity = scan->intensity[k];

                // Checking if mz, intensity are within specified ranges
                if (_maxMz && !isBetweenInclusive(mz,
                                                  _minMz,
                                                  _maxMz)) {
                    continue;
                }
                if (_maxIntensity && !isBetweenInclusive(intensity,
                                                         _minIntensity,
                                                         _maxIntensity)) {
                    continue;
                }

                // create new slice with the given bounds
                float cutoff = massCutoff->massCutoffValue(mz);
                mzSlice* s = new mzSlice(mz - cutoff,
                                         mz + cutoff,
                                         rt - rtWindow,
                                         rt + rtWindow);
                s->ionCount = intensity;
                s->rt = scan->rt;
                s->mz = mz;
                slices.push_back(s);
            }

            // progress update 
            if (mavenParameters->showProgressFlag ) {
                string progressText = to_string(i + 1)
                                      + num
                                      + " out of "
                                      + to_string(mavenParameters->samples.size())
                                      + " sample(s) processing…\n"
                                      + to_string(slices.size())
                                      + " slices created ";
                sendSignal(progressText,currentScans,totalScans);
            }
        }
    }

    cerr << "Found " << slices.size() << " slices" << endl;

    // before reduction sort by mz first then by rt
    sort(begin(slices),
         end(slices),
         [](const mzSlice* slice, const mzSlice* compSlice) {
             if (slice->mz == compSlice->mz) {
                 return slice->rt < compSlice->rt;
             }
             return slice->mz < compSlice->mz;
         });

    // reduction: merge slices sharing significant ROI
    for (auto first = begin(slices); first != end(slices); ++first) {
        if (mavenParameters->stop) {
            stopSlicing();
            break;
        }

        auto firstSlice = *first;
        if (mzUtils::almostEqual(firstSlice->ionCount, -1.0f))
            continue;

        for (auto second = next(first); second != end(slices); ++second) {
            auto secondSlice = *second;

            // stop iterating if the rest of the slices are too far
            if (firstSlice->mzmax < secondSlice->mzmin)
                break;

            if (mzUtils::almostEqual(secondSlice->ionCount, -1.0f))
                continue;

            // check if center of one of the slices lies in the other
            if ((firstSlice->mz > secondSlice->mzmin
                 && firstSlice->mz < secondSlice->mzmax
                 && firstSlice->rt > secondSlice->rtmin
                 && firstSlice->rt < secondSlice->rtmax)
                ||
                (secondSlice->mz > firstSlice->mzmin
                 && secondSlice->mz < firstSlice->mzmax
                 && secondSlice->rt > firstSlice->rtmin
                 && secondSlice->rt < firstSlice->rtmax)) {
                firstSlice->ionCount = std::max(firstSlice->ionCount,
                                                secondSlice->ionCount);
                firstSlice->rtmax = std::max(firstSlice->rtmax,
                                             secondSlice->rtmax);
                firstSlice->rtmin = std::min(firstSlice->rtmin,
                                             secondSlice->rtmin);
                firstSlice->mzmax = std::max(firstSlice->mzmax,
                                             secondSlice->mzmax);
                firstSlice->mzmin = std::min(firstSlice->mzmin,
                                             secondSlice->mzmin);

                firstSlice->mz = (firstSlice->mzmin + firstSlice->mzmax) / 2.0f;
                firstSlice->rt = (firstSlice->rtmin + firstSlice->rtmax) / 2.0f;
                float cutoff = massCutoff->massCutoffValue(firstSlice->mz);

                // make sure that mz window does not get out of control
                if (firstSlice->mzmin < firstSlice->mz - cutoff)
                    firstSlice->mzmin =  firstSlice->mz - cutoff;
                if (firstSlice->mzmax > firstSlice->mz + cutoff)
                    firstSlice->mzmax =  firstSlice->mz + cutoff;

                // recalculate center mz in case bounds changed
                firstSlice->mz = (firstSlice->mzmin + firstSlice->mzmax) / 2.0f;

                // flag this slice as already merged, and ignore henceforth
                secondSlice->ionCount = -1.0f;
            }
        }
        sendSignal("Reducing slices…", first - begin(slices), slices.size());
    }

    // remove merged slices
    slices.erase(remove_if(slices.begin(),
                           slices.end(),
                           [](mzSlice* slice) {
                               return (slice->ionCount == -1.0f);
                           }),
                 slices.end());

    cerr << "Reduced to " << slices.size() << " slices" << endl;

    float threshold = 100;
    removeDuplicateSlices(massCutoff, threshold);

    sort(slices.begin(), slices.end(), mzSlice::compMz);
    mergeNeighbouringSlices(massCutoff, rtWindow);
    adjustSlices();

    cerr << "After final merging and adjustments, "
         << slices.size()
         << " slices remain"
         << endl;
    sendSignal("Mass slicing done.", 1 , 1);
}

void MassSlices::algorithmC(float ppm, float minIntensity, float rtWindow) {
    delete_all(slices);
    slices.clear();
    cache.clear();

    for(unsigned int i=0; i < samples.size(); i++) {
        mzSample* s = samples[i];
        for(unsigned int j=0; j < s->scans.size(); j++) {
            Scan* scan = samples[i]->scans[j];
            if (scan->mslevel != 1 ) continue;
            vector<int> positions = scan->intensityOrderDesc();
            for(unsigned int k=0; k< positions.size() && k<10; k++ ) {
                int pos = positions[k];
                if (scan->intensity[pos] < minIntensity) continue;
                float rt = scan->rt;
                float mz = scan->mz[ pos ];
                float mzmax = mz + mz/1e6*ppm;
                float mzmin = mz - mz/1e6*ppm;
                if(! sliceExists(mzmin, mzmax, rt-2*rtWindow, rt+2*rtWindow) ) {
                    mzSlice* s = new mzSlice(mzmin,mzmax, rt-2*rtWindow, rt+2*rtWindow);
                    s->ionCount = scan->intensity[pos];
                    s->rt=scan->rt;
                    s->mz=mz;
                    slices.push_back(s);
                    int mzRange = mz*10;
                    cache.insert( pair<int,mzSlice*>(mzRange, s));
                }
            }
        }
    }
    cerr << "#algorithmC" << slices.size() << endl;
}

mzSlice* MassSlices::sliceExists(float mzMinBound,
                                 float mzMaxBound,
                                 float rtMinBound,
                                 float rtMaxBound)
{
    // calculate center for given bounds
    float mz = (mzMinBound + mzMaxBound) / 2.0f;
    float rt = (rtMinBound + rtMaxBound) / 2.0f;

    // putting all mz slices stored in cache within a particular range in ppp
    auto subcache = cache.equal_range(int(mz * 10));
    auto it = subcache.first;

    float bestDist = FLT_MAX;
    mzSlice* best = nullptr;

    // For loop to iterate till best MZ slice becomes second
    for (; it != subcache.second; ++it) {
        mzSlice* slice = (*it).second;
        float sliceMzMin = slice->mzmin;
        float sliceMzMax = slice->mzmax;
        float sliceRtMin = slice->rtmin;
        float sliceRtMax = slice->rtmax;

        // calculate center for the current slice
        float sliceMz = slice->mz;
        float sliceRt = slice->rt;

        // check if center of given bounds lies within slice, or whether center
        // of the slice lies within the given bounds
        if ((mz > sliceMzMin
             && mz < sliceMzMax
             && rt > sliceRtMin
             && rt < sliceRtMax)
            ||
            (sliceMz > mzMinBound
             && sliceMz < mzMaxBound
             && sliceRt > rtMinBound
             && sliceRt < rtMaxBound)) {
            float dist = hypotf(sliceMz - mz, sliceRt - rt);
            if (dist < bestDist) {
                best = slice;
                bestDist = dist;
            }
        }
    }
    return best;
}

void MassSlices::mergeNeighbouringSlices(MassCutoff* massCutoff,
                                         float rtTolerance)
{
    auto slicesInProximity = [&](vector<mzSample*>& samples,
                                 mzSlice* slice,
                                 mzSlice* comparisonSlice) {
        auto mz = slice->mz;
        auto mzMin = slice->mzmin;
        auto mzMax = slice->mzmax;
        auto rtMin = slice->rtmin;
        auto rtMax = slice->rtmax;
        auto comparisonMz = comparisonSlice->mz;
        auto comparisonMzMin = comparisonSlice->mzmin;
        auto comparisonMzMax = comparisonSlice->mzmax;
        auto comparisonRtMin = comparisonSlice->rtmin;
        auto comparisonRtMax = comparisonSlice->rtmax;
        auto mzCenter = (mz + comparisonMz) / 2.0f;

        // check to make sure slices are close to each other (or have some
        // overlap in mz domain); the tolerance is multiplied 10x so as to
        // include slices that may be further apart but should be merged
        float massTolerance = 10.0f * massCutoff->massCutoffValue(mzCenter);
        if (!(abs(mzCenter - mz) <= massTolerance
              && abs(mzCenter - comparisonMz) <= massTolerance)) {
            return make_pair(false, false);
        }

        // check if common RT regions exist between the slices being compared
        auto commonLowerRt = 0.0f;
        auto commonUpperRt = 0.0f;
        if (rtMin <= comparisonRtMin && rtMax >= comparisonRtMax) {
            commonLowerRt = comparisonRtMin;
            commonUpperRt = comparisonRtMax;
        } else if (rtMin >= comparisonRtMin && rtMax <= comparisonRtMax) {
            commonLowerRt = rtMin;
            commonUpperRt  = rtMax;
        } else if (rtMin >= comparisonRtMin && rtMin <= comparisonRtMax) {
            commonLowerRt = rtMin;
            commonUpperRt = min(rtMax, comparisonRtMax);
        } else if (rtMax >= comparisonRtMin && rtMax <= comparisonRtMax) {
            commonLowerRt = max(rtMin, comparisonRtMin);
            commonUpperRt = rtMax;
        }
        if (commonLowerRt == 0.0f && commonUpperRt == 0.0f)
            return make_pair(false, true);

        auto highestIntensity = 0.0f;
        auto mzAtHighestIntensity = 0.0f;
        auto rtAtHighestIntensity = 0.0f;
        auto highestCompIntensity = 0.0f;
        auto mzAtHighestCompIntensity = 0.0f;
        auto rtAtHighestCompIntensity = 0.0f;
        for (auto sample : samples) {
            // obtain EICs for the two slices
            auto eic = sample->getEIC(mzMin,
                                      mzMax,
                                      rtMin,
                                      rtMax,
                                      1,
                                      1,
                                      "");
            auto comparisonEic = sample->getEIC(comparisonMzMin,
                                                comparisonMzMax,
                                                comparisonRtMin,
                                                comparisonRtMax,
                                                1,
                                                1,
                                                "");

            // find out the highest intensity's mz and rt in the EICs
            for (size_t i = 0; i < eic->size(); ++i) {
                auto intensityAtIdx = eic->intensity[i];
                if (intensityAtIdx > highestIntensity) {
                    highestIntensity = intensityAtIdx;
                    mzAtHighestIntensity = eic->mz[i];;
                    rtAtHighestIntensity = eic->rt[i];
                }
            }
            for (size_t i = 0; i < comparisonEic->size(); ++i) {
                auto intensityAtIdx = comparisonEic->intensity[i];
                if (intensityAtIdx > highestCompIntensity) {
                    highestCompIntensity = intensityAtIdx;
                    mzAtHighestCompIntensity = comparisonEic->mz[i];
                    rtAtHighestCompIntensity = comparisonEic->rt[i];
                }
            }
            delete eic;
            delete comparisonEic;
        }

        if (highestIntensity == 0.0f && highestCompIntensity == 0.0f)
            return make_pair(false, true);

        // calculate and check for rt difference and mz difference, if
        // conditions are satisfied, mark the comparison slice to be merged
        auto rtDelta = abs(rtAtHighestIntensity - rtAtHighestCompIntensity);
        auto mzCenterForIntensity = (mzAtHighestIntensity
                                  + mzAtHighestCompIntensity) / 2.0f;
        auto massToleranceForIntensity =
            massCutoff->massCutoffValue(mzCenterForIntensity);
        auto mzDeltaNeg = abs(mzCenterForIntensity - mzAtHighestIntensity );
        auto mzDeltaPos = abs(mzAtHighestCompIntensity - mzCenterForIntensity );
        if (rtDelta <= rtTolerance
            && mzDeltaNeg <= massToleranceForIntensity
            && mzDeltaPos <= massToleranceForIntensity) {
            return make_pair(true, true);
        }

        return make_pair(false, true);
    };
    _mergeSlices(slicesInProximity,
                 massCutoff,
                 "Merging related slices in samples…");
}

void MassSlices::_mergeSlices(const function<pair<bool, bool>(vector<mzSample*>&,
                                                              mzSlice *,
                                                              mzSlice *)> &compareSlices,
                              MassCutoff* massCutoff,
                              const string &updateMessage)
{
    // lambda to help expand a given slice by merging a vector of slices into it
    auto expandSlice = [&](mzSlice* mergeInto, vector<mzSlice*> slices) {
        if (slices.empty())
            return;

        for (auto slice : slices) {
            mergeInto->ionCount = std::max(mergeInto->ionCount, slice->ionCount);
            mergeInto->rtmax = std::max(mergeInto->rtmax, slice->rtmax);
            mergeInto->rtmin = std::min(mergeInto->rtmin, slice->rtmin);
            mergeInto->mzmax = std::max(mergeInto->mzmax, slice->mzmax);
            mergeInto->mzmin = std::min(mergeInto->mzmin, slice->mzmin);
        }

        // calculate the new midpoints
        mergeInto->mz = (mergeInto->mzmin + mergeInto->mzmax) / 2.0f;
        mergeInto->rt = (mergeInto->rtmin + mergeInto->rtmax) / 2.0f;

        // make sure that mz window does not get out of control
        auto cutoff = massCutoff->massCutoffValue(mergeInto->mz);
        if (mergeInto->mzmin < mergeInto->mz - cutoff)
            mergeInto->mzmin =  mergeInto->mz - cutoff;
        if (mergeInto->mzmax > mergeInto->mz + cutoff)
            mergeInto->mzmax =  mergeInto->mz + cutoff;

        mergeInto->mz = (mergeInto->mzmin + mergeInto->mzmax) / 2.0f;
    };

    for(auto it = begin(slices); it != end(slices); ++it) {
        if (mavenParameters->stop) {
            stopSlicing();
            break;
        }

        sendSignal(updateMessage, it - begin(slices), slices.size());

        auto slice = *it;
        vector<mzSlice*> slicesToMerge;

        // search ahead
        for (auto ahead = next(it);
             ahead != end(slices) && it != end(slices);
             ++ahead) {
            auto comparisonSlice = *ahead;
            auto comparison = compareSlices(samples, slice, comparisonSlice);
            auto shouldMerge = comparison.first;
            auto continueIteration = comparison.second;
            if (shouldMerge)
                slicesToMerge.push_back(comparisonSlice);
            if (!continueIteration)
                break;
        }

        // search behind
        for (auto behind = prev(it);
             behind != begin(slices) && it != begin(slices);
             --behind) {
            auto comparisonSlice = *behind;
            auto comparison = compareSlices(samples, slice, comparisonSlice);
            auto shouldMerge = comparison.first;
            auto continueIteration = comparison.second;
            if (shouldMerge)
                slicesToMerge.push_back(comparisonSlice);
            if (!continueIteration)
                break;
        }

        // expand the current slice by merging all slices classified to be
        // part of the same, and then remove (and free) the slices already
        // merged
        expandSlice(slice, slicesToMerge);
        for (auto merged : slicesToMerge) {
            slices.erase(remove_if(begin(slices),
                                   end(slices),
                                   [&](mzSlice* s) { return s == merged; }),
                         slices.end());
            delete merged;
        }
        it = find_if(begin(slices),
                     end(slices),
                     [&](mzSlice* s) { return s == slice; });
    }
}

void MassSlices::adjustSlices()
{
    size_t progressCount = 0;
    for (auto slice : slices) {
        if (mavenParameters->stop) {
            stopSlicing();
            break;
        }

        auto eics = PeakDetector::pullEICs(slice,
                                           mavenParameters->samples,
                                           mavenParameters);
        float highestIntensity = 0.0f;
        float mzAtHighestIntensity = 0.0f;
        for (auto eic : eics) {
            for (size_t i = 0; i < eic->size(); ++i) {
                auto intensityAtIndex = eic->intensity[i];
                if (intensityAtIndex > highestIntensity) {
                    highestIntensity = intensityAtIndex;
                    mzAtHighestIntensity = eic->mz[i];
                }
            }
        }
        float cutoff = mavenParameters->massCutoffMerge
                                      ->massCutoffValue(mzAtHighestIntensity);
        slice->mzmin =  mzAtHighestIntensity - cutoff;
        slice->mzmax =  mzAtHighestIntensity + cutoff;
        slice->mz = (slice->mzmin + slice->mzmax) / 2.0f;

        delete_all(eics);

        ++progressCount;
        sendSignal("Adjusting slices…", progressCount, slices.size());
    }
}

void MassSlices::removeDuplicateSlices(MassCutoff *massCutoff, float threshold){

    vector<mzSlice*> returnSlices;
    mzSlice* slice;
    multimap<int,int>vectorCache;
    
   for(int i=0; i<slices.size(); i++) {
        slice = slices[i];
        float mz = slice->mz;
        multimap<int, int>::iterator start, end;
        
        start = vectorCache.equal_range( (int) (mz* 10 - 1) ).first;
        end = vectorCache.equal_range( (int) (mz* 10 + 1) ).second;
        float mzOverlap =  0.0;
        float rtOverlap = 0.0;
        float overlapArea, bestOverlapArea = 0.0;
        int bestSliceNum = -1;

        for(; start != end; start++) {
            int thisSliceNum = (*start).second;
            mzSlice *thisSlice = returnSlices[thisSliceNum];

            float low = thisSlice->mzmin > slice->mzmin ? thisSlice->mzmin : slice->mzmin;
            float high = thisSlice->mzmax < slice->mzmax ? thisSlice->mzmax : slice->mzmax;
            mzOverlap = high-low;

            low = thisSlice->rtmin > slice->rtmin ? thisSlice->rtmin : slice->rtmin;
            high = thisSlice->rtmax < slice->rtmax ? thisSlice->rtmax : slice->rtmax;
            rtOverlap = high-low;

            if(mzOverlap>0 && rtOverlap>0) overlapArea = mzOverlap * rtOverlap;
            else overlapArea = 0;

            float area1 = (thisSlice->mzmax-thisSlice->mzmin) * (thisSlice->rtmax-thisSlice->rtmin);
            float area2 = (slice->mzmax-slice->mzmin) * (slice->rtmax-slice->rtmin);
            float area = area1 < area2 ? area1 : area2;

            if (overlapArea/area >= threshold/100 && overlapArea > bestOverlapArea){
                bestOverlapArea = overlapArea;
                bestSliceNum = thisSliceNum;
            }
            
        }

        if(bestSliceNum >= 0){
            mzSlice* Z = returnSlices[bestSliceNum];
            
            Z->ionCount = std::max((float) Z->ionCount, (float ) slice->ionCount);
            Z->rtmax = std::max((float)Z->rtmax, (float)slice->rtmax);
            Z->rtmin = std::min((float)Z->rtmin, (float)slice->rtmin);
            Z->mzmax = std::max((float)Z->mzmax, (float)slice->mzmax);
            Z->mzmin = std::min((float)Z->mzmin, (float)slice->mzmin);

            //make sure that mz windown doesn't get out of control
            if (Z->mzmin < mz-massCutoff->massCutoffValue(mz)) Z->mzmin =  mz-massCutoff->massCutoffValue(mz);
            if (Z->mzmax > mz+massCutoff->massCutoffValue(mz)) Z->mzmax =  mz+massCutoff->massCutoffValue(mz);
            Z->mz = (Z->mzmin + Z->mzmax) / 2; Z->rt=(Z->rtmin + Z->rtmax) / 2;
        }
        else{
            vectorCache.insert(pair<int,int>(int (mz*10),returnSlices.size()));
            returnSlices.push_back(slice);
        }
    }
    slices = returnSlices;
}
