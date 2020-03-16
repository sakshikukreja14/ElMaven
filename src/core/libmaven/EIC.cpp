#include "EIC.h"
#include "Peak.h"
#include "PeakDetector.h"
#include "PeakGroup.h"
#include "SavGolSmoother.h"
#include "Scan.h"
#include "classifierNeuralNet.h"
#include "databases.h"
#include "datastructures/adduct.h"
#include "datastructures/mzSlice.h"
#include "doctest.h"
#include "masscutofftype.h"
#include "mavenparameters.h"
#include "mzPatterns.h"
#include "mzSample.h"

/**
 * @file EIC.cpp
 * @author Sabu George
 * @author Kiran
 * @author Sahil
 * @version 769
 */
EIC::EIC()
{
    _sample = NULL;
    spline = NULL;
    baseline = NULL;
    _mzmin = _mzmax = _rtmin = _rtmax = 0;
    _maxIntensity = _totalIntensity = 0;
    _rtAtMaxIntensity = 0.0f;
    _mzAtMaxIntensity = 0.0f;
    _maxAreaTopIntensity = 0;
    _maxAreaIntensity = 0;
    _maxAreaTopNotCorrectedIntensity = 0;
    _maxAreaNotCorrectedIntensity = 0;
    _filterSignalBaselineDiff = 0;
    _eic_noNoiseObs = 0;
    smootherType = GAUSSIAN;
    _baselineMode = BaselineMode::Threshold;
    baselineSmoothingWindow = 5;
    baselineDropTopX = 60;
    _aslsAsymmetry = 80;
    _aslsSmoothness = 2;
    for (unsigned int i = 0; i < 4; i++)
        color[i] = 0;
}

EIC::~EIC()
{
    if (spline != NULL)
        delete[] spline;
    spline = NULL;
    if (baseline != NULL)
        delete[] baseline;
    baseline = NULL;
    peaks.clear();
}

EIC* EIC::eicMerge(const vector<EIC*>& eics)
{
    // Merge to 776
    EIC* meic = new EIC();

    unsigned int maxlen = 0;
    float minRt = DBL_MAX;
    float maxRt = DBL_MIN;
    for (unsigned int i = 0; i < eics.size(); i++) {
        if (eics[i]->size() > maxlen)
            maxlen = eics[i]->size();
        if (eics[i]->_rtmin < minRt)
            minRt = eics[i]->_rtmin;
        if (eics[i]->_rtmax > maxRt)
            maxRt = eics[i]->_rtmax;
    }

    if (maxlen == 0)
        return meic;

    // create new EIC
    meic->_sample = NULL;
    vector<float> intensity(maxlen, 0);
    vector<float> rt(maxlen, 0);
    vector<int> scans(maxlen, 0);
    vector<float> mz(maxlen, 0);
    vector<int> mzcount(maxlen, 0);

    // smoothing   //initalize time array
    for (unsigned int i = 0; i < maxlen; i++) {
        rt[i] = minRt + i * ((maxRt - minRt) / maxlen);
        scans[i] = i;
    }

    // combine intensity data from all pulled eics
    for (unsigned int i = 0; i < eics.size(); i++) {
        EIC* e = eics[i];
        for (unsigned int j = 0; j < e->size(); j++) {
            unsigned int bin = ((e->rt[j] - minRt) / (maxRt - minRt) * maxlen);
            if (bin >= maxlen)
                bin = maxlen - 1;

            if (e->spline and e->spline[j] > 0) {
                intensity[bin] += e->spline[j];
            } else {
                intensity[bin] += e->intensity[j];
            }

            if (e->mz[j] > 0) {
                mz[bin] += e->mz[j];
                mzcount[bin]++;
            }
        }
    }

    unsigned int eicCount = eics.size();
    for (unsigned int i = 0; i < maxlen; i++) {
        intensity[i] /= eicCount;
        if (intensity[i] > meic->_maxIntensity) {
            meic->_maxIntensity = intensity[i];
            meic->_rtAtMaxIntensity = rt[i];
            meic->_mzAtMaxIntensity = mz[i];
        }
        if (mzcount[i])
            mz[i] /= mzcount[i];
        meic->_totalIntensity += intensity[i];
    }

    // copy to new EIC
    meic->_rtmin = minRt;
    meic->_rtmax = maxRt;
    meic->intensity = intensity;
    meic->rt = rt;
    meic->scannum = scans;
    meic->mz = mz;
    meic->_sampleName = eics[0]->_sampleName;
    meic->_sample = eics[0]->_sample;
    return meic;
}

bool EIC::_clearBaseline()
{
    if (baseline != nullptr) {  // delete previous baseline if exists
        delete[] baseline;
        baseline = nullptr;
        _eic_noNoiseObs = 0;
    }

    size_t n = intensity.size();
    if (!n)
        return false;

    baseline = new float[n];
    std::fill_n(baseline, n, 0.0f);

    return true;
}

void EIC::_computeAsLSBaseline(const float lambda,
                               const float p,
                               const int numIterations)
{
    // Eigen behaves better with double values compared to float
    vector<double> intensity;
    for (unsigned int i = 0; i < this->intensity.size(); ++i)
        intensity.push_back(static_cast<double>(this->intensity[i]));

    auto originalSize = intensity.size();

    // decimate the signal, if it is of very high-resolution
    auto resamplingFactor = mzUtils::approximateResamplingFactor(originalSize);
    intensity = mzUtils::resample(intensity, 1, resamplingFactor);

    using namespace Eigen;
    auto n = static_cast<unsigned int>(intensity.size());

    // Perform difference operation (pseudo-differentiation) on a sparse matrix
    auto diff = [](Eigen::SparseMatrix<double> mat) {
        Eigen::SparseMatrix<double> E1 =
            mat.block(0, 0, mat.rows() - 1, mat.cols());
        Eigen::SparseMatrix<double> E2 =
            mat.block(1, 0, mat.rows() - 1, mat.cols());
        return Eigen::SparseMatrix<double>(E2 - E1);
    };

    // create a sparse identity matrix of size `n`
    Eigen::SparseMatrix<double> ident(n, n);
    ident.setIdentity();

    // compute approximate second derivative of the identity matrix
    auto D = diff(diff(ident));

    // create a weights vector of length `n` initially containing only ones
    // for coefficients
    auto w = VectorXd(ArrayXd::Ones(n));

    // instantiate a Cholesky decomposition linear solver for sparse matrices
    SimplicialCholesky<SparseMatrix<double>> solver;

    // instantiate a vector for intensity and one to store baseline
    VectorXd intensityVec = Map<VectorXd>(intensity.data(), n);
    VectorXd baselineVec;

    // helper matrices for `select` operation used when creating binary vectors
    auto ones = VectorXd(MatrixXd::Ones(n, 1));
    auto zeros = VectorXd(MatrixXd::Zero(n, 1));

    // TODO: ideally this should converge to a point where the baseline does
    // not change anymnore, but since we do not have good float comparators
    // yet we can use a decent number of iterations to get as close to the true
    // baseline as possible
    for (int i = 0; i < numIterations; ++i) {
        // create a square matrix 'W' with the weights as its diagonal
        auto W = SparseMatrix<double>(w.asDiagonal());

        // compute 'A' and 'b', and then solve for 'x' that satisfies the
        // equation 'A·x = b', where x will be the iteratively estimated
        // baseline
        auto A = W + (lambda * (D.transpose() * D));
        solver.compute(A);
        auto b = VectorXd(w.array() * intensityVec.array());
        baselineVec = solver.solve(b);

        // calculate weights for the next iteration
        auto gtBin = VectorXd(
            ((intensityVec - baselineVec).array() > 0.0).select(ones, zeros));
        auto ltBin = VectorXd(
            ((intensityVec - baselineVec).array() < 0.0).select(ones, zeros));
        w = (p * gtBin) + ((1.0f - p) * ltBin);
    }

    // copy data from an Eigen vector to std::vector
    vector<double> tempVector(
        baselineVec.data(),
        baselineVec.data() + baselineVec.rows() * baselineVec.cols());

    // interpolate the signal after possible decimation
    tempVector = mzUtils::resample(tempVector, resamplingFactor, 1);

    // since the interpolated vector may not be of the same size as the original
    // intensity vector, we remove/pad (with zeros) until they are the same size
    while (tempVector.size() < originalSize)
        tempVector.push_back(0.0);
    while (tempVector.size() > originalSize)
        tempVector.pop_back();

    // clip negative values from the vector and switch back to float
    vector<float> clippedFloatBaseline;
    for (unsigned int i = 0; i < tempVector.size(); ++i) {
        auto val = static_cast<float>(tempVector[i]);
        if (val < 0.0f)
            clippedFloatBaseline.push_back(0.0f);
        else
            clippedFloatBaseline.push_back(val);
    }

    // since baseline right now is an array of float, the data from clipped
    // vector is copied to this array
    std::copy(begin(clippedFloatBaseline), end(clippedFloatBaseline), baseline);
}

void EIC::_computeThresholdBaseline(const int smoothingWindow,
                                    const int dropTopX)
{
    auto n = intensity.size();

    // sort intensity vector
    vector<float> tmpv = intensity;
    std::sort(tmpv.begin(), tmpv.end());

    // compute maximum intensity of baseline, any point above this value will
    // be dropped. User specifies quantile of points to keep, for example
    // drop 60% of highest intensities = cut at 40% value;

    float cutvalueF = (100.0 - (float)dropTopX) / 101;
    unsigned int pos = tmpv.size() * cutvalueF;

    float qcut = 0;
    pos < tmpv.size() ? qcut = tmpv[pos] : qcut = tmpv.back();

    // drop all points above maximum baseline value
    for (size_t i = 0; i < n; i++) {
        if (intensity[i] > qcut) {
            baseline[i] = qcut;
        } else {
            baseline[i] = intensity[i];
        }
    }

    // smooth baseline
    gaussian1d_smoothing(n, smoothingWindow, baseline);

    // count number of observation in EIC above baseline
    for (unsigned long int i = 0; i < n; i++) {
        if (intensity[i] > baseline[i])
            _eic_noNoiseObs++;
    }
}

void EIC::computeBaseline()
{
    if (!_clearBaseline())
        return;

    switch (_baselineMode) {
    case BaselineMode::Threshold:
        _computeThresholdBaseline(baselineSmoothingWindow, baselineDropTopX);
        break;
    case BaselineMode::AsLSSmoothing:
        _computeAsLSBaseline(pow(10.0f, static_cast<float>(_aslsSmoothness)),
                             static_cast<float>(_aslsAsymmetry) / 100.0f);
        break;
    }
}

void EIC::computeSpline(int smoothWindow)
{
    // Merged to 776
    int n = intensity.size();

    if (n == 0)
        return;
    if (this->spline != NULL) {
        delete[] spline;
        spline = NULL;
    }

    try {
        this->spline = new float[n];
        for (int i = 0; i < n; i++)
            spline[i] = 0;
    } catch (...) {
        cerr << "Exception caught while allocating memory " << n << "floats "
             << endl;
    }

    // initalize spline, set to intensity vector
    for (int i = 0; i < n; i++)
        spline[i] = intensity[i];

    // smoothing window is too large
    if (smoothWindow > n / 3)
        smoothWindow = n / 3;
    if (smoothWindow <= 1)
        return;  // nothing to smooth get out

    if (smootherType == SAVGOL) {
        // SAVGOL SMOOTHER
        mzUtils::SavGolSmoother smoother(smoothWindow, smoothWindow, 4);
        vector<float> smoothed = smoother.Smooth(intensity);
        for (int i = 0; i < n; i++)
            spline[i] = smoothed[i];
    } else if (smootherType == GAUSSIAN) {
        // GAUSSIAN SMOOTHER
        gaussian1d_smoothing(n, smoothWindow, spline);
    } else if (smootherType == AVG) {
        float* y = new float[n];
        for (int i = 0; i < n; i++)
            y[i] = intensity[i];
        smoothAverage(y, spline, smoothWindow, n);
        delete[] y;
    }
}

Peak* EIC::addPeak(int peakPos)
{
    peaks.push_back(Peak(this, peakPos));
    return &peaks[peaks.size() - 1];
}

void EIC::getPeakPositions(int smoothWindow)
{
    unsigned int N = intensity.size();
    if (N == 0)
        return;

    computeSpline(smoothWindow);
    if (spline == NULL)
        return;

    findPeaks();

    computeBaseline();
    getPeakStatistics();

    filterPeaks();
}

void EIC::findPeaks()
{
    peaks.clear();

    size_t N = intensity.size();

    for (size_t i = 1; i < N - 1; i++) {
        if (spline[i] > spline[i - 1] && spline[i] > spline[i + 1]) {
            addPeak(i);
        } else if (spline[i] > spline[i - 1] && spline[i] == spline[i + 1]) {
            float highpoint = spline[i];
            while (i < N - 1) {
                i++;
                if (spline[i + 1] == highpoint)
                    continue;
                if (spline[i + 1] > highpoint)
                    break;
                if (spline[i + 1] < highpoint) {
                    addPeak(i);
                    break;
                }
            }
        }
    }
}

void EIC::findPeakBounds(Peak& peak)
{
    int apex = peak.pos;

    int ii = apex - 1;
    int jj = apex + 1;
    int lb = ii;  // left bound
    int rb = jj;  // right bound

    // this will be used to calculate spline area
    int slb = ii;  // spline left bound
    int srb = jj;  // spline right bound

    size_t N = intensity.size();
    if (N == 0)
        return;
    if (!spline)
        return;
    if (!baseline)
        return;

    int directionality = 0;
    float lastValue = spline[apex];
    while (ii > 0 && ii < (int)N) {
        // walk left
        float relSlope = (spline[ii] - lastValue) / lastValue;
        relSlope > 0.01 ? directionality++ : directionality = 0;
        if (intensity[ii] <= intensity[lb])
            lb = ii;
        if (spline[ii] <= spline[slb])
            slb = ii;
        if (spline[ii] == 0)
            break;
        if (spline[ii] <= baseline[ii])
            break;
        if (spline[ii] <= spline[apex] * 0.01)
            break;

        if (directionality >= 2)
            break;
        lastValue = spline[ii];
        ii = ii - 1;
    }

    directionality = 0;
    lastValue = spline[apex];

    while (jj > 0 && jj < (int)N) {
        // walk right
        float relSlope = (spline[jj] - lastValue) / lastValue;
        relSlope > 0.01 ? directionality++ : directionality = 0;
        if (intensity[jj] <= intensity[rb])
            rb = jj;
        if (spline[jj] <= spline[srb])
            srb = jj;
        if (spline[jj] == 0)
            break;
        if (spline[jj] <= baseline[ii])
            break;
        if (spline[jj] <= spline[apex] * 0.01)
            break;

        if (directionality >= 2)
            break;
        lastValue = spline[jj];
        jj = jj + 1;
    }

    // find maximum point in the span from min to max position
    for (size_t k = lb; k < rb && k < N; k++) {
        if (intensity[k] > intensity[peak.pos] && mz[k] > 0)
            peak.pos = k;
    }

    // remove zero intensity points on the left
    for (size_t k = lb; k < peak.pos && k < N; k++) {
        if (intensity[k] > 0)
            break;
        lb = k;
    }

    // remove zero intensity points on the right
    for (unsigned int k = rb; k > peak.pos && k < N; k--) {
        if (intensity[k] > 0)
            break;
        rb = k;
    }

    // for rare cases where peak is a single observation
    if (lb == apex && lb - 1 >= 0)
        lb = apex - 1;
    if (rb == apex && rb + 1 < N)
        rb = apex + 1;

    peak.minpos = lb;
    peak.maxpos = rb;
    peak.splineminpos = slb;
    peak.splinemaxpos = srb;
}

void EIC::getPeakDetails(Peak& peak)
{
    unsigned int N = intensity.size();
    if (N == 0)
        return;
    if (baseline == NULL)
        return;
    if (peak.pos >= N)
        return;

    // intensity and mz at the apex of the peaks
    peak.peakIntensity = intensity[peak.pos];
    peak.noNoiseObs = 0;
    peak.peakAreaCorrected = 0;
    peak.peakArea = 0;
    peak.peakSplineArea = 0;
    float baselineArea = 0;

    if (_sample != NULL && _sample->isBlank) {
        peak.fromBlankSample = true;
    }

    StatisticsVector<float> allmzs;
    string bitstring;
    if (peak.maxpos >= N)
        peak.maxpos = N - 1;
    if (peak.minpos >= N)
        peak.minpos = peak.pos;  // unsigned number weirdness.

    for (unsigned int i = peak.splineminpos; i <= peak.splinemaxpos; i++) {
        if (peak.splineminpos == 0 && peak.splinemaxpos == 0)
            break;
        peak.peakSplineArea += spline[i];
    }

    float lastValue = intensity[peak.minpos];
    for (unsigned int j = peak.minpos; j <= peak.maxpos; j++) {
        peak.peakArea += intensity[j];
        baselineArea += baseline[j];
        if (intensity[j] > baseline[j])
            peak.noNoiseObs++;

        if (peak.peakIntensity < intensity[j]) {
            peak.peakIntensity = intensity[j];
            peak.pos = j;
        }

        if (mz.size() > 0 && mz[j] > 0)
            allmzs.push_back(mz[j]);

        if (intensity[j] <= baseline[j]) {
            bitstring += "0";
        } else if (intensity[j] > lastValue) {
            bitstring += "+";
        } else if (intensity[j] < lastValue) {
            bitstring += "-";
        } else if (intensity[j] == lastValue) {
            if (bitstring.length() > 1)
                bitstring += bitstring[bitstring.length() - 1];
            else
                bitstring += "0";
        }

        lastValue = intensity[j];
    }

    getPeakWidth(peak);

    if (rt.size() > 0 && rt.size() == N) {
        peak.rt = rt[peak.pos];
        peak.rtmin = rt[peak.minpos];
        peak.rtmax = rt[peak.maxpos];
    }

    if (scannum.size() && scannum.size() == N) {
        peak.scan = scannum[peak.pos];  // scan number at the apex of the peak
        peak.minscan = scannum[peak.minpos];  // scan number at left most bound
        peak.maxscan =
            scannum[peak.maxpos];  // scan number at the right most bound
    }

    int n = 1;
    peak.peakAreaTop = intensity[peak.pos];
    peak.peakAreaTopCorrected = intensity[peak.pos] - baseline[peak.pos];
    if (peak.pos - 1 < N) {
        peak.peakAreaTop += intensity[peak.pos - 1];
        peak.peakAreaTopCorrected +=
            intensity[peak.pos - 1] - baseline[peak.pos - 1];
        n++;
    }
    if (peak.pos + 1 < N) {
        peak.peakAreaTop += intensity[peak.pos + 1];
        peak.peakAreaTopCorrected +=
            intensity[peak.pos + 1] - baseline[peak.pos + 1];
        n++;
    }
    peak.peakAreaTop /= n;
    peak.peakAreaTopCorrected /= n;
    if (peak.peakAreaTopCorrected < 0)
        peak.peakAreaTopCorrected = 0;

    peak.peakMz = mz[peak.pos];

    float maxBaseLine =
        MAX(MAX(baseline[peak.pos], 10),
            MAX(intensity[peak.minpos], intensity[peak.maxpos]));
    peak.peakBaseLineLevel = baseline[peak.pos];
    peak.noNoiseFraction = (float)peak.noNoiseObs / (this->_eic_noNoiseObs + 1);
    peak.peakAreaCorrected = peak.peakArea - baselineArea;
    if (peak.peakAreaCorrected < 0)
        peak.peakAreaCorrected = 0;
    peak.peakAreaFractional = peak.peakAreaCorrected / (_totalIntensity + 1);
    peak.signalBaselineRatio = peak.peakIntensity / maxBaseLine;
    peak.signalBaselineDifference = peak.peakIntensity - maxBaseLine;

    if (allmzs.size() > 0) {
        peak.medianMz = allmzs.median();
        peak.baseMz = allmzs.mean();
        peak.mzmin = allmzs.minimum();
        peak.mzmax = allmzs.maximum();
    }

    if (peak.medianMz == 0) {
        peak.medianMz = peak.peakMz;
    }

    mzPattern p(bitstring);
    if (peak.width >= 5)
        peak.symmetry = p.longestSymmetry('+', '-');
    checkGaussianFit(peak);
}

void EIC::getPeakWidth(Peak& peak)
{
    int width = 1;
    int left = 0;
    int right = 0;
    size_t N = intensity.size();

    for (unsigned int i = peak.pos - 1; i > peak.minpos && i < N; i--) {
        if (intensity[i] < baseline[i]
            || mzUtils::almostEqual(baseline[i], intensity[i], 0.0000005f))
            break;
        else
            left++;
    }
    for (unsigned int j = peak.pos + 1; j < peak.maxpos && j < N; j++) {
        if (intensity[j] > baseline[j])
            right++;
        else
            break;
    }
    peak.width = width + left + right;
}

void EIC::filterPeaks()
{
    unsigned int i = 0;
    while (i < peaks.size()) {
        if (_filterSignalBaselineDiff > peaks[i].signalBaselineDifference) {
            peaks.erase(peaks.begin() + i);
        } else {
            ++i;
        }
    }
}

vector<mzPoint> EIC::getIntensityVector(Peak& peak)
{
    vector<mzPoint> y;

    if (intensity.size() > 0) {
        unsigned int maxi = peak.maxpos;
        unsigned int mini = peak.minpos;
        if (maxi >= intensity.size())
            maxi = intensity.size() - 1;

        for (unsigned int i = mini; i <= maxi; i++) {
            // TODO all intensity points are being pushed
            if (baseline and intensity[i] > baseline[i]) {
                y.push_back(mzPoint(rt[i], intensity[i], mz[i]));
            } else {
                y.push_back(mzPoint(rt[i], intensity[i], mz[i]));
            }
        }
    }
    return y;
}

void EIC::checkGaussianFit(Peak& peak)
{
    peak.gaussFitSigma = 0;
    peak.gaussFitR2 = 0.03;
    int left = peak.pos - peak.minpos;
    int right = peak.maxpos - peak.pos;
    if (left <= 0 || right <= 0)
        return;
    int moves = min(left, right);
    if (moves < 3)
        return;

    // copy intensities into separate vector
    // dim
    vector<float> pints(moves * 2 + 1);

    size_t j = peak.pos + moves;
    if (j >= intensity.size())
        j = intensity.size() - 1;
    if (j < 1)
        j = 1;
    size_t i = peak.pos - moves;
    if (i < 1)
        i = 1;

    int k = 0;
    for (; i <= j; i++) {
        pints[k] = intensity[i];
        k++;
    }
    mzUtils::gaussFit(pints, &(peak.gaussFitSigma), &(peak.gaussFitR2));
}

void EIC::getPeakStatistics()
{
    for (unsigned int i = 0; i < peaks.size(); i++) {
        findPeakBounds(peaks[i]);
        getPeakDetails(peaks[i]);

        if (peaks[i].peakAreaTopCorrected > _maxAreaTopIntensity)
            _maxAreaTopIntensity = peaks[i].peakAreaTopCorrected;

        if (peaks[i].peakAreaTop > _maxAreaTopNotCorrectedIntensity)
            _maxAreaTopNotCorrectedIntensity = peaks[i].peakAreaTop;

        if (peaks[i].peakAreaCorrected > _maxAreaIntensity)
            _maxAreaIntensity = peaks[i].peakAreaCorrected;

        if (peaks[i].peakArea > _maxAreaNotCorrectedIntensity)
            _maxAreaNotCorrectedIntensity = peaks[i].peakArea;
    }

    // assign peak ranks based on total area of the peak
    sort(peaks.begin(), peaks.end(), Peak::compArea);
    for (unsigned int i = 0; i < peaks.size(); i++)
        peaks[i].peakRank = i;
}

void EIC::summary()
{
    cerr << "EIC: mz=" << _mzmin << "-" << _mzmax << " rt=" << _rtmin << "-"
         << _rtmax << endl;
    cerr << "   : maxIntensity=" << _maxIntensity << endl;
    cerr << "   : peaks=" << peaks.size() << endl;
}

void EIC::removeLowRankGroups(vector<PeakGroup>& groups, unsigned int rankLimit)
{
    // Merged to 776
    if (groups.size() < rankLimit)
        return;
    std::sort(groups.begin(), groups.end(), PeakGroup::compIntensity);
    for (unsigned int i = 0; i < groups.size(); i++) {
        if (i > rankLimit) {
            groups.erase(groups.begin() + i);
            i--;
        }
    }
}

// TODO: Lots of parameters. Refactor this code - Sahil
vector<PeakGroup> EIC::groupPeaks(vector<EIC*>& eics,
                                  mzSlice* slice,
                                  int smoothingWindow,
                                  float maxRtDiff,
                                  double minQuality,
                                  double distXWeight,
                                  double distYWeight,
                                  double overlapWeight,
                                  bool useOverlap,
                                  double minSignalBaselineDifference,
                                  float productPpmTolerance,
                                  string scoringAlgo)
{
    vector<mzSample*> samples;
    for (size_t i = 0; i < eics.size(); ++i) {
        samples.push_back(
            eics[i]->_sample);  // collect all mzSample into vector samples
    }
    // list filled and return by this function
    vector<PeakGroup> pgroups;

    // create EIC compose from all sample eics
    EIC* m = EIC::eicMerge(eics);
    if (!m)
        return pgroups;

    // find peaks in merged eic
    m->setFilterSignalBaselineDiff(minSignalBaselineDifference);
    m->getPeakPositions(smoothingWindow);
    sort(m->peaks.begin(), m->peaks.end(), Peak::compRt);

    for (unsigned int i = 0; i < m->peaks.size(); i++) {
        PeakGroup grp;
        grp.groupId = i;
        if (slice) {
            grp.setSlice(*slice);
            grp.setAdduct(slice->adduct);
        }
        grp.setSelectedSamples(samples);
        pgroups.push_back(grp);
    }

    for (unsigned int i = 0; i < eics.size(); i++) {  // for every sample
        for (unsigned int j = 0; j < eics[i]->peaks.size();
             j++) {  // for every peak in the sample
            Peak& b = eics[i]->peaks[j];
            b.groupNum = -1;
            b.groupOverlap = FLT_MIN;

            vector<Peak>::iterator itr = lower_bound(
                m->peaks.begin(), m->peaks.end(), b, Peak::compRtMin);
            int lb = (itr - (m->peaks.begin())) - 1;
            if (lb < 0)
                lb = 0;
            // Find best matching group
            for (unsigned int k = 0; k < m->peaks.size(); k++) {
                Peak& a = m->peaks[k];

                float score;

                float overlap = checkOverlap(
                    a.rtmin, a.rtmax, b.rtmin, b.rtmax);  // check for overlap
                float distx = abs(b.rt - a.rt);
                float disty = abs(b.peakIntensity - a.peakIntensity);

                if (useOverlap) {
                    if (overlap == 0 and a.rtmax < b.rtmin)
                        continue;
                    if (overlap == 0 and a.rtmin > b.rtmax)
                        break;

                    if (distx > maxRtDiff && overlap < 0.2)
                        continue;

                    score = 1.0 / (distXWeight * distx + 0.01)
                            / (distYWeight * disty + 0.01)
                            * (overlapWeight * overlap);
                } else {
                    if (distx > maxRtDiff)
                        continue;

                    score = 1.0 / (distXWeight * distx + 0.01)
                            / (distYWeight * disty + 0.01);
                }

                if (score > b.groupOverlap) {
                    b.groupNum = k;
                    b.groupOverlap = score;
                }
            }

            if (b.groupNum != -1) {
                PeakGroup& bestPeakGroup = pgroups[b.groupNum];
                bestPeakGroup.addPeak(b);
            } else {
                PeakGroup grp;
                pgroups.push_back(grp);
                grp.groupId = pgroups.size() + 1;
                grp.setSlice(*slice);
                grp.addPeak(b);
                b.groupOverlap = 0;
            }
        }
    }

    // clean up peakgroup such that there is only one peak for each sample
    // does the same funtion of vector::erase(), but much faster
    pgroups.erase(std::remove_if(pgroups.begin(),
                                 pgroups.end(),
                                 [](const PeakGroup& grp) {
                                     return grp.peaks.size() <= 0;
                                 }),
                  pgroups.end());

    for (unsigned int i = 0; i < pgroups.size(); i++) {
        PeakGroup& grp = pgroups[i];
        grp.minQuality = minQuality;
        grp.reduce();
        // Feng note: fillInPeaks is unecessary
        grp.groupStatistics();
        grp.computeAvgBlankArea(eics);
        if (grp.getFragmentationEvents().size()) {
            grp.computeFragPattern(productPpmTolerance);
            grp.matchFragmentation(productPpmTolerance, scoringAlgo);
        }
    }

    if (m)
        delete (m);
    return (pgroups);
}

void EIC::getRTMinMaxPerScan()
{
    if (this->rt.size() > 0) {
        this->_rtmin = this->rt[0];
        this->_rtmax = this->rt[this->size() - 1];
    }
}

/**
 * This is the functon which gets the EIC of the given scan for the
 * given mzmin and mzmax. This function will go through the each scan
 * and find the the max intensity and mz corresponding to that max intensity.
 * Total intensity is calculated by adding the maxintensity from each scan.
 * @param[in] scan This is the
 */
bool EIC::makeEICSlice(mzSample* sample,
                       float mzmin,
                       float mzmax,
                       float rtmin,
                       float rtmax,
                       int mslevel,
                       int eicType,
                       string filterline)
{
    float eicMz = 0, eicIntensity = 0;
    int lb, scanNum;
    vector<float>::iterator mzItr;
    deque<Scan*>::iterator scanItr;
    deque<Scan*> scans;

    scans = sample->scans;
    // binary search rt domain iterator
    Scan tmpScan(sample, 0, 1, rtmin - 0.1, 0, -1);
    scanItr = lower_bound(scans.begin(), scans.end(), &tmpScan, Scan::compRt);
    if (scanItr >= scans.end()) {
        return false;
    }

    int estimatedScans = scans.size();

    // TODO: why is 10 added?
    if (sample->maxRt - sample->minRt > 0
        && (rtmax - rtmin) / (sample->maxRt - sample->minRt) <= 1) {
        estimatedScans = float(rtmax - rtmin) / (sample->maxRt - sample->minRt)
                             * scans.size()
                         + 10;
    }

    this->scannum.reserve(estimatedScans);
    this->rt.reserve(estimatedScans);
    this->intensity.reserve(estimatedScans);
    this->mz.reserve(estimatedScans);

    scanNum = scanItr - scans.begin() - 1;

    for (; scanItr != scans.end(); scanItr++) {
        Scan* scan = *(scanItr);
        scanNum++;

        if (!(scan->filterLine == filterline || filterline == ""))
            continue;
        if (scan->mslevel != mslevel)
            continue;
        if (scan->rt < rtmin)
            continue;
        if (scan->rt > rtmax)
            break;

        eicMz = 0;
        eicIntensity = 0;

        // binary search
        mzItr = lower_bound(scan->mz.begin(), scan->mz.end(), mzmin);
        lb = mzItr - scan->mz.begin();

        switch ((EIC::EicType)eicType) {
        // takes the maximum intensity for given m/z range in a scan
        case EIC::MAX: {
            for (unsigned int scanIdx = lb; scanIdx < scan->nobs(); scanIdx++) {
                if (scan->mz[scanIdx] < mzmin)
                    continue;
                if (scan->mz[scanIdx] > mzmax)
                    break;

                if (scan->intensity[scanIdx] > eicIntensity) {
                    eicIntensity = scan->intensity[scanIdx];
                    eicMz = scan->mz[scanIdx];
                }
            }
            break;
        }

        // takes the sum of all intensities for given m/z range in a scan
        // associated m/z is the weighted average(with intensities as weights)
        case EIC::SUM: {
            float n = 0;
            for (unsigned int scanIdx = lb; scanIdx < scan->nobs(); scanIdx++) {
                if (scan->mz[scanIdx] < mzmin)
                    continue;
                if (scan->mz[scanIdx] > mzmax)
                    break;

                eicIntensity += scan->intensity[scanIdx];
                eicMz += scan->mz[scanIdx] * scan->intensity[scanIdx];
                n += scan->intensity[scanIdx];
            }
            eicMz /= n;
            break;
        }

        default: {
            for (unsigned int scanIdx = lb; scanIdx < scan->nobs(); scanIdx++) {
                if (scan->mz[scanIdx] < mzmin)
                    continue;
                if (scan->mz[scanIdx] > mzmax)
                    break;

                if (scan->intensity[scanIdx] > eicIntensity) {
                    eicIntensity = scan->intensity[scanIdx];
                    eicMz = scan->mz[scanIdx];
                }
            }
            break;
        }
        }

        this->scannum.push_back(scanNum);
        this->rt.push_back(scan->rt);
        this->intensity.push_back(eicIntensity);
        this->mz.push_back(eicMz);
        this->_totalIntensity += eicIntensity;
        if (eicIntensity > this->_maxIntensity) {
            this->_maxIntensity = eicIntensity;
            this->_rtAtMaxIntensity = scan->rt;
            this->_mzAtMaxIntensity = eicMz;
        }
    }

    return true;
}

void EIC::normalizeIntensityPerScan(float scale)
{
    if (scale != 1.0) {
        for (unsigned int j = 0; j < this->size(); j++) {
            this->intensity[j] *= scale;
        }
    }
}

/////////////////////////Test case//////////////////////////////

class EICFixture
{
    private:
    vector<mzSample*> _samples;
    MavenParameters* _mavenparameters;
    vector<PeakGroup> _allgroup;
    Peak _peak;
    Databases _database;
    EIC* _eic;
    EIC* _eicTemp;

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
    EICFixture()
    {
        mzSample* sample1 = new mzSample();
        sample1->loadSample("bin/methods/091215_120i.mzXML");
        _samples.push_back(sample1);
        mzSample* sample2 = new mzSample();
        sample1->loadSample("bin/methods/091215_120M.mzXML");
        _samples.push_back(sample2);
        mzSample* sample3 = new mzSample();
        sample1->loadSample("bin/methods/091215_240i.mzXML");
        _samples.push_back(sample3);
        mzSample* sample4 = new mzSample();
        sample1->loadSample("bin/methods/091215_240M.mzXML");
        _samples.push_back(sample4);
        _mavenparameters = new MavenParameters();

        _getGroupsFromProcessCompounds("bin/methods/KNOWNS.csv");

        _allgroup = _mavenparameters->allgroups;
        sort(_allgroup.begin(), _allgroup.end(), PeakGroup::compIntensity);
        auto peakgroup = _allgroup[0];
        _peak = peakgroup.peaks[0];
        _eic = new EIC();
        _eic = _samples[0]->getEIC(402.9929f, 402.9969f, 12.0, 16.0, 1, 0, "");
        _eicTemp = new EIC();
        _eicTemp =
            _samples[1]->getEIC(402.9929f, 402.9969f, 12.0, 16.0, 1, 0, "");
    }

    Peak peak() const
    {
        return _peak;
    }

    vector<PeakGroup> allgroup()
    {
        return _allgroup;
    }
    vector<mzSample*> samples() const
    {
        return _samples;
    }

    EIC* eic()
    {
        return _eic;
    }

    Peak peak()
    {
        return _peak;
    }

    ~EICFixture()
    {
        delete _samples[0];
        delete _mavenparameters;
    }
    EIC* eicTemp()
    {
        return _eicTemp;
    }
};

TEST_CASE_FIXTURE(EICFixture, "Testing EIC")
{
    auto mzSamples = samples();
    auto Eic = eic();
    auto peaks = peak();
    auto allgroups = allgroup();
    auto secondEic = eicTemp();
    REQUIRE(Eic);
    REQUIRE(Eic->getSmootherType() == 1);

    SUBCASE("Testing bounds of peak")
    {
        Eic->findPeakBounds(peaks);
        REQUIRE(peaks.minpos == 67);
        REQUIRE(peaks.maxpos == 130);
        REQUIRE(peaks.splineminpos == 67);
        REQUIRE(peaks.splinemaxpos == 130);
    }

    SUBCASE("Testing normalization of intensity")
    {
        vector<float> oldIntensity = Eic->intensity;
        Eic->normalizeIntensityPerScan(2);
        for (size_t i = 0; i < Eic->intensity.size(); i++) {
            REQUIRE(doctest::Approx(Eic->intensity[i]) == oldIntensity[i] * 2);
        }
    }

    SUBCASE("Testing RTmin and Rtmax")
    {
        Eic->getRTMinMaxPerScan();
        REQUIRE(doctest::Approx(Eic->rtmin()) == 12.012);
        REQUIRE(doctest::Approx(Eic->rtmax()) == 15.9873);
    }

    SUBCASE("Removing low rank groups")
    {
        REQUIRE(allgroups.size() == 51);
        Eic->removeLowRankGroups(allgroups, 10);
        REQUIRE(allgroups.size() == 11);
    }

    SUBCASE("Testing peakPosition and related fucntions")
    {
        Eic->getPeakPositions(10);

        // Test for compute baseline
        REQUIRE(Eic->baseline);

        // Test for findPeaks();
        REQUIRE(Eic->peaks.size() == 6);

        // Test for peak statistics
        REQUIRE(doctest::Approx(Eic->maxAreaTopIntensity()) == 40502.6);
        REQUIRE(doctest::Approx(Eic->maxAreaTopNotCorrectedIntensity())
                == 40502.6);
        REQUIRE(doctest::Approx(Eic->maxAreaIntensity()) == 260251);
        REQUIRE(doctest::Approx(Eic->maxAreaNotCorrectedIntensity()) == 260251);
        for (size_t i = 0; i < Eic->peaks.size(); i++)
            REQUIRE(Eic->peaks[i].peakRank == i);

        // Test for compute spline
        int n = sizeof(Eic->spline) / sizeof(float);
        REQUIRE(n == 2);
        REQUIRE(doctest::Approx(Eic->spline[0]) == 0);
        REQUIRE(doctest::Approx(Eic->spline[1]) == 0);

        // Test for filter peaks
        REQUIRE(Eic->peaks.size() == 6);

        // Test for peak details
        Eic->getPeakDetails(peaks);
        REQUIRE(peaks.maxpos == 130);
        REQUIRE(peaks.minpos == 67);
        REQUIRE(peaks.pos == 115);
        REQUIRE(doctest::Approx(peaks.peakSplineArea) == 128891);
        REQUIRE(doctest::Approx(peaks.peakArea) == 129690);
        REQUIRE(peaks.noNoiseObs == 39);
        REQUIRE(doctest::Approx(peaks.peakIntensity) == 14238);
        REQUIRE(doctest::Approx(peaks.rt) == 13.7635);
        REQUIRE(doctest::Approx(peaks.rtmin) == 13.0324);
        REQUIRE(doctest::Approx(peaks.rtmax) == 13.992);
        REQUIRE(peaks.scan == 2549);
        REQUIRE(peaks.minscan == 2501);
        REQUIRE(peaks.maxscan == 2564);
        REQUIRE(doctest::Approx(peaks.peakAreaTop) == 12839.4);
        REQUIRE(doctest::Approx(peaks.peakAreaTopCorrected) == 12839.4);
        REQUIRE(peaks.peakBaseLineLevel == 0);
        REQUIRE(doctest::Approx(peaks.noNoiseFraction) == 0.378641);
        REQUIRE(doctest::Approx(peaks.peakAreaCorrected) == 129691);
        REQUIRE(doctest::Approx(peaks.peakAreaFractional) == 0.3245);
        REQUIRE(doctest::Approx(peaks.signalBaselineRatio) == 10.9127);
        REQUIRE(doctest::Approx(peaks.signalBaselineDifference) == 12933.3);
        REQUIRE(doctest::Approx(peaks.medianMz) == 402.995);
        REQUIRE(doctest::Approx(peaks.baseMz) == 402.995);
        REQUIRE(doctest::Approx(peaks.mzmin) == 402.995);
        REQUIRE(doctest::Approx(peaks.mzmin) == 402.995);
        REQUIRE(doctest::Approx(peaks.mzmax) == 402.997);
        REQUIRE(peaks.symmetry == 8);
    }

    SUBCASE("Testing make EICs slice")
    {
        auto status =
            Eic->makeEICSlice(mzSamples[0], 180.002, 180.004, 0, 2, 1, 0, "");
        REQUIRE(status == true);
    }

    SUBCASE("Merge EIC")
    {
        vector<EIC*> eics;
        eics.push_back(Eic);
        eics.push_back(secondEic);
        auto merged = Eic->eicMerge(eics);
        REQUIRE(merged);
    }

    SUBCASE("Testing group peaks")
    {
        vector<EIC*> eics;

        MavenParameters* mavenparameters = new MavenParameters();
        mavenparameters->compoundMassCutoffWindow->setMassCutoffAndType(10,
                                                                        "ppm");
        mavenparameters->samples = mzSamples;
        mavenparameters->eic_smoothingWindow = 10;
        mavenparameters->eic_smoothingAlgorithm = 1;
        mavenparameters->amuQ1 = 0.25;
        mavenparameters->amuQ3 = 0.30;
        mavenparameters->aslsBaselineMode = false;
        mavenparameters->baseline_smoothingWindow = 5;
        mavenparameters->baseline_dropTopX = 80;
        mavenparameters->grouping_maxRtWindow = 0.5;
        mavenparameters->distXWeight = 1;
        mavenparameters->distYWeight = 5;
        mavenparameters->overlapWeight = 2;
        mavenparameters->useOverlap = 0;

        mzSlice* slice = new mzSlice();
        Databases database;
        bool matchRtFlag = true;
        int ionizationMode = 1;
        float compoundRTWindow = 2;
        string loadCompoundDB = "bin/methods/KNOWNS.csv";
        database.loadCompoundCSVFile(loadCompoundDB.c_str());
        vector<Compound*> compounds = database.getCompoundsSubset("KNOWNS");
        slice->compound = compounds[4];
        slice->calculateRTMinMax(matchRtFlag, compoundRTWindow);
        slice->calculateMzMinMax(mavenparameters->compoundMassCutoffWindow,
                                 ionizationMode);

        eics = PeakDetector::pullEICs(
            slice, mavenparameters->samples, mavenparameters);
        vector<PeakGroup> peakgroups =
            Eic->groupPeaks(eics,
                            slice,
                            mavenparameters->eic_smoothingWindow,
                            mavenparameters->grouping_maxRtWindow,
                            mavenparameters->minQuality,
                            mavenparameters->distXWeight,
                            mavenparameters->distYWeight,
                            mavenparameters->overlapWeight,
                            mavenparameters->useOverlap,
                            mavenparameters->minSignalBaselineDifference,
                            mavenparameters->fragmentTolerance,
                            mavenparameters->scoringAlgo);
        REQUIRE(peakgroups.size() == 2);
        REQUIRE(doctest::Approx(peakgroups[0].meanRt) == 8.2331);
        REQUIRE(doctest::Approx(peakgroups[0].meanMz) == 89.0229);

        REQUIRE(doctest::Approx(peakgroups[1].meanRt) == 9.11757);
        REQUIRE(doctest::Approx(peakgroups[2].meanMz) == 0);

        REQUIRE(peakgroups[0].peaks.size() == 1);
        REQUIRE(peakgroups[1].peaks.size() == 1);
    }

    SUBCASE("Testing getters and setters")
    {
        EIC* eic = new EIC();
        // setter and getter for baseline mode.
        eic->setBaselineMode(EIC::BaselineMode::Threshold);
        REQUIRE(eic->baselineMode() == EIC::BaselineMode::Threshold);

        // setter and getter for smoother type.
        eic->setSmootherType(EIC::SmootherType::AVG);
        REQUIRE(eic->getSmootherType() == EIC::SmootherType::AVG);

        eic->setBaselineSmoothingWindow(10);
        eic->setBaselineDropTopX(10);
        eic->setAsLSSmoothness(1);
        eic->setAsLSAsymmetry(1);
        eic->setFilterSignalBaselineDiff(1.0);
        REQUIRE(true);

        eic->setSampleName("Sample");
        REQUIRE(eic->sampleName() == "Sample");

        eic->setSample(mzSamples[0]);
        //       REQUIRE(eic->getSample() == mzSamples[0]);

        eic->setMaxIntensity(1.1);
        REQUIRE(doctest::Approx(eic->maxIntensity()) == 1.1);

        eic->setRtAtMaxIntensity(1.2);
        REQUIRE(doctest::Approx(eic->rtAtMaxIntensity()) == 1.2);

        eic->setMzAtMaxIntensity(1.2);
        REQUIRE(doctest::Approx(eic->mzAtMaxIntensity()) == 1.2);

        eic->setMaxAreaTopIntensity(1.2);
        REQUIRE(doctest::Approx(eic->maxAreaTopIntensity()) == 1.2);

        eic->setMaxAreaIntensity(1.2);
        REQUIRE(doctest::Approx(eic->maxAreaIntensity()) == 1.2);

        eic->setMaxAreaNotCorrectedIntensity(1.2);
        REQUIRE(doctest::Approx(eic->maxAreaNotCorrectedIntensity()) == 1.2);

        eic->setMaxAreaTopNotCorrectedIntensity(1.2);
        REQUIRE(doctest::Approx(eic->maxAreaTopNotCorrectedIntensity()) == 1.2);

        eic->setTotalIntensity(1.2);
        REQUIRE(doctest::Approx(eic->totalIntensity()) == 1.2);

        eic->setEic_noNoiseObs(10);
        REQUIRE(eic->eic_noNoiseObs() == 10);

        eic->setMzmin(1.2);
        REQUIRE(doctest::Approx(eic->mzmin()) == 1.2);

        eic->setMzmax(1.2);
        REQUIRE(doctest::Approx(eic->mzmax()) == 1.2);

        eic->setRtmin(1.2);
        REQUIRE(doctest::Approx(eic->rtmin()) == 1.2);

        eic->setRtmax(1.2);
        REQUIRE(doctest::Approx(eic->rtmax()) == 1.2);
    }
}
