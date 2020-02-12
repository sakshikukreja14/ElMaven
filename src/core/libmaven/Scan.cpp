#include "doctest.h"
#include "Scan.h"
#include "masscutofftype.h"
#include "mzSample.h"
#include "mavenparameters.h"
#include "constants.h"
#include "SavGolSmoother.h"

Scan::Scan(mzSample* sample, int scannum, int mslevel, float rt,
           float precursorMz, int polarity)
{
    this->_sample = sample;
    this->_rt = rt;
    this->_originalRt = rt;
    this->_scannum = scannum;
    this->_precursorMz = precursorMz;
    this->_mslevel = mslevel;
    this->_polarity = polarity;
    this->_productMz = 0;
    this->_collisionEnergy = 0;
    this->_centroided = 0;
    this->_precursorCharge = 0;
    this->_precursorIntensity = 0;
    this->_isolationWindow = 1;
}

void Scan::deepcopy(Scan* b)
{
    this->_sample = b->_sample;
    this->_rt = b->_rt;
    this->_scannum = b->_scannum;
    this->_precursorMz = b->_precursorMz;
    this->_precursorIntensity = b->_precursorIntensity;
    this->_precursorCharge = b->_precursorCharge;
    this->_mslevel = b->_mslevel;
    this->_polarity = b->_polarity;
    this->_productMz = b->_productMz;
    this->_collisionEnergy = b->_collisionEnergy;
    this->_centroided = b->_centroided;
    this->intensity = b->intensity;
    this->mz    = b->mz;
    this->_scanType = b->_scanType;
    this->_filterLine = b->_filterLine;
    this->setPolarity( b->getPolarity() );
    this->_originalRt = b->_originalRt;
    this->_isolationWindow = b->_isolationWindow;

}

int Scan::findHighestIntensityPos(float _mz, MassCutoff *massCutoff)
{
    float mzmin = _mz - massCutoff->massCutoffValue(_mz);
    float mzmax = _mz + massCutoff->massCutoffValue(_mz);
    vector<float>::iterator itr = lower_bound(mz.begin(),
                                              mz.end(), mzmin-1);
    int lb = itr-mz.begin();
    int bestPos = -1;
    float highestIntensity = 0;
    for(unsigned int k = lb; k < nobs(); k++ ) {
            if (mz[k] < mzmin)
                continue;
            if (mz[k] > mzmax)
                break;
            if (intensity[k] > highestIntensity ) {
                    highestIntensity = intensity[k];
                    bestPos = k;
            }
    }
    return bestPos;
}

/*
    @author: Sahil
*/
//TODO: Sahil, Added while merging point
int Scan::findClosestHighestIntensityPos(float _mz, MassCutoff *massCutoff)
{
    float mzmin = _mz - massCutoff->getMassCutoff()-0.001;
    float mzmax = _mz + massCutoff->getMassCutoff()+0.001;

    vector<float>::iterator itr = lower_bound(mz.begin(),
                                              mz.end(),
                                              mzmin-0.1);
    int lb = itr - mz.begin();
    float highestIntensity=0;
    for(unsigned int k=lb; k < mz.size(); k++ ) {
            if (mz[k] < mzmin)
                continue;
            if (mz[k] > mzmax)
                break;
            if (intensity[k] > highestIntensity)
                highestIntensity = intensity[k];
    }

    int bestPos=-1; float bestScore=0;
    for(unsigned int k=lb; k < mz.size(); k++ ) {
            if (mz[k] < mzmin) continue;
            if (mz[k] > mzmax) break;
            float deltaMz = (mz[k] - _mz);
            float alignScore = sqrt(intensity[k] /
                               highestIntensity) -
                               (deltaMz*deltaMz);
            if (bestScore < alignScore){
                bestScore=alignScore;
                bestPos=k;
            }
    }

    return bestPos;
}

vector<int> Scan::findMatchingMzs(float mzmin, float mzmax)
{
    vector<int>matches;
    vector<float>::iterator itr = lower_bound(mz.begin(),
                                              mz.end(), mzmin-1);
    int lb = itr - mz.begin();
    for(unsigned int k=lb; k < nobs(); k++ ) {
            if (mz[k] < mzmin)
                continue;
            if (mz[k] > mzmax)
                break;
            matches.push_back(k);
    }
    return matches;
}

void Scan::quantileFilter(int minQuantile)
{
        if (intensity.size() == 0 ) return;
        if( minQuantile <= 0 || minQuantile >= 100 ) return;

        int vsize = intensity.size();
        vector<float>dist = quantileDistribution(this->intensity);
        vector<float>cMz;
        vector<float>cIntensity;
        for(int i = 0; i < vsize; i++ ) {
            if ( intensity[i] > dist[ minQuantile ]) {
                cMz.push_back(mz[i]);
                cIntensity.push_back(intensity[i]);
            }
        }
        vector<float>(cMz).swap(cMz);
        vector<float>(cIntensity).swap(cIntensity);
        mz.swap(cMz);
        intensity.swap(cIntensity);
}

void Scan::intensityFilter(int minIntensity)
{
        if (intensity.size() == 0 ) return;
        //first pass.. find local maxima in intensity space
        int vsize = intensity.size();
        vector<float>cMz;
        vector<float>cIntensity;
        for(int i = 0; i < vsize; i++ ) {
           if ( intensity[i] > minIntensity) { //local maxima
                cMz.push_back(mz[i]);
                cIntensity.push_back(intensity[i]);
            }
        }
        vector<float>(cMz).swap(cMz);
        vector<float>(cIntensity).swap(cIntensity);
        mz.swap(cMz);
        intensity.swap(cIntensity);
}

void Scan::simpleCentroid()
{
        if (intensity.size() < 5 ) return;
        vector<float> spline=smoothenIntensitites();
        //find local maxima in intensity space
        int vsize = spline.size();
        vector<float>cMz;
        vector<float>cIntensity;

        findLocalMaximaInIntensitySpace(vsize, &cMz, &cIntensity, &spline);

        updateIntensityWithTheLocalMaximas(&cMz, &cIntensity);

        _centroided = true;
}

vector<float> Scan::smoothenIntensitites()
{
    //pass zero smooth..
    int smoothWindow = intensity.size() / 20;
    int order=2;

    if (smoothWindow < 1 )  { smoothWindow = 2; }
    if (smoothWindow > 10 ) { smoothWindow = 10; }

    mzUtils::SavGolSmoother smoother(smoothWindow,smoothWindow,order);
    //smooth once
    vector<float>spline = smoother.Smooth(intensity);
    //smooth twice
    spline = smoother.Smooth(spline);

    return spline;
}

void Scan::findLocalMaximaInIntensitySpace(int vsize, vector<float> *cMz,
                                           vector<float> *cIntensity,
                                           vector<float> *spline)
{
    for(int i = 1; i < vsize-2; i++)
    {
        if ( (*spline)[i] > (*spline)[i-1] &&
            (*spline)[i] > (*spline)[i+1] ) { //local maxima in spline space
                //local maxima in real intensity space
                float maxMz = mz[i];
                float maxIntensity = intensity[i];
                for(int j = i-1; j < i+1; j++){
                        if (intensity[i] > maxIntensity) {
                            maxIntensity = intensity[i];
                            maxMz = mz[i]; }
                }
                cMz->push_back(maxMz);
                cIntensity->push_back(maxIntensity);
        }
    }
}

void Scan::updateIntensityWithTheLocalMaximas(vector<float> *cMz,
                                              vector<float> *cIntensity)
{
    vector<float>(*cMz).swap(*cMz);
    vector<float>(*cIntensity).swap(*cIntensity);
    mz.swap(*cMz);
    intensity.swap(*cIntensity);
}

bool Scan::hasMz(float _mz, MassCutoff *massCutoff)
{
    float mzmin = _mz - massCutoff->massCutoffValue(_mz);
    float mzmax = _mz + massCutoff->massCutoffValue(_mz);

    sort(mz.begin(), mz.end());

    vector<float>::iterator itr = lower_bound(mz.begin(), mz.end(), mzmin);
    for(unsigned int k = itr - mz.begin(); k < nobs(); k++ ) {
        if (mz[k] >= mzmin && mz[k] <= mzmax )
            return true;
        if (mz[k] > mzmax )
            return false;
    }
    return false;
}

void Scan::summary()
{
    cerr << "Polarity=" << getPolarity()
         << " msLevel="  << _mslevel
         << " rt=" << _rt
         << " m/z size=" << mz.size()
         << " ints size=" << intensity.size()
         << " precursorMz=" << _precursorMz
         << " productMz=" << _productMz
         << " srmID=" << _filterLine
         << " totalIntensty=" << this->totalIntensity()
         << endl;
}

bool Scan::setParentPeakData(float mzfocus,
                             float noiseLevel,
                             MassCutoff *massCutoffMerge,
                             float minSigNoiseRatio)
{
    bool flag = true;
    int mzfocus_pos = this->findHighestIntensityPos(mzfocus,massCutoffMerge);

    if (mzfocus_pos < 0 ){
        cerr << "ERROR: Can't find parent " << mzfocus << endl;
        flag=false;
        return flag;
    }

    _parentPeakIntensity = this->intensity[mzfocus_pos];
    float parentPeakSN = _parentPeakIntensity/noiseLevel;
    if(parentPeakSN <= minSigNoiseRatio){
        flag=false;
        return flag;
    }

    return flag;
}

void Scan::initialiseBrotherData(int z, float mzfocus)
{
    //predict what M ought to be
    brotherdata->expectedMass = (mzfocus*z)-z;
    brotherdata->countMatches=0;
    brotherdata->totalIntensity=0;
    brotherdata->upCount=0;
    brotherdata->downCount=0;
    brotherdata->minZ=z;
    brotherdata->maxZ=z;
}

void Scan::updateBrotherDataIfPeakFound(int loopdirection, int ii,
                                        bool *flag, bool *lastMatched,
                                        float *lastIntensity,
                                        float noiseLevel,
                                        MassCutoff *massCutoffMerge)
{
    float brotherMz = (brotherdata->expectedMass + ii) / ii;
    int pos = this->findHighestIntensityPos(brotherMz, massCutoffMerge);
    float brotherIntensity = pos >= 0 ? this->intensity[pos] : 0;
    float snRatio = brotherIntensity / noiseLevel;
    if (brotherIntensity < 1.1*(*lastIntensity) && snRatio > 2 &&
        withinXMassCutoff(this->mz[pos]*ii-ii,brotherdata->expectedMass,
                          massCutoffMerge))
    {
        if (loopdirection==1) {
            brotherdata->maxZ = ii;
            brotherdata->upCount++;
        } else if (loopdirection==-1) {
            brotherdata->minZ = ii;
            brotherdata->downCount++;
        }
        brotherdata->countMatches++;
        brotherdata->totalIntensity += brotherIntensity;
        *lastMatched=true;
        *lastIntensity=brotherIntensity;

    } else if (*lastMatched == true) {
        *flag=false;
        return;
    }

}

void Scan::findBrotherPeaks (ChargedSpecies* x,
                            float mzfocus,
                            float noiseLevel,
                            MassCutoff *massCutoffMerge,
                            int minDeconvolutionCharge,
                            int maxDeconvolutionCharge,
                            int minDeconvolutionMass,
                            int maxDeconvolutionMass,
                            int minChargedStates)
{
    brotherdata = &b;
    for(int z = minDeconvolutionCharge; z <= maxDeconvolutionCharge; z++ ) {

        initialiseBrotherData(z,mzfocus);
        if (brotherdata->expectedMass >= maxDeconvolutionMass ||
            brotherdata->expectedMass <= minDeconvolutionMass )
            continue;
        bool flag = true;
        bool lastMatched = false;
        int loopdirection;
        loopdirection = 1;
        float lastIntensity = _parentPeakIntensity;
        for(int ii=z; ii < z+50 && ii<maxDeconvolutionCharge; ii++ ) {
            updateBrotherDataIfPeakFound(loopdirection, ii, &flag,
                                         &lastMatched, &lastIntensity,
                                         noiseLevel, massCutoffMerge);
            if (flag==false)
               break;
        }
        flag = true;
        lastMatched = false;
        loopdirection = -1;
        lastIntensity = _parentPeakIntensity;
        for(int ii=z-1; ii > z-50 && ii>minDeconvolutionCharge; ii--) {
             updateBrotherDataIfPeakFound(loopdirection, ii, &flag,
                                         &lastMatched, &lastIntensity,
                                         noiseLevel, massCutoffMerge);
             if (flag==false)
                 break;
        }
        updateChargedSpeciesDataAndFindQScore(x, z, mzfocus, noiseLevel,
                                             massCutoffMerge, minChargedStates);

    }
}


void Scan::updateChargedSpeciesDataAndFindQScore(ChargedSpecies* x, int z,
                                                 float mzfocus,
                                                 float noiseLevel,
                                                 MassCutoff *massCutoffMerge,
                                                 int minChargedStates)
{
    if (x->totalIntensity < brotherdata->totalIntensity &&
        brotherdata->countMatches>minChargedStates &&
        brotherdata->upCount >= 2 &&
        brotherdata->downCount >= 2 ) {
            x->totalIntensity = brotherdata->totalIntensity;
            x->countMatches = brotherdata->countMatches;
            x->deconvolutedMass = (mzfocus*z)-z;
            x->minZ = brotherdata->minZ;
            x->maxZ = brotherdata->maxZ;
            x->scan = this;
            x->observedCharges.clear();
            x->observedMzs.clear();
            x->observedIntensities.clear();
            x->upCount = brotherdata->upCount;
            x->downCount = brotherdata->downCount;

            float qscore=0;
            for(int ii=brotherdata->minZ; ii <= brotherdata->maxZ; ii++ ) {
                    int pos = this->findHighestIntensityPos(
                                          (brotherdata->expectedMass + ii) / ii,
                                           massCutoffMerge );
                    if (pos > 0 ) {
                        x->observedCharges.push_back(ii);
                        x->observedMzs.push_back( this->mz[pos] );
                        x->observedIntensities.push_back(this->intensity[pos]);
                        float snRatio = this->intensity[pos]/noiseLevel;
                        qscore += log(pow(0.97,(int)snRatio));
                    }
            }
            x->qscore = -20*qscore;
    }
}

ChargedSpecies* Scan::deconvolute(float mzfocus, float noiseLevel,
                                  MassCutoff *massCutoffMerge,
                                  float minSigNoiseRatio,
                                  int minDeconvolutionCharge,
                                  int maxDeconvolutionCharge,
                                  int minDeconvolutionMass,
                                  int maxDeconvolutionMass,
                                  int minChargedStates )
{
    bool flag=setParentPeakData(mzfocus,noiseLevel,
                                  massCutoffMerge,minSigNoiseRatio);
    if (flag==false)
        return NULL;

    int scanTotalIntensity=0;
    for(unsigned int i= 0; i < this->nobs(); i++)
        scanTotalIntensity += this->intensity[i];

    ChargedSpecies* x = new ChargedSpecies();
    findBrotherPeaks (x, mzfocus, noiseLevel,
                     massCutoffMerge, minDeconvolutionCharge,
                     maxDeconvolutionCharge, minDeconvolutionMass,
                     maxDeconvolutionMass, minChargedStates);

    if ( x->countMatches > minChargedStates ) {
            findError(x);
            return x;
    } else {
            delete(x);
            x=NULL;
            return(x);
    }
}

void Scan::findError(ChargedSpecies* x)
{
    float totalError=0; brotherdata->totalIntensity=0;
    for(unsigned int i=0; i < x->observedCharges.size(); i++ ) {
        float My = (x->observedMzs[i] * x->observedCharges[i]) -
                   x->observedCharges[i];
        float deltaM = abs(x->deconvolutedMass - My);
        totalError += deltaM * deltaM;
        brotherdata->totalIntensity += x->observedIntensities[i];
    }

    x->error = sqrt(totalError / x->countMatches);
}

vector<int> Scan::intensityOrderDesc()
{
    vector<pair<float,int> > mzarray(nobs());
    vector<int>position(nobs());

    for(unsigned int pos=0; pos < nobs(); pos++ ) {
        mzarray[pos] = make_pair(intensity[pos],pos);
    }

   //reverse sort first key [ ie intensity ]
   sort(mzarray.rbegin(), mzarray.rend());
   //return positions in order from highest to lowest intenisty
   for(unsigned int i=0; i < mzarray.size(); i++){
       position[i] = mzarray[i].second;
   }
   return position;
}

vector <pair<float,float> > Scan::getTopPeaks(float minFracCutoff,
                                             float minSNRatio = 3,
                                             int dropTopX = 40)
{
    vector<pair<float,float>> selected;
    if (nobs() == 0)
        return selected;

    unsigned int N = nobs();

    vector<int> positions = this->intensityOrderDesc();
    float maxI = intensity[positions[0]];

   //compute baseline intensity.. 
   float cutvalueF = (100.0-(float) dropTopX)/101;
   float baseline=1; 
   unsigned int mid = N * cutvalueF;
   if(mid < N) baseline = intensity[positions[mid]];

   for(unsigned int i=0; i<N; i++) {
       int pos = positions[i];
       if (intensity[pos] / baseline > minSNRatio &&
           intensity[pos]/maxI > minFracCutoff) {
            selected.push_back( make_pair(intensity[pos], mz[pos]));
       } else {
            break;
       }
   }
    return selected;
}

string Scan::toMGF()
{
    //Merged with Maven776 - Kiran
    std::stringstream buffer;
    buffer << "BEGIN IONS" << endl;
    if (_sample){
        buffer << "TITLE="
               <<  _sample->sampleName
               << "."
               << _scannum
               << "."
               << _scannum
               << "."
               << _precursorCharge
               << endl;
    }
    buffer << "PEPMASS="
           << setprecision(8)
           << _precursorMz
           << " "
           << setprecision(3)
           << _precursorIntensity
           << endl;

    buffer << "RTINSECONDS="
           << setprecision(9)
           << _rt*60
           << endl;

    buffer << "CHARGE="
           << _precursorCharge;
    if(_polarity < 0){
        buffer << "-";
    } else
        buffer << "+";

    buffer << endl;

    for(unsigned int i=0; i < mz.size(); i++) {
        buffer << setprecision(8)
               << mz[i]
               << " "
               << setprecision(3)
               << intensity[i]
               << endl;
    }

    buffer << "END IONS" << endl;
    return buffer.str();
}

vector<int> Scan::assignCharges(MassCutoff *massCutoffTolr)
{
    if ( nobs() == 0) {
        vector<int>empty;
        return empty;
    }

    //current scan size
    int N = nobs();
    vector<int>chargeStates (N,0);
    vector<int>peakClusters = vector<int>(N,0);
    vector<int>parentPeaks = vector<int>(N,0);
    int clusterNumber=0;

    //order intensities from high to low
    vector<int>intensityOrder = intensityOrderDesc();
    double NMASS=C13_MASS-12.00;

    //a little silly, required number of peaks in a series in already to call a charge
                          //z=0,   z=1,    z=2,   z=3,    z=4,   z=5,    z=6,     z=7,   z=8,
    int minSeriesSize[9] = { 1,     2,     3,      3,      3,     4,      4,       4,     5  } ;

    //for every position in a scan
    for(int i = 0; i < N; i++ ) {
        int pos=intensityOrder[i];
        float centerMz = mz[pos];
        float centerInts = intensity[pos];
       MassCutoff massCutoff=*massCutoffTolr;
       massCutoff.setMassCutoffAndType(2*massCutoffTolr->getMassCutoff(),
                                       massCutoffTolr->getMassCutoffType());

        if (chargeStates[pos] != 0)
            continue;   //charge already assigned

        //check for charged peak groups
        int bestZ=0;
        int maxSeriesIntenisty = 0;
        vector<int>bestSeries;

        //determine most likely charge state
        for(int z = 5; z >= 1; z--) {
            float delta = NMASS/z;
            int zSeriesIntensity = centerInts;
            vector<int>series;
            for(int j=1; j<6; j++) {
                //forward
                float mz=centerMz+(j*delta);
                int matchedPos = findHighestIntensityPos(mz,&massCutoff);
                if (matchedPos>0 && intensity[matchedPos]<centerInts) {
                    series.push_back(matchedPos);
                    zSeriesIntensity += intensity[matchedPos];
                } else break;
            }

            for(int j=1; j<3; j++) {  //back
                float mz=centerMz-(j*delta);
                int matchedPos = findHighestIntensityPos(mz,&massCutoff);
                if (matchedPos>0 && intensity[matchedPos]<centerInts) {
                    series.push_back(matchedPos);
                    zSeriesIntensity += intensity[matchedPos];
                } else break;
            }
            if (zSeriesIntensity>maxSeriesIntenisty){
                bestZ=z;
                maxSeriesIntenisty = zSeriesIntensity;
                bestSeries=series;
            }
        }

        //series with highest intensity is taken to be be the right one
        if(bestZ > 0 and static_cast<int>(bestSeries.size()) >=
                              minSeriesSize[bestZ] ) {
            clusterNumber++;
            int parentPeakPos=pos;
            for(unsigned int j=0; j<bestSeries.size();j++) {
                int brother_pos =bestSeries[j];
                if(bestZ > 1 and mz[brother_pos] < mz[parentPeakPos]
                   and intensity[brother_pos] < intensity[parentPeakPos]
                   and intensity[brother_pos] > intensity[parentPeakPos]*0.25)
                   parentPeakPos=brother_pos;

                chargeStates[brother_pos]=bestZ;
                peakClusters[brother_pos]=clusterNumber;
             }

            peakClusters[parentPeakPos]=clusterNumber;
            parentPeaks[parentPeakPos]=bestZ;
        }
    }
    return parentPeaks;
}

Scan* Scan::getLastFullScan(int historySize)
{
        if (!this->_sample)
            return 0;
    int scanNum = this->_scannum;
    for(int i = scanNum; i > (scanNum - historySize); i--) {
        Scan* lscan = this->_sample->getScan(i);
        if (!lscan or lscan->_mslevel > 1)
            continue;
        return lscan; // found ms1 scan, all is good
    }
    return 0;
}

void Scan::recalculatePrecursorMz(float ppm)
{
    if (_mslevel != 2)
        return;
    
    Scan* fullScan = getLastFullScan(50);
    if (!fullScan)
        return;
    
    MassCutoff* massCutoff = new MassCutoff();
    massCutoff->setMassCutoffAndType(ppm, "ppm");

    //find highest intensity precursor for this ms2 scan
    //increase the error range till a precursor is found
    for (int i : {1, 2, 3, 4, 5}) {
        massCutoff->setMassCutoff(ppm * i);
        unsigned int pos = fullScan->findHighestIntensityPos(this->_precursorMz,
                                                    massCutoff);
        if (pos > 0 && pos < fullScan->nobs()) {
            this->_precursorMz = fullScan->mz[pos];
            break;
        }
    }
}

vector<mzPoint> Scan::getIsolatedRegion(float isolationWindowAmu)
{
    vector<mzPoint> isolatedSegment;
    if(! this->_sample)
        return isolatedSegment;

    //find last ms1 scan or get out
    Scan* lastFullScan = this->getLastFullScan();
    if (!lastFullScan)
        return isolatedSegment;

    //no precursor information
    if (this->_precursorMz <= 0)
        return isolatedSegment;

    //extract isolated region
    float minMz = this->_precursorMz - (isolationWindowAmu / 2.0f);
    float maxMz = this->_precursorMz + (isolationWindowAmu / 2.0f);

    for(unsigned int i = 0; i < lastFullScan->nobs(); i++ ) {
            if (lastFullScan->mz[i] < minMz) continue;
            if (lastFullScan->mz[i] > maxMz) break;
            isolatedSegment.push_back(mzPoint(lastFullScan->_rt,
                                      lastFullScan->intensity[i],
                                      lastFullScan->mz[i]));
    }
    return isolatedSegment;
}

double Scan::getPrecursorPurity(float ppm)
{
    if (this->_precursorMz <= 0 ) return 0;
    if (this->_sample == 0 ) return 0;

    //extract isolated window
    vector<mzPoint> isolatedSegment = this->getIsolatedRegion(
                                                        this->_isolationWindow);
    if (isolatedSegment.size() == 0) return 0;

    //get last full scan
    Scan* lastFullScan = this->getLastFullScan();
    if (!lastFullScan) return 0;

    //locate intensity of isolated mass
    MassCutoff* massCutoff = new MassCutoff();
    massCutoff->setMassCutoffAndType(ppm, "ppm");
    int pos = lastFullScan->findHighestIntensityPos(this->_precursorMz,
                                                    massCutoff);
    if (pos < 0) return 0;
    double targetInt = lastFullScan->intensity[pos];

    //calculate total intensity in isolated segment
    double totalInt = 0;
    for (mzPoint& point: isolatedSegment)
        totalInt += point.y;

    delete massCutoff;
    if (totalInt > 0) {
        return (targetInt / totalInt);
    } else {
        return 0;
    }
}



///////////////////////TestCases/////////////////////////////

TEST_SUITE("Testing Scan Class")
{
    TEST_CASE("Testing deepCopy")
    {
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        Scan* copy= new Scan();
        copy->deepcopy(scan);

        REQUIRE(copy->sample() == scan->sample());
        REQUIRE(copy->rt() == scan->rt());
        REQUIRE(copy->originalRt() == scan->originalRt());
        REQUIRE(copy->scannum() == scan->scannum());
        REQUIRE(copy->precursorMz() == scan->precursorMz());
        REQUIRE(copy->mslevel() == scan->mslevel());
        REQUIRE(copy->getPolarity() == scan->getPolarity());
        REQUIRE(copy->productMz() == scan->productMz());
        REQUIRE(copy->collisionEnergy() == scan->collisionEnergy());
        REQUIRE(copy->isCentroided() == scan->isCentroided());
        REQUIRE(copy->precursorCharge() == scan->precursorCharge());
        REQUIRE(copy->precursorIntensity() == scan->precursorIntensity());
        REQUIRE(copy->isolationWindow() == scan->isolationWindow());
    }

    TEST_CASE("Testing findHighestIntensityPosition"){
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(0.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(5.3);
        scan->intensity.push_back(6.2);
        scan->intensity.push_back(3.5);
        scan->intensity.push_back(4.9);
        int res = scan->findHighestIntensityPos(4.9, massCutoff);
        REQUIRE(res == 3);

    }

    TEST_CASE("Testing findClosestHighestIntensityPosition"){
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(0.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(5.3);
        scan->intensity.push_back(6.2);
        scan->intensity.push_back(3.5);
        scan->intensity.push_back(4.9);

        int res = scan->findHighestIntensityPos(4.9, massCutoff);
        REQUIRE(res == 3);
    }

    TEST_CASE("Testing findMatchingMzs"){

        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(5.3);
        scan->intensity.push_back(6.2);
        scan->intensity.push_back(3.5);
        scan->intensity.push_back(4.9);

        vector<int> res = scan->findMatchingMzs(3.5, 6.2);
        int cnt = 0;
        for(auto it = res.begin(); it != res.end(); it++)
            REQUIRE(*it == cnt++);
    }

    TEST_CASE("Testing QuantileFilter and IntensityFilter")
    {
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(2.1);
        scan->intensity.push_back(2.4);
        scan->intensity.push_back(9.8);
        scan->intensity.push_back(4.5);

        scan->quantileFilter(5);
        REQUIRE(scan->mz[0] == 6.2f);
        REQUIRE(scan->mz[1] == 3.5f);
        REQUIRE(scan->mz[2] == 4.9f);
        REQUIRE(scan->intensity[0] == 2.4f);
        REQUIRE(scan->intensity[1] == 9.8f);
        REQUIRE(scan->intensity[2] == 4.5f);

        scan->intensityFilter(4);
        REQUIRE(scan->mz[0] == 3.5f);
        REQUIRE(scan->mz[1] == 4.9f);
        REQUIRE(scan->intensity[0] == 9.8f);
        REQUIRE(scan->intensity[1] == 4.5f);

    }

    TEST_CASE("Testing hasMz"){
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(0.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(5.3);
        scan->intensity.push_back(6.2);
        scan->intensity.push_back(3.5);
        scan->intensity.push_back(4.9);

        bool res = scan->hasMz(3.5, massCutoff);
        REQUIRE(res == true);
    }

    TEST_CASE("Testing Number of mz and Min and Max Mz")
    {
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(5.3);
        scan->intensity.push_back(6.2);
        scan->intensity.push_back(3.5);
        scan->intensity.push_back(4.9);

        unsigned int res = scan->nobs();
        REQUIRE(res == 4);

        float min = scan->minMz();
        float max = scan->maxMz();

        REQUIRE(doctest::Approx(min) == 3.5);
        REQUIRE(doctest::Approx(max) == 6.2);
    }

    TEST_CASE("Testing Assign Charges")
    {
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(6.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(1.0);
        scan->intensity.push_back(2.0);
        scan->intensity.push_back(3.0);
        scan->intensity.push_back(0.0);

        vector<int> res = scan->assignCharges(massCutoff);
        REQUIRE(res[0] == 0);
        REQUIRE(res[1] == 0);
        REQUIRE(res[2] == 0);
        REQUIRE(res[3] == 0);
    }

    TEST_CASE("Testing Recalculate precursorMz")
    {
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(6.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(1.0);
        scan->intensity.push_back(2.0);
        scan->intensity.push_back(3.0);
        scan->intensity.push_back(0.0);

        scan->recalculatePrecursorMz(10);
        REQUIRE(doctest::Approx(scan->precursorMz()) == 2.086);
    }

    TEST_CASE("Testing toMgf")
    {
        MassCutoff* massCutoff = new MassCutoff();
        massCutoff->setMassCutoffAndType(6.3, "ppm");
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan= new Scan(sample, 5, 4, 0.01, 2.086, -1);
        scan->mz.push_back(5.3);
        scan->mz.push_back(6.2);
        scan->mz.push_back(3.5);
        scan->mz.push_back(4.9);
        scan->intensity.push_back(1.0);
        scan->intensity.push_back(2.0);
        scan->intensity.push_back(3.0);
        scan->intensity.push_back(0.0);

        string res = scan->toMGF();
        string check = "BEGIN IONS";
        check += "\n";
        check += "TITLE=091215_120i.5.5.0";
        check += "\n";
        check += "PEPMASS=2.086 0";
        check += "\n";
        check +="RTINSECONDS=0.599999964";
        check += "\n";
        check += "CHARGE=0-";
        check += "\n";
        check += "5.3000002 1";
        check += "\n";
        check += "6.1999998 2";
        check += "\n";
        check += "3.5 3";
        check += "\n";
        check += "4.9000001 0";
        check += "\n";
        check += "END IONS";
        check += "\n";

       REQUIRE(res == check);
    }

    TEST_CASE("Testing Deconvolute")
    {
        mzSample* sample  = new mzSample();
        sample->loadSample("bin/methods/091215_120i.mzXML");

        Scan* scan=new Scan (sample, 1, 2, 3.3, 4.4, 1);

        float intensityarr[34]={131.825, 190.206, 183.417, 0, 93.979, 62.354,
                                  65.6181, 0, 0, 23.841, 52.939, 33.7124,
                                  18.8579, 42.6058, 12.807, 89.702, 14.8700,
                                  14.56700, 14.58700, 13.4041, 9.2089, 15.139,
                                  10.980, 86.791, 0, 0, 0, 0, 0, 13.401, 92.589,
                                  15.139, 10.950, 86.11};
        scan->intensity.assign(intensityarr,intensityarr+34);
        float mzarr[34]={86.06158, 86.06161, 86.06158, 86.06159, 86.06158,
                           86.06161, 86.06163, 86.06163, 87.00926, 87.00927,
                           87.00927, 87.00926, 87.00926, 87.00926, 87.04557,
                           87.04557, 87.04557, 88.04099, 88.04099, 88.04099,
                           88.04099, 88.04099, 88.04096, 88.04095, 88.04093,
                           88.04098, 88.04099, 88.04099, 88.04099, 88.04099,
                           88.04099, 88.04099, 88.04096, 88.04095};

        scan->mz.assign(mzarr, mzarr+34);
        MavenParameters* mavenparameters = new MavenParameters();
        mavenparameters->massCutoffMerge->setMassCutoffAndType(100000, "ppm");
        ChargedSpecies* chargedspecies=scan->deconvolute(87, 4,
                                            mavenparameters->massCutoffMerge,
                                            2, 3, 100, 500, 2e5, 3);

        REQUIRE(doctest::Approx(chargedspecies->deconvolutedMass) == 7654);
        REQUIRE(chargedspecies->countMatches == 18);
        REQUIRE(doctest::Approx(chargedspecies->error) == 443.506);
    }

}




