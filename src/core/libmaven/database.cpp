#include "doctest.h"
#include "Compound.h"
#include "constants.h"
#include "database.h"
#include "masscutofftype.h"
#include "mgf/mgf.h"
#include "mzMassCalculator.h"
#include "mzSample.h"
#include "mzUtils.h"
#include "boost/algorithm/string.hpp"

void Database::removeDatabase(string dbName)
{
    auto iter = begin(compoundsDB);
    while (iter < end(compoundsDB)) {
        auto compound = *iter;
        if (compound->db() == dbName) {
            compoundIdenticalCount.erase(compound->id() + compound->name() + dbName);
            compoundIdNameDbMap.erase(compound->id() + compound->name() + dbName);
            iter = compoundsDB.erase(iter);
            delete compound;
        } else {
            ++iter;
        }
    }
}

bool Database::addCompound(Compound* newCompound)
{
    if(newCompound == nullptr)
        return false;

    // existing compound, change its name according to the number of
    // compounds with the same ID
    if (compoundIdenticalCount.count(newCompound->id()
                                     + newCompound->name()
                                     + newCompound->db())) {
        int loadOrder = compoundIdenticalCount.at(newCompound->id()
                                                  + newCompound->name()
                                                  + newCompound->db());

        // return false if any of the compounds having the same ID are the
        // exact same in all aspects.
        for (int i = 0; i < loadOrder; ++i) {
            string name = newCompound->name();
            if (i != 0)
                name = name + " (" + to_string(i) + ")";

            Compound* possibleCopy = compoundIdNameDbMap[newCompound->id()
                                                         + name
                                                         + newCompound->db()];
            if (possibleCopy != nullptr && *newCompound == *possibleCopy)
                return false;
        }

        auto originalName = newCompound->name();
        newCompound->setName (originalName
                            + " ("
                            + to_string(loadOrder)
                            + ")");
        compoundIdenticalCount[newCompound->id()
                               + originalName
                               + newCompound->db()] = ++loadOrder;
    } else {
        compoundIdenticalCount[newCompound->id()
                               + newCompound->name()
                               + newCompound->db()] = 1;
    }

    compoundIdNameDbMap[newCompound->id()
                        + newCompound->name()
                        + newCompound->db()] = newCompound;
    compoundsDB.push_back(newCompound);
    return true;
}


Compound* Database::findSpeciesByIdAndName(string id,
                                           string name,
                                           string dbName)
{
    if (compoundIdNameDbMap.count(id + name + dbName))
        return compoundIdNameDbMap[id + name + dbName];
    return NULL;
}

vector<Compound*> Database::findSpeciesById(string id, string dbName) {
    vector<Compound*> matches;
    for (auto compound : compoundsDB) {
        if (compound->id() == id && compound->db() == dbName) {
            matches.push_back(compound);
        }
    }

    return matches;
}

vector<Compound*> Database::findSpeciesByName(string name, string dbname) {
		vector<Compound*> set;
		for(unsigned int i=0; i < compoundsDB.size(); i++ ) {
                                if (compoundsDB[i]->name() == name && compoundsDB[i]->db() == dbname) {
					set.push_back(compoundsDB[i]);
				}
		}
		return set;
}

vector<Compound*> Database::getCompoundsSubset(string dbname) {
	vector<Compound*> subset;
	for (unsigned int i=0; i < compoundsDB.size(); i++ ) {
                        if (compoundsDB[i]->db() == dbname) {
					subset.push_back(compoundsDB[i]);
			}
	}
	return subset;
}

vector<Compound*> Database::getKnowns() {
    return getCompoundsSubset("KNOWNS");
}

map<string,int> Database::getDatabaseNames() {
	map<string,int>dbnames;
        for (unsigned int i=0; i < compoundsDB.size(); i++ ) dbnames[ compoundsDB[i]->db() ]++;
	return dbnames;
}


bool Database::_startsWith(string line, string text)
{
    transform(line.begin(), line.end(), line.begin(), ::tolower);
    transform(text.begin(), text.end(), text.begin(), ::tolower);
    for( size_t i = 0; i < text.size(); i++){
        if(line[i] != text[i])
            return false;
    }
    return true;
}

bool Database::_contain(string line, string text)
{
    transform(line.begin(), line.end(), line.begin(), ::tolower);
    transform(text.begin(), text.end(), text.begin(), ::tolower);
    if(boost::algorithm::contains(line, text))
        return true;
    else
        return false;
}


int Database::loadNISTLibrary(string filename,
                              bsignal::signal<void (string, int, int)>* signal)
{


    if (signal)
        (*signal)("Preprocessing database " + filename, 0, 0);

    cout << endl<< "Counting number of lines in NIST Libary file…" << filename;
    ifstream file(filename.c_str());
    if(!file.is_open())
        return 0;

    file.unsetf(ios_base::skipws); // do not skip newlines
    unsigned lineCount = std::count(istream_iterator<char>(file),
                                    istream_iterator<char>(),
                                    '\n');
    file.close();
    cout << endl<< "Loading NIST Libary: " << filename;

    regex whiteSpace("\\s+");
    regex formulaMatch("Formula\\=(C\\d+H\\d+\\S*)");
    regex retentionTimeMatch("AvgRt\\=(\\S+)");

    string dbName = mzUtils::cleanFilename(filename);
    Compound* currentCompound = nullptr;
    bool capturePeaks = false;
    int compoundCount = 0;
    int currentLine = 0;

    file.open(filename, ios::in);

    while (!file.eof()) {

        string line;
        getline(file,line);

        if (_startsWith(line, "NAME:")) {
            // before reading the next record or ending stream, save the
            // compound created from last record
            if (currentCompound and !currentCompound->name().empty()) {
                if (!currentCompound->formula().empty()) {
                    auto formula = currentCompound->formula();
                    auto exactMass = MassCalculator::computeMass(formula, 0);
                    currentCompound->setMass (exactMass);
                }
                bool flag = addCompound(currentCompound);

                if (flag){
                    ++compoundCount;
                }
            }

            // we need to check this again before creating a new compound,
            // otherwise it would create one at stream end as well
            if (_startsWith(line, "NAME:")) {
                // new compound
                string name = line.substr(6, line.length());
                currentCompound = new Compound(name, name, "", 0);
                currentCompound->setDb (dbName);
                capturePeaks = false;
            }
        }

        if(currentCompound == nullptr)
            continue;

        if (_startsWith(line, "MW:")) {
            currentCompound->setMass (string2float(line.substr(3, line.length())));
        } else if (_startsWith(line, "CE:")) {
            currentCompound->setCollisionEnergy (string2float(line.substr(3, line.length())));
        } else if (_startsWith(line, "ID:")){
            string id = line.substr(4, line.length());
            if (id.size() > 0)
                currentCompound->setId(id);
        } else if (_startsWith(line, "LOGP:")) {
            currentCompound->logP = string2float(line.substr(5, line.length()));
        } else if (_startsWith(line, "RT:")){
            currentCompound->setExpectedRt (string2float(line.substr(3, line.length())));
        } else if (_startsWith(line, "SMILE:")){
            string smileString = line.substr(7, line.length());
            if (smileString.size() > 0)
                currentCompound->smileString = smileString;
        } else if (_startsWith(line, "SMILES:")) {
            string smileString = line.substr(8, line.length());
            if (smileString.size() > 0)
                currentCompound->smileString = smileString;
        } else if (_startsWith(line, "PRECURSORMZ:")) {
            currentCompound->setPrecursorMz (string2float(line.substr(13, line.length())));
        } else if (_startsWith(line, "EXACTMASS:")) {
            currentCompound->setMass (string2float(line.substr(10, line.length())));
        } else if (_startsWith(line,"ADDUCT:")) {
            currentCompound->adductString = line.substr(8, line.length());
        } else if (_startsWith(line, "FORMULA:")){
            string formula = line.substr(9, line.length());
            boost::replace_all_copy(formula,"\"", "");
            if (formula.size() > 0)
                currentCompound->setFormula (formula);
        } else if (_startsWith(line, "MOLECULE FORMULA:")){
            string formula = line.substr(17, line.length());
            boost::replace_all_copy(formula, "\"", "");
            if (formula.size() > 0)
                currentCompound->setFormula(formula);
        } else if (_startsWith(line, "CATEGORY:")){
            currentCompound->category.push_back(line.substr(10, line.length()));
        } else if (_startsWith(line, "TAG:")){
            if (_contain(line, "VIRTUAL"))
                currentCompound->virtualFragmentation = true;
        } else if (_startsWith(line, "ION MODE:")
                   || _startsWith(line, "IONMODE:")
                   || _startsWith(line, "IONIZATION:")) {
            if (_contain(line, "NEG"))
                currentCompound->ionizationMode = -1;
            if (_contain(line, "POS"))
                currentCompound->ionizationMode = +1;
        } else if (_startsWith(line, "COMMENT:")) {
            smatch match;
            string comment = line.substr(8, line.length());
            if (regex_search(comment, match, formulaMatch) && match.size() > 1) {
                currentCompound->setFormula (match.str(1));
            }
            if (regex_search(comment, match, retentionTimeMatch) && match.size() > 1) {
                currentCompound->setExpectedRt (string2float(match.str(1)));
            }
        } else if (_startsWith(line, "NUM PEAKS:")
                   || _startsWith(line, "NUMPEAKS:")) {
            capturePeaks = true;
        } else if (capturePeaks) {
            vector<string> mzIntensityPair;
            boost::split(mzIntensityPair, line, boost::is_any_of(" "));
            if (mzIntensityPair.size() >= 2) {
                double mz = string2float(mzIntensityPair[0]);
                double in = string2float(mzIntensityPair[1]);
                if (mz >= 0.0 && in >= 0.0) {
                    currentCompound->fragmentMzValues.push_back(mz);
                    currentCompound->fragmentIntensities.push_back(in);

                    int fragIdx = currentCompound->fragmentMzValues.size() - 1;
                    if (mzIntensityPair.size() >= 3) {
                        currentCompound->fragmentIonTypes[fragIdx] =
                                        mzIntensityPair[2];
                    }
                }
            }
        }
        ++currentLine;
        if (signal) {
            (*signal)("Loading spectral library: " + filename,
                      currentLine,
                      lineCount);
        }
    }
    if (currentCompound and !currentCompound->name().empty()) {
        if (!currentCompound->formula().empty()) {
            auto formula = currentCompound->formula();
            auto exactMass = MassCalculator::computeMass(formula, 0);
            currentCompound->setMass (exactMass);
        }
        bool flag = addCompound(currentCompound);
        if (flag){
            ++compoundCount;
        }
    }
    return compoundCount;
}


bool Database::isSpectralLibrary(string dbName) {
    auto compounds = getCompoundsSubset(dbName);
    if (compounds.size() > 0) {
        return compounds.at(0)->type() == Compound::Type::PRM;
    }
    return false;
}



int Database::loadCompoundCSVFile(string filename) {

    ifstream myfile(filename.c_str());
    if (! myfile.is_open()) return 0;

    string line;
    int loadCount = 0, lineCount = 0;
    map<string, int> header;
    vector<string> headers;

    // reset the contents of the vector containing the names of invalid rows
    invalidRows.clear();

    //assume that files are tab delimited, unless matched ".csv", then comma delimited
    string sep="\t";
    if(filename.find(".csv") != -1 || filename.find(".CSV") != -1) sep=",";

    while (getline(myfile,line)) {
        //This is used to write commands
        if (!line.empty() && line[0] == '#') continue;

        //trim spaces on the left
        size_t found = line.find_last_not_of(" \n\r\t");
        if (found != string::npos)
            line.erase(found+1);
        else continue;

        lineCount++;

        vector<string> fields;
        mzUtils::splitNew(line, sep, fields);

        mzUtils::removeSpecialcharFromStartEnd(fields);

        //Getting the heading from the csv File
        if (lineCount == 1) {
            headers = fields;
            for(unsigned int i = 0; i < fields.size(); i++ ) {
                fields[i] = makeLowerCase(fields[i]);
                header[ fields[i] ] = i;
            }
            continue;
        }

        Compound* compound = extractCompoundfromEachLine(fields, header, loadCount, filename);

        if (compound) {
            if (addCompound(compound)) {
                loadCount++;
            }
        }
    }
    sort(compoundsDB.begin(),compoundsDB.end(), Compound::compMass);
    myfile.close();
    return loadCount;
}

Compound* Database::extractCompoundfromEachLine(vector<string>& fields, map<string, int> & header, int loadCount, string filename) {
    string id, name, formula, polarityString;
    string note;
    float rt = 0, mz = 0, charge = 0, collisionenergy = 0, precursormz = 0, productmz = 0;
    int NumOfFields = fields.size();
    vector<string> categorylist;

    string dbname = mzUtils::cleanFilename(filename);

    if (header.count("mz") && header["mz"] < NumOfFields)
        mz = string2float(fields[header["mz"]]);

    //Expected RT is given importance over RT. So if Expected RT is
    //present in the DB that will be taken to the RT field in compounds
    if (header.count("rt") && header["rt"] < NumOfFields)
        rt = string2float(fields[header["rt"]]);

    if (header.count("expectedrt") && header["expectedrt"] < NumOfFields)
        rt = string2float(fields[header["expectedrt"]]);

    //TODO: Not really a todo, just marking that I fixed a typo here
    //make sure we merge it in
    if (header.count("formula") && header["formula"] < NumOfFields)
        formula = fields[header["formula"]];

    if (header.count("id") && header["id"] < NumOfFields)
        id = fields[header["id"]];

    // Compound Field is given importance than the names field
    // compound is a better field to keep
    if (header.count("name") && header["name"] < NumOfFields)
        name = fields[header["name"]];

    if (header.count("compound") && header["compound"] < NumOfFields)
        name = fields[header["compound"]];

    if (header.count("precursormz") && header["precursormz"] < NumOfFields)
        precursormz = string2float(fields[ header["precursormz"]]);

    if (header.count("productmz") && header["productmz"] < NumOfFields)
        productmz = string2float(fields[header["productmz"]]);

    if (header.count("collisionenergy") && header["collisionenergy"] < NumOfFields)
        collisionenergy = string2float(fields[ header["collisionenergy"]]);

    if (header.count("Q1") && header["Q1"] < NumOfFields)
        precursormz = string2float(fields[ header["Q1"]]);

    if (header.count("Q3") && header["Q3"] < NumOfFields)
        productmz = string2float(fields[header["Q3"]]);

    if (header.count("CE") && header["CE"] < NumOfFields)
        collisionenergy=string2float(fields[header["CE"]]);

    if (header.count("note") && header["note"] < NumOfFields)
        note = fields[header["note"]];

    categorylist = getCategoryFromDB(fields, header);

    charge = getChargeFromDB(fields, header);
    //If Some of the imp fields are not present here is th way it will be
    //assigned

    if (id.empty() && !name.empty())
        id = name;

    if (id.empty() && name.empty())
        id = "cmpd:" + integer2string(loadCount);

    //The compound should atleast have formula so that
    //mass can be calculated from the formula
    if ( mz > 0 || !formula.empty() || precursormz > 0) {
        Compound* compound = new Compound(id,name,formula,charge);

        compound->setExpectedRt (rt);

        if (mz == 0)
            mz = MassCalculator::computeMass(formula, charge);


        compound->setMass(mz);
        compound->setDb  (dbname);
        compound->setExpectedRt(rt);
        compound->setPrecursorMz (precursormz);
        compound->setProductMz( productmz);
        compound->setCollisionEnergy (collisionenergy);
        compound->note = note;

        for(unsigned int i=0; i < categorylist.size(); i++)
            compound->category.push_back(categorylist[i]);

        return compound;
    }

    if (!name.empty())
        id = name;
    invalidRows.push_back(id);

    return NULL;

}

float Database::getChargeFromDB(vector<string>& fields, map<string, int> & header) {
    float charge = 0;
    int NumOfFields = fields.size();
    if (header.count("charge") && header["charge"] < NumOfFields)
        charge = string2float(fields[header["charge"]]);

    if ( header.count("polarity") && header["polarity"] < NumOfFields) {
        string polarityString = fields[header["polarity"]];
        if ( polarityString == "+" ) {
            charge = 1;
        } else if ( polarityString == "-" ) {
            charge = -1;
        } else  {
            charge = string2float(polarityString);
        }
    }
    return charge;
}

vector<string> Database::getCategoryFromDB(vector<string>& fields, map<string, int> & header)
{
    vector<string> categorylist;
    int NumOfFields = fields.size();
    //What is category?
    if ( header.count("category") && header["category"] < NumOfFields) {
        string catstring = fields[header["category"]];
        if (!catstring.empty()) {
            mzUtils::split(catstring,';', categorylist);
            if(categorylist.size() == 0) categorylist.push_back(catstring);
        }
    }
    return categorylist;
}



int Database::loadMascotLibrary(string filename,
                                bsignal::signal<void (string, int, int)> *signal)
{
    mgf::MgfFile mgfFile;
    mgf::Driver driver(mgfFile);
    driver.trace_parsing = false;
    driver.trace_scanning = false;

    ifstream ifs(filename.c_str());

    bool result = driver.parse_stream(ifs);

    if (signal)
        (*signal)("Reading file " + filename, 0, 0);

    if (!result) {
        std::cerr << "Error parsing data stream"
                  << std::endl;
        return 0;
    }

    for (auto specIter = begin(mgfFile); specIter != end(mgfFile); ++specIter) {
        auto charges = specIter->getCHARGE();
        int charge = 1;
        if (!charges.empty())
            charge = charges.front();

        Compound* compound = new Compound(specIter->getTITLE(),
                                          specIter->getTITLE(),
                                          "",
                                          charge);

        compound->setExpectedRt ( specIter->getRTINSECONDS().first / 60.0f);
        compound->setMass (specIter->getPEPMASS().first);
        compound->setPrecursorMz(compound->mass());
        compound->smileString = specIter->getSMILES();
        compound->ionizationMode = specIter->getIONMODE() == "Negative" ||
                                   specIter->getIONMODE() == "negative"? -1
                                                                        : 1;

        // create spectra
        vector<float> fragmentMzValues;
        vector<float> fragmentInValues;
        for (auto fragPair = specIter->begin();
             fragPair != specIter->end();
             ++fragPair) {
            fragmentMzValues.push_back(fragPair->first);
            fragmentInValues.push_back(fragPair->second);
        }
        compound->fragmentMzValues = fragmentMzValues;
        compound->fragmentIntensities = fragmentInValues;

        compound->setDb( mzUtils::cleanFilename(filename));
        addCompound(compound);

        if (signal) {
            (*signal)("Loading spectral library: " + filename,
                      (specIter - begin(mgfFile)),
                      mgfFile.size());
        }
    }
    if (signal)
        (*signal)("Finished loading " + filename, 0, 0);
    return mgfFile.size();
}



TEST_CASE("Testing database class")
{
    SUBCASE("Testing CSV Loading File")
    {
        Database db;
        string mgfFile = "tests/doctest/test_loadCSV.csv";
        int rescsv = db.loadCompoundCSVFile(mgfFile);
        REQUIRE(rescsv == 10);
        vector<Compound*> compounds = db.getCompoundsSubset("test_loadCSV");

        vector<Compound*> compoundInput;
        Compound* c1 = new Compound("HMDB00653", "cholesteryl sulfate",
                                    "C27H46O4S", 0, 17.25);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C05464", "Deoxycholic acid",
                          "C26H43NO5", 0, 16.79);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C15H28O7P2","trans_trans-farnesyl diphosphate",
                          "C00448",0, 16.74);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("HMDB00619", "Cholic acid",
                          "C24H40O5", 0, 16.69);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C00341", "Geranyl-PP",
                          "C10H20O7P2", 0, 16.46) ;
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C05463", "Taurodeoxycholic acid",
                          "C26H45NO6S", 0, 16.23);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C00725", "lipoate",
                          "C8H14O2S2", 0, 15.97);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C00356", "3-hydroxy-3-methylglutaryl-CoA-nega",
                          "C27H44N7O20P3S", 0, 15.72);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C00630", "butyryl-CoA",
                          "C25H42N7O17P3S", 0, 15.72);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        c1 = new Compound("C00083", "malonyl-CoA",
                          "C24H38N7O19P3S", 0, 15.7);
        c1->setDb("test_loadCSV");
        compoundInput.push_back(c1);
        delete(c1);

        for(size_t i = 0; i < compoundInput.size(); i++ ){
            for(size_t j = 0; j < compounds.size(); j++){
                if(compoundInput[i]->id() == compounds[j]->id())
                    REQUIRE(compoundInput[i] == compounds[j]);
            }
        }
    }





    SUBCASE("Testing Mascot Library")
    {
        Database db;
        string mgfFile = "tests/doctest/test_Mascot.mgf";
        int resMgf = db.loadMascotLibrary(mgfFile);
        REQUIRE(resMgf == 10);
        vector<Compound*> compounds = db.getCompoundsSubset("test_Mascot");

        Compound* c1 = new Compound("HMDB:HMDB04095-2361 5-Methoxytryptamine M-H",
                                    "HMDB:HMDB04095-2361 5-Methoxytryptamine M-H",
                                    "", 1);
        c1->setMass(189.103);
        float rt = -0.0166667f;
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = -1;
        c1->smileString = "COC1=CC2=C(NC=C2CCN)C=C1";
        c1->setPrecursorMz(189.103);
        c1->setDb("test_Mascot");
        double mz1[] = {143.0,144.0,173.0, 175.0,
                        188.0, 190.0};
        for(int i = 0; i < 6; i++)
            c1->fragmentMzValues.push_back(mz1[i]);
        double ion1[] = {7.605, 2.369, 100.0,
                         3.178, 87.625, 2.304};
        for(int i = 0; i < 6; i++)
            c1->fragmentIntensities.push_back(ion1[i]);
        REQUIRE(*c1 == *compounds[0]);
        delete (c1);


        c1 = new Compound("HMDB:HMDB00099-156 L-Cystathionine M+H",
                          "HMDB:HMDB00099-156 L-Cystathionine M+H",
                          "", 1);
        c1->setMass(223.075);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "N[C@@H](CCSC[C@H](N)C(O)=O)C(O)=O";
        c1->setPrecursorMz(223.075);
        c1->setDb("test_Mascot");
        double mz2[] = {56.261, 88.097,
                        133.997, 177.012,
                        222.983};
        for(int i = 0; i < 5; i++)
            c1->fragmentMzValues.push_back(mz2[i]);
        double ion2[] = {5.478,	9.877, 100.0,
                         3.318,91.975};
        for(int i = 0; i < 5; i++)
            c1->fragmentIntensities.push_back(ion2[i]);
        REQUIRE(*c1 == *compounds[1]);
        delete (c1);


        c1 = new Compound("HMDB:HMDB00021-38 Iodotyrosine M+H",
                          "HMDB:HMDB00021-38 Iodotyrosine M+H",
                          "", 1);
        c1->setMass(307.978);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "N[C@@H](CC1=CC=C(O)C(I)=C1)C(O)=O";
        c1->setPrecursorMz(307.978);
        c1->setDb("test_Mascot");
        double mz3[] = {89.975, 94.035, 106.854,
                        108.236, 118.003, 119.876,
                        134.711, 248.732, 261.798};
        for(int i = 0; i < 9; i++)
            c1->fragmentMzValues.push_back(mz3[i]);
        double ion3[] = {5.266, 47.917, 28.935, 11.285,
                         21.296, 14.525, 100.0, 8.912,
                         11.69};
        for(int i = 0; i < 9; i++)
            c1->fragmentIntensities.push_back(ion3[i]);
        REQUIRE(*c1 == *compounds[2]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00145-218 Estrone M+H",
                          "HMDB:HMDB00145-218 Estrone M+H",
                          "", 1);
        c1->setMass(271.17);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "[H][C@@]12CCC(=O)[C@@]1(C)CC[C@]1([H])C3=C(CC[C@@]21[H])C=C(O)C=C3";
        c1->setPrecursorMz(271.17);
        c1->setDb("test_Mascot");
        double mz4[] = {199.0, 271.0, 272.0};
        for(int i = 0; i < 3; i++)
            c1->fragmentMzValues.push_back(mz4[i]);
        double ion4[] = {15.265, 100.0, 23.44};
        for(int i = 0; i < 3; i++)
            c1->fragmentIntensities.push_back(ion4[i]);
        REQUIRE(*c1 == *compounds[3]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00258-449 Sucrose M+H",
                          "HMDB:HMDB00258-449 Sucrose M+H",
                          "", 1);
        c1->setMass(343.124);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O";
        c1->setPrecursorMz(343.124);
        c1->setDb("test_Mascot");
        double mz5[] = {55.239, 69.118, 73.123,
                        85.059, 97.026, 99.025,
                        108.994, 114.982, 126.974,
                        144.97, 162.881};
        for(int i = 0; i < 11; i++)
            c1->fragmentMzValues.push_back(mz5[i]);
        double ion5[] = {6.379, 11.724, 5.172,
                         100.0, 19.138, 6.983,
                         5.302, 5.129, 55.172,
                         31.552, 13.103};
        for(int i = 0; i < 11; i++)
            c1->fragmentIntensities.push_back(ion5[i]);
        REQUIRE(*c1 == *compounds[4]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00032-52 7-Dehydrocholesterol M+H",
                          "HMDB:HMDB00032-52 7-Dehydrocholesterol M+H",
                          "", 1);
        c1->setMass(385.347);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "[H][C@@]12CC[C@H]([C@H](C)CCCC(C)C)[C@@]1(C)CC[C@@]1([H])C2=CC=C2C[C@@H](O)CC[C@]12C";
        c1->setPrecursorMz(385.347);
        c1->setDb("test_Mascot");
        double mz6[] = {43.315, 55.297,
                        57.193, 69.167,
                        71.228, 81.111,
                        83.101, 85.037,
                        93.281, 95.201,
                        97.106, 99.081,
                        105.068, 107.083,
                        109.136, 110.869,
                        119.285, 121.033,
                        122.954,
                        131.261, 132.947,
                        135.008, 142.691,
                        145.095, 146.953,
                        148.851, 157.174,
                        158.173, 159.094,
                        160.211, 161.179,
                        163.264, 169.01,
                        170.939, 172.821,
                        175.124, 176.756,
                        182.901, 185.072,
                        187.18, 188.937,
                        195.309, 196.41,
                        197.409, 198.768,
                        199.65, 201.102,
                        203.109, 209.27,
                        211.246, 212.027,
                        212.855, 213.628,
                        214.932, 217.048,
                        218.774, 223.077,
                        225.154, 226.021,
                        226.934, 228.957,
                        231.354, 237.11,
                        239.015, 241.1,
                        254.22, 255.048,
                        325.076, 339.315,
                        349.47, 367.233,
                        367.999, 385.161,
                        386.021};
        for(int i = 0; i < 74; i++)
            c1->fragmentMzValues.push_back(mz6[i]);

        double ion6[] = {18.02, 21.827, 45.685,
                         70.558, 62.437, 83.756,
                         63.959, 38.579, 26.904,
                         70.558, 38.579, 17.893,
                         27.284, 40.102, 56.853,
                         25.127, 25.761, 43.655,
                         32.995, 19.543, 35.533,
                         31.091, 21.32, 36.041,
                         42.132, 26.904, 40.102,
                         22.589, 58.376, 23.223,
                         54.315, 25.888, 22.843,
                         38.071, 47.716, 30.838,
                         22.589, 22.462, 35.025,
                         36.041, 24.492, 23.604,
                         20.685, 26.777, 57.868,
                         28.68, 32.995, 20.051,
                         21.32, 34.518, 18.02,
                         56.853, 20.685, 32.995,
                         17.64, 19.67, 19.67,
                         17.513, 17.893, 51.777,
                         24.619, 22.97, 34.518,
                         22.97, 28.299, 22.462,
                         45.178, 19.797, 19.162,
                         22.081, 62.944, 17.893,
                         100.0, 32.995};
        for(int i = 0; i < 74; i++)
            c1->fragmentIntensities.push_back(ion6[i]);
        REQUIRE(*c1 == *compounds[5]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00536-758 Adenylsuccinic acid M+H",
                          "HMDB:HMDB00536-758 Adenylsuccinic acid M+H",
                          "", 1);
        c1->setMass(464.082);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)N1C=NC2=C1N=CN=C2NC(CC(O)=O)C(O)=O";
        c1->setPrecursorMz(464.082);
        c1->setDb("test_Mascot");
        double mz7[] = {251.0, 252.0,
                        253.0, 386.0,
                        387.0, 431.0,
                        463.0, 464.0,
                        465.0};
        for(int i = 0; i < 9; i++)
            c1->fragmentMzValues.push_back(mz7[i]);
        double ion7[] = {65.439, 36.631, 4.831,
                         3.04, 3.697, 2.281,
                         83.455, 100.0, 38.559};
        for(int i = 0; i < 9; i++)
            c1->fragmentIntensities.push_back(ion7[i]);
        REQUIRE(*c1 == *compounds[6]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00620-837 Glutaconic acid M-H",
                          "HMDB:HMDB00620-837 Glutaconic acid M-H",
                          "", 1);
        c1->setMass(129.019);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = -1;
        c1->smileString = "OC(=O)C\\C=C\\C(O)=O";
        c1->setPrecursorMz(129.019);
        c1->setDb("test_Mascot");
        double mz8[] = {41.317, 85.036, 128.921};
        for(int i = 0; i < 3; i++)
            c1->fragmentMzValues.push_back(mz8[i]);
        double ion8[] = {16.447, 100.0, 2.237};
        for(int i = 0; i < 3; i++)
            c1->fragmentIntensities.push_back(ion8[i]);
        REQUIRE(*c1 == *compounds[7]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00700-976 Hydroxypropionic acid M+H",
                          "HMDB:HMDB00700-976 Hydroxypropionic acid M+H",
                          "", 1);
        c1->setMass(91.0395);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "OCCC(O)=O";
        c1->setPrecursorMz(91.0395);
        c1->setDb("test_Mascot");
        double mz9[] = {54.152, 55.397,
                        55.522, 63.364,
                        72.451, 73.322,
                        73.571, 89.878,
                        91.496, 91.87};
        for(int i = 0; i < 10; i++)
            c1->fragmentMzValues.push_back(mz9[i]);
        double ion9[] = {10.812, 28.713, 9.206,
                         13.339, 25.173, 100.0,
                         88.755, 23.212, 52.287,
                         10.032};
        for(int i = 0; i < 10; i++)
            c1->fragmentIntensities.push_back(ion9[i]);
        REQUIRE(*c1 == *compounds[8]);
        delete (c1);

        c1 = new Compound("HMDB:HMDB00705-984 Hexanoylcarnitine M+H",
                          "HMDB:HMDB00705-984 Hexanoylcarnitine M+H",
                          "", 1);
        c1->setMass(260.186);
        c1->setExpectedRt(-1/60.0f);
        c1->ionizationMode = 1;
        c1->smileString = "CCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C";
        c1->setPrecursorMz(260.186);
        c1->setDb("test_Mascot");
        double mz10[] = {60.437, 85.247,
                        99.242, 144.142,
                        201.043, 260.119};
        for(int i = 0; i < 6; i++)
            c1->fragmentMzValues.push_back(mz10[i]);
        double ion10[] = {17.549, 74.118, 5.221,
                         1.863, 37.647, 100.0};
        for(int i = 0; i < 6; i++)
            c1->fragmentIntensities.push_back(ion10[i]);
        REQUIRE(*c1 == *compounds[9]);
        delete (c1);

    }


    SUBCASE("Testing MSP File")
    {
        Database db;
        string mspFile = "tests/doctest/test_NISTLibrary.msp";
        int resMsp = db.loadNISTLibrary(mspFile);
        REQUIRE(resMsp == 10);
        vector<Compound*> compounds = db.getCompoundsSubset("test_NISTLibrary");

        Compound* c1 = new Compound("HMDB00902", "NAD-20.0,50.0,100.0",
                                   "C21H27N7O14P2", 0);
        c1->setMass(663.109121804);
        c1->category.push_back("None");
        c1->smileString = "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -3.6479;
        c1->setDb("test_NISTLibrary");
        double mz1[] = {62.964214325, 78.9593391418,
                        92.0254847209, 96.9695861816,
                        107.036354065, 134.047064209,
                        138.980453491, 146.98580424,
                        150.980529785, 158.925656128,
                        174.980461121, 176.936406453,
                        192.991004944, 211.001739502,
                        254.946289062, 272.956855774,
                        328.045959473, 346.056945801,
                        408.012435913, 426.022750854,
                        540.051066081, 540.052429199};
        for(int i = 0; i < 22; i++)
            c1->fragmentMzValues.push_back(mz1[i]);

        double ion1[] = {248, 8816,97,821, 289, 2205, 56, 107,
                         56, 3351, 539, 95, 179, 130, 54, 1175,
                         523, 488, 427, 614, 7137, 7118};
        for(int i = 0; i < 22; i++)
            c1->fragmentIntensities.push_back(ion1[i]);

        c1->fragmentIonTypes[0] = "2";
        c1->fragmentIonTypes[1] = "4";
        c1->fragmentIonTypes[2] = "3";
        c1->fragmentIonTypes[3] = "5";
        c1->fragmentIonTypes[4] = "3";
        c1->fragmentIonTypes[5] = "5";
        c1->fragmentIonTypes[6] = "2";
        c1->fragmentIonTypes[7] = "3";
        c1->fragmentIonTypes[8] = "2";
        c1->fragmentIonTypes[9] = "5";
        c1->fragmentIonTypes[10] = "4";
        c1->fragmentIonTypes[11] = "3";
        c1->fragmentIonTypes[12] = "4";
        c1->fragmentIonTypes[13] = "4";
        c1->fragmentIonTypes[14] = "2";
        c1->fragmentIonTypes[15] = "4";
        c1->fragmentIonTypes[16] = "4";
        c1->fragmentIonTypes[17] = "4";
        c1->fragmentIonTypes[18] = "4";
        c1->fragmentIonTypes[19] = "4";
        c1->fragmentIonTypes[20] = "3";
        c1->fragmentIonTypes[21] = "2";
        REQUIRE(*c1 == *(compounds[0]));
        delete(c1);

        c1 = new Compound("HMDB00641","L-GLUTAMINE-20.0,50.0,100.0",
                          "C5H10N2O3", 0);
        c1->setMass(146.06914218);
        c1->category.push_back("None");
        c1->smileString = "NC(=O)CCC(N)C(=O)O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -1.3362;
        c1->setDb("test_NISTLibrary");
        double mz2[] = {52.0191116333, 58.0298690796,
                        66.0349807739, 67.0301628113,
                        70.0298519135, 71.0140228271,
                        71.0251178741, 72.0091705322,
                        72.0455932617, 72.9931793213,
                        74.0247955322, 81.0460205078,
                        82.0252304077, 82.0298728943,
                        84.0404510498, 84.0455093384,
                        86.0247936249, 97.0407562256,
                        98.0247980754, 99.0564365387,
                        101.071949005, 102.056114197,
                        107.02511851, 109.040693283,
                        125.035812378, 125.390640259,
                        127.042938232, 127.051460266,
                        128.035438538, 138.631530762,
                        142.038706462, 145.062143962 };
        for(int i = 0; i < 32; i++)
            c1->fragmentMzValues.push_back(mz2[i]);
        double ion2[]={51, 2961,57,225, 140, 34, 134,
                       254, 28, 42, 1813, 33,50, 925,
                       61, 1891, 253, 110, 117, 185,
                       548, 17, 108, 768, 397, 26, 42,
                       2527, 1147, 23, 34, 3244 };
        for(int i = 0; i < 32; i++)
            c1->fragmentIntensities.push_back(ion2[i]);
        c1->fragmentIonTypes[0] = "2";
        c1->fragmentIonTypes[1] = "4";
        c1->fragmentIonTypes[2] = "2";
        c1->fragmentIonTypes[3] = "4";
        c1->fragmentIonTypes[4] = "4";
        c1->fragmentIonTypes[5] = "4";
        c1->fragmentIonTypes[6] = "4";
        c1->fragmentIonTypes[7] = "4";
        c1->fragmentIonTypes[8] = "3";
        c1->fragmentIonTypes[9] = "3";
        c1->fragmentIonTypes[10] = "4";
        c1->fragmentIonTypes[11] = "3";
        c1->fragmentIonTypes[12] = "2";
        c1->fragmentIonTypes[13] = "4";
        c1->fragmentIonTypes[14] = "2";
        c1->fragmentIonTypes[15] = "4";
        c1->fragmentIonTypes[16] = "4";
        c1->fragmentIonTypes[17] = "4";
        c1->fragmentIonTypes[18] = "3";
        c1->fragmentIonTypes[19] = "4";
        c1->fragmentIonTypes[20] = "4";
        c1->fragmentIonTypes[21] = "3";
        c1->fragmentIonTypes[22] = "3";
        c1->fragmentIonTypes[23] = "4";
        c1->fragmentIonTypes[24] = "3";
        c1->fragmentIonTypes[25] = "2";
        c1->fragmentIonTypes[26] = "2";
        c1->fragmentIonTypes[27] = "4";
        c1->fragmentIonTypes[28] = "4";
        c1->fragmentIonTypes[29] = "2";
        c1->fragmentIonTypes[30] = "3";
        c1->fragmentIonTypes[31] = "3";
        REQUIRE(*c1 == *compounds[1]);
        delete(c1);

        c1 = new Compound("HMDB00965","HYPOTAURINE-20.0,50.0,100.0",
                          "C2H7NO2S",0);
        c1->setMass(109.019749464);
        c1->category.push_back("None");
        c1->smileString = "NCCS(=O)O";
        c1->setDb("test_NISTLibrary");
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -0.8332;
        double mz3[] = {63.9624799093, 64.9702987671, 108.01247406};
        for(int i=0; i < 3; i++)
            c1->fragmentMzValues.push_back(mz3[i]);
        double ion3[] = {9645, 981, 671};
        for(int i = 0; i < 3; i++)
            c1->fragmentIntensities.push_back(ion3[i]);
        c1->fragmentIonTypes[0] = "3";
        c1->fragmentIonTypes[1] = "3";
        c1->fragmentIonTypes[2] = "2";
        REQUIRE(*c1 == *compounds[2]);
        delete(c1);

        c1 = new Compound("HMDB00175","INOSINE 5'-PHOSPHATE-20.0,50.0,100.0",
                          "C10H13N4O8P", 0);
        c1->setMass(348.047100006);
        c1->setDb("test_NISTLibrary");
        c1->category.push_back("None");
        c1->smileString ="O=c1nc[nH]c2c1ncn2C1OC(COP(=O)(O)O)C(O)C1O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -2.1519;
        double mz4[]= {62.964050293, 65.0145721436, 66.0098419189,
                       78.9590606689, 92.0253245036, 96.9695968628,
                       135.031138102, 150.980484009, 192.991226196,
                       211.001693726, 347.040802002};
        for(int i = 0; i < 11; i++)
            c1->fragmentMzValues.push_back(mz4[i]);
        double ion4[] = {61, 355, 67, 7527, 1374,
                         871, 230, 113, 69, 523,
                         6304};
        for(int i = 0; i < 11; i++)
            c1->fragmentIntensities.push_back(ion4[i]);
        c1->fragmentIonTypes[0] = "3";
        c1->fragmentIonTypes[1] = "2";
        c1->fragmentIonTypes[2] = "2";
        c1->fragmentIonTypes[3] = "3";
        c1->fragmentIonTypes[4] = "3";
        c1->fragmentIonTypes[5] = "3";
        c1->fragmentIonTypes[6] = "3";
        c1->fragmentIonTypes[7] = "2";
        c1->fragmentIonTypes[8] = "2";
        c1->fragmentIonTypes[9] = "2";
        c1->fragmentIonTypes[10] = "2";
        REQUIRE(*c1 == *compounds[3]);
        delete(c1);

        c1 = new Compound("HMDB00094","CITRATE-20.0,50.0,100.0",
                          "C6H8O7",0);
        c1->setDb("test_NISTLibrary");
        c1->setMass(192.027002596);
        c1->category.push_back("None");
        c1->smileString = "O=C(O)CC(O)(CC(=O)O)C(=O)O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -1.2485;
        double mz5[]= {57.0346002579, 58.9589700699,
                       59.0139074326, 67.0189723969,
                       72.9933039347, 73.029571533,
                       83.0139490763, 85.0295581818,
                       87.0087833405, 87.9253025055,
                       102.948770523, 103.040133158,
                       103.920129776, 105.93572998,
                       111.008729935, 123.94648234,
                       129.019475301, 130.998474121,
                       146.939005534, 147.030131022,
                       154.998906453, 173.009328206,
                       191.020187378};
        for(int i = 0; i < 23; i++)
            c1->fragmentMzValues.push_back(mz5[i]);
        double ion5[]={2598, 96, 337, 971, 27,
                       30, 33, 1730, 2943, 377,
                       540, 26, 48, 41, 4779, 29,
                       384, 40, 291, 33, 36, 198,
                       1634};
        for(int i = 0; i < 23; i++)
            c1->fragmentIntensities.push_back(ion5[i]);
        c1->fragmentIonTypes[0] = "4";
        c1->fragmentIonTypes[1] = "4";
        c1->fragmentIonTypes[2] = "4";
        c1->fragmentIonTypes[3] = "4";
        c1->fragmentIonTypes[4] = "3";
        c1->fragmentIonTypes[5] = "3";
        c1->fragmentIonTypes[6] = "3";
        c1->fragmentIonTypes[7] = "4";
        c1->fragmentIonTypes[8] = "4";
        c1->fragmentIonTypes[9] = "4";
        c1->fragmentIonTypes[10] = "4";
        c1->fragmentIonTypes[11] = "3";
        c1->fragmentIonTypes[12] = "4";
        c1->fragmentIonTypes[13] = "4";
        c1->fragmentIonTypes[14] = "4";
        c1->fragmentIonTypes[15] = "3";
        c1->fragmentIonTypes[16] = "3";
        c1->fragmentIonTypes[17] = "3";
        c1->fragmentIonTypes[18] = "3";
        c1->fragmentIonTypes[19] = "3";
        c1->fragmentIonTypes[20] = "3";
        c1->fragmentIonTypes[21] = "3";
        c1->fragmentIonTypes[22] = "3";
        REQUIRE(*c1 == *compounds[4]);
        delete(c1);

        c1 = new Compound("HMDB00094", "CITRATE-20.0,50.0,100.0 (1)",
                          "C6H8O7", 0);
        c1->setMass(192.027002596);
        c1->setDb("test_NISTLibrary");
        c1->category.push_back("None");
        c1->smileString = "O=C(O)CC(O)(CC(=O)O)C(=O)O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -1.2485;
        double mz6[] = {57.0345865885, 58.9589468638,
                        59.0138645172, 67.0189628601,
                        83.0137527466, 85.0295346578,
                        87.0087547302, 87.9252548218,
                        102.948786418, 102.997827148,
                        103.919984182, 105.935714722,
                        106.943557739, 111.008686066,
                        123.94655482,129.019424438,
                        130.998524984, 146.938732147,
                        147.029987335, 154.998645782,
                        173.008839925, 191.020000458};
        for(int i =0; i < 22; i++)
            c1->fragmentMzValues.push_back(mz6[i]);
        double ion6[] = {2419, 265, 295, 509, 26,
                         1641, 2912, 1013, 703, 45,
                         218, 175, 48, 4798, 615,362,
                         31,296, 27,40,140, 1087};
        for(int i = 0; i < 22; i++)
            c1->fragmentIntensities.push_back(ion6[i]);
        c1->fragmentIonTypes[0] = "6";
        c1->fragmentIonTypes[1] = "6";
        c1->fragmentIonTypes[2] = "6";
        c1->fragmentIonTypes[3] = "6";
        c1->fragmentIonTypes[4] = "5";
        c1->fragmentIonTypes[5] = "6";
        c1->fragmentIonTypes[6] = "6";
        c1->fragmentIonTypes[7] = "6";
        c1->fragmentIonTypes[8] = "6";
        c1->fragmentIonTypes[9] = "5";
        c1->fragmentIonTypes[10] = "6";
        c1->fragmentIonTypes[11] = "6";
        c1->fragmentIonTypes[12] = "5";
        c1->fragmentIonTypes[13] = "6";
        c1->fragmentIonTypes[14] = "6";
        c1->fragmentIonTypes[15] = "4";
        c1->fragmentIonTypes[16] = "3";
        c1->fragmentIonTypes[17] = "4";
        c1->fragmentIonTypes[18] = "4";
        c1->fragmentIonTypes[19] = "4";
        c1->fragmentIonTypes[20] = "3";
        c1->fragmentIonTypes[21] = "4";
        REQUIRE(*c1 == *compounds[5]);
        delete(c1);

        c1 = new Compound("HMDB00094","CITRATE-20.0,50.0,100.0 (2)",
                          "C6H8O7",0);
        c1->setDb("test_NISTLibrary");
        c1->setMass(192.027002596);
        c1->category.push_back("None");
        c1->smileString ="O=C(O)CC(O)(CC(=O)O)C(=O)O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -1.2485;
        double mz7[] = {57.0345850405, 59.0138343464,
                        59.9852591621, 67.0189613674,
                        67.8221041361, 68.8053665161,
                        69.5959014893, 83.0137751653,
                        85.0295320594, 87.0087479301,
                        92.4431991577, 101.024281979,
                        102.948895264, 103.040128371,
                        108.914337158, 111.00138855,
                        111.008674622, 126.157676697,
                        129.019415463, 130.998455048,
                        147.030051736, 154.998764038,
                        173.009250641, 190.186141968,
                        191.020112879, 209.673370361};
        for(int i = 0; i < 26; i++)
            c1->fragmentMzValues.push_back(mz7[i]);
        double ion7[] = {3297, 372, 14, 1177, 12,
                         37, 39, 47, 1748, 3150,
                         41, 20, 9,  29,  38, 77,
                         5432, 37, 450, 46, 42,
                         46, 173, 36,1360, 36};
        for(int i = 0; i < 26; i++)
            c1->fragmentIntensities.push_back(ion7[i]);
        c1->fragmentIonTypes[0] = "23";
        c1->fragmentIonTypes[1] = "22";
        c1->fragmentIonTypes[2] = "9";
        c1->fragmentIonTypes[3] = "23";
        c1->fragmentIonTypes[4] = "6";
        c1->fragmentIonTypes[5] = "2";
        c1->fragmentIonTypes[6] = "2";
        c1->fragmentIonTypes[7] = "13";
        c1->fragmentIonTypes[8] = "23";
        c1->fragmentIonTypes[9] = "23";
        c1->fragmentIonTypes[10] = "2";
        c1->fragmentIonTypes[11] = "16";
        c1->fragmentIonTypes[12] = "15";
        c1->fragmentIonTypes[13] = "17";
        c1->fragmentIonTypes[14] = "2";
        c1->fragmentIonTypes[15] = "2";
        c1->fragmentIonTypes[16] = "23";
        c1->fragmentIonTypes[17] = "2";
        c1->fragmentIonTypes[18] = "17";
        c1->fragmentIonTypes[19] = "16";
        c1->fragmentIonTypes[20] = "17";
        c1->fragmentIonTypes[21] = "17";
        c1->fragmentIonTypes[22] = "16";
        c1->fragmentIonTypes[23] = "2";
        c1->fragmentIonTypes[24] = "17";
        c1->fragmentIonTypes[25] = "2";
        REQUIRE(*c1 == *compounds[6]);
        delete(c1);

        c1 = new Compound("HMDB00167","L-THREONINE-20.0,50.0,100.0",
                          "C4H9NO3",0);
        c1->setDb("test_NISTLibrary");
        c1->setMass(119.058243148);
        c1->category.push_back("None");
        c1->smileString ="CC(O)C(N)C(=O)O";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -1.2209;
        double mz8[]={72.0091552734, 74.0247751872, 118.050979614};
        for(int i = 0; i < 3; i++)
            c1->fragmentMzValues.push_back(mz8[i]);
        double ion8[] = {3389, 8436, 1811};
        for(int i = 0; i < 3; i++)
            c1->fragmentIntensities.push_back(ion8[i]);
        c1->fragmentIonTypes[0] = "3";
        c1->fragmentIonTypes[1] = "3";
        c1->fragmentIonTypes[2] = "2";
        REQUIRE(*c1 == *compounds[7]);
        delete(c1);

        c1 = new Compound("HMDB01366","PURINE-20.0,50.0,100.0",
                          "C5H4N4",0);

        c1->setDb("test_NISTLibrary");
        c1->setMass(120.043596128);
        c1->setCollisionEnergy(200);
        c1->ionizationMode = -1;
        c1->category.push_back("None");
        c1->smileString ="c1ncc2[nH]cnc2n1";
        c1->logP = 0.3529;
        double mz9[] = {64.0066299438, 65.0145670573,
                        67.0302454631,68.0254770915,
                        90.0099182129, 92.0254681905,
                        119.036338806};
        for(int i = 0; i < 7; i++)
            c1->fragmentMzValues.push_back(mz9[i]);
        double ion9[] = {62, 1786, 515, 204, 257,
                         1913,9786};
        for(int i = 0; i < 7; i++)
            c1->fragmentIntensities.push_back(ion9[i]);
        c1->fragmentIonTypes[0] = "6";
        c1->fragmentIonTypes[1] = "6";
        c1->fragmentIonTypes[2] = "6";
        c1->fragmentIonTypes[3] = "6";
        c1->fragmentIonTypes[4] = "6";
        c1->fragmentIonTypes[5] = "6";
        c1->fragmentIonTypes[6] = "6";
        REQUIRE(*c1 == *compounds[8]);
        delete(c1);

        c1 = new Compound("HMDB00230","N-ACETYLNEURAMINATE-20.0,50.0,100.0",
                          "C11H19NO9", 0);
        c1->setDb("test_NISTLibrary");
        c1->setMass(309.105981188);
        c1->category.push_back("None");
        c1->smileString = "CC(=O)NC1C(O)CC(O)(C(=O)O)OC1C(O)C(O)CO";
        c1->ionizationMode = -1;
        c1->setCollisionEnergy(200);
        c1->logP = -3.8718;
        double mz10[] = {52.0190925598, 55.0188713074,
                         57.0343933105, 58.0060195923,
                         58.0298500061, 59.0138956706,
                         59.9851913452, 66.0349349976,
                         71.0138778687, 71.0465774536,
                         71.0502929688, 72.0091629028,
                         72.9930648804, 82.0300013224,
                         84.0457458496, 87.0087865194,
                         89.0243759155, 98.0549468994,
                         98.0611139933, 100.040331523,
                         101.024398804, 119.035064697,
                         142.051315308, 170.046279907,
                         220.083190918, 290.088500977,
                         308.099243164};
        for(int i = 0; i < 27; i++)
            c1->fragmentMzValues.push_back(mz10[i]);
        double ion10[] = {140 , 46, 37, 37, 421, 1453,
                          34, 131, 227, 46, 1140, 49,
                          119, 277, 41, 6222, 60, 43,
                          1098, 127, 209, 409, 37, 1567,
                          71, 136, 907};
        for(int i = 0; i < 27; i++)
            c1->fragmentIntensities.push_back(ion10[i]);
        c1->fragmentIonTypes[0] = "2";
        c1->fragmentIonTypes[1] = "2";
        c1->fragmentIonTypes[2] = "2";
        c1->fragmentIonTypes[3] = "3";
        c1->fragmentIonTypes[4] = "3";
        c1->fragmentIonTypes[5] = "3";
        c1->fragmentIonTypes[6] = "2";
        c1->fragmentIonTypes[7] = "2";
        c1->fragmentIonTypes[8] = "3";
        c1->fragmentIonTypes[9] = "2";
        c1->fragmentIonTypes[10] = "2";
        c1->fragmentIonTypes[11] = "2";
        c1->fragmentIonTypes[12] = "2";
        c1->fragmentIonTypes[13] = "3";
        c1->fragmentIonTypes[14] = "2";
        c1->fragmentIonTypes[15] = "3";
        c1->fragmentIonTypes[16] = "2";
        c1->fragmentIonTypes[17] = "2";
        c1->fragmentIonTypes[18] = "3";
        c1->fragmentIonTypes[19] = "3";
        c1->fragmentIonTypes[20] = "2";
        c1->fragmentIonTypes[21] = "2";
        c1->fragmentIonTypes[22] = "2";
        c1->fragmentIonTypes[23] = "2";
        c1->fragmentIonTypes[24] = "2";
        c1->fragmentIonTypes[25] = "2";
        c1->fragmentIonTypes[26] = "2";
        REQUIRE(*c1 == *compounds[9]);
        delete(c1);
    }

}

