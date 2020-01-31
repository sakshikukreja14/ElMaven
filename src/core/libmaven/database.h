#ifndef DATABASE_H
#define DATABASE_H

#include <boost/signals2.hpp>
#include "standardincludes.h"

using namespace std;
class Adduct;
class Compound;
class MassCutoff;
class Pathway;
class Peak;
class Reaction;

namespace bsignal = boost::signals2;

/*
class Molecule2D {
       public:
	Molecule2D() {}
	QString id;
	QVector<QPointF> coord;
	QVector<QString> atoms;
};

*/

class Database {
       public:
            Database() { };
            Database(string filename) {	}
            ~Database() { }

            /**
             * @brief Remove an already loaded database (and all compounds part of it)
             * from the DB store.
             * @param dbName Name of the database to be removed.
             */
            void removeDatabase(string dbName);

            void loadKnowns();
            int loadCompoundCSVFile(string filename);
            bool addCompound(Compound* c);
            vector<Compound*> getCompoundsSubset(string database);
            vector<Compound*> getKnowns();

            /**
             * @brief Load metabolites from a file at a given path by treating it as
             * having NIST library format.
             * @param filepath The absolute path of the NIST library file.
             * @param signal Pointer to a boost signal object that can be called
             * with a string for update message, an integer for current steps of
             * progress and another integer for total steps to completion.
             * @return The number of compounds that were loaded into the database.
             */
            int loadNISTLibrary(string filepath,
                            bsignal::signal<void (string, int, int)>* signal=nullptr);

            int loadMascotLibrary(string filepath,
                              bsignal::signal<void (string, int, int)>* signal=nullptr);

            /**
             * @brief Checks whether the library with the given name is an NIST
             * library or not.
             * @details The first compound in the database (if it has any compounds)
             * is checked for PRM information and if found, the database is
             * regarded as an NIST library.
             * @param dbName String name of the database to be checked.
             * @return True if the database with given name is an NIST library,
             * false otherwise.
             */
            bool isSpectralLibrary(string dbName);

            map<string, int> getDatabaseNames();
            Compound* findSpeciesByIdAndName(string id, string name, string dbName);
            deque<Compound*> getCompoundsDB(){ 	return compoundsDB;}	set<Compound*> findSpeciesByMass(float mz, MassCutoff *massCutoff);
            vector<Compound*> findSpeciesByName(string name, string dbname);
            vector<Compound*> findSpeciesById(string id, string dbName);
            Compound* extractCompoundfromEachLine(vector<string>& fields, map<string, int> & header, int loadCount, string filename);
            float getChargeFromDB(vector<string>& fields, map<string, int> & header);
            vector<string> getCategoryFromDB(vector<string>& fields, map<string, int> & header);

            deque<Compound*> compoundsDB;
            map<string, Compound*> compoundIdNameDbMap;
            vector<string> invalidRows;
            vector<string> notFoundColumns;
            map<string, int> compoundIdenticalCount;
            //Added while merging with Maven776 - Kiran
            const std::string ANYDATABASE;

    private:
            bool _startsWith(string line, string text);
            bool _contain(string line, string text);

};

extern Database DB;

#endif
