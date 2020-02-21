#ifndef MZLINK_H
#define MZLINK_H

#include "standardincludes.h"

using namespace std;
/**
 * @brief Link two mass by charge ratios
 */
class mzLink
{
    public:
    /**
         * @brief Constructor for mzLink with no input arguments
         */
    mzLink();

    /**
         * @brief Constructor for mzLink with integer mzs
         */
    mzLink(int a, int b, string n);

    /**
         * @brief Constructor for mzLink with float mzs
         */
    mzLink(float a, float b, string n);

    /**
         * @brief mzLink Constructor for mzLink with
         * float mzs and values.
         */
    mzLink(float mz1, float mz2,
           float value1, float value2,
           string n);

    /**
         * @brief Destructor for mzLink
         */
    ~mzLink() {}

    string note;

    /**
         * @brief Compare m/z of two mzLinks
         * @param a mzLink 1
         * @param b mzLink 2
         * @return True if mzLink 1 has lower m/z than mzLink b, else false
         */
    static bool compMz(const mzLink &a, const mzLink &b)
    {
        return a._mz1 < b._mz1;
    }

    /**
         * @brief Compare correlation of two links
         * @param a mzLink 1
         * @param b mzLink 2
         * @return True if Link 1 has lower correlation than Link b, else false
         */
    static bool compCorrelation(const mzLink &a, const mzLink &b)
    {
        return a._correlation > b._correlation;
    }

    /**
         * @brief setCorrelation set the value for correlation
         * @param a
         */
    void setCorrelation(float a)
    {
        this->_correlation = a;
    }

    /**
         * @brief correlation returns correlation
         * @return
         */
    float correlation()
    {
        return this->_correlation;
    }

    /**
         * @brief mz1 returns mz1
         * @return
         */
    float mz1()
    {
        return this->_mz1;
    }

    /**
         * @brief mz2 returns mz2
         * @return
         */
    float mz2()
    {
        return this->_mz2;
    }

    /**
         * @brief value1 retruns value1
         * @return
         */
    float value1(){
        return this->_value1;
    }

    /**
         * @brief value2 returns value2
         * @return
         */
    float value2(){
        return this->_value2;
    }

    /**
         * @brief setValue1 sets value1
         * @param a
         */
    void setValue1(float a)
    {
        this->_value1 = a;
    }

    /**
         * @brief setValue2 sets value2
         * @param a
         */
    void setValue2(float a)
    {
        this->_value2 = a;
    }

    private:
    float _mz1;
    float _mz2;
    float _value1;
    float _value2;
    float _correlation;

};

#endif // MZLINK_H
