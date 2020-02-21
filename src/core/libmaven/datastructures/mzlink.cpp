#ifndef MZLINK_CPP
#define MZLINK_CPP

#include "mzlink.h"

mzLink::mzLink()
{
    _mz1 = _mz2 = 0.0;
    _value1 = _value2 = 0.0;
    //  data1 = data2 = NULL;
    _correlation = 0;
}

mzLink::mzLink(int a, int b, string n)
{
    _mz1 = a;
    _mz2 = b;
    note = n;
    _value1 = 0.0;
    _value2 = 0.0;
    _correlation = 0;
}
mzLink::mzLink(float a, float b, string n)
{
    _mz1 = a;
    _mz2 = b;
    note = n;
    _value1 = 0.0;
    _value2 = 0.0;
    _correlation = 0;
}

mzLink::mzLink(float mz1, float mz2,  float value1,
               float value2, string note)
{
    _mz1 = mz1;
    _mz2 = mz2;
    _value1 = value1;
    _value2 = value2;
    this->note  = note;
}
#endif // MZLINK_CPP
