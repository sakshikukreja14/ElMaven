#ifndef MZPOINT_H
#define MZPOINT_H

/**
 * @brief Store point in rt, intensity and mz space (3 Dimensional)
 *
 * @details The data we get from Mass Spectrometry has 3 datatypes, namely
 * rt (retention time), intensity and mz (Mass by charge ratio). Class mzPoint
 * stores a point in 3-D space where x, y and z are rt, intensity and mz axis
 *
 */
class mzPoint
{
    public:
    /**
         * @brief Constructor with default values
         */
    mzPoint()
    {
        _x = _y = _z = 0;
    }

    /**
         * @brief Constructor with input values
         * @param ix Value in x axis
         * @param iy Value in y axis
         * @param iz Value in z axis
         */
    mzPoint(double ix, double iy, double iz)
    {
        _x = ix;
        _y = iy;
        _z = iz;
    }

    /**
         * @brief operator =    Assignment operator overloaded
         * @param b
         * @return
         */
    mzPoint& operator=(const mzPoint& b)
    {
        _x = b._x;
        _y = b._y;
        _z = b._z;
        return *this;
    }

    /**
         * @brief Compare point in x plane
         * @param a object of class mzPoint
         * @param b object of class mzPoint
         * @return bool True if b is greater than a in x plane
         */
    static bool compX(const mzPoint& a, const mzPoint& b)
    {
        return a._x < b._x;
    }

    /**
         * @brief Compare point in y plane
         * @param a object of class mzPoint
         * @param b object of class mzPoint
         * @return bool True if b is greater than a in y plane
         */
    static bool compY(const mzPoint& a, const mzPoint& b)
    {
        return a._y < b._y;
    }

    /**
         * @brief Compare point in z plane
         * @param a object of class mzPoint
         * @param b object of class mzPoint
         * @return bool True if b is greater than a in z plane
         */
    static bool compZ(const mzPoint& a, const mzPoint& b)
    {
        return a._z < b._z;
    }

    /**
         * @brief x returns _x
         * @return
         */
    double x()
    {
        return _x;
    }

    /**
         * @brief y returns _y
         * @return
         */
    double y()
    {
        return _y;
    }

    /**
         * @brief z return _z
         * @return
         */
    double z()
    {
        return _z;
    }

    private:
    double _x, _y, _z;
};


#endif // MZPOINT_H
