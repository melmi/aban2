/**
 * vector.h
 *
 * by Mohammad Elmi
 * Created on 28 Aug 2013
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

namespace aban2
{
struct vector
{
    union
    {
        double components[3];
        struct
        {
            double x, y, z;
        };
    };

    vector(): vector(0, 0, 0)
    {
    }

    vector(double _x, double _y, double _z): x(_x), y(_y), z(_z)
    {
    }
    
    static vector from_data(double **data, size_t i)
    {
        return vector(data[0][i], data[1][i], data[2][i]);
    }
};
}

#endif