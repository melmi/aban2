/**
 * vector.h
 *
 * by Mohammad Elmi
 * Created on 28 Aug 2013
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <iostream>

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

    vector();

    vector(double _x, double _y, double _z);

    double operator*(vector v);
    vector operator+(vector v);
    vector operator-(vector v);

    double l2();
    double l();
    void normalize();

    static vector from_data(double **data, size_t i);
};

}

#endif