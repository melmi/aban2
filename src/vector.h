/**
 * vector.h
 *
 * by Mohammad Elmi
 * Created on 28 Aug 2013
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "common.h"
#include <ostream>

namespace aban2
{

struct vector
{
    union
    {
        double cmpnt[3];
        struct
        {
            double x, y, z;
        };
    };

    vector();
    vector(double _x, double _y, double _z);
    explicit vector(double c): vector(c, c, c) {}

    vector operator*(double r) const;
    double operator*(vector v) const;
    vector operator^(vector v) const;
    vector operator+(vector v) const;
    vector operator-(vector v) const;
    vector operator-() const;
    vector &operator*=(const double r);
    vector &operator+=(const vector v);
    vector &operator-=(const vector v);

    double l2() const;
    double l() const;
    void normalize();
    void normalize(double epsilon);

    void to_data(double **data, size_t i) const;
    static vector from_data(double **data, size_t i);
};

vector operator*(const double r, const vector v);
std::ostream& operator<< (std::ostream &out, vector v);

}

#endif