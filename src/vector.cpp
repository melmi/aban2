/**
 * vector.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "vector.h"
#include <cmath>

namespace aban2
{

vector::vector(): vector(0, 0, 0)
{
}

vector::vector(double _x, double _y, double _z): x(_x), y(_y), z(_z)
{
}

vector vector::from_data(double **data, size_t i)
{
    return vector(data[0][i], data[1][i], data[2][i]);
}

void vector::to_data(double **data, size_t i) const
{
    data[0][i] = x;
    data[1][i] = y;
    data[2][i] = z;
}

vector vector::operator*(double r) const
{
    return {r * x, r * y, r * z};
}

double vector::operator*(vector v) const
{
    return x * v.x + y * v.y + z * v.z;
}

vector vector::operator^(vector v) const
{
    return {y *v.z - z * v.y, z *v.x - x * v.z, x *v.y - y * v.x};
}

vector vector::operator+(vector v) const
{
    return {x + v.x, y + v.y, z + v.z};
}

vector vector::operator-(vector v) const
{
    return {x - v.x, y - v.y, z - v.z};
}

vector vector::operator-() const
{
    return -1.0 * (*this);
}

vector &vector::operator*=(double r)
{
    x *= r;
    y *= r;
    z *= r;
    return *this;
}

vector &vector::operator+=(vector v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

vector &vector::operator-=(vector v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

double vector::l2() const
{
    return x * x + y * y + z * z;
}

double vector::l() const
{
    return std::sqrt(l2());
}

void vector::normalize(double epsilon)
{
    double len = l() + epsilon;
    x /= len;
    y /= len;
    z /= len;
}

void vector::normalize()
{
    normalize(1e-8);
}

vector operator*(double r, const vector v)
{
    return v * r;
}

std::ostream& operator<< (std::ostream &out, vector v)
{
    out << "(" 
        << v.x << ", "
        << v.y << ", "
        << v.z << ")";
    return out;
}

}
