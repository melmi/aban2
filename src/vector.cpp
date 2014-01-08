/**
 * vector.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "vector.h"

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

double vector::operator*(vector v)
{
    return x * v.x + y * v.y + z * v.z;
}

vector vector::operator+(vector v)
{
    return {x + v.x, y + v.y, z + v.z};
}

vector vector::operator-(vector v)
{
    return {x + v.x, y + v.y, z + v.z};
}

}
