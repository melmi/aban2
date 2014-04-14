#include "ireconst.h"

#include <cmath>
#include <iostream>

namespace aban2
{

inline double x2(double x)
{
    return x * x;
}

inline double x3(double x)
{
    return x * x * x;
}

ireconst::vol_func_t ireconst::vol_funcs[]
{
    &ireconst::get_volume_1d,
    &ireconst::get_volume_2d,
    &ireconst::get_volume_3d
};

double ireconst::get_volume_3d(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    double result = x3(_alpha);

    if (_alpha - alpha_max + m.x * c.x > 0)result += x3(_alpha - alpha_max + m.x * c.x);
    if (_alpha - alpha_max + m.y * c.y > 0)result += x3(_alpha - alpha_max + m.y * c.y);
    if (_alpha - alpha_max + m.z * c.z > 0)result += x3(_alpha - alpha_max + m.z * c.z);

    if (_alpha - m.x * c.x > 0)result -= x3(_alpha - m.x * c.x);
    if (_alpha - m.y * c.y > 0)result -= x3(_alpha - m.y * c.y);
    if (_alpha - m.z * c.z > 0)result -= x3(_alpha - m.z * c.z);

    return result * base_vol;
}

double ireconst::get_volume_2d(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    double result = x2(_alpha);

    if (m.x > epsilon) if (_alpha - m.x * c.x > 0)result -= x2(_alpha - m.x * c.x);
    if (m.y > epsilon) if (_alpha - m.y * c.y > 0)result -= x2(_alpha - m.y * c.y);
    if (m.z > epsilon) if (_alpha - m.z * c.z > 0)result -= x2(_alpha - m.z * c.z);

    return result * base_vol;
}

double ireconst::get_volume_1d(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    return _alpha * base_vol;
}

void ireconst::init(vector _c, vector _m)
{
    c = _c;
    orig_m = _m;
    m = {std::abs(_m.x), std::abs(_m.y), std::abs(_m.z)};
    if (m.l2() < epsilon)m.x = 1;
    alpha_max = m * c;

    base_vol = 1;
    int ndirs = 0;
    if (m.x > epsilon) base_vol /= m.x * ++ndirs; else base_vol *= c.x;
    if (m.y > epsilon) base_vol /= m.y * ++ndirs; else base_vol *= c.y;
    if (m.z > epsilon) base_vol /= m.z * ++ndirs; else base_vol *= c.z;

    vol_func = vol_funcs[ndirs - 1];
}

void ireconst::set_alpha()
{
    double v0 = 0, v1 = c.x * c.y * c.z;
    double a0 = 0, a1 = alpha_max;
    while (std::abs(a0 - a1) > epsilon)
    {
        double anew = (a0 + a1) / 2.0;
        double vnew = (this->*vol_func)(anew);
        if (volume > vnew)
        {
            a0 = anew;
            v0 = vnew;
        }
        else
        {
            a1 = anew;
            v1 = vnew;
        }
    }
    alpha = a0;
}

void ireconst::set_volume()
{
    volume = (this->*vol_func)(alpha);
}

ireconst ireconst::get_remaining(size_t dir, double delta)
{
    bool from_start = delta * orig_m.components[dir] < 0.0;
    delta = std::abs(delta);

    vector new_c = c;
    new_c.components[dir] -= delta;

    double new_alpha;
    if (from_start)
        new_alpha = -(delta * m.components[dir] - alpha);
    else
        new_alpha = alpha;

    return ireconst::from_alpha(new_c, m, new_alpha);
}

void ireconst::split(size_t dir, double delta, vector orig_m, ireconst &remaining, ireconst &departing)
{
    remaining = get_remaining(dir, delta, orig_m);
    delta = std::copysign(c.components[dir] - std::abs(delta), -delta);
    departing = get_remaining(dir, delta, orig_m);
}

ireconst ireconst::from_volume(vector c, vector m, double volume)
{
    ireconst result;
    result.init(c, m);
    result.volume = volume;
    result.set_alpha();
    return result;
}

ireconst ireconst::from_alpha(vector c, vector m, double alpha)
{
    ireconst result;
    result.init(c, m);
    result.alpha = alpha;
    result.set_volume();
    return result;
}

ireconst ireconst::from_full  (vector c, vector m)
{
    ireconst result;
    result.init(c, m);
    result.alpha = result.alpha_max;
    result.volume = c.x * c.y * c.z;
    return result;
}

ireconst ireconst::from_empty (vector c, vector m)
{
    ireconst result;
    result.init(c, m);
    result.alpha = 0;
    result.volume = 0;
    return result;
}
}