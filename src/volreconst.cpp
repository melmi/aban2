#include "volreconst.h"

#include <cmath>
#include <tuple>
#include <algorithm>

namespace aban2
{

/******************** volreconst ********************/

volreconst *volreconst::from_base_data(vector _c, vector _m)
{
    bool has_elem[3];
    vector m = {std::abs(_m.x), std::abs(_m.y), std::abs(_m.z)};
    if (m.l2() < epsilon)m.x = 1;

    double base_vol = 1;
    size_t ndirs = 0;
    if (has_elem[0] = m.x > epsilon) base_vol /= m.x * ++ndirs; else base_vol *= _c.x;
    if (has_elem[1] = m.y > epsilon) base_vol /= m.y * ++ndirs; else base_vol *= _c.y;
    if (has_elem[2] = m.z > epsilon) base_vol /= m.z * ++ndirs; else base_vol *= _c.z;

    volreconst *result;
    switch (ndirs)
    {
    case 1: result = new volreconst1d(); break;
    case 2: result = new volreconst2d(); break;
    case 3: result = new volreconst3d(); break;
    }

    result->c = _c;
    result->orig_m = _m;
    result->m = m;
    result->alpha_max = m * _c;
    result->base_vol = base_vol;
    std::copy_n(has_elem, 3, result->has_elem);

    return result;
}

void volreconst::set_alpha()
{
    double v0 = 0, v1 = c.x * c.y * c.z;
    double a0 = 0, a1 = alpha_max;
    while (std::abs(a0 - a1) > epsilon)
    {
        double anew = (a0 + a1) / 2.0;
        double vnew = this->get_volume(anew);
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

void volreconst::set_volume()
{
    volume = this->get_volume(alpha);
}

double volreconst::get_correct_moment(size_t dir, double moment)
{
    // terms inside double prantheses are distances
    if (!std::signbit(moment)) return moment;

    return c.x * c.y * c.z * ((c.cmpnt[dir] / 2.0)) -
           volume * ((-c.cmpnt[dir] + moment / volume));
}

volreconst *volreconst::get_remaining(size_t dir, double delta)
{
    bool from_start = delta * orig_m.cmpnt[dir] < 0.0;
    delta = std::abs(delta);

    vector new_c = c;
    new_c.cmpnt[dir] -= delta;

    double new_alpha;
    if (from_start)
        new_alpha = -(delta * m.cmpnt[dir] - alpha);
    else
        new_alpha = alpha;

    return volreconst::from_alpha(new_c, m, new_alpha);
}

std::tuple<double, vector> volreconst::get_flux(size_t dir, double delta)
{
    auto remaining = get_remaining(dir, delta);
    double v0 = volume, v1 = remaining->volume;
    vector q0 { this->get_moment(0), this->get_moment(1), this->get_moment(2)};
    vector q1 { remaining->get_moment(0), remaining->get_moment(1), remaining->get_moment(2)};
    delete remaining;

    // moving q1 to the axis of q0
    if (delta > 0)
    {
        vector d;
        d.cmpnt[dir] = delta;
        q1 += v1 * d;
    }

    double vf = v0 - v1;
    vector qf = q0 - q1;

    // moving qf to the axis at cell center
    qf += (vf * 0.5) * c;

    return std::make_tuple(vf, qf);
}

volreconst *volreconst::from_volume(vector c, vector m, double volume)
{
    volreconst *result = volreconst::from_base_data(c, m);
    result->volume = volume;
    result->set_alpha();
    return result;
}

volreconst *volreconst::from_alpha(vector c, vector m, double alpha)
{
    volreconst *result = volreconst::from_base_data(c, m);
    result->alpha = alpha;
    result->set_volume();
    return result;
}

/******************** volreconst1d ********************/

double volreconst1d::get_volume(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    return _alpha * base_vol;
}

double volreconst1d::get_moment(size_t dir)
{
    // terms inside double prantheses are distances
    if (!has_elem[dir])return volume * ((c.cmpnt[dir] / 2.0));

    double result = volume * ((alpha / 2.0));
    return get_correct_moment(dir, result);
}

/******************** volreconst2d ********************/

double volreconst2d::x2(double x)
{
    if (x < 0.0)return 0.0;
    return x * x;
}

double volreconst2d::get_volume(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    double result = x2(_alpha);

    if (has_elem[0]) result -= x2(_alpha - m.x * c.x);
    if (has_elem[1]) result -= x2(_alpha - m.y * c.y);
    if (has_elem[2]) result -= x2(_alpha - m.z * c.z);

    return result * base_vol;
}

double volreconst2d::get_moment(size_t dir)
{
    // terms inside double prantheses are distances

    if (!has_elem[dir]) return volume * ((c.cmpnt[dir] / 2.0));

    double result = x2(alpha) * ((alpha / m.cmpnt[dir] / 3.0));

    for (int i = 0; i < 3; ++i)
        if (has_elem[i])
        {
            double d1 = i == dir ? c.cmpnt[i] : 0;
            double r1 = alpha - m.cmpnt[i] * c.cmpnt[i];
            result -= x2(r1) * ((d1 + r1 / m.cmpnt[i]  / 3.0));
        }
    result *= base_vol;
    return get_correct_moment(dir, result);
}

/******************** volreconst3d ********************/

double volreconst3d::x3(double x)
{
    if (x < 0)return 0;
    return x * x * x;
}

double volreconst3d::get_volume(double _alpha)
{
    if (_alpha < 0)return 0;
    if (_alpha > alpha_max) return c.x * c.y * c.z;

    double result = x3(_alpha);

    result -= x3(_alpha - m.x * c.x);
    result -= x3(_alpha - m.y * c.y);
    result -= x3(_alpha - m.z * c.z);

    result += x3(_alpha - alpha_max + m.x * c.x);
    result += x3(_alpha - alpha_max + m.y * c.y);
    result += x3(_alpha - alpha_max + m.z * c.z);

    return result * base_vol;
}

double volreconst3d::get_moment(size_t dir)
{
    // terms inside double prantheses are distances

    double result = x3(alpha) * ((alpha / m.cmpnt[dir] / 3.0));

    for (int i = 0; i < 3; ++i)
    {
        double d1 = i == dir ? c.cmpnt[i] : 0;
        double r1 = alpha - m.cmpnt[i] * c.cmpnt[i];
        result -= x3(r1) * ((d1 + r1 / m.cmpnt[i] / 3.0));

        double d2 = i == dir ? 0 : c.cmpnt[i];
        double r2 = alpha - alpha_max + m.cmpnt[i] * c.cmpnt[i];
        result += x3(r2) * ((d2 + r2 / m.cmpnt[i] / 3.0));
    }

    result *= base_vol;
    return get_correct_moment(dir, result);
}

}