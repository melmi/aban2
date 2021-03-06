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

    const double min_ratio = 0.001;
    has_elem[0] = std::abs(m.x / m.y) > min_ratio && std::abs(m.x / m.z) > min_ratio;
    has_elem[1] = std::abs(m.y / m.z) > min_ratio && std::abs(m.y / m.x) > min_ratio;
    has_elem[2] = std::abs(m.z / m.x) > min_ratio && std::abs(m.z / m.y) > min_ratio;

    double base_vol = 1;
    size_t ndirs = 0;
    for (int i = 0; i < 3; ++i)
        if (has_elem[i])
            base_vol /= m.cmpnt[i] * ++ndirs;
        else
        {
            base_vol *= _c.cmpnt[i];
            m.cmpnt[i] = _m.cmpnt[i] = 0;
        }
    m.normalize();
    _m.normalize();

    volreconst *result;
    switch (ndirs)
    {
    case 1:
        result = new volreconst1d();
        break;
    case 2:
        result = new volreconst2d();
        break;
    case 3:
        result = new volreconst3d();
        break;
    default:
        return nullptr;
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
    double a0 = 0, a1 = alpha_max;
    while (std::abs(a0 - a1) > epsilon)
    {
        double anew = (a0 + a1) / 2.0;
        double vnew = this->get_volume(anew);
        if (volume > vnew)
            a0 = anew;
        else
            a1 = anew;
    }
    alpha = a0;
}

void volreconst::set_volume()
{
    volume = this->get_volume(alpha);
}

volreconst *volreconst::get_cut(size_t dir, double delta)
{
    bool from_start = (delta * orig_m.cmpnt[dir]) < 0.0;
    delta = std::abs(delta);
    vector new_c = c;
    new_c.cmpnt[dir] = delta;

    double new_alpha = alpha;
    if (!from_start)
        new_alpha -= (c.cmpnt[dir] - delta) * m.cmpnt[dir];

    return volreconst::from_alpha(new_c, orig_m, new_alpha);
}

std::tuple<double, vector> volreconst::get_flux(size_t dir, double delta)
{
    auto cut = get_cut(dir, delta);
    double v = cut->volume;
    vector q = cut->get_moments();
    delete cut;

    // moving q to the axis of current volume
    if (delta > 0)q.cmpnt[dir] += v * (c.cmpnt[dir] - delta);
    // moving q to the axis at cell center
    q -= (v * 0.5) * c;
    // projecting q if necessary
    if (has_elem[0] && orig_m.cmpnt[0] < 0)q.cmpnt[0] *= -1;
    if (has_elem[1] && orig_m.cmpnt[1] < 0)q.cmpnt[1] *= -1;
    if (has_elem[2] && orig_m.cmpnt[2] < 0)q.cmpnt[2] *= -1;

    return std::make_tuple(v, q);
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
    result->alpha = std::min(alpha, result->alpha_max);
    result->set_volume();
    return result;
}

double volreconst::get_volume()
{
    return volume;
}

double volreconst::get_alpha()
{
    return alpha;
}

vector volreconst::get_moments()
{
    return {get_moment(0), get_moment(1), get_moment(2)};
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

    return volume * ((alpha / 2.0));
}

/******************** volreconst2d ********************/

double volreconst2d::x2(double x)
{
    return x < 0 ? 0 : x * x;
}

double volreconst2d::get_volume(double _alpha)
{
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

    for (size_t i = 0; i < 3; ++i)
        if (has_elem[i])
        {
            double d1 = i == dir ? c.cmpnt[i] : 0;
            double r1 = alpha - m.cmpnt[i] * c.cmpnt[i];
            result -= x2(r1) * ((d1 + r1 / m.cmpnt[dir]  / 3.0));
        }
    result *= base_vol;
    return result;
}

/******************** volreconst3d ********************/

double volreconst3d::x3(double x)
{
    return x < 0 ? 0 : x * x * x;
}

double volreconst3d::get_volume(double _alpha)
{
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

    for (size_t i = 0; i < 3; ++i)
    {
        double d1 = i == dir ? c.cmpnt[i] : 0;
        double r1 = alpha - m.cmpnt[i] * c.cmpnt[i];
        result -= x3(r1) * ((d1 + r1 / m.cmpnt[dir] / 3.0));

        double d2 = i == dir ? 0 : c.cmpnt[i];
        double r2 = alpha - alpha_max + m.cmpnt[i] * c.cmpnt[i];
        result += x3(r2) * ((d2 + r2 / m.cmpnt[dir] / 3.0));
    }

    result *= base_vol;
    return result;
}

}

