#include "ireconst.h"

#include <cmath>
#include <iostream>

namespace aban2
{
    inline double x2(double x) {return x * x;}
    inline double x3(double x) {return x * x * x;}

    ireconst::vol_func_t ireconst::vol_funcs[] {&ireconst::get_volume_1d, &ireconst::get_volume_2d, &ireconst::get_volume_3d};

    double ireconst::get_volume_3d(double _alpha)
    {
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
        double result = x2(_alpha);

        if (m.x > epsilon) if (_alpha - m.x * c.x > 0)result -= x2(_alpha - m.x * c.x);
        if (m.y > epsilon) if (_alpha - m.y * c.y > 0)result -= x2(_alpha - m.y * c.y);
        if (m.z > epsilon) if (_alpha - m.z * c.z > 0)result -= x2(_alpha - m.z * c.z);

        return result * base_vol;
    }

    double ireconst::get_volume_1d(double _alpha)
    {
        return _alpha * base_vol;
    }

    void ireconst::init(vector _c, vector _m)
    {
        c = _c;
        m = {std::abs(_m.x), std::abs(_m.y), std::abs(_m.z)};
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

    double ireconst::get_flux(size_t dir, double delta, bool from_start)
    {
        vector new_c = c;
        new_c.components[dir] -= delta;

        double new_alpha;
        if (from_start)
            new_alpha = alpha - delta * m.components[dir];
        else
            new_alpha = alpha;

        ireconst i2 = ireconst::from_alpha(new_c, m, new_alpha);
        return volume - i2.volume;
    }

    double ireconst::get_flux(size_t dir, double delta, vector orig_m)
    {
        if (delta * orig_m.components[dir] > 0)
            return get_flux(dir, std::abs(delta), true);
        else
            return get_flux(dir, std::abs(delta), false);
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
}