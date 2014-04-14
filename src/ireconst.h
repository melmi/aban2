/**
 * ireconst.h
 *
 * by Mohammad Elmi
 * Created on 2 Dec 2013
 */

/*
 * Reconstructs interface based on Gueyffier et al.
 *
 * Gueyffier, D., Li, J., Nadim, A., Scardovelli, R., and Zaleski, S. (1999),
 * "Volume-of-fluid interface tracking with smoothed surface stress methods for three-dimensional flows,"
 * Journal of Computational Physics, 152, 423–456.
 */

#ifndef _IRECONST_H_
#define _IRECONST_H_

#include "vector.h"

namespace aban2
{
class ireconst
{
public:
    typedef double(ireconst::*vol_func_t)(double x);

    double get_volume_3d(double _alpha);
    double get_volume_2d(double _alpha);
    double get_volume_1d(double _alpha);

    void init(vector _c, vector _m);

    void set_volume(); //sets volume assuming c, m and alpha  are known
    void set_alpha();  //sets alpha  assuming c, m and volume are known
    ireconst get_remaining(size_t dir, double delta);

    double alpha_max;
    double base_vol;
    static vol_func_t vol_funcs[];
    vol_func_t vol_func;
    vector orig_m;
public:
    constexpr static const double epsilon = 1e-7;

    vector c; // cell lengths
    vector m; // unit normal vector
    double alpha; // distance of freesurface from the most far corner of cell which is inside
    double volume;

    void split(size_t dir, double delta, ireconst &remaining, ireconst &departing);

    static ireconst from_full  (vector c, vector m);
    static ireconst from_empty (vector c, vector m);
    static ireconst from_volume(vector c, vector m, double volume);
    static ireconst from_alpha (vector c, vector m, double alpha );
};
}

#endif