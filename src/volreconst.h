/**
 * volreconst.h
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

#ifndef _VOLRECONST_H_
#define _VOLRECONST_H_

#include <tuple>
#include "vector.h"

namespace aban2
{

class volreconst
{
    constexpr static const double epsilon = 1e-7;

    static volreconst *from_base_data(vector _c, vector _m);
    void set_volume(); //sets volume assuming c, m and alpha  are known
    void set_alpha();  //sets alpha  assuming c, m and volume are known
protected:
    bool has_elem[3];
    double alpha_max;
    double base_vol;
    vector orig_m;

    vector c; // cell lengths
    vector m; // unit normal vector
    double alpha; // distance of freesurface from the most far corner of cell which is inside
    double volume;

    virtual double get_volume(double _alpha) = 0;
    volreconst *get_cut(size_t dir, double delta);
public:
    virtual ~volreconst() {}
    double get_volume();
    double get_alpha();
    vector get_moments();
    virtual double get_moment(size_t dir) = 0; // respect to cell cornet
    std::tuple<double, vector> get_flux(size_t dir, double delta);

    static volreconst *from_volume(vector c, vector m, double volume);
    static volreconst *from_alpha (vector c, vector m, double alpha );
};

class volreconst1d: public volreconst
{
protected:
    double get_volume(double _alpha);
public:
    double get_moment(size_t dir);
};

class volreconst2d: public volreconst
{
    static inline double x2(double x);
protected:
    double get_volume(double _alpha);
public:
    double get_moment(size_t dir);
};

class volreconst3d: public volreconst
{
    static inline double x3(double x);
protected:
    double get_volume(double _alpha);
public:
    double get_moment(size_t dir);
};

}

#endif
