/**
 * vof.h
 *
 * by Mohammad Elmi
 * Created on 10 Nov 2013
 */

#ifndef _VOF_H_
#define _VOF_H_

#include "common.h"
#include <array>
#include <tuple>
#include "domain.h"
#include "volreconst.h"
#include "fsnorm.h"

namespace aban2
{

class vof_t
{
    enum class fullness
    {
        full, empty, half
    };

    domain_t *d;
    fsnorm _fsnorm;
    size_t start_dir = 0;
    double *mass;
    double *rhou0[3], *rhou1[3];
    double *original_vof;
    fullness *fullnesses;
    bool *on_interface;
    volreconst **reconsts;

    void set_fullnesses();
    inline bool is_empty(size_t i, size_t j, size_t k, int di, int dj, int dk);
    bool is_on_interface(size_t i, size_t j, size_t k, size_t no);
    void detect_interfacial_cells();
    // void calculate_normals();
    void create_reconsts();
    void delete_reconsts();
    inline std::tuple<double, vector, vector> get_flux(row_t *row, size_t i, double udt, double ***grad_ustar);
    inline std::tuple<double, vector, vector> get_bc_flux(row_t *row, double ***grad_ustar, bcside side);
    void advect_row(row_t *row, double ***grad_ustar);
    void calculate_vof_masses_from_vars();
    void calculate_ustar_masses_from_vars(double ***grad_ustar);
    void calculate_vars_from_masses(size_t dir);
public:
    void calculate_normals();
    vof_t(domain_t *domain);
    ~vof_t();

    void advect();

    friend class vof_err_t;
};

class vof_err_t
{
    inline static int is_inside(vector n, double alpha, vector x);
    inline static double err(vector n1, double alpha1, vector n2, double alpha2, double h);
public:
    static double compare(vof_t &v);
};
}

#endif

