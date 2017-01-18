/**
 * domain.h
 *
 * by Mohammad Elmi
 * Created on 26 Aug 2013
 */

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "common.h"
#include <jsoncpp/json/json.h>
#include <list>
#include <fstream>

#include "bcondition.h"

namespace aban2
{

    struct row_t
    {
        size_t dir;
        size_t n;
        size_t start[3], end[3];
        char start_code, end_code;
    };

class domain_t
{
    struct varinfo
    {
        std::string name;
        char rank;
        bool show;
        union
        {
            void *generic;
            double **scalar;
            double ***vec;
        } data;

        varinfo(std::string _name, char _rank, bool _show, void *_data);
    };
public:
    //came from old mesh
    size_t n, ndir[3];
    double delta;
    double vcell, aface;
    size_t nrows[3];
    size_t *cellnos;
    row_t *rows[3];

    size_t idx(size_t i, size_t j, size_t k);
    size_t cellno(row_t *row, size_t i);
    bool is_inside(size_t i, size_t j, size_t k);
    bool is_inside(size_t i, size_t j, size_t k, size_t &no);
    bool is_inside(size_t i, size_t j, size_t k, int di, int dj, int dk, size_t &no);
    //

    double **uf, * *u, * *ustar, *p;
    double rho0, rho1, nu0, nu1;
    double *rho, *nu;
    double *vof, * *nb;
    double dt, tend;
    double t;
    vector g;
    int write_interval;
    flowbc **boundaries;
    std::list<varinfo> varlist;

    static domain_t *create_from_file(std::string file_name);
    domain_t(Json::Value *root);
    ~domain_t();
    void register_vars();
    void write_vtk(std::string file_name);
    double *extract_scalars(row_t *row, double *var);
    void insert_scalars(row_t *row, double *var, double *row_vals);
    vector *extract_vectors(row_t *row, double **var);
    void insert_vectors(row_t *row, double **var, vector *row_vals);
    size_t *get_row_cellnos(row_t *row);
    void *create_var(size_t rank);
    static void delete_var(size_t rank, void *v);

    double rho_bar(double _vof);
    double nu_bar(double _vof);

private:
    inline size_t idx(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk);
    void generate_rows(size_t dir);
    void set_cell_nos();

    void create_vars();
    void delete_vars();
};

}

#endif
