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
//#include <eigen3/Eigen/SparseCore>

#include "mesh.h"
#include "bcondition.h"

namespace aban2
{

class domain: public mesh
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

    static domain *create_from_file(std::string file_name);
    domain(Json::Value *root);
    ~domain();
    virtual void register_vars();
    void write_vtk(std::string file_name);
    double *extract_scalars(row *row, double *var);
    void insert_scalars(row *row, double *var, double *row_vals);
    vector *extract_vectors(row *row, double **var);
    void insert_vectors(row *row, double **var, vector *row_vals);
    size_t *get_row_cellnos(row *row);
    void *create_var(size_t rank);
    static void delete_var(size_t rank, void *v);

    double rho_bar(double _vof);
    double nu_bar(double _vof);

private:
    void create_vars();
    void delete_vars();
};

}

#endif
