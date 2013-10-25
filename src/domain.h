/**
 * domain.h
 *
 * by Mohammad Elmi
 * Created on 26 Aug 2013
 */

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "mesh.h"
#include <jsoncpp/json/json.h>
#include <list>
#include <fstream>
#include <eigen3/Eigen/SparseCore>

namespace aban2
{

enum class bctype
{
    neumann, dirichlet
};

struct bcond
{
    bctype type;
    double val;
};

struct flow_boundary
{
    bcond velbc[3];
    bcond pbc;

    void set_vel_bc(vector v);
    void set_pressure_bc(double val);
};

enum class vardim
{
    scalar, vector
};

struct varinfo
{
    std::string name;
    vardim dim;
    bool show;
    union
    {
        void *generic;
        double **scalar;
        double ** *vec;
    } data;

    varinfo(std::string _name, vardim _dim, bool _show, void *_data);
};

class domain: public mesh
{
public:
    double **u, * *ustar, *p;
    double dt, tend, rho, mu;
    int step_write;
    flow_boundary *boundaries;
    std::list<varinfo> varlist;

    static domain *create_from_file(std::string file_name);

    domain(Json::Value *root);

    virtual void register_vars();

    void write_vtk(std::string file_name);

    double *extract_scalars(mesh_row *row, double *var);

    void insert_scalars(mesh_row *row, double *var, double *row_vals);

    vector *extract_vectors(mesh_row *row, double **var);

    void insert_vectors(mesh_row *row, double **var, vector *row_vals);

    size_t *get_row_idxs(mesh_row *row);

    void *create_var(size_t dim);

    void delete_var(size_t dim, void *v);

private:
    void create_vars();

    void delete_vars();
};

}

#endif
