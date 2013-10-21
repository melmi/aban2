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

    void set_vel_bc(vector v)
    {
        pbc.type = bctype::neumann;
        pbc.val = 0;

        for (int i = 0; i < 3; ++i)
        {
            velbc[i].type = bctype::dirichlet;
            velbc[i].val = v.components[i];
        }
    }

    void set_pressure_bc(double val)
    {
        pbc.type = bctype::dirichlet;
        pbc.val = val;

        for (int i = 0; i < 3; ++i)
        {
            velbc[i].type = bctype::neumann;
            velbc[i].val = 0;
        }
    }
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

    varinfo(std::string _name, vardim _dim, bool _show, void *_data):
        name(_name), dim(_dim), show(_show)
    {
        data.generic = _data;
    }
};

class domain: public mesh
{
public:
    double **u, * *ustar, *p;
    double dt, tend, rho, mu;
    int step_write;
    flow_boundary *boundaries;
    std::list<varinfo> varlist;

    static domain *create_from_file(std::string file_name)
    {
        Json::Reader reader;
        Json::Value root;

        std::ifstream ifile(file_name);
        reader.parse(ifile, root, false);
        ifile.close();

        domain *result = new domain(&root);
        // msh = mesh::from_json(root);

        return result;
    }

    domain(Json::Value *root): mesh(root)
    {
        dt = root->get("dt", 1).asDouble();
        tend = root->get("tend", 1).asDouble();
        rho = root->get("rho", 1000).asDouble();
        mu = root->get("mu", 1e-3).asDouble();
        step_write = root->get("step_write", 1).asInt();

        boundaries = new flow_boundary[256];
        Json::Value boundaries_val = (*root)["boundaries"];
        auto bnames = boundaries_val.getMemberNames();
        for (auto & x : bnames)
        {
            Json::Value bc = boundaries_val[x];
            auto boundary = boundaries + x[0];
            if (bc["type"].asString() == "flux")
            {
                Json::Value bcval = bc["value"];
                boundary->set_vel_bc(vector(bcval[0].asDouble(), bcval[1].asDouble(), bcval[2].asDouble()));
            }
            if (bc["type"].asString() == "pressure")
                boundary->set_pressure_bc(bc.get("value", "0").asDouble());
        }

        register_vars();
        create_vars();
    }

    virtual void register_vars()
    {
        varlist.push_back(varinfo("p", vardim::scalar, true, &p));
        varlist.push_back(varinfo("u", vardim::vector, true, &u));
        varlist.push_back(varinfo("ustar", vardim::vector, true, &ustar));
    }

    void write_vtk(std::string file_name)
    {
        std::ofstream file(file_name);

        file << "# vtk DataFile Version 2.0" << std::endl;
        file << "aban2 output" << std::endl;
        file << "ASCII" << std::endl;

        file << "DATASET RECTILINEAR_GRID" << std::endl;
        file << "DIMENSIONS " << ndir[0] + 1 << " " << ndir[1] + 1 << " " << ndir[2] + 1 << std::endl;
        file << "X_COORDINATES " << ndir[0] + 1 << " float" << std::endl;
        for (int i = 0; i < ndir[0] + 1; ++i)file << (double)i *delta << " ";
        file << std::endl;
        file << "Y_COORDINATES " << ndir[1] + 1 << " float" << std::endl;
        for (int i = 0; i < ndir[1] + 1; ++i)file << (double)i *delta << " ";
        file << std::endl;
        file << "Z_COORDINATES " << ndir[2] + 1 << " float" << std::endl;
        for (int i = 0; i < ndir[2] + 1; ++i)file << (double)i *delta << " ";
        file << std::endl;

        file << "CELL_DATA " << ndir[0] *ndir[1] *ndir[2] << std::endl;
        file << "SCALARS existence int" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for (int k = 0; k < ndir[2]; ++k)
            for (int j = 0; j < ndir[1]; ++j)
                for (int i = 0; i < ndir[0]; ++i)
                    file << (exists(i, j, k) ? 1 : 0) << std::endl;

        for (auto & v : varlist)
            if (v.show)
                switch (v.dim)
                {
                case vardim::scalar:
                    file << "SCALARS " << v.name << " double" << std::endl;
                    file << "LOOKUP_TABLE default" << std::endl;
                    for (int k = 0; k < ndir[2]; ++k)
                        for (int j = 0; j < ndir[1]; ++j)
                            for (int i = 0; i < ndir[0]; ++i)
                            {
                                size_t x = idx(i, j, k);
                                double val = 0;
                                if (exists(i, j, k)) val = (*v.data.scalar)[idxs[x]];
                                file <<  val << std::endl;
                            }
                    break;
                case vardim::vector:
                    file << "VECTORS " << v.name << " double" << std::endl;
                    for (int k = 0; k < ndir[2]; ++k)
                        for (int j = 0; j < ndir[1]; ++j)
                            for (int i = 0; i < ndir[0]; ++i)
                            {
                                size_t x = idx(i, j, k);
                                vector val;
                                if (exists(i, j, k))
                                    val = vector::from_data(*v.data.vec, idxs[x]);
                                file <<  val.x << " " << val.y << " " << val.z << std::endl;
                            }
                    break;
                }
        file.close();
    }

    double *extract_scalars(mesh_row *row, double *var)
    {
        double *result = new double[row->n];
        size_t s = row->start[row->ii], e = row->end[row->ii];
        for (size_t i = s; i <= e; ++i)
        {
            size_t ix = idxs[idx(i, row->start[row->jj], row->start[row->kk], row->ii, row->jj, row->kk)];
            result[i - s] = var[ix];
        }
        return result;
    }

    void insert_scalars(mesh_row *row, double *var, double *row_vals)
    {
        size_t s = row->start[row->ii], e = row->end[row->ii];
        for (size_t i = s; i <= e; ++i)
        {
            size_t ix = idxs[idx(i, row->start[row->jj], row->start[row->kk], row->ii, row->jj, row->kk)];
            var[ix] = row_vals[i - s];
        }
    }

    vector *extract_vectors(mesh_row *row, double **var)
    {
        vector *result = new vector[row->n];
        size_t s = row->start[row->ii], e = row->end[row->ii];
        for (size_t i = s; i <= e; ++i)
        {
            size_t ix = idxs[idx(i, row->start[row->jj], row->start[row->kk], row->ii, row->jj, row->kk)];
            result[i - s] = vector::from_data(var, ix);
        } return result;
    }

    void insert_vectors(mesh_row *row, double **var, vector *row_vals)
    {
        size_t s = row->start[row->ii], e = row->end[row->ii];
        for (size_t i = s; i <= e; ++i)
        {
            size_t ix = idxs[idx(i, row->start[row->jj], row->start[row->kk], row->ii, row->jj, row->kk)];
            var[0][ix] = row_vals[i - row->start[row->ii]].x;
            var[1][ix] = row_vals[i - row->start[row->ii]].y;
            var[2][ix] = row_vals[i - row->start[row->ii]].z;
        }
    }

    size_t *get_row_idxs(mesh_row *row)
    {
        size_t *row_idxs = new size_t[row->n];
        size_t s = row->start[row->ii], e = row->end[row->ii];

        for (size_t i = s; i <= e; ++i)
        {
            size_t ix = idxs[idx(i, row->start[row->jj], row->start[row->kk], row->ii, row->jj, row->kk)];
            row_idxs[i - s] = ix;
        }

        return row_idxs;
    }

    void *create_var(size_t dim)
    {
        if (dim == 1)
        {
            double *result = new double[n];
            std::fill_n(result, n, 0);
            return result;
        }
        else
        {
            void **result = new void*[3];
            for (int i = 0; i < 3; ++i)
                result[i] = create_var(dim - 1);
            return result;
        }
    }

    void delete_var(size_t dim, void *v)
    {
        if (dim > 1)
            for (int i = 0; i < 3; ++i)
                delete_var(dim - 1, ((double **)v)[i]);
        delete[] v;
    }

private:
    void create_vars()
    {
        for (auto & v : varlist)
        {
            switch (v.dim)
            {
            case vardim::scalar:
                *v.data.scalar = (double *)create_var(1);
                break;
            case vardim::vector:
                *v.data.vec = (double **)create_var(2);
                break;
            }
        }
    }

    void delete_vars()
    {
        for (auto & v : varlist)
        {
            switch (v.dim)
            {
            case vardim::scalar:
                delete_var(1, *v.data.scalar);
                break;
            case vardim::vector:
                delete_var(2, *v.data.vec);
                break;
            }
        }
    }
};

}

#endif
