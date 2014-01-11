/**
 * domain.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "domain.h"

namespace aban2
{

varinfo::varinfo(std::string _name, vardim _dim, bool _show, void *_data):
    name(_name), dim(_dim), show(_show)
{
    data.generic = _data;
}

domain *domain::create_from_file(std::string file_name)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream ifile(file_name);
    reader.parse(ifile, root, false);
    ifile.close();

    domain *result = new domain(&root);

    return result;
}

domain::domain(Json::Value *root): mesh(root)
{
    t = 0;
    dt = root->get("dt", 1).asDouble();
    tend = root->get("tend", 1).asDouble();
    rho = root->get("rho", 1000).asDouble();
    mu = root->get("mu", 1e-3).asDouble();
    Json::Value gnode = (*root)["g"];
    g = vector(gnode[0].asDouble(), gnode[1].asDouble(), gnode[2].asDouble());
    step_write = root->get("step_write", 1).asInt();

    boundaries = new bcondition*[256];
    Json::Value boundaries_val = (*root)["boundaries"];
    bcondition::create_bcs(&boundaries_val, boundaries, this);

    register_vars();
    create_vars();
}

void domain::register_vars()
{
    varlist.push_back(varinfo("vof", vardim::scalar, true, &vof));
    varlist.push_back(varinfo("smooth_vof", vardim::scalar, false, &smooth_vof));
    varlist.push_back(varinfo("nb", vardim::vector, false, &nb));

    varlist.push_back(varinfo("p", vardim::scalar, true, &p));
    varlist.push_back(varinfo("u", vardim::vector, true, &u));
    varlist.push_back(varinfo("uf", vardim::vector, false, &uf));
    varlist.push_back(varinfo("ustar", vardim::vector, true, &ustar));
}

void domain::write_vtk(std::string file_name)
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
                            if (exists(i, j, k)) val = (*v.data.scalar)[cellnos[x]];
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
                                val = vector::from_data(*v.data.vec, cellnos[x]);
                            file <<  val.x << " " << val.y << " " << val.z << std::endl;
                        }
                break;
            }
    file.close();
}

double *domain::extract_scalars(mesh_row *row, double *var)
{
    double *result = new double[row->n];
    for (size_t i = 0; i < row->n; ++i)
        result[i] = var[cellno(row, i)];
    return result;
}

void domain::insert_scalars(mesh_row *row, double *var, double *row_vals)
{
    for (size_t i = 0; i < row->n; ++i)
        var[cellno(row, i)] = row_vals[i];
}

vector *domain::extract_vectors(mesh_row *row, double **var)
{
    vector *result = new vector[row->n];
    for (size_t i = 0; i < row->n; ++i)
        result[i] = vector::from_data(var, cellno(row, i));
    return result;
}

void domain::insert_vectors(mesh_row *row, double **var, vector *row_vals)
{
    for (size_t i = 0; i < row->n; ++i)
    {
        size_t ix = cellno(row, i);
        var[0][ix] = row_vals[i].x;
        var[1][ix] = row_vals[i].y;
        var[2][ix] = row_vals[i].z;
    }
}

size_t *domain::get_row_idxs(mesh_row *row)
{
    size_t *row_idxs = new size_t[row->n];
    for (size_t i = 0; i < row->n; ++i)
        row_idxs[i] = cellno(row, i);
    return row_idxs;
}

void *domain::create_var(size_t dim)
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

void domain::delete_var(size_t dim, void *v)
{
    if (dim > 1)
        for (int i = 0; i < 3; ++i)
            delete_var(dim - 1, ((double **)v)[i]);
    delete[] (double *)v;
}

void domain::create_vars()
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

void domain::delete_vars()
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

}

