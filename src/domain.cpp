/**
 * domain.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "domain.h"

namespace aban2
{

domain_t::varinfo::varinfo(std::string _name, char _rank, bool _show, void *_data)
    : name(_name), rank(_rank), show(_show)
{
    data.generic = _data;
}

domain_t *domain_t::create_from_file(std::string file_name)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream ifile(file_name);
    reader.parse(ifile, root, false);
    ifile.close();

    domain_t *result = new domain_t(&root);

    return result;
}

domain_t::domain_t(Json::Value *root)
{
    delta = root->get("delta", "1").asDouble();

    ndir[0] = root->get("ni", "1").asInt();
    ndir[1] = root->get("nj", "1").asInt();
    ndir[2] = root->get("nk", "1").asInt();

    nrows[0] = nrows[1] = nrows[2] = 0;

    for (size_t dir = 0; dir < NDIRS; ++dir)
        generate_rows(dir);

    set_cell_nos();

    vcell = delta * delta * delta;
    aface = delta * delta;

    t = 0;
    dt = root->get("dt", 1).asDouble();
    tend = root->get("tend", 1).asDouble();
    rho0 = root->get("rho0", 1000).asDouble();
    rho1 = root->get("rho1", 1000).asDouble();
    mu0 = root->get("mu0", 1e-3).asDouble();
    mu1 = root->get("mu1", 1e-3).asDouble();
    Json::Value gnode = (*root)["g"];
    g = vector(gnode[0].asDouble(), gnode[1].asDouble(), gnode[2].asDouble());
    write_interval = root->get("write_interval", 1).asInt();

    boundaries = new flowbc *[256];
    std::fill_n(boundaries, 256, nullptr);
    Json::Value boundaries_val = (*root)["boundaries"];
    flowbc::create_bcs(&boundaries_val, boundaries, this);

    register_vars();
    create_vars();
}

domain_t::~domain_t()
{
    for (int i = 0; i < 256; ++i)
        if (boundaries[i] != nullptr)
            delete boundaries[i];
    delete[] boundaries;
    delete_vars();

    delete[] cellnos;
    for (int i = 0; i < NDIRS; ++i)
        delete[] rows[i];
}

void domain_t::register_vars()
{
    varlist.push_back(varinfo("p", 1, true, &p));
    varlist.push_back(varinfo("u", 2, true, &u));
    varlist.push_back(varinfo("ustar", 2, false, &ustar));
    varlist.push_back(varinfo("rho", 1, false, &rho));
    varlist.push_back(varinfo("nu", 1, false, &nu));
    varlist.push_back(varinfo("uf", 2, false, &uf));

    varlist.push_back(varinfo("vof", 1, true, &vof));
    varlist.push_back(varinfo("nb", 2, false, &nb));
}

void domain_t::write_vtk(std::string file_name)
{
    std::ofstream file(file_name);

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "aban2 output" << std::endl;
    file << "ASCII" << std::endl;

    file << "DATASET RECTILINEAR_GRID" << std::endl;
    file << "DIMENSIONS " << ndir[0] + 1 << " " << ndir[1] + 1 << " " << ndir[2] + 1 << std::endl;
    file << "X_COORDINATES " << ndir[0] + 1 << " float" << std::endl;
    for (size_t i = 0; i < ndir[0] + 1; ++i)
        file << (double)i * delta << " ";
    file << std::endl;
    file << "Y_COORDINATES " << ndir[1] + 1 << " float" << std::endl;
    for (size_t i = 0; i < ndir[1] + 1; ++i)
        file << (double)i * delta << " ";
    file << std::endl;
    file << "Z_COORDINATES " << ndir[2] + 1 << " float" << std::endl;
    for (size_t i = 0; i < ndir[2] + 1; ++i)
        file << (double)i * delta << " ";
    file << std::endl;

    file << "CELL_DATA " << ndir[0] * ndir[1] * ndir[2] << std::endl;

    for (auto &v : varlist)
        if (v.show)
            switch (v.rank)
            {
            case 1:
                file << "SCALARS " << v.name << " double" << std::endl;
                file << "LOOKUP_TABLE default" << std::endl;
                for (size_t k = 0; k < ndir[2]; ++k)
                    for (size_t j = 0; j < ndir[1]; ++j)
                        for (size_t i = 0; i < ndir[0]; ++i)
                        {
                            size_t x = idx(i, j, k);
                            double val = 0;
                            if (is_inside(i, j, k))
                                val = (*v.data.scalar)[cellnos[x]];
                            file << val << std::endl;
                        }
                break;
            case 2:
                file << "VECTORS " << v.name << " double" << std::endl;
                for (size_t k = 0; k < ndir[2]; ++k)
                    for (size_t j = 0; j < ndir[1]; ++j)
                        for (size_t i = 0; i < ndir[0]; ++i)
                        {
                            size_t x = idx(i, j, k);
                            vector val;
                            if (is_inside(i, j, k))
                                val = vector::from_data(*v.data.vec, cellnos[x]);
                            file << val.x << " " << val.y << " " << val.z << std::endl;
                        }
                break;
            }
    file.close();
}

double *domain_t::extract_scalars(row_t *row, double *var)
{
    double *result = new double[row->n];
    for (size_t i = 0; i < row->n; ++i)
        result[i] = var[cellno(row, i)];
    return result;
}

void domain_t::insert_scalars(row_t *row, double *var, double *row_vals)
{
    for (size_t i = 0; i < row->n; ++i)
        var[cellno(row, i)] = row_vals[i];
}

vector *domain_t::extract_vectors(row_t *row, double **var)
{
    vector *result = new vector[row->n];
    for (size_t i = 0; i < row->n; ++i)
        result[i] = vector::from_data(var, cellno(row, i));
    return result;
}

void domain_t::insert_vectors(row_t *row, double **var, vector *row_vals)
{
    for (size_t i = 0; i < row->n; ++i)
    {
        size_t ix = cellno(row, i);
        var[0][ix] = row_vals[i].x;
        var[1][ix] = row_vals[i].y;
        var[2][ix] = row_vals[i].z;
    }
}

size_t *domain_t::get_row_cellnos(row_t *row)
{
    size_t *row_idxs = new size_t[row->n];
    for (size_t i = 0; i < row->n; ++i)
        row_idxs[i] = cellno(row, i);
    return row_idxs;
}

void *domain_t::create_var(size_t rank)
{
    if (rank == 1)
    {
        double *result = new double[n];
        std::fill_n(result, n, 0);
        return result;
    }
    else
    {
        void **result = new void *[3];
        for (int i = 0; i < 3; ++i)
            result[i] = create_var(rank - 1);
        return result;
    }
}

void domain_t::delete_var(size_t rank, void *v)
{
    if (v == nullptr)
        return;
    if (rank > 1)
        for (int i = 0; i < 3; ++i)
            delete_var(rank - 1, ((double **)v)[i]);
    delete[](double *) v;
}

void domain_t::create_vars()
{
    for (auto &v : varlist)
    {
        switch (v.rank)
        {
        case 1:
            *v.data.scalar = (double *)create_var(1);
            break;
        case 2:
            *v.data.vec = (double **)create_var(2);
            break;
        }
    }
}

void domain_t::delete_vars()
{
    for (auto &v : varlist)
    {
        switch (v.rank)
        {
        case 1:
            delete_var(1, *v.data.scalar);
            break;
        case 2:
            delete_var(2, *v.data.vec);
            break;
        }
    }
}

double domain_t::rho_bar(double _vof)
{
    return _vof * rho1 + (1.0 - _vof) * rho0;
}

double domain_t::nu_bar(double _vof, double rho_bar)
{
    double mu_bar=1.0/(_vof / mu1 + (1.0 - _vof) / mu0);
    return mu_bar/rho_bar;
}

/////////////////////////////////////////////////////////////

size_t domain_t::idx(size_t i, size_t j, size_t k)
{
    return k * ndir[1] + j * ndir[0] + i;
}

size_t domain_t::cellno(row_t *row, size_t i)
{
    row->start[row->dir] += i;
    size_t result = cellnos[idx(row->start[0], row->start[1], row->start[2])];
    row->start[row->dir] -= i;
    return result;
}

size_t domain_t::idx(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk)
{
    static size_t ijk[3];
    ijk[ii] = i;
    ijk[jj] = j;
    ijk[kk] = k;
    return idx(ijk[0], ijk[1], ijk[2]);
}

void get_dirs(size_t dir, size_t &ii, size_t &jj, size_t &kk)
{
    ii = dir;
    jj = (dir + 1) % 3;
    kk = (dir + 2) % 3;
}

bool domain_t::is_inside(size_t i, size_t j, size_t k)
{
    size_t dummy;
    return is_inside(i, j, k, 0, 0, 0, dummy);
}

bool domain_t::is_inside(size_t i, size_t j, size_t k, size_t &no)
{
    return is_inside(i, j, k, 0, 0, 0, no);
}

bool inside_bounds(size_t x, int dx, size_t max)
{
    long xx = x;
    xx += dx;
    return xx >= 0 && xx < (long)max;
}

bool domain_t::is_inside(size_t i, size_t j, size_t k, int di, int dj, int dk, size_t &no)
{
    bool inside = inside_bounds(i, di, ndir[0]) && inside_bounds(j, dj, ndir[1]) && inside_bounds(k, dk, ndir[2]);
    if (!inside)
        return false;

    i += di;
    j += dj;
    k += dk;

    size_t ix = idx(i, j, k);
    no = cellnos[ix];
    return true;
}

void domain_t::generate_rows(size_t dir)
{
    row_t row;
    std::vector<row_t> v;
    size_t ii, jj, kk;
    get_dirs(dir, ii, jj, kk);

    for (size_t k = 0; k < ndir[kk]; ++k)
        for (size_t j = 0; j < ndir[jj]; ++j)
        {
            row.dir = dir;
            row.start[ii] = 0;
            row.start[jj] = j;
            row.start[kk] = k;
            row.start_code = 1;
            row.end[ii] = ndir[ii] - 1;
            row.end[jj] = j;
            row.end[kk] = k;
            row.n = row.end[ii] - row.start[ii] + 1;
            row.end_code = 1;
            v.push_back(row);
        }

    nrows[dir] = v.size();
    rows[dir] = new struct row_t[nrows[dir]];
    for (size_t i = 0; i < nrows[dir]; ++i)
        rows[dir][i] = v[i];
}

void domain_t::set_cell_nos()
{
    size_t ii = 0;
    cellnos = new size_t[ndir[0] * ndir[1] * ndir[2]];
    for (size_t i = 0; i < ndir[0]; ++i)
        for (size_t j = 0; j < ndir[1]; ++j)
            for (size_t k = 0; k < ndir[2]; ++k)
            {
                size_t x = idx(i, j, k);
                cellnos[x] = ii++;
            }
    n = ii;
}
}
