/**
 * mesh.cpp
 *
 * by Mohammad Elmi
 * Created on 24 Oct 2013
 */

#include "mesh.h"

namespace aban2
{

size_t mesh::idx(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk)
{
    static size_t ijk[3];
    ijk[ii] = i;
    ijk[jj] = j;
    ijk[kk] = k;
    return idx(ijk[0], ijk[1], ijk[2]);
}

void mesh::get_dirs(size_t dir, size_t &ii, size_t &jj, size_t &kk)
{
    ii = dir;
    jj = (dir + 1) % 3;
    kk = (dir + 2) % 3;
}

bool mesh::exists(size_t i, size_t j, size_t k)
{
    size_t x = idx(i, j, k);
    return codes[x] == INSIDE;
}

void mesh::generate_rows(size_t dir)
{
    bool inside = false;
    mesh_row row;
    std::vector<mesh_row> v;
    size_t ii , jj, kk;
    get_dirs(dir, ii, jj, kk);

    for (size_t k = 0; k < ndir[kk]; ++k)
        for (size_t j = 0; j < ndir[jj]; ++j)
            for (size_t i = 0; i < ndir[ii]; ++i)
            {
                size_t x = idx(i, j, k, ii, jj, kk);
                if (!inside)
                {
                    if (codes[x] != INSIDE) continue;
                    inside = true;
                    row.dir = dir;
                    row.start[ii] = i;
                    row.start[jj] = j;
                    row.start[kk] = k;
                    row.start_code = codes[idx(i - 1, j, k, ii, jj, kk)];
                }
                else
                {
                    if (codes[x] == INSIDE) continue;
                    inside = false;
                    row.end[ii] = i - 1;
                    row.end[jj] = j;
                    row.end[kk] = k;
                    row.n = row.end[ii] - row.start[ii] + 1;
                    row.end_code = codes[x];
                    v.push_back(row);
                }
            }

    nrows[dir] = v.size();
    rows[dir] = new mesh_row[nrows[dir]];
    for (int i = 0; i < nrows[dir]; ++i)
        rows[dir][i] = v[i];
}

void mesh::set_cell_nos()
{
    size_t ii = 0;
    cellnos = new size_t[ndir[0] * ndir[1] * ndir[2]];
    for (int i = 0; i < ndir[0]; ++i)
        for (int j = 0; j < ndir[1]; ++j)
            for (int k = 0; k < ndir[2]; ++k)
            {
                size_t x = idx(i, j, k);
                if (codes[x] == INSIDE)
                    cellnos[x] = ii++;
            }
    n = ii;
}

void mesh::init()
{
    nrows[0] = nrows[1] = nrows[2] = 0;

    for (int dir = 0; dir < NDIRS; ++dir)
        generate_rows(dir);

    set_cell_nos();
}

mesh::mesh(Json::Value *root)
{
    delta = root->get("delta", "1").asDouble();

    ndir[0] = root->get("ni", "1").asInt();
    ndir[1] = root->get("nj", "1").asInt();
    ndir[2] = root->get("nk", "1").asInt();

    std::string str;
    Json::Value codes_val = (*root)["codes"];
    for (size_t k = 0; k < ndir[2]; ++k)
    {
        Json::Value codesk = codes_val[ndir[2] - k - 1];
        for (size_t j = 0; j < ndir[1]; ++j)
            str += codesk[ndir[1] - j - 1].asString();
    }

    codes = new char[str.length() + 1];
    strcpy(codes, str.c_str());

    init();
}

}