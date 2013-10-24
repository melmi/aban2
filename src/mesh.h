/**
 * mesh.h
 *
 * by Mohammad Elmi
 * Created on 26 Aug 2013
 */

#ifndef _MESH_H_
#define _MESH_H_

#include <jsoncpp/json/json.h>
#include <string>
#include <cstring>
#include <fstream>
#include <list>
#include "vector.h"

//#define THREE_D
#ifdef THREE_D
#define NDIRS 3
#else
#define NDIRS 2
#endif

namespace aban2
{

struct mesh_row
{
    size_t dir;
    size_t n;
    size_t start[3], end[3];
    char start_code, end_code;
};

class mesh
{
public:
    size_t n, ndir[3];
    double delta;
    char *codes;
    size_t *cellnos;
    size_t nrows[3];
    static const char INSIDE = ' ';
    mesh_row *rows[3];

    inline size_t idx(size_t i, size_t j, size_t k)
    {
        return k * ndir[1] + j * ndir[0] + i;
    }

    inline size_t idx(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk)
    {
        static size_t ijk[3];
        ijk[ii] = i;
        ijk[jj] = j;
        ijk[kk] = k;
        return idx(ijk[0], ijk[1], ijk[2]);
    }

    inline size_t cellno(mesh_row *r, size_t i)
    {
        r->start[r->dir] += i;
        size_t result = cellnos[idx(r->start[0], r->start[1], r->start[2])];
        r->start[r->dir] -= i;
        return result;
    }

    inline void get_dirs(size_t dir, size_t &ii, size_t &jj, size_t &kk)
    {
        ii = dir;
        jj = (dir + 1) % 3;
        kk = (dir + 2) % 3;
    }

    bool exists(size_t i, size_t j, size_t k)
    {
        size_t x = idx(i, j, k);
        return codes[x] == INSIDE;
    }

private:
    void generate_rows(size_t dir)
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

    void set_cell_nos()
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

    void init()
    {
        nrows[0] = nrows[1] = nrows[2] = 0;

        for (int dir = 0; dir < NDIRS; ++dir)
            generate_rows(dir);

        set_cell_nos();
    }

public:
    mesh(Json::Value *root)
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
};

}

#endif
