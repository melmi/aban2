/**
 * mesh.h
 *
 * by Mohammad Elmi
 * Created on 26 Aug 2013
 */

#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include <jsoncpp/json/json.h>
#include <string>
#include <cstring>
#include <fstream>
#include <list>
#include "vector.h"

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
    inline size_t cellno(mesh_row *r, size_t i)
    {
        r->start[r->dir] += i;
        size_t result = cellnos[idx(r->start[0], r->start[1], r->start[2])];
        r->start[r->dir] -= i;
        return result;
    }
    bool exists(size_t i, size_t j, size_t k);
    bool exists(size_t i, size_t j, size_t k, size_t &no);
    bool exists_and_inside(size_t i, size_t j, size_t k, size_t &no);

private:
    inline size_t idx(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk);
    inline void get_dirs(size_t dir, size_t &ii, size_t &jj, size_t &kk);
    void generate_rows(size_t dir);
    void set_cell_nos();
    void init();

public:
    mesh(Json::Value *root);
};

}

#endif
