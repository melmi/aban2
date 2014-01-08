/**
 * bcondition.h
 *
 * by Mohammad Elmi
 * Created on 25 Oct 2013
 */

#ifndef _BCONDITION_H_
#define _BCONDITION_H_

#include <jsoncpp/json/json.h>
#include "vector.h"

namespace aban2
{

class domain;
class mesh_row;

enum class bcside
{
    start, end
};

struct bcdesc
{
    // phi_face = coeff * phi[cellno] + val

    size_t cellno;
    double coeff;
    double val;
};

class bcondition
{
public:
    typedef bcdesc(bcondition::*func)(mesh_row *, bcside, size_t);

    static void create_bcs(Json::Value *bcroot, bcondition **boundaries, domain *_d);

    double face_val(func f, double *phi, mesh_row *r, bcside side, size_t cmpnt);
    double row_face_val(func f, double *phi, mesh_row *r, bcside side, size_t cmpnt);

    domain *d;

    virtual bcdesc p(mesh_row *r, bcside side, size_t cmpnt) = 0;
    virtual bcdesc u(mesh_row *r, bcside side, size_t cmpnt) = 0;
    virtual bcdesc vof(mesh_row *r, bcside side, size_t cmpnt) = 0;

    bcondition(domain *_d): d(_d) {}
};

class velocitybc: public bcondition
{
public:
    vector value;

    virtual bcdesc p(mesh_row *r, bcside side, size_t cmpnt);
    virtual bcdesc u(mesh_row *r, bcside side, size_t cmpnt);
    virtual bcdesc vof(mesh_row *r, bcside side, size_t cmpnt);

    velocitybc(Json::Value *bcdata, domain *_d);
};

class pressurebc: public bcondition
{
public:
    pressurebc(domain *_d): bcondition(_d) {}

    double value;

    virtual bcdesc p(mesh_row *r, bcside side, size_t cmpnt);
    virtual bcdesc u(mesh_row *r, bcside side, size_t cmpnt);
    virtual bcdesc vof(mesh_row *r, bcside side, size_t cmpnt);

    pressurebc(Json::Value *bcdata, domain *_d);
};

}

#endif