/**
 * bcondition.h
 *
 * by Mohammad Elmi
 * Created on 25 Oct 2013
 */

#ifndef _BCONDITION_H_
#define _BCONDITION_H_

#include <jsoncpp/json/json.h>

namespace aban2
{

class domain;
class mesh_row;

enum class bcside
{
    start, end
};

enum class bctype
{
    neumann, dirichlet
};

class bcondition
{
public:

    typedef double(bcondition::*func)(domain *, mesh_row *, bcside, size_t);
    typedef bctype bcondition::*type;

    bctype voftype, utype, ptype;

    virtual double p(domain *d, mesh_row *r, bcside side, size_t cmpnt) = 0;
    virtual double u(domain *d, mesh_row *r, bcside side, size_t cmpnt) = 0;
    virtual double vof(domain *d, mesh_row *r, bcside side, size_t cmpnt) = 0;

    static void create_bcs(Json::Value *bcroot, bcondition **boundaries);
};

class velocitybc: public bcondition
{
public:
    double value[3];

    double p(domain *d, mesh_row *r, bcside side, size_t cmpnt);
    double u(domain *d, mesh_row *r, bcside side, size_t cmpnt);
    double vof(domain *d, mesh_row *r, bcside side, size_t cmpnt);

    velocitybc(Json::Value *bcdata);
};

class pressurebc: public bcondition
{
public:
    double value;

    double p(domain *d, mesh_row *r, bcside side, size_t cmpnt);
    double u(domain *d, mesh_row *r, bcside side, size_t cmpnt);
    double vof(domain *d, mesh_row *r, bcside side, size_t cmpnt);

    pressurebc(Json::Value *bcdata);
};

}

#endif