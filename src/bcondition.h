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

struct bcdesc
{
    // phi_f = sw * (phi_P + (r_f-r_P) * grad) + cte
    // sw should be 0 or 1
    double sw, grad, cte;
    double face_val(double phi_P, double dx);
};

class bcondition
{
protected:
    domain *d;
public:
    bcondition(domain *_d): d(_d) {}
    virtual bcdesc desc(size_t cellno, size_t dir) = 0;
    double face_val(size_t cellno, size_t dir, double *phi, double dx);

    friend class flowbc;
};

class flowbc
{
public:
    typedef bcondition *flowbc::*member;

    union
    {
        struct
        {
            bcondition *u0, *u1, *u2;
        };
        bcondition *u[3];
    };
    bcondition *p, *vof;
    static member umembers[3];

    static void create_bcs(Json::Value *bcroot, flowbc **boundaries, domain *_d);
    flowbc(bcondition *_u0, bcondition *_u1, bcondition *_u2, bcondition *_p, bcondition *_vof)
        : u0(_u0), u1(_u1), u2(_u2), p(_p), vof(_vof) {}
    ~flowbc();

    static double fval_start(domain *d, mesh_row *row, double *phi, member bcmember);
    static double fval_end  (domain *d, mesh_row *row, double *phi, member bcmember);
};

}

#endif