/**
 * bcondition.h
 *
 * by Mohammad Elmi
 * Created on 25 Oct 2013
 */

#ifndef _BCONDITION_H_
#define _BCONDITION_H_

#include "common.h"
#include <jsoncpp/json/json.h>
#include "vector.h"
#include <functional>
#include "mesh.h"

namespace aban2
{
class domain;

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
    virtual ~bcondition() {}
    virtual bcdesc desc(size_t cellno, size_t dir) = 0;
    double face_val(size_t cellno, size_t dir, double *phi, double dx);

    friend class flowbc;
};

enum class bcside {start, end};

class flowbc
{
public:
    typedef bcondition *flowbc::*member;
    typedef std::function<double(domain *, mesh::row *, double *, bcside)> bc_val_getter;

    bcondition *u0, *u1, *u2, *p, *vof;

    flowbc(bcondition *_u0, bcondition *_u1, bcondition *_u2, bcondition *_p, bcondition *_vof)
        : u0(_u0), u1(_u1), u2(_u2), p(_p), vof(_vof) {}
    ~flowbc();

    static void create_bcs(Json::Value *bcroot, flowbc **boundaries, domain *_d);

    static member umembers[3];
    static bc_val_getter bc_p_getter, bc_vof_getter, bc_u_getter[3]; //, bc_rhou_getter[3];
};

}

#endif
