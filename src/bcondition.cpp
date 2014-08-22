/**
 * bcondition.cpp
 *
 * by Mohammad Elmi
 * Created on 26 Oct 2013
 */

#include "bcondition.h"

#include "domain.h"
#include "vector.h"

using namespace std::placeholders;

namespace aban2
{

//////////////////////////  bcdesc

double bcdesc::face_val(double phi_P, double dx)
{
    return sw * (phi_P + dx * grad) + cte;
}

//////////////////////////  bcondition


double bcondition::face_val(size_t cellno, size_t dir, double *phi, double dx)
{
    return desc(cellno, dir).face_val(phi[cellno], dx);
}

//////////////////////////  dirichlet

class dirichlet: public bcondition
{
    double value;
public:
    virtual bcdesc desc(size_t cellno, size_t dir)
    {
        return {0, 0, value};
    }
    dirichlet(domain *_d, double _value): bcondition(_d), value(_value) {}
};

//////////////////////////  neumann

class neumann: public bcondition
{
    double value;
public:
    virtual bcdesc desc(size_t cellno, size_t dir)
    {
        return {1, value, 0};
    }
    neumann(domain *_d, double _value): bcondition(_d), value(_value) {}
};

//////////////////////////  pressure_of_const_velocity

class pressure_of_const_velocity: public bcondition
{
protected:
    vector u_bc;
public:
    virtual bcdesc desc(size_t cellno, size_t dir)
    {
        return
        {
            1,
            - (u_bc - vector::from_data(d->ustar, cellno)).cmpnt[dir] *d->rho[cellno] / d->dt,
            0
        };
    }
    pressure_of_const_velocity(domain *_d, vector _u_bc): bcondition(_d), u_bc(_u_bc) {}
};

//////////////////////////  flow_bcondition

void flowbc::create_bcs(Json::Value *bcroot, flowbc **boundaries, domain *_d)
{
    auto bnames = bcroot->getMemberNames();
    for (auto &bname : bnames)
    {
        Json::Value bc_node = (*bcroot)[bname];
        Json::Value bc_type = bc_node["type"].asString();
        Json::Value bc_val = bc_node["value"];
        flowbc *bc;

        if (bc_type == "p")
        {
            bc = new flowbc
            {
                new neumann(_d, 0),
                new neumann(_d, 0),
                new neumann(_d, 0),
                new dirichlet(_d, bc_val.asDouble()),
                new neumann(_d, 0)
            };
        }
        if (bc_type == "u")
        {
            vector u {bc_val[0].asDouble(), bc_val[1].asDouble(), bc_val[2].asDouble()};
            bc = new flowbc
            {
                new dirichlet(_d, u.x),
                new dirichlet(_d, u.y),
                new dirichlet(_d, u.z),
                new pressure_of_const_velocity(_d, u),
                new neumann(_d, 0)
            };
        }

        boundaries[bname[0]] = bc;
    }
}

flowbc::~flowbc()
{
    delete p;
    delete u0;
    delete u1;
    delete u2;
    delete vof;
}

double face_val(domain *d, mesh_row *row, double *phi, flowbc::member mem, bcside side)
{
    if (side == bcside::start)
    {
        auto bc = d->boundaries[row->start_code]->*mem;
        auto cellno = d->cellno(row, 0);
        return bc->face_val(cellno, row->dir, phi, -d->delta / 2);
    }
    else
    {
        auto bc = d->boundaries[row->end_code]->*mem;
        auto cellno = d->cellno(row, row->n - 1);
        return bc->face_val(cellno, row->dir, phi, +d->delta / 2);
    }
}

double rho_u_face_val(domain *d, mesh_row *row, double *phi, bcside side, size_t dir)
{
    double uf = face_val(d, row, d->u[dir], flowbc::umembers[dir], side);
    double vof = face_val(d, row, d->vof, &flowbc::vof, side);

    return d->rho_bar(vof) * uf;
}

flowbc::member flowbc::umembers[3] {&flowbc::u0, &flowbc::u1, &flowbc::u2};
flowbc::bc_val_getter flowbc::bc_p_getter = std::bind(face_val, _1, _2, _3, &flowbc::p, _4);
flowbc::bc_val_getter flowbc::bc_vof_getter = std::bind(face_val, _1, _2, _3, &flowbc::vof, _4);
flowbc::bc_val_getter flowbc::bc_u_getter[3] =
{
    std::bind(face_val, _1, _2, _3, &flowbc::u0, _4),
    std::bind(face_val, _1, _2, _3, &flowbc::u1, _4),
    std::bind(face_val, _1, _2, _3, &flowbc::u2, _4)
};
flowbc::bc_val_getter flowbc::bc_rhou_getter[3] =
{
    std::bind(rho_u_face_val, _1, _2, _3, _4, 0),
    std::bind(rho_u_face_val, _1, _2, _3, _4, 1),
    std::bind(rho_u_face_val, _1, _2, _3, _4, 2)
};

}