/**
 * bcondition.cpp
 *
 * by Mohammad Elmi
 * Created on 26 Oct 2013
 */

#include "bcondition.h"

#include "domain.h"
#include "vector.h"

namespace aban2
{

//////////////////////////  bcondition

void bcondition::create_bcs(Json::Value *bcroot, bcondition **boundaries)
{
    auto bnames = bcroot->getMemberNames();
    for (auto &x : bnames)
    {
        Json::Value bc = (*bcroot)[x];
        if (bc["type"].asString() == "velocity")
            boundaries[x[0]] = new velocitybc(&bc);
        if (bc["type"].asString() == "pressure")
            boundaries[x[0]] = new pressurebc(&bc);
    }
}

//////////////////////////  velocitybc

double velocitybc::p(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno ;
    vector dx;
    if (side == bcside::start)
    {
        cellno =  d->cellno(r, 0);
        dx.components[r->dir] = -d->delta / 2.;
    }
    else
    {
        cellno = d->cellno(r, r->n - 1);
        dx.components[r->dir] = -d->delta / 2.;
    }
    double rho = d->rho;
    return d->p[cellno] + rho * (dx * d->g);
}

double velocitybc::u(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    return value[cmpnt];
}

double velocitybc::vof(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    return d->vof[cellno];
}

velocitybc::velocitybc(Json::Value *bcdata)
{
    voftype = bctype::neumann;
    utype = bctype::dirichlet;
    ptype = bctype::neumann;

    Json::Value bcval = (*bcdata)["value"];
    value[0] = bcval[0].asDouble();
    value[1] = bcval[1].asDouble();
    value[2] = bcval[2].asDouble();
}

//////////////////////////  pressurebc

double pressurebc::p(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    return value;
}

double pressurebc::u(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    return d->u[cmpnt][cellno];
}

double pressurebc::vof(domain *d, mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    return d->vof[cellno];
}

pressurebc::pressurebc(Json::Value *bcdata)
{
    voftype = bctype::neumann;
    utype = bctype::neumann;
    ptype = bctype::dirichlet;

    Json::Value bcval = (*bcdata)["value"];
    value = bcval.asDouble();
}

}
