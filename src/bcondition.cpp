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

void bcondition::create_bcs(Json::Value *bcroot, bcondition **boundaries, domain *_d)
{
    auto bnames = bcroot->getMemberNames();
    for (auto &x : bnames)
    {
        Json::Value bc = (*bcroot)[x];
        if (bc["type"].asString() == "velocity")
            boundaries[x[0]] = new velocitybc(&bc, _d);
        if (bc["type"].asString() == "pressure")
            boundaries[x[0]] = new pressurebc(&bc, _d);
    }
}

double bcondition::face_val(bcondition::func f, double *phi, mesh_row *r, bcside side, size_t cmpnt)
{
    bcdesc desc = (this->*f)(r, side, cmpnt);
    return phi[desc.cellno] * desc.coeff + desc.val;
}

double bcondition::row_face_val(bcondition::func f, double *phi, mesh_row *r, bcside side, size_t cmpnt)
{
    bcdesc desc = (this->*f)(r, side, cmpnt);
    double val = side == bcside::start ? phi[0] : phi[r->n - 1];
    return val * desc.coeff + desc.val;
}

//////////////////////////  velocitybc

bcdesc velocitybc::p(mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno;
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
    // return d->p[cellno] + rho * (dx * d->g);
    return {cellno, 1, ((value - vector::from_data(d->ustar, cellno)) * dx) *(-d->rho / d->dt)};
}

bcdesc velocitybc::u(mesh_row *r, bcside side, size_t cmpnt)
{
    return {0, 0, value.components[cmpnt]};
}

bcdesc velocitybc::vof(mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    return {cellno, 1, 0};
}

velocitybc::velocitybc(Json::Value *bcdata, domain *_d): bcondition(_d)
{
    Json::Value bcval = (*bcdata)["value"];
    value.components[0] = bcval[0].asDouble();
    value.components[1] = bcval[1].asDouble();
    value.components[2] = bcval[2].asDouble();
}

//////////////////////////  pressurebc

bcdesc pressurebc::p(mesh_row *r, bcside side, size_t cmpnt)
{
    return {0, 0, value};
}

bcdesc pressurebc::u(mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    // return d->u[cmpnt][cellno];
    return {cellno, 1, 0};
}

bcdesc pressurebc::vof(mesh_row *r, bcside side, size_t cmpnt)
{
    size_t cellno = side == bcside::start ? d->cellno(r, 0) : d->cellno(r, r->n - 1);
    // return d->vof[cellno];
    return {cellno, 1, 0};
}

pressurebc::pressurebc(Json::Value *bcdata, domain *_d): bcondition(_d)
{
    Json::Value bcval = (*bcdata)["value"];
    value = bcval.asDouble();
}

}
