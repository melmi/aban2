#include <iterator>

#include "lpw.h"

extern "C"
{
#include <laspack/rtc.h>
#include <laspack/vector.h>
#include <laspack/qmatrix.h>
#include <laspack/itersolv.h>
#include <laspack/precond.h>
}

namespace lpw
{
matrix_t::matrix_t(int _n): n(_n) {}

int matrix_t::idx(int i, int j)
{
    return j * n + i;
}

double& matrix_t::operator()(int i, int j)
{
    int ix = idx(i, j);
    auto ret = elements.insert(std::make_pair(ix, 0.0));
    return ret.first->second;
}

// structure of the algorith borrowed from here: http://en.cppreference.com/w/cpp/algorithm/count
template<class InputIt, class UnaryPredicate>
int count_while(InputIt first, InputIt last, UnaryPredicate p)
{
    int ret = 0;
    for (; first != last && p(*first); ++first)
    {
        if (p(*first))
        {
            ret++;
        }
    }
    return ret;
}

void matrix_t::construct_qmatrix (QMatrix* q, std::string s)
{
    Q_Constr(q, const_cast<char*>(s.c_str()), n, False, Rowws, Normal, True);

    auto iter = elements.begin();
    while (iter != elements.end())
    {
        int r = iter->first / n;
        int max_idx = idx(0,r + 1);
        int nr = count_while(iter, elements.end(), [max_idx](auto i)
        {
            return i.first < max_idx;
        });
        Q_SetLen(q, r+1, nr);
        for (int i = 0; i < nr; ++i)
        {
            Q_SetEntry(q, r + 1, i, iter->first % n + 1, iter->second);
            ++iter;
        }
    }
}

int matrix_t::get_n()
{
    return n;
}

void construct_vector(Vector* v, double*arr, int n, std::string s)
{
    V_Constr(v, const_cast<char*>(s.c_str()), n, Normal, True);
    for (int i = 0; i < n; ++i)
        V_SetCmp(v, i + 1, arr[i]);
}

void vector2arr(Vector* v, double* arr)
{
    int n=V_GetDim(v);
    for(int i=0;i<n;++i)
        arr[i]=V_GetCmp(v,i+1);
}

PrecondProcType get_precond_proc(precond_t p)
{
    switch (p)
    {
    case precond_t::ilu:
        return ILUPrecond;
    case precond_t::ssor:
        return SSORPrecond;
    case precond_t::jacobi:
        return JacobiPrecond;
    default:
        return NULL;
    }
}

int cg(matrix_t*a, double*x, double*b, int max_iter, precond_t p, double omega)
{
    QMatrix A;
    Vector xx, bb;

    a->construct_qmatrix(&A, "A");
    construct_vector(&xx, x, a->get_n(), "x");
    construct_vector(&bb, b, a->get_n(), "b");
    auto precond=get_precond_proc(p);
    CGIter(&A, &xx, &bb, max_iter, precond, omega);
    vector2arr(&xx, x);
    Q_Destr(&A);
    V_Destr(&xx);
    V_Destr(&bb);

    return GetLastNoIter();
}

}
