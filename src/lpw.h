#include <map>
#include <string>

struct QMatrixType;
typedef QMatrixType QMatrix;
namespace lpw
{

class matrix_t
{
    const double epsilon = 1e-10;
    int n;
    int idx(int, int);
    std::map<int, double> elements;

public:
    matrix_t(int);
    void construct_qmatrix(QMatrix*, std::string);
    double& operator()(int, int);
    int get_n();
};

enum class precond_t
{
    jacobi, ssor, ilu, none
};

int cg(matrix_t*a, double*x, double*b, int max_iter, precond_t p, double omega);

}
