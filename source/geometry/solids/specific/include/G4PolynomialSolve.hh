#ifndef G4POLYNOMIAL_SOLVE_HH
#define G4POLYNOMIAL_SOLVE_HH

#include <vector>

class G4PolynomialSolve
{
  public:

    G4PolynomialSolve() = default;
    ~G4PolynomialSolve() = default;

    static std::vector<double> SolvePolynomial(const std::vector<double>& coefficients);
};

#endif
