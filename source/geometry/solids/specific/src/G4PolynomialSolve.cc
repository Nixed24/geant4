#include <cmath>

#include "G4PolynomialSolve.hh"
#include "G4JTPolynomialSolver.hh"
#include "G4Types.hh"

std::vector<G4double> G4PolynomialSolve::SolvePolynomial(const std::vector<G4double>& coefficients)
{
  // begin copilot written
  std::vector<G4double> roots;

  // Check if the polynomial is of degree 0
  if (coefficients.size() == 1)
  {
    return roots; // No roots for a constant polynomial
  }

  // Check if the polynomial is of degree 1
  if (coefficients.size() == 2)
  {
    G4double root = -1 * coefficients[1] / coefficients[0];
    roots.push_back(root);
    return roots;
  }

  // Check if the polynomial is of degree 2
  if (coefficients.size() == 3)
  {
    G4double a = coefficients[0];
    G4double b = coefficients[1];
    G4double c = coefficients[2];

    G4double discriminant = b * b - 4 * a * c;

    if (discriminant > 0)
    {
      G4double root1 = (-b + std::sqrt(discriminant)) / (2 * a);
      G4double root2 = (-b - std::sqrt(discriminant)) / (2 * a); // should you really be doing this twice copilot?
      roots.push_back(root1);
      roots.push_back(root2);
    }
    else if (discriminant == 0)
    {
      G4double root = -b / (2 * a);
      roots.push_back(root);
    }
    // If discriminant < 0, no real roots

    return roots;
  }
  if (coefficients.size() == 4)
  {
    G4complex omega = std::polar(1.0, 2.0 * M_PI / 3.0);
    G4complex omega_squared = pow(omega, 2.0);
    G4double p = (3.0 * coefficients[0] * coefficients[2] - pow(coefficients[1], 2.0)) / (3.0 * pow(coefficients[0], 2.0));
    G4double q = (2.0 * pow(coefficients[1], 3) - 9 * coefficients[0] * coefficients[1] * coefficients[2] + 27 * pow(coefficients[0], 2) * coefficients[3]) / (27.0 * pow(coefficients[0], 3.0));
    G4complex s_1 = std::cbrt(-q / 2.0 + std::sqrt(pow(q / 2.0, 2.0) + pow(p / 3.0, 3.0)));
    G4complex s_2 = std::cbrt(-q / 2.0 - std::sqrt(pow(q / 2.0, 2.0) + pow(p / 3.0, 3.0)));
    G4complex x_0 = s_1 + s_2 - coefficients[1] / (3.0 * coefficients[0]);
    G4complex x_1 = omega * s_1 + omega_squared * s_2 - (coefficients[1] / (3.0 * coefficients[0]));
    G4complex x_2 = omega_squared * s_1 + omega * s_2 - (coefficients[1] / (3.0 * coefficients[0]));
    // I wrote this bit
    G4complex * pure_roots[3] = {&x_0, &x_1, &x_2};
    //
    for (auto it = std::begin(pure_roots); it != std::end(pure_roots); ++it)
    {
      if ((*it)->imag() == 0)
      {
        roots.push_back((*it)->real());
      }
    }
    return roots;
  }
  // end copilot written
  if (coefficients.size() == 5)
  {
    G4double three = 3;
    G4double two = 2;
    // ^^ Scalar multiplication between std::complex<G4double> and int is not allowed for some bizarre reason?
    G4double alpha = coefficients[1] * coefficients[1] \
    - 3 * coefficients[0] * coefficients[2] + 12 * coefficients[3];
    G4double beta = 2 * pow(coefficients[1], 3) \
    - 9 * coefficients[0] * coefficients[1] * coefficients[2]\
    + 27 * pow(coefficients[2], 2) \
    + 27 * pow(coefficients[0], 2) * coefficients[3] \
    - 72 * coefficients[1] * coefficients[3];
    G4complex gamma = sqrt(-4 * pow(alpha, 3) + (beta * beta));
    G4double delta = - 1 * pow(coefficients[0], 3) \
    + 4 * coefficients[0] * coefficients[1] \
    - 8 * coefficients[2];
    G4double eta = pow(coefficients[0], 2) - (2 * coefficients[1] / 3);
    G4complex phi = pow((beta + gamma), 1/3);
    G4complex kappa = (pow(2, 1/3) * alpha)/(three*phi);
    G4complex theta = phi/pow(54, 1/3);
    G4complex radical_1 = 0.5 * pow(eta + theta + kappa, 0.5);
    G4complex radical_2 = 0.5 * sqrt(2 * eta - theta - kappa - delta/pow((eta + theta + kappa), 0.25));
    G4complex x_0 = (-1 * coefficients[0]/4) + radical_1 + radical_2;
    G4complex x_1 = x_0 - (two * radical_1);
    G4complex x_2 = x_0 - (two * radical_2);
    G4complex x_3 = x_0 - (two * radical_1) - (two * radical_2);
    G4complex * pure_roots[4] = {&x_0, &x_1, &x_2, &x_3};
    for (auto it = std::begin(pure_roots); it != std::end(pure_roots); ++it)
    {
      if ((*it)->imag() == 0)
      {
        roots.push_back((*it)->real());
      }
    }
    return roots;
  }
  if (coefficients.size() > 5) // defer to numerical methods for polynomials of degree > 4
  {
    // begin copilot written
    G4JTPolynomialSolver solver;
    // i wrote this bit
    std::vector<G4double> coefficients_copy = coefficients;
    while (coefficients_copy.size() > 0 && coefficients_copy.front() == 0.0)
    {
      coefficients_copy.pop_back(); // Remove leading zeroes
    }
    if (coefficients_copy.size() == 0)
    {
      return roots; // No roots for a zero polynomial
    }
    // end
    G4int degree = coefficients_copy.size() - 1;
    std::vector<G4double> zeror(degree);
    std::vector<G4double> zeroi(degree);
    G4int numRoots = solver.FindRoots(const_cast<G4double*>(coefficients_copy.data()), degree, zeror.data(), zeroi.data());
    for (G4int i = 0; i < numRoots; ++i)
    {
      if (zeroi[i] == 0.0) // Only include real roots
      {
        roots.push_back(zeror[i]);
      }
    }
  }
  // end copilot written
  return roots;
}