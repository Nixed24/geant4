#include <algorithm>
#include <cmath>
#include "G4HalfQuadric.hh"
#include "G4PolynomialSolve.hh"

G4RotationMatrix G4HalfQuadric::_T(const G4RotationMatrix & M)
{
  G4RotationMatrix M_T;
  M_T.set(M.rowX(), M.rowY(), M.rowZ()); // set columns to rows
  return M_T;
}
std::vector<G4ThreeVector> G4HalfQuadric::GetLMContourIntersections(const G4ThreeVector& p) const
{
  const G4PolynomialSolve polysolve;
  const G4double r_x = p.x(), r_y = p.y(), r_z = p.z();
  const G4double M_0 = this->R;
  const G4double b_x = this->P.x() * 0.5, b_y = this->P.y() * 0.5, b_z = this->P.z() * 0.5;
  const G4double a = Q.xx(), b = Q.xy(), c = Q.xz(), e = Q.yy(), f = Q.yz(), i = Q.zz();
  // Edge cases: Determinant = 0 | Determinant is a cubic polynomial in lambda, could add a check for lambda being at one of its roots
  // I don't know what the determinant being zero would actually correspond to, it's not a plane because that would just give identity
  // but I think it could mean that there are infinitely many solutions for the closest point OR maybe the p is already on the surface
  // x of lambda will also div by zero so we'll know
  // Z -> lambda (lagrange multiplier)
  // const G4double characteristic_determinant = (-pow(Z,3)*a*e*i)+pow(Z,2)*e*i+pow(Z,3)*pow(b,2)*i+pow(Z,2)*a*i-Z*i+pow(Z,3)*a*pow(f,2)-pow(Z,2)*pow(f,2)-2*pow(Z,3)*b*c*f+pow(Z,3)*pow(c,2)*e+pow(Z,2)*a*e-Z*e-pow(Z,2)*pow(c,2)-pow(Z,2)*pow(b,2)-Z*a+1;
  // It's a LONG expression
  // (DONE) 07/11/26 these terms aren't correct; pull the correct ones from poly_redo_test_demo; that also uses b (P/2) instead of P
  const G4double power_6 =
(-a*pow(b_x,2)*pow(e,2)*pow(i,2))+M_0*pow(a,2)*pow(e,2)*pow(i,2)-pow(a,2)*pow(b_y,2)*e*pow(i,2)+2*a*b*b_x*b_y*e*pow(i,2)
                         +pow(b,2)*pow(b_x,2)*e*pow(i,2)-2*M_0*a*pow(b,2)*e*pow(i,2)+a*pow(b,2)*pow(b_y,2)*pow(i,2)
                         -2*pow(b,3)*b_x*b_y*pow(i,2)+M_0*pow(b,4)*pow(i,2)+2*a*pow(b_x,2)*e*pow(f,2)*i
                         -2*M_0*pow(a,2)*e*pow(f,2)*i+pow(a,2)*pow(b_y,2)*pow(f,2)*i
                         -2*a*b*b_x*b_y*pow(f,2)*i-pow(b,2)*pow(b_x,2)*pow(f,2)*i
                         +2*M_0*a*pow(b,2)*pow(f,2)*i-2*a*b_x*b_y*c*e*f*i
                         -2*b*pow(b_x,2)*c*e*f*i+4*M_0*a*b*c*e*f*i
                         +2*pow(a,2)*b_y*b_z*e*f*i-2*a*b*b_x*b_z*e*f*i
                         -2*a*b*pow(b_y,2)*c*f*i+6*pow(b,2)*b_x*b_y*c*f*i
                         -4*M_0*pow(b,3)*c*f*i-2*a*pow(b,2)*b_y*b_z*f*i
                         +2*pow(b,3)*b_x*b_z*f*i+pow(b_x,2)*pow(c,2)*pow(e,2)*i-2*M_0*a*pow(c,2)*pow(e,2)*i
                         +2*a*b_x*b_z*c*pow(e,2)*i-pow(a,2)*pow(b_z,2)*pow(e,2)*i
                         +2*a*pow(b_y,2)*pow(c,2)*e*i-2*b*b_x*b_y*pow(c,2)*e*i
                         +2*M_0*pow(b,2)*pow(c,2)*e*i-2*a*b*b_y*b_z*c*e*i
                         -2*pow(b,2)*b_x*b_z*c*e*i+2*a*pow(b,2)*pow(b_z,2)*e*i
                         -pow(b,2)*pow(b_y,2)*pow(c,2)*i+2*pow(b,3)*b_y*b_z*c*i-pow(b,4)*pow(b_z,2)*i
                         -a*pow(b_x,2)*pow(f,4)+M_0*pow(a,2)*pow(f,4)+2*a*b_x*b_y*c*pow(f,3)
                         +2*b*pow(b_x,2)*c*pow(f,3)-4*M_0*a*b*c*pow(f,3)-2*pow(a,2)*b_y*b_z*pow(f,3)
                         +2*a*b*b_x*b_z*pow(f,3)-pow(b_x,2)*pow(c,2)*e*pow(f,2)+2*M_0*a*pow(c,2)*e*pow(f,2)
                         -2*a*b_x*b_z*c*e*pow(f,2)+pow(a,2)*pow(b_z,2)*e*pow(f,2)-a*pow(b_y,2)*pow(c,2)*pow(f,2)
                         -4*b*b_x*b_y*pow(c,2)*pow(f,2)+4*M_0*pow(b,2)*pow(c,2)*pow(f,2)
                         +6*a*b*b_y*b_z*c*pow(f,2)-4*pow(b,2)*b_x*b_z*c*pow(f,2)
                         -a*pow(b,2)*pow(b_z,2)*pow(f,2)+2*b_x*b_y*pow(c,3)*e*f-4*M_0*b*pow(c,3)*e*f
                         -2*a*b_y*b_z*pow(c,2)*e*f+6*b*b_x*b_z*pow(c,2)*e*f
                         -2*a*b*pow(b_z,2)*c*e*f+2*b*pow(b_y,2)*pow(c,3)*f
                         -4*pow(b,2)*b_y*b_z*pow(c,2)*f+2*pow(b,3)*pow(b_z,2)*c*f+M_0*pow(c,4)*pow(e,2)
                         -2*b_x*b_z*pow(c,3)*pow(e,2)+a*pow(b_z,2)*pow(c,2)*pow(e,2)-pow(b_y,2)*pow(c,4)*e
                         +2*b*b_y*b_z*pow(c,3)*e-pow(b,2)*pow(b_z,2)*pow(c,2)*e;
  const G4double power_5 =
2*pow(b_x,2)*pow(e,2)*pow(i,2)-2*M_0*a*pow(e,2)*pow(i,2)+2*a*pow(b_y,2)*e*pow(i,2)-4*b*b_x*b_y*e*pow(i,2)
                      +2*a*pow(b_x,2)*e*pow(i,2)+2*M_0*pow(b,2)*e*pow(i,2)-2*M_0*pow(a,2)*e*pow(i,2)
                      +2*pow(a,2)*pow(b_y,2)*pow(i,2)-4*a*b*b_x*b_y*pow(i,2)+2*M_0*a*pow(b,2)*pow(i,2)
                      -4*pow(b_x,2)*e*pow(f,2)*i+4*M_0*a*e*pow(f,2)*i-2*a*pow(b_y,2)*pow(f,2)*i
                      +4*b*b_x*b_y*pow(f,2)*i-2*a*pow(b_x,2)*pow(f,2)*i-2*M_0*pow(b,2)*pow(f,2)*i
                      +2*M_0*pow(a,2)*pow(f,2)*i+4*b_x*b_y*c*e*f*i-4*M_0*b*c*e*f*i
                      -4*a*b_y*b_z*e*f*i+4*b*b_x*b_z*e*f*i+4*a*b_x*b_y*c*f*i
                      -4*M_0*a*b*c*f*i-4*pow(a,2)*b_y*b_z*f*i+4*a*b*b_x*b_z*f*i
                      +2*M_0*pow(c,2)*pow(e,2)*i-4*b_x*b_z*c*pow(e,2)*i+2*a*pow(b_z,2)*pow(e,2)*i
                      +2*a*pow(b_x,2)*pow(e,2)*i-2*M_0*pow(a,2)*pow(e,2)*i-2*pow(b_y,2)*pow(c,2)*e*i
                      -2*pow(b_x,2)*pow(c,2)*e*i+4*M_0*a*pow(c,2)*e*i+4*b*b_y*b_z*c*e*i
                      -4*a*b_x*b_z*c*e*i-2*pow(b,2)*pow(b_z,2)*e*i+2*pow(a,2)*pow(b_z,2)*e*i
                      +2*pow(a,2)*pow(b_y,2)*e*i-4*a*b*b_x*b_y*e*i-2*pow(b,2)*pow(b_x,2)*e*i
                      +4*M_0*a*pow(b,2)*e*i-4*a*pow(b_y,2)*pow(c,2)*i+4*b*b_x*b_y*pow(c,2)*i
                      -2*M_0*pow(b,2)*pow(c,2)*i+4*a*b*b_y*b_z*c*i-2*a*pow(b,2)*pow(b_z,2)*i
                      -2*a*pow(b,2)*pow(b_y,2)*i+4*pow(b,3)*b_x*b_y*i-2*M_0*pow(b,4)*i+2*pow(b_x,2)*pow(f,4)
                      -2*M_0*a*pow(f,4)-4*b_x*b_y*c*pow(f,3)+4*M_0*b*c*pow(f,3)
                      +4*a*b_y*b_z*pow(f,3)-4*b*b_x*b_z*pow(f,3)-2*M_0*pow(c,2)*e*pow(f,2)
                      +4*b_x*b_z*c*e*pow(f,2)-2*a*pow(b_z,2)*e*pow(f,2)-2*a*pow(b_x,2)*e*pow(f,2)
                      +2*M_0*pow(a,2)*e*pow(f,2)+2*pow(b_y,2)*pow(c,2)*pow(f,2)+2*pow(b_x,2)*pow(c,2)*pow(f,2)
                      -2*M_0*a*pow(c,2)*pow(f,2)-4*b*b_y*b_z*c*pow(f,2)+2*pow(b,2)*pow(b_z,2)*pow(f,2)
                      +2*pow(b,2)*pow(b_x,2)*pow(f,2)-2*M_0*a*pow(b,2)*pow(f,2)+4*a*b_x*b_y*c*e*f
                      -4*M_0*a*b*c*e*f-4*pow(a,2)*b_y*b_z*e*f+4*a*b*b_x*b_z*e*f
                      -4*b_x*b_y*pow(c,3)*f+4*M_0*b*pow(c,3)*f+4*a*b_y*b_z*pow(c,2)*f
                      -4*b*b_x*b_z*pow(c,2)*f-4*pow(b,2)*b_x*b_y*c*f+4*M_0*pow(b,3)*c*f
                      +4*a*pow(b,2)*b_y*b_z*f-4*pow(b,3)*b_x*b_z*f+2*M_0*a*pow(c,2)*pow(e,2)
                      -4*a*b_x*b_z*c*pow(e,2)+2*pow(a,2)*pow(b_z,2)*pow(e,2)-2*M_0*pow(c,4)*e
                      +4*b_x*b_z*pow(c,3)*e-2*a*pow(b_z,2)*pow(c,2)*e-2*a*pow(b_y,2)*pow(c,2)*e
                      -2*M_0*pow(b,2)*pow(c,2)*e+4*a*b*b_y*b_z*c*e+4*pow(b,2)*b_x*b_z*c*e
                      -4*a*pow(b,2)*pow(b_z,2)*e+2*pow(b_y,2)*pow(c,4)-4*b*b_y*b_z*pow(c,3)
                      +2*pow(b,2)*pow(b_z,2)*pow(c,2)+2*pow(b,2)*pow(b_y,2)*pow(c,2)-4*pow(b,3)*b_y*b_z*c
                      +2*pow(b,4)*pow(b_z,2);
  const G4double power_4 =
pow(a,2)*pow(e,2)*i*pow(r_z,2)-2*a*pow(b,2)*e*i*pow(r_z,2)+pow(b,4)*i*pow(r_z,2)-pow(a,2)*e*pow(f,2)*pow(r_z,2)
                      +a*pow(b,2)*pow(f,2)*pow(r_z,2)+2*a*b*c*e*f*pow(r_z,2)-2*pow(b,3)*c*f*pow(r_z,2)
                      -a*pow(c,2)*pow(e,2)*pow(r_z,2)+pow(b,2)*pow(c,2)*e*pow(r_z,2)-2*pow(a,2)*e*f*i*r_y*r_z
                      +2*a*pow(b,2)*f*i*r_y*r_z+2*a*b*c*e*i*r_y*r_z
                      -2*pow(b,3)*c*i*r_y*r_z+2*pow(a,2)*pow(f,3)*r_y*r_z-6*a*b*c*pow(f,2)*r_y*r_z
                      +2*a*pow(c,2)*e*f*r_y*r_z+4*pow(b,2)*pow(c,2)*f*r_y*r_z
                      -2*b*pow(c,3)*e*r_y*r_z+2*a*b*e*f*i*r_x*r_z-2*pow(b,3)*f*i*r_x*r_z
                      -2*a*c*pow(e,2)*i*r_x*r_z+2*pow(b,2)*c*e*i*r_x*r_z
                      -2*a*b*pow(f,3)*r_x*r_z+2*a*c*e*pow(f,2)*r_x*r_z
                      +4*pow(b,2)*c*pow(f,2)*r_x*r_z-6*b*pow(c,2)*e*f*r_x*r_z
                      +2*pow(c,3)*pow(e,2)*r_x*r_z+2*b*b_x*e*f*i*r_z-2*pow(b,2)*b_y*f*i*r_z
                      -2*pow(a,2)*b_y*f*i*r_z+2*a*b*b_x*f*i*r_z-2*b_x*c*pow(e,2)*i*r_z
                      +2*b*b_y*c*e*i*r_z+2*a*b*b_y*c*i*r_z-2*pow(b,2)*b_x*c*i*r_z
                      -2*b*b_x*pow(f,3)*r_z+2*b_x*c*e*pow(f,2)*r_z+2*b*b_y*c*pow(f,2)*r_z
                      -2*a*b_x*c*pow(f,2)*r_z+2*pow(b,2)*b_z*pow(f,2)*r_z+2*pow(a,2)*b_z*pow(f,2)*r_z
                      -2*b_y*pow(c,2)*e*f*r_z-4*b*b_z*c*e*f*r_z-2*pow(a,2)*b_y*e*f*r_z
                      +2*a*b*b_x*e*f*r_z+2*a*b_y*pow(c,2)*f*r_z+2*b*b_x*pow(c,2)*f*r_z
                      -4*a*b*b_z*c*f*r_z+2*a*pow(b,2)*b_y*f*r_z-2*pow(b,3)*b_x*f*r_z
                      +2*b_z*pow(c,2)*pow(e,2)*r_z-2*a*b_x*c*pow(e,2)*r_z+2*pow(a,2)*b_z*pow(e,2)*r_z
                      +2*a*b*b_y*c*e*r_z+2*pow(b,2)*b_x*c*e*r_z-4*a*pow(b,2)*b_z*e*r_z
                      -2*b*b_y*pow(c,3)*r_z+2*pow(b,2)*b_z*pow(c,2)*r_z-2*pow(b,3)*b_y*c*r_z
                      +2*pow(b,4)*b_z*r_z+pow(a,2)*e*pow(i,2)*pow(r_y,2)-a*pow(b,2)*pow(i,2)*pow(r_y,2)
                      -pow(a,2)*pow(f,2)*i*pow(r_y,2)+2*a*b*c*f*i*pow(r_y,2)-2*a*pow(c,2)*e*i*pow(r_y,2)
                      +pow(b,2)*pow(c,2)*i*pow(r_y,2)+a*pow(c,2)*pow(f,2)*pow(r_y,2)-2*b*pow(c,3)*f*pow(r_y,2)
                      +pow(c,4)*e*pow(r_y,2)-2*a*b*e*pow(i,2)*r_x*r_y+2*pow(b,3)*pow(i,2)*r_x*r_y
                      +2*a*b*pow(f,2)*i*r_x*r_y+2*a*c*e*f*i*r_x*r_y
                      -6*pow(b,2)*c*f*i*r_x*r_y+2*b*pow(c,2)*e*i*r_x*r_y
                      -2*a*c*pow(f,3)*r_x*r_y+4*b*pow(c,2)*pow(f,2)*r_x*r_y-2*pow(c,3)*e*f*r_x*r_y
                      -2*b*b_x*e*pow(i,2)*r_y+2*pow(b,2)*b_y*pow(i,2)*r_y+2*pow(a,2)*b_y*pow(i,2)*r_y
                      -2*a*b*b_x*pow(i,2)*r_y+2*b*b_x*pow(f,2)*i*r_y+2*b_x*c*e*f*i*r_y
                      -4*b*b_y*c*f*i*r_y+2*a*b_x*c*f*i*r_y-2*pow(b,2)*b_z*f*i*r_y
                      -2*pow(a,2)*b_z*f*i*r_y+2*b*b_z*c*e*i*r_y-4*a*b_y*pow(c,2)*i*r_y
                      +2*b*b_x*pow(c,2)*i*r_y+2*a*b*b_z*c*i*r_y-2*b_x*c*pow(f,3)*r_y
                      +2*b_y*pow(c,2)*pow(f,2)*r_y+2*b*b_z*c*pow(f,2)*r_y+2*pow(a,2)*b_y*pow(f,2)*r_y
                      -2*a*b*b_x*pow(f,2)*r_y-2*b_z*pow(c,2)*e*f*r_y+2*a*b_x*c*e*f*r_y
                      -2*pow(a,2)*b_z*e*f*r_y-2*b_x*pow(c,3)*f*r_y+2*a*b_z*pow(c,2)*f*r_y
                      -4*a*b*b_y*c*f*r_y+2*pow(b,2)*b_x*c*f*r_y+2*a*pow(b,2)*b_z*f*r_y
                      -2*b*b_x*pow(c,2)*e*r_y+2*a*b*b_z*c*e*r_y+2*b_y*pow(c,4)*r_y
                      -2*b*b_z*pow(c,3)*r_y+2*pow(b,2)*b_y*pow(c,2)*r_y-2*pow(b,3)*b_z*c*r_y
                      +a*pow(e,2)*pow(i,2)*pow(r_x,2)-pow(b,2)*e*pow(i,2)*pow(r_x,2)-2*a*e*pow(f,2)*i*pow(r_x,2)
                      +pow(b,2)*pow(f,2)*i*pow(r_x,2)+2*b*c*e*f*i*pow(r_x,2)-pow(c,2)*pow(e,2)*i*pow(r_x,2)
                      +a*pow(f,4)*pow(r_x,2)-2*b*c*pow(f,3)*pow(r_x,2)+pow(c,2)*e*pow(f,2)*pow(r_x,2)
                      +2*b_x*pow(e,2)*pow(i,2)*r_x-2*b*b_y*e*pow(i,2)*r_x-2*a*b*b_y*pow(i,2)*r_x
                      +2*pow(b,2)*b_x*pow(i,2)*r_x-4*b_x*e*pow(f,2)*i*r_x+2*b*b_y*pow(f,2)*i*r_x
                      +2*b_y*c*e*f*i*r_x+2*b*b_z*e*f*i*r_x+2*a*b_y*c*f*i*r_x
                      -4*b*b_x*c*f*i*r_x+2*a*b*b_z*f*i*r_x-2*b_z*c*pow(e,2)*i*r_x
                      +2*b*b_y*pow(c,2)*i*r_x-2*pow(b,2)*b_z*c*i*r_x+2*b_x*pow(f,4)*r_x
                      -2*b_y*c*pow(f,3)*r_x-2*b*b_z*pow(f,3)*r_x+2*b_z*c*e*pow(f,2)*r_x
                      +2*b_x*pow(c,2)*pow(f,2)*r_x-2*a*b_z*c*pow(f,2)*r_x-2*a*b*b_y*pow(f,2)*r_x
                      +2*pow(b,2)*b_x*pow(f,2)*r_x+2*a*b_y*c*e*f*r_x-4*b*b_x*c*e*f*r_x
                      +2*a*b*b_z*e*f*r_x-2*b_y*pow(c,3)*f*r_x+2*b*b_z*pow(c,2)*f*r_x
                      +2*pow(b,2)*b_y*c*f*r_x-2*pow(b,3)*b_z*f*r_x+2*b_x*pow(c,2)*pow(e,2)*r_x
                      -2*a*b_z*c*pow(e,2)*r_x-2*b*b_y*pow(c,2)*e*r_x+2*pow(b,2)*b_z*c*e*r_x
                      +M_0*pow(e,2)*pow(i,2)-pow(b_y,2)*e*pow(i,2)-4*pow(b_x,2)*e*pow(i,2)+4*M_0*a*e*pow(i,2)
                      -4*a*pow(b_y,2)*pow(i,2)+6*b*b_x*b_y*pow(i,2)-a*pow(b_x,2)*pow(i,2)-2*M_0*pow(b,2)*pow(i,2)
                      +M_0*pow(a,2)*pow(i,2)-2*M_0*e*pow(f,2)*i+pow(b_y,2)*pow(f,2)*i+4*pow(b_x,2)*pow(f,2)*i
                      -4*M_0*a*pow(f,2)*i+2*b_y*b_z*e*f*i-6*b_x*b_y*c*f*i
                      +4*M_0*b*c*f*i+8*a*b_y*b_z*f*i-6*b*b_x*b_z*f*i
                      -pow(b_z,2)*pow(e,2)*i-4*pow(b_x,2)*pow(e,2)*i+4*M_0*a*pow(e,2)*i-4*M_0*pow(c,2)*e*i
                      +8*b_x*b_z*c*e*i-4*a*pow(b_z,2)*e*i-4*a*pow(b_y,2)*e*i
                      +8*b*b_x*b_y*e*i-4*a*pow(b_x,2)*e*i-4*M_0*pow(b,2)*e*i
                      +4*M_0*pow(a,2)*e*i+4*pow(b_y,2)*pow(c,2)*i+pow(b_x,2)*pow(c,2)*i-2*M_0*a*pow(c,2)*i
                      -6*b*b_y*b_z*c*i+2*a*b_x*b_z*c*i+2*pow(b,2)*pow(b_z,2)*i
                      -pow(a,2)*pow(b_z,2)*i-4*pow(a,2)*pow(b_y,2)*i+8*a*b*b_x*b_y*i-4*M_0*a*pow(b,2)*i
                      +M_0*pow(f,4)-2*b_y*b_z*pow(f,3)+pow(b_z,2)*e*pow(f,2)+4*pow(b_x,2)*e*pow(f,2)
                      -4*M_0*a*e*pow(f,2)+2*M_0*pow(c,2)*pow(f,2)-2*b_x*b_z*c*pow(f,2)
                      -2*b*b_x*b_y*pow(f,2)+2*a*pow(b_x,2)*pow(f,2)+2*M_0*pow(b,2)*pow(f,2)
                      -2*M_0*pow(a,2)*pow(f,2)-6*b_x*b_y*c*e*f+4*M_0*b*c*e*f
                      +8*a*b_y*b_z*e*f-6*b*b_x*b_z*e*f-2*b_y*b_z*pow(c,2)*f
                      +2*b*pow(b_z,2)*c*f+2*b*pow(b_y,2)*c*f-6*a*b_x*b_y*c*f
                      +2*b*pow(b_x,2)*c*f+4*M_0*a*b*c*f-2*pow(b,2)*b_y*b_z*f
                      +6*pow(a,2)*b_y*b_z*f-6*a*b*b_x*b_z*f-2*M_0*pow(c,2)*pow(e,2)
                      +6*b_x*b_z*c*pow(e,2)-4*a*pow(b_z,2)*pow(e,2)-a*pow(b_x,2)*pow(e,2)+M_0*pow(a,2)*pow(e,2)
                      +2*pow(b_y,2)*pow(c,2)*e-4*M_0*a*pow(c,2)*e-6*b*b_y*b_z*c*e
                      +8*a*b_x*b_z*c*e+4*pow(b,2)*pow(b_z,2)*e-4*pow(a,2)*pow(b_z,2)*e-pow(a,2)*pow(b_y,2)*e
                      +2*a*b*b_x*b_y*e+pow(b,2)*pow(b_x,2)*e-2*M_0*a*pow(b,2)*e+M_0*pow(c,4)
                      -2*b_x*b_z*pow(c,3)+a*pow(b_z,2)*pow(c,2)+4*a*pow(b_y,2)*pow(c,2)-2*b*b_x*b_y*pow(c,2)
                      +2*M_0*pow(b,2)*pow(c,2)-6*a*b*b_y*b_z*c-2*pow(b,2)*b_x*b_z*c
                      +4*a*pow(b,2)*pow(b_z,2)+a*pow(b,2)*pow(b_y,2)-2*pow(b,3)*b_x*b_y+M_0*pow(b,4);
  const G4double power_3 =
    (-2*a*pow(e,2)*i*pow(r_z,2))+2*pow(b,2)*e*i*pow(r_z,2)-2*pow(a,2)*e*i*pow(r_z,2)+2*a*pow(b,2)*i*pow(r_z,2)
                         +2*a*e*pow(f,2)*pow(r_z,2)+2*pow(a,2)*pow(f,2)*pow(r_z,2)-4*b*c*e*f*pow(r_z,2)
                         -4*a*b*c*f*pow(r_z,2)+2*pow(c,2)*pow(e,2)*pow(r_z,2)+2*a*pow(c,2)*e*pow(r_z,2)
                         +4*a*e*f*i*r_y*r_z-4*pow(b,2)*f*i*r_y*r_z-4*a*pow(f,3)*r_y*r_z
                         +8*b*c*pow(f,2)*r_y*r_z-4*pow(c,2)*e*f*r_y*r_z
                         +4*a*c*e*i*r_x*r_z-4*pow(b,2)*c*i*r_x*r_z
                         -4*a*c*pow(f,2)*r_x*r_z+8*b*pow(c,2)*f*r_x*r_z-4*pow(c,3)*e*r_x*r_z
                         +4*a*b_y*f*i*r_z-4*b*b_x*f*i*r_z+4*b_x*c*e*i*r_z
                         -4*b*b_y*c*i*r_z-4*a*b_z*pow(f,2)*r_z+4*a*b_y*e*f*r_z
                         -4*b*b_x*e*f*r_z+8*b*b_z*c*f*r_z+4*pow(a,2)*b_y*f*r_z
                         -4*a*b*b_x*f*r_z+4*b_x*c*pow(e,2)*r_z-4*a*b_z*pow(e,2)*r_z
                         -4*b_z*pow(c,2)*e*r_z-4*b*b_y*c*e*r_z+4*a*b_x*c*e*r_z
                         +4*pow(b,2)*b_z*e*r_z-4*pow(a,2)*b_z*e*r_z-4*a*b*b_y*c*r_z
                         +4*a*pow(b,2)*b_z*r_z-2*a*e*pow(i,2)*pow(r_y,2)+2*pow(b,2)*pow(i,2)*pow(r_y,2)
                         +2*a*pow(f,2)*i*pow(r_y,2)-4*b*c*f*i*pow(r_y,2)+2*pow(c,2)*e*i*pow(r_y,2)
                         -2*pow(a,2)*e*i*pow(r_y,2)+2*a*pow(b,2)*i*pow(r_y,2)+2*pow(a,2)*pow(f,2)*pow(r_y,2)
                         -4*a*b*c*f*pow(r_y,2)+2*a*pow(c,2)*e*pow(r_y,2)+4*a*b*e*i*r_x*r_y
                         -4*pow(b,3)*i*r_x*r_y-4*a*b*pow(f,2)*r_x*r_y+8*pow(b,2)*c*f*r_x*r_y
                         -4*b*pow(c,2)*e*r_x*r_y-4*a*b_y*pow(i,2)*r_y+4*b*b_x*pow(i,2)*r_y
                         -4*b_x*c*f*i*r_y+4*a*b_z*f*i*r_y+4*b*b_x*e*i*r_y
                         +4*b_y*pow(c,2)*i*r_y-4*b*b_z*c*i*r_y-4*pow(b,2)*b_y*i*r_y
                         -4*pow(a,2)*b_y*i*r_y+4*a*b*b_x*i*r_y-4*a*b_y*pow(f,2)*r_y
                         -4*b_x*c*e*f*r_y+4*a*b_z*e*f*r_y+8*b*b_y*c*f*r_y
                         -4*a*b_x*c*f*r_y+4*pow(a,2)*b_z*f*r_y-4*b*b_z*c*e*r_y
                         +4*a*b_y*pow(c,2)*r_y-4*a*b*b_z*c*r_y-2*a*e*pow(i,2)*pow(r_x,2)
                         +2*pow(b,2)*pow(i,2)*pow(r_x,2)+2*a*pow(f,2)*i*pow(r_x,2)-4*b*c*f*i*pow(r_x,2)
                         -2*a*pow(e,2)*i*pow(r_x,2)+2*pow(c,2)*e*i*pow(r_x,2)+2*pow(b,2)*e*i*pow(r_x,2)
                         +2*a*e*pow(f,2)*pow(r_x,2)-4*b*c*e*f*pow(r_x,2)+2*pow(c,2)*pow(e,2)*pow(r_x,2)
                         -4*b_x*e*pow(i,2)*r_x+4*b*b_y*pow(i,2)*r_x+4*b_x*pow(f,2)*i*r_x
                         -4*b_y*c*f*i*r_x-4*b*b_z*f*i*r_x-4*b_x*pow(e,2)*i*r_x
                         +4*b_z*c*e*i*r_x+4*b*b_y*e*i*r_x+4*a*b*b_y*i*r_x
                         -4*pow(b,2)*b_x*i*r_x+4*b_x*e*pow(f,2)*r_x-4*b_y*c*e*f*r_x
                         -4*b*b_z*e*f*r_x-4*a*b_y*c*f*r_x+8*b*b_x*c*f*r_x
                         -4*a*b*b_z*f*r_x+4*b_z*c*pow(e,2)*r_x-4*b_x*pow(c,2)*e*r_x
                         +4*a*b_z*c*e*r_x-2*M_0*e*pow(i,2)+2*pow(b_y,2)*pow(i,2)+2*pow(b_x,2)*pow(i,2)
                         -2*M_0*a*pow(i,2)+2*M_0*pow(f,2)*i-4*b_y*b_z*f*i-2*M_0*pow(e,2)*i
                         +2*pow(b_z,2)*e*i+2*pow(b_y,2)*e*i+8*pow(b_x,2)*e*i-8*M_0*a*e*i
                         +2*M_0*pow(c,2)*i-4*b_x*b_z*c*i+2*a*pow(b_z,2)*i+8*a*pow(b_y,2)*i
                         -12*b*b_x*b_y*i+2*a*pow(b_x,2)*i+4*M_0*pow(b,2)*i-2*M_0*pow(a,2)*i
                         +2*M_0*e*pow(f,2)-4*pow(b_x,2)*pow(f,2)+4*M_0*a*pow(f,2)-4*b_y*b_z*e*f
                         +8*b_x*b_y*c*f-4*M_0*b*c*f-12*a*b_y*b_z*f
                         +8*b*b_x*b_z*f+2*pow(b_z,2)*pow(e,2)+2*pow(b_x,2)*pow(e,2)-2*M_0*a*pow(e,2)
                         +4*M_0*pow(c,2)*e-12*b_x*b_z*c*e+8*a*pow(b_z,2)*e+2*a*pow(b_y,2)*e
                         -4*b*b_x*b_y*e+2*a*pow(b_x,2)*e+2*M_0*pow(b,2)*e-2*M_0*pow(a,2)*e
                         -4*pow(b_y,2)*pow(c,2)+2*M_0*a*pow(c,2)+8*b*b_y*b_z*c-4*a*b_x*b_z*c
                         -4*pow(b,2)*pow(b_z,2)+2*pow(a,2)*pow(b_z,2)+2*pow(a,2)*pow(b_y,2)-4*a*b*b_x*b_y+2*M_0*a*pow(b,2);
  const G4double power_2 =
    pow(e,2)*i*pow(r_z,2)+4*a*e*i*pow(r_z,2)-2*pow(b,2)*i*pow(r_z,2)+pow(a,2)*i*pow(r_z,2)-e*pow(f,2)*pow(r_z,2)
                  -4*a*pow(f,2)*pow(r_z,2)+6*b*c*f*pow(r_z,2)-4*pow(c,2)*e*pow(r_z,2)-a*pow(c,2)*pow(r_z,2)
                  -2*e*f*i*r_y*r_z-2*b*c*i*r_y*r_z+2*pow(f,3)*r_y*r_z
                  +2*pow(c,2)*f*r_y*r_z+2*pow(b,2)*f*r_y*r_z+2*pow(a,2)*f*r_y*r_z
                  -2*b*c*e*r_y*r_z-2*a*b*c*r_y*r_z-2*b*f*i*r_x*r_z
                  -2*a*c*i*r_x*r_z+2*c*pow(f,2)*r_x*r_z-2*b*e*f*r_x*r_z
                  -2*a*b*f*r_x*r_z+2*c*pow(e,2)*r_x*r_z+2*pow(c,3)*r_x*r_z
                  +2*pow(b,2)*c*r_x*r_z-2*b_y*f*i*r_z-2*b_x*c*i*r_z+2*b_z*pow(f,2)*r_z
                  -2*b_y*e*f*r_z-8*a*b_y*f*r_z+6*b*b_x*f*r_z+2*b_z*pow(e,2)*r_z
                  -8*b_x*c*e*r_z+8*a*b_z*e*r_z+2*b_z*pow(c,2)*r_z+6*b*b_y*c*r_z
                  -2*a*b_x*c*r_z-4*pow(b,2)*b_z*r_z+2*pow(a,2)*b_z*r_z+e*pow(i,2)*pow(r_y,2)
                  -pow(f,2)*i*pow(r_y,2)+4*a*e*i*pow(r_y,2)-4*pow(b,2)*i*pow(r_y,2)-4*a*pow(f,2)*pow(r_y,2)
                  +6*b*c*f*pow(r_y,2)-2*pow(c,2)*e*pow(r_y,2)+pow(a,2)*e*pow(r_y,2)-a*pow(b,2)*pow(r_y,2)
                  +2*b*pow(i,2)*r_x*r_y-2*c*f*i*r_x*r_y+2*b*pow(f,2)*r_x*r_y
                  -2*c*e*f*r_x*r_y-2*a*c*f*r_x*r_y-2*a*b*e*r_x*r_y
                  +2*b*pow(c,2)*r_x*r_y+2*pow(b,3)*r_x*r_y+2*b_y*pow(i,2)*r_y-2*b_z*f*i*r_y
                  +8*a*b_y*i*r_y-8*b*b_x*i*r_y+2*b_y*pow(f,2)*r_y-2*b_z*e*f*r_y
                  +6*b_x*c*f*r_y-8*a*b_z*f*r_y-2*b*b_x*e*r_y-4*b_y*pow(c,2)*r_y
                  +6*b*b_z*c*r_y+2*pow(b,2)*b_y*r_y+2*pow(a,2)*b_y*r_y-2*a*b*b_x*r_y
                  +a*pow(i,2)*pow(r_x,2)+4*a*e*i*pow(r_x,2)-pow(c,2)*i*pow(r_x,2)-4*pow(b,2)*i*pow(r_x,2)
                  -2*a*pow(f,2)*pow(r_x,2)+6*b*c*f*pow(r_x,2)+a*pow(e,2)*pow(r_x,2)-4*pow(c,2)*e*pow(r_x,2)
                  -pow(b,2)*e*pow(r_x,2)+2*b_x*pow(i,2)*r_x+8*b_x*e*i*r_x-2*b_z*c*i*r_x
                  -8*b*b_y*i*r_x-4*b_x*pow(f,2)*r_x+6*b_y*c*f*r_x+6*b*b_z*f*r_x
                  +2*b_x*pow(e,2)*r_x-8*b_z*c*e*r_x-2*b*b_y*e*r_x+2*b_x*pow(c,2)*r_x
                  -2*a*b_z*c*r_x-2*a*b*b_y*r_x+2*pow(b,2)*b_x*r_x+M_0*pow(i,2)+4*M_0*e*i
                  -pow(b_z,2)*i-4*pow(b_y,2)*i-4*pow(b_x,2)*i+4*M_0*a*i-2*M_0*pow(f,2)+6*b_y*b_z*f
                  +M_0*pow(e,2)-4*pow(b_z,2)*e-pow(b_y,2)*e-4*pow(b_x,2)*e+4*M_0*a*e-2*M_0*pow(c,2)
                  +6*b_x*b_z*c-4*a*pow(b_z,2)-4*a*pow(b_y,2)+6*b*b_x*b_y-a*pow(b_x,2)
                  -2*M_0*pow(b,2)+M_0*pow(a,2);
  const G4double power_1 =
    (-2*e*i*pow(r_z,2))-2*a*i*pow(r_z,2)+2*pow(f,2)*pow(r_z,2)+2*pow(c,2)*pow(r_z,2)-4*a*f*r_y*r_z
                     +4*b*c*r_y*r_z+4*b*f*r_x*r_z-4*c*e*r_x*r_z+4*b_y*f*r_z
                     -4*b_z*e*r_z+4*b_x*c*r_z-4*a*b_z*r_z-2*e*i*pow(r_y,2)
                     +2*pow(f,2)*pow(r_y,2)-2*a*e*pow(r_y,2)+2*pow(b,2)*pow(r_y,2)-4*b*i*r_x*r_y
                     +4*c*f*r_x*r_y-4*b_y*i*r_y+4*b_z*f*r_y-4*a*b_y*r_y
                     +4*b*b_x*r_y-2*a*i*pow(r_x,2)-2*a*e*pow(r_x,2)+2*pow(c,2)*pow(r_x,2)
                     +2*pow(b,2)*pow(r_x,2)-4*b_x*i*r_x-4*b_x*e*r_x+4*b_z*c*r_x
                     +4*b*b_y*r_x-2*M_0*i-2*M_0*e+2*pow(b_z,2)+2*pow(b_y,2)+2*pow(b_x,2)
                     -2*M_0*a;
  const G4double power_0 =
    i*pow(r_z,2)+2*f*r_y*r_z+2*c*r_x*r_z+2*b_z*r_z+e*pow(r_y,2)+2*b*r_x*r_y+2*b_y*r_y
                        +a*pow(r_x,2)+2*b_x*r_x+M_0;
  const std::vector<G4double> coefficients {power_6, power_5, power_4, power_3, power_2, power_1, power_0};
  // With this method, you get also get roots at singular value of the determinant, can do an if to check
  const std::vector<G4double> lambdas = polysolve.SolvePolynomial(coefficients);
  std::vector<G4ThreeVector> candidate_points;
  G4ThreeVector candidate_point;
  std::vector<G4double> candidate_distances;
  if (lambdas.size() == 0) {
    // all roots are complex/not found => lambda contour doesn't exist
    // (think it happens if there are infinite solutions, e.g point is at the centre of a sphere)
    // or doesn't touch/intersect the quadric surface
    return std::vector<G4ThreeVector>();
  }
  for (unsigned long int i_loop = 0; i_loop < lambdas.size(); i_loop++) {
    // Process each lambda value
    // Will run into div by zero here if determinant is equal to zero
    G4double cur_lambda = lambdas[i_loop];
    G4double det_denom = ((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))
                     -cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)));
    if (det_denom == 0)
    {
      continue;
    }
    G4double x_component =
    ((pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e))*(r_z+cur_lambda*b_z))
                /det_denom
                +((cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)*(r_y+cur_lambda*b_y))
                 /det_denom
                +(((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))*(r_x+cur_lambda*b_x))
                 /det_denom;
    G4double y_component =
    ((cur_lambda*(1-cur_lambda*a)*f+pow(cur_lambda,2)*b*c)*(r_z+cur_lambda*b_z))
                  /det_denom
                  +(((1-cur_lambda*a)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(c,2))*(r_y+cur_lambda*b_y))
                   /det_denom
                  +((cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)*(r_x+cur_lambda*b_x))
                   /det_denom;
    G4double z_component = (((1-cur_lambda*a)*(1-cur_lambda*e)-pow(cur_lambda,2)*pow(b,2))*(r_z+cur_lambda*b_z))
                /det_denom
                +((cur_lambda*(1-cur_lambda*a)*f+pow(cur_lambda,2)*b*c)*(r_y+cur_lambda*b_y))
                 /det_denom
                +((pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e))*(r_x+cur_lambda*b_x))
                 /det_denom;

    candidate_point.set(x_component, y_component, z_component);
    candidate_points.push_back(candidate_point);
  }
  return candidate_points;
}

G4RotationMatrix G4HalfQuadric::get_Q() const
{
  return this->Q;
}
G4ThreeVector G4HalfQuadric::get_P() const
{
  return this->P;
}
G4double G4HalfQuadric::get_R() const
{
  return this->R;
}
void G4HalfQuadric::set_Q(const G4RotationMatrix& new_Q)
{
  this->Q = new_Q;
}
void G4HalfQuadric::set_P(const G4ThreeVector& new_P)
{
  this->P = new_P;
}
void G4HalfQuadric::set_R(G4double new_R)
{
  this->R = new_R;
}

G4double G4HalfQuadric::ImplicitFunction(const G4ThreeVector& p) const
{
  // x_T Q x + P dot x + R
  // P is a Hep3Vector, Q is a HepRotation/RotationMatrix
  G4double x_T_Q_x = this->Q.xx() * p.x() * p.x() + this->Q.xy() * p.x() * p.y() + this->Q.xz() * p.x() * p.z()
           + this->Q.yx() * p.y() * p.x() + this->Q.yy() * p.y() * p.y() + this->Q.yz() * p.y() * p.z()
           + this->Q.zx() * p.z() * p.x() + this->Q.zy() * p.z() * p.y() + this->Q.zz() * p.z() * p.z();
  G4double P_x = this->P * p;
  return x_T_Q_x + P_x + this->R;
}

std::vector<G4ThreeVector> G4HalfQuadric::RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  // a = v_T Q v
  // b = 2 v_T Q p + v * P
  // c = p_T Q p + R
  if (this->Inside(p) == kSurface)
  {
    return std::vector<G4ThreeVector> {p};
  }
    G4double a = this->Q.xx() * v.x() * v.x() + this->Q.xy() * v.x() * v.y() + this->Q.xz() * v.x() * v.z()
  + this->Q.yx() * v.y() * v.x() + this->Q.yy() * v.y() * v.y() + this->Q.yz() * v.y() * v.z()
  + this->Q.zx() * v.z() * v.x() + this->Q.zy() * v.z() * v.y() + this->Q.zz() * v.z() * v.z();
  if (a == 0) // avoid div by zero when solving -- probably an exotic case
  {
    return std::vector<G4ThreeVector> {p};
  }
  G4double b = 2 * (this->Q.xx() * v.x() * p.x() + this->Q.xy() * v.x() * p.y() + this->Q.xz() * v.x() * p.z()
+ this->Q.yx() * v.y() * p.x() + this->Q.yy() * v.y() * p.y() + this->Q.yz() * v.y() * p.z()
+ this->Q.zx() * v.z() * p.x() + this->Q.zy() * v.z() * p.y() + this->Q.zz() * v.z() * p.z())
  + (this->P * v);
  G4double c = this->Q.xx() * p.x() * p.x() + this->Q.xy() * p.x() * p.y() + this->Q.xz() * p.x() * p.z()
           + this->Q.yx() * p.y() * p.x() + this->Q.yy() * p.y() * p.y() + this->Q.yz() * p.y() * p.z()
           + this->Q.zx() * p.z() * p.x() + this->Q.zy() * p.z() * p.y() + this->Q.zz() * p.z() * p.z()
  + p * this->P + this->R;
  const std::vector<G4double> coefficients {a,b,c};
  std::vector<G4ThreeVector> intersections;
  const G4PolynomialSolve polysolve;
  const std::vector<G4double> v_multipliers = polysolve.SolvePolynomial(coefficients);
  for (auto it = v_multipliers.begin(); it != v_multipliers.end(); ++it)
  {
    intersections.push_back(p + (*it) * v);
  }
  return intersections;
}
G4double G4HalfQuadric::DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  std::vector<G4ThreeVector> intersections = this->RayIntersection(p,v);
  std::vector<G4double> distances;
  for (const auto& intersection : intersections)
  {
    distances.push_back((intersection - p).mag());
  }
  return *(std::min_element(distances.begin(), distances.end()));
}

std::ostream& G4HalfQuadric::StreamInfo(std::ostream& os) const
{
    os << "--- Half Surface Quadric ---" << std::endl;
    return os;
}

G4ThreeVector G4HalfQuadric::SurfaceNormal(const G4ThreeVector& p) const
{
  G4double grad_x = 2 * (Q.xx() * p[0] + Q.xy() * p[1] + Q.xz() * p[2]);
  G4double grad_y = 2 * (Q.yy() * p[0] + Q.yx() * p[1] + Q.yz() * p[2]);
  G4double grad_z = 2 * (Q.zz() * p[0] + Q.zx() * p[1] + Q.zy() * p[2]);
  G4ThreeVector grad {grad_x, grad_y, grad_z};
  return grad;
}
G4ThreeVector G4HalfQuadric::ClosestSurfacePoint(const G4ThreeVector& p) const
{
  const std::vector<G4ThreeVector> candidate_points = GetLMContourIntersections(p);
  std::vector<G4double> distances;
  for (const auto& candidate_point : candidate_points)
  {
    distances.push_back((candidate_point - p).mag());
  }
  const auto min_distance = std::min_element(distances.begin(), distances.end());
  const G4int element_index = std::distance(distances.begin(), min_distance);
  return candidate_points[element_index];
}

void G4HalfQuadric::Transform(const G4AffineTransform& t)
{
  G4RotationMatrix rotate = t.NetRotation(); // M_ij
  G4ThreeVector translate = t.NetTranslation(); // T_j
  // this is the vector that the ORIGIN moves by, if a minus sign is missed then that is because the translation is of the object and not its occupied space.
  G4RotationMatrix Q_prime = _T(rotate) * this->Q * rotate;
  G4ThreeVector P_prime = 2 * (_T(rotate) * this->Q * translate) + (_T(rotate) * (this->P));
  G4double R_prime = (this->P * translate) + (_T(this->Q) * translate) * translate + this->R; // was missing this->R
  this->Q = Q_prime;
  this->P = P_prime;
  this->R = R_prime;
  return;
}

G4double G4HalfQuadric::DistanceToSurface(const G4ThreeVector& p) const
{
  std::vector<G4ThreeVector> candidate_points = GetLMContourIntersections(p);
  std::vector<G4double> candidate_distances;
  for (const auto& candidate_point : candidate_points) // oh cool, copilot generated this
  {
    candidate_distances.push_back((candidate_point - p).mag());
  }
  const auto min_distance = std::min_element(candidate_distances.begin(), candidate_distances.end());
  return *min_distance;
}
