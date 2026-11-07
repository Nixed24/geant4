#include <algorithm>
#include <cmath>
#include "G4HalfQuadric.hh"
#include "G4PolynomialSolve.hh"

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
    G4double a = this->Q.xx() * v.x() * v.x() + this->Q.xy() * v.x() * v.y() + this->Q.xz() * v.x() * v.z()
  + this->Q.yx() * v.y() * v.x() + this->Q.yy() * v.y() * v.y() + this->Q.yz() * v.y() * v.z()
  + this->Q.zx() * v.z() * v.x() + this->Q.zy() * v.z() * v.y() + this->Q.zz() * v.z() * v.z();
  if (a == 0) // check a = 0 for parallel i think
  {
    return std::vector<G4ThreeVector> {p};
  }
  G4double b = this->Q.xx() * v.x() * p.x() + this->Q.xy() * v.x() * p.y() + this->Q.xz() * v.x() * p.z()
+ this->Q.yx() * v.y() * p.x() + this->Q.yy() * v.y() * p.y() + this->Q.yz() * v.y() * p.z()
+ this->Q.zx() * v.z() * p.x() + this->Q.zy() * v.z() * p.y() + this->Q.zz() * v.z() * p.z()
  + (this->P * v);
  G4double c = this->Q.xx() * p.x() * p.x() + this->Q.xy() * p.x() * p.y() + this->Q.xz() * p.x() * p.z()
           + this->Q.yx() * p.y() * p.x() + this->Q.yy() * p.y() * p.y() + this->Q.yz() * p.y() * p.z()
           + this->Q.zx() * p.z() * p.x() + this->Q.zy() * p.z() * p.y() + this->Q.zz() * p.z() * p.z()
  + this->R;

  G4double discriminant = b * b - 4 * a * c;

  if (discriminant < 0)
  {
    return std::vector<G4ThreeVector>();
  }

  G4double root_discriminant = pow(discriminant, 0.5);
  G4double lambda_1 = (-1 * b + root_discriminant) / (2*a);
  G4double lambda_2 = (-1 * b - root_discriminant) / (2*a);
  std::vector<G4ThreeVector> intersections;
  if (lambda_1 >= 0)
  {
    intersections.push_back(p + (lambda_1 * v));
  }
  if (lambda_2 >= 0)
  {
    intersections.push_back(p + (lambda_2 * v));
  }
  return intersections;
}
G4double G4HalfQuadric::DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  std::vector<G4ThreeVector> intersections = this->RayIntersection(p,v);
  // probably redundant
  return std::min(intersections[0].mag(), intersections[1].mag());
}

std::ostream& G4HalfQuadric::StreamInfo(std::ostream& os) const
{
    os << "--- Half Surface Quadric ---" << std::endl;
    return os;
}

G4double G4HalfQuadric::DistanceToSurface(const G4ThreeVector& p) const
{
  G4PolynomialSolve polysolve;
  const G4double r_x = p.x(), r_y = p.y(), r_z = p.z();
  const G4double M_0 = this->R;
  const G4double P_x = this->P.x() * 0.5, P_y = this->P.y() * 0.5, P_z = this->P.z() * 0.5;
  const G4double a = Q.xx(), b = Q.xy(), c = Q.xz(), e = Q.yy(), f = Q.yz(), i = Q.zz();
  // Edge cases: Determinant = 0 | Determinant is a cubic polynomial in lambda, could add a check for lambda being at one of its roots
  // I don't know what the determinant being zero would actually correspond to, it's not a plane because that would just give identity
  // but I think it could mean that there are infinitely many solutions for the closest point OR maybe the p is already on the surface
  // x of lambda will also div by zero so we'll know
  // Z -> lambda (lagrange multiplier)
  // const G4double characteristic_determinant = (-pow(Z,3)*a*e*i)+pow(Z,2)*e*i+pow(Z,3)*pow(b,2)*i+pow(Z,2)*a*i-Z*i+pow(Z,3)*a*pow(f,2)-pow(Z,2)*pow(f,2)-2*pow(Z,3)*b*c*f+pow(Z,3)*pow(c,2)*e+pow(Z,2)*a*e-Z*e-pow(Z,2)*pow(c,2)-pow(Z,2)*pow(b,2)-Z*a+1;
  // It's a LONG expression
  // 07/11/26 these terms aren't correct; pull the correct ones from poly_redo_test_demo; that also uses b (P/2) instead of P
  const G4double power_6 = M_0*pow(a,2)*pow(e,2)*pow(i,2)-2*M_0*a*pow(b,2)*e*pow(i,2)+M_0*pow(b,4)*pow(i,2)-2*M_0*pow(a,2)*e*pow(f,2)*i+2*M_0*a*pow(b,2)*pow(f,2)*i+4*M_0*a*b*c*e*f*i-4*M_0*pow(b,3)*c*f*i-2*M_0*a*pow(c,2)*pow(e,2)*i+2*M_0*pow(b,2)*pow(c,2)*e*i+M_0*pow(a,2)*pow(f,4)-4*M_0*a*b*c*pow(f,3)+2*M_0*a*pow(c,2)*e*pow(f,2)+4*M_0*pow(b,2)*pow(c,2)*pow(f,2)-4*M_0*b*pow(c,3)*e*f+M_0*pow(c,4)*pow(e,2);
  const G4double power_5 = P_x*a*b*e*f*i*r_z-P_y*pow(a,2)*e*f*i*r_z-P_x*pow(b,3)*f*i*r_z+P_y*a*pow(b,2)*f*i*r_z-P_x*a*c*pow(e,2)*i*r_z+P_z*pow(a,2)*pow(e,2)*i*r_z+P_x*pow(b,2)*c*e*i*r_z+P_y*a*b*c*e*i*r_z-2*P_z*a*pow(b,2)*e*i*r_z-P_y*pow(b,3)*c*i*r_z+P_z*pow(b,4)*i*r_z-P_x*a*b*pow(f,3)*r_z+P_y*pow(a,2)*pow(f,3)*r_z+P_x*a*c*e*pow(f,2)*r_z-P_z*pow(a,2)*e*pow(f,2)*r_z+2*P_x*pow(b,2)*c*pow(f,2)*r_z-3*P_y*a*b*c*pow(f,2)*r_z
                         +P_z*a*pow(b,2)*pow(f,2)*r_z-3*P_x*b*pow(c,2)*e*f*r_z+P_y*a*pow(c,2)*e*f*r_z+2*P_z*a*b*c*e*f*r_z+2*P_y*pow(b,2)*pow(c,2)*f*r_z-2*P_z*pow(b,3)*c*f*r_z+P_x*pow(c,3)*pow(e,2)*r_z-P_z*a*pow(c,2)*pow(e,2)*r_z-P_y*b*pow(c,3)*e*r_z+P_z*pow(b,2)*pow(c,2)*e*r_z-P_x*a*b*e*pow(i,2)*r_y+P_y*pow(a,2)*e*pow(i,2)*r_y+P_x*pow(b,3)*pow(i,2)*r_y-P_y*a*pow(b,2)*pow(i,2)*r_y+P_x*a*b*pow(f,2)*i*r_y
                         -P_y*pow(a,2)*pow(f,2)*i*r_y+P_x*a*c*e*f*i*r_y-P_z*pow(a,2)*e*f*i*r_y-3*P_x*pow(b,2)*c*f*i*r_y+2*P_y*a*b*c*f*i*r_y+P_z*a*pow(b,2)*f*i*r_y+P_x*b*pow(c,2)*e*i*r_y-2*P_y*a*pow(c,2)*e*i*r_y+P_z*a*b*c*e*i*r_y+P_y*pow(b,2)*pow(c,2)*i*r_y-P_z*pow(b,3)*c*i*r_y-P_x*a*c*pow(f,3)*r_y+P_z*pow(a,2)*pow(f,3)*r_y+2*P_x*b*pow(c,2)*pow(f,2)*r_y+P_y*a*pow(c,2)*pow(f,2)*r_y
                         -3*P_z*a*b*c*pow(f,2)*r_y-P_x*pow(c,3)*e*f*r_y+P_z*a*pow(c,2)*e*f*r_y-2*P_y*b*pow(c,3)*f*r_y+2*P_z*pow(b,2)*pow(c,2)*f*r_y+P_y*pow(c,4)*e*r_y-P_z*b*pow(c,3)*e*r_y+P_x*a*pow(e,2)*pow(i,2)*r_x-P_x*pow(b,2)*e*pow(i,2)*r_x-P_y*a*b*e*pow(i,2)*r_x+P_y*pow(b,3)*pow(i,2)*r_x-2*P_x*a*e*pow(f,2)*i*r_x+P_x*pow(b,2)*pow(f,2)*i*r_x+P_y*a*b*pow(f,2)*i*r_x+2*P_x*b*c*e*f*i*r_x
                         +P_y*a*c*e*f*i*r_x+P_z*a*b*e*f*i*r_x-3*P_y*pow(b,2)*c*f*i*r_x-P_z*pow(b,3)*f*i*r_x-P_x*pow(c,2)*pow(e,2)*i*r_x-P_z*a*c*pow(e,2)*i*r_x+P_y*b*pow(c,2)*e*i*r_x+P_z*pow(b,2)*c*e*i*r_x+P_x*a*pow(f,4)*r_x-2*P_x*b*c*pow(f,3)*r_x-P_y*a*c*pow(f,3)*r_x-P_z*a*b*pow(f,3)*r_x+P_x*pow(c,2)*e*pow(f,2)*r_x+P_z*a*c*e*pow(f,2)*r_x+2*P_y*b*pow(c,2)*pow(f,2)*r_x+2*P_z*pow(b,2)*c*pow(f,2)*r_x
                         -P_y*pow(c,3)*e*f*r_x-3*P_z*b*pow(c,2)*e*f*r_x+P_z*pow(c,3)*pow(e,2)*r_x-2*M_0*a*pow(e,2)*pow(i,2)+pow(P_x,2)*pow(e,2)*pow(i,2)+2*M_0*pow(b,2)*e*pow(i,2)-2*P_x*P_y*b*e*pow(i,2)-2*M_0*pow(a,2)*e*pow(i,2)+2*M_0*a*pow(b,2)*pow(i,2)+pow(P_y,2)*pow(b,2)*pow(i,2)+pow(P_x,2)*pow(b,2)*pow(i,2)-2*P_x*P_y*a*b*pow(i,2)+pow(P_y,2)*pow(a,2)*pow(i,2)+4*M_0*a*e*pow(f,2)*i-2*pow(P_x,2)*e*pow(f,2)*i-2*M_0*pow(b,2)*pow(f,2)*i+2*P_x*P_y*b*pow(f,2)*i
                         +2*M_0*pow(a,2)*pow(f,2)*i-4*M_0*b*c*e*f*i+2*P_x*P_y*c*e*f*i+2*P_x*P_z*b*e*f*i-4*M_0*a*b*c*f*i-2*pow(P_y,2)*b*c*f*i-2*pow(P_x,2)*b*c*f*i+2*P_x*P_y*a*c*f*i-2*P_y*P_z*pow(b,2)*f*i+2*P_x*P_z*a*b*f*i-2*P_y*P_z*pow(a,2)*f*i+2*M_0*pow(c,2)*pow(e,2)*i-2*P_x*P_z*c*pow(e,2)*i-2*M_0*pow(a,2)*pow(e,2)*i+4*M_0*a*pow(c,2)*e*i+2*P_y*P_z*b*c*e*i
                         +4*M_0*a*pow(b,2)*e*i-2*M_0*pow(b,2)*pow(c,2)*i+2*P_x*P_y*b*pow(c,2)*i-2*pow(P_y,2)*a*pow(c,2)*i-2*P_x*P_z*pow(b,2)*c*i+2*P_y*P_z*a*b*c*i-2*M_0*pow(b,4)*i-2*M_0*a*pow(f,4)+pow(P_x,2)*pow(f,4)+4*M_0*b*c*pow(f,3)-2*P_x*P_y*c*pow(f,3)-2*P_x*P_z*b*pow(f,3)-2*M_0*pow(c,2)*e*pow(f,2)+2*P_x*P_z*c*e*pow(f,2)+2*M_0*pow(a,2)*e*pow(f,2)-2*M_0*a*pow(c,2)*pow(f,2)+pow(P_y,2)*pow(c,2)*pow(f,2)+pow(P_x,2)*pow(c,2)*pow(f,2)
                         +2*P_y*P_z*b*c*pow(f,2)-2*P_x*P_z*a*c*pow(f,2)-2*M_0*a*pow(b,2)*pow(f,2)+pow(P_z,2)*pow(b,2)*pow(f,2)+pow(P_x,2)*pow(b,2)*pow(f,2)-2*P_x*P_y*a*b*pow(f,2)+pow(P_z,2)*pow(a,2)*pow(f,2)+pow(P_y,2)*pow(a,2)*pow(f,2)-2*P_y*P_z*pow(c,2)*e*f-4*M_0*a*b*c*e*f-2*pow(P_z,2)*b*c*e*f-2*pow(P_x,2)*b*c*e*f+2*P_x*P_y*a*c*e*f+2*P_x*P_z*a*b*e*f-2*P_y*P_z*pow(a,2)*e*f+4*M_0*b*pow(c,3)*f-2*P_x*P_y*pow(c,3)*f
                         +2*P_x*P_z*b*pow(c,2)*f+2*P_y*P_z*a*pow(c,2)*f+4*M_0*pow(b,3)*c*f+2*P_x*P_y*pow(b,2)*c*f-2*pow(P_z,2)*a*b*c*f-2*pow(P_y,2)*a*b*c*f-2*P_x*P_z*pow(b,3)*f+2*P_y*P_z*a*pow(b,2)*f+2*M_0*a*pow(c,2)*pow(e,2)+pow(P_z,2)*pow(c,2)*pow(e,2)+pow(P_x,2)*pow(c,2)*pow(e,2)-2*P_x*P_z*a*c*pow(e,2)+pow(P_z,2)*pow(a,2)*pow(e,2)-2*M_0*pow(c,4)*e-2*M_0*pow(b,2)*pow(c,2)*e-2*P_x*P_y*b*pow(c,2)*e+2*P_x*P_z*pow(b,2)*c*e
                         +2*P_y*P_z*a*b*c*e-2*pow(P_z,2)*a*pow(b,2)*e+pow(P_y,2)*pow(c,4)-2*P_y*P_z*b*pow(c,3)+pow(P_z,2)*pow(b,2)*pow(c,2)+pow(P_y,2)*pow(b,2)*pow(c,2)-2*P_y*P_z*pow(b,3)*c+pow(P_z,2)*pow(b,4);
  const G4double power_4 = pow(a,2)*pow(e,2)*i*pow(r_z,2)-2*a*pow(b,2)*e*i*pow(r_z,2)+pow(b,4)*i*pow(r_z,2)-pow(a,2)*e*pow(f,2)*pow(r_z,2)+a*pow(b,2)*pow(f,2)*pow(r_z,2)+2*a*b*c*e*f*pow(r_z,2)-2*pow(b,3)*c*f*pow(r_z,2)-a*pow(c,2)*pow(e,2)*pow(r_z,2)+pow(b,2)*pow(c,2)*e*pow(r_z,2)-2*pow(a,2)*e*f*i*r_y*r_z+2*a*pow(b,2)*f*i*r_y*r_z+2*a*b*c*e*i*r_y*r_z-2*pow(b,3)*c*i*r_y*r_z+2*pow(a,2)*pow(f,3)*r_y*r_z-6*a*b*c*pow(f,2)*r_y*r_z+2*a*pow(c,2)*e*f*r_y*r_z+4*pow(b,2)*pow(c,2)*f*r_y*r_z
                       -2*b*pow(c,3)*e*r_y*r_z+2*a*b*e*f*i*r_x*r_z-2*pow(b,3)*f*i*r_x*r_z-2*a*c*pow(e,2)*i*r_x*r_z+2*pow(b,2)*c*e*i*r_x*r_z-2*a*b*pow(f,3)*r_x*r_z+2*a*c*e*pow(f,2)*r_x*r_z+4*pow(b,2)*c*pow(f,2)*r_x*r_z-6*b*pow(c,2)*e*f*r_x*r_z+2*pow(c,3)*pow(e,2)*r_x*r_z+P_x*b*e*f*i*r_z+2*P_y*a*e*f*i*r_z-3*P_y*pow(b,2)*f*i*r_z+P_x*a*b*f*i*r_z-P_y*pow(a,2)*f*i*r_z
                       -P_x*c*pow(e,2)*i*r_z-2*P_z*a*pow(e,2)*i*r_z+P_y*b*c*e*i*r_z+2*P_x*a*c*e*i*r_z+2*P_z*pow(b,2)*e*i*r_z-2*P_z*pow(a,2)*e*i*r_z-3*P_x*pow(b,2)*c*i*r_z+P_y*a*b*c*i*r_z+2*P_z*a*pow(b,2)*i*r_z-P_x*b*pow(f,3)*r_z-2*P_y*a*pow(f,3)*r_z+P_x*c*e*pow(f,2)*r_z+2*P_z*a*e*pow(f,2)*r_z+5*P_y*b*c*pow(f,2)*r_z-3*P_x*a*c*pow(f,2)*r_z+P_z*pow(b,2)*pow(f,2)*r_z
                       +3*P_z*pow(a,2)*pow(f,2)*r_z-3*P_y*pow(c,2)*e*f*r_z-6*P_z*b*c*e*f*r_z+P_x*a*b*e*f*r_z-P_y*pow(a,2)*e*f*r_z+5*P_x*b*pow(c,2)*f*r_z+P_y*a*pow(c,2)*f*r_z-6*P_z*a*b*c*f*r_z-P_x*pow(b,3)*f*r_z+P_y*a*pow(b,2)*f*r_z+3*P_z*pow(c,2)*pow(e,2)*r_z-P_x*a*c*pow(e,2)*r_z+P_z*pow(a,2)*pow(e,2)*r_z-2*P_x*pow(c,3)*e*r_z+2*P_z*a*pow(c,2)*e*r_z+P_x*pow(b,2)*c*e*r_z+P_y*a*b*c*e*r_z
                       -2*P_z*a*pow(b,2)*e*r_z-P_y*b*pow(c,3)*r_z+P_z*pow(b,2)*pow(c,2)*r_z-P_y*pow(b,3)*c*r_z+P_z*pow(b,4)*r_z+pow(a,2)*e*pow(i,2)*pow(r_y,2)-a*pow(b,2)*pow(i,2)*pow(r_y,2)-pow(a,2)*pow(f,2)*i*pow(r_y,2)+2*a*b*c*f*i*pow(r_y,2)-2*a*pow(c,2)*e*i*pow(r_y,2)+pow(b,2)*pow(c,2)*i*pow(r_y,2)+a*pow(c,2)*pow(f,2)*pow(r_y,2)-2*b*pow(c,3)*f*pow(r_y,2)+pow(c,4)*e*pow(r_y,2)-2*a*b*e*pow(i,2)*r_x*r_y+2*pow(b,3)*pow(i,2)*r_x*r_y+2*a*b*pow(f,2)*i*r_x*r_y
                       +2*a*c*e*f*i*r_x*r_y-6*pow(b,2)*c*f*i*r_x*r_y+2*b*pow(c,2)*e*i*r_x*r_y-2*a*c*pow(f,3)*r_x*r_y+4*b*pow(c,2)*pow(f,2)*r_x*r_y-2*pow(c,3)*e*f*r_x*r_y-P_x*b*e*pow(i,2)*r_y-2*P_y*a*e*pow(i,2)*r_y+3*P_y*pow(b,2)*pow(i,2)*r_y-P_x*a*b*pow(i,2)*r_y+P_y*pow(a,2)*pow(i,2)*r_y+P_x*b*pow(f,2)*i*r_y+2*P_y*a*pow(f,2)*i*r_y+P_x*c*e*f*i*r_y+2*P_z*a*e*f*i*r_y-6*P_y*b*c*f*i*r_y
                       +P_x*a*c*f*i*r_y-3*P_z*pow(b,2)*f*i*r_y-P_z*pow(a,2)*f*i*r_y+2*P_y*pow(c,2)*e*i*r_y+P_z*b*c*e*i*r_y+2*P_x*a*b*e*i*r_y-2*P_y*pow(a,2)*e*i*r_y+P_x*b*pow(c,2)*i*r_y-2*P_y*a*pow(c,2)*i*r_y+P_z*a*b*c*i*r_y-2*P_x*pow(b,3)*i*r_y+2*P_y*a*pow(b,2)*i*r_y-P_x*c*pow(f,3)*r_y-2*P_z*a*pow(f,3)*r_y+P_y*pow(c,2)*pow(f,2)*r_y+5*P_z*b*c*pow(f,2)*r_y-3*P_x*a*b*pow(f,2)*r_y
                       +3*P_y*pow(a,2)*pow(f,2)*r_y-3*P_z*pow(c,2)*e*f*r_y+P_x*a*c*e*f*r_y-P_z*pow(a,2)*e*f*r_y-P_x*pow(c,3)*f*r_y+P_z*a*pow(c,2)*f*r_y+5*P_x*pow(b,2)*c*f*r_y-6*P_y*a*b*c*f*r_y+P_z*a*pow(b,2)*f*r_y-3*P_x*b*pow(c,2)*e*r_y+2*P_y*a*pow(c,2)*e*r_y+P_z*a*b*c*e*r_y+P_y*pow(c,4)*r_y-P_z*b*pow(c,3)*r_y+P_y*pow(b,2)*pow(c,2)*r_y-P_z*pow(b,3)*c*r_y+a*pow(e,2)*pow(i,2)*pow(r_x,2)
                       -pow(b,2)*e*pow(i,2)*pow(r_x,2)-2*a*e*pow(f,2)*i*pow(r_x,2)+pow(b,2)*pow(f,2)*i*pow(r_x,2)+2*b*c*e*f*i*pow(r_x,2)-pow(c,2)*pow(e,2)*i*pow(r_x,2)+a*pow(f,4)*pow(r_x,2)-2*b*c*pow(f,3)*pow(r_x,2)+pow(c,2)*e*pow(f,2)*pow(r_x,2)+P_x*pow(e,2)*pow(i,2)*r_x-P_y*b*e*pow(i,2)*r_x-2*P_x*a*e*pow(i,2)*r_x+3*P_x*pow(b,2)*pow(i,2)*r_x-P_y*a*b*pow(i,2)*r_x-2*P_x*e*pow(f,2)*i*r_x+P_y*b*pow(f,2)*i*r_x+2*P_x*a*pow(f,2)*i*r_x+P_y*c*e*f*i*r_x
                       +P_z*b*e*f*i*r_x-6*P_x*b*c*f*i*r_x+P_y*a*c*f*i*r_x+P_z*a*b*f*i*r_x-P_z*c*pow(e,2)*i*r_x-2*P_x*a*pow(e,2)*i*r_x+2*P_x*pow(c,2)*e*i*r_x+2*P_z*a*c*e*i*r_x+2*P_x*pow(b,2)*e*i*r_x+2*P_y*a*b*e*i*r_x+P_y*b*pow(c,2)*i*r_x-3*P_z*pow(b,2)*c*i*r_x-2*P_y*pow(b,3)*i*r_x+P_x*pow(f,4)*r_x-P_y*c*pow(f,3)*r_x-P_z*b*pow(f,3)*r_x+P_z*c*e*pow(f,2)*r_x
                       +2*P_x*a*e*pow(f,2)*r_x+P_x*pow(c,2)*pow(f,2)*r_x-3*P_z*a*c*pow(f,2)*r_x+P_x*pow(b,2)*pow(f,2)*r_x-3*P_y*a*b*pow(f,2)*r_x-6*P_x*b*c*e*f*r_x+P_y*a*c*e*f*r_x+P_z*a*b*e*f*r_x-P_y*pow(c,3)*f*r_x+5*P_z*b*pow(c,2)*f*r_x+5*P_y*pow(b,2)*c*f*r_x-P_z*pow(b,3)*f*r_x+3*P_x*pow(c,2)*pow(e,2)*r_x-P_z*a*c*pow(e,2)*r_x-2*P_z*pow(c,3)*e*r_x-3*P_y*b*pow(c,2)*e*r_x+P_z*pow(b,2)*c*e*r_x
                       +M_0*pow(e,2)*pow(i,2)+4*M_0*a*e*pow(i,2)-2*pow(P_x,2)*e*pow(i,2)-2*M_0*pow(b,2)*pow(i,2)+4*P_x*P_y*b*pow(i,2)+M_0*pow(a,2)*pow(i,2)-2*pow(P_y,2)*a*pow(i,2)-2*M_0*e*pow(f,2)*i-4*M_0*a*pow(f,2)*i+2*pow(P_x,2)*pow(f,2)*i+4*M_0*b*c*f*i-4*P_x*P_y*c*f*i-4*P_x*P_z*b*f*i+4*P_y*P_z*a*f*i+4*M_0*a*pow(e,2)*i-2*pow(P_x,2)*pow(e,2)*i-4*M_0*pow(c,2)*e*i+4*P_x*P_z*c*e*i-4*M_0*pow(b,2)*e*i
                       +4*P_x*P_y*b*e*i+4*M_0*pow(a,2)*e*i-2*M_0*a*pow(c,2)*i+2*pow(P_y,2)*pow(c,2)*i-4*P_y*P_z*b*c*i-4*M_0*a*pow(b,2)*i-2*pow(P_y,2)*pow(b,2)*i-2*pow(P_x,2)*pow(b,2)*i+4*P_x*P_y*a*b*i-2*pow(P_y,2)*pow(a,2)*i+M_0*pow(f,4)-4*M_0*a*e*pow(f,2)+2*pow(P_x,2)*e*pow(f,2)+2*M_0*pow(c,2)*pow(f,2)+2*M_0*pow(b,2)*pow(f,2)-2*M_0*pow(a,2)*pow(f,2)-2*pow(P_z,2)*a*pow(f,2)-2*pow(P_y,2)*a*pow(f,2)+4*M_0*b*c*e*f-4*P_x*P_y*c*e*f
                       -4*P_x*P_z*b*e*f+4*P_y*P_z*a*e*f+4*M_0*a*b*c*f+4*pow(P_z,2)*b*c*f+4*pow(P_y,2)*b*c*f+4*pow(P_x,2)*b*c*f-4*P_x*P_y*a*c*f-4*P_x*P_z*a*b*f+4*P_y*P_z*pow(a,2)*f-2*M_0*pow(c,2)*pow(e,2)+4*P_x*P_z*c*pow(e,2)+M_0*pow(a,2)*pow(e,2)-2*pow(P_z,2)*a*pow(e,2)-4*M_0*a*pow(c,2)*e-2*pow(P_z,2)*pow(c,2)*e-2*pow(P_x,2)*pow(c,2)*e-4*P_y*P_z*b*c*e+4*P_x*P_z*a*c*e-2*M_0*a*pow(b,2)*e
                       +2*pow(P_z,2)*pow(b,2)*e-2*pow(P_z,2)*pow(a,2)*e+M_0*pow(c,4)+2*M_0*pow(b,2)*pow(c,2)+2*pow(P_y,2)*a*pow(c,2)-4*P_y*P_z*a*b*c+M_0*pow(b,4)+2*pow(P_z,2)*a*pow(b,2);
  const G4double power_3 = (-2*a*pow(e,2)*i*pow(r_z,2))+2*pow(b,2)*e*i*pow(r_z,2)-2*pow(a,2)*e*i*pow(r_z,2)+2*a*pow(b,2)*i*pow(r_z,2)+2*a*e*pow(f,2)*pow(r_z,2)+2*pow(a,2)*pow(f,2)*pow(r_z,2)-4*b*c*e*f*pow(r_z,2)-4*a*b*c*f*pow(r_z,2)+2*pow(c,2)*pow(e,2)*pow(r_z,2)+2*a*pow(c,2)*e*pow(r_z,2)+4*a*e*f*i*r_y*r_z-4*pow(b,2)*f*i*r_y*r_z-4*a*pow(f,3)*r_y*r_z+8*b*c*pow(f,2)*r_y*r_z-4*pow(c,2)*e*f*r_y*r_z+4*a*c*e*i*r_x*r_z-4*pow(b,2)*c*i*r_x*r_z
                          -4*a*c*pow(f,2)*r_x*r_z+8*b*pow(c,2)*f*r_x*r_z-4*pow(c,3)*e*r_x*r_z-P_y*e*f*i*r_z-3*P_x*b*f*i*r_z+2*P_y*a*f*i*r_z+P_z*pow(e,2)*i*r_z+2*P_x*c*e*i*r_z+4*P_z*a*e*i*r_z-3*P_y*b*c*i*r_z-P_x*a*c*i*r_z-2*P_z*pow(b,2)*i*r_z+P_z*pow(a,2)*i*r_z+P_y*pow(f,3)*r_z-P_z*e*pow(f,2)*r_z+P_x*c*pow(f,2)*r_z-6*P_z*a*pow(f,2)*r_z-3*P_x*b*e*f*r_z
                          +2*P_y*a*e*f*r_z+P_y*pow(c,2)*f*r_z+10*P_z*b*c*f*r_z+P_y*pow(b,2)*f*r_z-3*P_x*a*b*f*r_z+3*P_y*pow(a,2)*f*r_z+3*P_x*c*pow(e,2)*r_z-2*P_z*a*pow(e,2)*r_z-6*P_z*pow(c,2)*e*r_z-3*P_y*b*c*e*r_z+2*P_x*a*c*e*r_z+2*P_z*pow(b,2)*e*r_z-2*P_z*pow(a,2)*e*r_z+P_x*pow(c,3)*r_z-P_z*a*pow(c,2)*r_z+P_x*pow(b,2)*c*r_z-3*P_y*a*b*c*r_z+2*P_z*a*pow(b,2)*r_z
                          -2*a*e*pow(i,2)*pow(r_y,2)+2*pow(b,2)*pow(i,2)*pow(r_y,2)+2*a*pow(f,2)*i*pow(r_y,2)-4*b*c*f*i*pow(r_y,2)+2*pow(c,2)*e*i*pow(r_y,2)-2*pow(a,2)*e*i*pow(r_y,2)+2*a*pow(b,2)*i*pow(r_y,2)+2*pow(a,2)*pow(f,2)*pow(r_y,2)-4*a*b*c*f*pow(r_y,2)+2*a*pow(c,2)*e*pow(r_y,2)+4*a*b*e*i*r_x*r_y-4*pow(b,3)*i*r_x*r_y-4*a*b*pow(f,2)*r_x*r_y+8*pow(b,2)*c*f*r_x*r_y-4*b*pow(c,2)*e*r_x*r_y+P_y*e*pow(i,2)*r_y+3*P_x*b*pow(i,2)*r_y
                          -2*P_y*a*pow(i,2)*r_y-P_y*pow(f,2)*i*r_y-P_z*e*f*i*r_y-3*P_x*c*f*i*r_y+2*P_z*a*f*i*r_y+2*P_x*b*e*i*r_y+4*P_y*a*e*i*r_y+2*P_y*pow(c,2)*i*r_y-3*P_z*b*c*i*r_y-6*P_y*pow(b,2)*i*r_y+2*P_x*a*b*i*r_y-2*P_y*pow(a,2)*i*r_y+P_z*pow(f,3)*r_y+P_x*b*pow(f,2)*r_y-6*P_y*a*pow(f,2)*r_y-3*P_x*c*e*f*r_y+2*P_z*a*e*f*r_y+P_z*pow(c,2)*f*r_y
                          +10*P_y*b*c*f*r_y-3*P_x*a*c*f*r_y+P_z*pow(b,2)*f*r_y+3*P_z*pow(a,2)*f*r_y-2*P_y*pow(c,2)*e*r_y-3*P_z*b*c*e*r_y-P_x*a*b*e*r_y+P_y*pow(a,2)*e*r_y+P_x*b*pow(c,2)*r_y+2*P_y*a*pow(c,2)*r_y-3*P_z*a*b*c*r_y+P_x*pow(b,3)*r_y-P_y*a*pow(b,2)*r_y-2*a*e*pow(i,2)*pow(r_x,2)+2*pow(b,2)*pow(i,2)*pow(r_x,2)+2*a*pow(f,2)*i*pow(r_x,2)-4*b*c*f*i*pow(r_x,2)-2*a*pow(e,2)*i*pow(r_x,2)
                          +2*pow(c,2)*e*i*pow(r_x,2)+2*pow(b,2)*e*i*pow(r_x,2)+2*a*e*pow(f,2)*pow(r_x,2)-4*b*c*e*f*pow(r_x,2)+2*pow(c,2)*pow(e,2)*pow(r_x,2)-2*P_x*e*pow(i,2)*r_x+3*P_y*b*pow(i,2)*r_x+P_x*a*pow(i,2)*r_x+2*P_x*pow(f,2)*i*r_x-3*P_y*c*f*i*r_x-3*P_z*b*f*i*r_x-2*P_x*pow(e,2)*i*r_x+2*P_z*c*e*i*r_x+2*P_y*b*e*i*r_x+4*P_x*a*e*i*r_x-P_x*pow(c,2)*i*r_x-P_z*a*c*i*r_x-6*P_x*pow(b,2)*i*r_x
                          +2*P_y*a*b*i*r_x+2*P_x*e*pow(f,2)*r_x+P_z*c*pow(f,2)*r_x+P_y*b*pow(f,2)*r_x-2*P_x*a*pow(f,2)*r_x-3*P_y*c*e*f*r_x-3*P_z*b*e*f*r_x+10*P_x*b*c*f*r_x-3*P_y*a*c*f*r_x-3*P_z*a*b*f*r_x+3*P_z*c*pow(e,2)*r_x+P_x*a*pow(e,2)*r_x-6*P_x*pow(c,2)*e*r_x+2*P_z*a*c*e*r_x-P_x*pow(b,2)*e*r_x-P_y*a*b*e*r_x+P_z*pow(c,3)*r_x+P_y*b*pow(c,2)*r_x
                          +P_z*pow(b,2)*c*r_x+P_y*pow(b,3)*r_x-2*M_0*e*pow(i,2)-2*M_0*a*pow(i,2)+pow(P_y,2)*pow(i,2)+pow(P_x,2)*pow(i,2)+2*M_0*pow(f,2)*i-2*P_y*P_z*f*i-2*M_0*pow(e,2)*i-8*M_0*a*e*i+4*pow(P_x,2)*e*i+2*M_0*pow(c,2)*i-2*P_x*P_z*c*i+4*M_0*pow(b,2)*i-8*P_x*P_y*b*i-2*M_0*pow(a,2)*i+4*pow(P_y,2)*a*i+2*M_0*e*pow(f,2)+4*M_0*a*pow(f,2)+pow(P_z,2)*pow(f,2)+pow(P_y,2)*pow(f,2)-2*pow(P_x,2)*pow(f,2)-2*P_y*P_z*e*f
                          -4*M_0*b*c*f+6*P_x*P_y*c*f+6*P_x*P_z*b*f-8*P_y*P_z*a*f-2*M_0*a*pow(e,2)+pow(P_z,2)*pow(e,2)+pow(P_x,2)*pow(e,2)+4*M_0*pow(c,2)*e-8*P_x*P_z*c*e+2*M_0*pow(b,2)*e-2*P_x*P_y*b*e-2*M_0*pow(a,2)*e+4*pow(P_z,2)*a*e+2*M_0*a*pow(c,2)+pow(P_z,2)*pow(c,2)-2*pow(P_y,2)*pow(c,2)+pow(P_x,2)*pow(c,2)+6*P_y*P_z*b*c-2*P_x*P_z*a*c+2*M_0*a*pow(b,2)-2*pow(P_z,2)*pow(b,2)+pow(P_y,2)*pow(b,2)+pow(P_x,2)*pow(b,2)
                          -2*P_x*P_y*a*b+pow(P_z,2)*pow(a,2)+pow(P_y,2)*pow(a,2);
  const G4double power_2 = pow(e,2)*i*pow(r_z,2)+4*a*e*i*pow(r_z,2)-2*pow(b,2)*i*pow(r_z,2)+pow(a,2)*i*pow(r_z,2)-e*pow(f,2)*pow(r_z,2)-4*a*pow(f,2)*pow(r_z,2)+6*b*c*f*pow(r_z,2)-4*pow(c,2)*e*pow(r_z,2)-a*pow(c,2)*pow(r_z,2)-2*e*f*i*r_y*r_z-2*b*c*i*r_y*r_z+2*pow(f,3)*r_y*r_z+2*pow(c,2)*f*r_y*r_z+2*pow(b,2)*f*r_y*r_z+2*pow(a,2)*f*r_y*r_z-2*b*c*e*r_y*r_z-2*a*b*c*r_y*r_z-2*b*f*i*r_x*r_z-2*a*c*i*r_x*r_z+2*c*pow(f,2)*r_x*r_z
                   -2*b*e*f*r_x*r_z-2*a*b*f*r_x*r_z+2*c*pow(e,2)*r_x*r_z+2*pow(c,3)*r_x*r_z+2*pow(b,2)*c*r_x*r_z-P_y*f*i*r_z-2*P_z*e*i*r_z-P_x*c*i*r_z-2*P_z*a*i*r_z+3*P_z*pow(f,2)*r_z-P_y*e*f*r_z+5*P_x*b*f*r_z-6*P_y*a*f*r_z+P_z*pow(e,2)*r_z-6*P_x*c*e*r_z+4*P_z*a*e*r_z+3*P_z*pow(c,2)*r_z+5*P_y*b*c*r_z-P_x*a*c*r_z-2*P_z*pow(b,2)*r_z+P_z*pow(a,2)*r_z
                   +e*pow(i,2)*pow(r_y,2)-pow(f,2)*i*pow(r_y,2)+4*a*e*i*pow(r_y,2)-4*pow(b,2)*i*pow(r_y,2)-4*a*pow(f,2)*pow(r_y,2)+6*b*c*f*pow(r_y,2)-2*pow(c,2)*e*pow(r_y,2)+pow(a,2)*e*pow(r_y,2)-a*pow(b,2)*pow(r_y,2)+2*b*pow(i,2)*r_x*r_y-2*c*f*i*r_x*r_y+2*b*pow(f,2)*r_x*r_y-2*c*e*f*r_x*r_y-2*a*c*f*r_x*r_y-2*a*b*e*r_x*r_y+2*b*pow(c,2)*r_x*r_y+2*pow(b,3)*r_x*r_y+P_y*pow(i,2)*r_y-P_z*f*i*r_y-2*P_y*e*i*r_y
                   -6*P_x*b*i*r_y+4*P_y*a*i*r_y+3*P_y*pow(f,2)*r_y-P_z*e*f*r_y+5*P_x*c*f*r_y-6*P_z*a*f*r_y-P_x*b*e*r_y-2*P_y*a*e*r_y-2*P_y*pow(c,2)*r_y+5*P_z*b*c*r_y+3*P_y*pow(b,2)*r_y-P_x*a*b*r_y+P_y*pow(a,2)*r_y+a*pow(i,2)*pow(r_x,2)+4*a*e*i*pow(r_x,2)-pow(c,2)*i*pow(r_x,2)-4*pow(b,2)*i*pow(r_x,2)-2*a*pow(f,2)*pow(r_x,2)+6*b*c*f*pow(r_x,2)+a*pow(e,2)*pow(r_x,2)-4*pow(c,2)*e*pow(r_x,2)
                   -pow(b,2)*e*pow(r_x,2)+P_x*pow(i,2)*r_x+4*P_x*e*i*r_x-P_z*c*i*r_x-6*P_y*b*i*r_x-2*P_x*a*i*r_x-2*P_x*pow(f,2)*r_x+5*P_y*c*f*r_x+5*P_z*b*f*r_x+P_x*pow(e,2)*r_x-6*P_z*c*e*r_x-P_y*b*e*r_x-2*P_x*a*e*r_x+3*P_x*pow(c,2)*r_x-P_z*a*c*r_x+3*P_x*pow(b,2)*r_x-P_y*a*b*r_x+M_0*pow(i,2)+4*M_0*e*i+4*M_0*a*i-2*pow(P_y,2)*i-2*pow(P_x,2)*i-2*M_0*pow(f,2)
                   +4*P_y*P_z*f+M_0*pow(e,2)+4*M_0*a*e-2*pow(P_z,2)*e-2*pow(P_x,2)*e-2*M_0*pow(c,2)+4*P_x*P_z*c-2*M_0*pow(b,2)+4*P_x*P_y*b+M_0*pow(a,2)-2*pow(P_z,2)*a-2*pow(P_y,2)*a;
  const G4double power_1 = (-2*e*i*pow(r_z,2))-2*a*i*pow(r_z,2)+2*pow(f,2)*pow(r_z,2)+2*pow(c,2)*pow(r_z,2)-4*a*f*r_y*r_z+4*b*c*r_y*r_z+4*b*f*r_x*r_z-4*c*e*r_x*r_z+P_z*i*r_z+3*P_y*f*r_z-2*P_z*e*r_z+3*P_x*c*r_z-2*P_z*a*r_z-2*e*i*pow(r_y,2)+2*pow(f,2)*pow(r_y,2)-2*a*e*pow(r_y,2)+2*pow(b,2)*pow(r_y,2)-4*b*i*r_x*r_y+4*c*f*r_x*r_y-2*P_y*i*r_y+3*P_z*f*r_y+P_y*e*r_y+3*P_x*b*r_y-2*P_y*a*r_y
                      -2*a*i*pow(r_x,2)-2*a*e*pow(r_x,2)+2*pow(c,2)*pow(r_x,2)+2*pow(b,2)*pow(r_x,2)-2*P_x*i*r_x-2*P_x*e*r_x+3*P_z*c*r_x+3*P_y*b*r_x+P_x*a*r_x-2*M_0*i-2*M_0*e-2*M_0*a+pow(P_z,2)+pow(P_y,2)+pow(P_x,2);
  const G4double power_0 = i*pow(r_z,2)+2*f*r_y*r_z+2*c*r_x*r_z+P_z*r_z+e*pow(r_y,2)+2*b*r_x*r_y+P_y*r_y+a*pow(r_x,2)+P_x*r_x+M_0;
  const std::vector<G4double> coefficients {power_6, power_5, power_4, power_3, power_2, power_1, power_0};
  // With this method, you get also get roots at singular value of the determinant, can do an if to check
  const std::vector<G4double> lambdas = polysolve.SolvePolynomial(coefficients);
  std::vector<G4ThreeVector> candidate_points;
  G4ThreeVector candidate_point;
  std::vector<G4double> candidate_distances;
  if (lambdas.size() == 0) {
    // all roots are complex/not found => lambda contour doesn't exist
    // (happens if there are infinite solutions, e.g point is at the centre of a sphere)
    // or doesn't touch/intersect the quadric surface
    return 0;
  }
  for (auto i_loop = 0; i_loop < lambdas.size(); i_loop++) {
    // Process each lambda value
    // Will run into div by zero here if determinant is equal to zero
    G4double cur_lambda = lambdas[i_loop];

    G4double x_component = ((pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e))*(r_z+cur_lambda*P_z))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))+((cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)*(r_y+cur_lambda*P_y))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))
    +(((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))*(r_x+cur_lambda*P_x))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)));
    G4double y_component = ((cur_lambda*(1-cur_lambda*a)*f+pow(cur_lambda,2)*b*c)*(r_z+cur_lambda*P_z))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))+(((1-cur_lambda*a)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(c,2))*(r_y+cur_lambda*P_y))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))
    +((cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)*(r_x+cur_lambda*P_x))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)));
    G4double z_component = (((1-cur_lambda*a)*(1-cur_lambda*e)-pow(cur_lambda,2)*pow(b,2))*(r_z+cur_lambda*P_z))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))+((cur_lambda*(1-cur_lambda*a)*f+pow(cur_lambda,2)*b*c)*(r_y+cur_lambda*P_y))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)))
    +((pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e))*(r_x+cur_lambda*P_x))/((1-cur_lambda*a)*((1-cur_lambda*e)*(1-cur_lambda*i)-pow(cur_lambda,2)*pow(f,2))-cur_lambda*b*(cur_lambda*b*(1-cur_lambda*i)+pow(cur_lambda,2)*c*f)-cur_lambda*c*(pow(cur_lambda,2)*b*f+cur_lambda*c*(1-cur_lambda*e)));

    candidate_point.set(x_component, y_component, z_component);
    candidate_points.push_back(candidate_point);
    candidate_distances.push_back((candidate_point - p).mag());
  }
  const auto min_distance = std::min_element(candidate_distances.begin(), candidate_distances.end());
  return *min_distance;
}
