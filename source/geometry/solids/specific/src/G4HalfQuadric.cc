#include "G4HalfQuadric.hh"
#include <algorithm>

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
+ this->Q.zx() * v.z() * p.x() + this->Q.zy() * v.z() * p.y() + this->Q.zz() * v.z() * p.z();
  + (this->P * v);
  G4double c = this->Q.xx() * p.x() * p.x() + this->Q.xy() * p.x() * p.y() + this->Q.xz() * p.x() * p.z()
           + this->Q.yx() * p.y() * p.x() + this->Q.yy() * p.y() * p.y() + this->Q.yz() * p.y() * p.z()
           + this->Q.zx() * p.z() * p.x() + this->Q.zy() * p.z() * p.y() + this->Q.zz() * p.z() * p.z();
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
G4double G4HalfQuadric::DistanceToSurface(const G4ThreeVector& p) const
{
  
  return 0;
}
std::ostream& G4HalfQuadric::StreamInfo(std::ostream& os) const
{
    os << "--- Half Surface Quadric ---" << std::endl;
    return os;
}