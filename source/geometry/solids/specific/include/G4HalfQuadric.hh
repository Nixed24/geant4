#ifndef G4HALF_QUADRIC_HH
#define G4HALF_QUADRIC_HH

#include <vector>

#include "G4HalfSurface.hh"
#include "G4Types.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
class G4HalfQuadric : public G4HalfSurface {
  private:
    G4RotationMatrix Q;
    G4ThreeVector P;
    G4double R;
  public:
    G4HalfQuadric() = default;
    ~G4HalfQuadric() override = default;

    G4double ImplicitFunction(const G4ThreeVector& p) const override;

    G4double Union(const G4ThreeVector& p, const G4HalfSurface& s) const override;
    std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double Intersection(const G4ThreeVector& p, const G4HalfSurface& s) const override;
    G4double Subtraction(const G4ThreeVector& p, const G4HalfSurface& s) const override;

    G4double DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double DistanceToSurface(const G4ThreeVector&p) const override;

    // This will involve taking a gradient of the implicit function and evaluating it at p
    G4ThreeVector SurfaceNormal(const G4ThreeVector& /*p*/) const override;

    std::ostream& StreamInfo(std::ostream&) const override;
};

#endif