#ifndef G4HALF_QUADRIC_HH
#define G4HALF_QUADRIC_HH

#include <vector>

#include "G4HalfSurface.hh"
#include "G4Types.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
class G4HalfQuadric : public G4HalfSurface {
  private:

    static G4RotationMatrix _T(const G4RotationMatrix & M);
    std::vector<G4ThreeVector> GetLMContourIntersections(const G4ThreeVector& p) const;

    G4RotationMatrix Q; // Default constructor of G4RotationMatrix gives 3x3 identity
    G4ThreeVector P;
    G4double R {-1}; // R = -1 gives a unit sphere at the origin
  public:

    G4RotationMatrix get_Q() const;
    G4ThreeVector get_P() const;
    G4double get_R() const;
    void set_Q(const G4RotationMatrix& new_Q);
    void set_P(const G4ThreeVector& new_P);
    void set_R(G4double new_R);

    G4HalfQuadric() = default;
    ~G4HalfQuadric() override = default;

    G4double ImplicitFunction(const G4ThreeVector& p) const override;

    std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

    G4double DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double DistanceToSurface(const G4ThreeVector&p) const override;

    // This will involve taking a gradient of the implicit function and evaluating it at p
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4ThreeVector ClosestSurfacePoint(const G4ThreeVector& p) const override;

    void Transform(const G4AffineTransform& t) override;

    std::ostream& StreamInfo(std::ostream&) const override;
};

#endif