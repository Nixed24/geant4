#ifndef G4HALF_SURFACE_HH
#define G4HALF_SURFACE_HH

#include "G4Types.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include <vector>


class G4HalfSurface
{
public:
    G4HalfSurface()= default;
    virtual ~G4HalfSurface() = default;

    EInside Inside(const G4ThreeVector& p) const;

    virtual G4double ImplicitFunction(const G4ThreeVector& p) const = 0;

    virtual std::vector<G4ThreeVector> RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const = 0;

    G4double Union(const G4ThreeVector& p, const G4HalfSurface& s) const;
    G4double Intersection(const G4ThreeVector& p, const G4HalfSurface& s) const;
    G4double Subtraction(const G4ThreeVector& p, const G4HalfSurface& s) const;

    virtual G4double DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const = 0;
    virtual G4double DistanceToSurface(const G4ThreeVector&p) const = 0;

    // This will involve taking a gradient of the implicit function and evaluating it at p
    virtual G4ThreeVector ClosestSurfacePoint(const G4ThreeVector& p) const = 0;
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& /*p*/) const = 0;

    virtual void Transform(const G4AffineTransform& t) = 0;

    virtual std::ostream& StreamInfo(std::ostream&) const;
};

#endif