#ifndef GEANT4_G4HALFSOLID_HH
#define GEANT4_G4HALFSOLID_HH

#include "G4VSolid.hh"
#include "G4HalfSurface.hh"

class G4HalfSolid : public G4VSolid
{
protected:
    struct RayHit
    {
        const G4HalfSurface* targetSurface;
        G4bool isSubtraction;
        G4ThreeVector intersection;
        G4double distance;
    };

public:
    G4HalfSolid(const G4String& name) : G4VSolid(name) {};

    virtual G4double ImplicitFunction(const G4ThreeVector&) const = 0;

    // Provides boolean operations on other half solids.
    G4double Intersection(const G4ThreeVector&, const G4HalfSolid&) const;
    G4double Union(const G4ThreeVector&, const G4HalfSolid&) const;
    G4double Subtraction(const G4ThreeVector&, const G4HalfSolid&) const;

    // Populates a vector of RayHit structs with data about where the ray intersected
    virtual void GetPotentialRayIntersections(const G4ThreeVector&, const G4ThreeVector&, std::vector<RayHit>&) const = 0;
    virtual void GetSolidRayIntersections(const G4ThreeVector&, const G4ThreeVector&, std::vector<RayHit>&, G4bool) const;

    EInside Inside(const G4ThreeVector&) const override;

    using G4VSolid::DistanceToOut;
    using G4VSolid::DistanceToIn;

    // Calculates shortest distance from a point, passes normal by reference
    virtual G4double DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const = 0;
    virtual G4double DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const = 0;

    virtual void GetOnSurfaceNormal(const G4ThreeVector&, G4ThreeVector&) const = 0;

    // Can be implemented in the base class if the derived classes have a defined GetRayIntersection
    G4double DistanceToIn(const G4ThreeVector&, const G4ThreeVector&) const override;
    G4double DistanceToOut(const G4ThreeVector&, const G4ThreeVector&, const G4bool, G4bool*, G4ThreeVector*) const override;

    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;

    virtual void Transform(const G4AffineTransform& t) const = 0;
};

#endif //GEANT4_G4HALFSOLID_HH
