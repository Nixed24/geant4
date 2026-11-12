#ifndef G4HALF_INTERSECTION_SOLID_HH
#define G4HALF_INTERSECTION_SOLID_HH

#include "G4HalfSurface.hh"
#include "G4HalfSolid.hh"
#include <vector>
#include "G4VSolid.hh"

class G4HalfIntersectionSolid : public G4HalfSolid
{
    std::vector<G4HalfSurface*> intersectionSurfaces_;
    std::vector<G4HalfSurface*> subtractionSurfaces_;

protected:
    void GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& normal) const override;

public:
    using G4HalfSolid::DistanceToIn;
    using G4HalfSolid::DistanceToOut;

    G4HalfIntersectionSolid(const std::vector<G4HalfSurface*>& intersectionSurfaces,
                            const std::vector<G4HalfSurface*>& subtractionSurfaces,
                            const G4String& name = "G4HalfIntersectionSolid");

    void AddIntersectionSurface(G4HalfSurface& s);
    void AddSubtractionSurface(G4HalfSurface& s);

    G4double ImplicitFunction(const G4ThreeVector &p) const override;

    // --- Inherited functions ---
    void GetPotentialRayIntersections(const G4ThreeVector&, const G4ThreeVector&,
        std::vector<RayHit>&, const G4bool) const override;


    G4double DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;

    G4double DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    void Transform(const G4AffineTransform& t) const override;

    G4GeometryType GetEntityType() const override;
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform,
                          G4double& pMin, G4double& pMax) const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    std::ostream& StreamInfo(std::ostream& os) const override;
};

#endif