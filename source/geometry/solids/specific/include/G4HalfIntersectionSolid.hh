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

public:
    G4HalfIntersectionSolid(const std::vector<G4HalfSurface*>& intersectionSurfaces,
                            const std::vector<G4HalfSurface*>& subtractionSurfaces,
                            const G4String& name = "G4HalfIntersectionSolid");

    void AddIntersectionSurface(G4HalfSurface& s);
    void AddSubtractionSurface(G4HalfSurface& s);

    void GetRayIntersections(const G4ThreeVector&, const G4ThreeVector&,
        std::vector<RayHit>&) const override;

    G4double ImplicitFunction(const G4ThreeVector &p) const override;

    // Provides a new implicit function for a subtraction with a solid
    G4double Subtraction(const G4ThreeVector& p, const G4HalfIntersectionSolid& subtraction) const;


    // --- Inherited functions ---
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const override = 0;

    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p, G4ThreeVector* n) const override;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v, const G4bool calcNorm, G4bool* validNorm,
                           G4ThreeVector* n) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, G4ThreeVector* n) const override;

    G4GeometryType GetEntityType() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    std::ostream& StreamInfo(std::ostream& os) const override;
};

#endif