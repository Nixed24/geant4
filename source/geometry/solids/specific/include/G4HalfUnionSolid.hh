#ifndef G4HALF_UNION_SOLID_H
#define G4HALF_UNION_SOLID_H

#include "G4VSolid.hh"
#include "G4HalfSolid.hh"
#include "G4HalfIntersectionSolid.hh"
#include "vector"

class G4HalfUnionSolid : public G4HalfSolid
{
    std::vector<std::unique_ptr<G4HalfSolid>> intersectionSolids_;
    std::vector<std::unique_ptr<G4HalfSolid>> subtractionSolids_;

protected:
    void GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const override;

public:
    using G4HalfSolid::DistanceToIn;
    using G4HalfSolid::DistanceToOut;

    G4HalfUnionSolid(std::vector<std::unique_ptr<G4HalfSolid>> intersectionSolids,
                            std::vector<std::unique_ptr<G4HalfSolid>> subtractionSolids,
                            const G4String& name = "G4HalfIntersectionSolid");

    void AddIntersectionSolid(std::unique_ptr<G4HalfIntersectionSolid> s);
    void AddSubtractionSolid(std::unique_ptr<G4HalfIntersectionSolid> s);


    // --- Inherited functions ---
    /*
     * Calculates the minimum and maximum length of the solid
     * along a given axis, clipped to pVoxelLimit.
     *
     * Used for rendering.
     *
     * Leaving this virtual for now until stewart gives me his renderer
     */
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const override {return true;}

    G4double ImplicitFunction(const G4ThreeVector &p) const override;

    void GetPotentialRayIntersections(const G4ThreeVector&, const G4ThreeVector&,
        std::vector<RayHit>&, const G4bool) const override;


    /*
     * Calculates shortest distance to the inside of the solid.
     */
    G4double DistanceToIn(const G4ThreeVector& p) const override;

    /*
     * Calculates shortest distance to the inside of the solid,
     * passing the closest surface normal by reference
     */
    G4double DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;

    G4double DistanceToOut(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const override;

    void Transform(const G4AffineTransform& t) const override;

    G4GeometryType GetEntityType() const override;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    std::ostream& StreamInfo(std::ostream& os) const override;
};

#endif
