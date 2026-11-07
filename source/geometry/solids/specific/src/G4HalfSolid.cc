#include "G4HalfSolid.hh"

G4double G4HalfSolid::Intersection(const G4ThreeVector& p, const G4HalfSolid& solid) const
{
    return std::max(this->ImplicitFunction(p), solid.ImplicitFunction(p));
}


G4double G4HalfSolid::Union(const G4ThreeVector& p, const G4HalfSolid& solid) const
{
    return std::min(this->ImplicitFunction(p), solid.ImplicitFunction(p));
}


G4double G4HalfSolid::Subtraction(const G4ThreeVector& p, const G4HalfSolid& solid) const
{
    return std::max(this->ImplicitFunction(p), -solid.ImplicitFunction(p));
}

EInside G4HalfSolid::Inside(const G4ThreeVector& p) const
{
    G4double implicitFunction = this->ImplicitFunction(p);

    return (implicitFunction > kCarTolerance) ? kOutside
          :(implicitFunction < -kCarTolerance) ? kInside
          : kSurface;
}

G4double G4HalfSolid::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    if (this->Inside(p) != kOutside) return 0;
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector.");

    G4ThreeVector vUnit = v.unit();

    std::vector<RayHit> rayHits;

    GetRayIntersections(p, vUnit, rayHits);

    for (const RayHit& hit : rayHits)
    {
        G4ThreeVector insideTestVector = hit.intersection + 2*kCarTolerance*vUnit;

        if (this->Inside(insideTestVector) == kInside) return hit.distance;
    }

    return kInfinity;
}

G4double G4HalfSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
    const G4bool calcNorm, G4bool* validNorm, G4ThreeVector* n) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude directional vector.");
    if (validNorm != nullptr) *validNorm = false;

    if (this->Inside(p) == kOutside) return 0;
    if (this->Inside(p) == kSurface)
    {
        if (calcNorm && n != nullptr) this->GetOnSurfaceNormal(p, *n);
        return 0;
    }

    const G4ThreeVector vUnit = v.unit();

    std::vector<RayHit> rayHits;

    GetRayIntersections(p, vUnit, rayHits);

    for (const RayHit& hit : rayHits)
    {
        const G4ThreeVector insideTestVector = hit.intersection + 2*kCarTolerance*vUnit;

        if (this->Inside(insideTestVector) == kOutside)
        {
            if (calcNorm && n != nullptr)
            {
                *n = (hit.isSubtraction) ? -hit.targetSurface->SurfaceNormal(hit.intersection)
                                         : hit.targetSurface->SurfaceNormal(hit.intersection);
                if (validNorm != nullptr) *validNorm = true;
            }

            return std::abs(hit.distance);
        }
    }

    return kInfinity;
}

G4ThreeVector G4HalfSolid::SurfaceNormal(const G4ThreeVector& p) const
{
    /*
     * Returns surface normal closest to a point
     */

    G4ThreeVector normal;
    if (this->Inside(p) == kInside) this->DistanceToOut(p, true, &normal);

    if (this->Inside(p) == kOutside) this->DistanceToIn(p, true, &normal);

    GetOnSurfaceNormal(p, normal);

    return normal;
}
