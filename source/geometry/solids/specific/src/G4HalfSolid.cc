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

void G4HalfSolid::GetSolidRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& rayHits, G4bool modifyRayHits) const
{
    auto isNotOnSurface = [&](const RayHit& hit) -> G4bool
    {
        return this->Inside(hit.intersection) != kSurface;
    };

    if (rayHits.empty() && !modifyRayHits)
        this->GetPotentialRayIntersections(p, v, rayHits, true);

    if (rayHits.empty() && modifyRayHits)
        throw std::logic_error("G4HalfSolid::GetSolidRayIntersections: rayHits is empty, but modifyRayHits was set to true.");

    if (!rayHits.empty() && modifyRayHits)
        rayHits.erase(std::remove_if(rayHits.begin(), rayHits.end(), isNotOnSurface),rayHits.end());
}

G4double G4HalfSolid::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    if (this->Inside(p) == kInside) return 0;

    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector.");

    G4ThreeVector vUnit = v.unit();

    std::vector<RayHit> rayHits;

    GetSolidRayIntersections(p, vUnit, rayHits, false);

    if (rayHits.empty()) return kInfinity;

    return rayHits[0].distance;
}

G4double G4HalfSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
    const G4bool calcNorm, G4bool* validNorm, G4ThreeVector* n) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude directional vector.");
    if (validNorm != nullptr) *validNorm = false;

    if (this->Inside(p) == kOutside) return 0;

    const G4ThreeVector vUnit = v.unit();

    std::vector<RayHit> rayHits;

    GetSolidRayIntersections(p, vUnit, rayHits, false);

    for (const RayHit& hit : rayHits)
    {
        if (this->Inside(hit.intersection) == kSurface)
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

    if (this->Inside(p) == kSurface) GetOnSurfaceNormal(p, normal);

    return normal;
}
