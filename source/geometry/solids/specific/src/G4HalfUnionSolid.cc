#include "G4HalfUnionSolid.hh"

G4HalfUnionSolid::G4HalfUnionSolid(const std::vector<G4HalfSolid*>& intersectionSolids,
                                   const std::vector<G4HalfSolid*>& subtractionSolids,
                                   const G4String& name)
                                      :  G4HalfSolid(name),
                                         intersectionSolids_(intersectionSolids),
                                         subtractionSolids_(subtractionSolids)

{
}

void G4HalfUnionSolid::AddIntersectionSolid(G4HalfIntersectionSolid& s)
{
    intersectionSolids_.push_back(&s);
}

void G4HalfUnionSolid::AddSubtractionSolid(G4HalfIntersectionSolid& s)
{
    subtractionSolids_.push_back(&s);
}

G4double G4HalfUnionSolid::ImplicitFunction(const G4ThreeVector& p) const
{
    G4double subtractionImplicitFunction = kInfinity;
    G4double intersectionImplicitFunction = kInfinity;

    for (const G4HalfSolid* s : intersectionSolids_)
    {
        intersectionImplicitFunction = std::min(intersectionImplicitFunction, s->ImplicitFunction(p));
    }

    for (const G4HalfSolid* s : subtractionSolids_)
    {
        subtractionImplicitFunction = std::min(subtractionImplicitFunction, s->ImplicitFunction(p));
    }

    G4double implicitFunction = std::max(intersectionImplicitFunction, -subtractionImplicitFunction);

    return implicitFunction;
}

void G4HalfUnionSolid::GetRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v, std::vector<RayHit>& hits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector");
    G4ThreeVector vUnit = v.unit();

    for (const G4HalfSolid* solid : intersectionSolids_) solid->GetRayIntersections(p, vUnit, hits);

    for (const G4HalfSolid* s : intersectionSolids_) s->GetRayIntersections(p, vUnit, hits);

    std::sort(hits.begin(), hits.end(),
            [](const RayHit& a, const RayHit& b)
            {
                return a.distance < b.distance;
            });
}

void G4HalfUnionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
    G4bool onSurface = this->Inside(p) == kSurface;

    for (const G4HalfSolid* solid : intersectionSolids_)
    {
        if (solid->Inside(p) == kSurface && this->Inside(p) == kSurface)
        {
            norm = solid->SurfaceNormal(p);
        }
    }
    for (const G4HalfSolid* solid : subtractionSolids_)
    {
        if (solid->Inside(p) == kSurface && this->Inside(p) == kSurface)
        {
            norm = -solid->SurfaceNormal(p);
        }
    }
}

G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    if (this->Inside(p) != kOutside) return 0;

    G4double distanceToInI = kInfinity;
    G4double distanceToInS = kInfinity;

    for (const G4HalfSolid* s : intersectionSolids_)
    {
        distanceToInI = std::min(distanceToInI, s->DistanceToIn(p));
    }
    for (const G4HalfSolid* s : subtractionSolids_)
    {
        distanceToInS = std::max(distanceToInS, s->DistanceToIn(p));
    }

    return std::max(distanceToInI, -distanceToInS);
}

G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        GetOnSurfaceNormal(p, *n);
        return 0;
    }

    G4double distanceToInI = kInfinity;
    G4double distanceToInS = kInfinity;

    G4ThreeVector normI;
    G4ThreeVector normS;

    G4ThreeVector candidateNormal;

    for (const G4HalfSolid* solid : intersectionSolids_)
    {
        G4double distance = solid->DistanceToIn(p, calcNorm,
            (calcNorm && n != nullptr) ? &candidateNormal : nullptr);

        if (distance < distanceToInI)
        {
            distanceToInI = distance;
            if (calcNorm && n != nullptr) normI = candidateNormal;
        }
    }

    for (const G4HalfSolid* solid : subtractionSolids_)
    {
       G4double distance = solid->DistanceToIn(p, calcNorm,
            (calcNorm && n != nullptr) ? &candidateNormal : nullptr);

        if (distance < distanceToInS)
        {
            distanceToInS = distance;
            candidateNormal *= -1;
            if (calcNorm && n != nullptr) normS = -candidateNormal;
        }
    }

    G4double distanceToIn = (distanceToInI > distanceToInS) ? distanceToInI : distanceToInS;
    n = (distanceToInI > distanceToInS) ? &normI : &normS;

    return distanceToIn;
}
