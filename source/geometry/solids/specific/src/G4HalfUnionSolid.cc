#include "G4HalfUnionSolid.hh"

#include "G4VGraphicsScene.hh"

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

void G4HalfUnionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v, std::vector<RayHit>& hits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector");
    G4ThreeVector vUnit = v.unit();

    for (const G4HalfSolid* solid : intersectionSolids_) solid->GetPotentialRayIntersections(p, vUnit, hits);

    for (const G4HalfSolid* s : subtractionSolids_) s->GetPotentialRayIntersections(p, vUnit, hits);

    std::sort(hits.begin(), hits.end(),
            [](const RayHit& a, const RayHit& b)
            {
                return a.distance < b.distance;
            });
}

void G4HalfUnionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
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

G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    struct CandidateSolid
    {
        G4double distance = kInfinity;
        const G4HalfSolid* solid = nullptr;
        G4int normSign = +1;
    };

    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        GetOnSurfaceNormal(p, *n);
        return 0;
    }

    auto setCandidate = [&](const std::vector<G4HalfSolid*>& solidList, const G4int sign, CandidateSolid& best)
    {
        for (G4HalfSolid* solid : solidList)
        {
            G4double candidateDistance = solid->DistanceToIn(p);

            if (candidateDistance < best.distance)
            {
                best = {candidateDistance, solid, sign};
            }
        }
    };

    CandidateSolid intersectionCandidate, subtractionCandidate;

    setCandidate(intersectionSolids_, +1, intersectionCandidate);
    setCandidate(subtractionSolids_, -1, subtractionCandidate);

    if (calcNorm && n != nullptr)
    {
        *n = (intersectionCandidate.distance > -subtractionCandidate.distance)
            ? intersectionCandidate.solid->SurfaceNormal(p)
            : -subtractionCandidate.solid->SurfaceNormal(p);
    }

    return std::max(intersectionCandidate.distance, -subtractionCandidate.distance);
}


G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    return DistanceToIn(p, false, nullptr);
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    struct CandidateSolid
    {
        G4double distance = kInfinity;
        const G4HalfSolid* solid = nullptr;
        G4int normSign = +1;
    };

    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        GetOnSurfaceNormal(p, *n);
        return 0;
    }


    auto setCandidate = [&](const std::vector<G4HalfSolid*>& solidList, const G4int sign, CandidateSolid& best)
    {
        for (G4HalfSolid* solid : solidList)
        {
            G4double candidateDistance = solid->DistanceToIn(p);

            if (candidateDistance < best.distance)
            {
                best = {candidateDistance, solid, sign};
            }
        }
    };

    CandidateSolid intersectionCandidate, subtractionCandidate;

    setCandidate(intersectionSolids_, +1, intersectionCandidate);
    setCandidate(subtractionSolids_, -1, subtractionCandidate);

    if (calcNorm && n != nullptr)
    {
        *n = (intersectionCandidate.distance > -subtractionCandidate.distance)
            ? intersectionCandidate.solid->SurfaceNormal(p)
            : -subtractionCandidate.solid->SurfaceNormal(p);
    }

    return std::max(intersectionCandidate.distance, -subtractionCandidate.distance);
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    return DistanceToOut(p, false, nullptr);
}

void G4HalfUnionSolid::Transform(const G4AffineTransform& t) const
{
    for (G4HalfSolid* solid : intersectionSolids_) solid->Transform(t);

    for (G4HalfSolid* solid : subtractionSolids_) solid->Transform(t);
}


G4GeometryType G4HalfUnionSolid::GetEntityType() const
{
    return "G4HalfUnionSolid";
}

void G4HalfUnionSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
    scene.AddSolid(*this);
}

std::ostream& G4HalfUnionSolid::StreamInfo(std::ostream& os) const
{
    os << "--- Half Union Solid: " << GetName() << " ---" << std::endl;

    return os;
}
