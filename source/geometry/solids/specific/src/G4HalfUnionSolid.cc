#include "G4HalfUnionSolid.hh"

#include "G4VGraphicsScene.hh"

G4HalfUnionSolid::G4HalfUnionSolid(std::vector<std::unique_ptr<G4HalfSolid>> intersectionSolids,
                                   std::vector<std::unique_ptr<G4HalfSolid>> subtractionSolids,
                                   const G4String& name)
                                      :  G4HalfSolid(name),
                                         intersectionSolids_(std::move(intersectionSolids)),
                                         subtractionSolids_(std::move(subtractionSolids))

{
}

void G4HalfUnionSolid::AddIntersectionSolid(std::unique_ptr<G4HalfIntersectionSolid> s)
{
    intersectionSolids_.push_back(std::move(s));
}

void G4HalfUnionSolid::AddSubtractionSolid(std::unique_ptr<G4HalfIntersectionSolid> s)
{
    subtractionSolids_.push_back(std::move(s));
}

G4double G4HalfUnionSolid::ImplicitFunction(const G4ThreeVector& p) const
{
    G4double subtractionImplicitFunction = kInfinity;
    G4double intersectionImplicitFunction = kInfinity;

    for (const auto &s : intersectionSolids_)
    {
        intersectionImplicitFunction = std::min(intersectionImplicitFunction, s->ImplicitFunction(p));
    }

    for (const auto &s : subtractionSolids_)
    {
        subtractionImplicitFunction = std::min(subtractionImplicitFunction, s->ImplicitFunction(p));
    }

    const G4double implicitFunction = std::max(intersectionImplicitFunction, -subtractionImplicitFunction);

    return implicitFunction;
}

void G4HalfUnionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits, const G4bool keepOnlySolidHits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude unit vector");
    const G4ThreeVector vUnit = v.unit();

    for (const auto &solid : intersectionSolids_)
        solid->GetPotentialRayIntersections(p, vUnit, hits, keepOnlySolidHits);

    for (const auto &solid : subtractionSolids_)
        solid->GetPotentialRayIntersections(p, vUnit, hits, keepOnlySolidHits);

    std::sort(hits.begin(), hits.end(),
            [](const RayHit& a, const RayHit& b)
            {
                return a.distance < b.distance;
            });
}

void G4HalfUnionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
    const G4bool onSolid = this->Inside(p) == kSurface;

    if (not onSolid)
        throw std::logic_error("G4HalfIntersectionSolid::GetOnSurfaceNormal: Point passed in is not on the surface.");

    for (const auto &solid : intersectionSolids_)
    {
        if (solid->Inside(p) == kSurface)
            norm = solid->SurfaceNormal(p);
    }
    for (const auto &solid : subtractionSolids_)
    {
        if (solid->Inside(p) == kSurface)
            norm = -solid->SurfaceNormal(p);
    }
}

G4double G4HalfUnionSolid::SignedDistance(const G4ThreeVector& p, G4ThreeVector* n, const G4bool calcNorm) const
{
    struct CandidateSolid
    {
        G4double signedDistance = kInfinity;
        G4ThreeVector norm;
        G4int normSign;
    };

    auto setCandidate = [&](const std::vector<std::unique_ptr<G4HalfSolid>>& solidList, const G4int sign, CandidateSolid& best)
    {
        for (const auto &solid : solidList)
        {
            CandidateSolid candidate;
            candidate.normSign = sign;
            candidate.signedDistance = solid->SignedDistance(p, &candidate.norm, calcNorm);

            if (candidate.signedDistance < best.signedDistance)
            {
                best = candidate;
            }
        }
    };

    CandidateSolid intersectionCandidate, subtractionCandidate;

    setCandidate(intersectionSolids_, +1, intersectionCandidate);
    setCandidate(subtractionSolids_, -1, subtractionCandidate);

    if (calcNorm && n != nullptr)
    {
        *n = (intersectionCandidate.signedDistance > -subtractionCandidate.signedDistance)
        ? intersectionCandidate.norm : subtractionCandidate.norm;
    }

    return std::max(intersectionCandidate.signedDistance, -subtractionCandidate.signedDistance);
}

G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    return std::min(0.0, this->SignedDistance(p, n, calcNorm));
}


G4double G4HalfUnionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    return DistanceToIn(p, false, nullptr);
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
   return DistanceToOut(p, calcNorm, n);
}


G4double G4HalfUnionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    return DistanceToOut(p, false, nullptr);
}

void G4HalfUnionSolid::Transform(const G4AffineTransform& t) const
{
    for (auto &solid : intersectionSolids_) solid->Transform(t);

    for (auto &solid : subtractionSolids_) solid->Transform(t);
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
