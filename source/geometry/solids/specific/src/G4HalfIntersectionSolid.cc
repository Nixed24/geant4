#include "G4HalfIntersectionSolid.hh"
#include <algorithm>

#include "G4HalfPlane.hh"
#include "G4VGraphicsScene.hh"
#include "G4VoxelLimits.hh"

const G4ThreeVector ORIGIN = {0,0,0};

G4HalfIntersectionSolid::G4HalfIntersectionSolid(std::vector<std::unique_ptr<G4HalfSurface>> intersectionSurfaces,
                                                 std::vector<std::unique_ptr<G4HalfSurface>> subtractionSurfaces,
                                                 const G4String& name)
                                                     :  G4HalfSolid(name),
                                                        intersectionSurfaces_(std::move(intersectionSurfaces)),
                                                        subtractionSurfaces_(std::move(subtractionSurfaces))

{
}


void G4HalfIntersectionSolid::AddIntersectionSurface(std::unique_ptr<G4HalfSurface> s)
{
    intersectionSurfaces_.push_back(std::move(s));
}


void G4HalfIntersectionSolid::AddSubtractionSurface(std::unique_ptr<G4HalfSurface> s)
{
    subtractionSurfaces_.push_back(std::move(s));
}

G4double G4HalfIntersectionSolid::ImplicitFunction(const G4ThreeVector &p) const
{
    G4double implicitFunction = -kInfinity;

    for (const auto &s : intersectionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, s->ImplicitFunction(p));
    }

    for (const auto &s : subtractionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, -s->ImplicitFunction(p));
    }

    return implicitFunction;
}


void G4HalfIntersectionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits, const G4bool keepOnlySolidHits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude directional vector.");

    G4ThreeVector vUnit = v.unit();

    auto getHits = [&](const std::unique_ptr<G4HalfSurface> &s, const G4bool subtraction)
    {
        std::vector<RayHit> surfaceHits;

        for (const G4ThreeVector& intersection : s->RayIntersection(p,vUnit))
        {
            if (keepOnlySolidHits && this->Inside(intersection) != kSurface) continue;

            RayHit hit;
            hit.intersection = intersection;
            hit.distance = (intersection - p).dot(vUnit);
            hit.targetSurface = s.get();
            hit.isSubtraction = subtraction;

            surfaceHits.push_back(hit);
        }

        return surfaceHits;
    };

    for (const auto &s : intersectionSurfaces_)
    {
        std::vector<RayHit> intersectionHits = getHits(s, false);
        hits.insert(hits.end(), intersectionHits.begin(), intersectionHits.end());
    }

    for (const auto &s : subtractionSurfaces_)
    {
        std::vector<RayHit> subtractionHits = getHits(s, true);
        hits.insert(hits.end(), subtractionHits.begin(), subtractionHits.end());
    }

    std::sort(hits.begin(), hits.end(),
        [](const RayHit& a, const RayHit& b)
        {
            return a.distance < b.distance;
        });
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        if (calcNorm && n != nullptr) GetOnSurfaceNormal(p, *n);
        return 0;
    }

    G4double distanceToIn = -kInfinity;
    const G4HalfSurface* candidateSurface = nullptr;
    G4int normSign;

    auto setDistance = [&](const std::unique_ptr<G4HalfSurface> &surface, const G4int sign)
    {
        G4double candidateDistance = surface->DistanceToSurface(p);
        if (sign*candidateDistance > distanceToIn)
        {
            distanceToIn = candidateDistance;
            candidateSurface = surface.get();
            normSign = sign;
        }
    };

    for (const auto &surface : intersectionSurfaces_)
    {
        setDistance(surface, +1);
    }

    for (const auto &surface : subtractionSurfaces_)
    {
        setDistance(surface, -1);
    }

    if (calcNorm &&  n != nullptr)
    {
        const G4ThreeVector surfacePoint = candidateSurface->ClosestSurfacePoint(p);
        *n = normSign*candidateSurface->SurfaceNormal(surfacePoint);
    }

    return distanceToIn;
}


G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    return DistanceToIn(p, false, nullptr);
}


G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p, const G4bool calcNorm, G4ThreeVector* n) const
{
    if (this->Inside(p) == kOutside) return 0;
    if (this->Inside(p) == kSurface)
    {
        if (calcNorm && n != nullptr) GetOnSurfaceNormal(p, *n);
        return 0;
    }

    G4double distanceToOut = -kInfinity;
    const G4HalfSurface* candidateSurface = nullptr;
    G4int normSign;

    auto setDistance = [&](const std::unique_ptr<G4HalfSurface> &surface, const G4int sign)
    {
        G4double candidateDistance = surface->DistanceToSurface(p);
        if (sign*candidateDistance > distanceToOut)
        {
            distanceToOut = candidateDistance;
            candidateSurface = surface.get();
            normSign = sign;
        }
    };

    for (const auto &surface : intersectionSurfaces_)
    {
        setDistance(surface, +1);
    }

    for (const auto &surface : subtractionSurfaces_)
    {
        setDistance(surface, -1);
    }

    if (calcNorm &&  n != nullptr)
    {
        const G4ThreeVector surfacePoint = candidateSurface->ClosestSurfacePoint(p);
        *n = normSign*candidateSurface->SurfaceNormal(surfacePoint);
    }

    return std::abs(distanceToOut);
}


G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    return DistanceToOut(p, false, nullptr);
}


void G4HalfIntersectionSolid::GetOnSurfaceNormal(const G4ThreeVector& p, G4ThreeVector& norm) const
{
    const G4bool onSolid = this->Inside(p) == kSurface;
    if (not onSolid)
        throw std::logic_error("G4HalfIntersectionSolid::GetOnSurfaceNormal: Point passed in is not on the surface.");

    for (const auto &s : intersectionSurfaces_)
    {
        if (s->Inside(p) == kSurface)
        {
            norm = (onSolid) ? s->SurfaceNormal(p) : s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
    for (const auto &s : subtractionSurfaces_)
    {
        if (s->Inside(p) == kSurface)
        {
            norm = (onSolid) ? -s->SurfaceNormal(p) : -s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
}


void G4HalfIntersectionSolid::Transform(const G4AffineTransform& t) const
{
    for (auto &s : intersectionSurfaces_) s->Transform(t);

    for (auto &s : subtractionSurfaces_) s->Transform(t);
}


G4GeometryType G4HalfIntersectionSolid::GetEntityType() const
{
    return "G4HalfIntersectionSolid";
}

void G4HalfIntersectionSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
    // Standard implementation as suggested by CoPilot
    scene.AddSolid(*this);
}

G4bool G4HalfIntersectionSolid::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const
{
    return true;
}

std::ostream& G4HalfIntersectionSolid::StreamInfo(std::ostream& os) const
{
    os << "--- Half Intersection Solid: " << GetName() << " ---" << std::endl;
    os << "Surface Intersections: " << std::endl;
    for (const auto &s : intersectionSurfaces_)
    {
        s->StreamInfo(os);
    }
    os << "Surface Subtractions: " << std::endl;
    for (const auto &s : subtractionSurfaces_)
    {
        s->StreamInfo(os);
    }
    return os;
}