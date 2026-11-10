#include "G4HalfIntersectionSolid.hh"
#include <algorithm>

#include "G4HalfPlane.hh"
#include "G4VGraphicsScene.hh"
#include "G4VoxelLimits.hh"

const G4ThreeVector ORIGIN = {0,0,0};

G4HalfIntersectionSolid::G4HalfIntersectionSolid(const std::vector<G4HalfSurface*>& intersectionSurfaces,
                                                 const std::vector<G4HalfSurface*>& subtractionSurfaces,
                                                 const G4String& name)
                                                     :  G4HalfSolid(name),
                                                        intersectionSurfaces_(intersectionSurfaces),
                                                        subtractionSurfaces_(subtractionSurfaces)

{
}


void G4HalfIntersectionSolid::AddIntersectionSurface(G4HalfSurface& s)
{
    intersectionSurfaces_.push_back(&s);
}


void G4HalfIntersectionSolid::AddSubtractionSurface(G4HalfSurface& s)
{
    subtractionSurfaces_.push_back(&s);
}


G4double G4HalfIntersectionSolid::ImplicitFunction(const G4ThreeVector& p) const
{
    G4double implicitFunction = -kInfinity;

    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, s->ImplicitFunction(p));
    }

    for (const G4HalfSurface* s : subtractionSurfaces_)
    {
        implicitFunction = std::max(implicitFunction, -s->ImplicitFunction(p));
    }

    return implicitFunction;
}


void G4HalfIntersectionSolid::GetPotentialRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits) const
{
    if (v.mag() == 0) throw std::logic_error("Cannot pass zero magnitude directional vector.");

    G4ThreeVector vUnit = v.unit();

    auto getHits = [&p, &vUnit](const G4HalfSurface* s, const G4bool subtraction)
    {
        std::vector<RayHit> surfaceHits;

        for (const G4ThreeVector& intersection : s->RayIntersection(p,vUnit))
        {
            RayHit hit;
            hit.intersection = intersection;
            hit.distance = (intersection - p).dot(vUnit);
            hit.targetSurface = s;
            hit.isSubtraction = subtraction;

            surfaceHits.push_back(hit);
        }

        return surfaceHits;
    };

    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        std::vector<RayHit> intersectionHits = getHits(s, false);
        hits.insert(hits.end(), intersectionHits.begin(), intersectionHits.end());
    }

    for (const G4HalfSurface* s : subtractionSurfaces_)
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
        GetOnSurfaceNormal(p, *n);
        return 0;
    }

    G4double distanceToIn = -kInfinity;
    const G4HalfSurface* candidateSurface = nullptr;
    G4int normSign;

    auto setDistance = [&](const G4HalfSurface* surface, const G4int sign)
    {
        G4double candidateDistance = surface->DistanceToSurface(p);
        if (sign*candidateDistance > distanceToIn)
        {
            distanceToIn = candidateDistance;
            candidateSurface = surface;
            normSign = sign;
        }
    };

    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
        setDistance(surface, +1);
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        setDistance(surface, -1);
    }

    if (calcNorm &&  n != nullptr)
    {
        G4ThreeVector surfacePoint = candidateSurface->ClosestSurfacePoint(p);
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

    auto setDistance = [&](const G4HalfSurface* surface, const G4int sign)
    {
        G4double candidateDistance = surface->DistanceToSurface(p);
        if (sign*candidateDistance > distanceToOut)
        {
            distanceToOut = candidateDistance;
            candidateSurface = surface;
            normSign = sign;
        }
    };

    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
        setDistance(surface, +1);
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        setDistance(surface, -1);
    }

    if (calcNorm &&  n != nullptr)
    {
        G4ThreeVector surfacePoint = candidateSurface->ClosestSurfacePoint(p);
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
    G4bool onSurface = this->Inside(p) == kSurface;

    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        if (s->Inside(p) == kSurface && this->Inside(p) == kSurface)
        {
            norm = (onSurface) ? s->SurfaceNormal(p) : s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
    for (const G4HalfSurface* s : subtractionSurfaces_)
    {
        if (s->Inside(p) == kSurface && this->Inside(p) == kSurface)
        {
            norm = (onSurface) ? -s->SurfaceNormal(p) : -s->SurfaceNormal(s->ClosestSurfacePoint(p));
        }
    }
}


void G4HalfIntersectionSolid::Transform(const G4AffineTransform& t) const
{
    for (G4HalfSurface* s : intersectionSurfaces_) s->Transform(t);

    for (G4HalfSurface* s : subtractionSurfaces_) s->Transform(t);
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

/*G4bool G4HalfIntersectionSolid::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax) const
{
    for (G4HalfSurface* s : intersectionSurfaces_) s->Transform(pTransform);
    for (G4HalfSurface* s : subtractionSurfaces_) s->Transform(pTransform);

    if (pAxis == kUndefined)
    {
        throw std::invalid_argument("G4HalfIntersectionSolid::CalculateExtent: pAxis is undefined.");
    }

    if (pAxis > 2)
    {
        throw std::invalid_argument("G4HalfIntersectionSolid::CalculateExtent: invalid pAxis.");
    }

    G4HalfPlane voxelLimitPlanePosX({+1,0,0}, pVoxelLimit.GetMaxXExtent());
    G4HalfPlane voxelLimitPlaneNegX({-1,0,0}, pVoxelLimit.GetMinXExtent());
    G4HalfPlane voxelLimitPlanePosY({0,+1,0}, pVoxelLimit.GetMaxYExtent());
    G4HalfPlane voxelLimitPlaneNegY({0,-1,0}, pVoxelLimit.GetMinYExtent());
    G4HalfPlane voxelLimitPlanePosZ({0,0,+1}, pVoxelLimit.GetMaxZExtent());
    G4HalfPlane voxelLimitPlaneNegZ({0,0,-1}, pVoxelLimit.GetMinZExtent());

    std::vector<G4HalfSurface*> voxelLimits =
        {
        &voxelLimitPlanePosX, &voxelLimitPlanePosY, &voxelLimitPlanePosZ,
        &voxelLimitPlaneNegX, &voxelLimitPlaneNegY, &voxelLimitPlaneNegZ
        };

    std::vector<G4HalfSurface*> nA;

    G4HalfIntersectionSolid voxelLimitCube(voxelLimits, nA, "VoxelLimitCube");

    G4ThreeVector axisRayVPos = ORIGIN;
    G4ThreeVector axisRayVNeg = ORIGIN;

    axisRayVPos[pAxis] = 1;
    axisRayVNeg[pAxis] = -1;

    std::vector<G4ThreeVector> axisIntersectionsPos;
    std::vector<G4ThreeVector> axisIntersectionsNeg;

    // Obtain intersections through all surfaces
    for (G4HalfSurface* s: intersectionSurfaces_)
    {
        std::vector<G4ThreeVector> surfaceIntersectionsPos = s->RayIntersection(ORIGIN, axisRayVPos);
        std::vector<G4ThreeVector> surfaceIntersectionsNeg = s->RayIntersection(ORIGIN, axisRayVNeg);

        axisIntersectionsPos.insert(axisIntersectionsPos.begin(),
            surfaceIntersectionsPos.begin(),surfaceIntersectionsPos.end());
        axisIntersectionsNeg.insert(axisIntersectionsNeg.begin(),
            surfaceIntersectionsNeg.begin(),surfaceIntersectionsNeg.end());
    }
    for (G4HalfSurface* s: subtractionSurfaces_)
    {
        std::vector<G4ThreeVector> surfaceIntersectionsPos = s->RayIntersection(ORIGIN, axisRayVPos);
        std::vector<G4ThreeVector> surfaceIntersectionsNeg = s->RayIntersection(ORIGIN, axisRayVNeg);

        axisIntersectionsPos.insert(axisIntersectionsPos.begin(),
            surfaceIntersectionsPos.begin(),surfaceIntersectionsPos.end());
        axisIntersectionsNeg.insert(axisIntersectionsNeg.begin(),
            surfaceIntersectionsNeg.begin(),surfaceIntersectionsNeg.end());
    }

    if (axisIntersectionsNeg.empty() && axisIntersectionsPos.empty()) return false;
    

    std::vector<G4ThreeVector> validIntersections;
    for (G4ThreeVector& intersectionPoint : axisIntersectionsPos)
    {
        if (voxelLimitCube.Inside(intersectionPoint) == kInside) continue;
        if (this->Inside(intersectionPoint) != kSurface) continue;
        validIntersections.push_back(intersectionPoint);
    }

    if (validIntersections.empty()) return false;

    pMin = pVoxelLimit.GetMinExtent(pAxis);
    pMax = pVoxelLimit.GetMaxExtent(pAxis);

    for (G4ThreeVector intersectionPoint : validIntersections)
    {
        pMin = std::min(pMin,intersectionPoint[pAxis]);
        pMax = std::max(pMax,intersectionPoint[pAxis]);
    }

    return true;
}*/

std::ostream& G4HalfIntersectionSolid::StreamInfo(std::ostream& os) const
{
    os << "--- Half Intersection Solid: " << GetName() << " ---" << std::endl;
    os << "Surface Intersections: " << std::endl;
    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        os << s << std::endl;
    }
    os << "Surface Subtractions: " << std::endl;
    for (const G4HalfSurface* s : subtractionSurfaces_)
    {
        os << s << std::endl;
    }
    return os;
}