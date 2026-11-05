#include "G4HalfIntersectionSolid.hh"
#include <algorithm>

#include "G4HalfPlane.hh"
#include "G4VGraphicsScene.hh"
#include "G4VoxelLimits.hh"

const G4ThreeVector ORIGIN = {0,0,0};

G4HalfIntersectionSolid::G4HalfIntersectionSolid(const std::vector<G4HalfSurface*>& intersectionSurfaces,
                                                 const std::vector<G4HalfSurface*>& subtractionSurfaces,
                                                 const G4String& name)
                                                     :  G4VSolid(name),
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

G4double G4HalfIntersectionSolid::Subtraction(const G4ThreeVector& p, const G4HalfIntersectionSolid& subtraction) const
{
    G4double implicitFunction = this->ImplicitFunction(p);

    implicitFunction = std::max(implicitFunction, -subtraction.ImplicitFunction(p));

    return implicitFunction;
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

void G4HalfIntersectionSolid::GetRayIntersections(const G4ThreeVector& p, const G4ThreeVector& v,
    std::vector<RayHit>& hits) const
{
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
        hits.insert(hits.begin(), intersectionHits.begin(), intersectionHits.end());
    }

   for (const G4HalfSurface* s : subtractionSurfaces_)
   {
        std::vector<RayHit> subtractionHits = getHits(s, true);
        hits.insert(hits.begin(), subtractionHits.begin(), subtractionHits.end());
   }
    std::sort(hits.begin(), hits.end(),
        [](const RayHit& a, const RayHit& b)
        {
            return a.distance < b.distance;
        });
}

G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p) const
{
    if (this->Inside(p) != kInside) return 0;

    G4double distanceToOut = -kInfinity;

    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        distanceToOut = std::max(distanceToOut, s->DistanceToSurface(p));
    }
    for (const G4HalfSurface* s : subtractionSurfaces_)
    {
        distanceToOut = std::max(distanceToOut, -s->DistanceToSurface(p));
    }

    return distanceToOut;
}
G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p, const G4HalfSurface*& s) const
{
    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        s = &this->WhatSurfaceAmIOn(p);
        return 0;
    }

    G4double distanceToOut = -kInfinity;

    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
       if (s->DistanceToSurface(p) > distanceToOut)
       {
           distanceToOut = s->DistanceToSurface(p);
           s = surface;
       }
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        if (-s->DistanceToSurface(p) > distanceToOut)
        {
            distanceToOut = -s->DistanceToSurface(p);
            s = surface;
        }
    }

    return distanceToOut;
}

G4double G4HalfIntersectionSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v, const G4bool calcNorm, G4bool* validNorm, G4ThreeVector* n) const
{
    G4ThreeVector vUnit = v.unit();
    if (this->Inside(p) != kInside) return 0;
    G4double distanceToOut = kInfinity;
    std::vector<RayHits> rayHits;

    GetRayIntersections(p, v, rayHits);

    for (const RayHits& hit : rayHits)
    {

    }

    return distanceToOut;
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p) const
{
    if (this -> Inside(p) == kInside) return 0;

    G4double distanceToIn = -kInfinity;

    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
        distanceToIn = std::max(distanceToIn, surface->DistanceToSurface(p));
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        distanceToIn = std::max(distanceToIn, -surface->DistanceToSurface(p));
    }

    return distanceToIn;
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p, const G4HalfSurface*& s) const
{
    if (this->Inside(p) == kInside) return 0;
    if (this->Inside(p) == kSurface)
    {
        s = &this->WhatSurfaceAmIOn(p);
        return 0;
    }

    G4double distanceToIn = -kInfinity;
    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
        G4double distanceToSurface = surface->DistanceToSurface(p);
        if (distanceToSurface > distanceToIn)
        {
            distanceToIn = distanceToSurface;
            s = surface;
        }
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        G4double distanceToSurface = -surface->DistanceToSurface(p);
        if (distanceToSurface > distanceToIn)
        {
            distanceToIn = distanceToSurface;
            s = surface;
        }
    }

    return distanceToIn;
}

G4double G4HalfIntersectionSolid::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4ThreeVector vUnit = v.unit();
    if (this -> Inside(p) != kOutside) return 0;

    G4double distanceToIn = kInfinity;
    for (const G4HalfSurface* surface : intersectionSurfaces_)
    {
        G4double distance = surface->DistanceToSurface(p,v);
        if (distance == kInfinity) continue;
        G4ThreeVector point = p + vUnit * (distance);

        if (this->Inside(point + vUnit*kCarTolerance) == kInside)
        {
            distanceToIn = std::min(distanceToIn, distance);
        }
    }

    for (const G4HalfSurface* surface : subtractionSurfaces_)
    {
        G4double distance = surface->DistanceToSurface(p,v);
        if (distance == kInfinity) continue;
        G4ThreeVector point = p + vUnit * (distance);

        if (this->Inside(point + vUnit*kCarTolerance) == kInside)
        {
            distanceToIn = std::min(distanceToIn, distance);
        }
    }

    return distanceToIn;
}


const G4HalfSurface& G4HalfIntersectionSolid::WhatSurfaceAmIOn(const G4ThreeVector& p) const
{
    if (this->Inside(p) != kSurface)
    {
        throw std::invalid_argument("G4HalfIntersectionSolid::WhatSurfaceAmIOn takes surface points as argument only.");
    }

    for (const G4HalfSurface* s : intersectionSurfaces_)
    {
        if (s->Inside(p) == kSurface) return *s;
    }
    for (const G4HalfSurface* s : subtractionSurfaces_)
    {
        if (s->Inside(p) == kSurface) return *s;
    }

    throw std::logic_error("G4HalfIntersectionSolid has no surfaces.");
}

G4ThreeVector G4HalfIntersectionSolid::SurfaceNormal(const G4ThreeVector& p) const
{
    /*
     * Returns surface normal closest to a point
     */

    const G4HalfSurface* surface = nullptr;
    if (this->Inside(p) == kInside)
    {
        this->DistanceToOut(p, surface);
        return surface->SurfaceNormal(p);
    }
    if (this->Inside(p) == kOutside)
    {
        this->DistanceToIn(p, surface);
        return surface->SurfaceNormal(p);
    }
    return this->WhatSurfaceAmIOn(p).SurfaceNormal(p);
}

EInside G4HalfIntersectionSolid::Inside(const G4ThreeVector& p) const
{
    EInside inside = (this->ImplicitFunction(p) > kCarTolerance) ? kOutside
                    :(this->ImplicitFunction(p) < -kCarTolerance) ? kInside
                    : kSurface;
    return inside;
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
}

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