#include "G4HalfPlane.hh"

G4HalfPlane::G4HalfPlane(const G4ThreeVector& normalIn, const G4double& distanceIn)
            :normal(normalIn.unit()),
             distance(distanceIn) {}

G4double G4HalfPlane::ImplicitFunction(const G4ThreeVector& p) const
{
    return normal.dot(p) - distance;
}

G4ThreeVector G4HalfPlane::closestSurfacePoint(const G4ThreeVector& p) const
{
    return this->RayIntersection(p, this->normal)[0];
}

std::vector<G4ThreeVector> G4HalfPlane::RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4ThreeVector vUnit = v.unit();

    const G4double distanceToSurface = this->DistanceToSurface(p,v);

    if (distanceToSurface == kInfinity && this->Inside(p) == kSurface) return {p};
    if (distanceToSurface == kInfinity) return {};

    return {vUnit*distanceToSurface + p};
}

// Returns the signed distance to surface
G4double G4HalfPlane::DistanceToSurface(const G4ThreeVector& p) const
{
    // For a plane, the implicit function is the same as the signed distance
    return this->ImplicitFunction(p);
}

G4double G4HalfPlane::DistanceToSurface(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

    const G4ThreeVector vUnit = v.unit();

    G4double projection = normal.dot(vUnit);
    G4double implicitFunction = this->ImplicitFunction(p);

    if (this->Inside(p) == kSurface) return 0;
    if (projection < -kCarTolerance && implicitFunction < -kCarTolerance) return kInfinity;
    if (projection > kCarTolerance && implicitFunction > kCarTolerance) return kInfinity;
    if (-kCarTolerance <= projection && projection <= kCarTolerance) return kInfinity; // Parallel to plane

    return -implicitFunction/projection;
}

void G4HalfPlane::Transform(const G4AffineTransform& t)
{
    const G4ThreeVector& normalPoint = this->distance*this->normal;
    const G4ThreeVector& transformedPoint = t.TransformPoint(normalPoint);

    this->normal = t.TransformAxis(normal).unit();
    this->distance = normal.dot(transformedPoint);
}
