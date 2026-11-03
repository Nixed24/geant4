#include "G4HalfSurface.hh"
#include <algorithm>

G4double G4HalfSurface::Union(const G4ThreeVector& p, const G4HalfSurface& surface) const
{
    return std::min(this->ImplicitFunction(p), surface.ImplicitFunction(p));
}

G4double G4HalfSurface::Intersection(const G4ThreeVector& p, const G4HalfSurface& surface) const
{
    return std::max(this->ImplicitFunction(p), surface.ImplicitFunction(p));
}

G4double G4HalfSurface::Subtraction(const G4ThreeVector& p, const G4HalfSurface& surface) const
{
    return std::max(this->ImplicitFunction(p), -1*surface.ImplicitFunction(p));
}

std::vector<G4ThreeVector> G4HalfSurface::RayIntersection(const G4ThreeVector& p, const G4ThreeVector& v) const
{
    const G4ThreeVector vUnit = v/v.mag();

    const G4double distanceToSurface = this->DistanceToSurface(p,v);

    if (distanceToSurface == kInfinity) return {};

    return {vUnit*distanceToSurface + p};
}

EInside G4HalfSurface::Inside(const G4ThreeVector& p) const
{
    if (this->ImplicitFunction(p) < 0) return kInside;
    if (this->ImplicitFunction(p) > 0) return kOutside;
    return kSurface;
}

std::ostream& G4HalfSurface::StreamInfo(std::ostream& os) const
{
    os << "--- Half Surface ---" << std::endl;
    return os;
}