#include "G4HalfUnionSolid.hh"

G4HalfUnionSolid::G4HalfUnionSolid(const std::vector<G4HalfIntersectionSolid*>& intersectionSolids,
                                   const std::vector<G4HalfIntersectionSolid*>& subtractionSolids,
                                   const G4String& name)
                                      :  G4VSolid(name),
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

    for (const G4HalfIntersectionSolid* s : intersectionSolids_)
    {
        intersectionImplicitFunction = std::min(intersectionImplicitFunction, s->ImplicitFunction(p));
    }

    for (const G4HalfIntersectionSolid* s : subtractionSolids_)
    {
        subtractionImplicitFunction = std::min(subtractionImplicitFunction, s->ImplicitFunction(p));
    }

    G4double implicitFunction = std::max(intersectionImplicitFunction, -subtractionImplicitFunction);

    return implicitFunction;
}

EInside G4HalfUnionSolid::Inside(const G4ThreeVector& p) const
{
    G4double implicitFunction = this->ImplicitFunction(p);

    return  (implicitFunction < -kCarTolerance) ? kInside
           :(implicitFunction > kCarTolerance)  ? kOutside
           :                                      kSurface;
}

