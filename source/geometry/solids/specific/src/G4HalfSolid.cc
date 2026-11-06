#include "G4HalfSolid.hh"

EInside G4HalfSolid::Inside(const G4ThreeVector& p) const
{
    G4double implicitFunction = this->ImplicitFunction(p);

    return (implicitFunction > kCarTolerance) ? kOutside
          :(implicitFunction < -kCarTolerance) ? kInside
          : kSurface;
}
