#ifndef GEANT4_G4HALFSOLID_HH
#define GEANT4_G4HALFSOLID_HH

#include "G4VSolid.hh"
#include "G4HalfSurface.hh"

class G4HalfSolid : public G4VSolid
{
protected:
    struct RayHit
    {
        const G4HalfSurface* targetSurface;
        G4bool isSubtraction;
        G4ThreeVector intersection;
        G4double distance;
    };

public:
    G4HalfSolid(const G4String& name) : G4VSolid(name) {};



};

#endif //GEANT4_G4HALFSOLID_HH
