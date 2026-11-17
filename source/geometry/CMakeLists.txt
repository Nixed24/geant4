# - G4geometry category builds
add_definitions(-DG4GEOM_ALLOC_EXPORT)

include(biasing/sources.cmake)
include(divisions/sources.cmake)
include(magneticfield/sources.cmake)
include(management/sources.cmake)
include(navigation/sources.cmake)
include(solids/Boolean/sources.cmake)
include(solids/CSG/sources.cmake)
include(solids/specific/sources.cmake)
include(volumes/sources.cmake)

geant4_add_category(G4geometry MODULES
  G4geombias
  G4geomdivision
  G4magneticfield
  G4geometrymng
  G4navigation
  G4geomBoolean
  G4csg
  G4specsolids
  G4volumes)
set(CGAL_DIR "/home/deck/Documents/mphys_project/2_my_geant4/g4/cgal/install/lib/cmake/CGAL")
find_package(CGAL 5.6 REQUIRED)
#set(Boost_DIR "/home/deck/Documents/mphys_project/2_my_geant4/g4/boost/install/lib/cmake/Boost-1.92.0")
#find_package(Boost 1.92.0 REQUIRED)
