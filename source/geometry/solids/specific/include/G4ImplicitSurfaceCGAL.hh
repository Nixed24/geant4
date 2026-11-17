//
// Created by deck on 11/14/26.
//

#ifndef GEANT4_G4IMPLICITSURFACECGAL_HH
#define GEANT4_G4IMPLICITSURFACECGAL_HH
#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3_GT;
typedef GT::Point_3 Point_3_GT;
typedef GT::FT FT_GT;
typedef FT_GT (*Function)(Point_3_GT);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3_GT;
typedef CGAL::Surface_mesh<Point_3_GT> Surface_mesh_GT;

namespace PMP = CGAL::Polygon_mesh_processing;

class G4HalfQuadric;
class G4SurfaceMeshCGAL;

G4SurfaceMeshCGAL* make_mesh(G4HalfQuadric &quadric, double sphere_size);
#endif  // GEANT4_G4IMPLICITSURFACECGAL_HH
