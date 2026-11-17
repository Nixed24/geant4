//
// Created by deck on 11/14/26.
//
#include "G4ImplicitSurfaceCGAL.hh"
#include "G4HalfQuadric.hh"
#include "G4SurfaceMeshCGAL.hh"

FT_GT sphere_function (Point_3_GT p) {
  const FT_GT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1000*1000;
}

G4SurfaceMeshCGAL* make_mesh(G4HalfQuadric &quadric, double sphere_size) {
  G4cout << "G4ImplicitSufaceCGAL make_mesh" << std::endl;
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  std::function<FT_GT(Point_3_GT)> surface_function = [&quadric](Point_3_GT p){
    G4ThreeVector g4_p(p.x(), p.y(), p.z());
    return quadric.ImplicitFunction(g4_p);
  };

  // defining the surface
  Surface_3_GT surface(surface_function,             // pointer to function
                       Sphere_3_GT(CGAL::ORIGIN, sphere_size)); // bounding sphere

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     75,  // radius bound
                                                     75); // distance bound
  // meshing surface
  //CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

  Surface_mesh_3 *sm = new Surface_mesh_3();
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, *sm);
  G4cout << "G4ImplicitSufaceCGAL Final number of points: " << tr.number_of_vertices() << "\n";

  // get mesh point, offset by normal and evaluate implicit function
  Surface_mesh_3::Face_index face_index;
  face_index = *(sm->faces().begin());
  auto vertex = sm->point(*CGAL::vertices_around_face(sm->halfedge(face_index), *sm).begin());
  auto normal = PMP::compute_face_normal(face_index, *sm);

  auto displaced = vertex + 3*normal;

  auto sgn = quadric.ImplicitFunction(G4ThreeVector(CGAL::to_double(displaced.x()),
                                       CGAL::to_double(displaced.y()),
                                       CGAL::to_double(displaced.z())));

  if (sgn < 0)
    CGAL::Polygon_mesh_processing::reverse_face_orientations(*sm);

  G4SurfaceMeshCGAL *g4_sm = new G4SurfaceMeshCGAL(sm);
  delete sm;

  return g4_sm;
}