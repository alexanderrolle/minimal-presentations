#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;

typedef Surface_mesh::Vertex_index Vertex_index;
typedef Surface_mesh::Edge_index Edge_index;

std::map<Edge_index,int> edge_map;
std::map<Vertex_index,int> vertex_map;

std::pair<double,double> get_max(std::vector<Point_3> vhs) {
  
  double max_x = CGAL::to_double(vhs[0].x());
  double max_y = CGAL::to_double(vhs[0].y());

  for(int i=1;i<vhs.size();i++) {
    double curr_x = CGAL::to_double(vhs[i].x());
    double curr_y = CGAL::to_double(vhs[i].y());
    if(curr_x>max_x) {
      max_x=curr_x;
    }
    if(curr_y>max_y) {
      max_y=curr_y;
    }
  }
  return std::make_pair(max_x,max_y);
}




int main(int argc, char* argv[])
{

  int n = (argc>1) ? atoi(argv[1]) : 100;
  std::vector<Point_3> points;
  CGAL::Random_points_on_sphere_3<Point_3> rnd;

  for(int i=0;i<n;i++) {
    Point_3 p = *rnd++;
    //std::cout << p.x() << " " << p.y() << " " << p.z() << " - " << p.x()*p.x() + p.y()*p.y()+p.z()*p.z() << std::endl;
    points.push_back(p);
  }
  Surface_mesh sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);
  std::cerr << "The convex hull contains " << sm.num_vertices() << " vertices, "
	    << sm.num_edges() << " edges, and " << sm.num_faces() << " faces"  << std::endl;
  std::cerr << "N=" << sm.num_edges()+sm.num_faces() << std::endl;

  int v_count=0;
  for(auto it = sm.vertices().begin();it!=sm.vertices().end();it++) {
    Point_3 p = sm.point(*it);
    //std::cout << "Found point " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    vertex_map[*it]=v_count++;
  }
  int e_count=0;
  for(auto it = sm.edges().begin();it!=sm.edges().end();it++) {
    edge_map[*it] = e_count++;

  }

  std::cout << "firep\nx-coordinate\ny-coordinate\n";
  std::cout << sm.number_of_faces() << " " << sm.number_of_edges() << " " << sm.number_of_vertices() << std::endl;

  
  for(auto it = sm.faces().begin();it!=sm.faces().end();it++) {
    std::vector<Point_3> endpoints;
    std::vector<int> bd;
    auto start_halfedge = sm.halfedge(*it);
    auto he = start_halfedge;
    do {
      endpoints.push_back(sm.point(sm.target(he)));
      bd.push_back(edge_map[sm.edge(he)]);
      he=sm.next(he);
    }
    while (he!=start_halfedge);
    assert(endpoints.size()==3);
    assert(bd.size()==3);
    std::pair<double,double> val = get_max(endpoints);
    std::sort(bd.begin(),bd.end());
    std::cout << val.first << " " << val.second << " ; " << bd[0] << " " << bd[1] << " " << bd[2] << " " << std::endl;
    
    
  }

  for(auto it = sm.edges().begin();it!=sm.edges().end();it++) {
    std::vector<Point_3> endpoints;
    auto halfedge = sm.halfedge(*it);
    Vertex_index v1 = sm.target(halfedge);
    Vertex_index v2 = sm.target(sm.opposite(halfedge));
    int i = vertex_map[v1];
    int j = vertex_map[v2];
    
    
    endpoints.push_back(sm.point(sm.target(halfedge)));
    endpoints.push_back(sm.point(sm.target(sm.opposite(halfedge))));
    Point_3& p = endpoints[0];
    //std::cout << "1st point " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    Point_3& q = endpoints[1];
    //std::cout << "2nd point " << q.x() << " " << q.y() << " " << q.z() << std::endl;
    std::pair<double,double> val = get_max(endpoints);
    
    std::cout << val.first << " " << val.second << " ; " << i << " " << j << " " << std::endl;
  }
  

  return 0;
}
