#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                    Coordinate;
typedef CGAL::Delaunay_triangulation_3<K>        Delaunay;
typedef Delaunay::Point                          Point;
typedef Delaunay::Vertex_handle                  Vertex_handle;
typedef Delaunay::Cell_handle                    Cell_handle;
typedef Delaunay::Facet                          Facet;
typedef Delaunay::Edge                           Edge; 



std::map<Vertex_handle,int> vertex_map;
std::map<std::pair<Vertex_handle,Vertex_handle>, int> edge_map;
std::map<Facet,int> facet_map;

std::pair<Vertex_handle,Vertex_handle> get_endpoints(const Edge& e) {
  Cell_handle cell = e.first;
  int i = e.second;
  int j = e.third;
  Vertex_handle vh1=cell->vertex(i);
  Vertex_handle vh2=cell->vertex(j);
  return std::make_pair(vh1,vh2);
}

std::vector<int> boundary_of(const Edge& e) {
  Vertex_handle vh1,vh2;
  boost::tie(vh1,vh2)=get_endpoints(e);
  int i = vertex_map[vh1];
  int j = vertex_map[vh2];
  std::vector<int> result;
  result.push_back(i);
  result.push_back(j);
  std::sort(result.begin(),result.end());
  return result;
}

std::vector<int> boundary_of(const Facet& f) {
  Cell_handle ch = f.first;
  int ind = f.second;
  Vertex_handle vh1 = ch->vertex((ind+1)%4);
  Vertex_handle vh2 = ch->vertex((ind+2)%4);
  Vertex_handle vh3 = ch->vertex((ind+3)%4);
  int i = edge_map[std::make_pair(vh1,vh2)];
  int j = edge_map[std::make_pair(vh1,vh3)];
  int k = edge_map[std::make_pair(vh2,vh3)];
  std::vector<int> result;
  result.push_back(i);
  result.push_back(j);
  result.push_back(k);
  std::sort(result.begin(),result.end());
  return result;
}

std::vector<int> boundary_of(const Cell_handle& ch) {
  int i = facet_map[Facet(ch,0)];
  int j = facet_map[Facet(ch,1)];
  int k = facet_map[Facet(ch,2)];
  int l = facet_map[Facet(ch,3)];
  std::vector<int> result;
  result.push_back(i);
  result.push_back(j);
  result.push_back(k);
  result.push_back(l);
  std::sort(result.begin(),result.end());
  return result;
}

std::pair<Coordinate,Coordinate> value_of(const Edge& e) {
  Vertex_handle vh1,vh2;
  boost::tie(vh1,vh2)=get_endpoints(e);  
  Coordinate x1=vh1->point().x();
  Coordinate x2=vh2->point().x();
  Coordinate y1=vh1->point().y();
  Coordinate y2=vh2->point().y();
  return std::make_pair(CGAL::max(x1,x2),CGAL::max(y1,y2));
}
  
std::pair<Coordinate,Coordinate> value_of(const Facet& f) {
  Cell_handle ch = f.first;
  int ind = f.second;
  Vertex_handle vh1 = ch->vertex((ind+1)%4);
  Vertex_handle vh2 = ch->vertex((ind+2)%4);
  Vertex_handle vh3 = ch->vertex((ind+3)%4);
  Coordinate x1=vh1->point().x();
  Coordinate x2=vh2->point().x();
  Coordinate x3=vh3->point().x();
  Coordinate y1=vh1->point().y();
  Coordinate y2=vh2->point().y();
  Coordinate y3=vh3->point().y();
  return std::make_pair(CGAL::max(CGAL::max(x1,x2),x3),
			CGAL::max(CGAL::max(y1,y2),y3));
}

std::pair<Coordinate,Coordinate> value_of(const Cell_handle& ch) {
  Vertex_handle vh1 = ch->vertex(0);
  Vertex_handle vh2 = ch->vertex(1);
  Vertex_handle vh3 = ch->vertex(2);
  Vertex_handle vh4 = ch->vertex(3);
  Coordinate x1=vh1->point().x();
  Coordinate x2=vh2->point().x();
  Coordinate x3=vh3->point().x();
  Coordinate x4=vh4->point().x();
  Coordinate y1=vh1->point().y();
  Coordinate y2=vh2->point().y();
  Coordinate y3=vh3->point().y();
  Coordinate y4=vh4->point().y();
  return std::make_pair(CGAL::max(CGAL::max(x1,x2),CGAL::max(x3,x4)),
			CGAL::max(CGAL::max(y1,y2),CGAL::max(y3,y4)));
}


int main(int argc, char** argv)
{
  int n = (argc>1) ? atoi(argv[1]) : 100;
  
  Delaunay T;
  std::vector<Point> pts;
  CGAL::Random_points_in_sphere_3<Point> rnd;
  for (int i = 0; i != n; ++i) {
    Point p = *rnd++;
    //std::cerr << "Inserting " << p.x() << " " << p.y() << std::endl;
    pts.push_back(p);
  }
  T.insert(pts.begin(),pts.end());
  std::cerr << "Final triangulation has " << T.number_of_vertices()
            << " vertices, " << T.number_of_finite_edges() << " edges, "
	    << T.number_of_finite_facets() << " facets, and "
	    << T.number_of_finite_cells() << " cells." << std::endl;


  int v_count=0;

  for(auto vit = T.finite_vertices_begin();
      vit!=T.finite_vertices_end();
      vit++) {
    Vertex_handle vh = vit;
    vertex_map[vh]=v_count++;
  }
  int e_count=0;
  for(auto eit = T.finite_edges_begin();
      eit!=T.finite_edges_end();
      eit++) {
    Vertex_handle vh1,vh2;
    boost::tie(vh1,vh2)=get_endpoints(*eit);
    edge_map[std::make_pair(vh1,vh2)]=e_count;
    edge_map[std::make_pair(vh2,vh1)]=e_count;
    e_count++;
  }

  int f_count=0;
  for(auto fit = T.finite_facets_begin();
      fit!=T.finite_facets_end();
      fit++) {
    facet_map[*fit]=f_count;
    facet_map[T.mirror_facet(*fit)]=f_count;
    f_count++;
  }

  

    
  
  std::cout << "scc2020\n2\n";


  std::cout << T.number_of_finite_cells() << " "
	    << T.number_of_finite_facets() << " " 
	    << T.number_of_finite_edges() << " "
	    << T.number_of_vertices() << std::endl;
  std::cerr << "N=" << T.number_of_finite_cells()+T.number_of_finite_facets()+T.number_of_finite_edges() << std::endl;

  for(auto cit = T.finite_cells_begin();
      cit!=T.finite_cells_end();
      cit++) {
    auto val = value_of(cit);
    std::vector<int> bd=boundary_of(cit);
    std::cout << std::setprecision(12) << std::fixed << CGAL::to_double(val.first) << " "
	      << std::setprecision(12) << std::fixed<< CGAL::to_double(val.second) << " ; ";
    for(auto it : bd) {
      std::cout << it <<  " ";
    }
    std::cout << std::endl;
  }
  for(auto fit = T.finite_facets_begin();
      fit!=T.finite_facets_end();
      fit++) {
    auto val = value_of(*fit);
    std::vector<int> bd=boundary_of(*fit);
    std::cout << std::setprecision(12) << std::fixed<< CGAL::to_double(val.first) << " "
	      << std::setprecision(12) << std::fixed<< CGAL::to_double(val.second) << " ; ";
    for(auto it : bd) {
      std::cout << it <<  " ";
    }
    std::cout << std::endl;
  }
  for(auto eit = T.finite_edges_begin();
      eit!=T.finite_edges_end();
      eit++) {
    auto val = value_of(*eit);
    std::vector<int> bd=boundary_of(*eit);
    std::cout << std::setprecision(12) << std::fixed<< CGAL::to_double(val.first) << " "
	      << std::setprecision(12) << std::fixed<< CGAL::to_double(val.second) << " ; ";
    for(auto it : bd) {
      std::cout << it <<  " ";
    }
    std::cout << std::endl;
  }
  

  return 0;
}
