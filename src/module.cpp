#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include<list>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
using namespace boost;
namespace py = pybind11;

// supports DELOCX and DELOCY without copying matrix A
// "Serializes" a map. That means, giving an index for every pixel.
// For example: [0123]
// 				[4567]
#define LOC(y,x,Arows) ((y)*(Arows+1)+(x))
// Deserialising. Returns the x or y coordinate for a given index.
#define DELOCY(i,Arows) ((i)/(Arows+1))
#define DELOCX(i,Arows) ((i)%(Arows+1))


// Typedefs for graph_t, vertex_descriptor, edge_descriptor, Edge,
// and out_edge_iterator.
typedef adjacency_list < listS, vecS, directedS,
    no_property, property < edge_weight_t, double > > graph_t;
typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
typedef std::pair<int, int> Edge;
typedef graph_traits < graph_t >::out_edge_iterator out_edge_iterator;


/*
  Global state, needs mlock!
*/


struct tState
{
    graph_t g;						                        // the graph
    property_map<graph_t, edge_weight_t>::type weightmap;   // the weights
    std::vector<vertex_descriptor> p;		                // predecessor map
    std::vector<double> d;				                    // distance map
    vertex_descriptor global_s;				      // cache the start of dijkstra
    int Arows;
} ;

// The graph

std::vector<tState> gState; 

int py_graphFromBitmap(py::array_t<uint8_t> pyo_image, uint8_t walkable=0)
{
    auto A  = pyo_image.mutable_unchecked<2>(); 
    tState state;

	// Create a matrix from the bitmap and store the number of rows for LOC/DELOC
    state.Arows = A.shape(0);
    std::cout << "Rows:" << state.Arows << ", Columns: " << A.shape(1) << "\n";
    
    // prepare the graph and search
    // This is the full graph. All points are vertices
    // hence, some vertices are disconnected and make their
    // own connected component, e.g., black pixels
    graph_t _g(A.shape(0)*A.shape(1));
    state.g = _g;
    // Define macro to check, whether an edge can be added from
    // i,j to i+k, j+l
    auto CHECK_EDGE = [&](int i,int j,int k,int l) { return (
	i+k >= 0 && i+k < A.shape(0) && 
	j+l >= 0 && j+l < A.shape(1) &&
	A(i,j) == walkable &&
	A(i+k,j+l) == walkable );};

	// Iterate through the map and add edges, if possible
    int e=0;
    for (int i=0; i< A.shape(0); i++) {
        for (int j=0; j< A.shape(1); j++) {
            // NOTE!!!
            // left      0 -1
            // right     0  1
            // up       -1  0
            // down      1  0
            if (CHECK_EDGE(i,j,-1,0))
                add_edge(LOC(i,j,state.Arows),LOC(i-1,j,state.Arows),1,state.g);
            if (CHECK_EDGE(i,j,+1,0))
                add_edge(LOC(i,j,state.Arows),LOC(i+1,j,state.Arows),1,state.g);
            if (CHECK_EDGE(i,j,0,-1))
                add_edge(LOC(i,j,state.Arows),LOC(i,j-1,state.Arows),1,state.g);
            if (CHECK_EDGE(i,j,0,1))
                add_edge(LOC(i,j,state.Arows),LOC(i,j+1,state.Arows),1,state.g);
            // DIAGONALS
            if (CHECK_EDGE(i,j,-1,-1))
                add_edge(LOC(i,j,state.Arows),LOC(i-1,j-1,state.Arows),sqrt(2),state.g);
            if (CHECK_EDGE(i,j,-1,1))
                add_edge(LOC(i,j,state.Arows),LOC(i-1,j+1,state.Arows),sqrt(2),state.g);
            if (CHECK_EDGE(i,j,1,-1))
                add_edge(LOC(i,j,state.Arows),LOC(i+1,j-1,state.Arows),sqrt(2),state.g);
            if (CHECK_EDGE(i,j,1,1))
                add_edge(LOC(i,j,state.Arows),LOC(i+1,j+1,state.Arows),sqrt(2),state.g);
        }
    }
    
    std::cout << "Added " << num_vertices(state.g) << " vertices and "
				  << num_edges(state.g)<<" edges.\n";
				  	
	// Store the state
    gState.push_back(state);
    int index = gState.size()-1;
    std::cout << "State is " << index << std::endl;
    return index;

}





// -------------
// pure C++ code
// -------------

std::vector<int> multiply(const std::vector<double>& input)
{
  std::vector<int> output(input.size());

  for ( size_t i = 0 ; i < input.size() ; ++i )
    output[i] = 10*static_cast<int>(input[i]);

  return output;
}

// ----------------
// Python interface
// ----------------


struct predicate_violated {};
/*This function actually performs our homotopy test*/
bool py_empty(py::array_t<uint8_t> pyo_image, py::array_t<double> pyo_points, uint8_t allowed=0)
{
    auto image  = pyo_image.mutable_unchecked<2>(); 
    auto points = pyo_points.unchecked<2>();

    const int polyCorners = points.shape(0);
    std::vector<double> nodeX;
    nodeX.reserve(polyCorners);
    try {
    for (size_t pixelY=0; pixelY < image.shape(1); pixelY++)
    {
	int i,j;
    
    auto polyX = [&](int row) {return points(row,0);};
    auto polyY = [&](int row) {return points(row,1);};
    
	nodeX.clear(); j=polyCorners-1;
	for (i=0; i<polyCorners; i++)
	{
	    if (polyY(i)<(double) pixelY && polyY(j)>=(double) pixelY
		||  polyY(j)<(double) pixelY && polyY(i)>=(double) pixelY)
	    {
		nodeX.push_back( (polyX(i)+(pixelY-polyY(i))/(polyY(j)-polyY(i))*(polyX(j)-polyX(i))));
	   }
	  j=i;
    } // for i ranging polyCorners
    std::sort(nodeX.begin(), nodeX.end());

    const auto IMAGE_LEFT=0;
    const auto IMAGE_RIGHT = image.shape(1);
    
     for (i=0; i<nodeX.size(); i+=2) {
	if   (nodeX[i  ]>=IMAGE_RIGHT) break;
	
	if   (nodeX[i+1]> IMAGE_LEFT ) {
	  if (nodeX[i  ]< IMAGE_LEFT ) nodeX[i  ]=IMAGE_LEFT ;
	  if (nodeX[i+1]> IMAGE_RIGHT) nodeX[i+1]=IMAGE_RIGHT;
	  for (auto pixelX=(int)nodeX[i]; pixelX<=(int) (nodeX[i+1]+0.5); pixelX++)
	      if (image(pixelX,pixelY) != allowed)
	         throw(predicate_violated());
	  }
	  }
   } // for y
   }catch(predicate_violated v)
    {
    return false;

   }
   return true;

 }


void py_dijkstra(int handle, py::array_t<uint32_t> pyo_p)
{
    tState *state = &gState[handle];  // @TODO: bound check? else crash
    auto p = pyo_p.unchecked<1>();
    
    vertex_descriptor s = LOC(p(0),p(1),state->Arows);
    state->global_s = s;
    
    std::cout << "Search: " << s << "\n"; 
    std::cout << "Using " << num_vertices(state->g) << " vertices and "
				  << num_edges(state->g) << " edges.\n";
    state->weightmap = get(edge_weight, state->g);	// weightmap
    state->p.resize(num_vertices(state->g));		// predecessor
    state->d.resize(num_vertices(state->g));		// distances

	// Perform dijkstra
    dijkstra_shortest_paths(state->g, s,
		predecessor_map(boost::make_iterator_property_map(state->p.begin(), get(boost::vertex_index, state->g))).
		distance_map(boost::make_iterator_property_map(state->d.begin(), get(boost::vertex_index, state->g))));
    std::cout << "Did it " << std::endl;
}


py::array_t<uint32_t> py_getpath(int handle, py::array_t<uint32_t> pyo_p, double penalty=1.0)
{
    tState *state = &gState[handle];  // @TODO: bound check? else crash
    auto p = pyo_p.unchecked<1>();
    
    vertex_descriptor e = LOC(p(0),p(1),state->Arows);
    if (state->p[e] == e) { // no path found
	  return (py::array_t<uint32_t>(0));
    }
    // Get the shortest path from s to e


    std::list<vertex_descriptor> shortest_path;
    for(vertex_descriptor v = e;;v=state->p[v] ) {
        shortest_path.push_front(v);
        if ( penalty != 1) {
			edge_descriptor e1; bool found;
			tie(e1, found) = edge(state->p[v], v, state->g);
			if(found) {
//				octave_stdout << e1 << " : " << get(state->weightmap, e1);
				put(state->weightmap, e1, get(state->weightmap, e1)*2);
//				octave_stdout << " -> " << get(state->weightmap, e1) << "\n";
			}
		}
        if(state->p[v] == v)
            break;
    }    
// allocate py::array (to pass the result of the C++ function to Python)
  // correctly transfer with a capsule
  uint32_t *result_ptr = new uint32_t[shortest_path.size()*2];
  uint32_t *ptr = result_ptr;
  for(const auto pt:shortest_path) {
        *ptr++ =  DELOCY(pt,state->Arows);
        *ptr++  = DELOCX(pt,state->Arows);
    }

      py::capsule free_when_done(result_ptr, [](void *f) {
            uint32_t *foo = reinterpret_cast<uint32_t *>(f);
            //std::cerr << "Element [0] = " << foo[0] << "\n";
            //std::cerr << "freeing memory @ " << f << "\n";
            delete[] foo;
        });

   return py::array_t<uint32_t>(
            shortest_path.size()*2, result_ptr,free_when_done); // numpy array references this parent
 }


 

void py_getweights(int handle, py::array_t<double> pyo_image)
{
    auto weights  = pyo_image.mutable_unchecked<2>(); 
    tState *state = &gState[handle];  // @TODO: bound check? else crash
    
    for (auto r =0; r < weights.shape(0); r++)
      for (auto c = 0; c < weights.shape(1); c++)
      {
	  vertex_descriptor v = LOC(r,c,state->Arows);
	  out_edge_iterator out_i, out_end;
	  edge_descriptor e;
	  double sum = 0;
	  double N = 0;
	  for (tie(out_i, out_end) = out_edges(v, state->g); out_i != out_end; ++out_i)
	  {
			sum += get(state->weightmap, *out_i);
			N++;
	  }
	  weights(r,c) = (N==0)?0:sum/N;
      }
}
/*This function is for debugging, but it is not needed*/


void py_fill(py::array_t<uint8_t> pyo_image, py::array_t<double> pyo_points)
{
    auto image  = pyo_image.mutable_unchecked<2>(); 
    auto points = pyo_points.unchecked<2>();

    const int polyCorners = points.shape(0);
    std::vector<double> nodeX;
    nodeX.reserve(polyCorners);
    
    for (size_t pixelY=0; pixelY < image.shape(1); pixelY++)
    {
//	std::cout << "scanning " << pixelY << std::endl;
    int nodes = 0; int i,j;

    
    auto polyX = [&](int row) {return points(row,0);};
    auto polyY = [&](int row) {return points(row,1);};
    
	nodeX.clear(); j=polyCorners-1;
	for (i=0; i<polyCorners; i++)
	{
	    if (polyY(i)<(double) pixelY && polyY(j)>=(double) pixelY
		||  polyY(j)<(double) pixelY && polyY(i)>=(double) pixelY)
	    {
		nodeX.push_back( (polyX(i)+(pixelY-polyY(i))/(polyY(j)-polyY(i))*(polyX(j)-polyX(i))));
	   }
	  j=i;
	  // now we have the X of the scanline
    } // for i ranging polyCorners
    std::sort(nodeX.begin(), nodeX.end());
    //std::cout << "On Scanline " << pixelY << " we have " << std::endl;
  //  Fill the pixels between node pairs.

    const auto IMAGE_LEFT=0;
    const auto IMAGE_RIGHT = image.shape(1);
    
    
     for (i=0; i<nodeX.size(); i+=2) {
	if   (nodeX[i  ]>=IMAGE_RIGHT) break;
	
	if   (nodeX[i+1]> IMAGE_LEFT ) {
	  if (nodeX[i  ]< IMAGE_LEFT ) nodeX[i  ]=IMAGE_LEFT ;
	  if (nodeX[i+1]> IMAGE_RIGHT) nodeX[i+1]=IMAGE_RIGHT;
	  for (auto pixelX=(int)nodeX[i]; pixelX<=(int) (nodeX[i+1]+0.5); pixelX++)
	      image(pixelX,pixelY) = 1;
	  }
	  }



    
    } // for y


    
}



/*// wrap C++ function with NumPy array IO
py::array_t<int> py_multiply(py::array_t<double, py::array::c_style | py::array::forcecast> array)
{
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> array_vec(array.size());

  // copy py::array -> std::vector
  std::memcpy(array_vec.data(),array.data(),array.size()*sizeof(double));

  // call pure C++ function
  std::vector<int> result_vec = multiply(array_vec);

  // allocate py::array (to pass the result of the C++ function to Python)
  auto result        = py::array_t<int>(array.size());
  auto result_buffer = result.request();
  int *result_ptr    = (int *) result_buffer.ptr;

  // copy std::vector -> py::array
  std::memcpy(result_ptr,result_vec.data(),result_vec.size()*sizeof(int));

  return result;
  }
  */ 

// wrap as Python module
PYBIND11_MODULE(cfsrouting,m)
{
  m.doc() = "pybind11 example plugin";

//  m.def("multiply", &py_multiply, "Convert all entries of an 1-D NumPy-array to int and multiply by 10");
//  m.def("test", &py_multiply, "Convert all entries of an 1-D NumPy-array to int and multiply by 10");
  m.def ("fill", &py_fill,py::arg().noconvert(), py::arg().none());
  m.def("empty", &py_empty);
  m.def("graphFromBitmap", &py_graphFromBitmap);
  m.def("dijkstra", &py_dijkstra);
  m.def("getpath", &py_getpath);
  m.def("getweights", &py_getweights, py::arg().noconvert(),py::arg().noconvert());
}


/*template<typename vtype>
py::array wrap(vtype v)
{
   return py::array(v.size(),v.data()); // does a copy
}


PYBIND11_MODULE(spatialfacet,m) {
    py::class_<SpatialFacetMiner>(m, "SpatialFacetMiner")
    .def(py::init<>())
    .def("add_database", &SpatialFacetMiner::add_database)
    .def("query", [](SpatialFacetMiner &m,std::string query_string, int first, int max, int check_at_least){return m.query(query_string,first,max,check_at_least, false);})
    .def("query_with_data", [](SpatialFacetMiner &m,std::string query_string, int first, int max, int check_at_least){return m.query(query_string,first,max,check_at_least, true);})
    .def("getSpyData", [](SpatialFacetMiner &m){
	auto &spy = m.getSpy();
	return py::make_tuple(
	    wrap(spy.coords[0]),
	    wrap(spy.coords[1]),
	    wrap(spy.docids),
	    wrap(spy.weights)
	 );
	 })
   .def("getSpyStringData",[](SpatialFacetMiner &m){
      auto &spy = m.getSpy();
      return py::make_tuple(spy.value1, spy.values);

      })

    .def ("augment",[](SpatialFacetMiner &m, std::string query_string,
		       py::array_t<double, py::array::c_style | py::array::forcecast> documents,
		       int n_terms)
    {
	    std::vector<int> stl_documents;
	    auto r = documents.unchecked<1>();
	    for (py::ssize_t i =0; i < r.shape(0); i++)
	      stl_documents.push_back(r(i));
	    
	    std::vector<string> terms; std::vector<double> weights; std::string query_out;	    
//    void augment_query_from_documents(std::string query_string, std::vector<int> documents, int n_terms,
//				      std::vector<std::string> &terms, std::vector<double> &weights, std::string &query_out)
	    
	    m.augment_query_from_documents(query_string, stl_documents,n_terms, terms, weights, query_out);
	   return py::make_tuple(terms, weights, query_out);
	    
    });



      ;
    
}
*/
