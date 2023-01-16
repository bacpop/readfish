#include "bindings.h"

PYBIND11_MODULE(query_cpp, m) {
  m.doc() = "Finds strings in bifrost DBG";

  py::class_<Graph, std::unique_ptr<Graph>>(m, "Graph")
        .def(py::init<>())
        .def("build", &Graph::build)
        .def("read", &Graph::read)
        .def("query", &Graph::query);
}