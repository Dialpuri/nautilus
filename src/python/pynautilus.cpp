//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "nautilus-include.h"
#include "nautilus-run.h"

namespace nb = nanobind;

using namespace nb::literals;

NB_MODULE(nautilus_module, m) {
nb::class_<NautilusInput>(m, "Input")
            .def(nb::init< const std::string&, // mtzin
             const std::string&, // seqin
             const std::string&, // pdbin
             const std::string&, // predin
             const std::string&, // colin_fo
             const std::string&, // colin_hl
             const std::string&, // colin_phifom
             const std::string&,  // colin_fc
             const std::string&  >()); // colin_free 

nb::class_<NautilusOutput>(m, "Output")
            .def(nb::init< const std::string& >()); // pdbout

m.def("run", &run);
}