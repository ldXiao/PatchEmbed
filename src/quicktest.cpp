//
// Created by Lind Xiao on 7/30/19.
//
#include <map>
#include <vector>
#include <iostream>
#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/eigen.h>
namespace py = pybind11;
int main() {
    py::scoped_interpreter guard{};

    py::exec(R"(
        kwargs = dict(name="World", number=42)
        message = "Hello, {name}! The answer is {number}".format(**kwargs)
        print(message)
    )");

    py::module sys = py::module::import("sys");
    py::module torch = py::module::import("torch");
    py::module os = py::module::import("os");
    sys.attr("path").attr("append")("../src");
    py::module utlis = py::module::import("utlis");
    utlis.attr("test")();
}