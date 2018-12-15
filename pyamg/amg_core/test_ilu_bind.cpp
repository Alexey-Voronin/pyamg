// DO NOT EDIT: this file is generated

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "test_ilu.h"

namespace py = pybind11;

template<class I, class F>
void _ilu_bsr_fact(
                const I n,
              const I lof,
    const I rows_in_block,
    const I cols_in_block,
py::array_t<I> & row_ptr_arr,
py::array_t<I> & cols_arr,
py::array_t<F> & val_arrIn,
py::array_t<I> & row_ptr_out,
py::array_t<I> & cols_out,
 py::array_t<F> & val_out,
 py::array_t<I> & out_len
                   )
{
    auto py_row_ptr_arr = row_ptr_arr.unchecked();
    auto py_cols_arr = cols_arr.unchecked();
    auto py_val_arrIn = val_arrIn.unchecked();
    auto py_row_ptr_out = row_ptr_out.mutable_unchecked();
    auto py_cols_out = cols_out.mutable_unchecked();
    auto py_val_out = val_out.mutable_unchecked();
    auto py_out_len = out_len.mutable_unchecked();
    const I *_row_ptr_arr = py_row_ptr_arr.data();
    const I *_cols_arr = py_cols_arr.data();
    const F *_val_arrIn = py_val_arrIn.data();
    I *_row_ptr_out = py_row_ptr_out.mutable_data();
    I *_cols_out = py_cols_out.mutable_data();
    F *_val_out = py_val_out.mutable_data();
    I *_out_len = py_out_len.mutable_data();

    return ilu_bsr_fact<I, F>(
                        n,
                      lof,
            rows_in_block,
            cols_in_block,
             _row_ptr_arr, row_ptr_arr.shape(0),
                _cols_arr, cols_arr.shape(0),
               _val_arrIn, val_arrIn.shape(0),
             _row_ptr_out, row_ptr_out.shape(0),
                _cols_out, cols_out.shape(0),
                 _val_out, val_out.shape(0),
                 _out_len, out_len.shape(0)
                              );
}

PYBIND11_MODULE(test_ilu, m) {
    m.doc() = R"pbdoc(
    Pybind11 bindings for test_ilu.h

    Methods
    -------
    ilu_bsr_fact
    )pbdoc";

    py::options options;
    options.disable_function_signatures();

    m.def("ilu_bsr_fact", &_ilu_bsr_fact<int, double>,
        py::arg("n"), py::arg("lof"), py::arg("rows_in_block"), py::arg("cols_in_block"), py::arg("row_ptr_arr").noconvert(), py::arg("cols_arr").noconvert(), py::arg("val_arrIn").noconvert(), py::arg("row_ptr_out").noconvert(), py::arg("cols_out").noconvert(), py::arg("val_out").noconvert(), py::arg("out_len").noconvert(),
R"pbdoc(
)pbdoc");

}

