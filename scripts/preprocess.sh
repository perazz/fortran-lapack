#!/bin/sh
#
#
# preprocess sources with fypp
fypp "../src/stdlib_linalg_solve.fypp" > ../src/stdlib_linalg_solve.f90
fypp "../src/stdlib_linalg_inverse.fypp" > ../src/stdlib_linalg_inverse.f90
fypp "../test/test_linalg_solve.fypp" > ../test/test_linalg_solve.f90
