"""
  GITFans.jl is a collection of Julia functions for computations
  like those shown in the paper
  [J. Boehm, S. Keicher, Y. Ren:
  Computing GIT-Fans with Symmetry and the Mori Chamber Decomposition of M06bar](https://arxiv.org/abs/1603.09241).

  The code relies on the package [Oscar.jl](https://github.com/oscar-system/Oscar.jl).
"""
module GITFans

# A first version of the code had been written with focus on
# combining the Julia interfaces of GAP, polymake, and Singular,
# without using Oscar objects.
# The functionality is still available, inside the submodule `OLD`.
include("code_interfaces.jl")

# This file contains the "oscarified" version of the code.
include("code_oscar.jl")

end
