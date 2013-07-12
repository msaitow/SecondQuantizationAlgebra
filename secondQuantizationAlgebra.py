#    file:  secondQuantizationAlgebra.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Collects the various functions associated with
#           the second quantization algebra code for easy
#           importing, as in
#           "import secondQuantizationAlgebra.py as sqa".
#
# (c) 2008-2009 Eric Neuscamman (eric.neuscamman@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)

from sqaOptions import \
  options

from sqaIndex import \
  index

from sqaSymmetry import \
  symmetry

from sqaTensor import \
  creDesTensor,       \
  creOp,              \
  desOp,              \
  kroneckerDelta,     \
  sfExOp,             \
  tensor

from sqaTerm import  \
  combineTerms,      \
  multiplyTerms,     \
  removeCoreOpPairs, \
  removeCoreOps_sf,  \
  removeVirtOps_sf,  \
  sortOps,           \
  term,              \
  termChop

from sqaNormalOrder import \
  normalOrder,             \
  contractCoreOps_sf

from sqaCommutator import \
  commutator

from sqaMisc import \
  allDifferent,     \
  get_num_perms,    \
  makePermutations, \
  makeTuples

from sqaMisc2 import      \
  assign_rdm_types,       \
  combine_transpose,      \
  convert_ops_to_rdms_so

from sqaDecomposition_sf import  \
  decomp_3op_to_2op_2rdm_sf,     \
  decomp_3ops_to_2ops_2rdms_sf,  \
  decomp_3rdm_to_2rdm_sf,        \
  decomp_3rdms_to_2rdms_sf,      \
  decomp_4op_to_2op_2rdm_sf,     \
  decomp_4ops_to_2ops_2rdms_sf,  \
  decomp_4rdm_to_2rdm_sf,        \
  decomp_4rdm_to_3rdm_sf,        \
  decomp_4rdms_to_2rdms_sf,      \
  decomp_4rdms_to_3rdms_sf

from sqaDecomposition_so import  \
  decomp_3op_to_2op_2rdm_so,     \
  decomp_3op_to_2op_3rdm_so,     \
  decomp_3ops_to_2ops_2rdms_so,  \
  decomp_3ops_to_2ops_3rdms_so,  \
  decomp_3rdm_to_2rdm_so,        \
  decomp_3rdms_to_2rdms_so,      \
  decomp_4op_to_2op_2rdm_so,     \
  decomp_4op_to_2op_3rdm_so,     \
  decomp_4ops_to_2ops_2rdms_so,  \
  decomp_4ops_to_2ops_3rdms_so,  \
  decomp_4rdm_to_2rdm_so,        \
  decomp_4rdm_to_3rdm_so,        \
  decomp_4rdms_to_2rdms_so,      \
  decomp_4rdms_to_3rdms_so

from sqaConvert import \
  convert2Mulliken,    \
  convertF90_0,        \
  convertF90_1,        \
  convertF90_ERI,      \
  convertTeX,          \
  factorize
