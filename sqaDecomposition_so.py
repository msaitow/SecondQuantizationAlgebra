#    file:  sqaDecomposition_so.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Contains functions for the various decompositions
#           used in Canonical Transformation theory with
#           respect to spin-orbital operators.
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


from sqaIndex import index
from sqaTensor import tensor, creOp, desOp
from sqaTerm import term, multiplyTerms, termChop
from sqaMisc import makePermutations, get_num_perms
from sqaMisc2 import assign_rdm_types
from sqaOptions import options


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3rdms_to_2rdms_so(inTerms, d3name, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Decomposes any 3-particle spin-orbital RDMs in inTerms into products of 1 and 2
  particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose 3RDMs and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the rdm variable
    rdm = False

    # Search for 3-particle RDMs in the term
    for j in range(len(t.tensors)-1, -1, -1):
      if t.tensors[j].name == d3name:
        rdm = t.tensors[j]
        pre_term = term(t.numConstant, t.constants, t.tensors[0:j])
        post_term = term(1.0, [], t.tensors[j+1:])
        break

    # If the term has no 3RDMs, incrment the index
    if rdm is False:
      i += 1

    # If the term has a 3RDM, decompose it
    else:
      opCount += 1
      decomp = decomp_3rdm_to_2rdm_so(rdm, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 3-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdms_to_2rdms_so(inTerms, d4name, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Decomposes any 4-particle spin-orbital RDMs in inTerms into products of 1 and 2
  particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose 4RDMs and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the rdm variable
    rdm = False

    # Search for 4-particle RDMs in the term
    for j in range(len(t.tensors)-1, -1, -1):
      if t.tensors[j].name == d4name:
        rdm = t.tensors[j]
        pre_term = term(t.numConstant, t.constants, t.tensors[0:j])
        post_term = term(1.0, [], t.tensors[j+1:])
        break

    # If the term has no 4RDMs, incrment the index
    if rdm is False:
      i += 1

    # If the term has a 4RDM, decompose it
    else:
      opCount += 1
      decomp = decomp_4rdm_to_2rdm_so(rdm, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdms_to_3rdms_so(inTerms, d4name, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Decomposes any 4-particle spin-orbital RDMs in inTerms into
  products of 1, 2, and 3 particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose 4RDMs and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the rdm variable
    rdm = False

    # Search for 4-particle RDMs in the term
    for j in range(len(t.tensors)-1, -1, -1):
      if t.tensors[j].name == d4name:
        rdm = t.tensors[j]
        pre_term = term(t.numConstant, t.constants, t.tensors[0:j])
        post_term = term(1.0, [], t.tensors[j+1:])
        break

    # If the term has no 4RDMs, incrment the index
    if rdm is False:
      i += 1

    # If the term has a 4RDM, decompose it
    else:
      opCount += 1
      decomp = decomp_4rdm_to_3rdm_so(rdm, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3rdm_to_2rdm_so(d3, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Returns a list of terms composing the decomposition of the 4 body reduced density matrix
  in terms of one and two body rdms"
  """

  # d3 -- Four body rdm to be decomposed
  # d1_aa -- One body rdm
  # d1_bb -- One body rdm
  # d2_aaaa -- Two body rdm
  # d2_bbbb -- Two body rdm
  # d2_abab -- Two body rdm

  # The types of each index on the input tensors must be the same except for the spin variable.
  # The spin variable must have exactly two possible values.

  # Check input
  if ( not isinstance(d3, tensor) ) or (len(d3.indices) != 6):
    raise TypeError, "d3 must be a tensor with 6 indices"
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d3.name, d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Initialize the return value
  decomp = []

  # Prepare terms from    sum (-1)^p %gamma^{i_1}_{i_4} %lambda^{i_2`i_3}_{i_5`i_6}
  for p0 in range(3):
    temp = range(3)
    del temp[p0]
    p1 = temp[0]
    p2 = temp[1]
    ti = [p0, p1, p2]
    for q0 in range(3):
      temp = range(3)
      del temp[q0]
      q1 = temp[0]
      q2 = temp[1]
      if q1 == p2 or q2 == p1:
        bi = [q0, q2, q1]
      else:
        bi = [q0, q1, q2]

      # determine the sign of the term based on the number of index permutations
      n_perm = get_num_perms(ti,bi)
      if n_perm > 1:
        raise ValueError, "n_perm = %i which is > 1.  This was unexpected." %n_perm
      sign = (-1)**n_perm

      # create the 1RDM 2RDM term
      coeff = sign
      bj = [i+3 for i in bi]
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [d3.indices[i] for i in (ti[0:1] + bj[0:1])]
      d2_copy_0 = d2.copy()
      d2_copy_0.indices = [d3.indices[i] for i in (ti[1:3] + bj[1:3])]
      decomp.append(term(coeff, [], [d1_copy_0, d2_copy_0]))

      # create the 1RDM 1RDM 1RDM terms
      coeff = sign * (-1)
      bj = [i+3 for i in bi]
      d1_copies = []
      for j in range(3):
        d1_copies.append(d1.copy())
        d1_copies[-1].indices = [d3.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
      decomp.append(term(coeff, [], d1_copies))
      coeff = sign
      bj = bi[0:1] + bi[2:3] + bi[1:2]
      bj = [i+3 for i in bj]
      d1_copies = []
      for j in range(3):
        d1_copies.append(d1.copy())
        d1_copies[-1].indices = [d3.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
      decomp.append(term(coeff, [], d1_copies))

  # Compute terms from the sum involving only one particle density matrices
  ti = range(3)
  for bi in makePermutations(3):
    
    # Determine the number of permutations
    n_perm = get_num_perms(ti,bi)

    # Determine the sign
    sign = (-1)**n_perm

    # Compute 1RDM 1RDM 1RDM 1RDM term
    bj = [i+3 for i in bi]
    coeff = sign
    d1_copies = []
    for j in range(3):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d3.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab])

  # Combine like terms
  for t in decomp:
    for ten in t.tensors:
      t.scale(ten.sortIndeces())
    t.tensors.sort()
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdm_to_2rdm_so(d4, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Returns a list of terms composing the decomposition of the 4 body reduced density matrix
  in terms of one and two body rdms"
  """

  # d4 -- Four body rdm to be decomposed
  # d1_aa -- One body rdm
  # d1_bb -- One body rdm
  # d2_aaaa -- Two body rdm
  # d2_bbbb -- Two body rdm
  # d2_abab -- Two body rdm

  # The types of each index on the input tensors must be the same except for the spin variable.
  # The spin variable must have exactly two possible values.

  # Check input
  if ( not isinstance(d4, tensor) ) or (len(d4.indices) != 8):
    raise TypeError, "d4 must be a tensor with 8 indices"
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d4.name, d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Initialize the return value
  decomp = []

  # Compute terms from the sum involving two 2-body cumulants
  for p in range(4):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(r+1,4):
          if s == p or r == q:
            bi = [s,r]
          else:
            bi = [r,s]
          temp = range(4)
          del temp[s]
          del temp[r]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp
          
          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Create 2RDM 2RDM term
          bj = [i+4 for i in bi]
          coeff = 0.5 * sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d2_copy_1 = d2.copy()
          d2_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))

          # Create 'un-permuted' 2RDM 1RDM 1RDM terms
          coeff = 0.5 * sign * (-1)
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Create 'permuted' 2RDM 1RDM 1RDM terms
          coeff = 0.5 * sign
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          bj = [i+4 for i in bj]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Create 1RDM 1RDM 1RDM 1RDM terms
          bj = [i for i in bi]
          bj = [i+4 for i in bj]
          coeff = 0.5 * sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = 0.5 * sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          bj = [i+4 for i in bj]
          coeff = 0.5 * sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = 0.5 * sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

  # Compute terms from the sum involving one 2-body cumulant
  for p in range(4):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(4):
          if r == s:
            continue
          bi = [r,s]
          temp = range(4)
          del temp[max(r,s)]
          del temp[min(r,s)]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp

          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Compute 1RDM 1RDM 2RDM term
          bj = [i for i in bi]
          bj = [i+4 for i in bj]
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d1_copy_0, d1_copy_1, d2_copy_0]))

          # Compute 'un-permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = [i for i in bi]
          bj = [i+4 for i in bj]
          coeff = sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

          # Compute 'permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

  # Compute terms from the sum involving only one particle density matrices
  ti = range(4)
  for bi in makePermutations(4):
    
    # Determine the number of permutations
    n_perm = get_num_perms(ti,bi)

    # Determine the sign
    sign = (-1)**n_perm

    # Compute 1RDM 1RDM 1RDM 1RDM term
    bj = [i+4 for i in bi]
    coeff = sign
    d1_copies = []
    for j in range(4):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab])

  # Combine like terms
  for t in decomp:
    t.tensors.sort()
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdm_to_3rdm_so(d4, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Returns a list of terms composing the decomposition of the 4 body
  reduced density matrix in terms of one, two, and three body rdms"
  """

  # d4 -- Four body rdm to be decomposed
  # d1_aa -- One body rdm
  # d1_bb -- One body rdm
  # d2_aaaa -- Two body rdm
  # d2_bbbb -- Two body rdm
  # d2_abab -- Two body rdm
  # d3_aaaaaa -- Three body rdm
  # d3_bbbbbb -- Three body rdm
  # d3_baabaa -- Three body rdm
  # d3_abbabb -- Three body rdm

  # The types of each index on the input tensors must be the same except for the spin variable.
  # The spin variable must have exactly two possible values.

  # Check input
  if ( not isinstance(d4, tensor) ) or (len(d4.indices) != 8):
    raise TypeError, "d4 must be a tensor with 8 indices"
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Create dummy tensors for unknown-spin rdms
  taken_names = [d4.name, d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name, \
                 d3_aaaaaa.name, d3_bbbbbb.name, d3_baabaa.name, d3_abbabb.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name
  d3 = d3_aaaaaa.copy()
  d3.name = dummy_name

  # Initialize the return value
  decomp = []

  # Compute terms from the sum involving the 3-body cumulant
  for p in range(4):
    ti = range(4)
    del(ti[p])
    ti.insert(0,p)
    for q in range(4):
      bi = range(4)
      del(bi[q])
      bi.insert(0,q)
      for i in range(1,4):
        for j in range(1,4):
          if i != j and ti[i] == bi[j]:
            temp = bi[i]
            bi[i] = bi[j]
            bi[j] = temp

      # Determine the number of index permutations
      n_perm = get_num_perms(ti,bi)

      # Determine the sign
      sign = (-1)**n_perm

      # Create the 1RDM 3RDM term
      bj = [i+4 for i in bi]
      coeff = sign
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
      d3_copy_0 = d3.copy()
      d3_copy_0.indices = [d4.indices[i] for i in (ti[1:4] + bj[1:4])]
      decomp.append(term(coeff, [], [d1_copy_0, d3_copy_0]))

      # Create terms not involving the 3RDM
      d3_decomp = decomp_3rdm_to_2rdm_so(d3_copy_0, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)
      for t in d3_decomp:
        coeff = (-1) * sign * t.numConstant
        decomp.append(term(coeff, [], [d1_copy_0] + t.tensors))
      
  # Compute terms from the sum involving two 2-body cumulants
  for p in range(1):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(r+1,4):
          if s == p or r == q:
            bi = [s,r]
          else:
            bi = [r,s]
          temp = range(4)
          del temp[s]
          del temp[r]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp
          
          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Create 2RDM 2RDM term
          bj = [i+4 for i in bi]
          coeff = sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d2_copy_1 = d2.copy()
          d2_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))

          # Create 'un-permuted' 2RDM 1RDM 1RDM terms
          coeff = sign * (-1)
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Create 'permuted' 2RDM 1RDM 1RDM terms
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          bj = [i+4 for i in bj]
          coeff = sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Create 1RDM 1RDM 1RDM 1RDM terms
          bj = [i for i in bi]
          bj = [i+4 for i in bj]
          coeff = sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          bj = [i+4 for i in bj]
          coeff = sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

  # Compute terms from the sum involving one 2-body cumulant
  for p in range(4):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(4):
          if r == s:
            continue
          bi = [r,s]
          temp = range(4)
          del temp[max(r,s)]
          del temp[min(r,s)]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp

          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Compute 1RDM 1RDM 2RDM term
          bj = [i+4 for i in bi]
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d1_copy_0, d1_copy_1, d2_copy_0]))

          # Compute 'un-permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = [i for i in bi]
          bj = [i+4 for i in bj]
          coeff = sign * (-1)
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

          # Compute 'permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          bj = [i+4 for i in bj]
          coeff = sign
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

  # Compute terms from the sum involving only one particle density matrices
  ti = range(4)
  for bi in makePermutations(4):

    # Determine the number of permutations
    n_perm = get_num_perms(ti,bi)

    # Determine the sign
    sign = (-1)**n_perm

    # Compute 1RDM 1RDM 1RDM 1RDM term
    bj = [i+4 for i in bi]
    coeff = sign
    d1_copies = []
    for j in range(4):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb])

  # Combine like terms
  for t in decomp:
    t.tensors.sort()
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3ops_to_2ops_2rdms_so(inTerms, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Decomposes any 3-particle excitation operators in inTerms into products of 1 and 2
  particle excitation operators and reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose three particle operators and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the operator variable
    op = False

    # Check whether the term has exactly 3 creation and 3 destruction operators
    creCount = 0
    desCount = 0
    for ten in t.tensors:
      if isinstance(ten, creOp):
        creCount += 1
      elif isinstance(ten, desOp):
        desCount += 1

    # Search for the 3-particle operator in the term
    if creCount == 3 and desCount == 3:
      for j in range(len(t.tensors)-1, 4, -1):
        if isinstance(t.tensors[j-0], desOp) and \
           isinstance(t.tensors[j-1], desOp) and \
           isinstance(t.tensors[j-2], desOp) and \
           isinstance(t.tensors[j-3], creOp) and \
           isinstance(t.tensors[j-4], creOp) and \
           isinstance(t.tensors[j-5], creOp):
          op = t.tensors[j-5:j+1]
          pre_term = term(t.numConstant, t.constants, t.tensors[0:j-5])
          post_term = term(1.0, [], t.tensors[j+1:])
          break

    # If the term has no three particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a three particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_3op_to_2op_2rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 3-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4ops_to_2ops_2rdms_so(inTerms, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Decomposes any 4-particle excitation operators in inTerms into products of 1 and 2
  particle excitation operators and reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose four particle operators and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the operator variable
    op = False

    # Check whether the term has exactly 4 creation and 4 destruction operators
    creCount = 0
    desCount = 0
    for ten in t.tensors:
      if isinstance(ten, creOp):
        creCount += 1
      elif isinstance(ten, desOp):
        desCount += 1

    # Search for a 4-particle operator in the term
    if creCount == 4 and desCount == 4:
      for j in range(len(t.tensors)-1, 6, -1):
        if isinstance(t.tensors[j-0], desOp) and \
           isinstance(t.tensors[j-1], desOp) and \
           isinstance(t.tensors[j-2], desOp) and \
           isinstance(t.tensors[j-3], desOp) and \
           isinstance(t.tensors[j-4], creOp) and \
           isinstance(t.tensors[j-5], creOp) and \
           isinstance(t.tensors[j-6], creOp) and \
           isinstance(t.tensors[j-7], creOp):
          op = t.tensors[j-7:j+1]
          pre_term = term(t.numConstant, t.constants, t.tensors[0:j-7])
          post_term = term(1.0, [], t.tensors[j+1:])
          break

    # If the term has no four particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a four particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_4op_to_2op_2rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3ops_to_2ops_3rdms_so(inTerms, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Decomposes any 3-particle excitation operators in inTerms into products of 1 and 2
  particle excitation operators and 1, 2, and 3 particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose three particle operators and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the operator variable
    op = False

    # Check whether the term has exactly 3 creation and 3 destruction operators
    creCount = 0
    desCount = 0
    for ten in t.tensors:
      if isinstance(ten, creOp):
        creCount += 1
      elif isinstance(ten, desOp):
        desCount += 1

    # Search for the 3-particle operator in the term
    if creCount == 3 and desCount == 3:
      for j in range(len(t.tensors)-1, 4, -1):
        if isinstance(t.tensors[j-0], desOp) and \
           isinstance(t.tensors[j-1], desOp) and \
           isinstance(t.tensors[j-2], desOp) and \
           isinstance(t.tensors[j-3], creOp) and \
           isinstance(t.tensors[j-4], creOp) and \
           isinstance(t.tensors[j-5], creOp):
          op = t.tensors[j-5:j+1]
          pre_term = term(t.numConstant, t.constants, t.tensors[0:j-5])
          post_term = term(1.0, [], t.tensors[j+1:])
          break

    # If the term has no three particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a three particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_3op_to_2op_3rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 3-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4ops_to_2ops_3rdms_so(inTerms, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Decomposes any 4-particle excitation operators in inTerms into products of 1 and 2
  particle excitation operators and 1, 2, and 3 particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Initialize index
  i = 0
  opCount = 0

  # One by one, decompose four particle operators and append the result
  # to the list of terms
  while i < len(inTerms):

    # Assign a short name for the ith term
    t = inTerms[i]

    # Check input
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

    # Initialize the operator variable
    op = False

    # Check whether the term has exactly 4 creation and 4 destruction operators
    creCount = 0
    desCount = 0
    for ten in t.tensors:
      if isinstance(ten, creOp):
        creCount += 1
      elif isinstance(ten, desOp):
        desCount += 1

    # Search for a 4-particle operator in the term
    if creCount == 4 and desCount == 4:
      for j in range(len(t.tensors)-1, 6, -1):
        if isinstance(t.tensors[j-0], desOp) and \
           isinstance(t.tensors[j-1], desOp) and \
           isinstance(t.tensors[j-2], desOp) and \
           isinstance(t.tensors[j-3], desOp) and \
           isinstance(t.tensors[j-4], creOp) and \
           isinstance(t.tensors[j-5], creOp) and \
           isinstance(t.tensors[j-6], creOp) and \
           isinstance(t.tensors[j-7], creOp):
          op = t.tensors[j-7:j+1]
          pre_term = term(t.numConstant, t.constants, t.tensors[0:j-7])
          post_term = term(1.0, [], t.tensors[j+1:])
          break

    # If the term has no four particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a four particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_4op_to_2op_3rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3op_to_2op_2rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Returns the list of terms composing the decomposition of the 3 body operator op
  in terms of one and two body operators"
  """

  # Prepare error message
  TypeErrorMessage = "op must be a list of three creOp objects followed by three desOp objects."

  # Check input
  if type(op) != type([]) or len(op) != 6 or \
     not (isinstance(op[0], creOp) and isinstance(op[1], creOp) and isinstance(op[2], creOp) and \
          isinstance(op[3], desOp) and isinstance(op[4], desOp) and isinstance(op[5], desOp)):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Convert to density matrix ordering
  (op[3],op[4],op[5]) = (op[5],op[4],op[3])

  # Initialize the return value
  decomp = []

  # Loop over the nine index permutations used in the decomposition
  for p0 in range(3):

    # Compute the top indices
    temp = range(3)
    del temp[p0]
    p1 = temp[0]
    p2 = temp[1]
    ti = [p0, p1, p2]

    for q0 in range(3):

      # Compute the bottom indices, remembering that the second and third index pairs
      # should retain original pairing as much as possible.
      temp = range(3)
      del temp[q0]
      q1 = temp[0]
      q2 = temp[1]
      if q1 == p2 or q2 == p1:
        bi = [q0, q2, q1]
      else:
        bi = [q0, q1, q2]

      # Determine the number of permutations from original pairing
      n_perms = get_num_perms(ti,bi)

      # Compute sign
      sign = (-1)**n_perms

      # Shift bi to match the bottom indices of op
      bj = [i+3 for i in bi]

      # operators below are restored to operator ordering

      # Create 1RDM 2op term
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [op[ti[0]].indices[0], op[bj[0]].indices[0]]
      opList = [op[i] for i in ti[1:3]+bj[2:3]+bj[1:2]]
      decomp.append( term( sign, [], [d1_copy_0] + opList ) )

      # Create un-permuted 1RDM 1RDM 1op terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      opList = [op[i] for i in ti[2:3]+bj[2:3]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1] + opList ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      opList = [op[i] for i in ti[1:2]+bj[1:2]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1] + opList ) )

      # Create permuted 1RDM 1RDM 1op terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      opList = [op[i] for i in ti[2:3]+bj[1:2]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1] + opList ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      opList = [op[i] for i in ti[1:2]+bj[2:3]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1] + opList ) )

      # Create un-permuted 1RDM 1RDM 1RDM terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )

      # Create permuted 1RDM 1RDM 1RDM terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )

      # Create 1RDM 2RDM term from extended-normal-ordered 2-body op
      d2_copy = d2.copy()
      d2_copy.indices = [op[i].indices[0] for i in ti[1:3]+bj[1:3]]
      decomp.append( term(-sign, [], [d1_copy_0, d2_copy] ) )

      # Create 2RDM 1op term
      d2_copy = d2.copy()
      d2_copy.indices = [op[i].indices[0] for i in ti[1:3] + bj[1:3]]
      opList = [op[i] for i in ti[0:1]+bj[0:1]]
      decomp.append( term( sign, [], [d2_copy] + opList ) )

      # Create 1RDM 2RDM term from extended-normal-ordered 2-body op
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [op[ti[0]].indices[0], op[bj[0]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d2_copy] ) )

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab])

  # Compute terms from 3RDM decomposition
  taken_names = [d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name]
  d3name = 'd3_temp'
  while d3name in taken_names:
    d3name += '_'
  d3 = tensor(d3name, [i.indices[0] for i in op], [])
  decomp.extend(decomp_3rdm_to_2rdm_so(d3, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab))

  # Combine like terms
  for t in decomp:
    i = 0
    while i < len(t.tensors) and \
          not isinstance(t.tensors[i], creOp) and \
          not isinstance(t.tensors[i], desOp):
      i += 1
    temp = t.tensors[0:i]
    for ten in temp:
      t.scale(ten.sortIndeces())
    temp.sort()
    temp.extend(t.tensors[i:])
    t.tensors = temp
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Convert back to operator ordering
  (op[3],op[4],op[5]) = (op[5],op[4],op[3])

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4op_to_2op_2rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab):
  """
  Returns the list of terms composing the decomposition of the 4 body operator op
  in terms of one and two body operators"
  """

  # Prepare error message
  TypeErrorMessage = "op must be a list of four creOp objects followed by four desOp objects."

  # Check input
  if type(op) != type([]) or len(op) != 8 or \
     not (isinstance(op[0], creOp) and isinstance(op[1], creOp) and \
          isinstance(op[2], creOp) and isinstance(op[3], creOp) and \
          isinstance(op[4], desOp) and isinstance(op[5], desOp) and \
          isinstance(op[6], desOp) and isinstance(op[7], desOp)):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Convert to density matrix ordering
  (op[4],op[5],op[6],op[7]) = (op[7],op[6],op[5],op[4])

  # Compute a list of the operator's indices
  indexList = [o.indices[0] for o in op]

  # Initialize the return value
  decomp = []

#  raise RuntimeError, "decomp_4op_to_2op_2rdm_so not yet implemented"

  # Compute terms from the extended-normal-ordered 2-body operator
  for p in range(4):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(r+1,4):
          if s == p or r == q:
            bi = [s,r]
          else:
            bi = [r,s]
          temp = range(4)
          del temp[s]
          del temp[r]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp
          
          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Adjust bottom indices to match indexing in op
          bj = [i+4 for i in bi]

          # Compute 2RDM 2op term
          coeff = sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [indexList[i] for i in (ti[0:2] + bj[0:2])]
          opList = [creOp(indexList[ti[2]]), creOp(indexList[ti[3]]), desOp(indexList[bj[3]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0] + opList))

          # Compute un-permuted 2RDM 1RDM 1op terms
          coeff = -sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          opList = [creOp(indexList[ti[3]]), desOp(indexList[bj[3]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          opList = [creOp(indexList[ti[2]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))

          # Compute permuted 2RDM 1RDM 1op terms
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          opList = [creOp(indexList[ti[3]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          opList = [creOp(indexList[ti[2]]), desOp(indexList[bj[3]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))

          # Compute un-permuted 2RDM 1RDM 1RDM terms
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Compute permuted 2RDM 1RDM 1RDM terms
          coeff = -sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Compute 2RDM 2RDM term
          coeff = -sign
          d2_copy_1 = d2.copy()
          d2_copy_1.indices = [indexList[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))

  # Compute terms from the extended-normal-ordered 1-body operator
  for p in range(4):
    ti = [p] + range(4)
    del(ti[p+1])
    for q in range(4):
      bi = [q] + range(4)
      del(bi[q+1])
      if bi[2] == ti[1] or bi[1] == ti[2]:
        (bi[1],bi[2]) = (bi[2],bi[1])
      if bi[3] == ti[1] or bi[1] == ti[3]:
        (bi[1],bi[3]) = (bi[3],bi[1])
      if bi[3] == ti[2] or bi[2] == ti[3]:
        (bi[2],bi[3]) = (bi[3],bi[2])

      # Get the number of permutations
      n_perm = get_num_perms(ti,bi)

      # Compute the sign
      sign = (-1)**n_perm

      # Adjust bottom indices to match indexing in op
      bj = [i+4 for i in bi]

      # Compute the decomposition of the 3rdm
      d3 = tensor('d3_dummy', [indexList[i] for i in (ti[1:4] + bj[1:4])], [])
      d3_decomp = decomp_3rdm_to_2rdm_so(d3, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab)

      # Compute terms involving 1op
      opList = [creOp(indexList[ti[0]]), desOp(indexList[bj[0]])]
      for t in d3_decomp:
        coeff = sign * t.numConstant
        decomp.append(term(coeff, [], t.tensors + opList))

      # Compute terms involving all RDMs
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [indexList[i] for i in (ti[0:1] + bj[0:1])]
      for t in d3_decomp:
        coeff = -sign * t.numConstant
        decomp.append(term(coeff, [], [d1_copy_0] + t.tensors))

  # Compute terms from the decomposed 4-body density matrix
  d4 = tensor('d4_dummy', indexList, [])
  decomp.extend(decomp_4rdm_to_2rdm_so(d4, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab])

  # Combine like terms
  for t in decomp:
    i = 0
    while i < len(t.tensors) and \
          not isinstance(t.tensors[i], creOp) and \
          not isinstance(t.tensors[i], desOp):
      i += 1
    temp = t.tensors[0:i]
    for ten in temp:
      t.scale(ten.sortIndeces())
    temp.sort()
    temp.extend(t.tensors[i:])
    t.tensors = temp
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Convert back to operator ordering
  (op[4],op[5],op[6],op[7]) = (op[7],op[6],op[5],op[4])

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3op_to_2op_3rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Returns the list of terms composing the decomposition of the 3 body operator op
  in terms of one and two body operators, using the exact 3rdm"
  """

  # Prepare error message
  TypeErrorMessage = "op must be a list of three creOp objects followed by three desOp objects."

  # Check input
  if type(op) != type([]) or len(op) != 6 or \
     not (isinstance(op[0], creOp) and isinstance(op[1], creOp) and isinstance(op[2], creOp) and \
          isinstance(op[3], desOp) and isinstance(op[4], desOp) and isinstance(op[5], desOp)):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Convert to density matrix ordering
  (op[3],op[4],op[5]) = (op[5],op[4],op[3])

  # Initialize the return value
  decomp = []

  # Loop over the nine index permutations used in the decomposition
  for p0 in range(3):

    # Compute the top indices
    temp = range(3)
    del temp[p0]
    p1 = temp[0]
    p2 = temp[1]
    ti = [p0, p1, p2]

    for q0 in range(3):

      # Compute the bottom indices, remembering that the second and third index pairs
      # should retain original pairing as much as possible.
      temp = range(3)
      del temp[q0]
      q1 = temp[0]
      q2 = temp[1]
      if q1 == p2 or q2 == p1:
        bi = [q0, q2, q1]
      else:
        bi = [q0, q1, q2]

      # Determine the number of permutations from original pairing
      n_perms = get_num_perms(ti,bi)

      # Compute sign
      sign = (-1)**n_perms

      # Shift bi to match the bottom indices of op
      bj = [i+3 for i in bi]

      # operators below are restored to operator ordering

      # Create 1RDM 2op term
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [op[ti[0]].indices[0], op[bj[0]].indices[0]]
      opList = [op[i] for i in ti[1:3]+bj[2:3]+bj[1:2]]
      decomp.append( term( sign, [], [d1_copy_0] + opList ) )

      # Create un-permuted 1RDM 1RDM 1op terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      opList = [op[i] for i in ti[2:3]+bj[2:3]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1] + opList ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      opList = [op[i] for i in ti[1:2]+bj[1:2]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1] + opList ) )

      # Create permuted 1RDM 1RDM 1op terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      opList = [op[i] for i in ti[2:3]+bj[1:2]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1] + opList ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      opList = [op[i] for i in ti[1:2]+bj[2:3]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1] + opList ) )

      # Create un-permuted 1RDM 1RDM 1RDM terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[2]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[1]].indices[0], op[bj[1]].indices[0]]
      decomp.append( term( sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )

      # Create permuted 1RDM 1RDM 1RDM terms
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )
      d1_copy_1 = d1.copy()
      d1_copy_1.indices = [op[ti[2]].indices[0], op[bj[1]].indices[0]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op[ti[1]].indices[0], op[bj[2]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d1_copy_1, d1_copy_2] ) )

      # Create 1RDM 2RDM term from extended-normal-ordered 2-body op
      d2_copy = d2.copy()
      d2_copy.indices = [op[i].indices[0] for i in ti[1:3]+bj[1:3]]
      decomp.append( term(-sign, [], [d1_copy_0, d2_copy] ) )

      # Create 2RDM 1op term
      d2_copy = d2.copy()
      d2_copy.indices = [op[i].indices[0] for i in ti[1:3] + bj[1:3]]
      opList = [op[i] for i in ti[0:1]+bj[0:1]]
      decomp.append( term( sign, [], [d2_copy] + opList ) )

      # Create 1RDM 2RDM term from extended-normal-ordered 2-body op
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [op[ti[0]].indices[0], op[bj[0]].indices[0]]
      decomp.append( term(-sign, [], [d1_copy_0, d2_copy] ) )

  # Compute the 3RDM term
  decomp.append(term(1.0, [], [tensor(dummy_name, [i.indices[0] for i in op], [])]))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb])

  # Combine like terms
  for t in decomp:
    i = 0
    while i < len(t.tensors) and \
          not isinstance(t.tensors[i], creOp) and \
          not isinstance(t.tensors[i], desOp):
      i += 1
    temp = t.tensors[0:i]
    for ten in temp:
      t.scale(ten.sortIndeces())
    temp.sort()
    temp.extend(t.tensors[i:])
    t.tensors = temp
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Convert back to operator ordering
  (op[3],op[4],op[5]) = (op[5],op[4],op[3])

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4op_to_2op_3rdm_so(op, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb):
  """
  Returns the list of terms composing the decomposition of the 4 body operator op
  in terms of one and two body operators, using the exact 3rdm"
  """

  # Prepare error message
  TypeErrorMessage = "op must be a list of four creOp objects followed by four desOp objects."

  # Check input
  if type(op) != type([]) or len(op) != 8 or \
     not (isinstance(op[0], creOp) and isinstance(op[1], creOp) and \
          isinstance(op[2], creOp) and isinstance(op[3], creOp) and \
          isinstance(op[4], desOp) and isinstance(op[5], desOp) and \
          isinstance(op[6], desOp) and isinstance(op[7], desOp)):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1_aa, tensor) ) or (len(d1_aa.indices) != 2):
    raise TypeError, "d1_aa must be a tensor with 2 indices"
  if ( not isinstance(d1_bb, tensor) ) or (len(d1_bb.indices) != 2):
    raise TypeError, "d1_bb must be a tensor with 2 indices"
  if ( not isinstance(d2_aaaa, tensor) ) or (len(d2_aaaa.indices) != 4):
    raise TypeError, "d2_aaaa must be a tensor with 4 indices"
  if ( not isinstance(d2_bbbb, tensor) ) or (len(d2_bbbb.indices) != 4):
    raise TypeError, "d2_bbbb must be a tensor with 4 indices"
  if ( not isinstance(d2_abab, tensor) ) or (len(d2_abab.indices) != 4):
    raise TypeError, "d2_abab must be a tensor with 4 indices"
  if ( not isinstance(d3_aaaaaa, tensor) ) or (len(d3_aaaaaa.indices) != 6):
    raise TypeError, "d3_aaaaaa must be a tensor with 6 indices"
  if ( not isinstance(d3_bbbbbb, tensor) ) or (len(d3_bbbbbb.indices) != 6):
    raise TypeError, "d3_bbbbbb must be a tensor with 6 indices"
  if ( not isinstance(d3_baabaa, tensor) ) or (len(d3_baabaa.indices) != 6):
    raise TypeError, "d3_baabaa must be a tensor with 6 indices"
  if ( not isinstance(d3_abbabb, tensor) ) or (len(d3_abbabb.indices) != 6):
    raise TypeError, "d3_abbabb must be a tensor with 6 indices"

  # Create dummy tensors for unknown-spin 1 and 2 rdms
  taken_names = [d1_aa.name, d1_bb.name, d2_aaaa.name, d2_bbbb.name, d2_abab.name] 
  dummy_name = 'rdm_dummy'
  while dummy_name in taken_names:
    dummy_name += '_'
  d1 = d1_aa.copy()
  d1.name = dummy_name
  d2 = d2_aaaa.copy()
  d2.name = dummy_name

  # Convert to density matrix ordering
  (op[4],op[5],op[6],op[7]) = (op[7],op[6],op[5],op[4])

  # Compute a list of the operator's indices
  indexList = [o.indices[0] for o in op]

  # Initialize the return value
  decomp = []

#  raise RuntimeError, "decomp_4op_to_2op_2rdm_so not yet implemented"

  # Compute terms from the extended-normal-ordered 2-body operator
  for p in range(4):
    for q in range(p+1,4):
      ti = [p,q]
      for i in range(4):
        if not (i in ti):
          ti.append(i)
      for r in range(4):
        for s in range(r+1,4):
          if s == p or r == q:
            bi = [s,r]
          else:
            bi = [r,s]
          temp = range(4)
          del temp[s]
          del temp[r]
          if temp[0] == ti[3] or temp[1] == ti[2]:
            temp = [temp[1], temp[0]]
          bi += temp
          
          # Determine the number of index permutations
          n_perm = get_num_perms(ti,bi)

          # Determine the sign
          sign = (-1)**n_perm

          # Adjust bottom indices to match indexing in op
          bj = [i+4 for i in bi]

          # Compute 2RDM 2op term
          coeff = sign
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [indexList[i] for i in (ti[0:2] + bj[0:2])]
          opList = [creOp(indexList[ti[2]]), creOp(indexList[ti[3]]), desOp(indexList[bj[3]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0] + opList))

          # Compute un-permuted 2RDM 1RDM 1op terms
          coeff = -sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          opList = [creOp(indexList[ti[3]]), desOp(indexList[bj[3]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          opList = [creOp(indexList[ti[2]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))

          # Compute permuted 2RDM 1RDM 1op terms
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          opList = [creOp(indexList[ti[3]]), desOp(indexList[bj[2]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          opList = [creOp(indexList[ti[2]]), desOp(indexList[bj[3]])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0] + opList))

          # Compute un-permuted 2RDM 1RDM 1RDM terms
          coeff = sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[3:4])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[2:3] + bj[2:3])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Compute permuted 2RDM 1RDM 1RDM terms
          coeff = -sign
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [indexList[i] for i in (ti[3:4] + bj[2:3])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [indexList[i] for i in (ti[2:3] + bj[3:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))

          # Compute 2RDM 2RDM term
          coeff = -sign
          d2_copy_1 = d2.copy()
          d2_copy_1.indices = [indexList[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))

  # Compute terms from the extended-normal-ordered 1-body operator
  for p in range(4):
    ti = [p] + range(4)
    del(ti[p+1])
    for q in range(4):
      bi = [q] + range(4)
      del(bi[q+1])
      if bi[2] == ti[1] or bi[1] == ti[2]:
        (bi[1],bi[2]) = (bi[2],bi[1])
      if bi[3] == ti[1] or bi[1] == ti[3]:
        (bi[1],bi[3]) = (bi[3],bi[1])
      if bi[3] == ti[2] or bi[2] == ti[3]:
        (bi[2],bi[3]) = (bi[3],bi[2])

      # Get the number of permutations
      n_perm = get_num_perms(ti,bi)

      # Compute the sign
      sign = (-1)**n_perm

      # Adjust bottom indices to match indexing in op
      bj = [i+4 for i in bi]

      # Compute term involving 1op
      coeff = sign
      d3 = tensor(dummy_name, [indexList[i] for i in (ti[1:4] + bj[1:4])], [])
      opList = [creOp(indexList[ti[0]]), desOp(indexList[bj[0]])]
      decomp.append(term(coeff, [], [d3] + opList))

      # Compute term involving 1RDM
      coeff = -sign
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [indexList[i] for i in (ti[0:1] + bj[0:1])]
      decomp.append(term(coeff, [], [d3, d1_copy_0]))

  # Compute terms from the decomposed 4-body density matrix
  d4 = tensor('d4_dummy', indexList, [])
  decomp.extend(decomp_4rdm_to_3rdm_so(d4, d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb))

  # Assign RDM types
  assign_rdm_types(decomp, dummy_name, [d1_aa, d1_bb, d2_aaaa, d2_bbbb, d2_abab, d3_aaaaaa, d3_bbbbbb, d3_baabaa, d3_abbabb])

  # Combine like terms
  for t in decomp:
    i = 0
    while i < len(t.tensors) and \
          not isinstance(t.tensors[i], creOp) and \
          not isinstance(t.tensors[i], desOp):
      i += 1
    temp = t.tensors[0:i]
    for ten in temp:
      t.scale(ten.sortIndeces())
    temp.sort()
    temp.extend(t.tensors[i:])
    t.tensors = temp
  decomp.sort()
  i = 0
  while i < len(decomp)-1:
    if decomp[i].tensors == decomp[i+1].tensors:
      decomp[i+1].numConstant += decomp[i].numConstant
      del(decomp[i])
    else:
      i += 1
  termChop(decomp)

  # Convert back to operator ordering
  (op[4],op[5],op[6],op[7]) = (op[7],op[6],op[5],op[4])

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


