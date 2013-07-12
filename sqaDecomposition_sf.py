#    file:  sqaDecomposition_sf.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Contains functions for the various decompositions
#           used in Canonical Transformation theory with
#           respect to spin-free operators.
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
from sqaTensor import tensor, sfExOp
from sqaTerm import term, multiplyTerms, termChop
from sqaMisc import makePermutations, get_num_perms
from sqaOptions import options


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3rdms_to_2rdms_sf(inTerms, d3name, d1, d2, indexOrder = '1212'):
  """
  Decomposes any 3-particle spin free RDMs in inTerms into products of 1 and 2
  particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"

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
      decomp = decomp_3rdm_to_2rdm_sf(rdm, d1, d2, indexOrder)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 3-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdms_to_2rdms_sf(inTerms, d4name, d1, d2, d2_hom, d2_het):
  """
  Decomposes any 4-particle spin free RDMs in inTerms into products of 1 and 2
  particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"
  if ( not isinstance(d2_hom, tensor) ) or (len(d2_hom.indices) != 4):
    raise TypeError, "d2_hom must be a tensor with 4 indices"
  if ( not isinstance(d2_het, tensor) ) or (len(d2_het.indices) != 4):
    raise TypeError, "d2_het must be a tensor with 4 indices"

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
      decomp = decomp_4rdm_to_2rdm_sf(rdm, d1, d2, d2_hom, d2_het)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdms_to_3rdms_sf(inTerms, d4name, d1, d2, d2_hom, d2_het, d3):
  """
  Decomposes any 4-particle spin free RDMs in inTerms into products of 1, 2, and 3
  particle reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"
  if ( not isinstance(d2_hom, tensor) ) or (len(d2_hom.indices) != 4):
    raise TypeError, "d2_hom must be a tensor with 4 indices"
  if ( not isinstance(d2_het, tensor) ) or (len(d2_het.indices) != 4):
    raise TypeError, "d2_het must be a tensor with 4 indices"
  if ( not isinstance(d3, tensor) ) or (len(d3.indices) != 6):
    raise TypeError, "d3 must be a tensor with 6 indices"

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
      decomp = decomp_4rdm_to_3rdm_sf(rdm, d1, d2, d2_hom, d2_het, d3)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body RDMs' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3rdm_to_2rdm_sf(d3, d1, d2, indexOrder = '1212'):
  """
  Returns the list of terms composing the decomposition of the 3 body spin free rdm d3
  in terms of one and two body rdms"
  """

  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"
  if ( not isinstance(d3, tensor) ) or (len(d3.indices) != 6):
    raise TypeError, "d3 must be a tensor with 6 indices"

  # sort the 3-body RDM's indices into 1212 ordering
  d3_indices = []
  if   indexOrder == '1212':
    d3_indices.extend(d3.indices)
  elif indexOrder == '1122':
    for i in range(3):
      d3_indices.append(d3.indices[2*i])
    for i in range(3):
      d3_indices.append(d3.indices[2*i+1])
  else:
    raise ValueError, "unexpected indexOrder parameter:  %s" %(indexOrder)

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
      bj = [i for i in bi]
      n_perm = get_num_perms(ti,bj)
      coeff = sign * (0.5)**n_perm
      bj = [i+3 for i in bj]
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [d3_indices[i] for i in (ti[0:1] + bj[0:1])]
      d2_copy_0 = d2.copy()
      d2_copy_0.indices = [d3_indices[i] for i in (ti[1:3] + bj[1:3])]
      decomp.append(term(coeff, [], [d1_copy_0, d2_copy_0]))

      # create the 1RDM 1RDM 1RDM terms
      bj = [i for i in bi]
      n_perm = get_num_perms(ti,bj)
      coeff = sign * (-1) * (0.5)**n_perm
      bj = [i+3 for i in bj]
      d1_copies = []
      for j in range(3):
        d1_copies.append(d1.copy())
        d1_copies[-1].indices = [d3_indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
      decomp.append(term(coeff, [], d1_copies))
      bj = bi[0:1] + bi[2:3] + bi[1:2]
      n_perm = get_num_perms(ti,bj)
      coeff = sign * (0.5)**n_perm
      bj = [i+3 for i in bj]
      d1_copies = []
      for j in range(3):
        d1_copies.append(d1.copy())
        d1_copies[-1].indices = [d3_indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
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
    coeff = sign * (0.5)**n_perm
    d1_copies = []
    for j in range(3):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d3_indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

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

  # if original ordering was 1122, convert back to it
  if indexOrder == '1122':
    for t in decomp:
      for ten in t.tensors:
        if ten.name == d2.name:
          ten.indices = [ten.indices[0], ten.indices[2], ten.indices[1], ten.indices[3]]

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4rdm_to_2rdm_sf(d4, d1, d2, d2_hom, d2_het):
  """
  Returns the list of terms composing the decomposition of the 4 body spin free rdm d4
  in terms of one and two body rdms"
  """

  # This decomposition assumes the density matrices correspond to a spin averaged ensemble state.

  # d4 -- Spin free 4rdm to be decomposed
  # d1 -- One body spin free rdm
  # d2 -- Two body spin free rdm
  # d2_hom -- Sum of all alpha and all beta two body rdms.  The "homogeneous spin" 2rdm.
  # d2_het -- Sum of abab and baba two body rdms.  The "heterogeneous spin" 2rdm.

  # Note that d2 = d2_hom + d2_het

  # Check input
  if ( not isinstance(d4, tensor) ) or (len(d4.indices) != 8):
    raise TypeError, "d4 must be a tensor with 8 indices"
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"
  if ( not isinstance(d2_hom, tensor) ) or (len(d2_hom.indices) != 4):
    raise TypeError, "d2_hom must be a tensor with 4 indices"
  if ( not isinstance(d2_het, tensor) ) or (len(d2_het.indices) != 4):
    raise TypeError, "d2_het must be a tensor with 4 indices"

  # Initialize the return value
  decomp = []

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
          if n_perm < 2:
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_copy_1 = d2.copy()
            d2_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))
          elif n_perm == 2:
            coeff = sign * 0.5
            d2_hom_copy_0 = d2_hom.copy()
            d2_hom_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_hom_copy_1 = d2_hom.copy()
            d2_hom_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_hom_copy_0, d2_hom_copy_1]))
            d2_het_copy_0 = d2_het.copy()
            d2_het_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_het_copy_1 = d2_het.copy()
            d2_het_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_het_copy_0, d2_het_copy_1]))
            d2_het_copy_0 = d2_het.copy()
            d2_het_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + [bj[1], bj[0]])]
            d2_het_copy_1 = d2_het.copy()
            d2_het_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + [bj[3], bj[2]])]
            decomp.append(term(coeff, [], [d2_het_copy_0, d2_het_copy_1]))
          else:
            raise RuntimeError, "Did not expect a term with %i permutations" %n_perm

          # Create 'un-permuted' 2RDM 1RDM 1RDM terms
          coeff = sign * (-1) * (0.5)**n_perm 
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
          if n_perm < 2:
            bj = bi[0:2] + bi[3:4] + bi[2:3]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d1_copy_0 = d1.copy()
            d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
            d1_copy_1 = d1.copy()
            d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
            decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
            bj = bi[1:2] + bi[0:1] + bi[2:4]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            d1_copy_0 = d1.copy()
            d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
            d1_copy_1 = d1.copy()
            d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
            decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          elif n_perm == 2:
            bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (-1) * (0.5)**n_perm
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
          else:
            raise RuntimeError, "Did not expect a term with %i permutations" %n_perm

          # Create 1RDM 1RDM 1RDM 1RDM terms
          bj = [i for i in bi]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
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
          coeff = 1
          bj = [i for i in bi]
          pairs_rdm1 = [ [ti[i],bj[i]] for i in range(2)]
          pairs_rdm2 = [ [bj[i],ti[i]] for i in range(2,4)]
          for pair in pairs_rdm2:
            if pair[0] == pair[1]:
              pairs_rdm2 = []
              break
          for pair in pairs_rdm2:
            if not (pair in pairs_rdm1):
              #print "swapping bj[2] and bj[3]"
              (bj[3],bj[2]) = (bj[2],bj[3])
              coeff *= -1
              break
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff *= sign * (0.5)**n_perm
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d1_copy_0, d1_copy_1, d2_copy_0]))

          # Compute 'un-permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = [i for i in bi]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

          # Compute 'permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
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
    coeff = sign * (0.5)**n_perm
    d1_copies = []
    for j in range(4):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

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


def decomp_4rdm_to_3rdm_sf(d4, d1, d2, d2_hom, d2_het, d3):
  """
  Returns the list of terms composing the decomposition of the 4 body spin free rdm d4
  in terms of one, two, and three body rdms"
  """

  # This decomposition assumes the density matrices correspond to a spin averaged ensemble state.

  # d4 -- Spin free 4rdm to be decomposed
  # d1 -- One body spin free rdm
  # d2 -- Two body spin free rdm
  # d2_hom -- Sum of all alpha and all beta two body rdms.  The "homogeneous spin" 2rdm.
  # d2_het -- Sum of abab and baba two body rdms.  The "heterogeneous spin" 2rdm.
  # d3 -- Three body spin free rdm

  # Note that d2 = d2_hom + d2_het

  # Check input
  if ( not isinstance(d4, tensor) ) or (len(d4.indices) != 8):
    raise TypeError, "d4 must be a tensor with 8 indices"
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"
  if ( not isinstance(d2_hom, tensor) ) or (len(d2_hom.indices) != 4):
    raise TypeError, "d2_hom must be a tensor with 4 indices"
  if ( not isinstance(d2_het, tensor) ) or (len(d2_het.indices) != 4):
    raise TypeError, "d2_het must be a tensor with 4 indices"
  if ( not isinstance(d3, tensor) ) or (len(d3.indices) != 6):
    raise TypeError, "d3 must be a tensor with 6 indices"

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
      coeff = sign * (0.5)**n_perm
      d1_copy_0 = d1.copy()
      d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
      d3_copy_0 = d3.copy()
      d3_copy_0.indices = [d4.indices[i] for i in (ti[1:4] + bj[1:4])]
      decomp.append(term(coeff, [], [d1_copy_0, d3_copy_0]))

      # Create terms not involving the 3RDM
      d3_decomp = decomp_3rdm_to_2rdm_sf(d3_copy_0, d1, d2)
      for t in d3_decomp:
        coeff = (-1) * sign * (0.5)**n_perm * t.numConstant
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
          if n_perm < 2:
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_copy_1 = d2.copy()
            d2_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_copy_0, d2_copy_1]))
          elif n_perm == 2:
            coeff = sign * 0.5
            d2_hom_copy_0 = d2_hom.copy()
            d2_hom_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_hom_copy_1 = d2_hom.copy()
            d2_hom_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_hom_copy_0, d2_hom_copy_1]))
            d2_het_copy_0 = d2_het.copy()
            d2_het_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d2_het_copy_1 = d2_het.copy()
            d2_het_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            decomp.append(term(coeff, [], [d2_het_copy_0, d2_het_copy_1]))
            d2_het_copy_0 = d2_het.copy()
            d2_het_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + [bj[1], bj[0]])]
            d2_het_copy_1 = d2_het.copy()
            d2_het_copy_1.indices = [d4.indices[i] for i in (ti[2:4] + [bj[3], bj[2]])]
            decomp.append(term(coeff, [], [d2_het_copy_0, d2_het_copy_1]))
          else:
            raise RuntimeError, "Did not expect a term with %i permutations" %n_perm

          # Create 'un-permuted' 2RDM 1RDM 1RDM terms
          coeff = sign * (-1) * (0.5)**n_perm 
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
          if n_perm < 2:
            bj = bi[0:2] + bi[3:4] + bi[2:3]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[0:2] + bj[0:2])]
            d1_copy_0 = d1.copy()
            d1_copy_0.indices = [d4.indices[i] for i in (ti[2:3] + bj[2:3])]
            d1_copy_1 = d1.copy()
            d1_copy_1.indices = [d4.indices[i] for i in (ti[3:4] + bj[3:4])]
            decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
            bj = bi[1:2] + bi[0:1] + bi[2:4]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (0.5)**n_perm
            d2_copy_0 = d2.copy()
            d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
            d1_copy_0 = d1.copy()
            d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
            d1_copy_1 = d1.copy()
            d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
            decomp.append(term(coeff, [], [d2_copy_0, d1_copy_0, d1_copy_1]))
          elif n_perm == 2:
            bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
            n_perm = get_num_perms(ti,bj)
            bj = [i+4 for i in bj]
            coeff = sign * (-1) * (0.5)**n_perm
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
          else:
            raise RuntimeError, "Did not expect a term with %i permutations" %n_perm

          # Create 1RDM 1RDM 1RDM 1RDM terms
          bj = [i for i in bi]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[2:4]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))
          bj = bi[1:2] + bi[0:1] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
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
          coeff = 1
          bj = [i for i in bi]
          pairs_rdm1 = [ [ti[i],bj[i]] for i in range(2)]
          pairs_rdm2 = [ [bj[i],ti[i]] for i in range(2,4)]
          for pair in pairs_rdm2:
            if pair[0] == pair[1]:
              pairs_rdm2 = []
              break
          for pair in pairs_rdm2:
            if not (pair in pairs_rdm1):
              #print "swapping bj[2] and bj[3]"
              (bj[3],bj[2]) = (bj[2],bj[3])
              coeff *= -1
              break
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff *= sign * (0.5)**n_perm
          d1_copy_0 = d1.copy()
          d1_copy_0.indices = [d4.indices[i] for i in (ti[0:1] + bj[0:1])]
          d1_copy_1 = d1.copy()
          d1_copy_1.indices = [d4.indices[i] for i in (ti[1:2] + bj[1:2])]
          d2_copy_0 = d2.copy()
          d2_copy_0.indices = [d4.indices[i] for i in (ti[2:4] + bj[2:4])]
          decomp.append(term(coeff, [], [d1_copy_0, d1_copy_1, d2_copy_0]))

          # Compute 'un-permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = [i for i in bi]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (-1) * (0.5)**n_perm
          d1_copies = []
          for j in range(4):
            d1_copies.append(d1.copy())
            d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
          decomp.append(term(coeff, [], d1_copies))

          # Compute 'permuted' 1RDM 1RDM 1RDM 1RDM term
          bj = bi[0:2] + bi[3:4] + bi[2:3]
          n_perm = get_num_perms(ti,bj)
          bj = [i+4 for i in bj]
          coeff = sign * (0.5)**n_perm
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
    coeff = sign * (0.5)**n_perm
    d1_copies = []
    for j in range(4):
      d1_copies.append(d1.copy())
      d1_copies[-1].indices = [d4.indices[i] for i in (ti[j:j+1] + bj[j:j+1])]
    decomp.append(term(coeff, [], d1_copies))

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


def decomp_3ops_to_2ops_2rdms_sf(inTerms, d1, d2, indexOrder = '1212'):
  """
  Decomposes any 3-particle sfExOps in inTerms into products of 1 and 2
  particle sfExOps and reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"

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

    # Search for a 3-particle operator in the term
    for j in range(len(t.tensors)-1, -1, -1):
      if isinstance(t.tensors[j], sfExOp) and t.tensors[j].order == 3:
        op = t.tensors[j]
        pre_term = term(t.numConstant, t.constants, t.tensors[0:j])
        post_term = term(1.0, [], t.tensors[j+1:])
        break

    # If the term has no three particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a three particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_3op_to_2op_2rdm_sf(op, d1, d2, indexOrder)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 3-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4ops_to_2ops_2rdms_sf(inTerms, d1, d2):
  """
  Decomposes any 4-particle sfExOps in inTerms into products of 1 and 2
  particle sfExOps and reduced density mattrices.
  """

  # Prepare error message
  TypeErrorMessage = "inTerms must be a list of term objects"

  # Check input
  if type(inTerms) != type([]):
    raise TypeError, TypeErrorMessage
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"

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

    # Search for a 4-particle operator in the term
    for j in range(len(t.tensors)-1, -1, -1):
      if isinstance(t.tensors[j], sfExOp) and t.tensors[j].order == 3:
        op = t.tensors[j]
        pre_term = term(t.numConstant, t.constants, t.tensors[0:j])
        post_term = term(1.0, [], t.tensors[j+1:])
        break

    # If the term has no four particle operators, incrment the index
    if op is False:
      i += 1

    # If the term has a four particle operator, decompose it
    else:
      opCount += 1
      decomp = decomp_4op_to_2op_2rdm_sf(op, d1, d2)
      for s in decomp:
        inTerms.append(pre_term.copy())
        inTerms[-1] = multiplyTerms(inTerms[-1], s)
        inTerms[-1] = multiplyTerms(inTerms[-1], post_term)
      del(inTerms[i])

  if options.verbose:
    print 'decomposed %i 4-body operators' %(opCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_3op_to_2op_2rdm_sf(op, d1, d2, indexOrder = '1212'):
  """
  Returns the list of terms composing the decomposition of the 3 body sfExOp op
  in terms of one and two body operators"
  """

  # Check input
  if ( not isinstance(op, sfExOp) ) or (op.order != 3):
    raise TypeError, "op must be a 3-particle sfExOp"
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"

  if indexOrder != '1212' and indexOrder != '1122':
    raise ValueError, "indexOrder must be '1212' or '1122'"

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

      # Shift bi to match the bottom indices of op
      bi = [i+op.order for i in bi]

      # Create 1RDM 2op term
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[0]], op.indices[bi[0]]]
      indexList = [op.indices[i] for i in ti[1:]] + [op.indices[i] for i in bi[1:]]
      decomp.append( term( (-0.5)**n_perms, [], [d1_copy, sfExOp(indexList)] ) )

      # Create 2RDM 1op term
      d2_copy = d2.copy()
      d2_copy.indices = [op.indices[i] for i in ti[1:]] + [op.indices[i] for i in bi[1:]]
      if indexOrder == '1122':
        d2_copy.indices = [d2_copy.indices[0], d2_copy.indices[2], d2_copy.indices[1], d2_copy.indices[3]]
      indexList = [op.indices[ti[0]], op.indices[bi[0]]]
      decomp.append( term( (-0.5)**n_perms, [], [d2_copy, sfExOp(indexList)] ) )

      # Create 1RDM 1RDM 1op terms
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[1]], op.indices[bi[1]]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op.indices[ti[2]], op.indices[bi[2]]]
      indexList = [op.indices[ti[0]], op.indices[bi[0]]]
      decomp.append( term( -2 * (-0.5)**n_perms, [], [d1_copy, d1_copy_2, sfExOp(indexList)] ) )
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[2]], op.indices[bi[1]]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op.indices[ti[1]], op.indices[bi[2]]]
      indexList = [op.indices[ti[0]], op.indices[bi[0]]]
      decomp.append( term( (-0.5)**n_perms, [], [d1_copy, d1_copy_2, sfExOp(indexList)] ) )

      # Create 1RDM 2RDM term
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[0]], op.indices[bi[0]]]
      d2_copy = d2.copy()
      d2_copy.indices = [op.indices[i] for i in ti[1:]] + [op.indices[i] for i in bi[1:]]
      if indexOrder == '1122':
        d2_copy.indices = [d2_copy.indices[0], d2_copy.indices[2], d2_copy.indices[1], d2_copy.indices[3]]
      decomp.append( term( - (-0.5)**n_perms, [], [d1_copy, d2_copy] ) )

      # Create 1RDM 1RDM 1RDM terms
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[0]], op.indices[bi[0]]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op.indices[ti[1]], op.indices[bi[1]]]
      d1_copy_3 = d1.copy()
      d1_copy_3.indices = [op.indices[ti[2]], op.indices[bi[2]]]
      decomp.append( term( (4.0 / 3.0) * (-0.5)**n_perms, [], [d1_copy, d1_copy_2, d1_copy_3] ) )
      d1_copy = d1.copy()
      d1_copy.indices = [op.indices[ti[0]], op.indices[bi[0]]]
      d1_copy_2 = d1.copy()
      d1_copy_2.indices = [op.indices[ti[2]], op.indices[bi[1]]]
      d1_copy_3 = d1.copy()
      d1_copy_3.indices = [op.indices[ti[1]], op.indices[bi[2]]]
      decomp.append( term( (-2.0 / 3.0) * (-0.5)**n_perms, [], [d1_copy, d1_copy_2, d1_copy_3] ) )

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def decomp_4op_to_2op_2rdm_sf(op, d1, d2):
  """
  Returns the list of terms composing the decomposition of the 4 body sfExOp op
  in terms of one and two body operators"
  """

  # Check input
  if ( not isinstance(op, sfExOp) ) or (op.order != 4):
    raise TypeError, "op must be a 4-particle sfExOp"
  if ( not isinstance(d1, tensor) ) or (len(d1.indices) != 2):
    raise TypeError, "d1 must be a tensor with 2 indices"
  if ( not isinstance(d2, tensor) ) or (len(d2.indices) != 4):
    raise TypeError, "d2 must be a tensor with 4 indices"

  # Initialize the return value
  decomp = []

  raise RuntimeError, "decomp_4op_to_2op_2rdm_sf not yet implemented"

  # Return the decomposition
  return decomp


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


