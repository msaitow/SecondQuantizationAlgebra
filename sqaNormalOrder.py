#    file:  sqaNormalOrder.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Converts a term into normal order.
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


from sqaTensor import tensor, kroneckerDelta, creOp, desOp, sfExOp
from sqaTerm import term, sortOps
from sqaIndex import index
from sqaMisc import makeTuples, allDifferent, makePermutations, get_num_perms
from sqaOptions import options


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def normalOrder(inTerm):
  "Returns a list of terms resulting from normal ordering the operators in inTerm."

#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps and has_sfExOps:
    raise RuntimeError, "Normal ordering not implemented when both creOp/desOp and sfExOp " + \
                        "tensors are present"

  # if the term is already normal ordered, return it unchanged
  elif inTerm.isNormalOrdered():
    outTerms = [inTerm.copy()]

  # Normal ordering for creOp/desOp
  elif has_creDesOps:

    # Separate the cre/des operators from other tensors
    ops = []
    nonOps = []
    for t in inTerm.tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        ops.append(t.copy())
      else:
        nonOps.append(t.copy())

    # Generate all contraction pairs
    contractionPairs = []
    for i in range(len(ops)):
      iTerm = ops[i]
      for j in range(i+1,len(ops)):
        jTerm = ops[j]
        if isinstance(iTerm, desOp) and isinstance(jTerm, creOp):
          contractionPairs.append((i,j));
    #print "contractionPairs\n", contractionPairs

    # Determine maximum contraction order
    creCount = 0
    maxConOrder = 0
    for i in range(len(ops)-1,-1,-1):
      iTerm = ops[i]
      if isinstance(iTerm, creOp):
        creCount +=1
      elif isinstance(iTerm, desOp) and creCount > 0:
        maxConOrder += 1
        creCount -= 1
    del(creCount,iTerm)

    # Generate all contractions
    contractions = []
    for i in range(maxConOrder+1):
      subCons = makeTuples(i,contractionPairs)
      j = 0
      while j < len(subCons):
        creOpTags = []
        desOpTags = []
        for k in range(i):
          creOpTags.append(subCons[j][k][1])
          desOpTags.append(subCons[j][k][0])
        if allDifferent(creOpTags) and allDifferent(desOpTags):
          j += 1
        else:
          del(subCons[j])
      for j in range(len(subCons)):
        contractions.append(subCons[j])
    del(subCons,creOpTags,desOpTags,contractionPairs)
    #print "contractions:\n", contractions

    # For each contraction, generate the resulting term
    outTerms = []
    for contraction in contractions:
      conSign = 1
      deltaFuncs = []
      conIndeces = []
      subOpString = []
      subOpString.extend(ops)
      for conPair in contraction:
        index1 = ops[conPair[0]].indices[0]
        index2 = ops[conPair[1]].indices[0]
        deltaFuncs.append(kroneckerDelta([index1,index2]))
        subOpString[conPair[0]] = 'contracted'
        subOpString[conPair[1]] = 'contracted'
        for q in subOpString[conPair[0]+1:conPair[1]]:
          if not (q is 'contracted'):
            conSign *= -1
      i = 0
      while i < len(subOpString):
        if subOpString[i] is 'contracted':
          del(subOpString[i])
        else:
          i += 1
      (sortSign,sortedOps) = sortOps(subOpString)
      totalSign = conSign * sortSign
      outTensors = []
      outTensors.extend(nonOps)
      outTensors.extend(deltaFuncs)
      outTensors.extend(sortedOps)
      outTerms.append( term(totalSign * inTerm.numConstant, inTerm.constants, outTensors) )

  # Normal ordering for sfExOps
  elif has_sfExOps:

    # Make separate lists of the spin free excitation operators and other tensors
    sfExOp_list = []
    other_list = []
    for t in inTerm.tensors:
      if isinstance(t, sfExOp):
        sfExOp_list.append(t.copy())
      else:
        other_list.append(t.copy())

    # Initialize n, the number of remaining spin free excitation operators
    n = len(sfExOp_list)

    # Set the original term, with all excitation operators moved to the end, as the
    # first iteration's input term
    iter_input_terms = [term(inTerm.numConstant, inTerm.constants, other_list + sfExOp_list)]

    # Successively normal order the last two excitation operators until each term
    # has only one exitation operator left (at which point the term is normal ordered)
    while n > 1:

      # Initialize the list to hold this iteration's output terms
      iter_output_terms = []

      # For each of this iteration's input terms, produce all terms resulting from normal ordering
      # the last two excitation operators
      for t in iter_input_terms:

        # Make a list of the term's tensors that excludes the last two excitation operators
        tensors_except_last_two = []
        tensors_except_last_two.extend(t.tensors[0:-2])

        # Give short names for the last two excitation operators and their orders
        e1 = t.tensors[-2]
        e2 = t.tensors[-1]
        o1 = e1.order
        o2 = e2.order

        # Loop over the number of contractions
        for nc in range(min(o1,o2)+1):
          
          # Compute the order of excitation operator for the current number of contractions
          newOrder = o1 + o2 - nc

          # Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
          # the order in which the tuples may be combined to form a contraction
          perms = [0]
          if nc > 0:
            perms = makePermutations(nc)
          tups1 = makeTuples(nc, range(o1))
          tups2 = makeTuples(nc, range(o2))

          # For each contraction, compute the resulting term
          for perm in perms:
            for tup1 in tups1:
              for tup2 in tups2:

                # Initialize the term's tensor list
                tensorList = []
                tensorList.extend(tensors_except_last_two)

                # Compute the pairs of indices to be contracted
                # Example: (conPairs[0][p], conPairs[1][p]) is the pth contraction pair.
                conPairs = [[],[]]
                for i in range(nc):
                  conPairs[0].append(tup1[perm[i]])
                  conPairs[1].append(tup2[i])

                # Initialize the index list for the new excitation operator
                indexList = [False] * (2*newOrder)

                # Populate the index list for the new excitation operator
                # Also, create a kronecker delta function for each contraction pair
                for i in range(o1):
                  indexList[i] = e1.indices[i]                     # The usual one
                  if i in conPairs[0]:
                    i1 = i+o1
                    i2 = conPairs[1][conPairs[0].index(i)]
                    ind1 = e1.indices[i1]
                    ind2 = e2.indices[i2]
                    tensorList.insert(0, kroneckerDelta([ind1,ind2]))
                    indexList[i+newOrder] = e2.indices[o2+i2]
                  else:
                    indexList[i+newOrder] = e1.indices[i+o1]       # Put the usual counterpart into indexList
                count = 0                                          # in case of that e1.indices[i] is not contracted
                for i in range(o2):
                  if not (i in conPairs[1]):
                    indexList[o1+count] = e2.indices[i]
                    indexList[o1+count+newOrder] = e2.indices[i+o2]
                    count += 1

                # Ensure that all slots in the index list have been filled
                for ind in indexList:
                  if ind is False:
                    raise RuntimeError, "There is at least one unassigned index in the new spin free operator."

                # Add the new excitation operator to the tensor list
                tensorList.append(sfExOp(indexList))

                # Add the resulting term to this iteration's list of output terms
                iter_output_terms.append(term(inTerm.numConstant, inTerm.constants, tensorList))

      # Set this iteration's list of output terms as the next iteration's input terms
      iter_input_terms = iter_output_terms

      # Decrement the counter for the number of excitation operators
      n -= 1; #print "MADE OK!!!"
      
    # Set the return value as the final iteration's output terms
    outTerms = iter_output_terms

  else:
    raise RuntimeError, "Normal ordering function failed to choose what to do."

#  print "Terms after normal ordering:"
#  for t in outTerms:
#    print t

  return outTerms


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def contractCoreOps_sf(inTerm):
  "Returns a list of terms resulting from normal ordering the operators in inTerm."

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

  # Make separate lists of the spin free excitation operators and other tensors
  sfExOp_list = []
  other_list = []
  for t in inTerm.tensors:
    if isinstance(t, sfExOp):
      sfExOp_list.append(t.copy())
    else:
      other_list.append(t.copy())

  # Initialize n, the number of remaining spin free excitation operators
  if len(sfExOp_list) != 1:
    raise RuntimeError, "terms should have single sfExOp"

  # Give short names for the excitation operator and their order
  e1 = sfExOp_list[0]
  o1 = e1.order

  # list of offset of "core" type indices
  cIndsCre = []
  cIndsDes = []
  for i in range(o1):
    if e1.indices[i].indType == (options.core_type,):
      cIndsCre.append(i)
    if e1.indices[i+o1].indType == (options.core_type,):
      cIndsDes.append(i+o1)

  # all the "core" type indices should be contracted
  nc = len(cIndsCre)
  outTerms = []
  if nc != len(cIndsDes):
    outTerms.append(term(0.0, [], []))
    return outTerms
  if nc == 0:
    outTerms = [inTerm.copy()]
    return outTerms
  newOrder = o1 - nc

  # For each contraction, compute the resulting term
  perms = makePermutations(nc)
  for perm in perms:

    # Compute the pairs of indices to be contracted
    # Example: (conPairs[0][p], conPairs[1][p]) is the pth contraction pair.
    conPairs = [[],[]]
    for i in range(nc):
      conPairs[0].append(cIndsCre[perm[i]])
      conPairs[1].append(cIndsDes[i])

    # evaluate constant prefactor
    prefactor = 1
    inds = []
    for i in range(nc):
      if i not in inds:
        inds.append(i)
        i0 = conPairs[0][i]
        i1 = conPairs[1][i] - o1
        if i1 == i0:
          prefactor *= 2
          continue
#old        if i1 in conPairs[0]:
#old          inds.append(conPairs[0].index(i1))
#old          i2 = conPairs[1][conPairs[0].index(i1)] - o1
#old          if i2 == i0:
#old            prefactor *= 2
#old            continue
#old          if i2 in conPairs[0]:
#old            inds.append(conPairs[0].index(i2))
#old            i3 = conPairs[1][conPairs[0].index(i2)] - o1
#old            if i3 == i0:
#old              prefactor *= 2
#old              continue
#old            if i3 in conPairs[0]:
#old              inds.append(conPairs[0].index(i3))
#old              i4 = conPairs[1][conPairs[0].index(i3)] - o1
#old              if i4 == i0:
#old                prefactor *= 2
#old                continue
#old              else:
#old                raise RuntimeError, "NYI:recursive algorithm1"
        while i1 in conPairs[0]:
          inds.append(conPairs[0].index(i1))
          i1 = conPairs[1][conPairs[0].index(i1)] - o1
          if i1 == i0:
            prefactor *= 2
            break
#
        
    # Initialize the term's tensor list
    tensorList = other_list[:]

    # Initialize the index list for the new excitation operator
    indexList = [False] * (2*newOrder)

    # Populate the index list for the new excitation operator
    # Also, create a kronecker delta function for each contraction pair
    count = 0
    for i1 in range(o1):
      if i1 in conPairs[0]:#create a kronecker delta
        i2 = conPairs[1][conPairs[0].index(i1)]
        ind1 = e1.indices[i1]
        ind2 = e1.indices[i2]
        tensorList.insert(0, kroneckerDelta([ind1,ind2]))
      else:               #insert a pair of indices
        indexList[count] = e1.indices[i1]
        i2 = i1+o1
#old        if i2 in conPairs[1]:
#old          i2 = conPairs[0][conPairs[1].index(i2)] + o1
#old          if i2 in conPairs[1]:
#old            i2 = conPairs[0][conPairs[1].index(i2)] + o1
#old            if i2 in conPairs[1]:
#old              i2 = conPairs[0][conPairs[1].index(i2)] + o1
#old              if i2 in conPairs[1]:
#old                i2 = conPairs[0][conPairs[1].index(i2)] + o1
#old                if i2 in conPairs[1]:
#old                  raise RuntimeError, "NYI:recursive algorithm2"
#old              else:
#old                indexList[count+newOrder] = e1.indices[i2]
#old            else:
#old              indexList[count+newOrder] = e1.indices[i2]
#old          else:
#old            indexList[count+newOrder] = e1.indices[i2]
#old        else:
#old          indexList[count+newOrder] = e1.indices[i2]
        while i2 in conPairs[1]:
          i2 = conPairs[0][conPairs[1].index(i2)] + o1
        else: indexList[count+newOrder] = e1.indices[i2]
#
        count += 1

    # Ensure that all slots in the index list have been filled
    for ind in indexList:
      if ind is False:
        raise RuntimeError, "There is at least one unassigned index in the new spin free operator."

    # Add the new excitation operator to the tensor list
    if len(indexList) != 0:
      tensorList.append(sfExOp(indexList))

    # Determine the sign
    indPairs = [[],[]]
    for i in range(o1):
      indPairs[0].append(i)
      if i in conPairs[0]:
        indPairs[1].append(conPairs[1][conPairs[0].index(i)]-o1)
      else:
        a  = e1.indices[i]
        b = indexList[indexList.index(a)+newOrder]
        #indPairs[1].append(e1.indices.index(b)-o1)#is not correct when there are duplicate indices
        indPairs[1].append(e1.indices[o1:2*o1].index(b))

    #for i in range(o1): indPairs[1].index(i)
    n_perm = get_num_perms(indPairs[0], indPairs[1])
    sign = (-1)**n_perm

    # Add the resulting term to this iteration's list of output terms
    #print term(sign*prefactor*inTerm.numConstant, inTerm.constants, tensorList)
    #print "-------------------------------------------------------------------"
    outTerms.append(term(sign*prefactor*inTerm.numConstant, inTerm.constants, tensorList))

  return outTerms
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


