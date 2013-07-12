#    file:  sqaMisc2.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Defines some miscilaneous functions.
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
from sqaTensor import tensor, sfExOp, creOp, desOp
from sqaTerm import term, combineTerms, multiplyTerms, termChop
from sqaMisc import makeTuples, allDifferent, makePermutations, get_num_perms
from sqaOptions import options
import time

alpha_type = options.alpha_type
beta_type = options.beta_type


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def combine_transpose(termList):
  """
  Combines any terms in termList that are the transpose of each other.
  """

  # record the starting time
  startTime = time.time()

  # ensure all terms are in canonical form
  for t in termList:
    if not t.isInCanonicalForm:
      if options.verbose:
        print "making canonical...  %s" %(str(t))
      t.makeCanonical()

  # if requested, print a greeting
  if options.verbose:
    print ""
    print "Checking for transpose equivalencies and combining..."

  # initialize counter variable
  count = 0

  # loop over the terms
  j = 0
  while j < len(termList):

    # if requested, print each term
    if options.verbose:
      print '%i %s' %(count, str(termList[j]))

    # increment the counter
    count += 1

    # create a temporary term that is the transpose of the current term
    temp = termList[j].copy()
    creDes = []
    for i in range(len(temp.tensors)-1,-1,-1):
      ten = temp.tensors[i]
      if isinstance(ten, creOp):
        creDes.append(desOp(temp.tensors.pop(i).indices[0]))
      if isinstance(ten, desOp):
        creDes.append(creOp(temp.tensors.pop(i).indices[0]))
      if isinstance(ten, sfExOp):
        ind_list = temp.tensors.pop(i).indices
        order = len(ind_list) / 2
        creDes.append( sfExOp(ind_list[order:] + ind_list[0:order]) )
    temp.tensors.extend(creDes)
    del(creDes)
    temp.isInCanonicalForm = False
    temp.makeCanonical()

    # Search the remaining terms for a term matching the transpose
    # of the current term.  If a match is found, add the current term's
    # transpose to the matching term and delete the current term.
    for k in range(j+1,len(termList)):
      if temp.sameForm(termList[k]):
        termList[k].numConstant += temp.numConstant
        del(termList[j])
        break
    else:
      j += 1

  # if requested, print the elapsed time
  if options.verbose:
    print 'Transpose combination complete in %.3f seconds' %(time.time() - startTime)
    print ''


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convert_ops_to_rdms_so(inTerms, name, ord = 0):
  """
  Converts normal ordered creation and destruction operators in each input term in to the
  corresponding reduced density matrix.  The RDMs are given the supplied name.
  If ord is positive, only RDMs of that order are converted.
  """

  # Check that the name is a string
  if type(name) != type('a'):
    raise TypeError, 'name must be a string'

  # Process each term
  for t in inTerms:

    # Check that the term is normal ordered
    if not t.isNormalOrdered():
      raise ValueError, "Cannot convert operators in '%s' to density matrix, they are not normal ordered" %(str(t))

    # Determine number of creation and destruction operators
    creCount = 0
    desCount = 0
    for ten in t.tensors:
      if isinstance(ten, creOp):
        creCount += 1
      elif isinstance(ten, desOp):
        desCount += 1

    # Skip the term if it will not have an RDM of the specified order.
    # Note that ord <= 0 implies that RDMs of all orders should be treated.
    if ord > 0 and (creCount != desCount or ord != creCount):
      continue

    # Initialize list of RDM's indices
    indexList = []

    # Loop through the term's tensors and compile the RDM's index list
    i = 0
    while i < len(t.tensors):

      # Check that the term has no spin free excitation operators
      if isinstance(t.tensors[i], sfExOp):
        raise TypeError, "Cannot convert spin orbital operators in '%s' to density matrix, a sfExOp is present" %(str(t))

      # Process a creation operator
      elif isinstance(t.tensors[i], creOp):
        indexList.insert(creCount, t.tensors.pop(i).indices[0])

      # Process a destruction operator
      elif isinstance(t.tensors[i], desOp):
        indexList.insert(creCount, t.tensors.pop(i).indices[0])

      # Increment the index
      else:
        i += 1
    
    # Check that the number of cre/des ops are the same
    if creCount != desCount:
      raise ValueError, "Cannot convert spin orbital operators in '%s' to density matrix, unequal number of cre/des operators" %(str(t))

    # If any cre/des ops were converted, add the rdm to the term's tensor list
    if creCount > 0:

      # Construct the density matrix and add it to the term's tensors
      t.tensors.append(tensor(name, indexList))

      # Mark the term as not being in canonical form
      t.isInCanonicalForm = False


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def getABCounts(indices):
    ACount = 0
    BCount = 0
    for ind in indices:
      hasA = False
      hasB = False
      if alpha_type in ind.indType:
        hasA = True
        ACount += 1
      if beta_type in ind.indType:
        hasB = True
        BCount += 1
      if hasA and hasB:
        raise ValueError, "Index '%s' has both alpha and beta type." %(ind.name)
      if not hasA and not hasB:
        raise ValueError, "Index '%s' has neither alpha nor beta type." %(ind.name)
    return (ACount, BCount)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def assign_rdm_types(inTerms, rdm_name, rdms):
  """
  Assigns the rdm type (aa, bb, aaaa, bbbb, abab, etc.) to each tensor in the list
  of terms inTerms with name rdm_name.  rdms is a list of tensors representing the
  rdm types available for assignment.
  """

  # Prepare lists for the number of top and bottom alpha and beta indices in the supplied density matrices
  rdm_topACounts = [0]*len(rdms)
  rdm_topBCounts = [0]*len(rdms)
  rdm_bottomACounts = [0]*len(rdms)
  rdm_bottomBCounts = [0]*len(rdms)

  # Process the input rdms
  for i in range(len(rdms)):

    # If the number of indices is not even, raise an error
    if len(rdms[i].indices) % 2 != 0:
      raise ValueError, "Input rdms must have an even number of indices. rdm '%s' does not." %(str(rdms[i]))

    # Compute the rdm's order
    order = len(rdms[i].indices) / 2

    # Determine the numbers of alpha and beta indices in the top rdm indices
    (rdm_topACounts[i], rdm_topBCounts[i]) = getABCounts(rdms[i].indices[:order])

    # Determine the numbers of alpha and beta indices in the bottom rdm indices
    (rdm_bottomACounts[i], rdm_bottomBCounts[i]) = getABCounts(rdms[i].indices[order:])

#  for i in range(len(rdms)):
#    print "rdm: %20s   counts: %2i %2i %2i %2i" %(str(rdms[i]), rdm_topACounts[i], rdm_topBCounts[i], rdm_bottomACounts[i], rdm_bottomBCounts[i])

  # For each input term, loop through the tensors to assign rdm types
  for t in inTerms:
    for i in range(len(t.tensors)):

      # Assign a short name to the current tensor
      ten = t.tensors[i]

#      # If the tensor is a creOp or desOp, skip it
#      if isinstance(ten, creOp) or isinstance(ten, desOp):
#        continue

      # If the tensor does not have the specified name, skip it
      if ten.name != rdm_name:
        continue

      # If the number of indices is not even, raise an error
      if len(ten.indices) % 2 != 0:
        raise ValueError, "Density matrices must have an even number of indices. Tensor '%s' does not." %(str(ten))

      # Compute the tensor's order
      order = len(ten.indices) / 2

      # Determine the numbers of alpha and beta indices in the top tensor indices
      (topACount, topBCount) = getABCounts(ten.indices[:order])

      # Determine the numbers of alpha and beta indices in the bottom tensor indices
      (bottomACount, bottomBCount) = getABCounts(ten.indices[order:])

      # If the number of alpha and beta indices in the top and bottom do not match, the density matrix is zero
      if topACount != bottomACount or topBCount != bottomBCount:
        t.numConstant = 0.0
        break

      # Find the rdm with the matching number of top and bottom alpha and beta indices
      for j in range(len(rdms)):
        if rdm_topACounts[j] == topACount and \
           rdm_topBCounts[j] == topBCount and \
           rdm_bottomACounts[j] == bottomACount and \
           rdm_bottomBCounts[j] == bottomBCount:
          matching_rdm = rdms[j]
          break

      # If no rdm matches the top and bottom alpha and beta pattern, leave the tensor as is
      else:
        continue
#        raise RuntimeError, "No rdm with the correct number of alpha and beta top and bottom indices"

      # Create an index list for the assigned rdm
      indexList = ten.indices #[ind for ind in ten.indices]

      # Sort the top indices to match the rdm's alpha/beta order
      for j in range(order-1):
        if alpha_type in matching_rdm.indices[j].indType:
          target_type = alpha_type
        else:
          target_type = beta_type
        if target_type not in indexList[j].indType:
          k = order-1
          while k > j:
            if target_type in indexList[k].indType:
              temp = indexList[k]
              indexList[k] = indexList[j]
              indexList[j] = temp
              t.numConstant *= -1
              break
            k -= 1
          if k == j:
#            print "rdm = '%s'" %(str(matching_rdm))
#            print "ten = '%s'" %(str(ten))
#            print "tensor indices:"
#            for ind in indexList:
#              print ind.name, ", ", ind.indType
#            print "counts: %2i %2i %2i %2i" %(topACount, topBCount, bottomACount, bottomBCount)
            raise RuntimeError, "Could not sort the top indices of '%s' to match the rdm's alpha/beta ordering." %(str(ten))

      # Sort the bottom indices to match the rdm's alpha/beta order
      for j in range(order, len(indexList)-1):
        if alpha_type in matching_rdm.indices[j].indType:
          target_type = alpha_type
        else:
          target_type = beta_type
        if target_type not in indexList[j].indType:
          k = len(indexList)-1
          while k > j:
            if target_type in indexList[k].indType:
              temp = indexList[k]
              indexList[k] = indexList[j]
              indexList[j] = temp
              t.numConstant *= -1
              break
            k -= 1
          if k == j:
            raise RuntimeError, "Could not sort the bottom indices of '%s' to match the rdm's alpha/beta ordering." %(str(ten))

      # Replace the tensor with the assigned rdm
      t.tensors[i] = tensor(matching_rdm.name, indexList, matching_rdm.symmetries)

      # Mark the term as not being in canonical form
      t.isInCanonicalForm = False

  # Remove any zero terms
  termChop(inTerms)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


