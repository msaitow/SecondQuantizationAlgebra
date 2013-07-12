#    file:  sqaTerm.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Defines the term class and some functions that use it.
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


import threading
from sqaIndex import index
from sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp, creDesTensor
from sqaMisc import makePermutations
from sqaOptions import options
import time

# The term class represents a set of constants and tensors that have been multiplied together.
# The class consists of a numerical constant factor, a list of named constants, and a list
# of tensors.
#
# It is common to encounter a list of term objects, which is often used to represent an expression.
# For example, the Hamiltonian in second quantization can be written as a list of terms.


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


class term:
  "A class for terms used in operator algebra.  Each term is a multiplicative string of constants, tensors, and operators."

  #------------------------------------------------------------------------------------------------

  def __init__(self, numConstant, constList, tensorList, isInCanonicalForm = False):
    self.constants = []
    self.tensors = []
    if type(numConstant) == type(1.0) or type(numConstant) == type(1):
      self.numConstant = float(numConstant)
    else:
      raise TypeError, "numConstant must be given as a float or an int."
    for c in constList:
      if type(c) != type('a'):
        raise TypeError, "constList must be a list of strings"
      self.constants.append(c)
    for t in tensorList:
      if not isinstance(t, tensor):
        raise TypeError, "tensorList must be a list of tensor objects"
      self.tensors.append(t.copy())
    if type(isInCanonicalForm) == type(True):
      self.isInCanonicalForm = isInCanonicalForm
    else:
      raise TypeError, "if specified, isInCanonicalForm must be True or False"

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):
    if not isinstance(other,term):
      raise TypeError, "term object can only be compared to other term objects."

    # sort by number of loose creation operators first
    retval = cmp(self.nCreOps(),other.nCreOps())
    if retval != 0:
      return retval

    # next sort by number of loose destruction operators
    retval = cmp(self.nDesOps(),other.nDesOps())
    if retval != 0:
      return retval

    # next sort by the orders of the spin free excitation operators
    retval = cmp(self.sfExOp_ranks(),other.sfExOp_ranks())
    if retval != 0:
      return retval

    # next sort by the number of constants
    retval = cmp(len(self.constants),len(other.constants))
    if retval != 0:
      return retval

    # next sort by the number of tensors
    retval = cmp(len(self.tensors),len(other.tensors))
    if retval != 0:
      return retval

    # next sort by the tensors' names
    retval = cmp([t.name for t in self.tensors], [t.name for t in other.tensors])
    if retval != 0:
      return retval

    # next sort by the tensors
    retval = cmp(self.tensors,other.tensors)
    if retval != 0:
      return retval

    # next sort by the constants
    retval = cmp(self.constants,other.constants)
    if retval != 0:
      return retval

    # finally compare the numerical constants
    numDiff = self.numConstant - other.numConstant
    if abs(numDiff) < 1e-6:
      return 0
    elif numDiff < 0:
      return -1
    elif numDiff > 0:
      return 1
    else:
      raise RuntimeError, "Failure in comparison of terms' numeric constants."
    return numDiff

  #------------------------------------------------------------------------------------------------

  def __str__(self):

    # numerical constant
    retval = " (%10.5f) " %self.numConstant

    # non-numerical constants
    for i in range(len(self.constants)):
      retval += self.constants[i] + " "

    # tensors
    for t in self.tensors:
      retval += str(t) + " "
    
    return retval

  #------------------------------------------------------------------------------------------------

  def __add__(self,other):
    "Adds the terms self and other.  Only works if the constants, tensors, and operators of two terms all match."
    if self.sameForm(other):
      retval = self.copy()
      retval.numConstant += other.numConstant
    else:
      raise RuntimeError, "self and other must have the same constants, tensors, and operators"
    return retval

  #------------------------------------------------------------------------------------------------

  def copy(self):
    "Returns a deep copy of the term."
    return term(self.numConstant, self.constants, self.tensors, self.isInCanonicalForm)

  #------------------------------------------------------------------------------------------------

  def nCreOps(self):
    "Returns the number of loose creation operators in the term"
    retval = 0
    for t in self.tensors:
      if isinstance(t, creOp):
        retval += 1
    return retval

  #------------------------------------------------------------------------------------------------

  def nDesOps(self):
    "Returns the number of loose destruction operators in the term"
    retval = 0
    for t in self.tensors:
      if isinstance(t, desOp):
        retval += 1
    return retval

  #------------------------------------------------------------------------------------------------

  def sfExOp_ranks(self):
    "Returns a list of the ranks of the spin free excitation operators in the term"
    retval = []
    for t in self.tensors:
      if isinstance(t, sfExOp):
        retval.append(len(t.indices)/2)
    return retval

  #------------------------------------------------------------------------------------------------

  def scale(self, factor):
    "Multiplies the term by factor."
    if type(factor) == type(1.0) or type(factor) == type(1):
      self.numConstant *= factor
    else:
      raise ValueError, "factor must an integer or float"

  #------------------------------------------------------------------------------------------------

  def sameForm(self,other):
    "Determines whether the terms are of the same form."
    self.makeCanonical()
    other.makeCanonical()
    if self.constants != other.constants or self.tensors != other.tensors:
      return False
    return True

  #------------------------------------------------------------------------------------------------

  def contractDeltaFuncs(self):
    "Contracts the indices within any kronecker delta functions in the term."

    i = 0
    while i < len(self.tensors):
      t = self.tensors[i]

#      print t
#      print "****************************"
#      if isinstance(t, kroneckerDelta):
#        i0 = t.indices[0]
#        i1 = t.indices[1]
#        print i0.indType
#        print i1.indType
      
      # If the term is a delta funciton with a repeated index, remove it
      if isinstance(t, kroneckerDelta) and (t.indices[0] == t.indices[1]):
#        print "2 ****************************"
        del(self.tensors[i])

      # If the term is a delta funciton with contractable indices, contract them and remove the delta func
      elif isinstance(t, kroneckerDelta) and (t.indices[0].isSummed or t.indices[1].isSummed):
#        print "3 ****************************"
        i0 = t.indices[0]
        i1 = t.indices[1]

        # Check that the indices have the same number of type groups
        if len(i0.indType) != len(i1.indType):
          raise RuntimeError, "Cannot contract indices %s, %s.  They have different numbers of type groups." %(str(i0), str(i1))

        # Determine the type of the new index based on the overlap between
        # the types of the two indices
        typeOverlap = []
        for j in range(len(i0.indType)):
          typeOverlap.append([])
          for typeString in i0.indType[j]:
            if typeString in i1.indType[j]:
              typeOverlap[-1].append(typeString)

#        for l0 in i0.indType:
#          typeOverlap.append([])
#          for s0 in l0:
#            for l1 in i1.indType:
#              if s0 in l1:
#                typeOverlap[-1].append(s0)

        # If there is no overlap between any of the type groups, the delta function is zero
        if ( len(i0.indType) > 0 or len(i1.indType) > 0 ) and [] in typeOverlap:
          self.numConstant = 0.0
          return

        # Create the new index
        if not i0.isSummed:
          newIndex = index(i0.name, typeOverlap, i0.isSummed)
        else:
          newIndex = index(i1.name, typeOverlap, i1.isSummed)

        # Remove the delta function
        del(self.tensors[i])

        # Excecute the index replacement
        for ten in self.tensors:
          for j in range(len(ten.indices)):
            if ten.indices[j] in [i0,i1]:
              ten.indices[j] = newIndex.copy()

      # Otherwise move on to the next tensor
      else:
        i += 1


#--------------------------------------------------------------------------------------------------

  def contractDeltaFuncs_new(self):
    "Contracts the indices within any kronecker delta functions in the term."

    i = 0
    while i < len(self.tensors):
      t = self.tensors[i]

      # If the term is a delta funciton with a repeated index, remove it
      if isinstance(t, kroneckerDelta) and (t.indices[0] == t.indices[1]):
        del(self.tensors[i])

      # If the term is a delta funciton with not-contractable indices
      elif isinstance(t, kroneckerDelta) and (not t.indices[0].isSummed) and (not t.indices[1].isSummed):
        i0 = t.indices[0]
        i1 = t.indices[1]
        # Check that the indices have the same number of type groups
        if len(i0.indType) != len(i1.indType):
          raise RuntimeError, "Cannot contract indices %s, %s.  They have different numbers of type groups." %(str(i0), str(i1))
        # Determine the type of the new index based on the overlap between
        # the types of the two indices
        typeOverlap = []
        for j in range(len(i0.indType)):
          typeOverlap.append([])
          for typeString in i0.indType[j]:
            if typeString in i1.indType[j]:
              typeOverlap[-1].append(typeString)
        # If there is no overlap between any of the type groups, the delta function is zero
        if ( len(i0.indType) > 0 or len(i1.indType) > 0 ) and [] in typeOverlap:
          self.numConstant = 0.0
          return
        else:
          i += 1

      # If the term is a delta funciton with contractable indices, contract them and remove the delta func
      elif isinstance(t, kroneckerDelta) and (t.indices[0].isSummed or t.indices[1].isSummed):
        i0 = t.indices[0]
        i1 = t.indices[1]

        # Check that the indices have the same number of type groups
        if len(i0.indType) != len(i1.indType):
          raise RuntimeError, "Cannot contract indices %s, %s.  They have different numbers of type groups." %(str(i0), str(i1))

        # Determine the type of the new index based on the overlap between
        # the types of the two indices
        typeOverlap = []
        for j in range(len(i0.indType)):
          typeOverlap.append([])
          for typeString in i0.indType[j]:
            if typeString in i1.indType[j]:
              typeOverlap[-1].append(typeString)

        # If there is no overlap between any of the type groups, the delta function is zero
        if ( len(i0.indType) > 0 or len(i1.indType) > 0 ) and [] in typeOverlap:
          self.numConstant = 0.0
          return

        # Create the new index
        if not i0.isSummed:
          newIndex = index(i0.name, typeOverlap, i0.isSummed)
        else:
          newIndex = index(i1.name, typeOverlap, i1.isSummed)

        # Remove the delta function
        del(self.tensors[i])

        # Excecute the index replacement
        for ten in self.tensors:
          for j in range(len(ten.indices)):
            if ten.indices[j] in [i0,i1]:
              ten.indices[j] = newIndex.copy()

      # Otherwise move on to the next tensor
      else:
        i += 1


  #------------------------------------------------------------------------------------------------

  def generateAlphabet(self):
    "Returns a list of strings that shares no elements with the term's index names."

    # Compile list of index names in self
    usedNames = []
    for t in self.tensors:
      for i in t.indices:
        if not (i.name in usedNames):
          usedNames.append(i.name)

    # Generate an alphabet with no elements overlapping with usedNames
    alphabet = []
    i = 0
    while len(alphabet) < len(usedNames):
      if not (str(i) in usedNames):
        alphabet.append(str(i))
      i += 1

    return alphabet

  #------------------------------------------------------------------------------------------------

  def isNormalOrdered(self):
    "Returns true if the term is in normal order and false otherwise"

    creFlag  = False 
    desFlag  = False
    sfExFlag = False
    for t in self.tensors:
      if ( isinstance(t,creOp) or isinstance(t,sfExOp) ) and ( desFlag or sfExFlag ):
        return False
      if isinstance(t, creOp):
        creFlag = True
      if isinstance(t, desOp):
        desFlag = True
      if isinstance(t, sfExOp):
        sfExFlag = True
    return True

  #------------------------------------------------------------------------------------------------

  def makeCanonical(self):
    "Converts the term to a unique canonical form."

    # Use the non recursive function
    self.makeCanonical_non_recursive()
    return

    # If the tensor is already in canonical form, do nothing
    if self.isInCanonicalForm:
      return

#    print "Converting to canonical form:"
#    print self

    # If the term is not normal ordered, raise an error
    if not self.isNormalOrdered():
      raise RuntimeError, "A term must have normal ordered operators to be converted to canonical form."

    # Sort the constants
    self.constants.sort()

    # If there are no tensors in the term then skip the tensor sorting.
    if len(self.tensors) == 0:
      self.isInCanonicalForm = True
      return

    # Create all tensor lists in which the tensors are sorted by type and name.
    # There can be more than one term here if there are multiple tensors with the same name.
    candidateTensorLists = self.getCandidateTensorLists()

    # Generate an alphabet for use in renaming indices
    # This is done to avoid renaming with an index name already in use.
    # Should try using numbers, i.e. '1', '2', '3', etc.
    alphabet = self.generateAlphabet()

    # For each candidate list, rename the indices canonically and record the score for that list.
    bestScore = [-1]
    nTopScore = 0
    for i in range(len(candidateTensorLists)):
      (score,map,sign,newTensorList) = getcim(candidateTensorLists[i],alphabet)
      if score > bestScore:
        (bestScore,bestMap,bestSign,bestTensorList) = (score,map,sign,newTensorList)
        nTopScore = 1
      elif score == bestScore:
        nTopScore += 1

    # Check to see that only one candidate achieved the top score
    if nTopScore > 1:
      raise RuntimeError, "%i candidates tied for the top score." %nTopScore

    # Set the tensors and their indices in the canonical order (the order with the highest score)
    self.tensors = bestTensorList

    # Apply the sign produced from the index sorting
    self.scale(bestSign)

    # Create an index mapping that converts to a canonical alphabet, i.e. a-z
    map = bestMap
    alphabet = list('abcdefghijklmnopqrstuvwxyz')
    if len(alphabet) < len(map):
      raise RuntimeError, "Alphabet smaller than number of indices, no more names left!"
    canonMap = {}
    while map.keys():
      minVal = min(map.values())
      for key in map.keys():
        if map[key] == minVal:
          canonMap[map[key].tup()] = map[key].copy()
          canonMap[map[key].tup()].name = alphabet.pop(0)
          del map[key]
          break
    map = canonMap

    # Rename the indices using the canonical mapping
    for t in self.tensors:
      for i in range(len(t.indices)):
        if t.indices[i].tup() in map.keys():
          t.indices[i] = map[t.indices[i].tup()].copy()

    # Turn on the canonical form flag to avoid calling this function again unnecessarily
    self.isInCanonicalForm = True

  #------------------------------------------------------------------------------------------------

  def makeCanonical_non_recursive(self):
    "Converts the term to a unique canonical form using a non-recursive algorithm."

    # If the tensor is already in canonical form, do nothing
    if self.isInCanonicalForm:
      return

    # If the term is not normal ordered, raise an error
    if not self.isNormalOrdered():
      raise RuntimeError, "A term must have normal ordered operators to be converted to canonical form."

    # Sort the constants
    self.constants.sort()

    # If there are no tensors in the term then skip the tensor sorting.
    if len(self.tensors) == 0:
      self.isInCanonicalForm = True
      return

    # Sort the freely commuting tensors by number and name
    fcList = []
    ncList = []
    for t in self.tensors:
      if t.freelyCommutes:
        fcList.append(t)
      else:
        ncList.append(t)
    fcList.sort(lambda x,y: cmp(x.name,y.name))
    nameGroups = []
    uniqueNames = []
    for t in fcList:
      if t.name in uniqueNames:
        nameGroups[-1].append(t)
      else:
        uniqueNames.append(t.name)
        nameGroups.append([t])
    nameGroups.sort(lambda x,y: cmp(len(x),len(y)))
    for t in ncList:
      nameGroups.append([t])
    del(uniqueNames,fcList,ncList,t)

    # Generate an alphabet for use in renaming indices
    # This is done to avoid renaming with an index name already in use.
    # Should try using numbers, i.e. '1', '2', '3', etc.
    alphabet = self.generateAlphabet()

    # Determine the best ordering and index mapping
    bestScore = [-1]
    nTopScore = 0
    # job format:  (map, gCount, tCount, aCount, gPerms)
    jobStack = [({},0,0,0,[])]
    while jobStack:

      # get the next job
      map,gCount,tCount,aCount,gPerms = jobStack.pop()

      # If there are no name groups remaining, compute the score
      if gCount == len(nameGroups):

        # Compute a new tensor list in which any dummy indices are given their
        # new names and all indices are sorted.
        # Also compute a list of these ordered indices.
        # Also compute the multiplicitive factor generated by sorting the indices.
        tenList = []
        indexList = []
        factor = 1
        for i in range(len(nameGroups)):
          for j in range(len(nameGroups[i])):
            tenList.append(nameGroups[i][gPerms[i][j]].copy())
            for k in range(len(tenList[-1].indices)):
              if tenList[-1].indices[k].isSummed:
                tenList[-1].indices[k] = map[tenList[-1].indices[k].tup()]
            factor *= tenList[-1].sortIndeces()
            for ind in tenList[-1].indices:
              indexList.append(ind)

        # Compute a score based on how alphabetical the indices are
        score = []
        for i in range(len(indexList)-1):
          score.append(0)
          for j in range(i+1,len(indexList)):
            if indexList[i] < indexList[j]:
              score[-1] += 1

        # If the current score is the best score, save the result
        if score > bestScore:
          nTopScore = 1
          bestScore = score
          bestMap = map
          best_factor = factor
          best_tensor_list = tenList

        # If the current score ties for the best, count the number of best scores
        elif score == bestScore:
          nTopScore += 1

      # If only cre/des operators remain, sort them and compute the score
      elif min([ (len(i) == 1 and (isinstance(i[0], creOp) or isinstance(i[0], desOp))) for i in nameGroups[gCount:] ]):

        # Compute a new tensor list in which any dummy indices are given their
        # new names and all indices are sorted.
        # Also compute a list of these ordered indices.
        # Also compute the multiplicitive factor generated by sorting the indices.
        tenList = []
        indexList = []
        factor = 1
        for i in range(gCount):
          for j in range(len(nameGroups[i])):
            tenList.append(nameGroups[i][gPerms[i][j]].copy())
            for k in range(len(tenList[-1].indices)):
              if tenList[-1].indices[k].isSummed:
                tenList[-1].indices[k] = map[tenList[-1].indices[k].tup()]
            factor *= tenList[-1].sortIndeces()
            for ind in tenList[-1].indices:
              indexList.append(ind)

        # Apply the input mapping to the creation/destruction operators
        # Create and apply a new mapping for any new dummy indices
        opList = [nameGroups[i][0].copy() for i in range(gCount,len(nameGroups))]
        nNewMaps = 0
        for op in opList:
          if op.indices[0].tup() in map:
            op.indices[0] = map[op.indices[0].tup()]
          elif op.indices[0].isSummed:
            map[op.indices[0].tup()] = index(alphabet[aCount+nNewMaps], op.indices[0].indType, op.indices[0].isSummed)
            nNewMaps += 1
            op.indices[0] = map[op.indices[0].tup()]

        # Sort the operators and apply the resulting sign
        (s,opList) = sortOps(opList)
        factor *= s

        # Add the operators' indices to the ordered list of indices.
        # Also add the sorted operators to the new tensor list.
        for op in opList:
          indexList.append(op.indices[0])
          tenList.append(op)

        # Compute a score based on how alphabetical the indices are
        score = []
        for i in range(len(indexList)-1):
          score.append(0)
          for j in range(i+1,len(indexList)):
            if indexList[i] < indexList[j]:
              score[-1] += 1

        # If the current score is the best score, save the result
        if score > bestScore:
          nTopScore = 1
          bestScore = score
          bestMap = map
          best_factor = factor
          best_tensor_list = tenList

        # If the current score ties for the best, count the number of best scores
        elif score == bestScore:
          nTopScore += 1

      # If only a sfExOp remains, sort its indices and compute the score
      elif (gCount == len(nameGroups)-1) and (len(nameGroups[gCount]) == 1) and isinstance(nameGroups[gCount][0], sfExOp):

        # Compute a new tensor list in which any dummy indices are given their
        # new names and all indices are sorted.
        # Also compute a list of these ordered indices.
        # Also compute the multiplicitive factor generated by sorting the indices.
        tenList = []
        indexList = []
        factor = 1
        for i in range(gCount):
          for j in range(len(nameGroups[i])):
            tenList.append(nameGroups[i][gPerms[i][j]].copy())
            for k in range(len(tenList[-1].indices)):
              if tenList[-1].indices[k].isSummed:
                tenList[-1].indices[k] = map[tenList[-1].indices[k].tup()]
            factor *= tenList[-1].sortIndeces()
            for ind in tenList[-1].indices:
              indexList.append(ind)

        # Apply the input mapping to the sfExOp
        # Create and apply a new mapping for any new dummy indices
        t = nameGroups[gCount][0].copy()
        nNewMaps = 0
        for i in range(len(t.indices)):
          if t.indices[i].tup() in map:
            t.indices[i] = map[t.indices[i].tup()]
          elif t.indices[i].isSummed:
            map[t.indices[i].tup()] = index(alphabet[aCount+nNewMaps], t.indices[i].indType, t.indices[i].isSummed)
            nNewMaps += 1
            t.indices[i] = map[t.indices[i].tup()]

        # Sort the indices of the sfExOp (go go gadget bubble sort!)
        i = 0
        while i < t.order-1:
          if t.indices[i] > t.indices[i+1]:
            temp = t.indices[i+1]
            t.indices[i+1] = t.indices[i]
            t.indices[i] = temp
            temp = t.indices[i+t.order+1]
            t.indices[i+t.order+1] = t.indices[i+t.order]
            t.indices[i+t.order] = temp
            i = 0
          else:
            i += 1

        # Add the sfExOp to the new tensor list
        tenList.append(t)

        # Add the sfExOp's indices to the ordered index list
        for i in t.indices:
          indexList.append(i)

        # Compute a score based on how alphabetical the indices are
        score = []
        for i in range(len(indexList)-1):
          score.append(0)
          for j in range(i+1,len(indexList)):
            if indexList[i] < indexList[j]:
              score[-1] += 1

        # If the current score is the best score, save the result
        if score > bestScore:
          nTopScore = 1
          bestScore = score
          bestMap = map
          best_factor = factor
          best_tensor_list = tenList

        # If the current score ties for the best, count the number of best scores
        elif score == bestScore:
          nTopScore += 1

      # If starting a new name group, sort the group's names based on the input mapping
      # and schedule a new job for each permutation of the tensors with no index assignments
      elif len(gPerms) <= gCount:
        withMapped = []
        withoutMapped = []
        for i in range(len(nameGroups[gCount])):
          leastMapped = False
          for ind in nameGroups[gCount][i].indices:
            if (ind.tup() in map) and ((leastMapped is False) or (map[ind.tup()] < leastMapped)):
              leastMapped = map[ind.tup()]
          if leastMapped is False:
            withoutMapped.append(i)
          else:
            withMapped.append((leastMapped,i))
        withMapped.sort(lambda x,y: cmp(x[0],y[0]))
        withMapped = [i[1] for i in withMapped]
        if len(withoutMapped) <= 1:
          jobStack.append((map,gCount,tCount,aCount,gPerms + [withMapped + withoutMapped]))
        else:
          for perm in makePermutations(len(withoutMapped)):
            new_gPerm = withMapped + [withoutMapped[i] for i in perm]
            jobStack.append((map,gCount,tCount,aCount,gPerms + [new_gPerm]))
          del(new_gPerm,perm)
        del(withMapped,withoutMapped,leastMapped)

      # If continuing an existing group, schedule a new job for each equivelent ordering
      # of the current tensor's indices
      else:

        # Give the current tensor a convenient name
        t = nameGroups[gCount][gPerms[gCount][tCount]]

        # Get the tensor's symmetry permutations
        (symPerms,factors) = t.symPermutes()

        # Compute gCount and tCount for the next job
        next_gCount = gCount
        next_tCount = tCount + 1
        if next_tCount == len(nameGroups[gCount]):
          next_gCount += 1
          next_tCount = 0

        # For each of the tensor's symmetry-equivelent index orderings, create mappings for any
        # un-mapped dummy indices and schedule a new job
        for perm in symPerms:
          nNewMaps = 0
          newMap = {}
          newMap.update(map)
          for ind in [t.indices[perm[i]] for i in range(len(t.indices))]:
            if ind.isSummed and ind.tup() not in newMap:
              newMap[ind.tup()] = index(alphabet[aCount+nNewMaps], ind.indType, ind.isSummed)
              nNewMaps += 1
          jobStack.append((newMap,next_gCount,next_tCount,aCount+nNewMaps,gPerms))

#    # Check to see that only one candidate achieved the top score
#    if nTopScore > 1:
#      #raise RuntimeError, "%i candidates tied for the top score." %nTopScore
#      print "WARNING: %i candidates tied for the top score." %nTopScore

    # Set the tensor list as the list with the 'best' index naming and ordering
    self.tensors = best_tensor_list 

    # Apply the factor produced from the canonical ordering
    self.scale(best_factor)

    # Create an index mapping that converts to a canonical alphabet, i.e. a-z
    alphabet = list('abcdefghijklmnopqrstuvwxyz')
    if len(alphabet) < len(bestMap):
      raise RuntimeError, "Alphabet smaller than number of indices, no more names left!"
    canonMap = {}
    while bestMap.keys():
      minVal = min(bestMap.values())
      for key in bestMap.keys():
        if bestMap[key] == minVal:
          canonMap[bestMap[key].tup()] = bestMap[key].copy()
          canonMap[bestMap[key].tup()].name = alphabet.pop(0)
          del bestMap[key]
          break

    # Rename the indices using the canonical mapping
    for t in self.tensors:
      for i in range(len(t.indices)):
        if t.indices[i].tup() in canonMap.keys():
          t.indices[i] = canonMap[t.indices[i].tup()].copy()

    # Turn on the canonical form flag to avoid calling this function again unnecessarily
    self.isInCanonicalForm = True

  #------------------------------------------------------------------------------------------------

  def getCandidateTensorLists(self):
    """
    Create all tensor lists in which the term's tensors are sorted by name.
    Multiple lists occur when there are multiple freely commuting tensors with the same name.
    """

    # Separate the tensors into two lists:  freely commuting and non-commuting
    fcList = []
    ncList = []
    for t in self.tensors:
      if t.freelyCommutes:
        fcList.append(t.copy())
      else:
        ncList.append(t.copy())

    if len(fcList) > 0:

      # Sort the freely commuting tensors by name
      fcList.sort(lambda x,y: cmp(x.name,y.name))

      # For the freely commuting tensors, compile a list of unique tensor names and the number of times they occur
      uniqueNames = []
      nameCounts = []
      for ten in fcList:
        if ten.name in uniqueNames:
          nameCounts[uniqueNames.index(ten.name)] += 1
        else:
          uniqueNames.append(ten.name)
          nameCounts.append(1)

      # For each unique name, generate ordering permutations
      permutes = []
      for i in nameCounts:
        permutes.append(makePermutations(i))

      # Combine the name perturbations above into all possible lists in which the freely commuting tensors
      # are ordered by name
      tensorsByName = []
      total = 0
      for i in nameCounts:
        tensorsByName.append(fcList[total:total+i])
        total += i
      candidateLists = []
      for perm in permutes[0]:
        candidateLists.append([])
        for p in perm:
          candidateLists[-1].append(tensorsByName[0][p].copy())
      for i in range(1,len(uniqueNames)):
        nOldVariants = len(candidateLists)
        for perm in permutes[i]:
          for j in range(nOldVariants):
            candidateLists.append([])
            candidateLists[-1].extend(candidateLists[j])
            for p in perm:
              candidateLists[-1].append(tensorsByName[i][p].copy())
        del(candidateLists[0:nOldVariants])

    else:
      candidateLists = [[]]

    # Append the non-commuting tensors to each candidate list
    for l in candidateLists:
      for t in ncList:
        l.append(t.copy())

    # Return the list of candidate tensor orders
    return candidateLists

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def combineTerms(termList, maxThreads = 1):
  "Combines any like terms in termList"

  if options.verbose:
    print ''
    print 'Combining like terms:'
    print 'Converting %i terms to canonical form...' %(len(termList))

  startTime = time.time()

  # Put the terms in termList into their canonical (unique) forms
  if maxThreads > 1:

    # Initialize counters and locks
    termCount = [0]
    printCount = [0]
    tLock = threading.Lock()
    pLock = threading.Lock()

    # Define function to use in threads
    def threadFunc(nTerms):

      batchSize = 100

      # Get first batch
      tLock.acquire()
      i = termCount[0]
      termCount[0] += batchSize
      tLock.release()
      k = i + batchSize

      while i < nTerms:

        termList[i].makeCanonical()

#        pLock.acquire()
#        print '%6i  %s' %(printCount[0],str(termList[i]))
#        printCount[0] += 1
#        pLock.release()

        i += 1

        if i == k:
          # Get next batch
          tLock.acquire()
          i = termCount[0]
          termCount[0] += batchSize
          tLock.release()
          k = i + batchSize

    # Start threads
    threads = []
    nTerms = len(termList)
    for i in range(maxThreads):
      threads.append(threading.Thread(target=threadFunc, args=(nTerms,)))
      threads[-1].start()

    # Wait for the threads to finish
    for thread in threads:
      thread.join()

  else:
    # Convert the terms to canonical form in the main thread
    for i in range(len(termList)):
      if options.verbose:
        print '%6i  %s' %(i,str(termList[i]))
      termList[i].makeCanonical()

  # Sort the terms
  termList.sort()

  # Combine any terms with the same canonical form
  i = 0
  while i < len(termList)-1:
    if termList[i].sameForm(termList[i+1]):
      termList[i+1] += termList[i]
      del(termList[i])
    else:
      i += 1
#  i = 0
#  while i < len(termList)-1:
#    j = i + 1
#    while j < len(termList):
#      if termList[j].sameForm(termList[i]):
#        termList[i] += termList[j]
#        del(termList[j])
#      else:
#        j += 1
#    i += 1

  # Remove terms with coefficients of zero
  termChop(termList)

  if options.verbose:
    print "Finished combining terms in %.3f seconds" %(time.time() - startTime)
    print ""

#  # Sort the terms
#  termList.sort()


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def multiplyTerms(t1,t2):
  if (not isinstance(t1,term)) or (not isinstance(t2,term)):
    raise TypeError, "t1 and t2 must be of type term"
  numConst = t1.numConstant * t2.numConstant
  constList = []
  constList.extend(t1.constants)
  constList.extend(t2.constants)
  tensorList = []
  tensorList.extend(t1.tensors)
  tensorList.extend(t2.tensors)
  return term(numConst,constList,tensorList)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def termChop(termList, tolerance = 1e-6):
  "Removes any terms with zero constant factors from termList."
  TypeErrorMessage = "termList must be a list of terms"
  if type(termList) != type([]):
    raise TypeError, TypeErrorMessage
  i = 0
  while i < len(termList):
    if not isinstance(termList[i],term):
      raise TypeError, TypeErrorMessage
    if abs(termList[i].numConstant) < tolerance:
      del(termList[i])
    else:
      i += 1


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def getcim(tenList, alphabet, tenCount = 0, alphaCount = 0, inputMaps = {}):
  """
  Determines the index mapping necessary to canonicalize the indices in tenList.
  Assumes any creation/destruction operators are in normal order.
  """

  # Determine what types of tensors are left to process
  no_tensors_left  = ( tenCount == len(tenList) )
  only_sfExOp_left = ( (tenCount == len(tenList)-1) and isinstance(tenList[-1], sfExOp) )
  only_creDes_left = ( tenCount < len(tenList) )
  for t in tenList[tenCount:]:
    if ( not isinstance(t, creOp) ) and ( not isinstance(t, desOp) ):
      only_creDes_left = False
      break

  if no_tensors_left or only_creDes_left or only_sfExOp_left:

    # Copy the input mapping into a new map that can be added to
    map = {}
    map.update(inputMaps)

    # For the preceeding tensors, create an ordered list of indeces
    # and a list of tensors with the canonical index sorting
    indexList = []
    newTensorList = []
    sign = 1
    for t in tenList[:tenCount]:
      tcopy = t.copy()
      for j in range(len(tcopy.indices)):
        if tcopy.indices[j].tup() in map.keys():
          tcopy.indices[j] = map[tcopy.indices[j].tup()]
      # Keep track of the sign produced by sorting the tensor's indices
      sign *= tcopy.sortIndeces()
      newTensorList.append(tcopy)
      for ind in tcopy.indices:
        indexList.append(ind)

  if only_creDes_left:

    # Apply the input mapping to the creation/destruction operators
    # Create and apply a new mapping for any new dummy indices
    opList = []
    for op in tenList[tenCount:]:
      opList.append(op.copy())
    nNewMaps = 0
    for op in opList:
      if op.indices[0].tup() in map.keys():
        # Apply input mapping
        op.indices[0] = map[op.indices[0].tup()].copy()
      elif op.indices[0].isSummed:
        # Create new mapping
        map[op.indices[0].tup()] = index(alphabet[alphaCount+nNewMaps], op.indices[0].indType, op.indices[0].isSummed)
        nNewMaps += 1
        op.indices[0] = map[op.indices[0].tup()].copy()

    # Sort the operators and apply the resulting sign
    (s,opList) = sortOps(opList)
    sign *= s

    # Add the operators' indices to the ordered list of indices.
    # Also add the sorted operators to the new tensor list.
    for op in opList:
      indexList.append(op.indices[0])
      newTensorList.append(op)

  if only_sfExOp_left:

    # Apply the input mapping to the sfExOp
    # Create and apply a new mapping for any new dummy indices
    t = tenList[-1].copy()
    nNewMaps = 0
    for i in range(len(t.indices)):
      if t.indices[i].tup() in map.keys():
        # Apply input mapping
        t.indices[i] = map[t.indices[i].tup()].copy()
      elif t.indices[i].isSummed:
        # Create new mapping
        map[t.indices[i].tup()] = index(alphabet[alphaCount+nNewMaps], t.indices[i].indType, t.indices[i].isSummed)
        nNewMaps += 1
        t.indices[i] = map[t.indices[i].tup()].copy()

    # Sort the indices of the sfExOp (go go gadget bubble sort!)
    i = 0
    while i < t.order-1:
      if t.indices[i] > t.indices[i+1]:
        temp = t.indices[i+1]
        t.indices[i+1] = t.indices[i]
        t.indices[i] = temp
        temp = t.indices[i+t.order+1]
        t.indices[i+t.order+1] = t.indices[i+t.order]
        t.indices[i+t.order] = temp
        i = 0
      else:
        i += 1

    # Add the sfExOp to the new tensor list
    newTensorList.append(t)

    # Add the sfExOp's indices to the ordered index list
    for i in t.indices:
      indexList.append(i)
        

  if no_tensors_left or only_creDes_left or only_sfExOp_left:

    # Compute a score based on how alphabetical the ordered list of indices is
    score = []
    for i in range(len(indexList)-1):
      score.append(0)
      for j in range(i+1,len(indexList)):
        if indexList[i] < indexList[j]:
          score[-1] += 1

    # Return the score, the mapping, the resulting sign, and the canonical tensor list produced by the mapping
    return (score,map,sign,newTensorList)

  # Otherwise, process the next tensor
  else:

    # Get the tensor's symmetry permutations
    (symPerms,factors) = tenList[tenCount].symPermutes()

    # Make a copy of the tensor to be processed
    t = tenList[tenCount].copy()

    # Apply all index maps from previous tensors to the current tensor
    for i in range(len(t.indices)):
      if t.indices[i].tup() in inputMaps:
        t.indices[i] = inputMaps[t.indices[i].tup()]

    # Determine the permutation that maximizes alphabetical order of the mapped indices
    bestPermScore = -1
    for perm in symPerms:
      permScore = 0
      indList = []
      for j in range(len(t.indices)):
        indList.append(t.indices[perm[j]])
      for i in range(len(indList)-1):
        for j in range(i+1,len(indList)):
          if indList[j] in inputMaps.values() and indList[i] in inputMaps.values() and indList[i] < indList[j]:
            permScore += 1
      if permScore > bestPermScore:
        bestPermScore = permScore
        bestInitialPerm = perm

    # Does it matter if there are more than one perm with max score?
    # I think it doesn't, one can select any of them

    # Reset the indices and sort them according to bestInitialPerm
    t = tenList[tenCount].copy()
    tcopy = t.copy()
    for i in range(len(t.indices)):
      tcopy.indices[bestInitialPerm[i]] = t.indices[i]
    t = tcopy

    # make a mapping of the next alphabet elements, in order, to the unassigned indices
    # also make any symmetry equivalent mappings
    bestScore = [-1]
    for k in range(len(symPerms)):
      nNewMaps = 0
      map = {}
      map.update(inputMaps)
      for i in range(len(t.indices)):
        p = symPerms[k][i]
        if t.indices[p].isSummed and ( not (t.indices[p].tup() in map) ):
          map[t.indices[p].tup()] = index(alphabet[alphaCount+nNewMaps], t.indices[p].indType, t.indices[p].isSummed)
          nNewMaps += 1
      (score,map,sign,newTensorList) = getcim(tenList, alphabet, tenCount+1, alphaCount+nNewMaps, map)
      if score > bestScore: #is it possible to have multiple top scores here? I think so.  Does it matter?
        bestScore = score
        bestMaps = map
        bestSign = sign
        bestNewTensorList = newTensorList

    # return the mapping that gives the best overall score
    return (bestScore,bestMaps,bestSign,bestNewTensorList)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def sortOps(unsortedOps, returnPermutation = False):
  """
  Sorts a list of creation/destruction operators into normal order and alphebetically.
  Performs no contractions.  Returns the overall sign resulting from the sort and the sorted operator list.
  """
  sortedOps = unsortedOps + []
  i = 0
  sign = 1
  if returnPermutation:
    perm = range(len(unsortedOps))
  while i < len(sortedOps)-1:
    if sortedOps[i] <= sortedOps[i+1]:
       i += 1
    else:
      temp = sortedOps[i]
      sortedOps[i] = sortedOps[i+1]
      sortedOps[i+1] = temp
      if returnPermutation:
        temp = perm[i]
        perm[i] = perm[i+1]
        perm[i+1] = temp
      i = 0
      sign *= -1
  if returnPermutation:
    return (sign,sortedOps,perm)
  return (sign,sortedOps)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeCoreOpPairs(inList):
  """
  Removes pairs of core creation and core destruction operators corresponding to the same core index.
  Does not remove a pair if it's creation or destruction operator is repeated.
  Input is a list of terms.
  The terms must be in normal order.
  """

  # prepare input argument
  if type(inList) != type([]):
    raise TypeError, "input must be a list of terms"

  # loop over input terms
  for t in inList:

    # Check that the term is indeed a term
    if not isinstance(t, term):
      raise TypeError, "input must be a term or list of terms"

    # Check that the term is normal ordered
    if not t.isNormalOrdered():
      raise ValueError, "core index removal function only works for normal ordered terms"

    # Initialize a counter for unremoved creation operators
    creCount = 0

    # Loop over the term's tensors
    i = 0
    while i < len(t.tensors):

      # Initialize flags
      operatorsRemoved = False
      repeatedCreOp = False
      repeatedDesOp = False

      # if the tensor is a core creation operator
      if isinstance(t.tensors[i], creOp) and options.core_type in t.tensors[i].indices[0].indType:

        # Check whether the creation operator is a repeat of an earlier creation operator
        for k in range(i):
          if t.tensors[k] == t.tensors[i]:
            repeatedCreOp = True

        # If a repeat, move to the next tensor
        if not repeatedCreOp:

          # otherwise...
          
          # create the matching destruction operator
          matchingDesOp = desOp(t.tensors[i].indices[0])

          # search for the matching destruction operator
          for j in range(i+1,len(t.tensors)):

            # if a repeat of the creation operator is found, move to the next tensor
            if t.tensors[j] == t.tensors[i]:
              break

            # if the matching destruction operator is found
            if t.tensors[j] == matchingDesOp:

              # initialize a counter for the number of operator commutations necessary to move the
              # matching creation and destruction operators to the begining and end of the term,
              # respectively
              commCount = creCount

              # count the number of destruction operators after the matching destruction operator
              for k in range(j+1, len(t.tensors)):
                if isinstance(t.tensors[k], desOp):
                  commCount += 1

                # while counting, check whether a repeat of the matching destruction operator is present
                if t.tensors[k] == t.tensors[j]:
                  repeatedDesOp = True

              # if there is a repeat of the matching destruction operator, move to the next tensor
              if repeatedDesOp:
                break

              # scale the term by the factor resulting from commuting the creation operator and
              # matching destruction operator to the beginning and end of the term, respectively
              t.scale((-1)**commCount)

              # delete the creation and matching destruction operators
              del t.tensors[j]
              del t.tensors[i]

              # set the operator removal flag to true
              operatorsRemoved = True

              # move to the next tensor
              break

      # If no operators were removed...
      if not operatorsRemoved:

        # If the current tensor is a creation operator, increase the creation operator count
        if isinstance(t.tensors[i], creOp):
          creCount += 1

        # increment the index
        i += 1


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeCoreOps_sf(inList):
  """
  Removes core indices from inList's terms' spin-free excitation operators.
  This function assumes that the spin-free operators will be converted to density matrices
  by taking their expectation value immediately after this function is finished.
  Terms that have a zero expectation value due to the nature of their spin-free operator's core
  indices are deleted from inList.
  """

  if options.verbose:
    print "removing core creation and destruction operators in preperation for conversion to RDMs by an expectation value..."
    print ""

  # loop repeatedly through the terms until no core indices are left
  hasCore = True
  while hasCore:

    hasCore = False

    # process each term
    t_num = 0
    while t_num < len(inList):

      # use a short name for the current term
      t = inList[t_num]

      # check that t is a term
      if not isinstance(t, term):
        raise TypeError, "inList must be a list of term objects"

      # check for normal ordering
      if not t.isNormalOrdered():
        raise ValueError, "input terms must be normal ordered"

      # check for spin-orbital creation and destruction operators
      for ten in t.tensors:
        if isinstance(ten, creOp) or isinstance(ten, desOp):
          raise TypeError, "input terms may not contain creOp or desOp objects"

      # find the spin-free excitation operator
      op = None
      for i in xrange(len(t.tensors)):
        if isinstance(t.tensors[i], sfExOp):
          op = t.tensors[i]
          opPos = i

      # if there is no spin-free excitation operator, skip the term
      if op is None:
        t_num += 1
        continue

      # ensure the sfExOp has no indices with multiple type groups or a type group with core and non-core types
      for ind in op.indices:
        if len(ind.indType) > 1:
          raise ValueError, "index %s in term (%s) has more than one type group:  %s" %(ind.name, str(t), str(typeGroup))
        for typeGroup in ind.indType:
          if options.core_type[0] in typeGroup and len(typeGroup) > 1:
            raise ValueError, "index %s in term (%s) has a type group including core and non-core types:  %s" %(ind.name, str(t), str(typeGroup))

      # compute the operator's order
      order = len(op.indices)/2

      # find a core index
      cInd = None
      for i in xrange(2*order):
        if op.indices[i].indType == (options.core_type,):
          cInd = op.indices[i]
          break

      # if there are no core indices, move to the next term
      if cInd is None:
        t_num += 1
        continue

      # if there is a core index, request another loop through the terms because all core indices have not been found
      hasCore = True

      # count the number of times the targeted core index appears among creation and destruction operators
      nCre = 0
      nDes = 0
      for i in xrange(order):
        if op.indices[i] == cInd:
          nCre += 1
        if op.indices[order+i] == cInd:
          nDes += 1

      # if the term is equal to zero, remove it and move to the next term
      if nCre != nDes or nCre > 2 or nDes > 2:
        del inList[t_num]
        continue

      # organize the operator's indices into vertical pairs of cre/des operator indices
      pairs = [ [op.indices[i],op.indices[order+i]] for i in xrange(order)]

      # print out the initial term
      if options.verbose:
        print "  initial term: ", t
#        print "verticle pairs: ",
#        for p in pairs:
#          print " [%s,%s]" %(p[0].name, p[1].name),
#        print ""

      # make sure the number of pairs is equal to the operator's order
      if len(pairs) != order:
        raise ValueError, "number of pairs not equal to operator's order"

      # determine the new operator's indices.
      # record how many pairs there were with both elements equal to the targeted core operator
      nMatch = 0
      topUnmatched = []
      botUnmatched = []
      i = order-1
      while i >= 0:
        if pairs[i][0] == cInd and pairs[i][1] == cInd:
          del pairs[i]
          nMatch += 1
        elif pairs[i][0] == cInd:
          botUnmatched.append(pairs.pop(i)[1])
        elif pairs[i][1] == cInd:
          topUnmatched.append(pairs.pop(i)[0])
        i -= 1
      newIndices = []
      newIndices.extend(topUnmatched)
      for p in pairs:
        newIndices.append(p[0])
      newIndices.extend(botUnmatched)
      for p in pairs:
        newIndices.append(p[1])

      # replace the old operator with the new operator in which the targeted core index has been removed
      if len(newIndices) > 0:
        t.tensors[opPos] = sfExOp(newIndices)
      # if there are no indices left after the core index's removal, remove the old operator
      else:
        del t.tensors[opPos]

      # apply the appropriate constant factor
      if   nCre == 1 and nMatch == 1:
        t.scale(2.0)
      elif nCre == 1 and nMatch == 0:
        t.scale(-1.0)
      elif nCre == 2 and nMatch == 2:
        t.scale(2.0)
      elif nCre == 2 and nMatch == 1:
        t.scale(-1.0)
      elif nCre == 2 and nMatch == 0:
        t.scale(1.0)
      else:
        raise ValueError, "unexpected values:  nCre = %i, nMatch = %i" %(nCre, nMatch)

      # print out the final term
      if options.verbose:
#        print "  topUnmatched: ",
#        for p in topUnmatched:
#          print " %s" %(p.name),
#        print ""
#        print "  botUnmatched: ",
#        for p in botUnmatched:
#          print " %s" %(p.name),
#        print ""
        print "    final term: ", t

      # for the special case of two unmatched pairs, the result is a sum of two different operators.
      # the first replaced the original operator, and the second is added here.
      if nCre == 2 and nMatch == 0:
        inList.append(t.copy())
        if len(newIndices) < 4:
          raise ValueError, "expected at least 4 remaining indices for nCre == 2 and nMatch == 0 case, but only %i are present" %len(newInices)
        (newIndices[0], newIndices[1]) = (newIndices[1], newIndices[0])
        inList[-1].tensors[opPos] = sfExOp(newIndices)
        # print out the additional final term
        if options.verbose:
          print "2nd final term: ", inList[-1]

      # print a blank line
      if options.verbose:
        print ""

      # increment the index to the next term
      t_num += 1


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeVirtOps_sf(inList):
  """
  Removes from inList any terms containing a spin-free operator with a virtual index.
  """

  if options.verbose:
    print "removing terms containing a spin-free operator with a virtual index..."
    print ""

  # loop over the terms in inList
  i = 0
  while i < len(inList):

    # ensure that each element of inList is a term object
    if not isinstance(inList[i], term):
      raise TypeError, "inList must be a list of term objects"

    # determine if the term's spin-free excitation operators have any virtual indices
    hasVirtual = False
    for ten in inList[i].tensors:
      if isinstance(ten, sfExOp):
        for ind in ten.indices:
          if options.virtual_type in ind.indType:
            hasVirtual = True

    # remove the term if a spin-free excitation operator had a virtual index
    if hasVirtual:
      if options.verbose:
        print " removing term: ", inList[i]
      del inList[i]

    # otherwise, move to the next term
    else:
      i += 1

  if options.verbose:
    print ""


#--------------------------------------------------------------------------------------------------


