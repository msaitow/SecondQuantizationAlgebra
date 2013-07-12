#    file:  sqaTensor.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Defines the tensor class and its children.
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
from sqaSymmetry import symmetry


# The tensor class represents an object consisting of a name, an ordered set of indices,
# and possibly a set of symmetry permutations among the indices.
#
# A tensor's name is a string.
#
# The indices are given by a list of objects of the index class.
#
# The symmetries are given by a list of objects of the symmetry class.  Note that this
# list does not have to be exhaustive, as the code will attempt to find all possible
# symmetry permutations that can be created from the supplied symmetries.
#
# For example, the list of symmetries
#   [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2),-1), sqa.symmetry((1,0,3,2),1)]
# is redundant, as the third permutation can and will be generated using the first two.
#
#
# There are 4 important children of the tensor class:  creOp, desOp, sfExOp, and kroneckerDelta.
#  - creOp and desOp are one-index tensors representing creation and destruction operators.
#  - sfExOp is a 2n-index tensor representing a spin-free excitation operator.
#  - kroneckerDelta is a two-index tensor representing the Kronecker delta function.

# Modified by M. Saitow May 16, 2012, hasExt is added in kroneckerDelta and tensor 
# Modified by M. Saitow May 16, 2012, purelyExt is added in kroneckerDelta 

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

class tensor:
  "A class to represent tensors in operator algebra.  Integrals and density matrices are examples."

  #------------------------------------------------------------------------------------------------

  freelyCommutes = True

  #------------------------------------------------------------------------------------------------

  def __init__(self, name, indices = [], symmetries = []):

    # Initialize data
    (self.permutations,self.factors) = (None,None)
    self.indices = []
    self.symmetries = []

    # Process name
    self.name = str(name)

    # Process indices
    indicesError = "indices must be a list of index objects"
    if not isinstance(indices, type([])):
      raise TypeError, indicesError
    for i in indices:
      if not isinstance(i, index):
        raise TypeError, indicesError
      self.indices.append( i.copy() )

    # Process symmetries
    symmetryError = "symmetries must be a list of symmetry objects"
    if not isinstance(symmetries, type([])):
      raise TypeError, symmetryError
    for sym in symmetries:
      if not isinstance(sym, symmetry):
        raise TypeError, symmetryError
      for s in self.symmetries:
        if s.pattern == sym.pattern:
          raise ValueError, "a tensor cannot have two symmetries with the same pattern"
      self.symmetries.append( sym.copy() )

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):
    if (not isinstance(other,tensor)):
      raise TypeError, "A tensor may only be compared to another tensor"

    # If other belongs to a tensor subclass, use the subclass's comparison method
    if isinstance(other,kroneckerDelta) or \
       isinstance(other,creOp) or \
       isinstance(other,desOp) or \
       isinstance(other,creDesTensor) or \
       isinstance(other,sfExOp):
      return  - cmp(other,self)

    # compare name next
    retval = cmp(self.name,other.name)
    if retval != 0:
      return retval

    # compare indices next
    retval = cmp(self.indices,other.indices)
    if retval != 0:
      return retval

    # compare symmetries next
    retval = cmp(self.symmetries,other.symmetries)
    return retval

  #------------------------------------------------------------------------------------------------

  def __str__(self):
    retval = self.name + "("
    for i in range(len(self.indices)):
      retval += self.indices[i].name #+ " " + str(self.indices[i].type)
      if i < len(self.indices)-1:
        retval += ","
    retval += ")"
    return retval

  #------------------------------------------------------------------------------------------------

  def copy(self):
    "Returns a copy of the tensor"
    retval = tensor(self.name, self.indices, self.symmetries)
    if self.permutations != None and self.factors != None:
      retval.permutations = [ perm + [] for perm in self.permutations ]
      retval.factors = self.factors + []
    return retval

  #------------------------------------------------------------------------------------------------

  def symPermutes(self, force = False):
    "Returns the index permutations and resulting factors allowed by the tensor's symmetry"

#    # If the result is already known, return it
#    if not force and self.permutations != None and self.factors != None:
#      return (self.permutations,self.factors)
    
    # Otherwise, compute the permutations and corresponding factors
    tuples = [range(len(self.indices))]
    factors = [1]
    allFound = False
    while not allFound:
      allFound = True
      newTuples = []
      newFactors = []
      for j in range(len(tuples)):
        for sym in self.symmetries:
          newTuples.append([])
          newFactors.append(sym.factor * factors[j])
          for i in sym.pattern:
            newTuples[-1].append(tuples[j][i])
      while len(newTuples) > 0:
        isNew = True
        for tup in tuples:
          if tup == newTuples[0]:
            isNew = False
            break
        if isNew:
          allFound = False
          tuples.append(newTuples[0])
          factors.append(newFactors[0])
          #print tuples[-1]
        del(newTuples[0],newFactors[0])

    # Save the results for later so they don't need to be computed again
    self.permutations, self.factors = tuples, factors

    return (self.permutations,self.factors)

  #------------------------------------------------------------------------------------------------

  def sortIndeces(self):
    "Sorts the indices alphabetically within the constraints of symmetry and returns the resulting factor."

    # If the tensor has no symmetry, do nothing
    if len(self.symmetries) == 0:
      return 1

    # Get the permutation tuples allowed by symmetry and the corresponding factors
    tuples,factors = [],[]
    (tup,fac) = self.symPermutes()
    tuples.extend(tup)
    factors.extend(fac)

    # Score the different permutations and select the winner
    scores = []
    for tup in tuples:
      scores.append(0)
    for i in range(len(self.indices)-1):
      for j in range(i+1,len(self.indices)):
        for k in range(len(tuples)):
          if self.indices[tuples[k][i]] < self.indices[tuples[k][j]]:
            scores[k] += 1
      # At each iteration, determine the max score and keep only the tuples with that score
      maxScore = max(scores)
      i = 0
      while i < len(scores):
        if scores[i] < maxScore:
          del(scores[i],tuples[i],factors[i])
        else:
          i += 1
      if len(scores) == 0:
        break

    # Raise an error if a unique winner was not found
    if len(scores) > 1:
      unique = True
      for i in range(1,len(scores)):
        for j in range(len(tuples[i])):
          if self.indices[tuples[i][j]] != self.indices[tuples[0][j]]:
            unique = False
            break
        if not unique:
          break
      if not unique:
        print "No unique winner produced when sorting the indices of tensor %s" %(str(self))
        raise RuntimeError, "Scoring system did not produce unique winner."

    # Sort indices in the uniquely determined order and return the resulting factor
    newIndeces = []
    for i in range(len(tuples[0])):
      newIndeces.append(self.indices[tuples[0][i]])
    self.indices = newIndeces
    return factors[0]

  #------------------------------------------------------------------------------------------------

  def hasIndex(self,i):
    "Returns True if i is one of the tensor's indices and False otherwise."
    if not isinstance(i,index):
      raise TypeError, "i must be of the index class"
    return (i in self.indices)

  #------------------------------------------------------------------------------------------------

  def hasExt(self):
    "Returns True if this tensor object contains the dummy index"
    retval = False
    for i in self.indices:
      if i.isExt == True: retval = True
    return retval

  #------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


class kroneckerDelta(tensor):
  "A tensor representation of the kronecker delta function."

  #------------------------------------------------------------------------------------------------

  freelyCommutes = True
  symmetries = [symmetry((1,0),1)]
  name = "kdelta"

  #------------------------------------------------------------------------------------------------

  def __init__(self,indices):
    if len(indices) != 2:
      raise ValueError, "The kronecker delta function takes exactly two indices"
    self.indices = []
    for i in indices:
      self.indices.append( i.copy() )
    (self.permutations,self.factors) = (None,None)

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):

    # comparison to another kroneckerDelta
    if isinstance(other, kroneckerDelta):
      retval = cmp(self.name,other.name)
      if retval != 0:
        return retval
      retval = cmp(self.indices,other.indices)
      if retval != 0:
        return retval
      retval = cmp(self.symmetries,other.symmetries)
      return retval

    # kroneckerDelta class is less than the creOp, desOp, creDesTensor and sfExOp sub classes
    elif isinstance(other, creOp):
      return -1
    elif isinstance(other, desOp):
      return -1
    elif isinstance(other, creDesTensor):
      return -1
    elif isinstance(other, sfExOp):
      return -1

    # kroneckerDelta class is greater than other tensor classes
    elif isinstance(other, tensor):
      return 1

    # Raise error if other is not a tensor
    else:
      raise TypeError, "A kronekerDelta may only be compared to another tensor"
    return 0

  #------------------------------------------------------------------------------------------------

  def copy(self):
    "Returns a copy of the knroneckerDelta object"
    return kroneckerDelta(self.indices)

  #------------------------------------------------------------------------------------------------

#   def hasExt(self):
#     "Returns whether this Kronecker delta object contains the external index"
#     retval = False
#     for i in self.indices:
#       if i.isExt == True: retval = True
#     return retvall

  #------------------------------------------------------------------------------------------------

  def purelyExt(self):
    "Returns whether this Kronecker delta object is composed of purely the external indices"
    retval = False
    for i in self.indices:
      if i.isExt == True: retval = True
      else:               retval= False
    return retvall

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


class sfExOp(tensor):
  """
  A tensor representation of a spin free excitation operator.
  e.g. E(i,j,k,l) = sum over ( sigma, tau )  of ( a+_(i sigma) a+_(j tau) a_(l tau) a_(k sigma) )
  Note that the indices i,j,k,l refer to spacial orbitals and sigma,tau refer to spins.
  """

  #------------------------------------------------------------------------------------------------

  freelyCommutes = False

  #------------------------------------------------------------------------------------------------

  def __init__(self, indices):
    # Check that there are an even number of indices
    if len(indices)/2 != (len(indices)+1)/2:
      raise ValueError, "A spin free excitation operator (the sfExOp class) must have an " + \
                        "even number of indices"

    # Initialize order
    self.order = len(indices)/2

    # Initialize name
    self.name = "E%i" %self.order

    # Initialize indices
    self.indices = []
    for i in indices:
      self.indices.append( i.copy() )

    # Initialize permutations and factors
    (self.permutations,self.factors) = (None,None)

    # Initialize symmetries
    self.symmetries = []
    for i in range(self.order-1):
      if i == 0:
        temp_tup = (1,)
      else:
        temp_tup = (0,)
      for j in range(1,2*self.order):
        if j == i:
          temp_tup = temp_tup + (i+1,)
        elif j == i+1:
          temp_tup = temp_tup + (i,)
        elif j == i+self.order:
          temp_tup = temp_tup + (i+1+self.order,)
        elif j == i+1+self.order:
          temp_tup = temp_tup + (i+self.order,)
        else:
          temp_tup = temp_tup + (j,)
      self.symmetries.append(symmetry(temp_tup, 1))
          

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):

    # comparison to another sfExOp
    if isinstance(other,sfExOp):
      retval = cmp(self.name,other.name)
      if retval != 0:
        return retval
      retval = cmp(self.indices,other.indices)
      if retval != 0:
        return retval
      retval = cmp(self.symmetries,other.symmetries)
      return retval

    # sfExOp class is less than the creOp, desOp, and creDesTensor classes
    elif isinstance(other,creOp):
      return -1
    elif isinstance(other,desOp):
      return -1
    elif isinstance(other,creDesTensor):
      return -1

    # sfExOp class is greater than other tensor subclasses
    elif isinstance(other,tensor):
      return 1

    # raise an error if other is not a tensor
    else:
      raise TypeError, "An sfExOp object may only be compared to another tensor"
    return 0

  #------------------------------------------------------------------------------------------------

  def copy(self):
    return sfExOp(self.indices)

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

class creDesTensor(tensor):
  """
  A tensor representation of a string of creation/destruction operators
  """

  #------------------------------------------------------------------------------------------------

  freelyCommutes = False
  name = "creDesTensor"

  #------------------------------------------------------------------------------------------------

  def __init__(self, ops):
    
    TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
    if not type(ops) == type([]):
      raise TypeError, TypeErrorMessage

    # Initialize permutations and factors
    (self.permutations,self.factors) = (None,None)

    # Build the index list and count the number of creation operators
    self.nCre = 0
    self.indices = []
    desFlag = False
    for op in ops:
      if (not isinstance(op, creOp)) and (not isinstance(op, desOp)):
        raise TypeError, TypeErrorMessage
      if isinstance(op, desOp):
        desFlag = True
      if isinstance(op, creOp):
        self.nCre += 1
        if desFlag:
          raise TypeError, TypeErrorMessage
      self.indices.append(op.indices[0].copy())

    # Initialize symmetries
    self.symmetries = []
    swapValues = range(len(self.indices)-1)
    if self.nCre > 0:
      del(swapValues[self.nCre-1])
    for i in swapValues:
      if i == 0:
        temp_tup = (1,)
      else:
        temp_tup = (0,)
      for j in range(1,len(self.indices)):
        if j == i:
          temp_tup = temp_tup + (i+1,)
        elif j == i+1:
          temp_tup = temp_tup + (i,)
        else:
          temp_tup = temp_tup + (j,)
      self.symmetries.append(symmetry(temp_tup, -1))
          
  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):

    # comparison to another creDesTensor
    if isinstance(other, creDesTensor):
      retval = cmp(self.name,other.name)
      if retval != 0:
        return retval
      retval = cmp(self.indices,other.indices)
      if retval != 0:
        return retval
      retval = cmp(self.symmetries,other.symmetries)
      return retval

    # creDesTensor class is less than the creOp and desOp classes
    elif isinstance(other,creOp):
      return -1
    elif isinstance(other,desOp):
      return -1

    # creDesTensor class is greater than other tensor subclasses
    elif isinstance(other,tensor):
      return 1

    # raise an error if other is not a tensor
    else:
      raise TypeError, "A creDesTensor object may only be compared to another tensor"
    return 0

  #------------------------------------------------------------------------------------------------

  def copy(self):
    ops = []
    for i in range(self.nCre):
      ops.append(creOp(self.indices[i]))
    for i in range(self.nCre,len(self.indices)):
      ops.append(desOp(self.indices[i]))
    return creDesTensor(ops)

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


class creOp(tensor):
  """
  A tensor representation for a creation operator
  """

  #------------------------------------------------------------------------------------------------

  freelyCommutes = False
  name = "cre"
  symmetries = []

  #------------------------------------------------------------------------------------------------

  def __init__(self, indices):

    # Initialize index
    if type(indices) == type([]) and len(indices) == 1 and isinstance(indices[0], index):
      inputIndex = indices[0].copy()
    elif isinstance(indices, index):
      inputIndex = indices.copy()
    else:
      raise TypeError, "indices must be an index or a list of indices with length 1"
    self.indices = [inputIndex]

    # Initialize permutations and factors
    (self.permutations,self.factors) = (None,None)

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):

    # comparison to another creOp
    if isinstance(other,creOp):
      retval = cmp(self.name,other.name)
      if retval != 0:
        return retval
      retval = cmp(self.indices,other.indices)
      if retval != 0:
        return retval
      retval = cmp(self.symmetries,other.symmetries)
      return retval

    # creOp class is less than the desOp class
    elif isinstance(other,desOp):
      return -1

    # creOp class is greater than other tensor subclasses
    elif isinstance(other,tensor):
      return 1

    # raise an error if other is not a tensor
    else:
      raise TypeError, "An creOp object may only be compared to another tensor"
    return 0

  #------------------------------------------------------------------------------------------------

  def copy(self):
    return creOp(self.indices)

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


class desOp(tensor):
  """
  A tensor representation for a destruction operator
  """

  #------------------------------------------------------------------------------------------------

  freelyCommutes = False
  name = "des"
  symmetries = []

  #------------------------------------------------------------------------------------------------

  def __init__(self, indices):

    # Initialize index
    if type(indices) == type([]) and len(indices) == 1 and isinstance(indices[0], index):
      inputIndex = indices[0].copy()
    elif isinstance(indices, index):
      inputIndex = indices.copy()
    else:
      raise TypeError, "indices must be an index or a list of indices with length 1"
    self.indices = [inputIndex]

    # Initialize permutations and factors
    (self.permutations,self.factors) = (None,None)

  #------------------------------------------------------------------------------------------------

  def __cmp__(self,other):

    # comparison to another desOp
    if isinstance(other,desOp):
      retval = cmp(self.name,other.name)
      if retval != 0:
        return retval
      retval = cmp(self.indices,other.indices)
      if retval != 0:
        return retval
      retval = cmp(self.symmetries,other.symmetries)
      return retval

    # desOp class is greater than other tensor subclasses
    elif isinstance(other,tensor):
      return 1

    # raise an error if other is not a tensor
    else:
      raise TypeError, "An desOp object may only be compared to another tensor"
    return 0

  #------------------------------------------------------------------------------------------------

  def copy(self):
    return desOp(self.indices)

  #------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

