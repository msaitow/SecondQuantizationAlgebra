#    file:  sqaIndex.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Defines the index class.
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




# The index class is used to represent a tensor index.
#
# An index consists of three parts:  a name, a list of types, and a bool.
#
# The index's name may be any string.
#
# The bool indicates whether the index is summed over or otherwise is a dummy
# index whose name may be changed.
# A value of True means that the various operator algebra functions are allowed
# the change the index's name as they wish.
# A value of False means that the index's name may not be changed.
#
# The list of types is actually a list of lists of strings.
# Each list of strings represents a type group.
# Examples of type groups are the index's spin type or whether the index is 
# core, active, virtual, etc.
# The reason for this format is to evaluate Kronecker delta functions, which will
# evaluate to zero if there is no overlap in any of the two indices' type groups.
# For example, if index1 has indexType = [['alpha'], ['core', 'active']] and
# index2 has indexType = [['beta'], ['active']], the Kronecker delta between them
# will be zero because the first type group has no matching strings.
# If the first type group had been omitted, then the delta function would not be
# zero because both indices have 'active' as one of the types in the second type group.
# Note that while the type groups can be inputted as a list of lists of strings, they
# are actually stored as a tuple of tuples of strings.

# Modified by M. Saitow May 16, 2012, isExt flag and __hash__ are added

class index(object):
  "A class for tensor and operator indices."

  def __init__(self, name, indexType = (), isSummed = False, isExt = False):
    # Initialize index name
    self.name = str(name)

    # Initialize index types
    indType = []
    for l in indexType:
      if not ( type(l) in [type([]), type(())] ):
        raise TypeError, "indexType must be a list or tuple of lists or tuples of strings"
      indType.append([])
      indType[-1].extend(l)
      indType[-1].sort()
#    indType.sort()
    self.indType = []
    for l in indType:
      self.indType.append(())
      for s in l:
        if type(s) != type('a'):
          raise TypeError, "indexType must be a list or tuple of lists or tuples of strings"
        self.indType[-1] = self.indType[-1] + (s,)
    self.indType = tuple(self.indType)

    # Initialize flag for whether the index is summed over
    if type(isSummed) != type(True):
      raise TypeError, "isSummed must be True or False"
    self.isSummed = isSummed

    # Initialize flag for whether the index is external one
    if type(isExt) != type(True):
      raise TypeError, "isExt must be True or False"
    self.isExt = isExt

  def __cmp__(self,other):
    if (not isinstance(other,index)):
      raise ValueError, "can only compare index class with other index class objects."
#     retval = cmp(self.isExt, other.isExt)
#     if retval != 0:
#       return retval
    retval = cmp(self.isSummed, other.isSummed)
    if retval != 0:
      return retval
    retval = cmp(self.name, other.name)
    if retval != 0:
      return retval
    return cmp(self.indType, other.indType)

  def __hash__(self):
    "Returns a hash tag. This implementation is a sort of ad hoc. So, this may cause some problem ..."
    retval = 0
    for c in self.name:
      retval = retval ^ ord(c)
    return retval

  def tup(self):
    "Returns a tuple representation of the index.  The return object in unmutable and thus can be used as a dictionary key."
    return (self.name, self.indType, self.isSummed, self.isExt)

  def copy(self):
    "Returns a deep copy of the index"
    return index(self.name, self.indType, self.isSummed, self.isExt)
