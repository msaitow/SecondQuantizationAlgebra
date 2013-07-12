#    file:  sqaCommutator.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Computes the commutator, including normal ordering,
#           between terms or lists of terms.
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

from sqaTerm import term, combineTerms, multiplyTerms, termChop
from sqaNormalOrder import normalOrder

def commutator(leftInput, rightInput, contract = True, combine = True):

  # Convert inputs that are terms into lists of terms
  if isinstance(leftInput,term):
    leftTerms = [leftInput]
  else:
    leftTerms = leftInput
  if isinstance(rightInput,term):
    rightTerms = [rightInput]
  else:
    rightTerms = rightInput

  # Check input integrity
  TypeErrorMessage = "commutator inputs must be terms or lists of terms"
  if type(leftTerms) != type([]) or type(rightTerms) != type([]):
      raise TypeError, TypeErrorMessage
  for t in rightTerms:
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage
  for t in leftTerms:
    if not isinstance(t, term):
      raise TypeError, TypeErrorMessage

  # Construct terms resulting from the commutator
  preNOTerms = [] #terms before normal ordering
  for lterm in leftTerms:
    for rterm in rightTerms:
      preNOTerms.append( multiplyTerms(lterm,rterm) )
      preNOTerms.append( multiplyTerms(rterm,lterm) )
      preNOTerms[-1].scale(-1)

  # For each term, apply Wick's theorem to convert it to normal order
  noTerms = [] #terms after normal ordering
  for t in preNOTerms:
    noTerms.extend(normalOrder(t))
  del(preNOTerms)

  # Contract any delta functions resulting from the normal ordering,
  # unless told not to
  if contract:
    for t in noTerms:
      t.contractDeltaFuncs()

  # Remove any terms that are zero 
  termChop(noTerms)

  # Combine any like terms unless told not to
  if combine:
    combineTerms(noTerms)

  # Return result
  return noTerms
