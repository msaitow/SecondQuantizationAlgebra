#    file:  sqaConvert.py
#  author:  Masaaki Saitow
#    date:  May 16, 2012
# summary:  Converts a stream of terms into normal Orz code.
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
#
# Also, 
#   Masaaki Saitow, Yuki Kurashige, and Takeshi Yanai.
#   J. Chem. Phys. ***, ****** (2013)
#


from math import fabs, log
import datetime
from sqaTensor import tensor, kroneckerDelta, creOp, desOp, sfExOp
from sqaTerm import term, sortOps
from sqaIndex import index
from sqaMisc import makeTuples, allDifferent, makePermutations, get_num_perms
from sqaOptions import options
from sqaConvert2 import CPloop, ReadBare, ReadERI, ReadD4, AccBare, LoopEnd, ReadRetval

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convert2Mulliken(inTerm, name_ofAmp = 'T'):
  "Return a list of terms in the Mulliken's notation, which is better suited for the Orz code"

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  retval = []
  # Sequentially process each term provided 
  for numTerm in range(len(inTerm)):

    terms = []
    for numTensor in range(len(inTerm[numTerm].tensors)):
      theTensor = inTerm[numTerm].tensors[numTensor]
      if theTensor.name == name_ofAmp: # Skip the conversion if thisTensor is the Amplitude ....
        terms.append(theTensor)
        continue 

      # In case of Kronecker delta or the one-body integrals ...
      if len(theTensor.indices) <= 2:
        terms.append(theTensor)
      else:
#       print "*** " + theTensor.name
#       print len(theTensor.indices)/2
        indices = []
        distance = len(theTensor.indices) / 2
        for numIndex in range(len(theTensor.indices)/2):
          indices.append(theTensor.indices[numIndex])
          indices.append(theTensor.indices[numIndex+distance])

        # Permutational symmetry is lost in the resultant tensor
        terms.append(tensor(theTensor.name, indices))
    retval.append(term(inTerm[numTerm].numConstant, inTerm[numTerm].constants, terms))

  return retval


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_0(LTensor, inTerm):
  "Returns a stream of Fortran code, which evaluate the tensorial contraction encountered in the Sigma equation."
  "In this function, existence of the outer iteration in C++ level is not allowed"
  
  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

    # If term has both creation/destruction operators and spin free excitation operators,
    # raise an error
    if has_creDesOps:
      raise RuntimeError, "creOp/desOp cannot be implemenented in ORz code as tensorial quantities" + \
                          "creOp/desOp are present " + "in term " + num

  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  print "!"
  print "!  Restate my assumptions:                                                   " 
  print "!  1, Mathematics is the language of nature.                                 " 
  print "!  2, Everything around us can be represented and understood through numbers." 
  print "!  3, If you graph the numbers of any system, patterns emerge.               " 
  print "!  Therefore: There are patterns everywhere in nature.                       " 
  print "!                                      ~~ Maximillian Cohen in Pi (film) ~~  " 
  print "!"

  # First of all, make the list of all the variables used in the subsequent tensorial contraction
  union_ofIndices = set(LTensor.indices)
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
      setIndices = set(thisTensor.indices)
      union_ofIndices = setIndices.union(union_ofIndices)

  list_ofIndices = list(union_ofIndices)
  print "! Indices used in the contractions below ... "
  for num in range(len(list_ofIndices)):
    thisIndex = list_ofIndices[num]
    if num % 5 != 0: print ", s_" + thisIndex.name + ", i_" + thisIndex.name, 
    else:
      print ""
      print "real(kind=8) :: ",
      print "s_" + thisIndex.name + ", i_" + thisIndex.name, 

  print ""
  print ""

  print ""
  # Search in each term ...
  num_loops    = []
  summed_Index = []
  for numTerm in range(len(inTerm)):
    thisTerm = inTerm[numTerm]
    print "! No.", numTerm
    print "!", LTensor, " <-- "
    print "!", thisTerm

    # Search inbetween each tensor ...
    loop_count = 0
    num_kdeltas  = 0
    summed_Index = LTensor.indices
    for numTensor in range(len(thisTerm.tensors)-1):
      thisTensor = thisTerm.tensors[numTensor]
      nextTensor = thisTerm.tensors[numTensor+1]

      set_i1       = set(thisTensor.indices)
      set_i2       = set(nextTensor.indices)
      set_i3       = set_i1.union(set_i2)
      summed_Index = set_i3.union(summed_Index)

    # Write down the loops over irreps ...
    for index in summed_Index:
      print "do s_" + index.name + " = 1,nir"
      loop_count += 1

    # Write down the constraints in terms of the irreps ...
    # Conditions for the left-hand tensor
    print "if( &"
    thisTensor = LTensor
    if len(thisTensor.indices) == 1:   # in case of unity
      print "S_" + thisTensor.indices[0].name + " == 0",
    elif len(thisTensor.indices) == 2: # in case of two
      print "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0",
    elif len(thisTensor.indices) == 3: # in case of three
      str  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      str += "s_" + thisTensor.indices[2].name + ") == 0"
      print str,
    elif len(thisTensor.indices) == 4: # in case of four
      str  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      str += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
      print str,
    elif len(thisTensor.indices) == 5: # in case of five
      str  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      str += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
      str += "s_" + thisTensor.indices[4].name + ")"
      print str,
    elif len(thisTensor.indices) == 6: # in case of six
      str  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      str += "s_" + thisTensor.indices[2].name + ") == "
      str += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
      str += "s_" + thisTensor.indices[5].name + ")"
      print str,
    else:
      raise TypeError, "Number of indices in tensor must be 1 - 6"
    print " & " 
    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if len(thisTensor.indices) == 1:   # in case of unity
        print "S_" + thisTensor.indices[0].name + " == 0",
      elif len(thisTensor.indices) == 2: # in case of two
        print "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0",
      elif len(thisTensor.indices) == 3: # in case of three
        str  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        str += "s_" + thisTensor.indices[2].name + ") == 0"
        print str,
      elif len(thisTensor.indices) == 4: # in case of four
        str  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        str += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
        print str,
      elif len(thisTensor.indices) == 5: # in case of five
        str  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        str += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
        str += "s_" + thisTensor.indices[4].name + ")"
        print str,
      elif len(thisTensor.indices) == 6: # in case of six
        str  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        str += "s_" + thisTensor.indices[2].name + ") == "
        str += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
        str += "s_" + thisTensor.indices[5].name + ")"
        print str,
      else:
        raise TypeError, "Number of indices in tensor must be 1 - 6"

      if not numTensor == len(thisTerm.tensors) - 1:
        print " .and. &"

    print ") then"
    
    # Write down the body of the tensorial contraction ...
    # *CAUTION* Only active and virtual orbitals are considered separately for now 
    count = 0
    for Index in summed_Index:
      iname = "i_" + Index.name
      sname = "s_" + Index.name
      if   Index.indType[0][0] == 'active':
        otype = "I_O"
      elif Index.indType[0][0] == 'virtual':
        otype = "I_V"
      print iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")"
      count += 1
    if count != loop_count:
      raise TypeError, "Algorithmic error occured .... "

    # As a first step of evalutation of the contraction, search the Kronecker deltas    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': num_kdeltas += 1

    if num_kdeltas >= 1:
      count = num_kdeltas
#      print "num_kdeltas ", num_kdeltas # *TEST*
      print "if(",
      for numTensor in range(len(thisTerm.tensors)):
        thisTensor = thisTerm.tensors[numTensor]
        if thisTensor.name == 'kdelta':
          str  = "s_" + thisTensor.indices[0].name + " == s_" + thisTensor.indices[1].name
          str += " .and. " + "i_" + thisTensor.indices[0].name + " == i_" + thisTensor.indices[1].name
          count -= 1
          if   count == 0: str += ") then"
          elif count >= 1: str += " .and. & "
          elif count >  0: raise TYpeError, "Algorithmic Error"
          print str

    print ""
    # First generate the left-hand side tensor
    nterm = LTensor.name + "_("
    for n_index in range(len(LTensor.indices)):
      nterm += "s_" + LTensor.indices[-(n_index+1)].name
      if not n_index == len(LTensor.indices)-1: nterm += ", "
      else: nterm += ")"
    nterm += "%array("
    for n_index in range(len(LTensor.indices)):
      nterm += "i_" + LTensor.indices[-(n_index+1)].name
      if not n_index == len(LTensor.indices)-1: nterm += ", "
      else: nterm += ")"
    print nterm + " = " + nterm + " &" 

    # Generate all the tensors ...
    if thisTerm.numConstant > 0.0: print "  + ", thisTerm.numConstant, " & "
    else:                          print "  - ", thisTerm.numConstant, " & " 

    if len(constList) > 0: 
      for thisConst in constList: print " * " + thisConst + " & "

    # Generate each tensor except Kronecker deltas ....
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': continue    
      nterm = "  * " + thisTensor.name + "_("
      for n_index in range(len(thisTensor.indices)):
        nterm += "s_" + thisTensor.indices[-(n_index+1)].name
        if not n_index == len(thisTensor.indices)-1: nterm += ", "
        else: nterm += ")"
      nterm += "%array("
      for n_index in range(len(thisTensor.indices)):
        nterm += "i_" + thisTensor.indices[-(n_index+1)].name
        if not n_index == len(thisTensor.indices)-1: nterm += ", "
        else: nterm += ")"
      if numTensor < len(thisTerm.tensors)-1: nterm += " & " 
      print nterm 
    
    print ""

    # Closing the iteration .... 
    if num_kdeltas >= 1: print "end if"
    for i in range(loop_count): print "end do"
    print "end if"
    for i in range(loop_count): print "end do"
    print "! Loop_count : ", loop_count
    print ""


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_ERI(LTensor, inTerm, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Separate given many-body equations in terms of the types of ERIs contained i.e. V[k*|**], V[a*|**], and 4-RDM2."
  "Then, interface to convertF90 function"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORZ code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  print ""
  print "! Searching in terms .... "
  print ""
  TermLeft    = []     # List composed of terms with neither ERI, nor 4-RDMs
  TermVc      = []     # List composed of terms that contain Vc type ERI, except those with 4-RDMs
  TermVi      = []     # List composed of terms that contain Vi type ERI, except those with 4-RDMs
  TermVa      = []     # List composed of terms that contain Va type ERI, except those with 4-RDMs
  TermVc_Rdm4 = []     # List composed if terms that contain both Vc and 4-RDMs 
  TermVi_Rdm4 = []     # List composed if terms that contain both Vi and 4-RDMs 
  TermVa_Rdm4 = []     # List composed if terms that contain both Va and 4-RDMs 
  isVc        = False  # Flag whether a term contains the Vi ERI
  isVi        = False  # Flag whether a term contains the Vi ERI
  isVa        = False  # Flag whether a term contains the Va ERI
  isRdm4      = False  # Flag wgether a term contains the 4-RDMs
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
      if   thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'core':    isVc = True
      elif thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'active':  isVi = True
      elif thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'virtual': isVa = True 
      elif thisTensor.name == name_ofRdm4: isRdm4 = True
    if   isVi and     isRdm4: TermVc_Rdm4.append(thisTerm)
    elif isVi and not isRdm4: TermVc.append(thisTerm) 
    elif isVi and     isRdm4: TermVi_Rdm4.append(thisTerm)
    elif isVi and not isRdm4: TermVi.append(thisTerm) 
    elif isVa and     isRdm4: TermVa_Rdm4.append(thisTerm)
    elif isVa and not isRdm4: TermVa.append(thisTerm);
    elif not isVi and not isVa and not isRdm4:                     
                              TermLeft.append(thisTerm)
    isVi   = False
    isVa   = False
    isRdm4 = False

  print "!  **** SUMMARY *****"
  print "! 0. Number of terms with neither ERI, nor 4-RDMs : " + str(len(TermLeft))
  print "! 1. Number of terms with V[c*|**] and no 4-RDM   : " + str(len(TermVc))
  print "! 2. Number of terms with V[c*|**] and 4-RDM      : " + str(len(TermVc_Rdm4))
  print "! 3. Number of terms with V[i*|**] and no 4-RDM   : " + str(len(TermVi))
  print "! 4. Number of terms with V[i*|**] and 4-RDM      : " + str(len(TermVi_Rdm4))
  print "! 5. Number of terms with V[a*|**] and no 4-RDM   : " + str(len(TermVa))
  print "! 6. Number of terms with V[a*|**] and 4-RDM      : " + str(len(TermVa_Rdm4))
  print "! ** Number of Total terms                        : " + str(len(inTerm))
  print ""
  if len(TermLeft) + len(TermVc) + len(TermVi) + len(TermVa) + len(TermVc_Rdm4) + len(TermVi_Rdm4) + len(TermVa_Rdm4) != len(inTerm):
    raise TypeError, "Algorithmic error .... "

  if len(TermLeft) > 0:
    print ""
    print "! Converting terms with neither ERI, nor 4-RDMs ..... "
    convertF90_1(LTensor, TermLeft, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVc) > 0:
    print ""
    print "! Converting terms with V[c*|**] and no 4-RDMs ..... "
    convertF90_1(LTensor, TermVc, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVc_Rdm4) > 0:
    print ""
    print "! Converting terms with V[c*|**] with 4-RDMs ..... "
    convertF90_1(LTensor, TermVc_Rdm4, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVi) > 0:
    print ""
    print "! Converting terms with V[i*|**] and no 4-RDMs ..... "
    convertF90_1(LTensor, TermVi, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVi_Rdm4) > 0:
    print ""
    print "! Converting terms with V[i*|**] with 4-RDMs ..... "
    convertF90_1(LTensor, TermVi_Rdm4, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVa) > 0:
    print ""
    print "! Converting terms with V[a*|**] and no 4-RDMs ..... "
    convertF90_1(LTensor, TermVa, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"
  if len(TermVa_Rdm4) > 0:
    print ""
    print "! Converting terms with V[a*|**] with 4-RDMs ..... "
    convertF90_1(LTensor, TermVa_Rdm4, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    print "! -----------------------------------------------------------------------"


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_1(LTensor, inTerm, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Returns a stream of Fortran code, which evaluate the tensorial contraction encountered in the Sigma equation."
  "In this function, existence of the outer iteration in C++ level is allowed but the names of the tensors such as "
  "ERI or BareAmp have to be specified, as well as the external index contained in the left-hand side tensor."
  "The Kronecker delta composed of purely the external indices are still evalusted at the most inner level"

  # In case of that the LHS is a Scalar 
  if len(inTerm) == 0:
    print "! RHS is ZERO .... "
    return

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORz code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # check whether the LTensor could be a BareAmp if specified so ....
  if len(LTensor.indices) != 4 and isBareLTS:
    raise TypeError, "LTensor could not be a BareAmpPack .... "

  print "!"
  print "!  Restate my assumptions:                                                   " 
  print "!  1, Mathematics is the language of nature.                                 " 
  print "!  2, Everything around us can be represented and understood through numbers." 
  print "!  3, If you graph the numbers of any system, patterns emerge.               " 
  print "!  Therefore: There are patterns everywhere in nature.                       " 
  print "!                                      ~~ Maximillian Cohen in Pi (film) ~~  " 
  print "!"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  print ""
  print "! ---------------------------- Parameters used ------------------------------"
  print "!"
  print "! Whether the LHS is a BareAmpPack ....... ", isBareLTS
  print "! LHS name ............................... ", lts_name
  print "! Name of ERI ............................ ", name_ofEri
  print "! Name of BareAmpPack appearing in RHS.... ", name_ofBareAmp
  print "!"
  print "!----------------------------------------------------------------------------"  
  print ""

  d = datetime.datetime.today()
  print "! Generated time: " + str(d.year) + "/" + str(d.month) + "/" + str(d.day),
  print " " + str(d.hour) + ":" + str(d.minute) + ":" + str(d.second)
  print ""

  print "! FEMTO BEGIN  **************************************************************" 
  print "use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock6"
  print ""
  print "implicit none"

  # Print all the external variables ....
  print ""
  isBare = False
  isEri  = False
  isRdm4 = False
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
      if   thisTensor.name == name_ofEri:     isEri  = True 
      elif thisTensor.name == name_ofBareAmp: isBare = True
      elif thisTensor.name == name_ofRdm4:    isRdm4 = True
  print ""

  if isBareLTS: 
    print "! Left tensor is BareAmp type ...."
    print "integer, intent(in) :: " + "i" + lts_name + ", s" + lts_name
  if isBare: 
    print "! Amplitude is found .... " 
    print "integer, intent(in) :: iamp, samp" # This can be different
  if isEri:
    print "! ERI is found .... "     
    print "integer, intent(in) :: ieri, seri"
  if isRdm4:
    print "! RDM4 is found .... "    
    print "integer, intent(in) :: irdm1, srdm1, irdm2, srdm2"

  print "! Information of the Irreps .... "
  print "integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)"  

  # Print all the name of parameters .... 
  list_ofAllConst = []
  for thisTerm in inTerm:
    for thisConst in thisTerm.constants:
      list_ofAllConst.append(thisConst)

  list_ofAllConst = list(set(list_ofAllConst))
  if len(list_ofAllConst) != 0:
    print "! Declaration of numerical constants .... "
  for thisConst in list_ofAllConst:
    print "real(kind=8), intent(inout) :: " + str(thisConst)

  # Print all the names of tensors .....
  if isBareLTS: 
    dict_ofAllTensors = {LTensor.name: 3}
  else:
    dict_ofAllTensors = {LTensor.name: len(LTensor.indices)}
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
#      print thisTensor.name # *TEST*
      if not dict_ofAllTensors.has_key(thisTensor.name):
        if thisTensor.name == name_ofBareAmp or thisTensor.name == name_ofEri:
          dict_ofAllTensors[thisTensor.name] = 3
        elif thisTensor.name == name_ofRdm4:
          dict_ofAllTensors[thisTensor.name] = 6
        else:
          dict_ofAllTensors[thisTensor.name] = len(thisTensor.indices)

#  global str
  print ""
  print "! Declare tensors used ..."
  print ""
  for name in dict_ofAllTensors.keys():
    if name == 'kdelta': continue
    dimTensor = dict_ofAllTensors[name]
    if dimTensor == 0:
      print "real(kind=8)                   :: " + name + "_"
      continue
    else:
      decTensor = "type(symblock" + str(dimTensor) + "), intent(inout) :: " + name + "_("    
    for i in range(dimTensor): 
      decTensor+= "0:nir-1"
      if i == dimTensor-1: decTensor += ")"
      else:                decTensor += ", "
    print decTensor
  print ""

#  print dict_ofAllTensors # *TEST*

  # First of all, make the list of all the variables used in the subsequent tensorial contraction
  list_ofIndices = []
  for thisIndex in LTensor.indices: list_ofIndices.append(thisIndex.name)
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
#      print "********* " + thisTensor.name # *TEST*
      for thisIndex in thisTensor.indices:
#        print thisIndex.name # *TEST*
        list_ofIndices.append(thisIndex.name)

  union_ofIndices = set(list_ofIndices)
  list_ofIndices = list(union_ofIndices)
  #print list_ofIndices # *TEST*
  print "! Indices used in the contractions as dummy ... "
  for num in range(len(list_ofIndices)):
    thisName = list_ofIndices[num]
    if num % 5 != 0: print ", s_" + thisName + ", i_" + thisName, 
    else:
      print ""
      print "integer :: ",
      print "s_" + thisName + ", i_" + thisName, 

  print ""
  print ""

  # If the LHS is a BareAmp, set so ...
  if isBareLTS == True: 
    LTensor.indices[extBare].isExt = True

  print ""
  # Search in each term ...
  num_loops    = []
#  summed_Index = []
  for numTerm in range(len(inTerm)):
    thisTerm = inTerm[numTerm]
    print "! No.", numTerm
    print "!", LTensor, " <-- "
    print "!", thisTerm

    # Dictionary that shows which index corresponds to ieri, isig, and so on ...
    extDict = {}
    # List of all external indices without any redundancy
    extIndices = []
    # If the LHS is a BareAmp, set so ...
    if isBareLTS == True: 
      extIndices.append(LTensor.indices[extBare]) # Record which indices are of external
      extDict["i"+lts_name] = LTensor.indices[extBare]

    for thisTensor in thisTerm.tensors:
      # If thisTensor is ERI, or BareAmp, the external indices have to be dropped off ...
      if   thisTensor.name == name_ofEri:     
        thisTensor.indices[extEri].isExt  = True
        extIndices.append(thisTensor.indices[extEri])
        extDict["ieri"] = thisTensor.indices[extEri]
      elif thisTensor.name == name_ofBareAmp: 
        thisTensor.indices[extBare].isExt = True
        extIndices.append(thisTensor.indices[extBare])
        extDict["iamp"] = thisTensor.indices[extBare] # This can be different one
      elif thisTensor.name == name_ofRdm4:
        thisTensor.indices[extRdm4[0]].isExt = True
        extIndices.append(thisTensor.indices[extRdm4[0]])
        extDict["irdm1"] = thisTensor.indices[extRdm4[0]]
        thisTensor.indices[extRdm4[1]].isExt = True
        extIndices.append(thisTensor.indices[extRdm4[1]])
        extDict["irdm2"] = thisTensor.indices[extRdm4[1]]
    # Avoid double counting     
    extIndices = list(set(extIndices)) 
    # Set all the external indices commonly beared by different tensors ....
    for extIndex in extIndices:
#      print extIndex.name, extIndex.isExt # *TEST*
      for thisTensor in thisTerm.tensors:
        for thisIndex in thisTensor.indices:
          if extIndex.name == thisIndex.name:
#            print thisIndex.name 
            thisIndex.isExt = True

#     # TESTS BEGIN ........
#     print "*********************"
#     for thisTensor in thisTerm.tensors:
#       print ""
#       for thisIndex in thisTensor.indices:
#         print thisIndex.name, thisIndex.isExt
#     print "*********************"
#     # TESTS END ..........

#     # Remove the internal loop if the term has Kronecker delta with the external and internal indices
#     for thisTensor in thisTerm.tensors:
#       if isinstance(thisTensor, kroneckerDelta) and thisTensor.hasExt and not thisTensorPurelyExt:      

    # Search inbetween each tensor ...
    loop_count  = 0
    num_kdeltas = 0
    summed_Index = LTensor.indices
    for numTensor in range(len(thisTerm.tensors)-1):
      thisTensor = thisTerm.tensors[numTensor]
      nextTensor = thisTerm.tensors[numTensor+1]

      set_i1       = set(thisTensor.indices)
      set_i2       = set(nextTensor.indices)
      set_i3       = set_i1.union(set_i2)
      summed_Index = set_i3.union(summed_Index)

    # Print whether the leg of ERI belongs to {occ} or {vir}
    # Core is not considered for now
    for thisTensor in thisTerm.tensors:
      if thisTensor.name == name_ofEri:
        print "! where i_" + thisTensor.indices[extEri].name + " in ",
        if   thisTensor.indices[extEri].indType[0][0] == 'core':
          print "{core}"
        elif thisTensor.indices[extEri].indType[0][0] == 'active':
          print "{occ}"
        elif thisTensor.indices[extEri].indType[0][0] == 'virtual':
          print "{vir}"
        else: raise TypeError, "Algorithmic Error"

    # Print guard block from the inappropriate call of external loop (lts_name is fixed to "sim")
    numGuard = 0

    # *NEW* I think this works more generally ..... 
    for numKey in range(len(extDict)):
      thisKey = extDict.keys()[numKey]
      for num in range(numKey+1,len(extDict)):
        thatKey = extDict.keys()[num]
        if extDict[thisKey].name == extDict[thatKey].name: 
          print "if(" + thisKey + " == " + thatKey + ") then"
          numGuard += 1
    # *NEW*

#     # *OLD*
#     if extDict.has_key('isig') and extDict.has_key('iamp'):    # in case isig == iamp 
#       if extDict["isig"].name == extDict["iamp"].name: 
#         print "if(isig == iamp .and. ssig == samp) then"
#         numGuard += 1
#     if extDict.has_key("isig") and extDict.has_key("ieri"):    # in case isig == ieri 
#       if extDict["isig"].name == extDict["ieri"].name:
#         print "if(isig == ieri .and. ssig == seri) then"
#         numGuard += 1
#     if extDict.has_key("isig") and extDict.has_key("irdm1"):   # in case isig == irdm1 
#       if extDict["isig"].name == extDict["irdm1"].name:
#         print "if(isig == irdm1 .and. ssig == srdm1) then"
#         numGuard += 1
#     if extDict.has_key("isig") and extDict.has_key("irdm2"):   # in case isig == irdm2 
#       if extDict["isig"].name == extDict["irdm2"].name:
#         print "if(isig == irdm2 .and. ssig == srdm2) then"
#         numGuard += 1
#     if extDict.has_key("iamp") and extDict.has_key("ieri"):    # in case iamp == ieri
#       if extDict["iamp"].name == extDict["ieri"].name:
#         print "if(iamp == ieri .and. samp == seri) then"
#         numGuard += 1
#     if extDict.has_key("iamp") and extDict.has_key("irdm1"):   # in case iamp == irdm1 
#       if extDict["iamp"].name == extDict["irdm1"].name:
#         print "if(iamp == irdm1 .and. samp == srdm1) then"
#         numGuard += 1
#     if extDict.has_key("iamp") and extDict.has_key("irdm2"):   # in case iamp == irdm2 
#       if extDict["iamp"].name == extDict["irdm2"].name:
#         print "if(iamp == irdm2 .and. samp == srdm2) then"
#         numGuard += 1
#     if extDict.has_key("irdm1") and extDict.has_key("ieri"):   # in case irdm1 == ieri
#       if extDict["irdm1"].name == extDict["ieri"].name:
#         print "if(irdm1 == ieri .and. srdm1 == seri) then"
#         numGuard += 1          
#     if extDict.has_key("irdm2") and extDict.has_key("ieri"):   # in case irdm2 == ieri    
#       if extDict["irdm2"].name == extDict["ieri"].name:
#         print "if(irdm2 == ieri .and. srdm2 == seri) then"
#         numGuard += 1
#     if extDict.has_key("irdm1") and extDict.has_key("irdm2"):  # in case irdm1 == irdm2
#       if extDict["irdm1"].name == extDict["irdm2"].name:
#         print "if(irdm1 == irdm2 .and. srdm1 == srdm2) then" 
#         numGuard += 1
#     # *OLD*

    # Declare the external indices ...
    if isBareLTS == True:
      print "   s_" + LTensor.indices[extBare].name + " = s" + lts_name # This can be different one
      print "   i_" + LTensor.indices[extBare].name + " = i" + lts_name # This can be different one

    countEri  = 0
    countBare = 0
    countRdm4 = 0
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == name_ofEri:  # in case of ERI
        print "   s_" + thisTensor.indices[extEri].name + " = seri"
        print "   i_" + thisTensor.indices[extEri].name + " = ieri"
        countEri += 1
      if thisTensor.name == name_ofBareAmp: # in case of BareAmp
        print "   s_" + thisTensor.indices[extBare].name + " = samp" # This can be different one
        print "   i_" + thisTensor.indices[extBare].name + " = iamp" # This can be different one
        countBare += 1 
      if thisTensor.name == name_ofRdm4:  # in case of RDM4
        print "   s_" + thisTensor.indices[extRdm4[0]].name + " = srdm1"
        print "   i_" + thisTensor.indices[extRdm4[0]].name + " = irdm1"
        print "   s_" + thisTensor.indices[extRdm4[1]].name + " = srdm2"
        print "   i_" + thisTensor.indices[extRdm4[1]].name + " = irdm2"
        countRdm4 += 1
      if countEri > 1 or countBare > 1 or countRdm4 > 1:
        raise TypeError, "Algorithmic error..."

    # Write down the loops over irreps ...
    for index in summed_Index:
      if index.isExt == True:
 #       print "**** " + index.name # *TEST* 
        continue      
      print "do s_" + index.name + " = 0, nir-1"
      loop_count += 1

    # Write down the constraints in terms of the irreps ...
    # Conditions for the left-hand tensor
    print "if( &"
    thisTensor = LTensor
    if len(thisTensor.indices) == 1:   # in case of unity
      print "S_" + thisTensor.indices[0].name + " == 0",
    elif len(thisTensor.indices) == 2: # in case of two
      print "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0",
    elif len(thisTensor.indices) == 3: # in case of three
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == 0"
      print cstr,
    elif len(thisTensor.indices) == 4: # in case of four
      cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      cstr += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
      print cstr,
    elif len(thisTensor.indices) == 5: # in case of five
      cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
      cstr += "s_" + thisTensor.indices[4].name + ")"
      print cstr,
    elif len(thisTensor.indices) == 6: # in case of six
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
      cstr += "s_" + thisTensor.indices[5].name + ")"
      print cstr,
    elif len(thisTensor.indices) == 7: # in case of seven
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
      cstr += "IEOR(s_" + thisTensor.indices[5].name + ", " + "s_" + thisTensor.indices[5].name + ")"
      print cstr,
    elif len(thisTensor.indices) == 8: # in case of eight (allowed only for 4-RDM)
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr +=      "IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + ")) == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[4].name + ", " + "s_" + thisTensor.indices[5].name + "),"
      cstr +=      "IEOR(s_" + thisTensor.indices[6].name + ", " + "s_" + thisTensor.indices[7].name + "))"
      print cstr,
    elif len(thisTensor.indices) != 0: # in case of scalar:
      raise TypeError, "Number of indices in tensor must be 0 - 8"
    if len(thisTensor.indices) != 0: print " .and. & " 
    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if len(thisTensor.indices) == 1:   # in case of unity
        print "S_" + thisTensor.indices[0].name + " == 0",
      elif len(thisTensor.indices) == 2: # in case of two
        print "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0",
      elif len(thisTensor.indices) == 3: # in case of three
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == 0"
        print cstr,
      elif len(thisTensor.indices) == 4: # in case of four
        cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        cstr += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
        print cstr,
      elif len(thisTensor.indices) == 5: # in case of five
        cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
        cstr += "s_" + thisTensor.indices[4].name + ")"
        print cstr,
      elif len(thisTensor.indices) == 6: # in case of six
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
        cstr += "s_" + thisTensor.indices[5].name + ")"
        print cstr,
      elif len(thisTensor.indices) == 7: # in case of seven
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
        cstr += "IEOR(s_" + thisTensor.indices[5].name + ", " + "s_" + thisTensor.indices[5].name + ")"
        print cstr,
      elif len(thisTensor.indices) == 8: # in case of eight (allowed only for 4-RDM)
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr +=      "IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + ")) == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[4].name + ", " + "s_" + thisTensor.indices[5].name + "),"
        cstr +=      "IEOR(s_" + thisTensor.indices[6].name + ", " + "s_" + thisTensor.indices[7].name + "))"
        print cstr,
      else:
        raise TypeError, "Number of indices in tensor must be 1 - 8"

      if not numTensor == len(thisTerm.tensors) - 1:
        print " .and. &"

    print ") then"
    
    # Write down the body of the tensorial contraction ...
    # *CAUTION* Only active and virtual orbitals are considered separately for now 
    count = 0
    for Index in summed_Index:
      if Index.isExt == True: 
#        print "**** " + Index.name # *TEST* 
        continue
      iname = "i_" + Index.name
      sname = "s_" + Index.name
      if   Index.indType[0][0] == 'core':
        otype = "I_C"
      elif Index.indType[0][0] == 'active':
        otype = "I_O"
      elif Index.indType[0][0] == 'virtual':
        otype = "I_V"
      print "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")"
      count += 1
    if count != loop_count:
      raise TypeError, "Algorithmic error occured .... "

    # As a first step of evalutation of the contraction, search the Kronecker deltas    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': num_kdeltas += 1

    if num_kdeltas >= 1:
      count = num_kdeltas
#      print "num_kdeltas ", num_kdeltas # *TEST*
      print "if(",
      for numTensor in range(len(thisTerm.tensors)):
        thisTensor = thisTerm.tensors[numTensor]
        if thisTensor.name == 'kdelta':
          cstr  = "s_" + thisTensor.indices[0].name + " == s_" + thisTensor.indices[1].name
          cstr += " .and. " + "i_" + thisTensor.indices[0].name + " == i_" + thisTensor.indices[1].name
          count -= 1
          if   count == 0: cstr += ") then"
          elif count >= 1: cstr += " .and. & "
          elif count >  0: raise TYpeError, "Algorithmic Error"
          print cstr

    print ""

    # First generate the left-hand side tensor
    # Make list of the indices which are of *NOT* external
    if not isBareLTS:
      LIndList = LTensor.indices
    else:                               # in case of BareAmP
      LIndList = []
      for n_index in range(len(LTensor.indices)):
        thisIndex = LTensor.indices[n_index]
        if n_index != extBare: LIndList.append(thisIndex)     
#     elif LTensor.name == name_ofEri:     # in case of ERI
#       LIndList = []
#       for n_index in range(len(LTensor.indices)):
#         thisIndex = LTensor.indices[n_index]
#         if thisIndex.isExt == False and n_index != extEri: LIndList.append(thisIndex)
#     elif LTensor.name == name_ofBareAmp: # in case of BareAmp
#       LIndList = []
#       for n_index in range(len(LTensor.indices)):
#         thisIndex = LTensor.indices[n_index]
#         if thisIndex.isExt == False and n_index != extBare: LIndList.append(thisIndex)
#     elif LTensor.name == name_ofRdm4:    # in case of 4-RDM
#       LIndList = []
#       for n_index in range(len(LTensor.indices)):
#         thisIndex = LTensor.indices[n_index]
#         # Dows it work correctly ??
#         if (thisIndex.isExt == False and n_index != extRdm4[0]) or (thisIndex.isExt == False and n_index != extRdm4[1]): 
#           LIndList.append(thisIndex)

    if len(LTensor.indices) != 0:   # in case of Tensor
      nterm = LTensor.name + "_("
      for n_index in range(len(LIndList)):
        nterm += "s_" + LIndList[-(n_index+1)].name
        if not n_index == len(LIndList)-1: nterm += ", "
        else: nterm += ")"
      nterm += "%array("
      for n_index in range(len(LIndList)):
        nterm += "i_" + LIndList[-(n_index+1)].name
        if not n_index == len(LIndList)-1: nterm += ", "
        else: nterm += ")"
      print nterm + " = " + nterm + " &"
    else: print LTensor.name + "_ = " + LTensor.name + "_ &"  # in case of Scalar  

    # Generate all the tensors ...
    if thisTerm.numConstant > 0.0: print "  + " + str(       thisTerm.numConstant) + " & "
    else:                          print "  - " + str((-1.0)*thisTerm.numConstant) + " & " 

    if len(thisTerm.constants) > 0: 
      for thisConst in thisTerm.constants: print "  * " + thisConst + " & "

    # Generate each tensor except Kronecker deltas ....
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': continue
     
      # Prepare the list of all indices, which are *NOT* of external, if the tensor is either of Eri or BareAmp
      if thisTensor.name != name_ofEri and thisTensor.name != name_ofBareAmp and thisTensor.name != name_ofRdm4:
        RIndList = thisTensor.indices
      elif thisTensor.name == name_ofEri:     # in case of ERI
        RIndList = []
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          if n_index != extEri: RIndList.append(thisIndex)
      elif thisTensor.name == name_ofBareAmp: # in case of BareAmp
        RIndList = []
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          if n_index != extBare: RIndList.append(thisIndex)
      elif thisTensor.name == name_ofRdm4:    # in case of 4-RDM
        RIndList = []
        dummy_count = 0
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          # Does it work correctly ??
          if (n_index != extRdm4[0]) and (n_index != extRdm4[1]): 
            RIndList.append(thisIndex)
            dummy_count += 1
        if(dummy_count > 6): raise TypeError, "4-RDMs must have at least 2 external indices" 

#       elif:
#         RIndTest = []
#         for index in thisTensor.indices: 
#           if index.isExt != True: RIndTest.append(index)

      nterm = "  * " + thisTensor.name + "_("
      for n_index in range(len(RIndList)):
        nterm += "s_" + RIndList[-(n_index+1)].name
        if not n_index == len(RIndList)-1: nterm += ", "
        else: nterm += ")"
      nterm += "%array("
      for n_index in range(len(RIndList)):
        nterm += "i_" + RIndList[-(n_index+1)].name
        if not n_index == len(RIndList)-1: nterm += ", "
        else: nterm += ")"
      if numTensor != len(thisTerm.tensors)-1: nterm += " & " 
      print nterm 
    
    print ""

    # Closing the iteration .... 
    if num_kdeltas >= 1: print "end if ! Kdeltas"                       # Kronecker deltas
    for i in range(loop_count): print "end do ! Irrep loop"             # Irrep loops 
    print "end if ! Irrep Cond"                                         # Irrep constraints
    for i in range(loop_count): print "end do ! Orbital loop"           # Orbital loops
    for i in range(numGuard) : print "end if ! Guard for external loop" # Guard for external loops
    print "! Loop_count     : ", loop_count
    print "! External_count : ", len(extIndices)
    print "! Total_count    : ", len(summed_Index)
    if loop_count + len(extIndices) != len(summed_Index):
      raise TypeError, "Algorithmic Error, Somethig is wrong !!!"
    print ""
    print ""
  print "! FEMTO END    **************************************************************" 


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertTeX(LTensor, inTerm, Term_per_row = 4):
  "Returns the TeX source of tenaorial equation given"

  # In case of that the LHS is a Scalar 
  if len(inTerm) == 0:
    print "! RHS is ZERO .... "
    return

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"

  nameRDMs = ["E2", "E3", "E4", "E5", "E6", "E7" ]

  print ""
  print "% Tensorial equation converted into TeX form ..... "
  print "%"
  print "% ", LTensor, "<-- "
  for thisTensor in inTerm:
    print "% ", thisTensor
  print "" 

  # Print name of Left tensor
  LTname = LTensor.name
  if len(LTensor.indices) != 0: LTname += "^{"
  for numIndex in range(len(LTensor.indices)):
    thisIndex = LTensor.indices[numIndex]
    LTname += thisIndex.name
    if   numIndex == len(LTensor.indices)/2 - 1: LTname += "}_{"
    elif numIndex == len(LTensor.indices)   - 1: LTname += "}"
    #else: LTname += ", "

  print LTname + " &+\!= ",

  # print all the terms on the RHS side ....
  isBroken = False
  for numTerm in range(len(inTerm)):
    thisTerm = inTerm[numTerm]
    thisLabel = ""
    if isBroken: thisLabel += "& "

    if   float(thisTerm.numConstant) > 0.0: thisLabel += "+ "
    else:                                   thisLabel += "- "
    if fabs(thisTerm.numConstant)%1.0 != 0 and fabs(thisTerm.numConstant) != 1.0:   
      thisLabel += str(fabs(thisTerm.numConstant)) 
    elif thisTerm.numConstant%1.0 == 0 and fabs(thisTerm.numConstant) != 1.0:      
      thisLabel += str(int(fabs(thisTerm.numConstant)))
    for thisConst in thisTerm.constants:
      thisLabel += thisConst     

    AllIndices = set([])
    # Relpace all the dummy indices ....
    for numTensor in range(len(thisTerm.tensors)-1): 
      thisTensor = thisTerm.tensors[numTensor  ]
      nextTensor = thisTerm.tensors[numTensor+1]

      set_i1 = set(thisTensor.indices)
      set_i2 = set(nextTensor.indices)
      set_i3 = set_i1.union(set_i2)
      AllIndices = set_i3.union(AllIndices)
    AllIndices = list(AllIndices)

    # *VERY* bithchy .....
    AllNames  = []
    AllSpaces = []
    for thisIndex in AllIndices:
      AllNames.append(thisIndex.name)
      AllSpaces.append(thisIndex.indType[0][0])

#    for ind in AllIndices: print ind.name, ind.indType[0][0], ind.isSummed # *TEST*

    C_count = 1
    O_Count = 1
    V_Count = 1
    for numIndex in range(len(AllIndices)):
      if   AllSpaces[numIndex] == 'core':    IndName = 'c_{' + str(C_Count) + '}'; C_Count += 1
      elif AllSpaces[numIndex] == 'active':  IndName = 'o_{' + str(O_Count) + '}'; O_Count += 1
      elif AllSpaces[numIndex] == 'virtual': IndName = 'v_{' + str(V_Count) + '}'; V_Count += 1
      else: raise TypeError, "Algorithmic Error"

      for thisTensor in thisTerm.tensors:
        for thisIndex in thisTensor.indices:
          if thisIndex.name == AllNames[numIndex] and thisIndex.indType[0][0] == AllSpaces[numIndex] and not thisIndex in LTensor.indices: 
            thisIndex.name = IndName
    
#     print ""                                                               # *TEST*    
#     for ind in AllIndices: print ind.name, ind.indType[0][0], ind.isSummed # *TEST*

    # Print each tensor ....
    for thisTensor in thisTerm.tensors:
      if thisTensor.name in nameRDMs:   thisLabel += "D^{"
      elif thisTensor.name == "kdelta": thisLabel += "\delta^{"
      else:                             thisLabel += thisTensor.name + "^{"
      for numIndex in range(len(thisTensor.indices)):
        thisIndex = thisTensor.indices[numIndex]
        thisLabel += thisIndex.name
        if   numIndex == len(thisTensor.indices)/2-1: thisLabel += "}_{"
        elif numIndex == len(thisTensor.indices)  -1: thisLabel += "} "

    if (numTerm+1)%Term_per_row == 0: print thisLabel + "\\\\"; isBroken = True
    else:                             print thisLabel,      ; isBroken = False     
  
  if numTerm != len(inTerm)-1: raise TypeError, "Algorithmic Error"


  
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def factorize(LTensor, inTerm, title, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2', priority = 'mem'):
  "Returns a stream of the tensorial formulae, which composes the overall tensorial equation such as the sigma equation."
  "Order of the memory requirement to keep the intermediates and polynomial scaling are optimized."
  "In this method, LHS is assumed to be a BareAmpPack quantity."
  "*WARNING* optimization condition is sometimes meaningless unless algorithm used in this code is not generalized."

  # In case of that the LHS is a Scalar 
  if len(inTerm) == 0:
    print "! RHS is ZERO .... "
    return
  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a string"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORz code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # check whether the LTensor could be a BareAmp if specified so ....
  if len(LTensor.indices) != 4 and isBareLTS:
    raise TypeError, "LTensor could not be a BareAmpPack .... "

  # check whether priority is string ....
  if not isinstance(priority, type('a')):
    raise TypeError, "Priority must given as a string"

  # Output for declarations for Fortaran codes ....
  CHname = "c_" + title + ".h"
  CHfile = open(CHname, 'w')

  CHguard = "c_" + title + "_h"
  CHguard = CHguard.upper()
  CHbegin = \
"""
// #include <tensor/tensor.h>
// #include <sci/hint/hintmo/hintmo.h>
// #include <sci/ctnew2/ctclass_input.h>
// #include <sci/ctnew2/ctclass_symblock.h>
// #include <sci/ctnew2/ctclass_rdmpack.h>
// #include <sci/ctnew2/ctclass_bareamppack.h>
 
extern "C"{

"""

  CHend = \
"""

} 


#endif

"""

  CHfile.write( "\n" )
  CHfile.write( "#ifndef " + CHguard + "\n" )
  CHfile.write( "#define " + CHguard + "\n" )
  CHfile.write( "\n" )
  CHfile.write( CHbegin )

  # Output for tensorial contractions ....
  Fname = "f_" + title + ".F90"
  F90file = open(Fname, 'w')
  F90file.write('#include \"../f_ct.fh\"\n\n')

  # Output for rest of tensorial contractions ....
  BareFlag = False
  for thisTerm in inTerm:
    for t in thisTerm.tensors:    
      if t.name == name_ofBareAmp: BareFlag = True; break

  CPfile = "c_" + title + ".cpp"
  CPfile = open(CPfile, "w")
  CPfile.write( \
"""

#include <orz/orz.h>
#include <orz/openmp.h>
#include <orz/cblas.h>
#include <orz/clapack.h>
#include <tensor/tensor.h>
#include <sci/hint/para_disttools.h>
#include <sci/ctnew2/ct.h>
#include <sci/ctnew2/ct_f.h>
#include <sci/ctnew2/ctclass_input.h>
#include <sci/ctnew2/ctclass_symblock.h>
#include <sci/ctnew2/ctclass_hintmo.h>
#include <sci/ctnew2/ctclass_rdmpack.h>
#include <sci/ctnew2/ctclass_bareamppack.h>
#include <sci/ctnew2/ctclass_orthamppack.h>
#include <sci/ctnew2/diaghessian.h>
#include <sci/ctnew2/symamp2.h>
#include <sci/ctnew2/mrci.h>

using std::cout;
using std::endl;

#define FLOPCOUNT


// ***************************************************************************
// orz::ct::mrci
// ***************************************************************************
									    /*!
   @brief CT input
									     */

""")

  if isBareLTS and not BareFlag:
    CPfile.write(\
"""
orz::ct::BareAmpPack orz::ct::%s(const orz::ct::Input &ctinp,
					 const orz::ct::SymBlockInfo &symblockinfo,
					 const orz::ct::HintMO &hintmo,
					 const orz::ct::RdmPack &rdmPack_sym,
					 const orz::DTensor &rdm4,
					 const double zeroth_energy) {

  // set up nmo nclosed, nocc
  const FC_INT nclosed = ctinp.nclosed();
  const FC_INT nocc    = ctinp.nocc();
  const FC_INT nvir    = ctinp.nvir();
  const FC_INT nmo     = nclosed + nocc + nvir;

  const FC_INT nir = symblockinfo.nir();
  const FC_INT * const nsym = symblockinfo.nsym().cptr();
  const FC_INT * const psym = symblockinfo.psym().cptr();
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();

  orz::DTensor moint1 = hintmo.int1(); // Setting up the one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  orz::ct::BareAmpPack retval;
  retval = orz::ct::BareAmpPack(ctinp, symblockinfo, \"%s\");
  orz::DTensor %s(nmo,nmo,nmo);
  orz::DTensor rdm4_sym;
  double * const %s_ptr = %s.cptr();

""" % (title, LTensor.name, name_ofEri, name_ofEri, name_ofEri))

  elif isBareLTS and BareFlag:
    CPfile.write(\
"""
orz::ct::BareAmpPack orz::ct::%s(const orz::ct::Input &ctinp,
					 const orz::ct::SymBlockInfo &symblockinfo,
					 const orz::ct::HintMO &hintmo,
					 const orz::ct::RdmPack &rdmPack_sym,
					 const orz::DTensor &rdm4,
					 const orz::ct::BareAmpPack &%s,
					 const int num_sigma,
					 const double T0,
					 const double zeroth_energy) {


  // set up nmo nclosed, nocc
  const FC_INT nclosed = ctinp.nclosed();
  const FC_INT nocc    = ctinp.nocc();
  const FC_INT nvir    = ctinp.nvir();
  const FC_INT nmo     = nclosed + nocc + nvir;
  const FC_INT nir     = symblockinfo.nir();
  const FC_INT * const nsym    = symblockinfo.nsym().cptr();
  const FC_INT * const psym    = symblockinfo.psym().cptr();
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();

  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  std::ostringstream stm;
  stm << num_sigma;
  std::string name_of_sigma = \"%s[\" + stm.str() + \"]\"; // Name of the Sigma vector
  orz::ct::BareAmpPack retval
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma); // Sigma(a, a', e, e') tensor

  orz::DTensor %sb; // Container of S2_aae,[b] tensor
  orz::DTensor %sb; // Container of T2_aae,[b] tensor

  orz::DTensor %s(nmo,nmo,nmo);
  orz::DTensor rdm4_sym;
  double * const %s_ptr = %s.cptr();

""" % (title, name_ofBareAmp, LTensor.name, LTensor.name, name_ofBareAmp, name_ofEri, name_ofEri, name_ofEri))

  elif not len(LTensor.indices) and BareFlag:
    CPfile.write(\
"""
double orz::ct::%s(const orz::ct::Input &ctinp,
			     const orz::ct::SymBlockInfo &symblockinfo,
			     const orz::ct::HintMO &hintmo,
			     const orz::ct::RdmPack &rdmPack_sym,
			     const double init_value,
			     const orz::ct::BareAmpPack &%s,
			     const double zeroth_energy) {

  // set up nmo nclosed, nocc
  const FC_INT nclosed = ctinp.nclosed();
  const FC_INT nocc    = ctinp.nocc();
  const FC_INT nvir    = ctinp.nvir();
  const FC_INT nmo     = nclosed + nocc + nvir;

  const FC_INT nir = symblockinfo.nir();
  const FC_INT * const nsym = symblockinfo.nsym().cptr();
  const FC_INT * const psym = symblockinfo.psym().cptr();
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();

  orz::DTensor moint1 = hintmo.int1(); // Setitng up one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  double %s = init_value;
  orz::DTensor %s(nmo,nmo,nmo); // V(p,nmo,nmo,nmo) for an orbital index p
  orz::DTensor %sb;             // Container of T2_ija,[b] tensor

  double * const %s_ptr = %s.cptr();
""" % (title, name_ofBareAmp, LTensor.name, name_ofEri, name_ofBareAmp, name_ofEri, name_ofEri))

  else:
    raise RunTimeError, " I don't know how to handle this case ....  "

  print "!"
  print "!  Restate my assumptions:                                                   " 
  print "!  1, Mathematics is the language of nature.                                 " 
  print "!  2, Everything around us can be represented and understood through numbers." 
  print "!  3, If you graph the numbers of any system, patterns emerge.               " 
  print "!  Therefore: There are patterns everywhere in nature.                       " 
  print "!                                      ~~ Maximillian Cohen in Pi (film) ~~  " 
  print "!"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  LegTensors = [name_ofEri, name_ofRdm4, name_ofBareAmp]

  # Size of the moel system ....
  numCore = 0    
  numOcc  = 30
  numVirt = 500 

  # The maximum size of intermediate allowed ....
  MaxSize = numOcc ** 6 # 3-RDM 

  # If the LHS is a BareAmp, set so ...
  if isBareLTS == True: 
    LTensor.indices[extBare].isExt = True

  # whether X is processed in on-the-fly fashion if V and S(or Hdiag) have no index in common
  isOF = True

  # Pointer which stands for the end point (so don't change it ..... )
  SigmaCount = False
  D4Count    = False

  # Whether the indices in Sigma is rotated or not (so don't change it ..... )
  isRotated = False

  # Whether a kdelta is killed (and as a consequences, indices in LTensor are altered) or not
  isKDkilled = False

  print ""
  # Search in each term ...
  num_loops    = []
  for numTerm in range(len(inTerm)):
    thisTerm = inTerm[numTerm]

    # As a preliminary step, kill the kdeltas if they have same indices
    kill_k = False
    for num_t1 in range(len(thisTerm.tensors)):
      t1 = thisTerm.tensors[num_t1]
      if t1.name == 'kdelta':
        for num_t2 in range(num_t1+1, len(thisTerm.tensors)):
          t2 = thisTerm.tensors[num_t2]
          if t2.name == 'kdelta':
            temp_i1 = t1.indices[:]
            temp_i2 = t2.indices[:]
            temp_i1.sort()
            temp_i2.sort()
            if temp_i1 == temp_i2: 
              thisTerm.tensors.remove(t2)
#              print "!!! WARNING: one kdelta is killed !!!"   
              kill_k = True
              break
        if kill_k: break

    # Kill kdelta if this is composed of purely dummy indices
    for t in thisTerm.tensors:
      if t.name == 'kdelta':
        if (not t.indices[0].isSummed) and (not t.indices[1].isSummed):
#          print "!!! WARNING: Indices in LHS is changed by kdelta !!!"
          isKDkilled = True
          LTcopy = LTensor.copy()
          tempInd = t.indices[:]
          tempInd.sort()
          for i in LTensor.indices:
            if i == tempInd[1]: 
              ind_i = LTensor.indices.index(i)
              LTensor.indices[ind_i] = tempInd[0].copy()
          for t2 in thisTerm.tensors:
            if t2.name == 'kdelta': continue
            else:
              for i in t2.indices:
                if i == tempInd[1]:
                  ind_i = t2.indices.index(i)
                  t2.indices[ind_i] = tempInd[0].copy()
          if isBareLTS:
            LTensor.indices[extBare].isExt = True
            temp_i = LTensor.indices[extBare].copy()
            for i in LTensor.indices:
              if i.name == temp_i.name: i.isExt = True
          thisTerm.tensors.remove(t)
              
        else: 
          raise TypeError, "I cannot handle Kronecker delta with purely dummy indices"

    set_V  = False; num_V  = False
    set_T  = False; num_T  = False
    set_D4 = False; num_D4 = False 
    if isBareLTS:
      # Check the BUG !!!!!!!!!!!!
      # First of all, rotate the leg index in ERI to be the same to that of the sigma
      thisName = LTensor.indices[extBare].name
      # Consistency between S(or Hdiag) and V, which is prior to all the other ones
      for t in thisTerm.tensors:
        if t.name == name_ofEri and LTensor.indices[extBare] in t.indices: 
          pLeg = t.indices.index(LTensor.indices[extBare])
          if   pLeg == 0:
            newIndices = [t.indices[0], t.indices[1], t.indices[2], t.indices[3]]
          elif pLeg == 1:
            newIndices = [t.indices[1], t.indices[0], t.indices[2], t.indices[3]]
          elif pLeg == 2:
            newIndices = [t.indices[2], t.indices[3], t.indices[0], t.indices[1]]
          elif pLeg == 3:
            newIndices = [t.indices[3], t.indices[2], t.indices[1], t.indices[0]]
          t.indices = newIndices[:]
          set_V = True
          num_V = thisTerm.tensors.index(t)
      # Consistency between S(or Hdiag) and T
      for t in thisTerm.tensors:
        if t.name == name_ofBareAmp and LTensor.indices[extBare] in [t.indices[2],t.indices[3]] \
          and LTensor.indices[2].indType[0][0] == LTensor.indices[3].indType[0][0]:
          if   t.indices.index(LTensor.indices[extBare]) == 3: set_T = True; num_T = thisTerm.tensors.index(t)
          elif t.indices.index(LTensor.indices[extBare]) == 2 and t.indices[2].indType[0][0] == t.indices[3].indType[0][0] and not set_V:
            newIndices = [t.indices[1], t.indices[0], t.indices[3], t.indices[2]]
            t.indices = newIndices[:]
            set_T = True

    # *VERY* bitchy and suspicious ......
    # Consisitency between V and D4
    # MS: I think this works (2012/06/14)
    Ocount = False
    is_V2  = False
    is_D4  = False
    if (not set_V):
      Ocount = 0
      for t in thisTerm.tensors:
        if   t.name == name_ofEri:
          is_V2 = True
          num_V  = thisTerm.tensors.index(t)
          for i in t.indices:
            if not i.indType[0][0] == "virtual": Ocount += 1
        elif t.name == name_ofRdm4: 
          num_D4 = thisTerm.tensors.index(t)
          is_D4  = True

    if is_V2 and is_D4 and Ocount == 4:
      for i1 in thisTerm.tensors[num_V].indices:
        for num_i2 in range(len(thisTerm.tensors[num_D4].indices)):
          if num_i2 % 2: continue
          list_i2 = [thisTerm.tensors[num_D4].indices[num_i2], thisTerm.tensors[num_D4].indices[num_i2+1]]
          list_other = thisTerm.tensors[num_D4].indices[:]
          list_other.remove(list_i2[0])
          list_other.remove(list_i2[1])
          if i1 in list_i2:
            thisTerm.tensors[num_D4].indices = list_i2 + list_other
            V2 = thisTerm.tensors[num_V]
            pLeg = V2.indices.index(i1)
            if   pLeg == 0:
              newIndices = [V2.indices[0], V2.indices[1], V2.indices[2], V2.indices[3]]
            elif pLeg == 1:
              newIndices = [V2.indices[1], V2.indices[0], V2.indices[2], V2.indices[3]]
            elif pLeg == 2:
              newIndices = [V2.indices[2], V2.indices[3], V2.indices[0], V2.indices[1]]
            elif pLeg == 3:
              newIndices = [V2.indices[3], V2.indices[2], V2.indices[1], V2.indices[0]]
            V2.indices = newIndices[:]
            set_D4 = True

          if set_D4: break    
        if set_D4: break          

    print "! No.", numTerm
    print "!", LTensor, " <-- "
    print "!", thisTerm

    # Dictionary that shows which index corresponds to ieri, isig, and so on ...
    extDict = {}
    # List of all external indices without any redundancy
    extIndices = []
    # If the LHS is a BareAmp, set so ...
    if isBareLTS == True: 
      extIndices.append(LTensor.indices[extBare]) # Record which indices are of external
      extDict["i"+lts_name] = LTensor.indices[extBare]

    for thisTensor in thisTerm.tensors:
      # If thisTensor is ERI, or BareAmp, the external indices have to be dropped off ...
      if   thisTensor.name == name_ofEri:     
        thisTensor.indices[extEri].isExt  = True
        extIndices.append(thisTensor.indices[extEri])
        extDict["ieri"] = thisTensor.indices[extEri]
      elif thisTensor.name == name_ofBareAmp: 
        thisTensor.indices[extBare].isExt = True
        extIndices.append(thisTensor.indices[extBare])
        extDict["iamp"] = thisTensor.indices[extBare] # This can be different one
      elif thisTensor.name == name_ofRdm4:
        thisTensor.indices[extRdm4[0]].isExt = True
        extIndices.append(thisTensor.indices[extRdm4[0]])
        extDict["irdm1"] = thisTensor.indices[extRdm4[0]]
        thisTensor.indices[extRdm4[1]].isExt = True
        extIndices.append(thisTensor.indices[extRdm4[1]])
        extDict["irdm2"] = thisTensor.indices[extRdm4[1]]

    # Avoid double counting     
    extIndices = list(set(extIndices)) 
    # Set all the external indices commonly beared by different tensors ....
    for extIndex in extIndices:
      for thisTensor in thisTerm.tensors:
        for thisIndex in thisTensor.indices:
          if extIndex.name == thisIndex.name:
            thisIndex.isExt = True

    # Search inbetween each tensor ...
    allIndices = LTensor.indices
    for numTensor in range(len(thisTerm.tensors)-1):
      thisTensor = thisTerm.tensors[numTensor]
      nextTensor = thisTerm.tensors[numTensor+1]

      set_i1     = set(thisTensor.indices)
      set_i2     = set(nextTensor.indices)
      set_i3     = set_i1.union(set_i2)
      allIndices = set_i3.union(allIndices)

    # Factorize into intermediates
    numCase = 0
    caseTensor = []
    for numTensor1 in range(len(thisTerm.tensors)):
      Tensor1 = thisTerm.tensors[numTensor1]
      for numTensor2 in range(numTensor1+1, len(thisTerm.tensors)):
        Tensor2 = thisTerm.tensors[numTensor2]
        # ERI should be the first. Due to the length of ERI tensor
        if   Tensor1.name == name_ofEri:
          Pattern = term(1.0, [], [Tensor1, Tensor2])
        elif Tensor2.name == name_ofEri:
          Pattern = term(1.0, [], [Tensor2, Tensor1])
        else:
          Pattern = term(1.0, [], [Tensor1, Tensor2])
#        Pattern = term(thisTerm.numConstant, thisTerm.constants, [Tensor1, Tensor2])

        summed_Index = LTensor.indices
        set_i1       = set(Tensor1.indices)
        set_i2       = set(Tensor2.indices)
        set_i3       = set_i1.union(set_i2)
        summed_Index = set_i3.union(set_i3)

        Interm_Index = []
        for i in summed_Index:
          if (i in Tensor1.indices) and (i in Tensor2.indices) and i.isSummed:
            Interm_Index.append(i)
        Interm_Index = summed_Index.difference(Interm_Index)

#         print "***** 1 IntermIndex *****" # *TEST*
#         for i in Interm_Index:            # *TEST*
#           print i.name, i.isExt           # *TEST*
#         print                             # *TEST*   

        # Reorder Interm_Index as [O, O, .... ,V, V, ....]
        Interm_Index = list(Interm_Index)
        Occ_Index = []
        Vir_Index = []
        for i in Interm_Index:
          if   i.indType[0][0] == 'core' or i.indType[0][0] == 'active':
            Occ_Index.append(i)
          elif i.indType[0][0] == 'virtual':
            Vir_Index.append(i)
        Interm_Index = Occ_Index + Vir_Index

#         print "***** 2 IntermIndex *****" # *TEST*
#         for i in Interm_Index:            # *TEST*
#           print i.name, i.isExt           # *TEST*
#         print                             # *TEST*   
        
        Interm = tensor('X', Interm_Index[:])        
        print "Case " + str(numCase) + " ..... ", Interm, " <----- ", Pattern

#         print "***** 3 IntermIndex *****" # *TEST*
#         for i in Interm.indices:          # *TEST*
#           print i.name, i.isExt           # *TEST*
#         print                             # *TEST*   

        # Estimate the polynomial order of the first contraction ...
        nCore = 0
        nOcc  = 0
        nVir  = 0
        for i in summed_Index:
          if   i.indType[0][0] == 'core':    nCore += 1
          elif i.indType[0][0] == 'active':  nOcc  += 1
          elif i.indType[0][0] == 'virtual': nVir  += 1

        # Estimate the polynomial order of the second contraction ...
        nCtemp = 0
        nOtemp = 0
        nVtemp = 0
        leftIndex = set(allIndices).difference(set(summed_Index))
        leftOrder = leftIndex.union(set(Interm_Index))
        for i in leftOrder:
          if   i.indType[0][0] == 'core':    nCtemp += 1
          elif i.indType[0][0] == 'active':  nOtemp += 1
          elif i.indType[0][0] == 'virtual': nVtemp += 1

#        if (log(numCore)*nCtemp+log(numOcc)*nOtemp+log(numVirt)*nVtemp > log(numCore)*nCore+log(numOcc)*nOcc+log(numVirt)*nVir):
        if (log(numOcc)*nOtemp+log(numVirt)*nVtemp > log(numOcc)*nOcc+log(numVirt)*nVir):
          nCore = nCtemp
          nOcc  = nOtemp
          nVir  = nVtemp
        print "! Polynomial order   is O(c^" + str(nCore) + "o^" + str(nOcc) + "v^" + str(nVir) + ")"

        thisCase = {}
        # The heaviest one shoul be stored as polynomial order ....
        thisCase['Poly_C']       = nCore         
        thisCase['Poly_O']       = nOcc
        thisCase['Poly_V']       = nVir

        # Estimate the memory requirement of the first contraction
        nCore = 0
        nOcc  = 0
        nVir  = 0
        for i in Interm_Index:
          if   i.indType[0][0] == 'core':    nCore += 1
          elif i.indType[0][0] == 'active':  nOcc  += 1
          elif i.indType[0][0] == 'virtual': nVir  += 1
        print "! Memory requirement is O(c^" + str(nCore) + "o^" + str(nOcc) + "v^" + str(nVir) + ")"
        print ""

        # Store the information on the pattern of formation of the intermediate .... 
        thisCase['Intermediate'] = Interm
        thisCase['Tensors']      = Pattern
        thisCase['Int_C']       = nCore         
        thisCase['Int_O']       = nOcc
        thisCase['Int_V']       = nVir
        caseTensor.append(thisCase)
        numCase += 1

      # Select which one is the best way ....
      # MEM is prior to INT
      memReq_C  = 1E+70
      memReq_O  = 1E+70
      memReq_V  = 1E+70
      intReq_C  = 1E+70
      intReq_O  = 1E+70
      intReq_V  = 1E+70
      mem_prior = 0
      int_prior = 0

      init_VFlag =True
      for t in thisTerm.tensors:
        if t.name == name_ofEri: init_VFlag = False

      for theCase in range(numCase):
        VFlag = init_VFlag
        for t in caseTensor[theCase]['Tensors'].tensors:
          if t.name == name_ofEri: VFlag = True
        if not VFlag: continue
#        if ((memReq_C*log(numCore))+(memReq_O*log(numOcc))+(memReq_V*log(numVirt))) > ((caseTensor[theCase]['Int_C']*log(numCore))+(caseTensor[theCase]['Int_O']*log(numOcc))+(caseTensor[theCase]['Int_V']*log(numVirt))):
        if ((memReq_O*log(numOcc))+(memReq_V*log(numVirt))) > ((caseTensor[theCase]['Int_O']*log(numOcc))+(caseTensor[theCase]['Int_V']*log(numVirt))):
          memReq_C = caseTensor[theCase]['Int_C']
          memReq_O = caseTensor[theCase]['Int_O']
          memReq_V = caseTensor[theCase]['Int_V']
          mem_prior = theCase
#        if ((intReq_C*log(numCore))+(intReq_O*log(numOcc))+(intReq_V*log(numVirt))) > ((caseTensor[theCase]['Poly_C']*log(numCore))+(caseTensor[theCase]['Poly_O']*log(numOcc))+(caseTensor[theCase]['Poly_V']*log(numVirt))):
        if ((intReq_O*log(numOcc))+(intReq_V*log(numVirt))) > ((caseTensor[theCase]['Poly_O']*log(numOcc))+(caseTensor[theCase]['Poly_V']*log(numVirt))):
          intReq_C = caseTensor[theCase]['Poly_C']
          intReq_O = caseTensor[theCase]['Poly_O']
          intReq_V = caseTensor[theCase]['Poly_V']
          int_prior = theCase

    if mem_prior != int_prior: print "\n! Conflict between optimal memory requirement and polynomial order ...."

    if   priority == 'mem':
      opt_prior = int_prior # This can be switched ..... 
    elif priority == 'int':
      opt_prior = mem_prior # This can be switched ..... 
    else:
      raise TypeError, "priority should be given either mem or int"

# *NO_GOOD*     # Check the bug!!!
# *NO_GOOD*     # If the leg of ERI cannot be matched with that of Sigma, rotate leg of Sigma if possible ....
# *NO_GOOD*     if isBareLTS and (not set_V) and LTensor.indices[2].indType[0][0] == LTensor.indices[3].indType[0][0]:
# *NO_GOOD*       thisName = LTensor.indices[extBare-1].name
# *NO_GOOD*       for t in thisTerm.tensors:
# *NO_GOOD*         if t.name == name_ofEri and LTensor.indices[extBare-1] in t.indices: 
# *NO_GOOD*           pLeg = t.indices.index(LTensor.indices[extBare-1])
# *NO_GOOD*           t.indices[extEri].isExt = False
# *NO_GOOD*           if   pLeg == 0:
# *NO_GOOD*             newIndices = [t.indices[0], t.indices[1], t.indices[2], t.indices[3]]
# *NO_GOOD*           elif pLeg == 1:
# *NO_GOOD*             newIndices = [t.indices[1], t.indices[0], t.indices[2], t.indices[3]]
# *NO_GOOD*           elif pLeg == 2:
# *NO_GOOD*             newIndices = [t.indices[2], t.indices[3], t.indices[0], t.indices[1]]
# *NO_GOOD*           elif pLeg == 3:
# *NO_GOOD*             newIndices = [t.indices[3], t.indices[2], t.indices[1], t.indices[0]]
# *NO_GOOD*           t.indices = newIndices[:]
# *NO_GOOD*           t.indices[extEri].isExt = True
# *NO_GOOD*           set_V = True
# *NO_GOOD*           LTensor.indices[extBare].isExt = False
# *NO_GOOD*           newIndices = [LTensor.indices[1], LTensor.indices[0], LTensor.indices[3], LTensor.indices[2]]
# *NO_GOOD*           LTensor.indices = newIndices[:]
# *NO_GOOD*           LTensor.indices[extBare].isExt = True
# *NO_GOOD*           isRotated = True
# *NO_GOOD* 
# *NO_GOOD*           for t2 in thisTerm.tensors + caseTensor[opt_prior]['Tensors'].tensors:
# *NO_GOOD*             if   t2.name == name_ofEri     and LTensor.indices[extBare-1] == t2.indices[extEri]: 
# *NO_GOOD*               LTensor.indices[extBare-1].isExt = True
# *NO_GOOD*             elif t2.name == name_ofBareAmp and LTensor.indices[extBare-1] == t2.indices[extBare]:
# *NO_GOOD*               LTensor.indices[extBare-1].isExt = True
# *NO_GOOD* 
# *NO_GOOD*       for t in caseTensor[opt_prior]['Tensors'].tensors:
# *NO_GOOD*         if t.name == name_ofEri and LTensor.indices[extBare-1] in t.indices: 
# *NO_GOOD*           pLeg = t.indices.index(LTensor.indices[extBare-1])
# *NO_GOOD*           t.indices[extEri].isExt = False
# *NO_GOOD*           if   pLeg == 0:
# *NO_GOOD*             newIndices = [t.indices[0], t.indices[1], t.indices[2], t.indices[3]]
# *NO_GOOD*           elif pLeg == 1:
# *NO_GOOD*             newIndices = [t.indices[1], t.indices[0], t.indices[2], t.indices[3]]
# *NO_GOOD*           elif pLeg == 2:
# *NO_GOOD*             newIndices = [t.indices[2], t.indices[3], t.indices[0], t.indices[1]]
# *NO_GOOD*           elif pLeg == 3:
# *NO_GOOD*             newIndices = [t.indices[3], t.indices[2], t.indices[1], t.indices[0]]
# *NO_GOOD*           t.indices = newIndices[:]
# *NO_GOOD*           t.indices[extEri].isExt = True
# *NO_GOOD*           LTensor.indices[extBare].isExt = False 
# *NO_GOOD*           set_V = True
# *NO_GOOD*           newIndices = [LTensor.indices[1], LTensor.indices[0], LTensor.indices[3], LTensor.indices[2]]
# *NO_GOOD*           LTensor.indices = newIndices[:]
# *NO_GOOD*           LTensor.indices[extBare].isExt = True
# *NO_GOOD*           isRotated = True
# *NO_GOOD* 
# *NO_GOOD*           for t2 in thisTerm.tensors + caseTensor[opt_prior]['Tensors'].tensors:
# *NO_GOOD*             if   t2.name == name_ofEri     and LTensor.indices[extBare-1] == t2.indices[extEri]: 
# *NO_GOOD*               LTensor.indices[extBare-1].isExt = True
# *NO_GOOD*             elif t2.name == name_ofBareAmp and LTensor.indices[extBare-1] == t2.indices[extBare]:
# *NO_GOOD*               LTensor.indices[extBare-1].isExt = True
        
    print "! ----------------------------------"
    print "! set_V       : ", set_V      # *TEST*
    print "! set_T       : ", set_T      # *TEST*
    print "! set_D4      : ", set_D4     # *TEST*
    print "! isRotated   : ", isRotated  # *TEST*
    print "! isKD_killed : ", isKDkilled # *TEST* 
    print "! ----------------------------------"

    CPfile.write("\n\n  {\n  // No.%s\n" % numTerm)

    if set_V and set_D4:
      raise RunTimeError, "Algorithmic Error"

    print ""
    print "! The optimal choice is .... "
    if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices):
      caseTensor[opt_prior]['Intermediate'] = LTensor
      caseTensor[opt_prior]['Tensors']      = thisTerm
    print caseTensor[opt_prior]['Intermediate'], " <-- ", caseTensor[opt_prior]['Tensors'],
    print ""

    # Substitute the intermediate into the tensorial equation ....
    if not set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices):
      for Tensor1 in caseTensor[opt_prior]['Tensors'].tensors:
        thisTerm.tensors.remove(Tensor1)
      thisTerm.tensors.append(caseTensor[opt_prior]['Intermediate'].copy())
      print LTensor, " <-- ", thisTerm
    print ""

#     print "***** 4 IntermIndex *****" # *TEST*
#     for t in thisTerm.tensors:
#       print t
#       for i in t.indices:
#         print i.name, i.isExt 

    # Print all the information on the factorization ....
    print "! Scaling       : O(c^" + str(caseTensor[opt_prior]['Poly_C']) + " o^" + str(caseTensor[opt_prior]['Poly_O']) + " v^" + str(caseTensor[opt_prior]['Poly_V']) + ")"
    print "! Max size of X : c^"   + str(caseTensor[opt_prior]['Int_C'])  + " o^" + str(caseTensor[opt_prior]['Int_O'])  + " v^" + str(caseTensor[opt_prior]['Int_V'])

    Indent = ""
    print "\n! * Begin scaling analysis .... *\n"

    # Starts the first contraction .... 
    LoopCount = 0
    Indent = []
    Declared = []
    Cs     = "" 
    for i in range(5000): # this number can be varied
      Indent.append(Cs)
      Cs += "  "

    # ********* First Contraction ************
    title1if = "g_if_" + title + "_no0_x" + str(numTerm) 
    title1   = "g_"    + title + "_no0_x" + str(numTerm)

    if not isOF and not set_V and num_V:
      for i in caseTensor[opt_prior]['Intermediate'].indices:
        i.isExt = False
    
    Set_i  = []
    ExtInd = []
    NameT  = []
    Consts = caseTensor[opt_prior]['Tensors'].constants
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      if (not 'E' in t.name) and (not t.name == 'kdelta'): NameT.append(t.name)
      Set_i.extend(t.indices[:])
    for i in Set_i:
      if i.isExt == True: ExtInd.append(i.copy())
    ExtInd = list(set(ExtInd))
    NameT.append(caseTensor[opt_prior]['Intermediate'].name)

    Set_i = list(Set_i)
    DecFirst = {"DecInd":ExtInd, "DecConst":Consts, "DecTensor":NameT}    
 
    if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices):
      caseTensor[opt_prior]['Tensors'].constants   = thisTerm.constants
      caseTensor[opt_prior]['Tensors'].numConstant = thisTerm.numConstant
      # Declarations for Fortran code included from C++ 
      convertCpp_header(LTensor, [caseTensor[opt_prior]['Tensors']], title1if, CHfile, DecFirst, isBareLTS, \
        lts_name, name_ofEri, name_ofBareAmp) # *TEST*      
      # Interface inbetween C++ and Fortran parts
      convertF90_interface(LTensor, [caseTensor[opt_prior]['Tensors']], title1if, F90file, DecFirst, isBareLTS, \
        lts_name, name_ofEri, name_ofBareAmp) # *TEST*
      # Body of the contraction
      convertF90_ERI2(LTensor, [caseTensor[opt_prior]['Tensors']], title1, F90file, DecFirst, isBareLTS, \
        lts_name, name_ofEri, name_ofBareAmp) # *TEST*
    else:
      # Declarations for Fortran code included from C++ 
      convertCpp_header(caseTensor[opt_prior]['Intermediate'], [caseTensor[opt_prior]['Tensors']], title1if, \
        CHfile, DecFirst, False, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
      # Interface inbetween C++ and Fortran parts
      convertF90_interface(caseTensor[opt_prior]['Intermediate'], [caseTensor[opt_prior]['Tensors']], title1if, \
        F90file, DecFirst, False, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
      # Body of the contraction
      convertF90_ERI2(caseTensor[opt_prior]['Intermediate'], [caseTensor[opt_prior]['Tensors']], title1, \
        F90file, DecFirst, False, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
#    print ""

    O1 = []
    for i in caseTensor[opt_prior]['Intermediate'].indices:
      if   i.isExt == True: O1.append(i) 
    O1 = list(set(O1))

    # Newly added
    # Check the bug!!
    if not isOF and not set_V and num_V:
      O1 = []
#       for t in caseTensor[opt_prior]['Tensors'].tensors:
#         if t.name == name_ofEri: 
#           num_V = caseTensor[opt_prior]['Tensors'].tensors.index(t)
#       V2 = caseTensor[opt_prior]['Tensors'].tensors[num_V]
# #      print V2
# #      print num_V
#       for Ind in LTensor.indices:
# #        print type(Ind)
# #        print type(V2.indices)
#         if Ind in V2.indices:
#           O1.append(Ind)
#           pLeg = V2.indices.index(Ind)
# #          print pLeg
#           V2.indices[extEri].isExt = False
#           if   pLeg == 0:
#             newIndices = [V2.indices[0], V2.indices[1], V2.indices[2], V2.indices[3]]
#           elif pLeg == 1:
#             newIndices = [V2.indices[1], V2.indices[0], V2.indices[2], V2.indices[3]]
#           elif pLeg == 2:
#             newIndices = [V2.indices[2], V2.indices[3], V2.indices[0], V2.indices[1]]
#           elif pLeg == 3:
#             newIndices = [V2.indices[3], V2.indices[2], V2.indices[1], V2.indices[0]]
#           V2.indices = newIndices[:]
#           V2.indices[extEri].isExt = True
#           MaximizeEri =True
# #          print V2
#         if MaximizeEri: break
      print "! Intermediate is not processed in ad hoc fashion .... "
      print ""
#      print chr(7), 

    # Order of the leg index of LHS has to be top
    if isBareLTS:
      if len(O1) > 0 and LTensor.indices[extBare] in O1:
        O1.remove(LTensor.indices[extBare])
        O1.insert(0, LTensor.indices[extBare])  
  
    # Check the bug!!!!!
    # if O1 has both irdm1 and irdm2, irdm1 must be prior to irdm2
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      if t.name == name_ofRdm4:
        irdm1Flag = False
        irdm1     = False
        irdm2Flag = False
        irdm2     = False
        for i in O1: 
          if i == t.indices[0]: irdm1Flag = True; irdm1 = O1.index(i)
          if i == t.indices[1]: irdm2Flag = True; irdm2 = O1.index(i)
        if irdm1Flag and irdm2Flag and irdm1 > irdm2:
          temp_i = O1[irdm2].copy()
          del O1[irdm2]
          O1.insert(irdm2, O1[irdm1-1])
          del O1[irdm1]
          O1.insert(irdm1, temp_i)

    #Order of Eri has to be the TOP!!!!
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      if t.name == name_ofEri and t.indices[extEri] in O1:
        O1.remove(t.indices[extEri])
        O1.insert(0, t.indices[extEri])
        
    t_list = [] # List of tensors already read from disk, or GA
    for i in O1:
      if   i.indType[0][0] == 'core':    Label = "{core}"
      elif i.indType[0][0] == 'active':  Label = "{occ}"
      elif i.indType[0][0] == 'virtual': Label = "{vir}"
      print Indent[LoopCount] + "for " + i.name + " in " + Label + ":"
      CPloop(Indent[LoopCount], i, CPfile)
      LoopCount += 1
      Declared.append(i)
      if isBareLTS:
        if   i == LTensor.indices[extBare]:
          print Indent[LoopCount] + "Read " + LTensor.name + " from Disk or GA for " + LTensor.indices[extBare].name 
          ReadRetval(Indent[LoopCount], LTensor, extBare, CPfile)
          SigmaCount = LoopCount
      for t in caseTensor[opt_prior]['Tensors'].tensors:
        if   t.name == name_ofEri     and i == t.indices[extEri]:
          print Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[extEri].name
          ReadERI(Indent[LoopCount], t, extEri, CPfile)
          t_list.append(t.name)
        elif t.name == name_ofBareAmp and i == t.indices[extBare]:         
          print Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[extBare].name
          ReadBare(Indent[LoopCount], t, extBare, CPfile)
          t_list.append(t.name)

    # In case of 4-RDM
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      if t.name == name_ofRdm4 and t.indices[0] in O1 and t.indices[1] in O1:
        print Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[0].name + "," + t.indices[1].name
        ReadD4(Indent[LoopCount], t, extRdm4, CPfile)
        D4Count = LoopCount
        t_list.append(t.name)         

    Ccount = 0
    Ocount = 0
    Vcount = 0
    printVar = Indent[LoopCount] + "Declare "
    for i in caseTensor[opt_prior]['Intermediate'].indices:
      if   not i in O1 and i.indType[0][0] == 'core'   : Ccount += 1
      elif not i in O1 and i.indType[0][0] == 'active' : Ocount += 1
      elif not i in O1 and i.indType[0][0] == 'virtual': Vcount += 1

    printVar += caseTensor[opt_prior]['Intermediate'].name + " as a "
    if Ccount: printVar += "c^" + str(Ccount)
    if Ocount: printVar += "o^" + str(Ocount)
    if Vcount: printVar += "v^" + str(Vcount)
    if not Ccount and not Ocount and not Vcount: printVar += "scalar"
    else: printVar += " tensor"    
    if not set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices):
      print printVar

      Var = []
      for c in range(Ccount): Var.append("nclosed")
      for o in range(Ocount): Var.append("nocc")
      for v in range(Vcount): Var.append("nvir") 
      DecX = "  %sorz::DTensor X(" % Indent[LoopCount]
      for num in range(len(Var)):
        DecX += Var[num]
        if not num == len(Var)-1: DecX += ", "
        else:  DecX += ");"
      extI = []
      for i in caseTensor[opt_prior]['Intermediate'].indices:
        if i.isExt: extI.append(i.copy())
      if not len(extI): symmLabel = "0"
      else:
        symmLabel = ""
        for i in extI:
          symmLabel += "s" + i.name
          if not extI.index(i) == len(extI)-1: symmLabel += "^"
      DecSymm = "  %sorz::DTensor X" % Indent[LoopCount] + "c" * Ccount + "a" * Ocount + "v" * Vcount + " = orz::ct::sympack_X" + "c" * Ccount + "a" * Ocount + "v" * Vcount  + "(" + "symblockinfo, " + symmLabel + ", X);"
      CPfile.write( DecX    + "\n" )
      CPfile.write( DecSymm + "\n" )      

    intSize = (numCore**Ccount) * (numOcc**Ocount) * (numVirt**Vcount)
    if intSize > MaxSize:
      print ""
      print "!!!! *WARNING* Size of intermediate exceeds the MaxSize !!!!"
      print chr(7)*100,
      print ""

    O2     = []
    printVar = ""
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      if   t.name == name_ofEri and not t.name in t_list:
        i = t.indices[extEri]

        if not t.indices[extEri] in Declared:
          O2.append(i)
          if   i.indType[0][0] == 'core':    Label = "{core}"
          elif i.indType[0][0] == 'active':  Label = "{occ}"
          elif i.indType[0][0] == 'virtual': Label = "{vir}"
          printVar += Indent[LoopCount] + "for " + i.name + " in " + Label + ":\n"
          CPloop(Indent[LoopCount], i, CPfile)
          LoopCount += 1
          Declared.append(i) 

        printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + i.name + "\n" 
        ReadERI(Indent[LoopCount], t, extEri, CPfile)
        t_list.append(t.name)
      elif t.name == name_ofBareAmp and not t.name in t_list:      
        i = t.indices[extBare]

        if not t.indices[extBare] in Declared:
          O2.append(i)
          if   i.indType[0][0] == 'core':    Label = "{core}"
          elif i.indType[0][0] == 'active':  Label = "{occ}"
          elif i.indType[0][0] == 'virtual': Label = "{vir}"
          printVar += Indent[LoopCount] + "for " + i.name + " in " + Label + ":\n"
          CPloop(Indent[LoopCount], i, CPfile)
          LoopCount += 1
          Declared.append(i) 

        printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + i.name + "\n"
        ReadBare(Indent[LoopCount], t, extBare, CPfile)
        t_list.append(t.name)
      elif t.name == name_ofRdm4 and not t.name in t_list:      
        i = [t.indices[0], t.indices[1]]
        O2 += i
        for j in i:
          if   j.indType[0][0] == 'core':    Label = "{core}"
          elif j.indType[0][0] == 'active':  Label = "{occ}"
          elif j.indType[0][0] == 'virtual': Label = "{vir}"
          if not j in Declared:
            printVar += Indent[LoopCount] + "for " + j.name + " in " + Label + ":\n"
            CPloop(Indent[LoopCount], j, CPfile)
            LoopCount += 1
            Declared.append(j)
          if i[0] in Declared and i[1] in Declared and not t.name in t_list:  
            printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + i[0].name + "," + i[1].name + "\n"
            ReadD4(Indent[LoopCount], t, extRdm4, CPfile)
            D4Count = LoopCount
            t_list.append(t.name)
      
    print printVar,

    iCount = 0
    T1name = caseTensor[opt_prior]['Intermediate'].name + "_("
    if not (caseTensor[opt_prior]['Intermediate'].name == LTensor.name and isBareLTS):
      for num_i in range(len(O1)):
        i = O1[num_i] 
        T1name += i.name
        if num_i != len(O1) - 1: T1name += "," 
      T1name += ")("
      for num_i in range(len(caseTensor[opt_prior]['Intermediate'].indices)):
        i = caseTensor[opt_prior]['Intermediate'].indices[num_i]
        if not i in O1: 
          T1name += i.name
          iCount += 1
          if iCount != (len(caseTensor[opt_prior]['Intermediate'].indices)-len(O1)): T1name += ","
          else: T1name += ")"
    else:
      T1name += "(" + LTensor.indices[extBare].name + ")(" + LTensor.indices[0].name + ", " + LTensor.indices[1].name \
              + ", " + LTensor.indices[2].name + ", )"

    Sumname = "sum("
    O3 = []
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      for i in t.indices: 
        if not i in O1 and not i in O2 and not i in caseTensor[opt_prior]['Intermediate'].indices: O3.append(i)

    O3 = list(set(O3))

    for num_i in range(len(O3)):
      i = O3[num_i]
      Sumname += i.name
      if    num_i != len(O3)-1: Sumname += ","
    Sumname += ")"

    T2names = " "
    for t in caseTensor[opt_prior]['Tensors'].tensors:
      T2names += t.name + "("
      for num_i in range(len(t.indices)):
        i = t.indices[num_i]
        if   t.name == name_ofEri     and not num_i == extEri:               T2names += i.name
        elif t.name == name_ofRdm4    and not num_i == 0 and not num_i == 1: T2names += i.name
        elif t.name == name_ofBareAmp and not num_i == extBare:              T2names += i.name
        elif not t.name == name_ofEri and not t.name == name_ofRdm4 and not t.name == name_ofBareAmp:
          T2names += i.name

        if num_i != len(t.indices)-1: T2names += ","

      T2names += ")"
      if t != caseTensor[opt_prior]['Tensors'].tensors[-1]: T2names += " * "

    if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices): 
      Constnames = str(thisTerm.numConstant) + " "
      for i in thisTerm.constants: Constnames += i
      Constnames += " "
    else: Constnames = str(1.0) + " "

    # Call the Fortran function
    if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices): 
      convertCpp_body(LTensor, [caseTensor[opt_prior]['Tensors']], title1if, CPfile, DecFirst, Indent[LoopCount], isBareLTS, \
        lts_name, name_ofEri, name_ofBareAmp) # *TEST*        
    else:
      convertCpp_body(caseTensor[opt_prior]['Intermediate'], [caseTensor[opt_prior]['Tensors']], title1if, \
        CPfile, DecFirst, Indent[LoopCount], False, lts_name, name_ofEri, name_ofBareAmp) # *TEST*

    print Indent[LoopCount] + T1name + " += " + Constnames + Sumname + T2names    

    if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices): 
      if isBareLTS:

        LTname  = LTensor.name + "_(" 
        if (isBareLTS): LTname += LTensor.indices[extBare].name
        LTname += ")("
        for num_i in range(len(LTensor.indices)-1):
          i = LTensor.indices[num_i]
          LTname += i.name
          if num_i != len(LTensor.indices)-2: LTname += ","
        LTname += ")"

        AccBare(Indent[SigmaCount], LTensor, extBare, CPfile)
        print 
        print Indent[SigmaCount] + "Accumulate " + LTname + " for " + LTensor.indices[extBare].name 

      print ""
      print "! ----------------------------------------------------------------------"
      print "! ----------------------------------------------------------------------"
      print ""

      if SigmaCount == False:
        raise RunTimeError, "Algorithmic Error"

      if isRotated and isBareLTS:
        LTensor.indices[extBare].isExt = False
        newIndices = [LTensor.indices[1], LTensor.indices[0], LTensor.indices[3], LTensor.indices[2]]
        LTensor.indices = newIndices[:]
        for i in LTensor.indices: i.isExt = False
        LTensor.indices[extBare].isExt = True
        isRotated = False

      if isKDkilled: 
        LTensor.indices = LTcopy.indices[:]
        isKDkilled = False

      #print D4Count
      for i in range(LoopCount+1, len(O1), -1):
        if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices): 
          if isBareLTS and i == SigmaCount:  
            AccBare(Indent[SigmaCount], LTensor, extBare, CPfile)
        if D4Count and i == D4Count:
          D4Count = False
          CPfile.write("  %sFC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();\n" % Indent[i])
        LoopEnd(Indent[i-2], CPfile)
  
      CPfile.write("  }\n")
      continue
    else:     print ""

    #print D4Count
    for i in range(LoopCount, len(O1), -1):
      if set(LTensor.indices) == set(caseTensor[opt_prior]['Intermediate'].indices): 
        if isBareLTS and i == SigmaCount:  
          AccBare(Indent[SigmaCount], LTensor, extBare, CPfile)
      if D4Count and i == D4Count:
        D4Count = False
        CPfile.write("  %sFC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();\n" % Indent[i])
      LoopEnd(Indent[i-1], CPfile)

    # ********* Second Contraction ************
    title2if = "g_if_" + title + "_no1_x" + str(numTerm) 
    title2   = "g_"    + title + "_no1_x" + str(numTerm) 

#     print "***** 5 IntermIndex *****" # *TEST*
#     for t in thisTerm.tensors:
#       print t
#       for i in t.indices:
#         print i.name, i.isExt 

    Set_i  = LTensor.indices[:]
    ExtInd = []
    NameT  = []
    Consts = thisTerm.constants
#    if not len(LTensor.indices): Consts.append(LTensor.name)
    for t in thisTerm.tensors:
      if (not 'E' in t.name) and (not t.name == 'kdelta'): NameT.append(t.name)
      Set_i.extend(t.indices[:])
    for i in Set_i:
#      print i.name + ", ", i.isExt # *TEST*
      if i.isExt == True: ExtInd.append(i.copy());
    ExtInd = list(set(ExtInd))
    NameT.append(LTensor.name)
    Set_i = list(Set_i)
    DecSecond = {"DecInd":ExtInd, "DecConst":Consts, "DecTensor":NameT}    
    # Declarations for Fortran code included from C++ 
    convertCpp_header(LTensor, [thisTerm], title2if, CHfile, DecSecond, isBareLTS, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
    # Interface inbetween C++ and Fortran parts
    convertF90_interface(LTensor, [thisTerm], title2if, F90file, DecSecond, isBareLTS, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
    # Body of the tensorisl contraction
    convertF90_ERI2(LTensor, [thisTerm], title2, F90file, DecSecond, isBareLTS, lts_name, name_ofEri, name_ofBareAmp) # *TEST*

    if not isOF and not set_V and num_V:
      temp_i = []
      for t in thisTerm.tensors:
        if t.name == 'X':
          for i in t.indices:
            if i.isExt == True: i.isExt = False; temp_i.append(i.copy())
      for i in temp_i:
        for t in thisTerm.tensors:
          if i in t.indices:
            for j in t.indices:
              if j == i: j.isExt = False

    # Starts the second contraction ....
    Declared = []
    for i in O1:
      Declared.append(i.copy())

    LoopCount = len(O1)
    printVar = ""
    if isBareLTS and not LTensor.indices[extBare] in O1:
      if   LTensor.indices[extBare].indType[0][0] == 'core':    Label = "{core}"
      elif LTensor.indices[extBare].indType[0][0] == 'active':  Label = "{occ}"
      elif LTensor.indices[extBare].indType[0][0] == 'virtual': Label = "{vir}"
      printVar += Indent[LoopCount] + "for " + LTensor.indices[extBare].name + " in " + Label + ":\n"
      CPloop(Indent[LoopCount], LTensor.indices[extBare], CPfile)
      LoopCount += 1
      Declared.append(LTensor.indices[extBare])

      printVar += Indent[LoopCount] + "Read " + LTensor.name + " from Disk or GA for " + LTensor.indices[extBare].name       
      ReadRetval(Indent[LoopCount], LTensor, extBare, CPfile)
      print printVar
      SigmaCount = LoopCount

    O4 = []
    t_list = []
    printVar = ""
    for t in thisTerm.tensors:
      if   t.name == name_ofEri:     O4.append(t.indices[extEri])
      elif t.name == name_ofBareAmp: O4.append(t.indices[extBare])
      elif t.name == name_ofRdm4:
        O4.append(t.indices[0])
        O4.append(t.indices[1])
    for i in O4:
      if   i.indType[0][0] == 'core':    Label = "{core}"
      elif i.indType[0][0] == 'active':  Label = "{occ}"
      elif i.indType[0][0] == 'virtual': Label = "{vir}"

      if not i in Declared: 
        print Indent[LoopCount] + "for " + i.name + " in " + Label
        CPloop(Indent[LoopCount], i, CPfile)
        LoopCount += 1
        Declared.append(i)

      for t in thisTerm.tensors:
        if   t.name == name_ofEri     and not t.name in t_list:        
          printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[extEri].name       
          ReadERI(Indent[LoopCount], t, extEri, CPfile)
          t_list.append(t.name)
          print printVar
        elif t.name == name_ofBareAmp and not t.name in t_list:        
          printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[extBare].name       
          ReadBare(Indent[LoopCount], t, extBare, CPfile)
          t_list.append(t.name)
          print printVar
        elif t.name == name_ofRdm4    and not t.name in t_list and t.indices[0] in Declared and t.indices[1] in Declared:        
          printVar += Indent[LoopCount] + "Read " + t.name + " from Disk or GA for " + t.indices[0].name + "," + t.indices[1].name       
          ReadD4(Indent[LoopCount], t, extRdm4, CPfile)
          D4Count = LoopCount
          t_list.append(t.name)
          print printVar
        printVar = ""
      
    IndList = set(LTensor.indices)
    for num_t in range(len(thisTerm.tensors)-1):
      this_t = thisTerm.tensors[num_t  ]
      next_t = thisTerm.tensors[num_t+1]

      set_i   = set(this_t.indices).union(set(next_t.indices))
      IndList = set(IndList).union(set_i)

    UnDeclared = set([])
    UnDeclared = IndList.difference(set(Declared))

    tempList = set(LTensor.indices)
    if (isBareLTS): tempList.remove(LTensor.indices[extBare])
    UnDeclared = UnDeclared.difference(set(tempList))
    UnDeclared = list(UnDeclared)

    LTname  = LTensor.name + "_(" 
    if (isBareLTS): LTname += LTensor.indices[extBare].name
    LTname += ")("
    for num_i in range(len(LTensor.indices)-1):
      i = LTensor.indices[num_i]
      LTname += i.name
      if num_i != len(LTensor.indices)-2: LTname += ","
    LTname += ")"

    Sumname = "sum("
    for num_i in range(len(UnDeclared)):
      i = UnDeclared[num_i]
      Sumname += i.name
      if num_i != len(UnDeclared)-1: Sumname += ","
    Sumname += ")"

    T2names = " "
    for t in thisTerm.tensors:
      if   t.name != 'X': T2names += t.name + "(" # This can be different
      else:
        T2names += t.name + "_("
        for i in O1:
          T2names += i.name
          if i != O1[-1]: T2names += ","
        T2names += ")("
      for num_i in range(len(t.indices)):
        i = t.indices[num_i]
        if   t.name == name_ofEri     and not num_i == extEri:               T2names += i.name
        elif t.name == name_ofRdm4    and not num_i == 0 and not num_i == 1: T2names += i.name
        elif t.name == name_ofBareAmp and not num_i == extBare:              T2names += i.name
        elif t.name == 'X'            and not i in O1:                       T2names += i.name # This can be different
        elif not t.name == name_ofEri and not t.name == name_ofRdm4 and not t.name == name_ofBareAmp and not t.name == 'X':
          T2names += i.name

        if num_i != len(t.indices)-1: T2names += ","

      T2names += ")"
      if t != thisTerm.tensors[-1]: T2names += " * "

    Constnames = str(thisTerm.numConstant)
    for i in thisTerm.constants: Constnames += " * " + i + " * "
    Constnames += " "

    # Call the Fortran function
    convertCpp_body(LTensor, [thisTerm], title2if, CPfile, DecSecond, Indent[LoopCount], isBareLTS, lts_name, name_ofEri, name_ofBareAmp) # *TEST*
    print Indent[LoopCount] + LTname + " += " + Constnames + Sumname + T2names
    if isBareLTS:
      print 
      print Indent[SigmaCount] + "Accumulate " + LTname + " for " + LTensor.indices[extBare].name 

    for i in range(LoopCount, 0, -1):
      if isBareLTS and i == SigmaCount: AccBare(Indent[SigmaCount], LTensor, extBare, CPfile)
      if D4Count and i == D4Count:
        CPfile.write("  %sFC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();\n" % Indent[i])
      LoopEnd(Indent[i-1], CPfile)

    CPfile.write("  }\n")

    print ""
    print "! ----------------------------------------------------------------------"
    print "! ----------------------------------------------------------------------"
    print ""

    if SigmaCount == False and isBareLTS:
      raise RunTimeError, "Algorithmic Error"

    if isRotated and isBareLTS:
      LTensor.indices[extBare].isExt = False
      newIndices = [LTensor.indices[1], LTensor.indices[0], LTensor.indices[3], LTensor.indices[2]]
      LTensor.indices = newIndices[:]
      for i in LTensor.indices: i.isExt = False
      LTensor.indices[extBare].isExt = True
      isRotated = False

    if isKDkilled: 
      LTensor.indices = LTcopy.indices[:]
      isKDkilled = False

  if not len(LTensor.indices) and BareFlag:
    CPfile.write(\
"""

  return %s;
}

""" % LTensor.name)
  else:
    CPfile.write (\
"""

  return retval;
}

""")

  CHfile.write( CHend )
  CHfile.close()
  F90file.close()
  CPfile.close()



#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_interface(LTensor, inTerm, title, Ffile, Declare, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Print out the interface inbetween C++ and Fortran parts"

  # Check is Declare is of a dictionary
  if not type(Declare) == type({}):
    raise TypeError, "Declare should be given as a diactionary"

  # Check if Ffile is a file object
  aaafile = open('aaa.out', 'w')
  if not type(Ffile) == type(aaafile):
    raise TypeError, 'Ffile shoud be given as a file object'

  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a atring"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORZ code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  # Declaration of the external variables
  ExtIndices = Declare["DecInd"]
  Consts     = Declare["DecConst"]
  NameT      = Declare["DecTensor"]

  printGuard  = "\n\n"
  printGuard += "! **********************************************************************\n"
  printGuard += "!                                                                       \n"
  printGuard += "! **********************************************************************\n"

  printName = "subroutine " + title + "("  

  for i in ExtIndices:
    printName += "s" + i.name + ", i" + i.name + ", "
  for c in Consts:
    printName += c + ", "
  for t in NameT:
    if t == 'kdelta': continue
    printName += t + ", "

  printName += "nir, nsym, psym)\n\n"
  DeclareComm8  = "use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6\n"
  DeclareComm8 += "use mod_g_if_comm08\n"
  DeclareComm8 += "implicit none\n"
  DeclareComm8 += "\n"


  Ffile.write( printGuard   )
  Ffile.write( printName    )
  Ffile.write( DeclareComm8 )

  if len(ExtIndices): 
    DeclareInd = "integer, intent(inout) :: "
    for i in ExtIndices:
      DeclareInd += "s" + i.name + ", i" + i.name
      if ExtIndices.index(i) != len(ExtIndices)-1: DeclareInd += ", "
    DeclareInd += "\n"
    Ffile.write(DeclareInd)
  if len(Consts):
    DeclareConst = "real(kind=8), intent(inout) :: "
    for i in Consts:
      DeclareConst += i
      if Consts.index(i) != len(Consts)-1: DeclareConst += ", "
    DeclareConst += "\n"
    Ffile.write( DeclareConst )
  if len(NameT):
    DeclareTensor = "real(kind=8), intent(inout) :: "
    for i in NameT:
      if len(LTensor.indices) == 0 and i == LTensor.name:
        DeclareTensor += i
      else:
        DeclareTensor += i + "(*)"
      if NameT.index(i) != len(NameT)-1: DeclareTensor += ", "
    DeclareTensor += "\n"
    Ffile.write( DeclareTensor )   

  DeclareSymm  = "! Information of the Irreps ....\n" 
  DeclareSymm += "integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)\n"  
  DeclareSymm += "! Some extra stuff\n"
  DeclareSymm += "integer :: sleft\n\n" 
  Ffile.write( DeclareSymm )

  thisTerm = inTerm[0]
  extTdict = {}
  for t in thisTerm.tensors+[LTensor]:
    if t.name == 'X':
      SetX =  ""
      exCore = []
      inCore = []
      exAct  = []
      inAct  = []
      exVir  = []
      inVir  = []
      for i in t.indices:
        if       i.isExt and i.indType[0][0] == 'core'   : exCore.append(i.copy())
        elif     i.isExt and i.indType[0][0] == 'active' : exAct.append(i.copy())
        elif     i.isExt and i.indType[0][0] == 'virtual': exVir.append(i.copy())
        elif not i.isExt and i.indType[0][0] == 'core'   : inCore.append(i.copy())
        elif not i.isExt and i.indType[0][0] == 'active' : inAct.append(i.copy())
        elif not i.isExt and i.indType[0][0] == 'virtual': inVir.append(i.copy())
      tempList = exCore + exAct + exVir
      if   len(tempList) == 0: SetX += "sleft = 0"
      elif len(tempList) == 1:
        SetX += "sleft = s" + tempList[0].name + "\n"
      elif len(tempList) == 2:
        SetX += "sleft = IEOR(s" + tempList[0].name + ",s" + tempList[1].name + ")\n"
      elif len(tempList) == 3:
        SetX += "sleft = IEOR(s" + tempList[0].name + ",IEOR(s" + tempList[1].name + ",s" + tempList[2].name + "))\n"
      elif len(tempList) == 4:
        SetX += "sleft = IEOR(IEOR(" + tempList[0].name + "," + tempList[1].name + "),IEOR(" \
             + tempList[2].name + "," + tempList[3].name + ")"

      SetX += "\n"
      Ffile.write( SetX )

      nameX = "x" + "c" * len(inCore) + "a" * len(inAct) + "v" * len(inVir)
      SetX = "call set_symblock_" + nameX + "(sleft, x, nir, nsym, psym) ! -> " + nameX + " (allocate) \n" 
      Ffile.write( SetX )
      extTdict['X'] = nameX

    elif t.name == name_ofBareAmp:
      iname = t.indices[extBare].name
      SetT = "call set_symblock_av2(s" + iname + ", " + t.name + ", nir, nsym, psym) ! -> av2_i (allocate)\n"
      extTdict[name_ofBareAmp] = "av2_i"
      Ffile.write( SetT )

    elif t.name == 'h':
      Seth1 = "call set_symblock_h1(" + t.name + ", nir, nsym, psym) ! -> h1 (allocate)\n"
      extTdict[t.name] = "h1"
      Ffile.write( Seth1 )

    elif t.name == name_ofEri:
      iname = t.indices[extEri].name
      SetV2 = "call set_symblock_h2(s" + iname + ", " + t.name + ", nir, nsym, psym) ! -> h2_i (allocate)\n"
      extTdict[name_ofEri] = "h2_i"
      Ffile.write( SetV2 )

  if isBareLTS: 
    extTdict[LTensor.name] = "av2_i2"
    SetLTensor = "call set_symblock_av2_2(s" + LTensor.indices[extBare].name + ", " + LTensor.name + ", nir, nsym, psym) ! -> av2_i2 (allocate)\n"
    Ffile.write( SetLTensor )

  FuncName = list(title)
  FuncName.pop(1)
  FuncName.pop(1)
  FuncName.pop(1)
  FuncName = "".join(FuncName)
  SetFunc = "call " + FuncName + "("
  
  # Declaration of the external variables
#   ExtIndices = Declare["DecInd"]
#   Consts     = Declare["DecConst"]
#   NameT      = Declare["DecTensor"]
  for i in ExtIndices:
    SetFunc += "s" + i.name + ", i" + i.name + ", "
  for c in Consts:
    SetFunc += c + ", "
  for t in NameT:
    if t == 'kdelta': continue
    if t in extTdict.keys():
      name_t = extTdict[t]
      SetFunc += name_t + ", "
    else: 
      SetFunc += t + ", "

  E_list = []
  for t in thisTerm.tensors:
    if 'E' in t.name: E_list.append(t.name)

  for i in E_list:
    if   i == "E1": SetFunc += "d1, "
    elif i == "E2": SetFunc += "d2, "
    elif i == "E3": SetFunc += "d3, "
    elif i == "E4": SetFunc += "d4_ij, "
    else:
      raise TypeError, "Algorithmic Error"

  SetFunc += "nir, nsym, psym)\n\n"
  Ffile.write( SetFunc )

  for t in extTdict.keys():
    SetDealloc = "deallocate(" + extTdict[t] + ")\n"
    Ffile.write( SetDealloc )
  Ffile.write( "\n" )

  nameEnd = "end subroutine " + title + "\n\n\n"
  Ffile.write( nameEnd )


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertCpp_header(LTensor, inTerm, title, CHfile, Declare, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Print out the interface inbetween C++ and Fortran parts"

  # Check is Declare is of a dictionary
  if not type(Declare) == type({}):
    raise TypeError, "Declare should be given as a diactionary"

  # Check if Ffile is a file object
  aaafile = open('aaa.out', 'w')
  if not type(CHfile) == type(aaafile):
    raise TypeError, 'Ffile shoud be given as a file object'

  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a atring"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORZ code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  # Declaration of the external variables
  ExtIndices = Declare["DecInd"]
  Consts     = Declare["DecConst"]
  NameT      = Declare["DecTensor"]

  FuncName = "void FC_FUNC(" + title.lower() + ", " + title.upper() + ")\n"
  FuncBody = "  ("

  num = 0 
  for i in ExtIndices:
    FuncBody += "const FC_INT &s" + i.name + ", const FC_INT &i" + i.name + ", "
    num += 1
  if len(ExtIndices): FuncBody += "\n   "
  for c in Consts:
   FuncBody += "const double * const " + c + ", "
  num += 1
  if len(Consts): FuncBody += "\n   "
  for t in NameT:
    if t == 'kdelta': continue
    FuncBody += "const double * const " + t + ", "
    num += 1
  if len(NameT): FuncBody += "\n   "
  FuncBody += "const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);\n\n"

  CHfile.write( FuncName )
  CHfile.write( FuncBody )


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertCpp_body(LTensor, inTerm, title, CPfile, Declare, Indents, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Print out the calling of the Fortran function from C++ code"

  # Check is Declare is of a dictionary
  if not type(Declare) == type({}):
    raise TypeError, "Declare should be given as a diactionary"

  # Check if Ffile is a file object
  aaafile = open('aaa.out', 'w')
  if not type(CPfile) == type(aaafile):
    raise TypeError, 'CPfile shoud be given as a file object'

  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a atring"

  # Check if Indents is given correctly
  if not type(Indents) == type('aaa'):
    raise TypeError, "Indents should be given as a atring"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORZ code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  # Declaration of the external variables
  ExtIndices = Declare["DecInd"]
  Consts     = Declare["DecConst"]
  NameT      = Declare["DecTensor"]

  FuncName = ("  %sFC_FUNC(" % Indents) + title.lower() + ", " + title.upper() + ")\n" 
  FuncBody =  "    %s("        % Indents

  num = 0 
  for i in ExtIndices:
    FuncBody += "s" + i.name + ", i" + i.name + ", "
    num += 1
#  if len(ExtIndices): FuncBody += "\n   "
  for c in Consts:
   FuncBody += "&" + c + ", "
  num += 1
#  if len(Consts): FuncBody += "\n   "
  for t in NameT:
    if   t == 'kdelta': continue  # Skip kdeltas                        
    elif t == name_ofEri: # in case of ERI
      FuncBody += t + "_sym.cptr(), "
    elif t == name_ofBareAmp: # in case of Amplitude 
      FuncBody += t + "b.cptr(), "
    elif t == 'X': # in case of the intermediate
      Ccount = 0
      Ocount = 0
      Vcount = 0
      if LTensor.name != 'X':
        for ten in inTerm[0].tensors:
          if ten.name == 'X':
            for i in ten.indices:
              if   i.indType[0][0] == 'core'    and not i.isExt: Ccount += 1
              elif i.indType[0][0] == 'active'  and not i.isExt: Ocount += 1
              elif i.indType[0][0] == 'virtual' and not i.isExt: Vcount += 1
      else: 
        for i in LTensor.indices:
          if   i.indType[0][0] == 'core'    and not i.isExt: Ccount += 1
          elif i.indType[0][0] == 'active'  and not i.isExt: Ocount += 1
          elif i.indType[0][0] == 'virtual' and not i.isExt: Vcount += 1
      FuncBody += t + 'c' * Ccount + 'a' * Ocount + 'v' * Vcount + '.cptr(), '
    elif t == 'h': # in case of one-body int
      FuncBody += 'moint1_sym.cptr(), '
    elif 'E' in t: continue # in case of RDM's
    elif t == LTensor.name and LTensor.name != 'X':
      if    isBareLTS: FuncBody += t + "b.cptr(), "
      elif  not len(LTensor.indices): FuncBody += "&" + t + ", "
      else: FuncBody += t + ".cptr(), " 
    else:
      raise TypeError, "I cannot handle this .... %s" % t
    num += 1
#  if len(NameT): FuncBody += "\n   "
  FuncBody += "nir, nsym, psym);\n"

  CPfile.write( FuncName )
  CPfile.write( FuncBody )


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_ERI2(LTensor, inTerm, title, Ffile, Declare, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "Separate given many-body equations in terms of the types of ERIs contained i.e. V[k*|**], V[a*|**], and 4-RDM2."
  "Then, interface to convertF90 function"

  # Check is Declare is of a dictionary
  if not type(Declare) == type({}):
    raise TypeError, "Declare should be given as a diactionary"

  # Check if Ffile is a file object
  aaafile = open('aaa.out', 'w')
  if not type(Ffile) == type(aaafile):
    raise TypeError, 'Ffile shoud be given as a file object'

  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a atring"

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORZ code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

#   print ""
#   print "! Searching in terms .... "
#   print ""
  TermLeft    = []     # List composed of terms with neither ERI, nor 4-RDMs
  TermVc      = []     # List composed of terms that contain Vi type ERI, except those with 4-RDMs
  TermVi      = []     # List composed of terms that contain Vi type ERI, except those with 4-RDMs
  TermVa      = []     # List composed of terms that contain Va type ERI, except those with 4-RDMs
  Term_Rdm4   = []     # List composed if terms that contain no ERIs and 4-RDMs 
  TermVc_Rdm4 = []     # List composed if terms that contain both Vi and 4-RDMs 
  TermVi_Rdm4 = []     # List composed if terms that contain both Vi and 4-RDMs 
  TermVa_Rdm4 = []     # List composed if terms that contain both Va and 4-RDMs 
  isVc        = False  # Flag whether a term contains the Vi ERI
  isVi        = False  # Flag whether a term contains the Vi ERI
  isVa        = False  # Flag whether a term contains the Va ERI
  isRdm4      = False  # Flag wgether a term contains the 4-RDMs
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
      if   thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'core':    isVc = True
      elif thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'active':  isVi = True
      elif thisTensor.name == name_ofEri and thisTensor.indices[extEri].indType[0][0] == 'virtual': isVa = True 
      elif thisTensor.name == name_ofRdm4: isRdm4 = True
    if   isVc and     isRdm4: TermVc.append(thisTerm)
    elif isVc and not isRdm4: TermVc_Rdm4.append(thisTerm)
    elif isVi and     isRdm4: TermVi_Rdm4.append(thisTerm)
    elif isVi and not isRdm4: TermVi.append(thisTerm) 
    elif isVa and     isRdm4: TermVa_Rdm4.append(thisTerm)
    elif isVa and not isRdm4: TermVa.append(thisTerm);
    elif not isVc and not isVi and not isVa and isRdm4:
                              Term_Rdm4.append(thisTerm)
    elif not isVc and not isVi and not isVa and not isRdm4:                     
                              TermLeft.append(thisTerm)
    isVc   = False
    isVi   = False
    isVa   = False
    isRdm4 = False

#   print "!  **** SUMMARY *****"
#   print "! 0. Number of terms with neither ERI, nor 4-RDMs : " + str(len(TermLeft))
#   print "! 1. Number of terms with V[i*|**] and no 4-RDM   : " + str(len(TermVi))
#   print "! 2. Number of terms with V[i*|**] and 4-RDM      : " + str(len(TermVi_Rdm4))
#   print "! 3. Number of terms with V[a*|**] and no 4-RDM   : " + str(len(TermVa))
#   print "! 4. Number of terms with V[a*|**] and 4-RDM      : " + str(len(TermVa_Rdm4))
#   print "! 4. Number of terms with no ERIs  and 4-RDM      : " + str(len(Term_Rdm4))
#   print "! ** Number of Total terms                        : " + str(len(inTerm))
#   print ""
  if len(TermLeft) + len(TermVc) + len(TermVi) + len(TermVa) + len(TermVc_Rdm4) + len(TermVi_Rdm4) + len(TermVa_Rdm4) + len(Term_Rdm4)!= len(inTerm):
    raise TypeError, "Algorithmic error .... "

  if len(TermLeft) > 0:
    Ffile.write( "! Converting terms with neither ERI, nor 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermLeft, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVc) > 0:
    Ffile.write( "! Converting terms with V[c*|**] and no 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVc, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVc_Rdm4) > 0:
    Ffile.write( "! Converting terms with V[c*|**] and with 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVc_Rdm4, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVi) > 0:
    Ffile.write( "! Converting terms with V[i*|**] and no 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVi, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVi_Rdm4) > 0:
    Ffile.write( "! Converting terms with V[i*|**] with 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVi_Rdm4, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVa) > 0:
    Ffile.write( "! Converting terms with V[a*|**] and no 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVa, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(TermVa_Rdm4) > 0:
    Ffile.write( "! Converting terms with V[a*|**] with 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, TermVa_Rdm4, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )
  if len(Term_Rdm4) > 0:
    Ffile.write( "! Converting terms with no ERIs with 4-RDMs ..... \n\n" )
    convertF90_2(LTensor, Term_Rdm4, title, Ffile, Declare, isBareLTS, lts_name, name_ofEri, name_ofBareAmp)
    Ffile.write( "! -----------------------------------------------------------------------\n" )


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def convertF90_2(LTensor, inTerm, title, F90file, Declare, isBareLTS = False, lts_name = 'sig', name_ofEri = 'V', name_ofBareAmp = 'S2_i'):
  "A simplified version of convertF90_1, designed only for sqa.factorize"

  # Check is Declare is of a dictionary
  if not type(Declare) == type({}):
    raise TypeError, "Declare should be given as a diactionary"

  # Check if Ffile is a file object
  aaafile = open('aaa.out', 'w')
  if not type(F90file) == type(aaafile):
    raise TypeError, 'Ffile shoud be given as a file object'

  # Check if title is given correctly
  if not type(title) == type('aaa'):
    raise TypeError, "title should be given as a atring"

  # In case of that the LHS is a Scalar 
  if len(inTerm) == 0:
    print "! RHS is ZERO .... "
    return

  # Check that LTensor is of tensor class
  if not isinstance(LTensor, tensor):
    raise TypeError, "LTensor must be a tensor object"

  # Check that isBareLTS is of boolean
  if not isinstance(isBareLTS, type(True)):
    raise TypeError, "isBareLTS must be True, or False"

  # Check that name_ofEri is of string
  if not isinstance(name_ofEri, type('a')):
    raise TypeError, "name_ofEri must be a string" 

  # Check that name_ofBareAmp is of string
  if not isinstance(name_ofBareAmp, type('a')):
    raise TypeError, "name_ofBareAmp must be a string" 

  # check that inTerm is a list of terms
  if not isinstance(inTerm, type([])):
    raise TypeError, "inTerm must be provided as a list"
  if not len(inTerm):
    raise TypeError, "Length of inTerm must be 1 (for now)"
  if len(inTerm[0].tensors) > 2:
    print len(inTerm[0].tensors)
    raise TypeError, "Number of tensors must be 2 in this code. If you wanna treat more, use convertF90_1"
  for num in range(len(inTerm)):
    if not isinstance(inTerm[num], term):
      raise TypeError, "Each element of inTerm must be class term"
  
  # determine what types of operators the term contains
  for num in range(len(inTerm)):
    has_creDesOps = False
    has_sfExOps = False
    for t in inTerm[num].tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        has_creDesOps = True
      elif isinstance(t, sfExOp):
        has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps:
    raise RuntimeError, "creOp/desOp cannot be implemenented in ORz code as tensorial quantities" + \
                        "creOp/desOp are present " + "in term " + num

  # name_ofEri must be provided as a name of electron repulsion integral
  if not isinstance(name_ofEri, type("a")):
    raise TypeError, "name_ofEri has to be provided as a charactar"

  # names_ofBareAmp has to be given as a name of tensors in the BareAmpPack
  if not isinstance(name_ofBareAmp, type("a")):
    raise TypeError, "names_ofAmp have to be given as a charactar"      
      
  # check that LTensor, which is supposed to represent a sigma vector, is a class tensor
  if not isinstance(inTerm[num], term):
    raise TypeError, "LTensor must be class term"

  # check whether the LTensor could be a BareAmp if specified so ....
  if len(LTensor.indices) != 4 and isBareLTS:
    raise TypeError, "LTensor could not be a BareAmpPack .... "

  # Order of the externally summed indices
  extEri  = 0        # in case of Eri
  extBare = 3        # in case of the BareAmpPack 
  extRdm4 = [0, 1]   # in case of the 4-RDM
  name_ofRdm4 = 'E4' # Name of 4-RDM

  F90file.write( "\n" )                                                                                  
  F90file.write( "! ---------------------------- Parameters used ------------------------------\n" )
  F90file.write( "!                                                                            \n" )    
  F90file.write( "! Whether the LHS is a BareAmpPack ....... " + str(isBareLTS)             + "\n" )
  F90file.write( "! LHS name ............................... " + lts_name                   + "\n" )
  F90file.write( "! Name of ERI ............................ " + name_ofEri                 + "\n" )
  F90file.write( "! Name of BareAmpPack appearing in RHS.... " + name_ofBareAmp             + "\n" )
  F90file.write( "!                                                                            \n" )
  F90file.write( "!----------------------------------------------------------------------------\n" )  
#  F90file.write( "\n" )

  printName  = "! subroutine " + title + "("
  printName2 = "subroutine " + title + "("

  printEnd  = "end subroutine " + title + "\n\n"

  d = datetime.datetime.today()
  printData  =  "! Generated time: " + str(d.year) + "/" + str(d.month) + "/" + str(d.day)
  printData +=  " " + str(d.hour) + ":" + str(d.minute) + ":" + str(d.second) 
  printData += "\n\n"

  DeclareComm8  = "! FEMTO BEGIN  **************************************************************\n"
  DeclareComm8 += "use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6\n"
  DeclareComm8 += "implicit none\n\n"
  DeclareExtInd = ""

  tempSet = set([])
  Tensor1 = inTerm[0].tensors[0]
  if len(inTerm[0].tensors) == 2: 
    Tensor2 = inTerm[0].tensors[1]
    tempSet = set(Tensor1.indices).union(set(Tensor2.indices))
  else:                           
    tempSet = set(Tensor1.indices[:]) 
  tempSet = list(tempSet)
  for i in tempSet:
    if i.isExt == True:
      DeclareExtInd += "integer, intent(in) :: i_" +  i.name + ", s_" + i.name + "\n"
      printName     += "i_" + i.name + ", s_" + i.name + ", "
  
  DeclareSymm  = "! Information of the Irreps ....\n" 
  DeclareSymm += "integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)\n"  

  # Print all the name of parameters .... 
  list_ofAllConst = []
  for thisTerm in inTerm:
    for thisConst in thisTerm.constants:
      list_ofAllConst.append(thisConst)

  DeclareConst = ""
  list_ofAllConst = list(set(list_ofAllConst))
  if len(list_ofAllConst) != 0:
    DeclareConst = "! Declaration of numerical constants .... \n"
  for thisConst in list_ofAllConst:
    DeclareConst += "real(kind=8), intent(inout) :: " + str(thisConst)
    printName += str(thisConst) + ", "

  # Print all the names of tensors .....
  if isBareLTS: 
    dict_ofAllTensors = {LTensor.name: 3}
  elif not isBareLTS and LTensor.name == 'X': # This can be different
    NDcount = 0
    for i in LTensor.indices:
      if i.isExt == False: NDcount += 1
    dict_ofAllTensors = {LTensor.name: NDcount}
  else:
    dict_ofAllTensors = {LTensor.name: len(LTensor.indices)}
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
#      print thisTensor.name # *TEST*
      if not dict_ofAllTensors.has_key(thisTensor.name):
        if thisTensor.name == name_ofBareAmp or thisTensor.name == name_ofEri:
          dict_ofAllTensors[thisTensor.name] = 3
        elif thisTensor.name == name_ofRdm4:
          dict_ofAllTensors[thisTensor.name] = 6
        elif thisTensor.name == 'X': # This can be different
          NDcount = 0
          for i in thisTensor.indices:
            if i.isExt == False: NDcount += 1      
          dict_ofAllTensors[thisTensor.name] = NDcount
        else:
          dict_ofAllTensors[thisTensor.name] = len(thisTensor.indices)

  # Declaration of the external variables
  ExtIndices = Declare["DecInd"]
  Consts     = Declare["DecConst"]
  NameT      = Declare["DecTensor"]
  for i in ExtIndices:
    printName2 += "s_" + i.name + ", i_" + i.name + ", "
  for c in Consts:
    printName2 += c + ", "
  for t in NameT:
    if t == 'kdelta': continue
    printName2 += t + "_, "

#  global str
  DeclareTensor = """\n! Declare tensors used ...\n"""
  for name in dict_ofAllTensors.keys():
    if 'E' in name: printName2 += name + "_, "
    if name == 'kdelta': continue
    dimTensor = dict_ofAllTensors[name]
    printName += name + "_, "
    if dimTensor == 0:
      DeclareTensor += "real(kind=8)                   :: " + name + "_\n"
      continue
    else:
      decTensor = "type(symblock" + str(dimTensor) + "), intent(inout) :: " + name + "_("    
    for i in range(dimTensor): 
      decTensor+= "0:nir-1"
      if i == dimTensor-1: decTensor += ")"; DeclareTensor += decTensor + "\n"
      else:                decTensor += ", "
  DeclareTensor += "\n"
  F90file.write( "\n" )

  printName  += "nir, nsym, psym)\n\n"
  printName2 += "nir, nsym, psym)\n\n"
  
  F90file.write( printData     ) 
#  F90file.write( printName     )
  F90file.write( printName2    )
  F90file.write( DeclareComm8  )
  F90file.write( DeclareExtInd )
  F90file.write( DeclareSymm   )
  F90file.write( DeclareConst  )
  F90file.write( DeclareTensor )

  # First of all, make the list of all the variables used in the subsequent tensorial contraction
  list_ofIndices = []
  for thisIndex in LTensor.indices: list_ofIndices.append(thisIndex.name)
  for thisTerm in inTerm:
    for thisTensor in thisTerm.tensors:
      for thisIndex in thisTensor.indices:
        list_ofIndices.append(thisIndex.name)

  union_ofIndices = set(list_ofIndices)
  list_ofIndices = list(union_ofIndices)

  NDindices = []
  for t in thisTerm.tensors:
    for i in t.indices:
      if not i.isExt: NDindices.append(i)
  NDindices = list(set(NDindices))

  F90file.write( "! Indices used in the contractions as dummy ... \n" )
  for num in range(len(NDindices)):
    thisName = NDindices[num].name
    if num % 5 != 0: F90file.write( ", s_" + thisName + ", i_" + thisName ) 
    else:
      F90file.write( "\n" )          
      F90file.write( "integer :: "                       )
      F90file.write( "s_" + thisName + ", i_" + thisName )

#  F90file.write( "\n" )
  F90file.write( "\n" )

  # If the LHS is a BareAmp, set so ...
  if isBareLTS == True: 
    LTensor.indices[extBare].isExt = True

  F90file.write( "\n" )
  # Search in each term ...
  num_loops    = []
#  summed_Index = []
  for numTerm in range(len(inTerm)):
    thisTerm = inTerm[numTerm]
    F90file.write( "!" + LTensor.__str__()  + " <-- \n" )
    F90file.write( "!" + thisTerm.__str__() + "\n" )

    # Dictionary that shows which index corresponds to ieri, isig, and so on ...
    extDict = {}
    # List of all external indices without any redundancy
    extIndices = []
    # If the LHS is a BareAmp, set so ...
    if isBareLTS == True: 
      extIndices.append(LTensor.indices[extBare]) # Record which indices are of external
      extDict["i"+lts_name] = LTensor.indices[extBare]

    for thisTensor in thisTerm.tensors:
      # If thisTensor is ERI, or BareAmp, the external indices have to be dropped off ...
      for thisIndex in thisTensor.indices:
        if thisIndex.isExt == True: extIndices.append(thisIndex)        
    # Avoid double counting     
    extIndices = list(set(extIndices)) 

    # Search inbetween each tensor ...
    loop_count  = 0
    num_kdeltas = 0
    summed_Index = set(LTensor.indices)
    for numTensor in range(len(thisTerm.tensors)-1):
      thisTensor = thisTerm.tensors[numTensor]
      nextTensor = thisTerm.tensors[numTensor+1]

      set_i1       = set(thisTensor.indices)
      set_i2       = set(nextTensor.indices)
      set_i3       = set_i1.union(set_i2)
      summed_Index = set_i3.union(summed_Index)

    # Print whether the leg of ERI belongs to {occ} or {vir}
    # Core is not considered for now
    for thisTensor in thisTerm.tensors:
      if thisTensor.name == name_ofEri:
        F90file.write( "! where i_" + thisTensor.indices[extEri].name + " in " )
        if   thisTensor.indices[extEri].indType[0][0] == 'active':
          F90file.write( "{core}" + "\n" )
        elif thisTensor.indices[extEri].indType[0][0] == 'core':
          F90file.write( "{occ}" + "\n" )
        elif thisTensor.indices[extEri].indType[0][0] == 'virtual':
          F90file.write( "{vir}" + "\n" )
        else: raise TypeError, "Algorithmic Error"

    # Write down the loops over irreps ...
    for index in summed_Index:
      if index.isExt == True:
        continue      
      F90file.write( "do s_" + index.name + " = 0, nir-1" + "\n" )
      loop_count += 1

    # Write down the constraints in terms of the irreps ...
    # Conditions for the left-hand tensor
    F90file.write( "if( &" + "\n" )
    thisTensor = LTensor
    if len(thisTensor.indices) == 1:   # in case of unity
      F90file.write( "S_" + thisTensor.indices[0].name + " == 0" )
    elif len(thisTensor.indices) == 2: # in case of two
      F90file.write( "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0" )
    elif len(thisTensor.indices) == 3: # in case of three
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == 0"
      F90file.write( cstr )
    elif len(thisTensor.indices) == 4: # in case of four
      cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      cstr += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
      F90file.write( cstr )
    elif len(thisTensor.indices) == 5: # in case of five
      cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
      cstr += "s_" + thisTensor.indices[4].name + ")"
      F90file.write( cstr )
    elif len(thisTensor.indices) == 6: # in case of six
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
      cstr += "s_" + thisTensor.indices[5].name + ")"
      F90file.write( cstr )
    elif len(thisTensor.indices) == 7: # in case of seven
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr += "s_" + thisTensor.indices[2].name + ") == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
      cstr += "IEOR(s_" + thisTensor.indices[5].name + ", " + "s_" + thisTensor.indices[5].name + ")"
      F90file.write( cstr )
    elif len(thisTensor.indices) == 8: # in case of eight (allowed only for 4-RDM)
      cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
      cstr +=      "IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + ")) == "
      cstr += "IEOR(IEOR(s_" + thisTensor.indices[4].name + ", " + "s_" + thisTensor.indices[5].name + "),"
      cstr +=      "IEOR(s_" + thisTensor.indices[6].name + ", " + "s_" + thisTensor.indices[7].name + "))"
      F90file.write( cstr )
    elif len(thisTensor.indices) != 0: # in case of scalar:
      raise TypeError, "Number of indices in tensor must be 0 - 8"
    if len(thisTensor.indices) != 0: F90file.write( " .and. & " + "\n" ) 
    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if len(thisTensor.indices) == 1:   # in case of unity
        F90file.write( "S_" + thisTensor.indices[0].name + " == 0" )
      elif len(thisTensor.indices) == 2: # in case of two
        F90file.write( "IEOR(" + "S_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + ") == 0" )
      elif len(thisTensor.indices) == 3: # in case of three
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == 0"
        F90file.write( cstr )
      elif len(thisTensor.indices) == 4: # in case of four
        cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        cstr += "IEOR(s_" + thisTensor.indices[2].name + ", s_" + thisTensor.indices[3].name + ")"
        F90file.write( cstr )
      elif len(thisTensor.indices) == 5: # in case of five
        cstr  = "IEOR(s_" + thisTensor.indices[0].name + ", s_" + thisTensor.indices[1].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + "),"
        cstr += "s_" + thisTensor.indices[4].name + ")"
        F90file.write( cstr )
      elif len(thisTensor.indices) == 6: # in case of six
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
        cstr += "s_" + thisTensor.indices[5].name + ")"
        F90file.write( cstr )
      elif len(thisTensor.indices) == 7: # in case of seven
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr += "s_" + thisTensor.indices[2].name + ") == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[3].name + ", " + "s_" + thisTensor.indices[4].name + "),"
        cstr += "IEOR(s_" + thisTensor.indices[5].name + ", " + "s_" + thisTensor.indices[5].name + ")"
        F90file.write( cstr )
      elif len(thisTensor.indices) == 8: # in case of eight (allowed only for 4-RDM)
        cstr  = "IEOR(IEOR(s_" + thisTensor.indices[0].name + ", " + "s_" + thisTensor.indices[1].name + "),"
        cstr +=      "IEOR(s_" + thisTensor.indices[2].name + ", " + "s_" + thisTensor.indices[3].name + ")) == "
        cstr += "IEOR(IEOR(s_" + thisTensor.indices[4].name + ", " + "s_" + thisTensor.indices[5].name + "),"
        cstr +=      "IEOR(s_" + thisTensor.indices[6].name + ", " + "s_" + thisTensor.indices[7].name + "))"
        F90file.write( cstr )
      else:
        raise TypeError, "Number of indices in tensor must be 1 - 8"

      if not numTensor == len(thisTerm.tensors) - 1:
        F90file.write( " .and. &" + "\n" )

    F90file.write( ") then" + "\n" )
    
    # Write down the body of the tensorial contraction ...
    # *CAUTION* Only active and virtual orbitals are considered separately for now 
    count = 0
    for Index in summed_Index:
      if Index.isExt == True: continue
      iname = "i_" + Index.name
      sname = "s_" + Index.name
      if   Index.indType[0][0] == 'core':
        otype = "I_C"
      elif Index.indType[0][0] == 'active':
        otype = "I_O"
      elif Index.indType[0][0] == 'virtual':
        otype = "I_V"
      F90file.write( "do " + iname + " = psym(I_BEGIN, " + otype + ", " + sname + "), psym(I_END, " + otype + ", " + sname+ ")" + "\n" )
      count += 1
    if count != loop_count:
      raise TypeError, "Algorithmic error occured .... "

    # As a first step of evalutation of the contraction, search the Kronecker deltas    
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': num_kdeltas += 1

    if num_kdeltas >= 1:
      count = num_kdeltas
      F90file.write( "if( &" + "\n" )
      for numTensor in range(len(thisTerm.tensors)):
        thisTensor = thisTerm.tensors[numTensor]
        if thisTensor.name == 'kdelta':
          cstr  = "s_" + thisTensor.indices[0].name + " == s_" + thisTensor.indices[1].name
          cstr += " .and. " + "i_" + thisTensor.indices[0].name + " == i_" + thisTensor.indices[1].name
          count -= 1
          if   count == 0: cstr += ") then"
          elif count >= 1: cstr += " .and. & "
          elif count >  0: raise TYpeError, "Algorithmic Error"
          F90file.write( cstr + "\n" )

    F90file.write( "" + "\n" )

    # First generate the left-hand side tensor
    # Make list of the indices which are of *NOT* external
    if   not isBareLTS and LTensor.name != 'X':
      LIndList = LTensor.indices
    elif not isBareLTS and LTensor.name == 'X':           # This can be varied .....
      LIndList = []
      for i in LTensor.indices:
        if i.isExt == False: LIndList.append(i) 
    else:                               # in case of BareAmP
      LIndList = []
      for n_index in range(len(LTensor.indices)):
        thisIndex = LTensor.indices[n_index]
        if n_index != extBare: LIndList.append(thisIndex)     

    if len(LTensor.indices) != 0:   # in case of Tensor
      nterm = LTensor.name + "_("
      for n_index in range(len(LIndList)):
        nterm += "s_" + LIndList[-(n_index+1)].name
        if not n_index == len(LIndList)-1: nterm += ", "
        else: nterm += ")"
      nterm += "%array("
      for n_index in range(len(LIndList)):
        nterm += "i_" + LIndList[-(n_index+1)].name
        if not n_index == len(LIndList)-1: nterm += ", "
        else: nterm += ")"
      F90file.write( nterm + " = " + nterm + " &" + "\n" )
    else: F90file.write( LTensor.name + "_ = " + LTensor.name + "_ &" + "\n" ) # in case of Scalar  

    # Generate all the tensors ...
    if thisTerm.numConstant > 0.0: F90file.write( "  + " + str(       thisTerm.numConstant) + "d+00" + " & " + "\n" )
    else:                          F90file.write( "  - " + str((-1.0)*thisTerm.numConstant) + "d+00" + " & " + "\n" )

    if len(thisTerm.constants) > 0: 
      for thisConst in thisTerm.constants: F90file.write( "  * " + thisConst + " & " + "\n" )

    # Generate each tensor except Kronecker deltas ....
    for numTensor in range(len(thisTerm.tensors)):
      thisTensor = thisTerm.tensors[numTensor]
      if thisTensor.name == 'kdelta': continue
     
      # Prepare the list of all indices, which are *NOT* of external, if the tensor is either of Eri or BareAmp
      if not thisTensor.name in [name_ofEri, name_ofBareAmp, name_ofRdm4, 'X']: # This can be different
        RIndList = thisTensor.indices
      elif thisTensor.name == name_ofEri:     # in case of ERI
        RIndList = []
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          if n_index != extEri: RIndList.append(thisIndex)
      elif thisTensor.name == name_ofBareAmp: # in case of BareAmp
        RIndList = []
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          if n_index != extBare: RIndList.append(thisIndex)
      elif thisTensor.name == name_ofRdm4:    # in case of 4-RDM
        RIndList = []
        dummy_count = 0
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          # Does it work correctly ??
          if (n_index != extRdm4[0]) and (n_index != extRdm4[1]): 
            RIndList.append(thisIndex)
            dummy_count += 1
        if(dummy_count > 6): raise TypeError, "4-RDMs must have at least 2 external indices" 
      elif thisTensor.name == 'X': # in case the intermediate
        RIndList = []
        for n_index in range(len(thisTensor.indices)):
          thisIndex = thisTensor.indices[n_index]
          if thisIndex.isExt == False: RIndList.append(thisIndex)

      nterm = "  * " + thisTensor.name + "_("
      for n_index in range(len(RIndList)):
        nterm += "s_" + RIndList[-(n_index+1)].name
        if not n_index == len(RIndList)-1: nterm += ", "
        else: nterm += ")"
      nterm += "%array("
      for n_index in range(len(RIndList)):
        nterm += "i_" + RIndList[-(n_index+1)].name
        if not n_index == len(RIndList)-1: nterm += ", "
        else: nterm += ")"
      if numTensor != len(thisTerm.tensors)-1: nterm += " & " 
      F90file.write( nterm + "\n" )
    
    F90file.write( "\n" )

    # Closing the iteration .... 
    if num_kdeltas >= 1: F90file.write( "end if ! Kdeltas\n" )                       # Kronecker deltas
    for i in range(loop_count): F90file.write( "end do ! Irrep loop\n" )             # Irrep loops 
    F90file.write( "end if ! Irrep Cond\n" )                                         # Irrep constraints
    for i in range(loop_count): F90file.write( "end do ! Orbital loop\n" )           # Orbital loops
    F90file.write( "! Loop_count     : " + str(loop_count)        + "\n" )
    F90file.write( "! External_count : " + str(len(extIndices))   + "\n" )
    F90file.write( "! Total_count    : " + str(len(summed_Index)) + "\n" )
#OK?#    if loop_count + len(extIndices) != len(summed_Index):
#OK?#      print "!!! WARNING: Sum of external and internal indices is not equal to that of all indices !!!\n"
#OK?#      print LTensor
#OK?#      for i in LTensor.indices: print i.isExt
#OK?#      print 
#OK?#      for t in thisTerm.tensors:
#OK?#        print t
#OK?#        for i in t.indices: print i.isExt
#OK?#        print
#OK?#      for i in summed_Index: print i.name, ", ", i.isSummed, ", ", i.isExt, ", " + i.indType[0][0]# *TEST*
#OK?#      raise TypeError, "Algorithmic Error, Somethig is wrong !!!"
    F90file.write( "\n" )
    F90file.write( "\n" )
  F90file.write( "! FEMTO END    **************************************************************\n\n" ) 
  F90file.write( printEnd + "\n" )

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
