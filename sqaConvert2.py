#    file:  sqaConvert2.py
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

from sqaTensor import tensor, kroneckerDelta, creOp, desOp, sfExOp
from sqaTerm import term, sortOps
from sqaIndex import index
from sqaMisc import makeTuples, allDifferent, makePermutations, get_num_perms
from sqaOptions import options


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def CPloop(Spaces, Index, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(Index, index):
    raise TypeError, "Index should be a sqaIndex object"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  Itype = ""
  if   Index.indType[0][0] == 'core':
    Itype = "I_C"
  elif Index.indType[0][0] == 'active':
    Itype = "I_O"
  elif Index.indType[0][0] == 'virtual':
    Itype = "I_V"

  File.write(\
"""  %sfor(int s%s = 0;s%s < nir;++s%s){
  %sfor(int i%s = symblockinfo.psym()(s%s,%s,I_BEGIN);i%s <= symblockinfo.psym()(s%s,%s,I_END);++i%s){
""" % (Spaces, Index.name, Index.name, Index.name, \
       Spaces, Index.name, Index.name, Itype,      Index.name, Index.name, Itype, Index.name))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def ReadBare(Spaces, T2, extBare, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(T2, tensor):
    raise TypeError, "T2 should be a tensor object"

  if not type(extBare) == type(1):
    raise TypeError, "extBare should be an integer"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  File.write("  %s%sb = %s.get_amp2(i%s);\n" % (Spaces, T2.name, T2.name, T2.indices[extBare].name))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def ReadRetval(Spaces, T2, extBare, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(T2, tensor):
    raise TypeError, "T2 should be a tensor object"

  if not type(extBare) == type(1):
    raise TypeError, "extBare should be an integer"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  File.write("  %s%sb = orz::DTensor(retval.namps_iamp()[i%s]);\n" % (Spaces, T2.name, T2.indices[extBare].name))



#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def ReadERI(Spaces, V2, extEri, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(V2, tensor):
    raise TypeError, "V2 should be a tensor object"

  if not type(extEri) == type(1):
    raise TypeError, "extEri should be an integer"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  i = V2.indices[extEri].name

  File.write(\
"""  %s// Load ERIs from somewhere, e.g. disk, GA, etc..
  %s%s <<= 0.0;
  %sshared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(i%s);
  %sfor(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
    %s// Load a signle record of integals
    %sconst int &imo2 = loadbuf_ptr->i0;
    %sconst int &imo3 = loadbuf_ptr->i1;
    %sconst int &imo4 = loadbuf_ptr->i2;
    %sconst double &v = loadbuf_ptr->v;
    %s%s_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
  %s}
  %sconst orz::DTensor %s_sym = orz::ct::sympack_int2(symblockinfo, i%s, s%s, %s); // V2=(IR-COV index)
""" % (Spaces, \
       Spaces, V2.name, \
       Spaces, i, \
       Spaces, \
       Spaces, \
       Spaces, \
       Spaces, \
       Spaces, \
       Spaces, \
       Spaces, V2.name, \
       Spaces, \
       Spaces, V2.name, i, i, V2.name))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def ReadD4(Spaces, D4, extRdm4, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(D4, tensor):
    raise TypeError, "D4 should be a tensor object"

  if not type(extRdm4) == type([]):
    raise TypeError, "extRdm4 should be a list"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  i1 = D4.indices[extRdm4[0]].name
  i2 = D4.indices[extRdm4[1]].name

  File.write(\
"""  %s// Load D4 from disk, or GA .... 
  %sint imoi = amo2imo[i%s];
  %sint imoj = amo2imo[i%s];

  %sorz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
  %s                                   orz::Slice(),            orz::Slice(),
  %s                                   orz::Slice(),            orz::Slice(),
  %s                                   orz::Slice(),            orz::Slice()).copy();
  %srdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, i%s, s%s, i%s, s%s, rdm4_ij_sliced);
  %sFC_FUNC(g_if_set_d4,G_IF_SET_D4)(s%s, s%s, i%s, i%s, rdm4_sym.cptr(), nir, nsym, psym);
""" % (Spaces,     \
       Spaces, i1, \
       Spaces, i2, \
                   \
       Spaces,     \
       Spaces,     \
       Spaces,     \
       Spaces,     \
       Spaces, i1, i1, i2, i2, \
       Spaces, i1, i2, i1, i2))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def AccBare(Spaces, T2, extBare, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  if not isinstance(T2, tensor):
    raise TypeError, "T2 should be a tensor object"

  if not type(extBare) == type(1):
    raise TypeError, "extBare should be an integer"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  File.write("  %sretval.acc_amp2(i%s, %sb);\n" % (Spaces, T2.indices[extBare].name, T2.name))


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def LoopEnd(Spaces, File):

  if not type(Spaces) == type("   "):
    raise TypeError, "Spaces should be a string"

  a = open("a.out", "w")
  if not type(File) == type(a):
    raise TypeError, "File should be a file object"

  File.write("  %s}\n" % Spaces)
  File.write("  %s}\n" % Spaces)



