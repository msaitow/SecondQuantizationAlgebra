#    file:  sqaTest.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Performs some simple tests of the second quantization algebra program.
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



import secondQuantizationAlgebra as sqa
import time

# Define density matrix symmetries
d1sym = [sqa.symmetry((1,0),1)]
d2sym_aaaa = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((2,3,0,1),1)]
d2sym_abab = [sqa.symmetry((2,3,0,1),1)]
d3sym_aaaaaa = [sqa.symmetry((1,0,2,3,4,5),-1), sqa.symmetry((0,2,1,3,4,5),-1), sqa.symmetry((3,4,5,0,1,2),1)]
d3sym_aabaab = [sqa.symmetry((1,0,2,3,4,5),-1), sqa.symmetry((3,4,5,0,1,2),1)]

# Define Amplitude symmetries
a1sym = []
a2sym_aaaa = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2),-1)]
a2sym_abab = []
a2sym_sf = [sqa.symmetry((1,0,3,2),1)]

# Define integral symmetries
h1sym = d1sym
h2sym_aaaa = d2sym_aaaa
h2sym_abab = d2sym_abab
h2sym_sf = [sqa.symmetry((2,3,0,1),1), sqa.symmetry((1,0,3,2),1)]

# Define dummy indices
dummyInd = []
for i in range(20):
  dummyInd.append(sqa.index('i%i' %i, [], True))

# Define non-dummy indices
fixedInd = []
for i in range(20):
  fixedInd.append(sqa.index('j%i' %i, [], False))

# Test 1: Normal ordering of a(j1) a(j2) a+(j3) a+(j4)
print ""
print "Starting test 1"
print ""
startTime = time.time()
test1_tensors = []
test1_tensors.append(sqa.desOp(fixedInd[1]))
test1_tensors.append(sqa.desOp(fixedInd[2]))
test1_tensors.append(sqa.creOp(fixedInd[3]))
test1_tensors.append(sqa.creOp(fixedInd[4]))
test1_term = sqa.term(1.0, [], test1_tensors)
test1_output = sqa.normalOrder(test1_term)
sqa.combineTerms(test1_output)
test1_string_output = ''
for t in test1_output:
  test1_string_output += str(t) + '\n'
test1_correct_answer = " (  -1.00000) kdelta(j1,j3) kdelta(j2,j4) \n" + \
                       " (   1.00000) kdelta(j1,j4) kdelta(j2,j3) \n" + \
                       " (   1.00000) kdelta(j1,j3) cre(j4) des(j2) \n" + \
                       " (  -1.00000) kdelta(j1,j4) cre(j3) des(j2) \n" + \
                       " (  -1.00000) kdelta(j2,j3) cre(j4) des(j1) \n" + \
                       " (   1.00000) kdelta(j2,j4) cre(j3) des(j1) \n" + \
                       " (   1.00000) cre(j3) cre(j4) des(j1) des(j2) \n"
print "Test 1 output:"
print test1_string_output
if test1_string_output == test1_correct_answer:
  print "Test 1 passed",
else:
  print "Test 1 failed",
print " (%.3f seconds)\n" %(time.time()-startTime)

# Test 2: Commutator of all alpha 1e hamiltonian with all alpha 2 electron amplitude (without index types)
#         The result is transpose equated, so that symmetrization is required in the resulting Hamiltonian
print ""
print "Starting test 2"
print ""
startTime = time.time()
h1 = sqa.tensor('h1', dummyInd[0:2], h1sym)
a2 = sqa.tensor('a2', dummyInd[2:6], a2sym_aaaa)
h1_term = sqa.term(1.0, [], [h1, sqa.creOp(dummyInd[0]), sqa.desOp(dummyInd[1])])
a2_terms = []
a2_terms.append( sqa.term( 1.0, [], [a2, sqa.creOp(dummyInd[2]), sqa.creOp(dummyInd[3]), sqa.desOp(dummyInd[5]), sqa.desOp(dummyInd[4])]) )
a2_terms.append( sqa.term(-1.0, [], [a2, sqa.creOp(dummyInd[4]), sqa.creOp(dummyInd[5]), sqa.desOp(dummyInd[3]), sqa.desOp(dummyInd[2])]) )
test2_output = sqa.commutator(h1_term, a2_terms)
sqa.combine_transpose(test2_output)
test2_string_output = ''
for t in test2_output:
  test2_string_output += str(t) + '\n'
test2_correct_answer = \
                       " (   4.00000) a2(a,b,c,d) h1(a,e) cre(c) cre(d) des(b) des(e) \n" + \
                       " (  -4.00000) a2(a,b,c,d) h1(c,e) cre(d) cre(e) des(a) des(b) \n"
#                       " (   4.00000) a2(a,b,c,d) h1(a,e) cre(b) cre(e) des(c) des(d) \n" + \
#                       " (  -4.00000) a2(a,b,c,d) h1(c,e) cre(a) cre(b) des(d) des(e) \n"
print "Test 2 output:"
print test2_string_output
if test2_string_output == test2_correct_answer:
  print "Test 2 passed",
else:
  print "Test 2 failed",
print " (%.3f seconds)\n" %(time.time()-startTime)

# Test 3: Commutator of all alpha 2e hamiltonian with alpha/beta 2 electron amplitude, index types included
print ""
print "Starting test 3"
print ""
startTime = time.time()
test3_indices = []
test3_indices.append(sqa.index('i0', [['alpha'], ['core', 'active', 'virtual']], True))
test3_indices.append(sqa.index('i1', [['alpha'], ['core', 'active', 'virtual']], True))
test3_indices.append(sqa.index('i2', [['alpha'], ['core', 'active', 'virtual']], True))
test3_indices.append(sqa.index('i3', [['alpha'], ['core', 'active', 'virtual']], True))
test3_indices.append(sqa.index('i4', [['alpha'], ['active', 'virtual']], True))
test3_indices.append(sqa.index('i5', [[ 'beta'], ['active', 'virtual']], True))
test3_indices.append(sqa.index('i6', [['alpha'], ['core', 'active']], True))
test3_indices.append(sqa.index('i7', [[ 'beta'], ['core', 'active']], True))
h2_aaaa = sqa.tensor('h2_aaaa', test3_indices[0:4], h2sym_aaaa)
a2_abab = sqa.tensor('a2_abab', test3_indices[4:8], a2sym_abab)
h2_aaaa_term = sqa.term(1.0, [], [h2_aaaa, sqa.creOp(test3_indices[0]),
                                           sqa.creOp(test3_indices[1]),
                                           sqa.desOp(test3_indices[3]),
                                           sqa.desOp(test3_indices[2])])
a2_abab_terms = []
a2_abab_terms.append( sqa.term( 1.0, [], [a2_abab, sqa.creOp(test3_indices[4]),
                                                   sqa.creOp(test3_indices[5]),
                                                   sqa.desOp(test3_indices[7]),
                                                   sqa.desOp(test3_indices[6])]) )
a2_abab_terms.append( sqa.term(-1.0, [], [a2_abab, sqa.creOp(test3_indices[6]),
                                                   sqa.creOp(test3_indices[7]),
                                                   sqa.desOp(test3_indices[5]),
                                                   sqa.desOp(test3_indices[4])]) )
test3_output = sqa.commutator(h2_aaaa_term, a2_abab_terms)
test3_string_output = ''
for t in test3_output:
  test3_string_output += str(t) + '\n'
  test3_ind_list = []
  for ten in t.tensors:
    for ind in ten.indices:
      if not ( ind in test3_ind_list ):
        test3_ind_list.append(ind.copy())
        test3_string_output += 'index %s types = %s\n' %(ind.name, str(ind.indType))
test3_correct_answer = " (   2.00000) a2_abab(a,b,c,d) h2_aaaa(a,e,f,g) cre(b) cre(f) cre(g) des(c) des(d) des(e) \n" + \
                       "index a types = (('alpha',), ('active', 'virtual'))\n" + \
                       "index b types = (('beta',), ('active', 'virtual'))\n" + \
                       "index c types = (('alpha',), ('active', 'core'))\n" + \
                       "index d types = (('beta',), ('active', 'core'))\n" + \
                       "index e types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index f types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index g types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       " (   2.00000) a2_abab(a,b,c,d) h2_aaaa(a,e,f,g) cre(c) cre(d) cre(e) des(b) des(f) des(g) \n" + \
                       "index a types = (('alpha',), ('active', 'virtual'))\n" + \
                       "index b types = (('beta',), ('active', 'virtual'))\n" + \
                       "index c types = (('alpha',), ('active', 'core'))\n" + \
                       "index d types = (('beta',), ('active', 'core'))\n" + \
                       "index e types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index f types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index g types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       " (  -2.00000) a2_abab(a,b,c,d) h2_aaaa(c,e,f,g) cre(a) cre(b) cre(e) des(d) des(f) des(g) \n" + \
                       "index a types = (('alpha',), ('active', 'virtual'))\n" + \
                       "index b types = (('beta',), ('active', 'virtual'))\n" + \
                       "index c types = (('alpha',), ('active', 'core'))\n" + \
                       "index d types = (('beta',), ('active', 'core'))\n" + \
                       "index e types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index f types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index g types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       " (  -2.00000) a2_abab(a,b,c,d) h2_aaaa(c,e,f,g) cre(d) cre(f) cre(g) des(a) des(b) des(e) \n" + \
                       "index a types = (('alpha',), ('active', 'virtual'))\n" + \
                       "index b types = (('beta',), ('active', 'virtual'))\n" + \
                       "index c types = (('alpha',), ('active', 'core'))\n" + \
                       "index d types = (('beta',), ('active', 'core'))\n" + \
                       "index e types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index f types = (('alpha',), ('active', 'core', 'virtual'))\n" + \
                       "index g types = (('alpha',), ('active', 'core', 'virtual'))\n"
print "Test 3 output:"
print test3_string_output
if test3_string_output == test3_correct_answer:
  print "Test 3 passed",
else:
  print "Test 3 failed",
print " (%.3f seconds)\n" %(time.time()-startTime)

# Test 4: Normal ordering of the spin free operator string E(i0,i1) E(i2,i3,i4,i5) E(i6,i7,i8,i9)
print ""
print "Starting test 4"
print ""
startTime = time.time()
test4_indices = [sqa.index('i%i' %i) for i in range(10)]
test4_term = sqa.term(1.0, [], [sqa.sfExOp(test4_indices[0:2]), sqa.sfExOp(test4_indices[2:6]), sqa.sfExOp(test4_indices[6:10])])
test4_output = sqa.normalOrder(test4_term)
sqa.combineTerms(test4_output)
test4_string_output = ""
for t in test4_output:
  test4_string_output += str(t) + '\n'
test4_correct_answer = " (   1.00000) kdelta(i1,i2) kdelta(i4,i6) kdelta(i5,i7) E2(i0,i3,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i2) kdelta(i4,i7) kdelta(i5,i6) E2(i0,i3,i9,i8) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i4,i6) kdelta(i5,i7) E2(i0,i2,i9,i8) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i4,i7) kdelta(i5,i6) E2(i0,i2,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i2) kdelta(i4,i6) E3(i0,i3,i7,i8,i5,i9) \n" + \
                       " (   1.00000) kdelta(i1,i2) kdelta(i4,i7) E3(i0,i3,i6,i9,i5,i8) \n" + \
                       " (   1.00000) kdelta(i1,i2) kdelta(i5,i6) E3(i0,i3,i7,i4,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i2) kdelta(i5,i7) E3(i0,i3,i6,i4,i9,i8) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i4,i6) E3(i0,i2,i7,i5,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i4,i7) E3(i0,i2,i6,i5,i9,i8) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i5,i6) E3(i0,i2,i7,i8,i4,i9) \n" + \
                       " (   1.00000) kdelta(i1,i3) kdelta(i5,i7) E3(i0,i2,i6,i9,i4,i8) \n" + \
                       " (   1.00000) kdelta(i1,i6) kdelta(i4,i7) E3(i0,i2,i3,i8,i9,i5) \n" + \
                       " (   1.00000) kdelta(i1,i6) kdelta(i5,i7) E3(i0,i2,i3,i8,i4,i9) \n" + \
                       " (   1.00000) kdelta(i1,i7) kdelta(i4,i6) E3(i0,i2,i3,i9,i8,i5) \n" + \
                       " (   1.00000) kdelta(i1,i7) kdelta(i5,i6) E3(i0,i2,i3,i9,i4,i8) \n" + \
                       " (   1.00000) kdelta(i4,i6) kdelta(i5,i7) E3(i0,i2,i3,i1,i8,i9) \n" + \
                       " (   1.00000) kdelta(i4,i7) kdelta(i5,i6) E3(i0,i2,i3,i1,i9,i8) \n" + \
                       " (   1.00000) kdelta(i1,i2) E4(i0,i3,i6,i7,i4,i5,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i3) E4(i0,i2,i6,i7,i5,i4,i8,i9) \n" + \
                       " (   1.00000) kdelta(i1,i6) E4(i0,i2,i3,i7,i8,i4,i5,i9) \n" + \
                       " (   1.00000) kdelta(i1,i7) E4(i0,i2,i3,i6,i9,i4,i5,i8) \n" + \
                       " (   1.00000) kdelta(i4,i6) E4(i0,i2,i3,i7,i1,i8,i5,i9) \n" + \
                       " (   1.00000) kdelta(i4,i7) E4(i0,i2,i3,i6,i1,i9,i5,i8) \n" + \
                       " (   1.00000) kdelta(i5,i6) E4(i0,i2,i3,i7,i1,i4,i8,i9) \n" + \
                       " (   1.00000) kdelta(i5,i7) E4(i0,i2,i3,i6,i1,i4,i9,i8) \n" + \
                       " (   1.00000) E5(i0,i2,i3,i6,i7,i1,i4,i5,i8,i9) \n"
print "Test 4 output:"
print test4_string_output
if test4_string_output == test4_correct_answer:
  print "Test 4 passed",
else:
  print "Test 4 failed",
print " (%.3f seconds)\n" %(time.time()-startTime)

# Test 5: Commutator of 2e spin free hamiltonian with 2e spin free amplitudes
print ""
print "Starting test 5"
print ""
startTime = time.time()
test5_indices = []
test5_indices.append( sqa.index('i0', [['core', 'active', 'virtual']], True) )
test5_indices.append( sqa.index('i1', [['core', 'active', 'virtual']], True) )
test5_indices.append( sqa.index('i2', [['core', 'active', 'virtual']], True) )
test5_indices.append( sqa.index('i3', [['core', 'active', 'virtual']], True) )
test5_indices.append( sqa.index('i4', [['active', 'virtual']], True) )
test5_indices.append( sqa.index('i5', [['active', 'virtual']], True) )
test5_indices.append( sqa.index('i6', [['core', 'active']], True) )
test5_indices.append( sqa.index('i7', [['core', 'active']], True) )
test5_h2 = sqa.tensor('h2', test5_indices[0:4], h2sym_sf)
test5_hterm = sqa.term(1.0, [], [test5_h2, sqa.sfExOp(test5_indices[0:4])])
test5_a2 = sqa.tensor('a2', test5_indices[4:8], a2sym_sf)
test5_aterms = []
test5_aterms.append( sqa.term( 1.0, [], [test5_a2, sqa.sfExOp(test5_indices[4:6] + test5_indices[6:8])]) )
test5_aterms.append( sqa.term(-1.0, [], [test5_a2, sqa.sfExOp(test5_indices[6:8] + test5_indices[4:6])]) )
test5_output = sqa.commutator(test5_hterm, test5_aterms)
test5_string_output = ""
for t in test5_output:
  test5_string_output += str(t) + '\n'
  test5_ind_list = []
  for ten in t.tensors:
    for ind in ten.indices:
      if not ( ind in test5_ind_list ):
        test5_ind_list.append(ind.copy())
        test5_string_output += 'index %s types = %s\n' %(ind.name, str(ind.indType))
test5_correct_answer = " (   2.00000) a2(a,b,c,d) h2(a,b,e,f) E2(c,d,e,f) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       " (   2.00000) a2(a,b,c,d) h2(a,b,e,f) E2(e,f,c,d) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       " (  -2.00000) a2(a,b,c,d) h2(c,d,e,f) E2(a,b,e,f) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       " (  -2.00000) a2(a,b,c,d) h2(c,d,e,f) E2(e,f,a,b) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       " (   4.00000) a2(a,b,c,d) h2(a,e,f,g) E3(b,f,g,d,c,e) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       "index g types = (('active', 'core', 'virtual'),)\n" + \
                       " (   4.00000) a2(a,b,c,d) h2(a,e,f,g) E3(c,d,e,f,b,g) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       "index g types = (('active', 'core', 'virtual'),)\n" + \
                       " (  -4.00000) a2(a,b,c,d) h2(c,e,f,g) E3(a,b,e,f,d,g) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       "index g types = (('active', 'core', 'virtual'),)\n" + \
                       " (  -4.00000) a2(a,b,c,d) h2(c,e,f,g) E3(d,f,g,b,a,e) \n" + \
                       "index a types = (('active', 'virtual'),)\n" + \
                       "index b types = (('active', 'virtual'),)\n" + \
                       "index c types = (('active', 'core'),)\n" + \
                       "index d types = (('active', 'core'),)\n" + \
                       "index e types = (('active', 'core', 'virtual'),)\n" + \
                       "index f types = (('active', 'core', 'virtual'),)\n" + \
                       "index g types = (('active', 'core', 'virtual'),)\n"
print "Test 5 output:"
print test5_string_output
if test5_string_output == test5_correct_answer:
  print "Test 5 passed",
else:
  print "Test 5 failed",
print " (%.3f seconds)\n" %(time.time()-startTime)
