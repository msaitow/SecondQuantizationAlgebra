#    file:  sqaTest_removeCoreOpPairs.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Performs a test on the removeCoreOpPairs function.
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

# Define short names for type groups
ta = sqa.options.alpha_type
tb = sqa.options.beta_type
tc = sqa.options.core_type
tt = sqa.options.active_type
tv = sqa.options.virtual_type

# Define core indices
ac = [sqa.index('c%i' %i, [ta,tc], False) for i in range(10)]
bc = [sqa.index('c%i' %i, [tb,tc], False) for i in range(10)]

# Define active indices
at = [sqa.index('a%i' %i, [ta,tt], True) for i in range(10)]
bt = [sqa.index('a%i' %i, [tb,tt], True) for i in range(10)]

terms = []
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.desOp(ac[0])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]), sqa.tensor('r0', [at[0], at[1]], []), sqa.creOp(ac[1]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]), sqa.creOp(ac[1]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]), sqa.creOp(ac[0]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[1]), sqa.creOp(ac[0]), \
                                sqa.desOp(ac[1]), sqa.desOp(ac[0]), sqa.desOp(ac[1])]))

print ""
print "Testing removeCoreOpPairs function."
print ""
print "Core indices labeled with c."
print "Active indices labeled with a."

print ""
print "initial terms:"
for t in terms:
  print t

sqa.removeCoreOpPairs(terms)

output = ""
for t in terms:
  output += str(t) + "\n"

correctOutput = \
" (   1.00000) \n" + \
" (  -1.00000) cre(c0) \n" + \
" (   1.00000) cre(a0) r0(a0,a1) cre(c1) des(a1) \n" + \
" (   1.00000) cre(a0) des(a1) \n" + \
" (  -1.00000) cre(c0) cre(c0) des(c0) des(a1) \n" + \
" (  -1.00000) cre(a0) cre(c1) des(c1) des(c1) \n"

print ""
print "terms after core operator pair removal:"
print output

if output == correctOutput:
  print "Test passed!"
else:
  print "Test failed.  Correct output is:"
  print correctOutput
print ""
