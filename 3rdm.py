import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions
tg_core    = sqa.options.core_type
tg_active  = sqa.options.active_type
tg_virtual = sqa.options.virtual_type

a1 = sqa.index('a1', [tg_active], False)
b1 = sqa.index('b1', [tg_active], False)
a2 = sqa.index('a2', [tg_active], False)
b2 = sqa.index('b2', [tg_active], False)
a3 = sqa.index('a3', [tg_active], False)
b3 = sqa.index('b3', [tg_active], False)
a4 = sqa.index('a4', [tg_active], False)
b4 = sqa.index('b4', [tg_active], False)

a = sqa.index('a', [tg_active], True)

result = []
term = sqa.term(1.0, [], [sqa.sfExOp([a1,b1])] + [sqa.sfExOp([a2,b2])] + [sqa.sfExOp([a3,b3])])
result += sqa.normalOrder(term)

for t in result:
    t.contractDeltaFuncs()

sqa.combineTerms(result)

print "results"
for t in result:
    print t
