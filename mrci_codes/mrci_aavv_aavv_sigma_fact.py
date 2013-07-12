import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions

tag_core    = sqa.options.core_type
tag_active  = sqa.options.active_type
tag_virtual = sqa.options.virtual_type

i = sqa.index('Ai', [tag_active],  False)
j = sqa.index('Aj', [tag_active],  False)
a = sqa.index('Va', [tag_virtual], False)
b = sqa.index('Vb', [tag_virtual], False)

k = sqa.index('Ak', [tag_active],  True)
l = sqa.index('Al', [tag_active],  True)
c = sqa.index('Vc', [tag_virtual], True)
d = sqa.index('Vd', [tag_virtual], True)

hsym   = sqa.symmetry((1,0), 1)
Dsym_a = sqa.symmetry((2,1, 0,3), 1)
Dsym_b = sqa.symmetry((0,3, 2,1), 1)
Dsym_c = sqa.symmetry((1,0, 3,2), 1)
Tsym   = sqa.symmetry((0,1, 2,3), 1)

E_ijab = [sqa.sfExOp([i, j, a, b])]
E_dclk = [sqa.sfExOp([d, c, l, k])]

T_dclk = [sqa.tensor('T2', [k, l, c, d], [Tsym])]

print ""
print "Evaluation of the IC-MRCI Hamiltonian elements"
print "based on the spin-free generator"
print ""
print "!!!!! The operator space is L: <Psi|", E_ijab[0], E_dclk[0], "|Psi>: R !!!!!"
print "@@@ SIGMA(Ai, Aj, Va, Vb) <-- sum_{klcd} H(Ai, Aj, Va, Vb; Vd, Vc, Al, Ak) T(Ak, Al, Vc, Vd)"
print ""

result = []
for tag_h1_p in [tag_active, tag_virtual]:
  for tag_h1_q in [tag_active, tag_virtual]:

    p = sqa.index('p', [tag_h1_p], True)
    q = sqa.index('q', [tag_h1_q], True)
    E_pq = [sqa.sfExOp([p, q])]
    h_pq = sqa.tensor('h', [p, q], [hsym])
        
    term1 = sqa.term( 1.0, [], E_ijab + [h_pq] + E_pq + T_dclk + E_dclk)
        
    batch1 = sqa.normalOrder(term1)

    sqa.removeVirtOps_sf(batch1)

    for t in batch1:
      t.contractDeltaFuncs_new()

    sqa.combineTerms(batch1)
    sqa.termChop(batch1)
                
    result += batch1

    sqa.combineTerms(result)
    sqa.termChop(result)                

for tag_h2_p in [tag_active, tag_virtual]:
  for tag_h2_q in [tag_active, tag_virtual]:
    for tag_h2_r in [tag_active, tag_virtual]:
      for tag_h2_s in [tag_active, tag_virtual]: 
        
        p = sqa.index('p', [tag_h2_p], True)
        q = sqa.index('q', [tag_h2_q], True)
        r = sqa.index('r', [tag_h2_r], True)
        s = sqa.index('s', [tag_h2_s], True)
        E_pqrs = [sqa.sfExOp([p, q, r, s])]
        V_pqrs = sqa.tensor('V2', [p, q, r, s], [Dsym_a, Dsym_b, Dsym_c])
        
        term1 = sqa.term( 0.5, [], E_ijab + [V_pqrs] + E_pqrs + T_dclk + E_dclk)
        
        batch1 = sqa.normalOrder(term1)
        
        sqa.removeVirtOps_sf(batch1)
        
        for t in batch1:
          t.contractDeltaFuncs_new()
        
        sqa.combineTerms(batch1)
        sqa.termChop(batch1)
        
        result += batch1
                    
        sqa.combineTerms(result)
        sqa.termChop(result)                


print ""
print "* The one- and two-body parts ....."
print ""

result = sqa.convert2Mulliken(result, 'T2')

num = 0
for t in result:
  print num, t
  num += 1

print ""

sigma = sqa.tensor('S2', [i, j, a, b], [])
sqa.factorize(sigma, result, 'sigma_oovv_oovv', True, 'sig', 'V2', 'T2', 'int')


        
