import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions

tag_core    = sqa.options.core_type
tag_active  = sqa.options.active_type
tag_virtual = sqa.options.virtual_type

i = sqa.index('Ai', [tag_active],  False)
j = sqa.index('Aj', [tag_active],  False)
k = sqa.index('Ak', [tag_active],  False)
a = sqa.index('Va', [tag_virtual], False)

hsym   = sqa.symmetry((1,0), 1)
Dsym_a = sqa.symmetry((2,1, 0,3), 1)
Dsym_b = sqa.symmetry((0,3, 2,1), 1)
Dsym_c = sqa.symmetry((1,0, 3,2), 1)
Tsym   = sqa.symmetry((0,1, 2,3), 1)

E_ijka = [sqa.sfExOp([i, j, k, a])]

print ""
print "Evaluation of the IC-MRCI Hamiltonian elements"
print "based on the spin-free generator"
print ""
print "!!!!! The operator space is L: <Psi|", E_ijka , "|Psi>: R !!!!!"
print "@@@ SIGMA(Ai, Aj, Ak, Va) <--  H(Ai, Aj, Ak, Va;0) T(0)"
print ""

result = []
for tag_h1_p in [tag_active, tag_virtual]:
    for tag_h1_q in [tag_active, tag_virtual]:

        p = sqa.index('p', [tag_h1_p], True)
        q = sqa.index('q', [tag_h1_q], True)
        E_pq = [sqa.sfExOp([p, q])]
        h_pq = sqa.tensor('h', [p, q], [hsym])
        
        term1 = sqa.term( 1.0, ["T0"], E_ijka + [h_pq] + E_pq)
        
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
        
                term1 = sqa.term( 0.5, ["T0"], E_ijka + [V_pqrs] + E_pqrs)
        
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
print "* The one- and two-body part....."
print ""

result = sqa.convert2Mulliken(result, 'T2')

num = 0
for t in result:
    print num, t
    num += 1

print ""

sigma = sqa.tensor('S2', [i, j, k, a], [])
sqa.factorize(sigma, result, 'sigma_ooov_G', True, 'sig', 'V2', 'T2', 'int')
