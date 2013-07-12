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

E_ijka = [sqa.sfExOp([i, j, k, a])]
E_bnml = [sqa.sfExOp([a, k, j, i])]

print ""
print "Evaluation of the IC-MRCI Hamiltonian elements based on"
print "the spin-free generator and the commutator technique"
print ""
print "!!!!! The operator space is L: <Psi|", E_ijka[0],", ",E_bnml[0], "|Psi>: R !!!!!"
print "@@@ SIGMA(Ai, Aj, Ak, Va) <-- H(Ai, Aj, Ak, Va; Vb, Am, Al, Ak) T(Ak, Al, Am, Vb)"
print ""

result = []
# One-body part ....
for tag_h1_p in [tag_active, tag_virtual]:
    for tag_h1_q in [tag_active, tag_virtual]:

        p = sqa.index('p', [tag_h1_p], True)
        q = sqa.index('q', [tag_h1_q], True)
        E_pq = [sqa.sfExOp([p, q])]
        h_pq = sqa.tensor('h', [p, q], [hsym])
        
        term1 = sqa.term( 1.0, [], E_ijka + [h_pq] + E_pq + E_bnml)
        
        batch1 = sqa.normalOrder(term1)

        sqa.removeVirtOps_sf(batch1)

        for t in batch1:
            t.contractDeltaFuncs_new()

        sqa.combineTerms(batch1)
        sqa.termChop(batch1)
                
        result += batch1

        term2 = sqa.term(-1.0, [], E_ijka + E_bnml + [h_pq] + E_pq)
    
        batch2 = sqa.normalOrder(term2)

        for t in batch2:
            t.contractDeltaFuncs_new()
                
        sqa.removeVirtOps_sf(batch2)
        sqa.combineTerms(batch2)
        sqa.termChop(batch2)
                
        result += batch2
    
        sqa.combineTerms(result)
        sqa.termChop(result)

# Two-body part ....
for tag_h2_p in [tag_active, tag_virtual]:
    for tag_h2_q in [tag_active, tag_virtual]:
        for tag_h2_r in [tag_active, tag_virtual]:
            for tag_h2_s in [tag_active, tag_virtual]: 
        
                p = sqa.index('p', [tag_h2_p], True)
                q = sqa.index('q', [tag_h2_q], True)
                r = sqa.index('r', [tag_h2_r], True)
                s = sqa.index('s', [tag_h2_s], True)
                E_pqrs = [sqa.sfExOp([p, q, r, s])]
                V_pqrs = sqa.tensor('V', [p, q, r, s], [Dsym_a, Dsym_b, Dsym_c])
        
                term1 = sqa.term( 0.5, [], E_ijka + [V_pqrs] +  E_pqrs  + E_bnml)
                term2 = sqa.term(-0.5, [], E_ijka +  E_bnml  + [V_pqrs] + E_pqrs)
        
                batch1 = sqa.normalOrder(term1)
        
                sqa.removeVirtOps_sf(batch1)
        
                for t in batch1:
                    t.contractDeltaFuncs_new()
        
                sqa.combineTerms(batch1)
                sqa.termChop(batch1)
        
                result += batch1
        
                batch2 = sqa.normalOrder(term2)
        
                for t in batch2:
                    t.contractDeltaFuncs_new()
                    
                sqa.removeVirtOps_sf(batch2)

                sqa.combineTerms(batch2)
                sqa.termChop(batch2)
        
                result += batch2

                sqa.combineTerms(result)
                sqa.termChop(result)                

term3 = sqa.term( 1.0, ["Ecas"], E_ijka + E_bnml)
result2 = sqa.normalOrder(term3)
sqa.removeVirtOps_sf(result2)

for t in result2:
    t.contractDeltaFuncs_new()

sqa.combineTerms(result2)
sqa.termChop(result2)

result += result2

print ""
print "* The zero, one and two-body parts ....."
print ""

result = sqa.convert2Mulliken(result)

num = 0
for t in result:
    print num, t
    num += 1

Hdiag = sqa.tensor('Hdiag', [i, j, k, a], [])
sqa.convertF90_ERI(Hdiag, result, True, 'ham', 'V', 'T')



        
