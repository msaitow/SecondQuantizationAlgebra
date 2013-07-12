import secondQuantizationAlgebra as sqa

sqa.options.verbose = False
# definitions

tag_core    = sqa.options.core_type
tag_active  = sqa.options.active_type
tag_virtual = sqa.options.virtual_type

i = sqa.index('Ai', [tag_active],  False)
j = sqa.index('Aj', [tag_active],  False)
k = sqa.index('Ak', [tag_active],  False)
l = sqa.index('Al', [tag_active],  False)
a = sqa.index('Va', [tag_virtual], False)
b = sqa.index('Vb', [tag_virtual], False)
c = sqa.index('Vc', [tag_virtual], False)
d = sqa.index('Vd', [tag_virtual], False)

hsym   = sqa.symmetry((1,0), 1)
Dsym_a = sqa.symmetry((2,1, 0,3), 1)
Dsym_b = sqa.symmetry((0,3, 2,1), 1)
Dsym_c = sqa.symmetry((1,0, 3,2), 1)

E_iajb = [sqa.sfExOp([i, j, a, b])]
E_dlck = [sqa.sfExOp([b, a, j, i])]

print ""
print "Evaluation of the IC-MRCI Hamiltonian elements"
print "based on the spin-free generator"
print ""
print "!!!!! The operator space is L: <Psi|", E_iajb[0], E_dlck[0], "|Psi>: R !!!!!"
print ""

result = []
for tag_h1_p in [tag_active, tag_virtual]:
    for tag_h1_q in [tag_active, tag_virtual]:

        p = sqa.index('p', [tag_h1_p], True)
        q = sqa.index('q', [tag_h1_q], True)
        E_pq = [sqa.sfExOp([p, q])]
        h_pq = sqa.tensor('h', [p, q], [hsym])
        
        term1 = sqa.term( 1.0, [], E_iajb + [h_pq] + E_pq + E_dlck)
        
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
                V_pqrs = sqa.tensor('V', [p, q, r, s], [Dsym_a, Dsym_b, Dsym_c])
        
                term1 = sqa.term( 0.5, [], E_iajb + [V_pqrs] + E_pqrs + E_dlck)
        
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
print "* The one- and two-body parts....."
print ""

result = sqa.convert2Mulliken(result)

num = 0
for t in result:
    print num, t
    num += 1

print ""

Hdiag = sqa.tensor('Hdiag', [i, j, a, b], [])
sqa.factorize(Hdiag, result, 'diag_oovv_oovv', True, 'ham', 'V', 'T', 'int')

# # Calculation of the overlap elements .....
# term1 = sqa.term( 1.0, [], E_iajb + E_dlck)
        
# batch1 = sqa.normalOrder(term1)

# sqa.removeVirtOps_sf(batch1)

# for t in batch1:
#     t.contractDeltaFuncs_new()

# sqa.combineTerms(batch1)
# sqa.termChop(batch1)
                
# result = batch1

# sqa.combineTerms(result)
# sqa.termChop(result)                

# print ""
# print "* The overlap element....."
# print ""

# result = sqa.convert2Mulliken(result)

# num = 0
# for t in result:
#     print num, t
#     num += 1


        
