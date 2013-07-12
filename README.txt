#    file:  README.txt
#  author:  Eric Neuscamman
#    date:  March 31, 2009
# summary:  Provides instruction for using the operator algebra program.
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



This program is intended to automate many of the tedious manipulations
that one encounters when working in second quantization.
To use it, be sure that your Python path includes this directory.
Then, in the python interpreter or in a python script, use a command
such as:

  >> import secondQuantizationAlgebra as sqa

Here we've given the module a conveniently short name.  Note that in
this tutorial, the >> symbol indicates the prompt in the Python
interpreter.  What follows is a set of examples intended to introduce
the capabilities of the module.


#####################################################################
# 1. Normal ordering
#####################################################################

One of the most common tasks in second quantization is to convert
a set of creation and destruction operators in to normal order.
Before we do, we need to define the operators, and to do that we
need to define some indices.  Let's use p, q, r, and s.

>> p = sqa.index('p')
>> q = sqa.index('q')
>> r = sqa.index('r')
>> s = sqa.index('s')

Now that we have our indices, let's define two creation and two
destruction operators with them.

>> p_cre = sqa.creOp(p)
>> q_cre = sqa.creOp(q)
>> r_des = sqa.desOp(r)
>> s_des = sqa.desOp(s)

Note that the creOp and desOp classes represent the creation and
destruction operators, and that each takes an index as its argument.
We are now ready to define a term, which will represent a non-
normal-ordered multiplication of the four operators.

>> unordered = sqa.term(1.0, [], [r_des, s_des, p_cre, q_cre])

In order to build the term, we have grouped the operators into
a list and passed them as the third argument.  The first argument is
the numerical constant, which we just set to 1 for now.
The second argument is a list of any constant variables multiplying
the term, given as as list of strings.  In this case we have no
constant variables, so we supply an empty list.
Before we normal order the operators, let's see what the term
looks like by printing it.

>> print unordered

You should see the following output:

 (   1.00000) des(r) des(s) cre(p) cre(q)

Now let's put the operators in normal order.

>> ordered_terms = sqa.normalOrder(unordered)

The result of normal ordering, which consists of more than one term,
is stored as a list of term objects.  To see the result, just print
out the terms one by one.

>> for t in ordered_terms:
>>   print t

The result that prints should look like this.

 (   1.00000) cre(p) cre(q) des(r) des(s) 
 (   1.00000) kdelta(r,p) cre(q) des(s) 
 (  -1.00000) kdelta(r,q) cre(p) des(s) 
 (  -1.00000) kdelta(s,p) cre(q) des(r) 
 (   1.00000) kdelta(s,q) cre(p) des(r) 
 (  -1.00000) kdelta(r,p) kdelta(s,q) 
 (   1.00000) kdelta(r,q) kdelta(s,p)

As expected, there are two double contractions, four single
contractions, and the uncontracted set of normal ordered operators.
kdelta is an abbreviation for the Kronecker delta function, which in
this program is represented by the kroneckerDelta class.


#####################################################################
# 2. Canonical Form and Combining Like Terms
#####################################################################

An important aspect of this automated algebra program is its ability
to combine like terms.  Term comparison and combination is carried out
by converting terms to unique canonical forms, which are based on a
lexicographical ordering of tensors and indices with some special rules
for when a tensor occurs multiple times in a term.  This canonical form
is useful because terms with dummy indices and high-symmetry tensors
can be written in many ways, making it possible for equivalent terms to
appear very different.

The first obstacle to term comparison is that dummy indices can take
on any name.  By default, indices in the program are created as non-
dummy indices, i.e. their name is fixed.  If you wish to create a 
dummy index instead, additional arguments are used as shown.

>> my_dummy = sqa.index('d', [], True)

As before, the first argument is the index's name, d in this case.
The second argument is a list of index type groups which will be
explained later and for now is left empty.  The third argument is a
bool that controls whether the index is a dummy index (True for dummy,
False for non-dummy).

The second obstacle to term comparison is that tensors often have
symmetries that allow one to rearrange their indices.  These symmetries
are represented in the program by the symmetry class.  For example:

>> my_symmetry = sqa.symmetry((1,0),1)

This object represents a symmetry in a 2-index tensor.  The first
argument is the index permutation corresponding to the symmetry (in
this case swapping the first and second indices), while the second
argument is the numerical factor that is applied when one executes the
permutation.  Thus this symmetry represents that of a symmetric matrix.
We may swap the two indices and the value of the element is the same.
A more complicated symmetry is:

>> my_symmetry_2 = sqa.symmetry((0,1,3,2),-1)

This symmetry corresponds to a four index tensor in which swapping the
last two indices changes the sign of the element.  An example of such
a tensor is the all-alpha-spin spin-orbital 2-body reduced density
matrix (RDM).  Let's use this density matrix as an example of how to
construct a tensor with symmetry.  First we'll need its four indices,
for which we'll create a list of non-dummy indices.

>> rdm_index_list = [sqa.index('p%i' %j, [], False) for j in range(4)]

This command uses a handy python syntax for constructing lists to create
a list of non-dummy indices with names 'p0', 'p1', 'p2', and 'p3'.
Now lets create the symmetry objects to describe the RDM's symmetry.
Even though the 2-body all-alpha-spin RDM has 8-fold symmetry, we only
need the following two symmetry objects to create it:

>> rdm_sym_1 = sqa.symmetry((0,1,3,2),-1)
>> rdm_sym_2 = sqa.symmetry((2,3,0,1),1)

The reason we don't need the others is that they can all be created by
repeated applications of the above two, and the program will do exactly
that when exploring the possible arrangements of the tensor's indices.
Now we are ready to build our first tensor.

>> rdm_2_aaaa = sqa.tensor('d2aaaa', rdm_index_list, [rdm_sym_1, rdm_sym_2])

This creates a tensor with name 'd2aaaa', the indices in the list we created
earlier, and the two symmetries.  Note that both indices and symmetries are
passed as lists.  To see what we've made lets create a term that is just
this tensor and print it out.

>> rdm_term = sqa.term(1.0, [], [rdm_2_aaaa])
>> print rdm_term

You should see the output:
 (   1.00000) d2aaaa(p0,p1,p2,p3)

When written on paper, it is common to differentiate between the RDM's top
and bottom indices.  In this program, when such a differentiation makes sense,
the first half of the indices are the 'top' indices and the second half are
the 'bottom' indices.  Indices pair vertically, so in this example p0 and p2
 are a pair and p1 and p3 are a pair.  Care must be taken when combining
tensors and creation/destruction operators, as these operators use a different
order for pairing.  The RDM above results from the expectation value of the
following string of normal-ordered creation and destruction operators:

 (   1.00000) cre(p0) cre(p1) des(p3) des(p2)

Notice that the outermost operators form a pair, as do the innermost.  This
convention of index pairing in tensors and operators is arbitrary, but
certain aspects of the program rely on this assumption.

So what about canonical forms and term combination?  Now that we know how to
define dummy indices and tensors with symmetry, we can perform an example.
Lets consider a term involving the all-alpha-spin 2-body RDM and the
all-alpha-spin 2-body portion of the Hamiltonian.  Note that the symmetries
of the two tensors involved happen to be the same.

>> p0 = sqa.index('p0', [], True)   # dummy
>> p1 = sqa.index('p1', [], True)   # dummy
>> p2 = sqa.index('p2', [], True)   # dummy
>> p3 = sqa.index('p3', [], True)   # dummy
>> q0 = sqa.index('q0', [], False)  # non-dummy
>> q1 = sqa.index('q1', [], False)  # non-dummy
>> p0_cre = sqa.creOp(p0)
>> p1_cre = sqa.creOp(p1)
>> p2_des = sqa.desOp(p2)
>> p3_des = sqa.desOp(p3)
>> sym_1 = sqa.symmetry((0,1,3,2),-1)
>> sym_2 = sqa.symmetry((2,3,0,1),1)
>> rdm_2_aaaa = sqa.tensor('d2aaaa', [p0, q0, p1, q1], [sym_1, sym_2])
>> h_2_aaaa   = sqa.tensor('h2aaaa', [p0, p1, p2, p3], [sym_1, sym_2])
>> my_term = sqa.term(1.0, [], [rdm_2_aaaa, h_2_aaaa, p0_cre, p1_cre, p3_des, p2_des])
>> print my_term

The output should be:
 (   1.00000) d2aaaa(p0,q0,p1,q1) h2aaaa(p0,p1,p2,p3) cre(p0) cre(p1) des(p3) des(p2)

Here the dummy indices p0, p1, p2, and p3 are summed over, while the non-dummy
indices q0 and q1 are not.  The 2-body portion of the Hamiltonian has been
written as the sum of the integrals (h2aaaa) and corresponding creation and
destruction operators.  Note that the order of destruction operators is
opposite the order of the last two indices in h2aaaa, as per our convention.
The presence of the RDM acts as a weighting on the integrals in the sum, which
may not have any physical meaning but acts as a good example for us.
In addition to being able to rename the dummy indices, the tensors' indices
can be reordered in many ways using their symmetries.  Further, commutation
relations give the creation and destruction operators the freedom to reorder
as well.  Thus it may be very difficult to determine whether this term is the
same as another term.  This is where canonical form is helpful.  It will
re-write the term in an equivalent but unique way.  Lets try it out.

>> my_term.makeCanonical()
>> print my_term

The output should be:
 (  -1.00000) d2aaaa(q0,a,q1,b) h2aaaa(a,b,c,d) cre(a) cre(b) des(c) des(d)

Notice here that the dummy indices p0, p1, p2, and p3 have been renamed as
a, b, c, and d.  The non-dummy indices q0 and q1 have been left alone, because
they are not summed over and so changing their names would change the term.
The indices in the d2aaaa tensor have also been rearranged through the
permutation (1,0,3,2).  Even though this permutation symmetry was not supplied
specifically, it can be obtained by repeated applications of the two that we
did supply.  This permutation should not have changed the term's sign, however,
so what did?  After renaming the indices p2 and p3 as c and d, the destruction
operators were not in lexicographical order.  The program therefore swapped
them, which produced the sign.  Finally, why did the program put q0 ahead of
a and q1 ahead of b, when both these arrangements are not alphabetical?  The
reason is that before arranging alphabetically, it orders non-dummy indices
ahead of dummy indices.

Now that we've seen how a term is converted to canonical form, let's put this
to practical use by combining two terms that look different but are actually
the same.  First we'll create a copy of the term we just used, but with
different dummy index names and index orderings.

>> i0 = sqa.index('i0', [], True)
>> i1 = sqa.index('i1', [], True)
>> i2 = sqa.index('i2', [], True)
>> i3 = sqa.index('i3', [], True)
>> i0_cre = sqa.creOp(i0)
>> i1_cre = sqa.creOp(i1)
>> i2_des = sqa.desOp(i2)
>> i3_des = sqa.desOp(i3)
>> rdm_2_aaaa_other = sqa.tensor('d2aaaa', [q0, i0, q1, i1], [sym_1, sym_2])
>> h_2_aaaa_other   = sqa.tensor('h2aaaa', [i2, i3, i0, i1], [sym_1, sym_2])
>> my_other_term = sqa.term(1.0, [], [rdm_2_aaaa_other, h_2_aaaa_other, i0_cre, i1_cre, i3_des, i2_des])

Now lets create a list containing our two terms and print them out.

>> term_list = [my_term, my_other_term]
>> for t in term_list:
>>   print t

The result should be:
 (  -1.00000) d2aaaa(q0,a,q1,b) h2aaaa(a,b,c,d) cre(a) cre(b) des(c) des(d) 
 (   1.00000) d2aaaa(q0,i0,q1,i1) h2aaaa(i2,i3,i0,i1) cre(i0) cre(i1) des(i3) des(i2) 

These two terms certainly don't look the same at first glance, but we have
constructed them so that they are.  To attempt to combine them, we use the
following command.

>> sqa.combineTerms(term_list)

The combineTerms function will print out its progress in converting terms to
canonical form if the variable sqa.options.verbose is set to True.
This is the default, so the printout should show as:

Combining like terms:
Converting 2 terms to canonical form...
     0   (  -1.00000) d2aaaa(q0,a,q1,b) h2aaaa(a,b,c,d) cre(a) cre(b) des(c) des(d) 
     1   (   1.00000) d2aaaa(q0,i0,q1,i1) h2aaaa(i2,i3,i0,i1) cre(i0) cre(i1) des(i3) des(i2) 
Finished combining terms in 0.175 seconds

We can check the result by printing out the terms in our list again.

>> for t in term_list:
>>   print t

The printed result is:
 (  -2.00000) d2aaaa(q0,a,q1,b) h2aaaa(a,b,c,d) cre(a) cre(b) des(c) des(d) 

Here we see that the terms were indeed the same and have been added together
to form a single term.  Thus to combine any like terms, we simply put them
in a list together and use the combineTerms function.  This is why a
canonical form is so useful.


#####################################################################
# 3. Index Type Groups
#####################################################################

When normal ordering terms and contracting delta functions, there are many
cases in which a delta function may be identified as zero.  For example, the
delta function between indices relating to spin-orbitals of different spin is
zero.  We can teach the program to take advantage of these cases by using
index type groups.  An index type group is a tuple of strings that define what
range the index refers to.  For example, the type group ('alpha', 'beta')
implies that an index ranges over both alpha and beta orbitals.  Another
example could be ('core', 'active', 'virtual'), implying the index ranges over
the core, active, and virtual orbitals.

As an example, let's consider the operation of a set of 1-body excitation
operators (which we may think of as a 1-body alpha-spin cluster operator from
coupled cluster theory) on the 2-electron alpha/beta spin portion of the
Hamiltonian.  To start, we define some short names for the type groups we'll
need using the built in types that the program has.  One can certainly define
new type groups, but the defaults are provided for consistency and because
some of the program's advanced functions rely on them.

>> tg_a = sqa.options.alpha_type
>> tg_b = sqa.options.beta_type
>> tg_c = sqa.options.core_type
>> tg_v = sqa.options.virtual_type
>> tg_h = sqa.options.core_type + sqa.options.virtual_type

The last type group, tg_h, is formed using python's ability to add tuples.
Now let's define the indices we'll need, using p as the base name for the
Hamiltonian's indices and q as the base name for the cluster operator's
indices.

>> p0 = sqa.index('p0', [tg_a, tg_h], True)
>> p1 = sqa.index('p1', [tg_b, tg_h], True)
>> p2 = sqa.index('p2', [tg_a, tg_h], True)
>> p3 = sqa.index('p3', [tg_b, tg_h], True)
>> q0 = sqa.index('q0', [tg_a, tg_v], True)
>> q1 = sqa.index('q1', [tg_a, tg_c], True)

Note that the Hamiltonian's indices range over both core and virtual orbitals,
while the cluster operator has one virtual index and one core index.  Now
let's define the tensors, remembering the Hamiltonian's transpose symmetry.

>> sym_h2abab = sqa.symmetry((2,3,0,1),1)
>> h2abab = sqa.tensor('h2abab', [p0, p1, p2, p3], [sym_h2abab])
>> p0_c = sqa.creOp(p0)
>> p1_c = sqa.creOp(p1)
>> p2_d = sqa.desOp(p2)
>> p3_d = sqa.desOp(p3)
>> t1aa = sqa.tensor('t1aa', [q0, q1], [])
>> q0_c = sqa.creOp(q0)
>> q1_d = sqa.desOp(q1)

Now let's define a term that represents the multiplication of the abab part of
the Hamiltonian by the cluster operator.  Remember that all the indices are
dummy and summed over.

>> term = sqa.term(1.0, [], [t1aa, q0_c, q1_d, h2abab, p0_c, p1_c, p3_d, p2_d])
>> print term

We should see the term:
 (   1.00000) t1aa(q0,q1) cre(q0) des(q1) h2abab(p0,p1,p2,p3) cre(p0) cre(p1) des(p3) des(p2)

Let's now put the term in normal order, which should produce a couple of delta
functions, one of which will be zero due to alpha-beta orthogonality.

>> result = sqa.normalOrder(term)
>> for t in result:
>>   print t

The result of normal order is:
 (  -1.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) cre(p0) cre(p1) cre(q0) des(p2) des(p3) des(q1) 
 (   1.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) kdelta(q1,p0) cre(p1) cre(q0) des(p2) des(p3) 
 (  -1.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) kdelta(q1,p1) cre(p0) cre(q0) des(p2) des(p3)

As expected there are two terms containing delta functions.  The kdelta(q1,p1)
tensor is identically zero because q1 is an alpha index while p1 is a beta
index.  kdelta(q1,p0) is not zero because both its indices are alpha and both
range over (at least) the core orbitals.  We can now use the program to
contract the delta functions.

>> for t in result:
>>   t.contractDeltaFuncs()
>>   print t

After contraction our terms are:
 (  -1.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) cre(p0) cre(p1) cre(q0) des(p2) des(p3) des(q1) 
 (   1.00000) t1aa(q0,p0) h2abab(p0,p1,p2,p3) cre(p1) cre(q0) des(p2) des(p3) 
 (   0.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) kdelta(q1,p1) cre(p0) cre(q0) des(p2) des(p3)

While the program did not delete the term with kdelta(q1,p1), it did set the
numerical coefficient to zero, which is effectively the same thing.  Notice
that in the term that had the kdelta(q1,p0) delta function, the q1 index was
replaced everywhere by p0 and the delta function was removed.  The program
also did something more subtle.  If we inspect the type groups on the 'new' p0
in the second term we will see that they are different from what we initially
set them to.

>> print result[1].tensors[1].indices[0].indType

This prints:
(('alpha',), ('core',))

We see that the 'new' p0 still has two type groups and is still an alpha
index, but now it only ranges over core indices.  This is because the terms in
the summation in which p0 was a virtual index would have caused kdelta(q1,p0)
to evaluate to zero, as q1 ranged over only core orbitals.  In contrast, the
p0 index in the term that had no delta function still ranges over both core
and virtual orbitals.

>> print result[0].tensors[1].indices[0].indType

This prints
(('alpha',), ('core', 'virtual'))
which is what we expected.

Finally, we can automatically remove the zero term using the termChop
function, which deletes terms with very small numerical coefficients.

>> sqa.termChop(result)
>> for t in result:
>>   print t

We are left with only the two non-zero terms:
 (  -1.00000) t1aa(q0,q1) h2abab(p0,p1,p2,p3) cre(p0) cre(p1) cre(q0) des(p2) des(p3) des(q1) 
 (   1.00000) t1aa(q0,p0) h2abab(p0,p1,p2,p3) cre(p1) cre(q0) des(p2) des(p3)


#####################################################################
# 4. Commutators
#####################################################################

Now let's put everything together to compute something relevant to Canonical
Transformation theory:  the commutator of the all alpha portions of the 2-body
anti-symmetric amplitudes and the 2-body Hamiltonian.  As usual, the first
step is to define some indices, which we do by creating a list of indices
called p.

>> p = [sqa.index('p%i' %j, [sqa.options.alpha_type], True) for j in range(8)]

Next define the amplitude and Hamiltonian tensors (remembering symmetries).

>> a2aaaa_syms = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
>> a2aaaa = sqa.tensor('a2aaaa', p[0:4], a2aaaa_syms)
>> h2aaaa_syms = [sqa.symmetry((2,3,0,1),1), sqa.symmetry((0,1,3,2), -1)]
>> h2aaaa = sqa.tensor('h2aaaa', p[4:8], h2aaaa_syms)

Now define the creation and destruction operator strings for the amplitude
and Hamiltonian operators.

>> a_ops1 = [sqa.creOp(p[0]), sqa.creOp(p[1]), sqa.desOp(p[3]), sqa.desOp(p[2])]
>> a_ops2 = [sqa.creOp(p[2]), sqa.creOp(p[3]), sqa.desOp(p[1]), sqa.desOp(p[0])]
>> h_ops  = [sqa.creOp(p[4]), sqa.creOp(p[5]), sqa.desOp(p[7]), sqa.desOp(p[6])]

In CT the amplitude operator is anti-symmetric, so we must write it as the sum
of two terms, which are the negative transposes of each other.  In the
program, sums of terms are stored as a list of the terms.

>> ampOp = [sqa.term(1.0, [], [a2aaaa] + a_ops1), sqa.term(-1.0, [], [a2aaaa] + a_ops2)]

Finally, we define the sum of terms (actually just a single term, but in list
form) for the Hamiltonian.

>> hamOp = [sqa.term(1.0, [], [h2aaaa] + h_ops)]

Taking the commutator of the amplitude and Hamiltonian operators can now be
done automatically using the commutator function.  This function expands all
the terms of the commutator, normal orders them, contracts all delta
functions, converts the terms to canonical form, and then combines any like
terms.  Now that's automation!

>> result = sqa.commutator(hamOp, ampOp)

If sqa.options.verbose is set to True, we will see a printout as the program
converts the resulting terms to canonical form:

Combining like terms:
Converting 28 terms to canonical form...
     0   (   1.00000) h2aaaa(p4,p5,p6,p7) a2aaaa(p0,p1,p2,p3) cre(p0) cre(p1) cre(p4) cre(p5) des(p2) des(p3) des(p6) des(p7) 
     1   (  -1.00000) h2aaaa(p4,p5,p6,p0) a2aaaa(p0,p1,p2,p3) cre(p1) cre(p4) cre(p5) des(p2) des(p3) des(p6) 
     2   (   1.00000) h2aaaa(p4,p5,p6,p1) a2aaaa(p0,p1,p2,p3) cre(p0) cre(p4) cre(p5) des(p2) des(p3) des(p6) 
     3   (   1.00000) h2aaaa(p4,p5,p0,p7) a2aaaa(p0,p1,p2,p3) cre(p1) cre(p4) cre(p5) des(p2) des(p3) des(p7) 
     4   (  -1.00000) h2aaaa(p4,p5,p1,p7) a2aaaa(p0,p1,p2,p3) cre(p0) cre(p4) cre(p5) des(p2) des(p3) des(p7) 
     5   (   1.00000) h2aaaa(p4,p5,p1,p0) a2aaaa(p0,p1,p2,p3) cre(p4) cre(p5) des(p2) des(p3) 
     6   (  -1.00000) h2aaaa(p4,p5,p0,p1) a2aaaa(p0,p1,p2,p3) cre(p4) cre(p5) des(p2) des(p3) 
     7   (  -1.00000) a2aaaa(p0,p1,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) cre(p4) cre(p5) des(p2) des(p3) des(p6) des(p7) 
     8   (   1.00000) a2aaaa(p0,p1,p2,p4) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) cre(p5) des(p2) des(p6) des(p7) 
     9   (  -1.00000) a2aaaa(p0,p1,p2,p5) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) cre(p4) des(p2) des(p6) des(p7) 
    10   (  -1.00000) a2aaaa(p0,p1,p4,p3) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) cre(p5) des(p3) des(p6) des(p7) 
    11   (   1.00000) a2aaaa(p0,p1,p5,p3) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) cre(p4) des(p3) des(p6) des(p7) 
    12   (  -1.00000) a2aaaa(p0,p1,p5,p4) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) des(p6) des(p7) 
    13   (   1.00000) a2aaaa(p0,p1,p4,p5) h2aaaa(p4,p5,p6,p7) cre(p0) cre(p1) des(p6) des(p7) 
    14   (  -1.00000) h2aaaa(p4,p5,p6,p7) a2aaaa(p0,p1,p2,p3) cre(p2) cre(p3) cre(p4) cre(p5) des(p0) des(p1) des(p6) des(p7) 
    15   (   1.00000) h2aaaa(p4,p5,p6,p2) a2aaaa(p0,p1,p2,p3) cre(p3) cre(p4) cre(p5) des(p0) des(p1) des(p6) 
    16   (  -1.00000) h2aaaa(p4,p5,p6,p3) a2aaaa(p0,p1,p2,p3) cre(p2) cre(p4) cre(p5) des(p0) des(p1) des(p6) 
    17   (  -1.00000) h2aaaa(p4,p5,p2,p7) a2aaaa(p0,p1,p2,p3) cre(p3) cre(p4) cre(p5) des(p0) des(p1) des(p7) 
    18   (   1.00000) h2aaaa(p4,p5,p3,p7) a2aaaa(p0,p1,p2,p3) cre(p2) cre(p4) cre(p5) des(p0) des(p1) des(p7) 
    19   (  -1.00000) h2aaaa(p4,p5,p3,p2) a2aaaa(p0,p1,p2,p3) cre(p4) cre(p5) des(p0) des(p1) 
    20   (   1.00000) h2aaaa(p4,p5,p2,p3) a2aaaa(p0,p1,p2,p3) cre(p4) cre(p5) des(p0) des(p1) 
    21   (   1.00000) a2aaaa(p0,p1,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) cre(p4) cre(p5) des(p0) des(p1) des(p6) des(p7) 
    22   (  -1.00000) a2aaaa(p0,p4,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) cre(p5) des(p0) des(p6) des(p7) 
    23   (   1.00000) a2aaaa(p0,p5,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) cre(p4) des(p0) des(p6) des(p7) 
    24   (   1.00000) a2aaaa(p4,p1,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) cre(p5) des(p1) des(p6) des(p7) 
    25   (  -1.00000) a2aaaa(p5,p1,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) cre(p4) des(p1) des(p6) des(p7) 
    26   (   1.00000) a2aaaa(p5,p4,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) des(p6) des(p7) 
    27   (  -1.00000) a2aaaa(p4,p5,p2,p3) h2aaaa(p4,p5,p6,p7) cre(p2) cre(p3) des(p6) des(p7) 
Finished combining terms in 3.876 seconds

Now we can examine the result by printing its terms.

>> for t in result:
>>   print t

This should give:
 (  -2.00000) a2aaaa(a,b,c,d) h2aaaa(a,b,e,f) cre(c) cre(d) des(e) des(f) 
 (  -2.00000) a2aaaa(a,b,c,d) h2aaaa(a,b,e,f) cre(e) cre(f) des(c) des(d) 
 (   2.00000) a2aaaa(a,b,c,d) h2aaaa(c,d,e,f) cre(a) cre(b) des(e) des(f) 
 (   2.00000) a2aaaa(a,b,c,d) h2aaaa(c,d,e,f) cre(e) cre(f) des(a) des(b) 
 (   4.00000) a2aaaa(a,b,c,d) h2aaaa(a,e,f,g) cre(b) cre(f) cre(g) des(c) des(d) des(e) 
 (   4.00000) a2aaaa(a,b,c,d) h2aaaa(a,e,f,g) cre(c) cre(d) cre(e) des(b) des(f) des(g) 
 (  -4.00000) a2aaaa(a,b,c,d) h2aaaa(c,e,f,g) cre(a) cre(b) cre(e) des(d) des(f) des(g) 
 (  -4.00000) a2aaaa(a,b,c,d) h2aaaa(c,e,f,g) cre(d) cre(f) cre(g) des(a) des(b) des(e) 

Note that our original index names have been replaced by their canonical
counterparts, which is not a problem as all indices are summed over and so
their names are irrelevant.  Of course in CT theory we also have to worry
about the other portions of the amplitude and Hamiltonian operators, but this
example shows how one can symbolically evaluate commutators automatically.


#####################################################################
# 5. Decompositions
#####################################################################

Operator and density matrix decompositions, described in J. Chem. Phys.
130, 124102 (2009), are a central part of Canonical Transformation
theory and have been implemented in this program.

The available decompositions for spin-orbital reduced density matrices
(RDMs) are:

decomp_3rdms_to_2rdms_so
decomp_4rdms_to_2rdms_so
decomp_4rdms_to_3rdms_so

These functions take a list of terms, a string, and various RDM tensors
(depending on the function) as arguments.  They search the list of terms for
any tensor whose name matches the provided string and decompose that tensor
(which is presumably a density matrix of the correct order) into lower body
RDMs.  This amounts to replacing the term containing the original RDM with a
group of terms, one for each term in the decomposition.

The analogous functions for decomposing operators are:

decomp_3ops_to_2ops_2rdms_so
decomp_3ops_to_2ops_3rdms_so
decomp_4ops_to_2ops_2rdms_so
decomp_4ops_to_2ops_3rdms_so

These take a list of terms and various RDM tensors (again depending on the
function) as inputs.  These functions act similarly to the RDM decompositions,
although here the operators are detected automatically rather than by name.
As a consequence, only normal ordered operators will be decomposed.

The decomposition functions require all indices of the RDMs, creation
operators, and destruction operators to have either the sqa.options.alpha_type
or sqa.options.beta_type as one of their type groups.  This is so that the
correct RDM tensors can be chosen during decomposition.

As an example of using the decomposition, lets create a term that involves the
3-body all-alpha-spin RDM and then decompose it.  Note that we use the same
indices for the 1- and 2-body RDMs as we do for the 3-body RDM.  This does not
matter, because those of the 1- and 2-body RDMs will be overwritten during the
decomposition.  What is important is that the 1- and 2-body RDMs have the
correct alpha and beta index combinations.

>> p = [sqa.index('p%i' %j, [sqa.options.alpha_type], False) for j in range(8)]
>> q = [sqa.index('q%i' %j, [sqa.options.beta_type],  False) for j in range(8)]
>> d1_sym = [sqa.symmetry((1,0),1)]
>> d2_sym_aaaa = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((2,3,0,1),1)]
>> d2_sym_abab = [sqa.symmetry((2,3,0,1),1)]
>> d3_sym_aaaaaa = [sqa.symmetry((1,0,2,3,4,5),-1), sqa.symmetry((0,2,1,3,4,5),-1), sqa.symmetry((3,4,5,0,1,2),1)]
>> d1aa = sqa.tensor('d1aa', p[0:2], d1_sym)
>> d1bb = sqa.tensor('d1bb', q[0:2], d1_sym)
>> d2aaaa = sqa.tensor('d2aaaa', p[0:4], d2_sym_aaaa)
>> d2bbbb = sqa.tensor('d2bbbb', q[0:4], d2_sym_aaaa)
>> d2abab = sqa.tensor('d2abab', [p[0],q[0],p[1],q[1]], d2_sym_abab)
>> d3aaaaaa = sqa.tensor('d3aaaaaa', p[0:6], d3_sym_aaaaaa)
>> term = sqa.term(1.0, [], [d3aaaaaa, sqa.creOp(p[6]), sqa.desOp(p[7])])
>> print term

This should print out our initial term, which is just a 1-body excitation
operator multiplied by the 3-body RDM:
 (   1.00000) d3aaaaaa(p0,p1,p2,p3,p4,p5) cre(p6) des(p7) 

Now put the term in a list and pass the list to the decomposition function to
replace the 3-body RDM by its 1- and 2-body decomposition.

>> term_list = [term]
>> sqa.decomp_3rdms_to_2rdms_so(term_list, 'd3aaaaaa', d1aa, d1bb, d2aaaa, d2bbbb, d2abab)

If sqa.options.verbose is set to true, you should see a printout telling you
how many decompositions were performed.  In this case there should be just 1:
decomposed 1 3-body RDMs

Our original term in the list has now been replaced by the terms making up its
decomposition.  To see them, just print them out.

>> for t in term_list:
>>   print t

This prints:
 (   1.00000) d1aa(p0,p3) d2aaaa(p1,p2,p4,p5) cre(p6) des(p7) 
 (  -1.00000) d1aa(p0,p4) d2aaaa(p1,p2,p3,p5) cre(p6) des(p7) 
 (   1.00000) d1aa(p0,p5) d2aaaa(p1,p2,p3,p4) cre(p6) des(p7) 
 (  -1.00000) d1aa(p1,p3) d2aaaa(p0,p2,p4,p5) cre(p6) des(p7) 
 (   1.00000) d1aa(p1,p4) d2aaaa(p0,p2,p3,p5) cre(p6) des(p7) 
 (  -1.00000) d1aa(p1,p5) d2aaaa(p0,p2,p3,p4) cre(p6) des(p7) 
 (   1.00000) d1aa(p2,p3) d2aaaa(p0,p1,p4,p5) cre(p6) des(p7) 
 (  -1.00000) d1aa(p2,p4) d2aaaa(p0,p1,p3,p5) cre(p6) des(p7) 
 (   1.00000) d1aa(p2,p5) d2aaaa(p0,p1,p3,p4) cre(p6) des(p7) 
 (  -2.00000) d1aa(p0,p3) d1aa(p1,p4) d1aa(p2,p5) cre(p6) des(p7) 
 (   2.00000) d1aa(p0,p3) d1aa(p1,p5) d1aa(p2,p4) cre(p6) des(p7) 
 (   2.00000) d1aa(p0,p4) d1aa(p1,p3) d1aa(p2,p5) cre(p6) des(p7) 
 (  -2.00000) d1aa(p0,p4) d1aa(p1,p5) d1aa(p2,p3) cre(p6) des(p7) 
 (  -2.00000) d1aa(p0,p5) d1aa(p1,p3) d1aa(p2,p4) cre(p6) des(p7) 
 (   2.00000) d1aa(p0,p5) d1aa(p1,p4) d1aa(p2,p3) cre(p6) des(p7) 

That's all for this tutorial, at least for now.
