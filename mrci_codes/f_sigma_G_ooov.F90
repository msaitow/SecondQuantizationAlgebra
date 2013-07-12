#include "../f_ct.fh"



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no0_x0(sc, ic, sd, id, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sd, id
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_G_ooov_no0_x0(sc, ic, sd, id, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_G_ooov_no0_x0


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no0_x0(s_c, i_c, s_d, i_d, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_d, s_d
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_e, i_e, s_f, i_f

!X(a,b,e,f) <-- 
! (   1.00000) V2(c,e,d,f) T2(a,b,c,d) 
! where i_c in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_e, s_f) .and. & 
IEOR(s_c, s_e) == IEOR(s_d, s_f) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)

X_(s_f, s_e, s_b, s_a)%array(i_f, i_e, i_b, i_a) = X_(s_f, s_e, s_b, s_a)%array(i_f, i_e, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_e)%array(i_f, i_d, i_e) & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 4
! External_count : 2
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no0_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no1_x0(X, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: X(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_G_ooov_no1_x0(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_G_ooov_no1_x0


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no1_x0(X_, S0_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: S0_
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_e, i_e, s_f, i_f

!S0() <-- 
! (   1.00000) E2(e,a,f,b) X(a,b,e,f) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
if( &
IEOR(s_e, s_a) == IEOR(s_f, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_e, s_f)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)

S0_ = S0_ &
  + 1.0d+00 & 
  * E2_(s_b, s_f, s_a, s_e)%array(i_b, i_f, i_a, i_e) & 
  * X_(s_f, s_e, s_b, s_a)%array(i_f, i_e, i_b, i_a)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 4
! External_count : 0
! Total_count    : 4


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no1_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no0_x1(sd, id, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id
real(kind=8), intent(inout) :: T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sd

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_G_ooov_no0_x1(sd, id, av2_i, xa, d2, nir, nsym, psym)

deallocate(xa)
deallocate(av2_i)

end subroutine g_if_sigma_G_ooov_no0_x1


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no0_x1(s_d, i_d, T2_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock1), intent(inout) :: X_(0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e

!X(e,d) <-- 
! (   1.00000) T2(a,b,c,d) E2(c,a,e,b) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(S_e, s_d) == 0 .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_d) .and. &
IEOR(s_c, s_a) == IEOR(s_e, s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e)%array(i_e) = X_(s_e)%array(i_e) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * E2_(s_b, s_e, s_a, s_c)%array(i_b, i_e, i_a, i_c)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 4
! External_count : 1
! Total_count    : 5


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no0_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no1_x1(sd, id, h, X, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id
real(kind=8), intent(inout) :: h(*), X(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sd

call set_symblock_xa(sleft, x, nir, nsym, psym) ! -> xa (allocate) 
call g_sigma_G_ooov_no1_x1(sd, id, h1, xa, S0, nir, nsym, psym)

deallocate(h1)
deallocate(xa)

end subroutine g_if_sigma_G_ooov_no1_x1


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no1_x1(s_d, i_d, h_, X_, S0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)
real(kind=8)                   :: S0_
type(symblock1), intent(inout) :: X_(0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_e, i_e

!S0() <-- 
! (   1.00000) h(d,e) X(e,d) 
do s_e = 0, nir-1
if( &
IEOR(S_d, s_e) == 0 .and. &
IEOR(S_e, s_d) == 0) then
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

S0_ = S0_ &
  + 1.0d+00 & 
  * h_(s_e, s_d)%array(i_e, i_d) & 
  * X_(s_e)%array(i_e)

end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
! Loop_count     : 1
! External_count : 1
! Total_count    : 2


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no1_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no0_x2(sd, id, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sd, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sd

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_G_ooov_no0_x2(sd, id, h2_i, xaaa, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)

end subroutine g_if_sigma_G_ooov_no0_x2


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no0_x2(s_d, i_d, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_f, i_f
integer :: s_g, i_g

!X(a,b,c,d) <-- 
! (   1.00000) V2(d,f,e,g) E3(c,a,e,g,f,b) 
! where i_d in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
do s_g = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_c, s_d) .and. & 
IEOR(s_d, s_f) == IEOR(s_e, s_g) .and. &
IEOR(IEOR(s_c, s_a),s_e) == IEOR(IEOR(s_g, s_f),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_g = psym(I_BEGIN, I_O, s_g), psym(I_END, I_O, s_g)

X_(s_c, s_b, s_a)%array(i_c, i_b, i_a) = X_(s_c, s_b, s_a)%array(i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_g, s_e, s_f)%array(i_g, i_e, i_f) & 
  * E3_(s_b, s_f, s_g, s_e, s_a, s_c)%array(i_b, i_f, i_g, i_e, i_a, i_c)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 6
! External_count : 1
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no0_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_ooov_no1_x2(sd, id, T2, X, S0, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id
real(kind=8), intent(inout) :: T2(*), X(*), S0
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sd

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_G_ooov_no1_x2(sd, id, av2_i, xaaa, S0, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i)

end subroutine g_if_sigma_G_ooov_no1_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:5:3

subroutine g_sigma_G_ooov_no1_x2(s_d, i_d, T2_, X_, S0_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
real(kind=8)                   :: S0_
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c

!S0() <-- 
! (   1.00000) T2(a,b,c,d) X(a,b,c,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_c, s_d) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)

S0_ = S0_ &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 1
! Total_count    : 4


! FEMTO END    **************************************************************

end subroutine g_sigma_G_ooov_no1_x2


! -----------------------------------------------------------------------
