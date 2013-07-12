#include "../f_ct.fh"



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x0(sa, ia, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x0(sa, ia, h2_i, xaaa, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x0


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x0(s_a, i_a, V2_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_b, i_b, s_c, i_c, s_d, i_d

!X(a,b,Ai,Aj) <-- 
! (   1.00000) V2(a,c,b,d) E2(Ai,c,Aj,d) 
! where i_a in {core}
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(s_a, s_c) == IEOR(s_b, s_d) .and. &
IEOR(s_Ai, s_c) == IEOR(s_Aj, s_d)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b) = X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b) &
  + 1.0d+00 & 
  * V2_(s_d, s_b, s_c)%array(i_d, i_b, i_c) & 
  * E2_(s_d, s_Aj, s_c, s_Ai)%array(i_d, i_Aj, i_c, i_Ai)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x0(sa, ia, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x0(sa, ia, sVa, iVa, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x0


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x0(s_a, i_a, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,Ak,Va) X(a,b,Ai,Aj) 
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b)

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

end subroutine g_sigma_ooov_ooov_no1_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x1(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x1(sd, id, sVa, iVa, h2_i, av2_i, xaaa, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x1


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x1(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_Ak, i_Ak, s_c, i_c, s_b, i_b

!X(a,b,Ak,Va) <-- 
! (   1.00000) V2(Va,d,Ak,c) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Va, s_d) == IEOR(s_Ak, s_c) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) = X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_c, s_Ak, s_d)%array(i_c, i_Ak, i_d) & 
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

end subroutine g_sigma_ooov_ooov_no0_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x1(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x1(sVa, iVa, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x1


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x1(s_Va, i_Va, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x2(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x2(sd, id, sVa, iVa, h2_i, av2_i, xaaa, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x2


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x2(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_Ak, i_Ak, s_c, i_c, s_b, i_b

!X(a,b,Ak,Va) <-- 
! (   1.00000) V2(Va,c,Ak,d) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Va, s_c) == IEOR(s_Ak, s_d) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) = X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_d, s_Ak, s_c)%array(i_d, i_Ak, i_c) & 
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

end subroutine g_sigma_ooov_ooov_no0_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x2(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x2(sVa, iVa, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x2(s_Va, i_Va, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E2(Ai,b,Aj,a) X(a,b,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Ai, s_b) == IEOR(s_Aj, s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_a, s_Aj, s_b, s_Ai)%array(i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x3(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x3(h1, xaaaa, d2, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x3(h_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!X(Ai,a,b,Aj) <-- 
! (   1.00000) h(a,c) E2(Ai,c,Aj,b) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_a) == IEOR(s_b, s_Aj) .and. & 
IEOR(S_a, s_c) == 0 .and. &
IEOR(s_Ai, s_c) == IEOR(s_Aj, s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai) = X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai) &
  + 1.0d+00 & 
  * h_(s_c, s_a)%array(i_c, i_a) & 
  * E2_(s_b, s_Aj, s_c, s_Ai)%array(i_b, i_Aj, i_c, i_Ai)

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
! Loop_count     : 5
! External_count : 0
! Total_count    : 5


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x3(sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x3(sVa, iVa, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x3(s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,Ak,Va) X(Ai,a,b,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_Ai, s_a) == IEOR(s_b, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x4(h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

real(kind=8), intent(inout) :: h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x4(h1, xaaaa, d2, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x4


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x4(h_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!X(Ai,a,b,Aj) <-- 
! (   1.00000) h(b,c) E2(Ai,a,Aj,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_a) == IEOR(s_b, s_Aj) .and. & 
IEOR(S_b, s_c) == 0 .and. &
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_c)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai) = X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai) &
  + 1.0d+00 & 
  * h_(s_c, s_b)%array(i_c, i_b) & 
  * E2_(s_c, s_Aj, s_a, s_Ai)%array(i_c, i_Aj, i_a, i_Ai)

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
! Loop_count     : 5
! External_count : 0
! Total_count    : 5


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x4


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x4(sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x4(sVa, iVa, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x4


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x4(s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,Ak,Va) X(Ai,a,b,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_Ai, s_a) == IEOR(s_b, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * X_(s_Aj, s_b, s_a, s_Ai)%array(i_Aj, i_b, i_a, i_Ai)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x4


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x5(sc, ic, sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x5(sc, ic, sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x5


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x5(s_c, i_c, s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_Ak, i_Ak

!X(a,b,Ak,Va) <-- 
! (   1.00000) T2(a,b,Ak,c) h(Va,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_c) .and. &
IEOR(S_Va, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) = X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * h_(s_c, s_Va)%array(i_c, i_Va)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 2
! Total_count    : 5


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x5(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x5(sVa, iVa, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x5


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x5(s_Va, i_Va, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x6(sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x6(sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x6


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x6(s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Ak, i_Ak

!X(a,b,Ak,Va) <-- 
! (   1.00000) T2(a,b,c,Va) h(Ak,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(S_Ak, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) = X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_c, s_Ak)%array(i_c, i_Ak)

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

end subroutine g_sigma_ooov_ooov_no0_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x6(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x6(sVa, iVa, xaaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x6


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x6(s_Va, i_Va, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x7(sa, ia, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x7(sa, ia, h2_i, xaaa, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x7


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x7(s_a, i_a, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) V2(a,d,c,e) E3(Ai,d,Aj,b,c,e) 
! where i_a in {core}
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(s_a, s_d) == IEOR(s_c, s_e) .and. &
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_b, s_c),s_e)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b) = X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * E3_(s_e, s_c, s_b, s_Aj, s_d, s_Ai)%array(i_e, i_c, i_b, i_Aj, i_d, i_Ai)

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

end subroutine g_sigma_ooov_ooov_no0_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x7(sa, ia, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x7(sa, ia, sVa, iVa, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x7


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x7(s_a, i_a, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,Ak,Va) X(a,b,Ai,Aj) 
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_b)%array(i_Aj, i_Ai, i_b)

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

end subroutine g_sigma_ooov_ooov_no1_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x8(sb, ib, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sb, ib
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sb, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sb

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x8(sb, ib, h2_i, xaaa, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x8


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x8(s_b, i_b, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_b, s_b
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) V2(b,d,c,e) E3(Ai,a,Aj,d,c,e) 
! where i_b in {core}
do s_a = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(s_b, s_d) == IEOR(s_c, s_e) .and. &
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_d, s_c),s_e)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_a)%array(i_Aj, i_Ai, i_a) = X_(s_Aj, s_Ai, s_a)%array(i_Aj, i_Ai, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * E3_(s_e, s_c, s_d, s_Aj, s_a, s_Ai)%array(i_e, i_c, i_d, i_Aj, i_a, i_Ai)

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

end subroutine g_sigma_ooov_ooov_no0_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x8(sb, ib, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sb, ib, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sb

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x8(sb, ib, sVa, iVa, av2_i, xaaa, av2_i2, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x8


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x8(s_b, i_b, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_b, s_b
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_Ak, i_Ak, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,Ak,Va) X(a,b,Ai,Aj) 
do s_a = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_a)%array(i_Aj, i_Ai, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x9(sc, ic, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x9(sc, ic, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x9


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x9(s_c, i_c, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_Ak, i_Ak, s_d, i_d, s_e, i_e

!X(a,b,d,e,Ak,Va) <-- 
! (   1.00000) V2(Va,d,c,e) T2(a,b,Ak,c) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Ak),s_Va) .and. & 
IEOR(s_Va, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_c)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a) = X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x9(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x9(sVa, iVa, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x9


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x9(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,a,Aj,d,e,b) X(a,b,d,e,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_d, s_e),s_b) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Ak),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_b, s_e, s_d, s_Aj, s_a, s_Ai)%array(i_b, i_e, i_d, i_Aj, i_a, i_Ai) & 
  * X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x10(sc, ic, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x10(sc, ic, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x10


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x10(s_c, i_c, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_Ak, i_Ak, s_d, i_d, s_e, i_e

!X(a,b,d,e,Ak,Va) <-- 
! (   1.00000) V2(Va,c,d,e) T2(a,b,Ak,c) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Ak),s_Va) .and. & 
IEOR(s_Va, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_Ak, s_c)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a) = X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x10


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x10(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x10(sVa, iVa, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x10


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x10(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,a,Aj,b,d,e) X(a,b,d,e,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_b, s_d),s_e) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Ak),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_e, s_d, s_b, s_Aj, s_a, s_Ai)%array(i_e, i_d, i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Ak, s_e, s_d, s_b, s_a)%array(i_Ak, i_e, i_d, i_b, i_a)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x10


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x11(sAk, iAk, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sAk, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x11(sAk, iAk, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x11


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x11(s_Ak, i_Ak, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e
integer :: s_Ai, i_Ai, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(Ak,d,a,e) E3(Ai,d,Aj,b,c,e) 
! where i_Ak in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_Ak, s_d) == IEOR(s_a, s_e) .and. &
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_b, s_c),s_e)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) = X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_a, s_d)%array(i_e, i_a, i_d) & 
  * E3_(s_e, s_c, s_b, s_Aj, s_d, s_Ai)%array(i_e, i_c, i_b, i_Aj, i_d, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x11


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x11(sAk, iAk, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x11(sAk, iAk, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x11


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x11(s_Ak, i_Ak, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x11


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x12(sAk, iAk, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sAk, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x12(sAk, iAk, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x12


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x12(s_Ak, i_Ak, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e
integer :: s_Ai, i_Ai, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(Ak,d,b,e) E3(Ai,d,Aj,e,c,a) 
! where i_Ak in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_Ak, s_d) == IEOR(s_b, s_e) .and. &
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_e, s_c),s_a)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) = X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_b, s_d)%array(i_e, i_b, i_d) & 
  * E3_(s_a, s_c, s_e, s_Aj, s_d, s_Ai)%array(i_a, i_c, i_e, i_Aj, i_d, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x12


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x12(sAk, iAk, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x12(sAk, iAk, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x12


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x12(s_Ak, i_Ak, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x12


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x13(sAk, iAk, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sAk, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x13(sAk, iAk, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x13


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x13(s_Ak, i_Ak, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e
integer :: s_Ai, i_Ai, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(Ak,d,c,e) E3(Ai,d,Aj,b,e,a) 
! where i_Ak in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_Ak, s_d) == IEOR(s_c, s_e) .and. &
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_b, s_e),s_a)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) = X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * E3_(s_a, s_e, s_b, s_Aj, s_d, s_Ai)%array(i_a, i_e, i_b, i_Aj, i_d, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x13(sAk, iAk, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x13(sAk, iAk, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x13


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x13(s_Ak, i_Ak, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x14(sAk, iAk, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sAk, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x14(sAk, iAk, h2_i, xaaaaa, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x14


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x14(s_Ak, i_Ak, V2_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e
integer :: s_Ai, i_Ai, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(Ak,c,d,e) E3(Ai,a,Aj,b,d,e) 
! where i_Ak in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_Ak, s_c) == IEOR(s_d, s_e) .and. &
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_b, s_d),s_e)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) = X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
  * E3_(s_e, s_d, s_b, s_Aj, s_a, s_Ai)%array(i_e, i_d, i_b, i_Aj, i_a, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x14(sAk, iAk, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sAk

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x14(sAk, iAk, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x14


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x14(s_Ak, i_Ak, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_c, s_b, s_a)%array(i_Aj, i_Ai, i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x15(sa, ia, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x15(sa, ia, sVa, iVa, h2_i, av2_i, xaaa, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x15


! Converting terms with V[i*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x15(s_a, i_a, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(c,d,e,Va) <-- 
! (   1.00000) V2(a,d,b,e) T2(a,b,c,Va) 
! where i_a in {core}
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(s_c, s_d) == IEOR(s_e, s_Va) .and. & 
IEOR(s_a, s_d) == IEOR(s_b, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Va)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_c)%array(i_e, i_d, i_c) = X_(s_e, s_d, s_c)%array(i_e, i_d, i_c) &
  + 1.0d+00 & 
  * V2_(s_e, s_b, s_d)%array(i_e, i_b, i_d) & 
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

end subroutine g_sigma_ooov_ooov_no0_x15


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x15(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x15(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x15


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x15(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c, i_c, s_d, i_d, s_e, i_e, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) E3(Ai,Ak,Aj,e,c,d) X(c,d,e,Va) 
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_e, s_c),s_d) .and. &
IEOR(s_c, s_d) == IEOR(s_e, s_Va)) then
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * E3_(s_d, s_c, s_e, s_Aj, s_Ak, s_Ai)%array(i_d, i_c, i_e, i_Aj, i_Ak, i_Ai) & 
  * X_(s_e, s_d, s_c)%array(i_e, i_d, i_c)

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

end subroutine g_sigma_ooov_ooov_no1_x15


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x16(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x16(sd, id, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x16


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x16(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c, i_c, s_a, i_a, s_Ak, i_Ak, s_b, i_b, s_e, i_e

!X(a,b,c,e,Ak,Va) <-- 
! (   1.00000) V2(Va,e,Ak,d) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_Ak),s_Va) .and. & 
IEOR(s_Va, s_e) == IEOR(s_Ak, s_d) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a) = X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_d, s_Ak, s_e)%array(i_d, i_Ak, i_e) & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x16


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x16(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x16(sVa, iVa, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x16


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x16(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,b,Aj,e,c,a) X(a,b,c,e,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_e, s_c),s_a) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_Ak),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_c, s_e, s_Aj, s_b, s_Ai)%array(i_a, i_c, i_e, i_Aj, i_b, i_Ai) & 
  * X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x16


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x17(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x17(sd, id, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x17


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no0_x17(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_c, i_c, s_a, i_a, s_Ak, i_Ak, s_b, i_b, s_e, i_e

!X(a,b,c,e,Ak,Va) <-- 
! (   1.00000) V2(Va,d,Ak,e) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_Ak),s_Va) .and. & 
IEOR(s_Va, s_d) == IEOR(s_Ak, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a) = X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_Ak, s_d)%array(i_e, i_Ak, i_d) & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x17


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x17(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x17(sVa, iVa, xaaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x17


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:43

subroutine g_sigma_ooov_ooov_no1_x17(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,e,Aj,b,c,a) X(a,b,c,e,Ak,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_e),s_Aj) == IEOR(IEOR(s_b, s_c),s_a) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_Ak),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_c, s_b, s_Aj, s_e, s_Ai)%array(i_a, i_c, i_b, i_Aj, i_e, i_Ai) & 
  * X_(s_Ak, s_e, s_c, s_b, s_a)%array(i_Ak, i_e, i_c, i_b, i_a)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 1
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x17


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x18(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x18(sd, id, sVa, iVa, h2_i, av2_i, xaaa, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x18


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x18(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e

!X(a,b,e,Va) <-- 
! (   1.00000) V2(Va,d,c,e) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_e, s_Va) .and. & 
IEOR(s_Va, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_b, s_a)%array(i_e, i_b, i_a) = X_(s_e, s_b, s_a)%array(i_e, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
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

end subroutine g_sigma_ooov_ooov_no0_x18


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x18(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x18(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x18


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x18(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_e, i_e, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,Ak,Aj,b,e,a) X(a,b,e,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_b, s_e),s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_e, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_e, s_b, s_Aj, s_Ak, s_Ai)%array(i_a, i_e, i_b, i_Aj, i_Ak, i_Ai) & 
  * X_(s_e, s_b, s_a)%array(i_e, i_b, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x18


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x19(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x19(sd, id, sVa, iVa, h2_i, av2_i, xaaa, nir, nsym, psym)

deallocate(xaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x19


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x19(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e

!X(a,b,e,Va) <-- 
! (   1.00000) V2(Va,c,d,e) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_e, s_Va) .and. & 
IEOR(s_Va, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_b, s_a)%array(i_e, i_b, i_a) = X_(s_e, s_b, s_a)%array(i_e, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
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

end subroutine g_sigma_ooov_ooov_no0_x19


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x19(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x19(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x19


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x19(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_e, i_e, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,Ak,Aj,a,e,b) X(a,b,e,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_a, s_e),s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_e, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_b, s_e, s_a, s_Aj, s_Ak, s_Ai)%array(i_b, i_e, i_a, i_Aj, i_Ak, i_Ai) & 
  * X_(s_e, s_b, s_a)%array(i_e, i_b, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x19


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x20(sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x20(sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x20


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x20(s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d

!X(b,c,d,Va) <-- 
! (   1.00000) T2(a,b,c,Va) h(a,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
if( &
IEOR(s_b, s_c) == IEOR(s_d, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(S_a, s_d) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)

X_(s_d, s_c, s_b)%array(i_d, i_c, i_b) = X_(s_d, s_c, s_b)%array(i_d, i_c, i_b) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_d, s_a)%array(i_d, i_a)

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

end subroutine g_sigma_ooov_ooov_no0_x20


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x20(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x20(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x20


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x20(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_b, i_b, s_c, i_c, s_d, i_d, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) E3(Ai,Ak,Aj,b,c,d) X(b,c,d,Va) 
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_b, s_c),s_d) .and. &
IEOR(s_b, s_c) == IEOR(s_d, s_Va)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * E3_(s_d, s_c, s_b, s_Aj, s_Ak, s_Ai)%array(i_d, i_c, i_b, i_Aj, i_Ak, i_Ai) & 
  * X_(s_d, s_c, s_b)%array(i_d, i_c, i_b)

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

end subroutine g_sigma_ooov_ooov_no1_x20


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x21(sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x21(sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x21


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x21(s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d

!X(a,c,d,Va) <-- 
! (   1.00000) T2(a,b,c,Va) h(b,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
if( &
IEOR(s_a, s_c) == IEOR(s_d, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(S_b, s_d) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)

X_(s_d, s_c, s_a)%array(i_d, i_c, i_a) = X_(s_d, s_c, s_a)%array(i_d, i_c, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_d, s_b)%array(i_d, i_b)

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

end subroutine g_sigma_ooov_ooov_no0_x21


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x21(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x21(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x21


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x21(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_c, i_c, s_d, i_d, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) E3(Ai,Ak,Aj,d,c,a) X(a,c,d,Va) 
do s_a = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_d, s_c),s_a) .and. &
IEOR(s_a, s_c) == IEOR(s_d, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * E3_(s_a, s_c, s_d, s_Aj, s_Ak, s_Ai)%array(i_a, i_c, i_d, i_Aj, i_Ak, i_Ai) & 
  * X_(s_d, s_c, s_a)%array(i_d, i_c, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x21


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x22(sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x22(sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x22


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x22(s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d

!X(a,b,d,Va) <-- 
! (   1.00000) T2(a,b,c,Va) h(c,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_d, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(S_c, s_d) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)

X_(s_d, s_b, s_a)%array(i_d, i_b, i_a) = X_(s_d, s_b, s_a)%array(i_d, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_d, s_c)%array(i_d, i_c)

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

end subroutine g_sigma_ooov_ooov_no0_x22


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x22(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x22(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x22


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x22(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,Ak,Aj,b,d,a) X(a,b,d,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_b, s_d),s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_d, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_d, s_b, s_Aj, s_Ak, s_Ai)%array(i_a, i_d, i_b, i_Aj, i_Ak, i_Ai) & 
  * X_(s_d, s_b, s_a)%array(i_d, i_b, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x22


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x23(sd, id, sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call g_sigma_ooov_ooov_no0_x23(sd, id, sVa, iVa, av2_i, h1, xaaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaaa)

end subroutine g_if_sigma_ooov_ooov_no0_x23


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x23(s_d, i_d, s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c

!X(a,b,c,Va) <-- 
! (   1.00000) T2(a,b,c,d) h(Va,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_d) .and. &
IEOR(S_Va, s_d) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)

X_(s_c, s_b, s_a)%array(i_c, i_b, i_a) = X_(s_c, s_b, s_a)%array(i_c, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_d, s_Va)%array(i_d, i_Va)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 2
! Total_count    : 5


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x23


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x23(sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaa(sleft, x, nir, nsym, psym) ! -> xaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x23(sVa, iVa, xaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x23


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x23(s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E3(Ai,Ak,Aj,b,c,a) X(a,b,c,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_b, s_c),s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_c, s_b, s_Aj, s_Ak, s_Ai)%array(i_a, i_c, i_b, i_Aj, i_Ak, i_Ai) & 
  * X_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

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

end subroutine g_sigma_ooov_ooov_no1_x23


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x24(sc, ic, se, ie, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, se, ie
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(se, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sc

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x24(sc, ic, se, ie, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x24


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x24(s_c, i_c, s_e, i_e, V2_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_e, s_e
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_f, i_f, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(e,a,d,f) E4(c,e,Ai,Ak,Aj,b,d,f) 
! where i_e in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_f = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_e, s_a) == IEOR(s_d, s_f) .and. &
IEOR(IEOR(s_c, s_e),IEOR(s_Ai, s_Ak)) == IEOR(IEOR(s_Aj, s_b),IEOR(s_d, s_f))) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ak, s_Ai, s_b, s_a)%array(i_Aj, i_Ak, i_Ai, i_b, i_a) = X_(s_Aj, s_Ak, s_Ai, s_b, s_a)%array(i_Aj, i_Ak, i_Ai, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_a)%array(i_f, i_d, i_a) & 
  * E4_(s_f, s_d, s_b, s_Aj, s_Ak, s_Ai)%array(i_f, i_d, i_b, i_Aj, i_Ak, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 2
! Total_count    : 9


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x24


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x24(sc, ic, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sc

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x24(sc, ic, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x24


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x24(s_c, i_c, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ak, s_Ai, s_b, s_a)%array(i_Aj, i_Ak, i_Ai, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x24


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x25(sAj, iAj, se, ie, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAj, iAj, se, ie
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(se, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sAj

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x25(sAj, iAj, se, ie, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x25


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x25(s_Aj, i_Aj, s_e, i_e, V2_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_e, s_e
integer, intent(in) :: i_Aj, s_Aj
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_f, i_f
integer :: s_Ai, i_Ai, s_Ak, i_Ak

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(e,b,d,f) E4(Aj,e,Ai,Ak,c,a,d,f) 
! where i_e in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_f = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_e, s_b) == IEOR(s_d, s_f) .and. &
IEOR(IEOR(s_Aj, s_e),IEOR(s_Ai, s_Ak)) == IEOR(IEOR(s_c, s_a),IEOR(s_d, s_f))) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

X_(s_Ak, s_Ai, s_c, s_b, s_a)%array(i_Ak, i_Ai, i_c, i_b, i_a) = X_(s_Ak, s_Ai, s_c, s_b, s_a)%array(i_Ak, i_Ai, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_b)%array(i_f, i_d, i_b) & 
  * E4_(s_f, s_d, s_a, s_c, s_Ak, s_Ai)%array(i_f, i_d, i_a, i_c, i_Ak, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 2
! Total_count    : 9


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x25


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x25(sAj, iAj, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAj, iAj, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sAj

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x25(sAj, iAj, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x25


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x25(s_Aj, i_Aj, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Aj, s_Aj
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_c, i_c, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (  -1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  - 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Ak, s_Ai, s_c, s_b, s_a)%array(i_Ak, i_Ai, i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x25


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x26(sa, ia, se, ie, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, se, ie
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(se, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x26(sa, ia, se, ie, h2_i, xaaaaa, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x26


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x26(s_a, i_a, s_e, i_e, V2_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_e, s_e
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_b, i_b, s_c, i_c, s_d, i_d, s_f, i_f, s_Ai, i_Ai
integer :: s_Ak, i_Ak, s_Aj, i_Aj

!X(a,b,c,Ai,Ak,Aj) <-- 
! (   1.00000) V2(e,c,d,f) E4(e,a,Ai,Ak,Aj,b,d,f) 
! where i_e in {core}
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_f = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj) .and. & 
IEOR(s_e, s_c) == IEOR(s_d, s_f) .and. &
IEOR(IEOR(s_e, s_a),IEOR(s_Ai, s_Ak)) == IEOR(IEOR(s_Aj, s_b),IEOR(s_d, s_f))) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ak, s_Ai, s_c, s_b)%array(i_Aj, i_Ak, i_Ai, i_c, i_b) = X_(s_Aj, s_Ak, s_Ai, s_c, s_b)%array(i_Aj, i_Ak, i_Ai, i_c, i_b) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_c)%array(i_f, i_d, i_c) & 
  * E4_(s_f, s_d, s_b, s_Aj, s_Ak, s_Ai)%array(i_f, i_d, i_b, i_Aj, i_Ak, i_Ai)

end do ! Irrep loop
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
end do ! Orbital loop
! Loop_count     : 7
! External_count : 2
! Total_count    : 9


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x26


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x26(sa, ia, sVa, iVa, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sVa, iVa
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x26(sa, ia, sVa, iVa, av2_i, xaaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x26


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x26(s_a, i_a, s_Va, i_Va, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_b, i_b, s_c, i_c, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) T2(a,b,c,Va) X(a,b,c,Ai,Ak,Aj) 
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_Ai, s_Ak),s_Aj)) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * X_(s_Aj, s_Ak, s_Ai, s_c, s_b)%array(i_Aj, i_Ak, i_Ai, i_c, i_b)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x26


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x27(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x27(sd, id, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x27


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x27(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_f, i_f

!X(a,b,c,e,f,Va) <-- 
! (   1.00000) V2(Va,e,d,f) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_f),s_Va) .and. & 
IEOR(s_Va, s_e) == IEOR(s_d, s_f) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)

X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a) = X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_e)%array(i_f, i_d, i_e) & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x27


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x27(sAi, iAi, sAk, iAk, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAi, iAi, sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x27(sAi, iAi, sAk, iAk, sVa, iVa, xaaaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x27


! Converting terms with no ERIs with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x27(s_Ai, i_Ai, s_Ak, i_Ak, s_Va, i_Va, X_, S2_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ai, s_Ai
integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_f, i_f
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E4(Ai,Ak,Aj,e,c,a,f,b) X(a,b,c,e,f,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),IEOR(s_Aj, s_e)) == IEOR(IEOR(s_c, s_a),IEOR(s_f, s_b)) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_f),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E4_(s_b, s_f, s_a, s_c, s_e, s_Aj)%array(i_b, i_f, i_a, i_c, i_e, i_Aj) & 
  * X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a)

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
! External_count : 3
! Total_count    : 9


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x27


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x28(sd, id, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call g_sigma_ooov_ooov_no0_x28(sd, id, sVa, iVa, h2_i, av2_i, xaaaaa, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x28


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x28(s_d, i_d, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_f, i_f

!X(a,b,c,e,f,Va) <-- 
! (   1.00000) V2(Va,d,e,f) T2(a,b,c,d) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_f),s_Va) .and. & 
IEOR(s_Va, s_d) == IEOR(s_e, s_f) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)

X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a) = X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_e, s_d)%array(i_f, i_e, i_d) & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a)

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
! Loop_count     : 5
! External_count : 2
! Total_count    : 7


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x28


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no1_x28(sAi, iAi, sAk, iAk, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sAi, iAi, sAk, iAk, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaaaaa(sleft, x, nir, nsym, psym) ! -> xaaaaa (allocate) 
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no1_x28(sAi, iAi, sAk, iAk, sVa, iVa, xaaaaa, av2_i2, d4_ij, nir, nsym, psym)

deallocate(xaaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_ooov_ooov_no1_x28


! Converting terms with no ERIs with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no1_x28(s_Ai, i_Ai, s_Ak, i_Ak, s_Va, i_Va, X_, S2_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Ai, s_Ai
integer, intent(in) :: i_Ak, s_Ak
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_e, i_e, s_f, i_f
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) E4(Ai,Ak,Aj,b,c,a,e,f) X(a,b,c,e,f,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_e = 0, nir-1
do s_f = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(IEOR(s_Ai, s_Ak),IEOR(s_Aj, s_b)) == IEOR(IEOR(s_c, s_a),IEOR(s_e, s_f)) .and. &
IEOR(IEOR(s_a, s_b),s_c) == IEOR(IEOR(s_e, s_f),s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E4_(s_f, s_e, s_a, s_c, s_b, s_Aj)%array(i_f, i_e, i_a, i_c, i_b, i_Aj) & 
  * X_(s_f, s_e, s_c, s_b, s_a)%array(i_f, i_e, i_c, i_b, i_a)

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
! External_count : 3
! Total_count    : 9


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no1_x28


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x29(sVa, iVa, Ecas, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no0_x29(sVa, iVa, Ecas, av2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x29


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x29(s_Va, i_Va, Ecas, T2_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Ak, i_Ak

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) Ecas T2(a,b,Ak,Va) E2(Ai,a,Aj,b) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Ak, s_Va) .and. &
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * Ecas & 
  * T2_(s_Ak, s_b, s_a)%array(i_Ak, i_b, i_a) & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai)

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
! Loop_count     : 5
! External_count : 1
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_ooov_ooov_no0_x29


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_ooov_ooov_no0_x30(sVa, iVa, Ecas, T2, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: Ecas
real(kind=8), intent(inout) :: T2(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_av2_2(sVa, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_ooov_ooov_no0_x30(sVa, iVa, Ecas, av2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i2)
deallocate(av2_i)

end subroutine g_if_sigma_ooov_ooov_no0_x30


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:44

subroutine g_sigma_ooov_ooov_no0_x30(s_Va, i_Va, Ecas, T2_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Declaration of numerical constants .... 
real(kind=8), intent(inout) :: Ecas
! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Ai, i_Ai, s_Ak, i_Ak
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Ak,Va) <-- 
! (   1.00000) Ecas T2(a,b,c,Va) E3(Ai,Ak,Aj,b,c,a) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Ai = 0, nir-1
do s_Ak = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Ak, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(IEOR(s_Ai, s_Ak),s_Aj) == IEOR(IEOR(s_b, s_c),s_a)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Ak = psym(I_BEGIN, I_O, s_Ak), psym(I_END, I_O, s_Ak)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) = S2_(s_Ak, s_Aj, s_Ai)%array(i_Ak, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * Ecas & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * E3_(s_a, s_c, s_b, s_Aj, s_Ak, s_Ai)%array(i_a, i_c, i_b, i_Aj, i_Ak, i_Ai)

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

end subroutine g_sigma_ooov_ooov_no0_x30


! -----------------------------------------------------------------------
