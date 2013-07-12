#include "../f_ct.fh"



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x0(sd, id, sVb, iVb, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVb, iVb
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x0(sd, id, sVb, iVb, h2_i, av2_i, xaav, nir, nsym, psym)

deallocate(xaav)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x0


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x0(s_d, i_d, s_Vb, i_Vb, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Va, i_Va

!X(a,b,Va,Vb) <-- 
! (   1.00000) V2(Vb,d,Va,c) T2(a,b,c,d) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Vb, s_d) == IEOR(s_Va, s_c) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) = X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_c, s_Va, s_d)%array(i_c, i_Va, i_d) & 
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

end subroutine g_sigma_oovv_oovv_no0_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x0(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x0(sVb, iVb, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x0


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x0(s_Vb, i_Vb, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x1(sd, id, sVb, iVb, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sd, id, sVb, iVb
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sd, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x1(sd, id, sVb, iVb, h2_i, av2_i, xaav, nir, nsym, psym)

deallocate(xaav)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x1


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x1(s_d, i_d, s_Vb, i_Vb, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_d, s_d
integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Va, i_Va

!X(a,b,Va,Vb) <-- 
! (   1.00000) V2(Vb,c,Va,d) T2(a,b,c,d) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Vb, s_c) == IEOR(s_Va, s_d) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) = X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_d, s_Va, s_c)%array(i_d, i_Va, i_c) & 
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

end subroutine g_sigma_oovv_oovv_no0_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x1(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x1(sVb, iVb, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x1


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x1(s_Vb, i_Vb, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,b,Aj,a) X(a,b,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_b) == IEOR(s_Aj, s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_a, s_Aj, s_b, s_Ai)%array(i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x2(sc, ic, sVb, iVb, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVb, iVb
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x2(sc, ic, sVb, iVb, av2_i, h1, xaav, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x2(s_c, i_c, s_Vb, i_Vb, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_Va, i_Va

!X(a,b,Vb,Va) <-- 
! (   1.00000) T2(a,b,Va,c) h(Vb,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_Va, s_c) .and. &
IEOR(S_Vb, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) = X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) & 
  * h_(s_c, s_Vb)%array(i_c, i_Vb)

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

end subroutine g_sigma_oovv_oovv_no0_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x2(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x2(sVb, iVb, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x2(s_Vb, i_Vb, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Vb,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x3(sVb, iVb, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x3(sVb, iVb, av2_i, h1, xaav, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x3(s_Vb, i_Vb, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Va, i_Va

!X(a,b,Vb,Va) <-- 
! (   1.00000) T2(b,a,c,Vb) h(Va,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va) .and. & 
IEOR(s_b, s_a) == IEOR(s_c, s_Vb) .and. &
IEOR(S_Va, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) = X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_a, s_b)%array(i_c, i_a, i_b) & 
  * h_(s_c, s_Va)%array(i_c, i_Va)

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

end subroutine g_sigma_oovv_oovv_no0_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x3(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x3(sVb, iVb, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x3(s_Vb, i_Vb, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,b,Aj,a) X(a,b,Vb,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_b) == IEOR(s_Aj, s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_a, s_Aj, s_b, s_Ai)%array(i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x4(sVb, iVb, sVa, iVa, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = IEOR(sVb,sVa)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_sigma_oovv_oovv_no0_x4(sVb, iVb, sVa, iVa, av2_i, h1, xaa, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaa)

end subroutine g_if_sigma_oovv_oovv_no0_x4


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x4(s_Vb, i_Vb, s_Va, i_Va, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c

!X(a,b,Vb,Va) <-- 
! (   1.00000) T2(a,b,c,Va) h(Vb,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Va) .and. &
IEOR(S_Vb, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)

X_(s_b, s_a)%array(i_b, i_a) = X_(s_b, s_a)%array(i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_c, s_Vb)%array(i_c, i_Vb)

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

end subroutine g_sigma_oovv_oovv_no0_x4


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x4(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVb,sVa)

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x4(sVb, iVb, sVa, iVa, xaa, av2_i2, d2, nir, nsym, psym)

deallocate(xaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x4


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x4(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,b,Aj,a) X(a,b,Vb,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_b) == IEOR(s_Aj, s_a) .and. &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_a, s_Aj, s_b, s_Ai)%array(i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_b, s_a)%array(i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x4


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x5(sVb, iVb, T2, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_sigma_oovv_oovv_no0_x5(sVb, iVb, av2_i, h1, xaav, nir, nsym, psym)

deallocate(h1)
deallocate(av2_i)
deallocate(xaav)

end subroutine g_if_sigma_oovv_oovv_no0_x5


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x5(s_Vb, i_Vb, T2_, h_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_Va, i_Va

!X(a,b,Vb,Va) <-- 
! (   1.00000) T2(a,b,c,Vb) h(Va,c) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va) .and. & 
IEOR(s_a, s_b) == IEOR(s_c, s_Vb) .and. &
IEOR(S_Va, s_c) == 0) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) = X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) &
  + 1.0d+00 & 
  * T2_(s_c, s_b, s_a)%array(i_c, i_b, i_a) & 
  * h_(s_c, s_Va)%array(i_c, i_Va)

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

end subroutine g_sigma_oovv_oovv_no0_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x5(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x5(sVb, iVb, xaav, av2_i2, d2, nir, nsym, psym)

deallocate(xaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x5


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x5(s_Vb, i_Vb, X_, S2_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E2(Ai,a,Aj,b) X(a,b,Vb,Va) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Ai, s_a) == IEOR(s_Aj, s_b) .and. &
IEOR(s_a, s_b) == IEOR(s_Vb, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E2_(s_b, s_Aj, s_a, s_Ai)%array(i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x6(sc, ic, sVb, iVb, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVb, iVb
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVb

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x6(sc, ic, sVb, iVb, h2_i, av2_i, xaaaav, nir, nsym, psym)

deallocate(xaaaav)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x6


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x6(s_c, i_c, s_Vb, i_Vb, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Va, i_Va

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Vb,d,c,e) T2(a,b,Va,c) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Vb, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_Va, s_c)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a) = X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * T2_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no0_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x6(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x6(sVb, iVb, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x6


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x6(s_Vb, i_Vb, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,a,Aj,d,e,b) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_d, s_e),s_b) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_b, s_e, s_d, s_Aj, s_a, s_Ai)%array(i_b, i_e, i_d, i_Aj, i_a, i_Ai) & 
  * X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x7(sc, ic, sVb, iVb, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, sVb, iVb
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sc, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = sVb

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call g_sigma_oovv_oovv_no0_x7(sc, ic, sVb, iVb, h2_i, av2_i, xaaaav, nir, nsym, psym)

deallocate(xaaaav)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x7


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x7(s_c, i_c, s_Vb, i_Vb, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Va, i_Va

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Vb,c,d,e) T2(a,b,Va,c) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Vb, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_Va, s_c)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a) = X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
  * T2_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no0_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x7(sVb, iVb, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaaaav(sleft, x, nir, nsym, psym) ! -> xaaaav (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x7(sVb, iVb, xaaaav, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaav)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x7


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x7(s_Vb, i_Vb, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock5), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,a,Aj,b,d,e) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_b, s_d),s_e) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_e, s_d, s_b, s_Aj, s_a, s_Ai)%array(i_e, i_d, i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_Va, s_e, s_d, s_b, s_a)%array(i_Va, i_e, i_d, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x8(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x8(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x8


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x8(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Va,d,c,e) T2(b,a,c,Vb) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Va, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_b, s_a) == IEOR(s_c, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
  * T2_(s_c, s_a, s_b)%array(i_c, i_a, i_b)

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

end subroutine g_sigma_oovv_oovv_no0_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x8(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x8(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x8


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x8(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,d,Aj,a,e,b) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_a, s_e),s_b) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_b, s_e, s_a, s_Aj, s_d, s_Ai)%array(i_b, i_e, i_a, i_Aj, i_d, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x9(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x9(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x9


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x9(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Va,c,d,e) T2(b,a,c,Vb) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Va, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_b, s_a) == IEOR(s_c, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
  * T2_(s_c, s_a, s_b)%array(i_c, i_a, i_b)

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

end subroutine g_sigma_oovv_oovv_no0_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x9(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x9(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x9


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x9(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,b,Aj,a,d,e) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_a, s_d),s_e) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_e, s_d, s_a, s_Aj, s_b, s_Ai)%array(i_e, i_d, i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x10(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x10(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x10


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x10(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Vb,d,c,e) T2(a,b,c,Va) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Vb, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
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

end subroutine g_sigma_oovv_oovv_no0_x10


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x10(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x10(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x10


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x10(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,b,Aj,d,e,a) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_d, s_e),s_a) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_e, s_d, s_Aj, s_b, s_Ai)%array(i_a, i_e, i_d, i_Aj, i_b, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x10


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x11(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVa, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x11(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x11


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x11(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Vb,c,d,e) T2(a,b,c,Va) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Vb, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Va)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
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

end subroutine g_sigma_oovv_oovv_no0_x11


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x11(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x11(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x11


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x11(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,b,Aj,a,d,e) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_a, s_d),s_e) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_e, s_d, s_a, s_Aj, s_b, s_Ai)%array(i_e, i_d, i_a, i_Aj, i_b, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x11


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x12(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x12(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x12


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x12(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Va,d,c,e) T2(a,b,c,Vb) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Va, s_d) == IEOR(s_c, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_c, s_d)%array(i_e, i_c, i_d) & 
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

end subroutine g_sigma_oovv_oovv_no0_x12


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x12(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x12(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x12


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x12(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,d,Aj,b,e,a) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_d),s_Aj) == IEOR(IEOR(s_b, s_e),s_a) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_a, s_e, s_b, s_Aj, s_d, s_Ai)%array(i_a, i_e, i_b, i_Aj, i_d, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x12


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x13(sVb, iVb, sVa, iVa, V2, T2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: V2(*), T2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V2, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x13(sVb, iVb, sVa, iVa, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x13


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no0_x13(s_Vb, i_Vb, s_Va, i_Va, V2_, T2_, X_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_e, i_e

!X(a,b,d,e,Va,Vb) <-- 
! (   1.00000) V2(Va,c,d,e) T2(a,b,c,Vb) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
if( &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb) .and. & 
IEOR(s_Va, s_c) == IEOR(s_d, s_e) .and. &
IEOR(s_a, s_b) == IEOR(s_c, s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_V, s_c), psym(I_END, I_V, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)

X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) = X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_e, s_d, s_c)%array(i_e, i_d, i_c) & 
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

end subroutine g_sigma_oovv_oovv_no0_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x13(sVb, iVb, sVa, iVa, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = IEOR(sVa,sVb)

call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x13(sVb, iVb, sVa, iVa, xaaaa, av2_i2, d3, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x13


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:3

subroutine g_sigma_oovv_oovv_no1_x13(s_Vb, i_Vb, s_Va, i_Va, X_, S2_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_e, i_e, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) E3(Ai,a,Aj,b,d,e) X(a,b,d,e,Va,Vb) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_e = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_b, s_d),s_e) .and. &
IEOR(IEOR(s_a, s_b),s_d) == IEOR(IEOR(s_e, s_Va),s_Vb)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_e = psym(I_BEGIN, I_O, s_e), psym(I_END, I_O, s_e)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * E3_(s_e, s_d, s_b, s_Aj, s_a, s_Ai)%array(i_e, i_d, i_b, i_Aj, i_a, i_Ai) & 
  * X_(s_e, s_d, s_b, s_a)%array(i_e, i_d, i_b, i_a)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no1_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x14(h, X, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x14(h1, xaaaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x14


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no0_x14(h_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) h(c,d) E3(Ai,a,Aj,b,c,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(S_c, s_d) == 0 .and. &
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_b, s_c),s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) = X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) &
  + 1.0d+00 & 
  * h_(s_d, s_c)%array(i_d, i_c) & 
  * E3_(s_d, s_c, s_b, s_Aj, s_a, s_Ai)%array(i_d, i_c, i_b, i_Aj, i_a, i_Ai)

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
! External_count : 0
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no0_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x14(sVb, iVb, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x14(sVb, iVb, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x14


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no1_x14(s_Vb, i_Vb, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) T2(a,b,Va,Vb) X(a,b,Ai,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * T2_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x15(h, X, nir, nsym, psym)

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
call g_sigma_oovv_oovv_no0_x15(h1, xaaaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaaaa)

end subroutine g_if_sigma_oovv_oovv_no0_x15


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no0_x15(h_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_c, i_c, s_d, i_d, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) h(c,d) E3(Ai,b,Aj,a,c,d) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_c = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(S_c, s_d) == 0 .and. &
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_a, s_c),s_d)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_c = psym(I_BEGIN, I_O, s_c), psym(I_END, I_O, s_c)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) = X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) &
  + 1.0d+00 & 
  * h_(s_d, s_c)%array(i_d, i_c) & 
  * E3_(s_d, s_c, s_a, s_Aj, s_b, s_Ai)%array(i_d, i_c, i_a, i_Aj, i_b, i_Ai)

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
! External_count : 0
! Total_count    : 6


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no0_x15


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x15(sVb, iVb, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x15(sVb, iVb, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x15


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no1_x15(s_Vb, i_Vb, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   1.00000) T2(b,a,Va,Vb) X(a,b,Ai,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_b, s_a) == IEOR(s_Va, s_Vb) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * T2_(s_Va, s_a, s_b)%array(i_Va, i_a, i_b) & 
  * X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x15


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x16(sc, ic, se, ie, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, se, ie
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x16(sc, ic, se, ie, h2_i, xaaaa, d4_ij, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x16


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no0_x16(s_c, i_c, s_e, i_e, V2_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_e, s_e
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_f, i_f, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) V2(c,e,d,f) E4(c,e,Ai,a,Aj,b,d,f) 
! where i_c in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_f = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(s_c, s_e) == IEOR(s_d, s_f) .and. &
IEOR(IEOR(s_c, s_e),IEOR(s_Ai, s_a)) == IEOR(IEOR(s_Aj, s_b),IEOR(s_d, s_f))) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) = X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_e)%array(i_f, i_d, i_e) & 
  * E4_(s_f, s_d, s_b, s_Aj, s_a, s_Ai)%array(i_f, i_d, i_b, i_Aj, i_a, i_Ai)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no0_x16


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x16(sVb, iVb, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x16(sVb, iVb, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x16


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no1_x16(s_Vb, i_Vb, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   0.50000) T2(a,b,Va,Vb) X(a,b,Ai,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_a, s_b) == IEOR(s_Va, s_Vb) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 0.5d+00 & 
  * T2_(s_Va, s_b, s_a)%array(i_Va, i_b, i_a) & 
  * X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x16


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no0_x17(sc, ic, se, ie, V2, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sc, ic, se, ie
real(kind=8), intent(inout) :: V2(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sc, V2, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call g_sigma_oovv_oovv_no0_x17(sc, ic, se, ie, h2_i, xaaaa, d4_ij, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)

end subroutine g_if_sigma_oovv_oovv_no0_x17


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no0_x17(s_c, i_c, s_e, i_e, V2_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_c, s_c
integer, intent(in) :: i_e, s_e
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_a, i_a, s_b, i_b, s_d, i_d, s_f, i_f, s_Ai, i_Ai
integer :: s_Aj, i_Aj

!X(a,b,Ai,Aj) <-- 
! (   1.00000) V2(c,e,d,f) E4(c,e,Ai,b,Aj,a,d,f) 
! where i_c in {core}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_f = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj) .and. & 
IEOR(s_c, s_e) == IEOR(s_d, s_f) .and. &
IEOR(IEOR(s_c, s_e),IEOR(s_Ai, s_b)) == IEOR(IEOR(s_Aj, s_a),IEOR(s_d, s_f))) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_f = psym(I_BEGIN, I_O, s_f), psym(I_END, I_O, s_f)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) = X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a) &
  + 1.0d+00 & 
  * V2_(s_f, s_d, s_e)%array(i_f, i_d, i_e) & 
  * E4_(s_f, s_d, s_a, s_Aj, s_b, s_Ai)%array(i_f, i_d, i_a, i_Aj, i_b, i_Ai)

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
! External_count : 2
! Total_count    : 8


! FEMTO END    **************************************************************

end subroutine g_sigma_oovv_oovv_no0_x17


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_oovv_oovv_no1_x17(sVb, iVb, T2, X, S2, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: T2(*), X(*), S2(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_av2(sVb, T2, nir, nsym, psym) ! -> av2_i (allocate)
sleft = 0
call set_symblock_xaaaa(sleft, x, nir, nsym, psym) ! -> xaaaa (allocate) 
call set_symblock_av2_2(sVb, S2, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_sigma_oovv_oovv_no1_x17(sVb, iVb, av2_i, xaaaa, av2_i2, nir, nsym, psym)

deallocate(xaaaa)
deallocate(av2_i)
deallocate(av2_i2)

end subroutine g_if_sigma_oovv_oovv_no1_x17


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 14:53:4

subroutine g_sigma_oovv_oovv_no1_x17(s_Vb, i_Vb, T2_, X_, S2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: S2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: T2_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj, s_Va, i_Va

!S2(Ai,Aj,Va,Vb) <-- 
! (   0.50000) T2(b,a,Va,Vb) X(a,b,Ai,Aj) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_b, s_a) == IEOR(s_Va, s_Vb) .and. &
IEOR(s_a, s_b) == IEOR(s_Ai, s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = S2_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 0.5d+00 & 
  * T2_(s_Va, s_a, s_b)%array(i_Va, i_a, i_b) & 
  * X_(s_Aj, s_Ai, s_b, s_a)%array(i_Aj, i_Ai, i_b, i_a)

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

end subroutine g_sigma_oovv_oovv_no1_x17


! -----------------------------------------------------------------------
