#include "../f_ct.fh"



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x0(sVb, iVb, V, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: V(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x0(sVb, iVb, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x0


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x0(s_Vb, i_Vb, V_, Hdiag_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_Va, i_Va

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) V(Vb,Va,Va,Vb) E2(Ai,Aj,Aj,Ai) 
! where i_Vb in {vir}
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Vb, s_Va) == IEOR(s_Va, s_Vb) .and. &
IEOR(s_Ai, s_Aj) == IEOR(s_Aj, s_Ai)) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_Vb, s_Va, s_Va)%array(i_Vb, i_Va, i_Va) & 
  * E2_(s_Ai, s_Aj, s_Aj, s_Ai)%array(i_Ai, i_Aj, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x1(sVb, iVb, V, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: V(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x1(sVb, iVb, h2_i, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x1


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x1(s_Vb, i_Vb, V_, Hdiag_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_Va, i_Va

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) V(Vb,Vb,Va,Va) E2(Ai,Ai,Aj,Aj) 
! where i_Vb in {vir}
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(s_Vb, s_Vb) == IEOR(s_Va, s_Va) .and. &
IEOR(s_Ai, s_Ai) == IEOR(s_Aj, s_Aj)) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_Va, s_Va, s_Vb)%array(i_Va, i_Va, i_Vb) & 
  * E2_(s_Aj, s_Aj, s_Ai, s_Ai)%array(i_Aj, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x1


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x2(h, X, nir, nsym, psym)

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
call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call g_diag_oovv_oovv_no0_x2(h1, xaav, d2, nir, nsym, psym)

deallocate(h1)
deallocate(xaav)

end subroutine g_if_diag_oovv_oovv_no0_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x2(h_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_Va, i_Va

!X(Ai,Aj,Va) <-- 
! (   1.00000) h(Va,Va) E2(Ai,Ai,Aj,Aj) 
do s_Va = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0 .and. & 
IEOR(S_Va, s_Va) == 0 .and. &
IEOR(s_Ai, s_Ai) == IEOR(s_Aj, s_Aj)) then
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = X_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * h_(s_Va, s_Va)%array(i_Va, i_Va) & 
  * E2_(s_Aj, s_Aj, s_Ai, s_Ai)%array(i_Aj, i_Aj, i_Ai, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 0
! Total_count    : 3


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no0_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x2(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaav(sleft, x, nir, nsym, psym) ! -> xaav (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x2(sVb, iVb, xaav, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaav)

end subroutine g_if_diag_oovv_oovv_no1_x2


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no1_x2(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: X_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj, s_Va, i_Va

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Va) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x2


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x3(sVb, iVb, h, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: h(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x3(sVb, iVb, h1, xaa, d2, nir, nsym, psym)

deallocate(h1)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no0_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x3(s_Vb, i_Vb, h_, X_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!X(Ai,Aj,Vb) <-- 
! (   1.00000) h(Vb,Vb) E2(Ai,Ai,Aj,Aj) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0 .and. & 
IEOR(S_Vb, s_Vb) == 0 .and. &
IEOR(s_Ai, s_Ai) == IEOR(s_Aj, s_Aj)) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * h_(s_Vb, s_Vb)%array(i_Vb, i_Vb) & 
  * E2_(s_Aj, s_Aj, s_Ai, s_Ai)%array(i_Aj, i_Aj, i_Ai, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 2
! External_count : 1
! Total_count    : 3


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no0_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x3(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x3(sVb, iVb, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x3


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no1_x3(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Vb) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x3


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x4(sVa, iVa, h, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: h(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h1(h, nir, nsym, psym) ! -> h1 (allocate)
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x4(sVa, iVa, h1, av2_i2, d2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h1)

end subroutine g_if_diag_oovv_oovv_no0_x4


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x4(s_Va, i_Va, h_, Hdiag_, E2_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)
type(symblock4), intent(inout) :: E2_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   2.00000) h(Va,Va) E2(Ai,Aj,Aj,Ai) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(S_Va, s_Va) == 0 .and. &
IEOR(s_Ai, s_Aj) == IEOR(s_Aj, s_Ai)) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 2.0d+00 & 
  * h_(s_Va, s_Va)%array(i_Va, i_Va) & 
  * E2_(s_Ai, s_Aj, s_Aj, s_Ai)%array(i_Ai, i_Aj, i_Aj, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 2
! External_count : 1
! Total_count    : 3


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no0_x4


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x5(sVa, iVa, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sVa

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x5(sVa, iVa, h2_i, xaa, d3, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x5


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x5(s_Va, i_Va, V_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj,Va) <-- 
! (   1.00000) V(Va,a,Va,b) E3(Ai,b,Aj,Aj,a,Ai) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0 .and. & 
IEOR(s_Va, s_a) == IEOR(s_Va, s_b) .and. &
IEOR(IEOR(s_Ai, s_b),s_Aj) == IEOR(IEOR(s_Aj, s_a),s_Ai)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_Va, s_a)%array(i_b, i_Va, i_a) & 
  * E3_(s_Ai, s_a, s_Aj, s_Aj, s_b, s_Ai)%array(i_Ai, i_a, i_Aj, i_Aj, i_b, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x5(sVb, iVb, sVa, iVa, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x5(sVb, iVb, sVa, iVa, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x5


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no1_x5(s_Vb, i_Vb, s_Va, i_Va, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Va) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 2
! Total_count    : 4


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no1_x5


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x6(sVa, iVa, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sVa

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x6(sVa, iVa, h2_i, xaa, d3, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x6


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no0_x6(s_Va, i_Va, V_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj,Va) <-- 
! (   1.00000) V(Va,Va,a,b) E3(Ai,Ai,Aj,Aj,a,b) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0 .and. & 
IEOR(s_Va, s_Va) == IEOR(s_a, s_b) .and. &
IEOR(IEOR(s_Ai, s_Ai),s_Aj) == IEOR(IEOR(s_Aj, s_a),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_a, s_Va)%array(i_b, i_a, i_Va) & 
  * E3_(s_b, s_a, s_Aj, s_Aj, s_Ai, s_Ai)%array(i_b, i_a, i_Aj, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x6(sVb, iVb, sVa, iVa, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb, sVa, iVa
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVa

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x6(sVb, iVb, sVa, iVa, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x6


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:57

subroutine g_diag_oovv_oovv_no1_x6(s_Vb, i_Vb, s_Va, i_Va, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Va) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Va) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 3
! External_count : 2
! Total_count    : 4


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no1_x6


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x7(sVb, iVb, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x7(sVb, iVb, h2_i, xaa, d3, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x7


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x7(s_Vb, i_Vb, V_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj,Vb) <-- 
! (   1.00000) V(Vb,a,Vb,b) E3(Ai,Ai,Aj,b,a,Aj) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0 .and. & 
IEOR(s_Vb, s_a) == IEOR(s_Vb, s_b) .and. &
IEOR(IEOR(s_Ai, s_Ai),s_Aj) == IEOR(IEOR(s_b, s_a),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_Vb, s_a)%array(i_b, i_Vb, i_a) & 
  * E3_(s_Aj, s_a, s_b, s_Aj, s_Ai, s_Ai)%array(i_Aj, i_a, i_b, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x7(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x7(sVb, iVb, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x7


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x7(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Vb) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x7


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x8(sVb, iVb, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVb, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x8(sVb, iVb, h2_i, xaa, d3, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x8


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x8(s_Vb, i_Vb, V_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj,Vb) <-- 
! (   1.00000) V(Vb,Vb,a,b) E3(Ai,Ai,Aj,Aj,a,b) 
! where i_Vb in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0 .and. & 
IEOR(s_Vb, s_Vb) == IEOR(s_a, s_b) .and. &
IEOR(IEOR(s_Ai, s_Ai),s_Aj) == IEOR(IEOR(s_Aj, s_a),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_a, s_Vb)%array(i_b, i_a, i_Vb) & 
  * E3_(s_b, s_a, s_Aj, s_Aj, s_Ai, s_Ai)%array(i_b, i_a, i_Aj, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x8(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = sVb

call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x8(sVb, iVb, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x8


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x8(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Vb, s_Vb
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj,Vb) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(IEOR(s_Ai, s_Aj),s_Vb) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x8


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x9(h, X, nir, nsym, psym)

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
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x9(h1, xaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no0_x9


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x9(h_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj) <-- 
! (   1.00000) h(a,b) E3(Ai,Ai,Aj,Aj,a,b) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(S_Ai, s_Aj) == 0 .and. & 
IEOR(S_a, s_b) == 0 .and. &
IEOR(IEOR(s_Ai, s_Ai),s_Aj) == IEOR(IEOR(s_Aj, s_a),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * h_(s_b, s_a)%array(i_b, i_a) & 
  * E3_(s_b, s_a, s_Aj, s_Aj, s_Ai, s_Ai)%array(i_b, i_a, i_Aj, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x9(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x9(sVb, iVb, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x9


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x9(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   1.00000) X(Ai,Aj) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(S_Ai, s_Aj) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x9


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x10(sVa, iVa, V, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: V(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x10(sVa, iVa, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x10


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x10(s_Va, i_Va, V_, Hdiag_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   1.00000) V(Va,a,Va,b) E3(Ai,Aj,Aj,b,a,Ai) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(s_Va, s_a) == IEOR(s_Va, s_b) .and. &
IEOR(IEOR(s_Ai, s_Aj),s_Aj) == IEOR(IEOR(s_b, s_a),s_Ai)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_Va, s_a)%array(i_b, i_Va, i_a) & 
  * E3_(s_Ai, s_a, s_b, s_Aj, s_Aj, s_Ai)%array(i_Ai, i_a, i_b, i_Aj, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x10


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x11(sVa, iVa, V, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: V(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x11(sVa, iVa, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x11


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x11(s_Va, i_Va, V_, Hdiag_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   1.00000) V(Va,a,Va,b) E3(Ai,a,Aj,Ai,b,Aj) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(s_Va, s_a) == IEOR(s_Va, s_b) .and. &
IEOR(IEOR(s_Ai, s_a),s_Aj) == IEOR(IEOR(s_Ai, s_b),s_Aj)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_b, s_Va, s_a)%array(i_b, i_Va, i_a) & 
  * E3_(s_Aj, s_b, s_Ai, s_Aj, s_a, s_Ai)%array(i_Aj, i_b, i_Ai, i_Aj, i_a, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x11


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x12(sVa, iVa, V, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: V(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sVa, V, nir, nsym, psym) ! -> h2_i (allocate)
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no0_x12(sVa, iVa, h2_i, av2_i2, d3, nir, nsym, psym)

deallocate(av2_i2)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x12


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x12(s_Va, i_Va, V_, Hdiag_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_Va, s_Va
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   2.00000) V(Va,Va,a,b) E3(Ai,Aj,Aj,Ai,a,b) 
! where i_Va in {vir}
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(s_Va, s_Va) == IEOR(s_a, s_b) .and. &
IEOR(IEOR(s_Ai, s_Aj),s_Aj) == IEOR(IEOR(s_Ai, s_a),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 2.0d+00 & 
  * V_(s_b, s_a, s_Va)%array(i_b, i_a, i_Va) & 
  * E3_(s_b, s_a, s_Ai, s_Aj, s_Aj, s_Ai)%array(i_b, i_a, i_Ai, i_Aj, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x12


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x13(h, X, nir, nsym, psym)

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
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x13(h1, xaa, d3, nir, nsym, psym)

deallocate(h1)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no0_x13


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x13(h_, X_, E3_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E3_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: h_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_a, i_a, s_b, i_b, s_Aj, i_Aj

!X(Ai,Aj) <-- 
! (   1.00000) h(a,b) E3(Ai,Aj,Aj,Ai,a,b) 
do s_a = 0, nir-1
do s_b = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(S_Ai, s_Aj) == 0 .and. & 
IEOR(S_a, s_b) == 0 .and. &
IEOR(IEOR(s_Ai, s_Aj),s_Aj) == IEOR(IEOR(s_Ai, s_a),s_b)) then
do i_a = psym(I_BEGIN, I_O, s_a), psym(I_END, I_O, s_a)
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * h_(s_b, s_a)%array(i_b, i_a) & 
  * E3_(s_b, s_a, s_Ai, s_Aj, s_Aj, s_Ai)%array(i_b, i_a, i_Ai, i_Aj, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x13(sVa, iVa, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x13(sVa, iVa, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x13


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x13(s_Va, i_Va, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   1.00000) X(Ai,Aj) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(S_Ai, s_Aj) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 1.0d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 2
! External_count : 1
! Total_count    : 3


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no1_x13


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x14(sa, ia, sc, ic, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x14(sa, ia, sc, ic, h2_i, xaa, d4_ij, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x14


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x14(s_a, i_a, s_c, i_c, V_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_b, i_b, s_Aj, i_Aj, s_d, i_d

!X(Ai,Aj) <-- 
! (   1.00000) V(a,c,b,d) E4(a,c,Ai,Ai,Aj,Aj,b,d) 
! where i_a in {core}
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(S_Ai, s_Aj) == 0 .and. & 
IEOR(s_a, s_c) == IEOR(s_b, s_d) .and. &
IEOR(IEOR(s_a, s_c),IEOR(s_Ai, s_Ai)) == IEOR(IEOR(s_Aj, s_Aj),IEOR(s_b, s_d))) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_d, s_b, s_c)%array(i_d, i_b, i_c) & 
  * E4_(s_d, s_b, s_Aj, s_Aj, s_Ai, s_Ai)%array(i_d, i_b, i_Aj, i_Aj, i_Ai, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x14(sVb, iVb, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVb, iVb
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVb, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x14(sVb, iVb, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x14


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x14(s_Vb, i_Vb, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Vb) <-- 
! (   0.50000) X(Ai,Aj) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
do s_Va = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Vb) .and. & 
IEOR(S_Ai, s_Aj) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)
do i_Va = psym(I_BEGIN, I_V, s_Va), psym(I_END, I_V, s_Va)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 0.5d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no1_x14


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no0_x15(sa, ia, sc, ic, V, X, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sa, ia, sc, ic
real(kind=8), intent(inout) :: V(*), X(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

call set_symblock_h2(sa, V, nir, nsym, psym) ! -> h2_i (allocate)
sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call g_diag_oovv_oovv_no0_x15(sa, ia, sc, ic, h2_i, xaa, d4_ij, nir, nsym, psym)

deallocate(xaa)
deallocate(h2_i)

end subroutine g_if_diag_oovv_oovv_no0_x15


! Converting terms with V[i*|**] with 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no0_x15(s_a, i_a, s_c, i_c, V_, X_, E4_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

integer, intent(in) :: i_a, s_a
integer, intent(in) :: i_c, s_c
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)
type(symblock6), intent(inout) :: E4_(0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1, 0:nir-1)
type(symblock3), intent(inout) :: V_(0:nir-1, 0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_b, i_b, s_Aj, i_Aj, s_d, i_d

!X(Ai,Aj) <-- 
! (   1.00000) V(a,c,b,d) E4(a,c,Ai,Aj,Aj,Ai,b,d) 
! where i_a in {core}
do s_b = 0, nir-1
do s_d = 0, nir-1
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(S_Ai, s_Aj) == 0 .and. & 
IEOR(s_a, s_c) == IEOR(s_b, s_d) .and. &
IEOR(IEOR(s_a, s_c),IEOR(s_Ai, s_Aj)) == IEOR(IEOR(s_Aj, s_Ai),IEOR(s_b, s_d))) then
do i_b = psym(I_BEGIN, I_O, s_b), psym(I_END, I_O, s_b)
do i_d = psym(I_BEGIN, I_O, s_d), psym(I_END, I_O, s_d)
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) = X_(s_Aj, s_Ai)%array(i_Aj, i_Ai) &
  + 1.0d+00 & 
  * V_(s_d, s_b, s_c)%array(i_d, i_b, i_c) & 
  * E4_(s_d, s_b, s_Ai, s_Aj, s_Aj, s_Ai)%array(i_d, i_b, i_Ai, i_Aj, i_Aj, i_Ai)

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

end subroutine g_diag_oovv_oovv_no0_x15


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_diag_oovv_oovv_no1_x15(sVa, iVa, X, Hdiag, nir, nsym, psym)

use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
use mod_g_if_comm08
implicit none

integer, intent(inout) :: sVa, iVa
real(kind=8), intent(inout) :: X(*), Hdiag(*)
! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)
! Some extra stuff
integer :: sleft

sleft = 0
call set_symblock_xaa(sleft, x, nir, nsym, psym) ! -> xaa (allocate) 
call set_symblock_av2_2(sVa, Hdiag, nir, nsym, psym) ! -> av2_i2 (allocate)
call g_diag_oovv_oovv_no1_x15(sVa, iVa, xaa, av2_i2, nir, nsym, psym)

deallocate(av2_i2)
deallocate(xaa)

end subroutine g_if_diag_oovv_oovv_no1_x15


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... True
! LHS name ............................... ham
! Name of ERI ............................ V
! Name of BareAmpPack appearing in RHS.... T
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:2:58

subroutine g_diag_oovv_oovv_no1_x15(s_Va, i_Va, X_, Hdiag_, nir, nsym, psym)

! FEMTO BEGIN  **************************************************************
use comm8_module, only : symblock1, symblock2, symblock3, symblock4, symblock5, symblock6
implicit none

! Information of the Irreps ....
integer, intent(inout) :: nir, nsym(3,0:nir-1), psym(3,6,0:nir-1)

! Declare tensors used ...
type(symblock3), intent(inout) :: Hdiag_(0:nir-1, 0:nir-1, 0:nir-1)
type(symblock2), intent(inout) :: X_(0:nir-1, 0:nir-1)

! Indices used in the contractions as dummy ... 

integer :: s_Ai, i_Ai, s_Aj, i_Aj

!Hdiag(Ai,Aj,Va,Va) <-- 
! (   0.50000) X(Ai,Aj) 
do s_Ai = 0, nir-1
do s_Aj = 0, nir-1
if( &
IEOR(s_Ai, s_Aj) == IEOR(s_Va, s_Va) .and. & 
IEOR(S_Ai, s_Aj) == 0) then
do i_Ai = psym(I_BEGIN, I_O, s_Ai), psym(I_END, I_O, s_Ai)
do i_Aj = psym(I_BEGIN, I_O, s_Aj), psym(I_END, I_O, s_Aj)

Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) = Hdiag_(s_Va, s_Aj, s_Ai)%array(i_Va, i_Aj, i_Ai) &
  + 0.5d+00 & 
  * X_(s_Aj, s_Ai)%array(i_Aj, i_Ai)

end do ! Irrep loop
end do ! Irrep loop
end if ! Irrep Cond
end do ! Orbital loop
end do ! Orbital loop
! Loop_count     : 2
! External_count : 1
! Total_count    : 3


! FEMTO END    **************************************************************

end subroutine g_diag_oovv_oovv_no1_x15


! -----------------------------------------------------------------------
