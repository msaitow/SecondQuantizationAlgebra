#include "../f_ct.fh"



! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_oovv_no0_x0(sc, ic, sd, id, V2, T2, X, nir, nsym, psym)

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
call g_sigma_G_oovv_no0_x0(sc, ic, sd, id, h2_i, av2_i, xaaaa, nir, nsym, psym)

deallocate(xaaaa)
deallocate(h2_i)
deallocate(av2_i)

end subroutine g_if_sigma_G_oovv_no0_x0


! Converting terms with V[a*|**] and no 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:4:58

subroutine g_sigma_G_oovv_no0_x0(s_c, i_c, s_d, i_d, V2_, T2_, X_, nir, nsym, psym)

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
! where i_c in {vir}
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

end subroutine g_sigma_G_oovv_no0_x0


! -----------------------------------------------------------------------


! **********************************************************************
!                                                                       
! **********************************************************************
subroutine g_if_sigma_G_oovv_no1_x0(X, S0, nir, nsym, psym)

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
call g_sigma_G_oovv_no1_x0(xaaaa, S0, d2, nir, nsym, psym)

deallocate(xaaaa)

end subroutine g_if_sigma_G_oovv_no1_x0


! Converting terms with neither ERI, nor 4-RDMs ..... 


! ---------------------------- Parameters used ------------------------------
!                                                                            
! Whether the LHS is a BareAmpPack ....... False
! LHS name ............................... sig
! Name of ERI ............................ V2
! Name of BareAmpPack appearing in RHS.... T2
!                                                                            
!----------------------------------------------------------------------------

! Generated time: 2013/7/12 15:4:58

subroutine g_sigma_G_oovv_no1_x0(X_, S0_, E2_, nir, nsym, psym)

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

end subroutine g_sigma_G_oovv_no1_x0


! -----------------------------------------------------------------------
