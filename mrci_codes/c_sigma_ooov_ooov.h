
#ifndef C_SIGMA_OOOV_OOOV_H
#define C_SIGMA_OOOV_OOOV_H


// #include <tensor/tensor.h>
// #include <sci/hint/hintmo/hintmo.h>
// #include <sci/ctnew2/ctclass_input.h>
// #include <sci/ctnew2/ctclass_symblock.h>
// #include <sci/ctnew2/ctclass_rdmpack.h>
// #include <sci/ctnew2/ctclass_bareamppack.h>
 
extern "C"{

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x0, G_IF_SIGMA_OOOV_OOOV_NO0_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x0, G_IF_SIGMA_OOOV_OOOV_NO1_X0)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x1, G_IF_SIGMA_OOOV_OOOV_NO0_X1)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x1, G_IF_SIGMA_OOOV_OOOV_NO1_X1)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x2, G_IF_SIGMA_OOOV_OOOV_NO0_X2)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x2, G_IF_SIGMA_OOOV_OOOV_NO1_X2)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x3, G_IF_SIGMA_OOOV_OOOV_NO0_X3)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x3, G_IF_SIGMA_OOOV_OOOV_NO1_X3)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x4, G_IF_SIGMA_OOOV_OOOV_NO0_X4)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x4, G_IF_SIGMA_OOOV_OOOV_NO1_X4)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x5, G_IF_SIGMA_OOOV_OOOV_NO0_X5)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x5, G_IF_SIGMA_OOOV_OOOV_NO1_X5)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x6, G_IF_SIGMA_OOOV_OOOV_NO0_X6)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x6, G_IF_SIGMA_OOOV_OOOV_NO1_X6)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x7, G_IF_SIGMA_OOOV_OOOV_NO0_X7)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x7, G_IF_SIGMA_OOOV_OOOV_NO1_X7)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x8, G_IF_SIGMA_OOOV_OOOV_NO0_X8)
  (const FC_INT &sb, const FC_INT &ib, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x8, G_IF_SIGMA_OOOV_OOOV_NO1_X8)
  (const FC_INT &sb, const FC_INT &ib, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x9, G_IF_SIGMA_OOOV_OOOV_NO0_X9)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x9, G_IF_SIGMA_OOOV_OOOV_NO1_X9)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x10, G_IF_SIGMA_OOOV_OOOV_NO0_X10)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x10, G_IF_SIGMA_OOOV_OOOV_NO1_X10)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x11, G_IF_SIGMA_OOOV_OOOV_NO0_X11)
  (const FC_INT &sAk, const FC_INT &iAk, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x11, G_IF_SIGMA_OOOV_OOOV_NO1_X11)
  (const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x12, G_IF_SIGMA_OOOV_OOOV_NO0_X12)
  (const FC_INT &sAk, const FC_INT &iAk, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x12, G_IF_SIGMA_OOOV_OOOV_NO1_X12)
  (const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x13, G_IF_SIGMA_OOOV_OOOV_NO0_X13)
  (const FC_INT &sAk, const FC_INT &iAk, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x13, G_IF_SIGMA_OOOV_OOOV_NO1_X13)
  (const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x14, G_IF_SIGMA_OOOV_OOOV_NO0_X14)
  (const FC_INT &sAk, const FC_INT &iAk, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x14, G_IF_SIGMA_OOOV_OOOV_NO1_X14)
  (const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x15, G_IF_SIGMA_OOOV_OOOV_NO0_X15)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x15, G_IF_SIGMA_OOOV_OOOV_NO1_X15)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x16, G_IF_SIGMA_OOOV_OOOV_NO0_X16)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x16, G_IF_SIGMA_OOOV_OOOV_NO1_X16)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x17, G_IF_SIGMA_OOOV_OOOV_NO0_X17)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x17, G_IF_SIGMA_OOOV_OOOV_NO1_X17)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x18, G_IF_SIGMA_OOOV_OOOV_NO0_X18)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x18, G_IF_SIGMA_OOOV_OOOV_NO1_X18)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x19, G_IF_SIGMA_OOOV_OOOV_NO0_X19)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x19, G_IF_SIGMA_OOOV_OOOV_NO1_X19)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x20, G_IF_SIGMA_OOOV_OOOV_NO0_X20)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x20, G_IF_SIGMA_OOOV_OOOV_NO1_X20)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x21, G_IF_SIGMA_OOOV_OOOV_NO0_X21)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x21, G_IF_SIGMA_OOOV_OOOV_NO1_X21)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x22, G_IF_SIGMA_OOOV_OOOV_NO0_X22)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x22, G_IF_SIGMA_OOOV_OOOV_NO1_X22)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x23, G_IF_SIGMA_OOOV_OOOV_NO0_X23)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x23, G_IF_SIGMA_OOOV_OOOV_NO1_X23)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x24, G_IF_SIGMA_OOOV_OOOV_NO0_X24)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &se, const FC_INT &ie, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x24, G_IF_SIGMA_OOOV_OOOV_NO1_X24)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x25, G_IF_SIGMA_OOOV_OOOV_NO0_X25)
  (const FC_INT &sAj, const FC_INT &iAj, const FC_INT &se, const FC_INT &ie, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x25, G_IF_SIGMA_OOOV_OOOV_NO1_X25)
  (const FC_INT &sAj, const FC_INT &iAj, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x26, G_IF_SIGMA_OOOV_OOOV_NO0_X26)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &se, const FC_INT &ie, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x26, G_IF_SIGMA_OOOV_OOOV_NO1_X26)
  (const FC_INT &sa, const FC_INT &ia, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x27, G_IF_SIGMA_OOOV_OOOV_NO0_X27)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x27, G_IF_SIGMA_OOOV_OOOV_NO1_X27)
  (const FC_INT &sAi, const FC_INT &iAi, const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x28, G_IF_SIGMA_OOOV_OOOV_NO0_X28)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no1_x28, G_IF_SIGMA_OOOV_OOOV_NO1_X28)
  (const FC_INT &sAi, const FC_INT &iAi, const FC_INT &sAk, const FC_INT &iAk, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x29, G_IF_SIGMA_OOOV_OOOV_NO0_X29)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const Ecas, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_ooov_ooov_no0_x30, G_IF_SIGMA_OOOV_OOOV_NO0_X30)
  (const FC_INT &sVa, const FC_INT &iVa, 
   const double * const Ecas, 
   const double * const T2, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);



} 


#endif

