
#ifndef C_SIGMA_OOVV_OOVV_H
#define C_SIGMA_OOVV_OOVV_H


// #include <tensor/tensor.h>
// #include <sci/hint/hintmo/hintmo.h>
// #include <sci/ctnew2/ctclass_input.h>
// #include <sci/ctnew2/ctclass_symblock.h>
// #include <sci/ctnew2/ctclass_rdmpack.h>
// #include <sci/ctnew2/ctclass_bareamppack.h>
 
extern "C"{

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x0, G_IF_SIGMA_OOVV_OOVV_NO0_X0)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVb, const FC_INT &iVb, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x0, G_IF_SIGMA_OOVV_OOVV_NO1_X0)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x1, G_IF_SIGMA_OOVV_OOVV_NO0_X1)
  (const FC_INT &sd, const FC_INT &id, const FC_INT &sVb, const FC_INT &iVb, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x1, G_IF_SIGMA_OOVV_OOVV_NO1_X1)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x2, G_IF_SIGMA_OOVV_OOVV_NO0_X2)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x2, G_IF_SIGMA_OOVV_OOVV_NO1_X2)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x3, G_IF_SIGMA_OOVV_OOVV_NO0_X3)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x3, G_IF_SIGMA_OOVV_OOVV_NO1_X3)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x4, G_IF_SIGMA_OOVV_OOVV_NO0_X4)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x4, G_IF_SIGMA_OOVV_OOVV_NO1_X4)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x5, G_IF_SIGMA_OOVV_OOVV_NO0_X5)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x5, G_IF_SIGMA_OOVV_OOVV_NO1_X5)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x6, G_IF_SIGMA_OOVV_OOVV_NO0_X6)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVb, const FC_INT &iVb, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x6, G_IF_SIGMA_OOVV_OOVV_NO1_X6)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x7, G_IF_SIGMA_OOVV_OOVV_NO0_X7)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &sVb, const FC_INT &iVb, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x7, G_IF_SIGMA_OOVV_OOVV_NO1_X7)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x8, G_IF_SIGMA_OOVV_OOVV_NO0_X8)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x8, G_IF_SIGMA_OOVV_OOVV_NO1_X8)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x9, G_IF_SIGMA_OOVV_OOVV_NO0_X9)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x9, G_IF_SIGMA_OOVV_OOVV_NO1_X9)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x10, G_IF_SIGMA_OOVV_OOVV_NO0_X10)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x10, G_IF_SIGMA_OOVV_OOVV_NO1_X10)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x11, G_IF_SIGMA_OOVV_OOVV_NO0_X11)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x11, G_IF_SIGMA_OOVV_OOVV_NO1_X11)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x12, G_IF_SIGMA_OOVV_OOVV_NO0_X12)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x12, G_IF_SIGMA_OOVV_OOVV_NO1_X12)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x13, G_IF_SIGMA_OOVV_OOVV_NO0_X13)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const V2, const double * const T2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x13, G_IF_SIGMA_OOVV_OOVV_NO1_X13)
  (const FC_INT &sVb, const FC_INT &iVb, const FC_INT &sVa, const FC_INT &iVa, 
   const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x14, G_IF_SIGMA_OOVV_OOVV_NO0_X14)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x14, G_IF_SIGMA_OOVV_OOVV_NO1_X14)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x15, G_IF_SIGMA_OOVV_OOVV_NO0_X15)
  (const double * const h, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x15, G_IF_SIGMA_OOVV_OOVV_NO1_X15)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x16, G_IF_SIGMA_OOVV_OOVV_NO0_X16)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &se, const FC_INT &ie, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x16, G_IF_SIGMA_OOVV_OOVV_NO1_X16)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no0_x17, G_IF_SIGMA_OOVV_OOVV_NO0_X17)
  (const FC_INT &sc, const FC_INT &ic, const FC_INT &se, const FC_INT &ie, 
   const double * const V2, const double * const X, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_sigma_oovv_oovv_no1_x17, G_IF_SIGMA_OOVV_OOVV_NO1_X17)
  (const FC_INT &sVb, const FC_INT &iVb, 
   const double * const T2, const double * const X, const double * const S2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);



} 


#endif

