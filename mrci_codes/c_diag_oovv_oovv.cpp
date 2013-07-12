

#include <orz/orz.h>
#include <orz/openmp.h>
#include <orz/cblas.h>
#include <orz/clapack.h>
#include <tensor/tensor.h>
#include <sci/hint/para_disttools.h>
#include <sci/ctnew2/ct.h>
#include <sci/ctnew2/ct_f.h>
#include <sci/ctnew2/ctclass_input.h>
#include <sci/ctnew2/ctclass_symblock.h>
#include <sci/ctnew2/ctclass_hintmo.h>
#include <sci/ctnew2/ctclass_rdmpack.h>
#include <sci/ctnew2/ctclass_bareamppack.h>
#include <sci/ctnew2/ctclass_orthamppack.h>
#include <sci/ctnew2/diaghessian.h>
#include <sci/ctnew2/symamp2.h>
#include <sci/ctnew2/mrci.h>

using std::cout;
using std::endl;

#define FLOPCOUNT


// ***************************************************************************
// orz::ct::mrci
// ***************************************************************************
									    /*!
   @brief CT input
									     */


orz::ct::BareAmpPack orz::ct::diag_oovv_oovv(const orz::ct::Input &ctinp,
					 const orz::ct::SymBlockInfo &symblockinfo,
					 const orz::ct::HintMO &hintmo,
					 const orz::ct::RdmPack &rdmPack_sym,
					 const orz::DTensor &rdm4,
					 const double zeroth_energy) {

  // set up nmo nclosed, nocc
  const FC_INT nclosed = ctinp.nclosed();
  const FC_INT nocc    = ctinp.nocc();
  const FC_INT nvir    = ctinp.nvir();
  const FC_INT nmo     = nclosed + nocc + nvir;

  const FC_INT nir = symblockinfo.nir();
  const FC_INT * const nsym = symblockinfo.nsym().cptr();
  const FC_INT * const psym = symblockinfo.psym().cptr();
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();

  orz::DTensor moint1 = hintmo.int1(); // Setting up the one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  orz::ct::BareAmpPack retval;
  retval = orz::ct::BareAmpPack(ctinp, symblockinfo, "Hdiag");
  orz::DTensor V(nmo,nmo,nmo);
  orz::DTensor rdm4_sym;
  double * const V_ptr = V.cptr();



  {
  // No.0
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V); // V2=(IR-COV index)
    FC_FUNC(g_if_diag_oovv_oovv_no0_x0, G_IF_DIAG_OOVV_OOVV_NO0_X0)
      (sVb, iVb, V_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.1
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V); // V2=(IR-COV index)
    FC_FUNC(g_if_diag_oovv_oovv_no0_x1, G_IF_DIAG_OOVV_OOVV_NO0_X1)
      (sVb, iVb, V_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.2
  orz::DTensor X(nocc, nocc, nvir);
  orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, 0, X);
  FC_FUNC(g_if_diag_oovv_oovv_no0_x2, G_IF_DIAG_OOVV_OOVV_NO0_X2)
    (moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x2, G_IF_DIAG_OOVV_OOVV_NO1_X2)
      (sVb, iVb, Xaav.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.3
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    orz::DTensor X(nocc, nocc);
    orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVb, X);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x3, G_IF_DIAG_OOVV_OOVV_NO0_X3)
      (sVb, iVb, moint1_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x3, G_IF_DIAG_OOVV_OOVV_NO1_X3)
      (sVb, iVb, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.4
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x4, G_IF_DIAG_OOVV_OOVV_NO0_X4)
      (sVa, iVa, moint1_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  {
  // No.5
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc);
    orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x5, G_IF_DIAG_OOVV_OOVV_NO0_X5)
      (sVa, iVa, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
      FC_FUNC(g_if_diag_oovv_oovv_no1_x5, G_IF_DIAG_OOVV_OOVV_NO1_X5)
        (sVb, iVb, sVa, iVa, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.6
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc);
    orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x6, G_IF_DIAG_OOVV_OOVV_NO0_X6)
      (sVa, iVa, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
      FC_FUNC(g_if_diag_oovv_oovv_no1_x6, G_IF_DIAG_OOVV_OOVV_NO1_X6)
        (sVb, iVb, sVa, iVa, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, Hdiagb);
    }
    }
  }
  }
  }


  {
  // No.7
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc);
    orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVb, X);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x7, G_IF_DIAG_OOVV_OOVV_NO0_X7)
      (sVb, iVb, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x7, G_IF_DIAG_OOVV_OOVV_NO1_X7)
      (sVb, iVb, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.8
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc);
    orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVb, X);
    FC_FUNC(g_if_diag_oovv_oovv_no0_x8, G_IF_DIAG_OOVV_OOVV_NO0_X8)
      (sVb, iVb, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x8, G_IF_DIAG_OOVV_OOVV_NO1_X8)
      (sVb, iVb, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.9
  orz::DTensor X(nocc, nocc);
  orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, 0, X);
  FC_FUNC(g_if_diag_oovv_oovv_no0_x9, G_IF_DIAG_OOVV_OOVV_NO0_X9)
    (moint1_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x9, G_IF_DIAG_OOVV_OOVV_NO1_X9)
      (sVb, iVb, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.10
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V); // V2=(IR-COV index)
    FC_FUNC(g_if_diag_oovv_oovv_no0_x10, G_IF_DIAG_OOVV_OOVV_NO0_X10)
      (sVa, iVa, V_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  {
  // No.11
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V); // V2=(IR-COV index)
    FC_FUNC(g_if_diag_oovv_oovv_no0_x11, G_IF_DIAG_OOVV_OOVV_NO0_X11)
      (sVa, iVa, V_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  {
  // No.12
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V); // V2=(IR-COV index)
    FC_FUNC(g_if_diag_oovv_oovv_no0_x12, G_IF_DIAG_OOVV_OOVV_NO0_X12)
      (sVa, iVa, V_sym.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  {
  // No.13
  orz::DTensor X(nocc, nocc);
  orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, 0, X);
  FC_FUNC(g_if_diag_oovv_oovv_no0_x13, G_IF_DIAG_OOVV_OOVV_NO0_X13)
    (moint1_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x13, G_IF_DIAG_OOVV_OOVV_NO1_X13)
      (sVa, iVa, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  {
  // No.14
  orz::DTensor X(nocc, nocc);
  orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, 0, X);
  for(int sa = 0;sa < nir;++sa){
  for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V); // V2=(IR-COV index)
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_O,I_BEGIN);ic <= symblockinfo.psym()(sc,I_O,I_END);++ic){
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ia];
      int imoj = amo2imo[ic];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ia, sa, ic, sc, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sa, sc, ia, ic, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_diag_oovv_oovv_no0_x14, G_IF_DIAG_OOVV_OOVV_NO0_X14)
        (sa, ia, sc, ic, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVb]);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x14, G_IF_DIAG_OOVV_OOVV_NO1_X14)
      (sVb, iVb, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, Hdiagb);
  }
  }
  }


  {
  // No.15
  orz::DTensor X(nocc, nocc);
  orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, 0, X);
  for(int sa = 0;sa < nir;++sa){
  for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V); // V2=(IR-COV index)
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_O,I_BEGIN);ic <= symblockinfo.psym()(sc,I_O,I_END);++ic){
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ia];
      int imoj = amo2imo[ic];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ia, sa, ic, sc, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sa, sc, ia, ic, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_diag_oovv_oovv_no0_x15, G_IF_DIAG_OOVV_OOVV_NO0_X15)
        (sa, ia, sc, ic, V_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    Hdiagb = orz::DTensor(retval.namps_iamp()[iVa]);
    FC_FUNC(g_if_diag_oovv_oovv_no1_x15, G_IF_DIAG_OOVV_OOVV_NO1_X15)
      (sVa, iVa, Xaa.cptr(), Hdiagb.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, Hdiagb);
  }
  }
  }


  return retval;
}

