

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


orz::ct::BareAmpPack orz::ct::sigma_ooov_ooov(const orz::ct::Input &ctinp,
					 const orz::ct::SymBlockInfo &symblockinfo,
					 const orz::ct::HintMO &hintmo,
					 const orz::ct::RdmPack &rdmPack_sym,
					 const orz::DTensor &rdm4,
					 const orz::ct::BareAmpPack &T2,
					 const int num_sigma,
					 const double T0,
					 const double zeroth_energy) {


  // set up nmo nclosed, nocc
  const FC_INT nclosed = ctinp.nclosed();
  const FC_INT nocc    = ctinp.nocc();
  const FC_INT nvir    = ctinp.nvir();
  const FC_INT nmo     = nclosed + nocc + nvir;
  const FC_INT nir     = symblockinfo.nir();
  const FC_INT * const nsym    = symblockinfo.nsym().cptr();
  const FC_INT * const psym    = symblockinfo.psym().cptr();
  const FC_INT * const amo2imo = symblockinfo.amo2imo().cptr();

  orz::DTensor moint1 = hintmo.int1(); // Setting up one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  std::ostringstream stm;
  stm << num_sigma;
  std::string name_of_sigma = "S2[" + stm.str() + "]"; // Name of the Sigma vector
  orz::ct::BareAmpPack retval
    = orz::ct::BareAmpPack(ctinp, symblockinfo, name_of_sigma); // Sigma(a, a', e, e') tensor

  orz::DTensor S2b; // Container of S2_aae,[b] tensor
  orz::DTensor T2b; // Container of T2_aae,[b] tensor

  orz::DTensor V2(nmo,nmo,nmo);
  orz::DTensor rdm4_sym;
  double * const V2_ptr = V2.cptr();



  {
  // No.0
  for(int sa = 0;sa < nir;++sa){
  for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x0, G_IF_SIGMA_OOOV_OOOV_NO0_X0)
      (sa, ia, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x0, G_IF_SIGMA_OOOV_OOOV_NO1_X0)
        (sa, ia, sVa, iVa, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.1
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x1, G_IF_SIGMA_OOOV_OOOV_NO0_X1)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x1, G_IF_SIGMA_OOOV_OOOV_NO1_X1)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.2
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x2, G_IF_SIGMA_OOOV_OOOV_NO0_X2)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x2, G_IF_SIGMA_OOOV_OOOV_NO1_X2)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.3
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_ooov_ooov_no0_x3, G_IF_SIGMA_OOOV_OOOV_NO0_X3)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x3, G_IF_SIGMA_OOOV_OOOV_NO1_X3)
      (sVa, iVa, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.4
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_ooov_ooov_no0_x4, G_IF_SIGMA_OOOV_OOOV_NO0_X4)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x4, G_IF_SIGMA_OOOV_OOOV_NO1_X4)
      (sVa, iVa, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.5
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x5, G_IF_SIGMA_OOOV_OOOV_NO0_X5)
        (sc, ic, sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x5, G_IF_SIGMA_OOOV_OOOV_NO1_X5)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.6
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x6, G_IF_SIGMA_OOOV_OOOV_NO0_X6)
      (sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x6, G_IF_SIGMA_OOOV_OOOV_NO1_X6)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.7
  for(int sa = 0;sa < nir;++sa){
  for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x7, G_IF_SIGMA_OOOV_OOOV_NO0_X7)
      (sa, ia, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x7, G_IF_SIGMA_OOOV_OOOV_NO1_X7)
        (sa, ia, sVa, iVa, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.8
  for(int sb = 0;sb < nir;++sb){
  for(int ib = symblockinfo.psym()(sb,I_O,I_BEGIN);ib <= symblockinfo.psym()(sb,I_O,I_END);++ib){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ib);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ib, sb, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sb, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x8, G_IF_SIGMA_OOOV_OOOV_NO0_X8)
      (sb, ib, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x8, G_IF_SIGMA_OOOV_OOOV_NO1_X8)
        (sb, ib, sVa, iVa, T2b.cptr(), Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.9
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x9, G_IF_SIGMA_OOOV_OOOV_NO0_X9)
        (sc, ic, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x9, G_IF_SIGMA_OOOV_OOOV_NO1_X9)
      (sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.10
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x10, G_IF_SIGMA_OOOV_OOOV_NO0_X10)
        (sc, ic, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x10, G_IF_SIGMA_OOOV_OOOV_NO1_X10)
      (sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.11
  for(int sAk = 0;sAk < nir;++sAk){
  for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iAk);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iAk, sAk, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sAk, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x11, G_IF_SIGMA_OOOV_OOOV_NO0_X11)
      (sAk, iAk, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x11, G_IF_SIGMA_OOOV_OOOV_NO1_X11)
        (sAk, iAk, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.12
  for(int sAk = 0;sAk < nir;++sAk){
  for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iAk);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iAk, sAk, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sAk, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x12, G_IF_SIGMA_OOOV_OOOV_NO0_X12)
      (sAk, iAk, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x12, G_IF_SIGMA_OOOV_OOOV_NO1_X12)
        (sAk, iAk, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.13
  for(int sAk = 0;sAk < nir;++sAk){
  for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iAk);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iAk, sAk, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sAk, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x13, G_IF_SIGMA_OOOV_OOOV_NO0_X13)
      (sAk, iAk, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x13, G_IF_SIGMA_OOOV_OOOV_NO1_X13)
        (sAk, iAk, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.14
  for(int sAk = 0;sAk < nir;++sAk){
  for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iAk);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iAk, sAk, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sAk, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x14, G_IF_SIGMA_OOOV_OOOV_NO0_X14)
      (sAk, iAk, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x14, G_IF_SIGMA_OOOV_OOOV_NO1_X14)
        (sAk, iAk, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.15
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sa = 0;sa < nir;++sa){
    for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
      // Load ERIs from somewhere, e.g. disk, GA, etc..
      V2 <<= 0.0;
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ia);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
        // Load a signle record of integals
        const int &imo2 = loadbuf_ptr->i0;
        const int &imo3 = loadbuf_ptr->i1;
        const int &imo4 = loadbuf_ptr->i2;
        const double &v = loadbuf_ptr->v;
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
      }
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ia, sa, V2); // V2=(IR-COV index)
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x15, G_IF_SIGMA_OOOV_OOOV_NO0_X15)
        (sa, ia, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x15, G_IF_SIGMA_OOOV_OOOV_NO1_X15)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.16
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x16, G_IF_SIGMA_OOOV_OOOV_NO0_X16)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x16, G_IF_SIGMA_OOOV_OOOV_NO1_X16)
      (sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.17
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x17, G_IF_SIGMA_OOOV_OOOV_NO0_X17)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x17, G_IF_SIGMA_OOOV_OOOV_NO1_X17)
      (sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.18
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x18, G_IF_SIGMA_OOOV_OOOV_NO0_X18)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x18, G_IF_SIGMA_OOOV_OOOV_NO1_X18)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.19
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x19, G_IF_SIGMA_OOOV_OOOV_NO0_X19)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x19, G_IF_SIGMA_OOOV_OOOV_NO1_X19)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.20
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x20, G_IF_SIGMA_OOOV_OOOV_NO0_X20)
      (sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x20, G_IF_SIGMA_OOOV_OOOV_NO1_X20)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.21
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x21, G_IF_SIGMA_OOOV_OOOV_NO0_X21)
      (sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x21, G_IF_SIGMA_OOOV_OOOV_NO1_X21)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.22
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x22, G_IF_SIGMA_OOOV_OOOV_NO0_X22)
      (sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x22, G_IF_SIGMA_OOOV_OOOV_NO1_X22)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.23
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x23, G_IF_SIGMA_OOOV_OOOV_NO0_X23)
        (sd, id, sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_ooov_ooov_no1_x23, G_IF_SIGMA_OOOV_OOOV_NO1_X23)
      (sVa, iVa, Xaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.24
  for(int sc = 0;sc < nir;++sc){
  for(int ic = symblockinfo.psym()(sc,I_O,I_BEGIN);ic <= symblockinfo.psym()(sc,I_O,I_END);++ic){
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sc, X);
    for(int se = 0;se < nir;++se){
    for(int ie = symblockinfo.psym()(se,I_O,I_BEGIN);ie <= symblockinfo.psym()(se,I_O,I_END);++ie){
      // Load ERIs from somewhere, e.g. disk, GA, etc..
      V2 <<= 0.0;
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ie);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
        // Load a signle record of integals
        const int &imo2 = loadbuf_ptr->i0;
        const int &imo3 = loadbuf_ptr->i1;
        const int &imo4 = loadbuf_ptr->i2;
        const double &v = loadbuf_ptr->v;
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
      }
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ie, se, V2); // V2=(IR-COV index)
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ic];
      int imoj = amo2imo[ie];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ic, sc, ie, se, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sc, se, ic, ie, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x24, G_IF_SIGMA_OOOV_OOOV_NO0_X24)
        (sc, ic, se, ie, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x24, G_IF_SIGMA_OOOV_OOOV_NO1_X24)
        (sc, ic, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.25
  for(int sAj = 0;sAj < nir;++sAj){
  for(int iAj = symblockinfo.psym()(sAj,I_O,I_BEGIN);iAj <= symblockinfo.psym()(sAj,I_O,I_END);++iAj){
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sAj, X);
    for(int se = 0;se < nir;++se){
    for(int ie = symblockinfo.psym()(se,I_O,I_BEGIN);ie <= symblockinfo.psym()(se,I_O,I_END);++ie){
      // Load ERIs from somewhere, e.g. disk, GA, etc..
      V2 <<= 0.0;
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ie);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
        // Load a signle record of integals
        const int &imo2 = loadbuf_ptr->i0;
        const int &imo3 = loadbuf_ptr->i1;
        const int &imo4 = loadbuf_ptr->i2;
        const double &v = loadbuf_ptr->v;
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
      }
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ie, se, V2); // V2=(IR-COV index)
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[iAj];
      int imoj = amo2imo[ie];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, iAj, sAj, ie, se, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sAj, se, iAj, ie, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x25, G_IF_SIGMA_OOOV_OOOV_NO0_X25)
        (sAj, iAj, se, ie, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x25, G_IF_SIGMA_OOOV_OOOV_NO1_X25)
        (sAj, iAj, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.26
  for(int sa = 0;sa < nir;++sa){
  for(int ia = symblockinfo.psym()(sa,I_O,I_BEGIN);ia <= symblockinfo.psym()(sa,I_O,I_END);++ia){
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sa, X);
    for(int se = 0;se < nir;++se){
    for(int ie = symblockinfo.psym()(se,I_O,I_BEGIN);ie <= symblockinfo.psym()(se,I_O,I_END);++ie){
      // Load ERIs from somewhere, e.g. disk, GA, etc..
      V2 <<= 0.0;
      shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ie);
      for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
        // Load a signle record of integals
        const int &imo2 = loadbuf_ptr->i0;
        const int &imo3 = loadbuf_ptr->i1;
        const int &imo4 = loadbuf_ptr->i2;
        const double &v = loadbuf_ptr->v;
        V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
      }
      const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ie, se, V2); // V2=(IR-COV index)
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ie];
      int imoj = amo2imo[ia];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ie, se, ia, sa, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(se, sa, ie, ia, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x26, G_IF_SIGMA_OOOV_OOOV_NO0_X26)
        (sa, ia, se, ie, V2_sym.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      S2b = orz::DTensor(retval.namps_iamp()[iVa]);
      T2b = T2.get_amp2(iVa);
      FC_FUNC(g_if_sigma_ooov_ooov_no1_x26, G_IF_SIGMA_OOOV_OOOV_NO1_X26)
        (sa, ia, sVa, iVa, T2b.cptr(), Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVa, S2b);
    }
    }
  }
  }
  }


  {
  // No.27
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x27, G_IF_SIGMA_OOOV_OOOV_NO0_X27)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    for(int sAi = 0;sAi < nir;++sAi){
    for(int iAi = symblockinfo.psym()(sAi,I_O,I_BEGIN);iAi <= symblockinfo.psym()(sAi,I_O,I_END);++iAi){
      for(int sAk = 0;sAk < nir;++sAk){
      for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
        // Load D4 from disk, or GA .... 
        int imoi = amo2imo[iAi];
        int imoj = amo2imo[iAk];

        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                           orz::Slice(),            orz::Slice(),
                                           orz::Slice(),            orz::Slice(),
                                           orz::Slice(),            orz::Slice()).copy();
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, iAi, sAi, iAk, sAk, rdm4_ij_sliced);
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sAi, sAk, iAi, iAk, rdm4_sym.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_sigma_ooov_ooov_no1_x27, G_IF_SIGMA_OOOV_OOOV_NO1_X27)
          (sAi, iAi, sAk, iAk, sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.28
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVa);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVa, sVa, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nocc);
    orz::DTensor Xaaaaa = orz::ct::sympack_Xaaaaa(symblockinfo, sVa, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_ooov_ooov_no0_x28, G_IF_SIGMA_OOOV_OOOV_NO0_X28)
        (sd, id, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaaa.cptr(), nir, nsym, psym);
    }
    }
    for(int sAi = 0;sAi < nir;++sAi){
    for(int iAi = symblockinfo.psym()(sAi,I_O,I_BEGIN);iAi <= symblockinfo.psym()(sAi,I_O,I_END);++iAi){
      for(int sAk = 0;sAk < nir;++sAk){
      for(int iAk = symblockinfo.psym()(sAk,I_O,I_BEGIN);iAk <= symblockinfo.psym()(sAk,I_O,I_END);++iAk){
        // Load D4 from disk, or GA .... 
        int imoi = amo2imo[iAi];
        int imoj = amo2imo[iAk];

        orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                           orz::Slice(),            orz::Slice(),
                                           orz::Slice(),            orz::Slice(),
                                           orz::Slice(),            orz::Slice()).copy();
        rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, iAi, sAi, iAk, sAk, rdm4_ij_sliced);
        FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sAi, sAk, iAi, iAk, rdm4_sym.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_sigma_ooov_ooov_no1_x28, G_IF_SIGMA_OOOV_OOOV_NO1_X28)
          (sAi, iAi, sAk, iAk, sVa, iVa, Xaaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
        FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
      }
      }
    }
    }
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.29
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x29, G_IF_SIGMA_OOOV_OOOV_NO0_X29)
      (sVa, iVa, &Ecas, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  {
  // No.30
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
    S2b = orz::DTensor(retval.namps_iamp()[iVa]);
    T2b = T2.get_amp2(iVa);
    FC_FUNC(g_if_sigma_ooov_ooov_no0_x30, G_IF_SIGMA_OOOV_OOOV_NO0_X30)
      (sVa, iVa, &Ecas, T2b.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVa, S2b);
  }
  }
  }


  return retval;
}

