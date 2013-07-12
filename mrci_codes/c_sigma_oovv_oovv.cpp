

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


orz::ct::BareAmpPack orz::ct::sigma_oovv_oovv(const orz::ct::Input &ctinp,
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
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sVb, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x0, G_IF_SIGMA_OOVV_OOVV_NO0_X0)
        (sd, id, sVb, iVb, V2_sym.cptr(), T2b.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x0, G_IF_SIGMA_OOVV_OOVV_NO1_X0)
      (sVb, iVb, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.1
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sVb, X);
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x1, G_IF_SIGMA_OOVV_OOVV_NO0_X1)
        (sd, id, sVb, iVb, V2_sym.cptr(), T2b.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x1, G_IF_SIGMA_OOVV_OOVV_NO1_X1)
      (sVb, iVb, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.2
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sVb, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x2, G_IF_SIGMA_OOVV_OOVV_NO0_X2)
        (sc, ic, sVb, iVb, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x2, G_IF_SIGMA_OOVV_OOVV_NO1_X2)
      (sVb, iVb, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.3
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sVb, X);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x3, G_IF_SIGMA_OOVV_OOVV_NO0_X3)
      (sVb, iVb, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x3, G_IF_SIGMA_OOVV_OOVV_NO1_X3)
      (sVb, iVb, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.4
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      T2b = T2.get_amp2(iVa);
      orz::DTensor X(nocc, nocc);
      orz::DTensor Xaa = orz::ct::sympack_Xaa(symblockinfo, sVb^sVa, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x4, G_IF_SIGMA_OOVV_OOVV_NO0_X4)
        (sVb, iVb, sVa, iVa, T2b.cptr(), moint1_sym.cptr(), Xaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x4, G_IF_SIGMA_OOVV_OOVV_NO1_X4)
        (sVb, iVb, sVa, iVa, Xaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.5
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    orz::DTensor X(nocc, nocc, nvir);
    orz::DTensor Xaav = orz::ct::sympack_Xaav(symblockinfo, sVb, X);
    FC_FUNC(g_if_sigma_oovv_oovv_no0_x5, G_IF_SIGMA_OOVV_OOVV_NO0_X5)
      (sVb, iVb, T2b.cptr(), moint1_sym.cptr(), Xaav.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x5, G_IF_SIGMA_OOVV_OOVV_NO1_X5)
      (sVb, iVb, Xaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.6
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sVb, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x6, G_IF_SIGMA_OOVV_OOVV_NO0_X6)
        (sc, ic, sVb, iVb, V2_sym.cptr(), T2b.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x6, G_IF_SIGMA_OOVV_OOVV_NO1_X6)
      (sVb, iVb, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.7
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc, nocc, nvir);
    orz::DTensor Xaaaav = orz::ct::sympack_Xaaaav(symblockinfo, sVb, X);
    for(int sc = 0;sc < nir;++sc){
    for(int ic = symblockinfo.psym()(sc,I_V,I_BEGIN);ic <= symblockinfo.psym()(sc,I_V,I_END);++ic){
      T2b = T2.get_amp2(ic);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x7, G_IF_SIGMA_OOVV_OOVV_NO0_X7)
        (sc, ic, sVb, iVb, V2_sym.cptr(), T2b.cptr(), Xaaaav.cptr(), nir, nsym, psym);
    }
    }
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x7, G_IF_SIGMA_OOVV_OOVV_NO1_X7)
      (sVb, iVb, Xaaaav.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.8
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
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
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      S2b = orz::DTensor(retval.namps_iamp()[iVb]);
      T2b = T2.get_amp2(iVb);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x8, G_IF_SIGMA_OOVV_OOVV_NO0_X8)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x8, G_IF_SIGMA_OOVV_OOVV_NO1_X8)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, S2b);
    }
    }
  }
  }
  }


  {
  // No.9
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
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
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      S2b = orz::DTensor(retval.namps_iamp()[iVb]);
      T2b = T2.get_amp2(iVb);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x9, G_IF_SIGMA_OOVV_OOVV_NO0_X9)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x9, G_IF_SIGMA_OOVV_OOVV_NO1_X9)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, S2b);
    }
    }
  }
  }
  }


  {
  // No.10
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      T2b = T2.get_amp2(iVa);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x10, G_IF_SIGMA_OOVV_OOVV_NO0_X10)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x10, G_IF_SIGMA_OOVV_OOVV_NO1_X10)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.11
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(iVb);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, iVb, sVb, V2); // V2=(IR-COV index)
    for(int sVa = 0;sVa < nir;++sVa){
    for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
      T2b = T2.get_amp2(iVa);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x11, G_IF_SIGMA_OOVV_OOVV_NO0_X11)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x11, G_IF_SIGMA_OOVV_OOVV_NO1_X11)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    }
    }
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.12
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
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
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      S2b = orz::DTensor(retval.namps_iamp()[iVb]);
      T2b = T2.get_amp2(iVb);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x12, G_IF_SIGMA_OOVV_OOVV_NO0_X12)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x12, G_IF_SIGMA_OOVV_OOVV_NO1_X12)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, S2b);
    }
    }
  }
  }
  }


  {
  // No.13
  for(int sVa = 0;sVa < nir;++sVa){
  for(int iVa = symblockinfo.psym()(sVa,I_V,I_BEGIN);iVa <= symblockinfo.psym()(sVa,I_V,I_END);++iVa){
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
    for(int sVb = 0;sVb < nir;++sVb){
    for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
      S2b = orz::DTensor(retval.namps_iamp()[iVb]);
      T2b = T2.get_amp2(iVb);
      orz::DTensor X(nocc, nocc, nocc, nocc);
      orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, sVa^sVb, X);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x13, G_IF_SIGMA_OOVV_OOVV_NO0_X13)
        (sVb, iVb, sVa, iVa, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no1_x13, G_IF_SIGMA_OOVV_OOVV_NO1_X13)
        (sVb, iVb, sVa, iVa, Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
      retval.acc_amp2(iVb, S2b);
    }
    }
  }
  }
  }


  {
  // No.14
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x14, G_IF_SIGMA_OOVV_OOVV_NO0_X14)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x14, G_IF_SIGMA_OOVV_OOVV_NO1_X14)
      (sVb, iVb, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.15
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  FC_FUNC(g_if_sigma_oovv_oovv_no0_x15, G_IF_SIGMA_OOVV_OOVV_NO0_X15)
    (moint1_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x15, G_IF_SIGMA_OOVV_OOVV_NO1_X15)
      (sVb, iVb, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.16
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sc = 0;sc < nir;++sc){
  for(int ic = symblockinfo.psym()(sc,I_O,I_BEGIN);ic <= symblockinfo.psym()(sc,I_O,I_END);++ic){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index)
    for(int se = 0;se < nir;++se){
    for(int ie = symblockinfo.psym()(se,I_O,I_BEGIN);ie <= symblockinfo.psym()(se,I_O,I_END);++ie){
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ic];
      int imoj = amo2imo[ie];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ic, sc, ie, se, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sc, se, ic, ie, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x16, G_IF_SIGMA_OOVV_OOVV_NO0_X16)
        (sc, ic, se, ie, V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x16, G_IF_SIGMA_OOVV_OOVV_NO1_X16)
      (sVb, iVb, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  {
  // No.17
  orz::DTensor X(nocc, nocc, nocc, nocc);
  orz::DTensor Xaaaa = orz::ct::sympack_Xaaaa(symblockinfo, 0, X);
  for(int sc = 0;sc < nir;++sc){
  for(int ic = symblockinfo.psym()(sc,I_O,I_BEGIN);ic <= symblockinfo.psym()(sc,I_O,I_END);++ic){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(ic);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, ic, sc, V2); // V2=(IR-COV index)
    for(int se = 0;se < nir;++se){
    for(int ie = symblockinfo.psym()(se,I_O,I_BEGIN);ie <= symblockinfo.psym()(se,I_O,I_END);++ie){
      // Load D4 from disk, or GA .... 
      int imoi = amo2imo[ic];
      int imoj = amo2imo[ie];

      orz::DTensor rdm4_ij_sliced = rdm4(orz::Slice(imoi,imoi+1), orz::Slice(imoj,imoj+1),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice(),
                                         orz::Slice(),            orz::Slice()).copy();
      rdm4_sym = orz::ct::sympack_rdm4_2(symblockinfo, ic, sc, ie, se, rdm4_ij_sliced);
      FC_FUNC(g_if_set_d4,G_IF_SET_D4)(sc, se, ic, ie, rdm4_sym.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_sigma_oovv_oovv_no0_x17, G_IF_SIGMA_OOVV_OOVV_NO0_X17)
        (sc, ic, se, ie, V2_sym.cptr(), Xaaaa.cptr(), nir, nsym, psym);
      FC_FUNC(g_if_unset_d4,G_IF_UNSET_D4)();
    }
    }
  }
  }
  for(int sVb = 0;sVb < nir;++sVb){
  for(int iVb = symblockinfo.psym()(sVb,I_V,I_BEGIN);iVb <= symblockinfo.psym()(sVb,I_V,I_END);++iVb){
    S2b = orz::DTensor(retval.namps_iamp()[iVb]);
    T2b = T2.get_amp2(iVb);
    FC_FUNC(g_if_sigma_oovv_oovv_no1_x17, G_IF_SIGMA_OOVV_OOVV_NO1_X17)
      (sVb, iVb, T2b.cptr(), Xaaaa.cptr(), S2b.cptr(), nir, nsym, psym);
    retval.acc_amp2(iVb, S2b);
  }
  }
  }


  return retval;
}

