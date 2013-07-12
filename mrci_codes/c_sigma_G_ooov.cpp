

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


double orz::ct::sigma_G_ooov(const orz::ct::Input &ctinp,
			     const orz::ct::SymBlockInfo &symblockinfo,
			     const orz::ct::HintMO &hintmo,
			     const orz::ct::RdmPack &rdmPack_sym,
			     const double init_value,
			     const orz::ct::BareAmpPack &T2,
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

  orz::DTensor moint1 = hintmo.int1(); // Setitng up one-body integrals
  const orz::DTensor moint1_sym = orz::ct::sympack_int1(symblockinfo, moint1); // moint1=(IR-COV index)

  double S0 = init_value;
  orz::DTensor V2(nmo,nmo,nmo); // V(p,nmo,nmo,nmo) for an orbital index p
  orz::DTensor T2b;             // Container of T2_ija,[b] tensor

  double * const V2_ptr = V2.cptr();


  {
  // No.0
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
    for(int sd = 0;sd < nir;++sd){
    for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
      T2b = T2.get_amp2(id);
      FC_FUNC(g_if_sigma_g_ooov_no0_x0, G_IF_SIGMA_G_OOOV_NO0_X0)
        (sc, ic, sd, id, V2_sym.cptr(), T2b.cptr(), Xaaaa.cptr(), nir, nsym, psym);
    }
    }
  }
  }
  FC_FUNC(g_if_sigma_g_ooov_no1_x0, G_IF_SIGMA_G_OOOV_NO1_X0)
    (Xaaaa.cptr(), &S0, nir, nsym, psym);
  }


  {
  // No.1
  for(int sd = 0;sd < nir;++sd){
  for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
    T2b = T2.get_amp2(id);
    orz::DTensor X(nocc);
    orz::DTensor Xa = orz::ct::sympack_Xa(symblockinfo, sd, X);
    FC_FUNC(g_if_sigma_g_ooov_no0_x1, G_IF_SIGMA_G_OOOV_NO0_X1)
      (sd, id, T2b.cptr(), Xa.cptr(), nir, nsym, psym);
    FC_FUNC(g_if_sigma_g_ooov_no1_x1, G_IF_SIGMA_G_OOOV_NO1_X1)
      (sd, id, moint1_sym.cptr(), Xa.cptr(), &S0, nir, nsym, psym);
  }
  }
  }


  {
  // No.2
  for(int sd = 0;sd < nir;++sd){
  for(int id = symblockinfo.psym()(sd,I_V,I_BEGIN);id <= symblockinfo.psym()(sd,I_V,I_END);++id){
    // Load ERIs from somewhere, e.g. disk, GA, etc..
    V2 <<= 0.0;
    shared_ptr<orz::hint::LoadSortBuffer> loadbuf_ptr = hintmo.get_LoadSortBuffer(id);
    for(; !loadbuf_ptr->out_of_range(); ++(*loadbuf_ptr)) {
      // Load a signle record of integals
      const int &imo2 = loadbuf_ptr->i0;
      const int &imo3 = loadbuf_ptr->i1;
      const int &imo4 = loadbuf_ptr->i2;
      const double &v = loadbuf_ptr->v;
      V2_ptr[((imo2)*nmo+imo3)*nmo+imo4] = v;
    }
    const orz::DTensor V2_sym = orz::ct::sympack_int2(symblockinfo, id, sd, V2); // V2=(IR-COV index)
    orz::DTensor X(nocc, nocc, nocc);
    orz::DTensor Xaaa = orz::ct::sympack_Xaaa(symblockinfo, sd, X);
    FC_FUNC(g_if_sigma_g_ooov_no0_x2, G_IF_SIGMA_G_OOOV_NO0_X2)
      (sd, id, V2_sym.cptr(), Xaaa.cptr(), nir, nsym, psym);
    T2b = T2.get_amp2(id);
    FC_FUNC(g_if_sigma_g_ooov_no1_x2, G_IF_SIGMA_G_OOOV_NO1_X2)
      (sd, id, T2b.cptr(), Xaaa.cptr(), &S0, nir, nsym, psym);
  }
  }
  }


  return S0;
}

