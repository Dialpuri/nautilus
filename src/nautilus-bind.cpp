#include <emscripten/bind.h>
#include <emscripten.h>

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper.h>

#include "nucleicacid_db.h"
#include "nautilus-tools.h"
#include "nautilus-ss-find.h"
#include "nautilus-target.h"
#include "nautilus-join.h"
#include "nautilus-sequence.h"
#include "nautilus-rebuild-bases.h"
#include "nautilus-tidy.h"
#include "nautilus-util.h"
#include "nautilus-mlfind.h"


using namespace emscripten; 

clipper::Xmap<float> load_mtz_file(const std::string& file_name, const std::string& fobs, const std::string& fcalc) { 
  
  int ncyc = 3;
  bool doanis = false;
  int nhit = 100;
  double res_in = 2.0;         // Resolution limit
  double srchst = 18.0;        // Search angle step
  int verbose = 0;

  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  using clipper::data32::Compute_fphi_from_fsigf_phifom;
  using clipper::data32::Compute_scale_u_aniso_fphi;
  clipper::Resolution resol;
  clipper::CCP4MTZfile mtzfile;
  mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // Get work reflection data
  clipper::HKL_info hkls;
  mtzfile.open_read( file_name );
  double res = clipper::Util::max( mtzfile.resolution().limit(), res_in );
  resol = clipper::Resolution( res );
  hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );
  clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
  clipper::HKL_data<clipper::data32::ABCD>    wrk_hl( hkls );
  clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw( hkls );
  clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );
  clipper::HKL_data<clipper::data32::Flag>    flag( hkls );
  mtzfile.import_hkl_data( wrk_f , fobs );
  mtzfile.import_hkl_data( fphi,  fcalc );
  mtzfile.close_read();

  // do anisotropy correction
  clipper::U_aniso_orth uaniso( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

    // scale obs data
    typedef clipper::SFscale_aniso<float> SFscale;
    SFscale sfscl( 3.0, SFscale::SHARPEN );
    sfscl( wrk_f );
    uaniso = sfscl.u_aniso_orth( SFscale::F );
    // scale map coeffs
    Compute_scale_u_aniso_fphi compute_aniso( 1.0, -uaniso );
    fphi.compute( fphi, compute_aniso );
    // output
    std::cout << std::endl << "Applying anisotropy correction:"
              << std::endl << uaniso.format() << std::endl << std::endl;
  

  // apply free flag
  clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;
  //wrk_f1.mask( flag != 0 );
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_f1[ih] = clipper::data32::F_sigF();  //ugly hack for broken SGI compilers

  // work map
  wrk_hl.compute( wrk_pw, clipper::data32::Compute_abcd_from_phifom() );
  wrk_pw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );
  clipper::Spacegroup cspg = hkls.spacegroup();
  clipper::Cell       cxtl = hkls.cell();
  clipper::Grid_sampling grid( cspg, cxtl, hkls.resolution() );
  clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
  xwrk.fft_from( fphi );

  return xwrk;
}

void write_file(const clipper::MiniMol& minimol, const std::string& name) { 
  clipper::MMDBfile mfile;
  mfile.export_minimol(minimol);
  mfile.write_file(name);
}

extern "C" void run(const std::string& file_name, const std::string& fobs, const std::string& fcalc) { 
  clipper::Xmap<float> xwrk = load_mtz_file(file_name, fobs, fcalc);

  std::cout << xwrk.cell().format() << std::endl;

  clipper::MiniMol mol_wrk( xwrk.spacegroup(), xwrk.cell() );
  NucleicAcidTargets natools;
  natools.add_pdb( "/nautilus_lib.pdb" );
  natools.init_stats( xwrk );

  mol_wrk = NucleicAcidTools::flag_chains( mol_wrk );
  mol_wrk = natools.find( xwrk, mol_wrk, 50, 50, 20 );

  write_file(mol_wrk, "find.pdb");

  mol_wrk = natools.grow( xwrk, mol_wrk, 25, 0.001 );
  write_file(mol_wrk, "grow.pdb");

  NucleicAcidJoin na_join;
  mol_wrk = na_join.join( mol_wrk );
  write_file(mol_wrk, "join.pdb");


  mol_wrk = natools.link( xwrk, mol_wrk );
  write_file(mol_wrk, "link.pdb");


  mol_wrk = natools.prune( mol_wrk );
  write_file(mol_wrk, "prune.pdb");


  mol_wrk = natools.rebuild_chain( xwrk, mol_wrk );
  write_file(mol_wrk, "rebuild_chain.pdb");


  // NucleicAcidSequence na_seqnc;
  // mol_wrk = na_seqnc.sequence( xwrk, mol_wrk, seq_wrk );
  // write_file(mol_wrk, "sequence.pdb");


  NucleicAcidRebuildBases na_bases;
  mol_wrk = na_bases.rebuild_bases( xwrk, mol_wrk );
  write_file(mol_wrk, "rebuild_bases.pdb");


}

EMSCRIPTEN_BINDINGS(nautilus_module) { 
  function("load_mtz", &run);
}