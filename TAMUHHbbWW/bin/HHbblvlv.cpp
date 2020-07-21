#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/TAMUHHbbWW/interface/HttSystematics_MSSMRun2.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"
#include "TF1.h"

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template <typename T>
void To1Bin(T* proc)
{
    std::unique_ptr<TH1> originalHist = proc->ClonedScaledShape();
    TH1F *hist = new TH1F("hist","hist",1,0,1);
    double err = 0;
    double rate =
        originalHist->IntegralAndError(0, originalHist->GetNbinsX() + 1, err);
    hist->SetDirectory(0);
    hist->SetBinContent(1, rate);
    hist->SetBinError(1, err);
    proc->set_shape(*hist, true);  // True means adjust the process rate to the
                                   // integral of the hist
}

bool BinIsControlRegion(ch::Object const* obj)
{
    return boost::regex_search(obj->bin(),boost::regex{"_cr$"});
}

// Useful to have the inverse sometimes too
bool BinIsNotControlRegion(ch::Object const* obj)
{
    return !BinIsControlRegion(obj);
}



int main(int argc, char** argv) {
  // First define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  //string SM125= "";
  string mass = "300";
  //string input_folder = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/HH_WWbb/GGToX0ToHHTo2B2L2Nu_NoMjjbinsMjjcut_nnout_MTandMT2_MJJ_nnstep0p04_nncut0p0_limits/";
  string input_folder = "GGToX0ToHHTo2B2L2Nu_NoMjjbinsMjjcut_nnout_MTandMT2_MJJ_nnstep0p04_nncut0p0_limits";
  string output_folder = "shapes/";
  string postfix="th1shapes";
  bool auto_rebin = false;
  bool manual_rebin = false;
  bool real_data = false;
  int control_region = 0;
  bool check_neg_bins = false;
  bool poisson_bbb = false;
  bool do_w_weighting = true;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
    ("mass,m", po::value<string>(&mass)->default_value(mass))
    ("input_folder", po::value<string>(&input_folder)->default_value("GGToX0ToHHTo2B2L2Nu_NoMjjbinsMjjcut_nnout_MTandMT2_MJJ_nnstep0p04_nncut0p0_limits"))
    ("postfix", po::value<string>(&postfix)->default_value("th1shapes"))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("real_data", po::value<bool>(&real_data)->default_value(false))
    ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(false))
    ("output_folder", po::value<string>(&output_folder)->default_value("HHbblvlv_output/"))
    ("control_region", po::value<int>(&control_region)->default_value(0))
    ("check_neg_bins", po::value<bool>(&check_neg_bins)->default_value(false))
    ("poisson_bbb", po::value<bool>(&poisson_bbb)->default_value(false))
    ("w_weighting", po::value<bool>(&do_w_weighting)->default_value(true));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  typedef vector<string> VString;
  //typedef vector<pair<int, string>> Categories;
  string input_dir;
  input_dir  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/TAMUHHbbWW/shapes/"+input_folder+"/";

  VString chns =
  //    {"tt"};
 //     {"mt"};
      {"MuMu","ElEl","MuEl"};


  map<string, VString> bkg_procs;
  bkg_procs["MuMu"] = {"TTbar","SingleTop","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV"}; // DY:data-driven for MuMu, ElEl
  bkg_procs["ElEl"] = {"TTbar","SingleTop","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV"};
  bkg_procs["MuEl"] = {"TTbar","SingleTop","Drell_Yan", "ttV","VV"}; //DY: MC


  //Example - could fill this map with hardcoded binning for different
  //categories if manual_rebin is turned on
  //map<string, vector<double> > binning;
  //binning["et_nobtag"] = {500, 700, 900, 3900};
  //binning["et_btag"] = {500,3900};
  //binning["mt_nobtag"] = {500,700,900,1300,1700,1900,3900};
  //binning["mt_btag"] = {500,1300,3900};
  //binning["tt_nobtag"] = {500,3900};
  //binning["tt_btag"] = {500,3900};
  //binning["em_nobtag"] = {500,3900};
  //binning["em_btag"] = {500,3900};

  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;


  // in future, it can extend to one-btag, two-btag, Run2016, Run2017, Run2018
   ch::Categories cats = {
       {1, "2Mbtag"}//only one cat: >= 2 medium btag
   };
  // Uncomment this next line to see a *lot* of debug information
  // cb.SetVerbosity(3);

  // Here we will just define two categories for an 8TeV analysis. Each entry in
  // the vector below specifies a bin name and corresponding bin_id.
  //

  //vector<string> masses = {"90","100","110","120","130","140","160","180", "200", "250", "300", "350", "400", "450", "500", "600", "700", "800", "900","1000","1200","1400","1500","1600","1800","2000","2300","2600","2900","3200"};

  //map<string, VString> signal_types = {
  //  {"ggH", {"ggh_htautau", "ggH_Htautau", "ggA_Atautau"}},
  //  {"bbH", {"bbh_htautau", "bbH_Htautau", "bbA_Atautau"}}
  //};
  //if(mass=="MH"){
  //  signal_types = {
  //    {"ggH", {"ggH"}},
  //    {"bbH", {"bbH"}}
  //  };
  //}
  //  vector<string> sig_procs = {"ggH","bbH"};
  for(auto chn : chns){
    cb.AddObservations({mass}, {"GGToX0ToHHTo2B2L2Nu"}, {"13TeV"}, {chn}, cats);

    // add background
    cb.AddProcesses({mass}, {"GGToX0ToHHTo2B2L2Nu"}, {"13TeV"}, {chn}, bkg_procs[chn], cats, false);

    // add signal
    cb.AddProcesses({mass}, {"GGToX0ToHHTo2B2L2Nu"}, {"13TeV"}, {chn}, {"Signal"}, cats, true);


    }

  //add systematics 
  using ch::syst::SystMap;
  using ch::syst::era;
  using ch::syst::bin_id;
  using ch::syst::process;
  using ch::syst::channel;


  cb.cp().signals().AddSyst(cb, "lumi_13TeV_2016", "lnN", SystMap<>::init(1.025));
  for(auto chn : chns){
    for (auto bkg : bkg_procs[chn]){
	if (bkg != "data_untagged"){
	    cb.cp().channel({chn}).process({bkg}).AddSyst(cb, "lumi_13TeV_2016", "lnN", SystMap<>::init(1.025));
	    cb.cp().channel({chn}).process({bkg}).AddSyst(cb, "CMS_eff_b_heavy", "shape", SystMap<>::init(1.00));
	    cb.cp().channel({chn}).process({bkg}).AddSyst(cb, "CMS_eff_b_light", "shape", SystMap<>::init(1.00));
	    cb.cp().channel({chn}).process({bkg}).AddSyst(cb, "CMS_pdf", "shape", SystMap<>::init(1.00));
	    cb.cp().channel({chn}).process({bkg}).AddSyst(cb, "CMS_pu", "shape", SystMap<>::init(1.00));
	}
    
    }
  }

  cb.cp().process({"SingleTop"}).AddSyst(cb, "singleTop_xsec", "lnN", SystMap<>::init(1.072));
  cb.cp().process({"SingleTop_untagged"}).AddSyst(cb, "singleTop_xsec", "lnN", SystMap<>::init(1.072));

  cb.cp().process({"TTbar"}).AddSyst(cb, "ttbar_xsec", "lnN", SystMap<>::init(1.053));
  cb.cp().process({"TTbar_untagged"}).AddSyst(cb, "ttbar_xsec", "lnN", SystMap<>::init(1.053));

  cb.cp().channel({"MuEl"}).process({"Drell_Yan"}).AddSyst(cb, "dy_mc_xsec", "lnN", SystMap<>::init(1.05));


  cb.cp().AddSyst(cb,      "dy_rwgt_norm_$channel", "lnN", SystMap<channel, process>::init
	  ({"MuMu"}, {"data_untagged"},      1.05)
	  ({"MuMu"}, {"TTbar_untagged"},     1.05)
	  ({"MuMu"}, {"SingleTop_untagged"}, 1.05)
	  ({"ElEl"}, {"data_untagged"},      1.05)
	  ({"ElEl"}, {"TTbar_untagged"},     1.05)
	  ({"ElEl"}, {"SingleTop_untagged"}, 1.05)
	  );
 
  cb.cp().AddSyst(cb, "CMS_eff_mu", "shape", SystMap<channel>::init
	  ({"MuMu"}, 1.00)
	  ({"MuEl"}, 1.00));
  cb.cp().AddSyst(cb, "CMS_iso_mu", "shape", SystMap<channel>::init
	  ({"MuMu"}, 1.00)
	  ({"MuEl"}, 1.00));
  cb.cp().AddSyst(cb, "CMS_eff_e",  "shape", SystMap<channel>::init
	  ({"ElEl"}, 1.00)
	  ({"ElEl"}, 1.00));
  cb.cp().AddSyst(cb, "CMS_iso_e",  "shape", SystMap<channel>::init
	  ({"ElEl"}, 1.00)
	  ({"ElEl"}, 1.00));
  cb.cp().channel({"MuEl"}).AddSyst(cb, "CMS_eff_trigger_MuEl", "shape", SystMap<>::init(1.00));
  cb.cp().channel({"MuMu"}).AddSyst(cb, "CMS_eff_trigger_MuMu", "shape", SystMap<>::init(1.00));
  cb.cp().channel({"ElEl"}).AddSyst(cb, "CMS_eff_trigger_ElEl", "shape", SystMap<>::init(1.00));

  //QCD for signal, TTbar, SingleTop, VV, DY
  cb.cp().process({"Signal"}).AddSyst(cb, "QCDscaleSignal", "shape", SystMap<>::init(1.00));
  cb.cp().process({"VV"}).AddSyst(cb, "QCDscaleVV", "shape", SystMap<>::init(1.00));
  cb.cp().process({"ttV"}).AddSyst(cb, "QCDscalettV", "shape", SystMap<>::init(1.00));
  cb.cp().channel({"MuEl"}).process({"Drell_Yan"}).AddSyst(cb, "QCDscaleDrell_Yan", "shape", SystMap<>::init(1.00));

  cb.cp().AddSyst(cb, "QCDscaleSingleTop", "shape", SystMap<channel, process>::init
	  ({"MuMu"}, {"SingleTop"}, 1.00)
	  ({"MuMu"}, {"SingleTop_untagged"}, 1.00)
	  ({"MuEl"}, {"SingleTop"}, 1.00)
	  ({"ElEl"}, {"SingleTop"}, 1.00)
	  ({"ElEl"}, {"SingleTop_untagged"}, 1.00)
	  );
  cb.cp().AddSyst(cb, "TTbar", "shape", SystMap<channel, process>::init
	  ({"MuMu"}, {"TTbar"}, 1.00)
	  ({"MuMu"}, {"TTbar_untagged"}, 1.00)
	  ({"MuEl"}, {"TTbar"}, 1.00)
	  ({"ElEl"}, {"TTbar"}, 1.00)
	  ({"ElEl"}, {"TTbar_untagged"}, 1.00)
	  );


  //! [part7]
  //extract the shapes from root file 
  for (string chn:chns){
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
	     //GGToX0ToHHTo2B2L2Nu_M800_ElEl_th1shapes.root
        input_dir+"GGToX0ToHHTo2B2L2Nu_M"+mass+"_"+chn+""+postfix+".root",
        "$PROCESS",
        "$PROCESS_$SYSTEMATIC");
   //cb.cp().channel({chn}).process("Signal").ExtractShapes(
    cb.cp().channel({chn}).signals().ExtractShapes(
        input_dir+"GGToX0ToHHTo2B2L2Nu_M"+mass+"_"+chn+""+postfix+".root",
        "$PROCESS",
        "$PROCESS_$SYSTEMATIC");
  }


 //Now delete processes with 0 yield
 cb.FilterProcs([&](ch::Process *p) {
  bool null_yield = !(p->rate() > 0. || BinIsControlRegion(p));
  if (null_yield){
     std::cout << "[Null yield] Removing process with null yield: \n ";
     std::cout << ch::Process::PrintHeader << *p << "\n"; 
     cb.FilterSysts([&](ch::Systematic *s){
       bool remove_syst = (MatchingProcess(*p,*s));
       return remove_syst;
    });
  }
  return null_yield;
 });


  // And convert any shapes in the CRs to lnN:
  // Convert all shapes to lnN at this stage

   //Replacing observation with the sum of the backgrounds (asimov) - nice to ensure blinding
    auto bins = cb.cp().bin_set();
    // For control region bins data (should) = sum of bkgs already
    // useful to be able to check this, so don't do the replacement
    // for these
  if(!real_data){
      for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
          std::cout << " - Replacing data with asimov in bin " << b << "\n";
          cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
            obs->set_shape(cb.cp().bin({b}).backgrounds().GetShape(), true);
          });
        }
  }

  //// and for checking the effect of negative bin contents
  //std::map<std::string, TH1F> before_rebin;
  //std::map<std::string, TH1F> after_rebin;
  //std::map<std::string, TH1F> after_rebin_neg;
  //if (check_neg_bins) {
  //  for (auto b : bins) {
  //    before_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
  //  }
  //}


  //auto rebin = ch::AutoRebin()
  //  .SetBinThreshold(0.)
  //  // .SetBinUncertFraction(0.5)
  //  .SetRebinMode(1)
  //  .SetPerformRebin(true)
  //  .SetVerbosity(1);
  //if(auto_rebin) rebin.Rebin(cb, cb);

  //if(manual_rebin) {
  //  for(auto b : bins) {
  //    std::cout << "Rebinning by hand for bin: " << b <<  std::endl;
  //    cb.cp().bin({b}).VariableRebin(binning[b]);
  //  }
  //}

  //if (check_neg_bins) {
  //  for (auto b : bins) {
  //    after_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
  //    // std::cout << "Bin: " << b << " (before)\n";
  //    // before_rebin[b].Print("range");
  //    // std::cout << "Bin: " << b << " (after)\n";
  //    // after_rebin[b].Print("range");
  //    // Build a sum-of-bkgs TH1 that doesn't truncate the negative yields
  //    // like the CH GetShape does
  //    for (auto p : cb.cp().bin({b}).backgrounds().process_set()) {
  //      TH1F proc_hist;
  //      cb.cp().bin({b}).process({p}).ForEachProc([&](ch::Process *proc) {
  //        proc_hist = proc->ShapeAsTH1F();
  //        proc_hist.Scale(proc->no_norm_rate());
  //        for (int i = 1; i <= proc_hist.GetNbinsX(); ++i) {
  //          if (proc_hist.GetBinContent(i) < 0.) {
  //            std::cout << p << " bin " << i << ": " << proc_hist.GetBinContent(i) << "\n";
  //          }
  //        }
  //      });
  //      if (after_rebin_neg.count(b)) {
  //        after_rebin_neg[b].Add(&proc_hist);
  //      } else {
  //        after_rebin_neg[b] = proc_hist;
  //      }
  //    }
  //    std::cout << "Bin: " << b << "\n";
  //    for (int i = 1; i <= after_rebin[b].GetNbinsX(); ++i) {
  //      double offset = after_rebin[b].GetBinContent(i) - after_rebin_neg[b].GetBinContent(i);
  //      double offset_by_yield = offset / after_rebin[b].GetBinContent(i);
  //      double offset_by_err = offset / after_rebin[b].GetBinError(i);
  //      printf("%-2i offset %-10.4f tot %-10.4f err %-10.4f off/tot %-10.4f off/err %-10.4f\n", i , offset, after_rebin[b].GetBinContent(i), after_rebin[b].GetBinError(i), offset_by_yield, offset_by_err);
  //    }
  //  }
  //}

  // Uncomment this to inject 1 obs event in the last bin of every signal-region
  // category
  // if(!real_data){
  //     for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
  //       std::cout << " - Adjusting data in bin " << b << "\n";
  //         cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
  //           TH1F new_obs = cb.cp().bin({b}).GetObservedShape();
  //           new_obs.SetBinContent(new_obs.GetNbinsX(), 1.);
  //           new_obs.Print("range");
  //           obs->set_shape(new_obs, true);
  //         });
  //       }
  // }

  // At this point we can fix the negative bins
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      std::cout << "[Negative bins] Fixing negative bins for " << p->bin()
                << "," << p->process() << "\n";
      // std::cout << "[Negative bins] Before:\n";
      // p->shape()->Print("range");
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      p->set_shape(std::move(newhist), false);
      // std::cout << "[Negative bins] After:\n";
      // p->shape()->Print("range");
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos) return;
    if (ch::HasNegativeBins(s->shape_u()) || ch::HasNegativeBins(s->shape_d())) {
      std::cout << "[Negative bins] Fixing negative bins for syst" << s->bin()
                << "," << s->process() << "," << s->name() << "\n";
      // std::cout << "[Negative bins] Before:\n";
      // s->shape_u()->Print("range");
      // s->shape_d()->Print("range");
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
      // std::cout << "[Negative bins] After:\n";
      // s->shape_u()->Print("range");
      // s->shape_d()->Print("range");
    }
  });

  cout << "Generating bbb uncertainties...";
  auto bbb = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(true)
    // .SetMergeZeroBins(false)
    .SetPoissonErrors(poisson_bbb);
  for (auto chn : chns) {
    std::cout << " - Doing bbb for channel " << chn << "\n";
    bbb.MergeAndAdd(cb.cp().channel({chn}).process({"ZTT", "QCD", "W", "ZJ", "ZL", "TT", "VV", "Ztt", "ttbar", "EWK", "Fakes", "ZMM", "TTJ", "WJets", "Dibosons"}).FilterAll([](ch::Object const* obj) {
                return BinIsControlRegion(obj);
                }), cb);
  }

  //// And now do bbb for the control region with a slightly different config:
  //auto bbb_ctl = ch::BinByBinFactory()
  //  .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
  //  .SetAddThreshold(0.)
  //  .SetMergeThreshold(0.4)
  //  .SetFixNorm(false)  // contrary to signal region, bbb *should* change yield here
  //  .SetVerbosity(1);
  //// Will merge but only for non W and QCD processes, to be on the safe side
  //bbb_ctl.MergeBinErrors(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion));
  //bbb_ctl.AddBinByBin(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion), cb);
  cout << " done\n";

  //Switch JES over to lnN:
  //cb.cp().syst_name({"CMS_scale_j_13TeV"}).ForEachSyst([](ch::Systematic *sys) { sys->set_type("lnN");});

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  // which is commonly used in the htt analyses
  ch::SetStandardBinNames(cb);
  //! [part8]


  //! [part9]
  // First we generate a set of bin names:


 //Write out datacards. Naming convention important for rest of workflow. We
 //make one directory per chn-cat, one per chn and cmb. In this code we only
 //store the individual datacards for each directory to be combined later, but
 //note that it's also possible to write out the full combined card with CH
  ch::CardWriter writer("shapes/" + output_folder + "/$TAG/$BIN.txt",
                        "shapes/" + output_folder + "/$TAG/$BIN_input.root");
  // We're not using mass as an identifier - which we need to tell the CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});
  writer.SetVerbosity(1);

  writer.WriteCards("cmb", cb);
  for (auto chn : chns) {
    // per-channel
    writer.WriteCards(chn, cb.cp().channel({chn}));
    // And per-channel-category
    writer.WriteCards("GGToX0ToHHTo2B2L2Nu_"+chn, cb.cp().channel({chn}));
  }

  cb.PrintAll();
  cout << " done\n";
}
