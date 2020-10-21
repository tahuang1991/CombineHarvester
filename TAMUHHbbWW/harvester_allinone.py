import ROOT
import os
import sys
import Datacards

masslist = 	[260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
#masslist = 	[260]
channellist =	["ElEl", "MuEl", "MuMu"]
trainlist =	["MTonly", "MTandMT2", "MTandMT2_MJJ"]
trainlist =     ["MTandMT2_MJJ"]
processlist = 	["TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"]
statslist =     ["CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale"]


create_shapes = False
create_datacards = False
auto_rebin = False #Requires create_datacards to be on
plot_sig_bkg = True
create_limits_scripts = False
create_impacts_scripts = False
raw_files_dir = '/afs/cern.ch/user/t/tahuang/public/HHbbWW/'
pwd = '/afs/cern.ch/work/d/daebi/diHiggs/CMSSW_8_1_0/src/CombineHarvester/TAMUHHbbWW/shapes/'
my_raw_files_dir = '/eos/user/d/daebi/20200331_oldNanoAOD/'

directory_dict =	[
	{
		"raw_input"		:	raw_files_dir+"HHbbWW_20201006_NNvsHME_linearized1D_MjjCR_HMEbinsv4/",
                "raw_file_nnstep"       :       "0p1",
		"shapes_output"		:	pwd+"HHbbWW_20201006_NNvsHME_linearized1D_MjjCR_HMEbinsv4_nnstep0p1/",
		"datacard_output"	:	"2D_MjjCR_0p1_autoRebinTrue_th0_fra0p3_emptybins_errorfix_1020"
	},
        {
                "raw_input"             :	raw_files_dir+"HHbbWW_20200814_NNvsHME_linearized1D_MjjCR_HMEbinsv4_1p15/",
                "raw_file_nnstep"       :       "0p1",
                "shapes_output"         :	pwd+"HHbbWW_20200814_NNvsHME_linearized1D_MjjCR_HMEbinsv4_1p15_nnstep0p1/",
                "datacard_output"       :	"2D_MjjCR_0p1_autoRebinTrue_threshold0_fraction0p3_test",
        }




        #{
        #        "raw_input"             :
        #        "raw_file_nnstep"       :
        #        "shapes_output"         :
        #        "datacard_output"       :
        #}
			]

directory = directory_dict[0]

out_prefix = "har_test"  ##Outfile names are {prefix}_{channel}_M{mass}_shapes.root

indir = directory["raw_input"]
out_folder = directory["shapes_output"]
print out_folder

condor_out_folder = pwd+directory["datacard_output"]

os.system("mkdir -p {condor_out_folder}".format(condor_out_folder = condor_out_folder))
os.system("mkdir -p {out_folder}".format(out_folder = out_folder))

os.system("mkdir -p "+condor_out_folder+"/condor")
os.system("mkdir -p "+condor_out_folder+"/condor/output")
os.system("mkdir -p "+condor_out_folder+"/condor/error")
os.system("mkdir -p "+condor_out_folder+"/condor/log")

if create_shapes:
  print "Starting to create shapes"
  ############# Creates the shapes files for CombineHarvester
  for train in trainlist:
    os.system("mkdir -p {out_folder}/{tr}".format(out_folder = out_folder, tr = train))
    for mass in masslist:
      processnames = ["TT","sT","DY","data_untagged","TT_untagged","sT_untagged","ttV","VV", "RadionM{}".format(mass)]
      ################ Might need to change format of infile_name
      #infile_name = indir+"Hhh_FinalBGYield_xsec1pb_NNvsHME_nnout_{tr}_nnstep{nnstep}_nncut0p0_SignalM{m}.root".format(m = mass, tr = train, nnstep = directory["raw_file_nnstep"])
      infile_name = indir+"Hhh_FinalBGYield_xsec1pb_NNvsHME_nnout_{tr}_nnstep{nnstep}_nncut0p0_HME1p15_SignalM{m}.root".format(m = mass, tr = train, nnstep = directory["raw_file_nnstep"])
      #infile_name = indir+"Hhh_FinalBGYield_xsec1pb_NN_nnout_{tr}_nnstep{nnstep}_nncut0p0_SignalM{m}.root".format(m = mass, tr = train, nnstep = directory["raw_file_nnstep"])
      #infile_name = indir+"Hhh_FinalBGYield_xsec1pb_HME_nnout_{tr}cut{nncut}_SignalM{m}.root".format(tr = train, nncut = directory["raw_file_nnstep"], m = mass)
      #infile_name = indir+"Hhh_FinalBGYield_xsec1pb_NN_nnout_{tr}_nnstep{nnstep}_nncut0p0_SignalM{m}.root".format(m = mass, tr = train, nnstep = directory["raw_file_nnstep"])
      infile = ROOT.TFile.Open(infile_name)
      for channel in channellist:
        os.system("mkdir -p {out_folder}/{tr}/{m}".format(tr = train, m = mass, out_folder = out_folder))
        outfile_name = "{out_folder}/{tr}/{m}/{out_prefix}_{ch}_M{m}_shapes.root".format(out_prefix = out_prefix, tr = train, ch = channel, m = mass, out_folder = out_folder)
  
        outfile = ROOT.TFile.Open(outfile_name, "recreate")
        histname_data_obs = "data_obs_{ch}_M{m}".format(ch = channel, m = mass)
        hist_data_obs = infile.Get(histname_data_obs)
        hist_clone_data_obs = hist_data_obs.Clone("data_obs")
        outfile.Write()
        del hist_clone_data_obs
        for i, process in enumerate(processlist):
          if (channel == "MuMu" or channel == "ElEl") and process == "Drell_Yan":
            continue
          if channel == "MuEl" and ("untagged" in process):
            continue
          histname_nominal = processnames[i]+"_"+channel
          hist_nominal = infile.Get(histname_nominal)
          if "Signal" in process:
            hist_nominal.Scale(1e-3/5.0)
          hist_clone_nominal = hist_nominal.Clone(processlist[i])
          outfile.Write()
          del hist_clone_nominal
          if "data" in process:
            continue
          for stat in statslist:
            histname_up = processnames[i]+"_"+channel+"_"+stat+"_up"
            histname_down = processnames[i]+"_"+channel+"_"+stat+"_down"
            histname_sys = process+"_"+stat
            hist_up = infile.Get(histname_up)
            hist_down = infile.Get(histname_down)
            if "Signal" in process:
              hist_up.Scale(1e-3/5.0)
              hist_down.Scale(1e-3/5.0)
            if "QCDscale" in histname_sys and "untagged" in histname_sys:
              histname_sys = process+"_"+stat+process.replace("_untagged", "")
            elif "QCDscale" in histname_sys and "untagged" not in histname_sys:
              histname_sys = process+"_"+stat+process
            if "CMS_eff_trigger" in histname_sys:
              histname_sys = process+"_"+stat+"_"+channel
            hist_clone_up = hist_up.Clone(histname_sys+"Up")
            hist_clone_down = hist_down.Clone(histname_sys+"Down")
            outfile.Write()
            del hist_clone_up
            del hist_clone_down
        outfile.Close()
      infile.Close()


if create_datacards:
  print "Starting to create datacards"
  ################### Create datacards with CombineHarvester
  for tr in trainlist:
    for m in masslist:
      inputdir = directory["shapes_output"]
      outputdir = directory["datacard_output"]
      print "Running M", m, " tr ", tr, " from ", inputdir, " to ", outputdir
      os.system("HHbblvlv --input_folder={inputdir} --output_folder={outputdir} --mass={mass} --training={train} --file_prefix={prefix} --auto_rebin={auto_rebin}".format(inputdir = inputdir, mass = m, train = tr, outputdir = outputdir, prefix = out_prefix, auto_rebin = auto_rebin))





if plot_sig_bkg:
  print "Begin plotting background/signal"
  for tr in trainlist:
    print "Training ", tr
    for mass in masslist:
      print "Mass ", mass
      signal = "RadionM{m}".format(m = mass)

      bgtypes = ["TT", "sT", "VV", "ttV"]
      bgtypes = ["TTbar", "SingleTop", "ttV", "VV"]

      untagged = ["data_unntagged", "TT_untagged", "sT_untagged"]

      branchlist = ["_MC_statistical_", "_CMS_pu_", "_CMS_eff_b_light", "_CMS_eff_trigger_", "_CMS_pdf_", "_CMS_eff_e", "_CMS_iso_mu_", "_CMS_eff_mu_", "_QCDscale_", "_CMS_eff_b_heavy_"]

      channellist = ["MuMu", "MuEl", "ElEl"]

      colorlist = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]

      for channel in channellist:
        filedir = "shapes/"+directory["datacard_output"]+"/"+tr+"/{m}/".format(m = mass)
        print channel
        filename = filedir+'GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV_input.root'.format(tr = tr, m = mass, ch = channel)
        F = ROOT.TFile.Open(filename)
        f = F.Get("GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV".format(tr = tr, m = mass, ch = channel))

        legend = ROOT.TLegend(0.78,0.66,0.88,0.69+len(bgtypes)*.04)
        legend.SetTextSize(0.045)
        legend.SetTextFont(42)
        legend.SetBorderSize(0)

        stack = ROOT.THStack("bgs", "")
        for i in range(len(bgtypes)):
          print bgtypes[i]
          h = f.Get(bgtypes[i])
          h.SetFillColor(colorlist[i])
          h.SetLineColor(colorlist[i])
          stack.Add(h)
          legend.AddEntry(h, bgtypes[i], "f")


        if channel != "MuEl":
          h = f.Get("data_untagged") - f.Get("TTbar_untagged") - f.Get("SingleTop_untagged")
          h.SetFillColor(colorlist[4])
          h.SetLineColor(colorlist[4])
          stack.Add(h)
          legend.AddEntry(h, "DY", "f")
        else:
          h = f.Get("Drell_Yan")
          h.SetFillColor(colorlist[4])
          h.SetLineColor(colorlist[4])
          stack.Add(h)
          legend.AddEntry(h, "DY", "f")


        h_data = f.Get("Signal{m}".format(m = mass))
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerColor(1)
        h_data.SetLineColor(1)
        legend.AddEntry(h_data, "Signal", "p")


        ROOT.gStyle.SetPadLeftMargin(0.13)
        c1 = ROOT.TCanvas("c", "c", 1600, 600)
        c1.Clear()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0, 1, 1.0)
        #pad1.SetBottomMargin(.02)
        pad1.Draw()
        pad1.cd()
        stack.Draw("hist")
        stack.SetTitle(tr+" "+channel+" M{m}".format(m = mass))
        stack.GetXaxis().SetTitle("Linearized DNN output")
        stack.GetYaxis().SetTitle("Entries")
        h_data.Draw("epsame")
        pad1.SetLogy()
        legend.Draw("same")
        os.system("mkdir -p {filedir}plots/".format(filedir = filedir))
        c1.SaveAs(filedir+"plots/{ch}.pdf".format(ch = channel))
        print "Saved to ", filedir, "plots/", channel, ".pdf"








if create_limits_scripts:
  print "Starting to create scripts"
  ################## Create scripts to run combine jobs
  out_folder = pwd+directory["datacard_output"]
  channellist = ["ElEl", "MuEl", "MuMu", "ElEl_MuEl_MuMu"]

  fname_all = out_folder+"/Run_all_scripts.sh"
  cname_all = out_folder+"/Run_all_condors.sh"

  script_all = open(fname_all, "write")
  script_all.write("#!/bin/bash\n")
  script_all.write("cd "+out_folder+"\n")
  script_all.write("eval `scramv1 runtime -sh`\n")

  condor_all = open(cname_all, "write")
  condor_all.write("#!/bin/bash\n")
  condor_all.write("cd "+out_folder+"/condor\n")
  condor_all.write("eval `scramv1 runtime -sh`\n")

  for train in trainlist:
    for mass in masslist:
      for channel in channellist:
        massdir = "{out_folder}/{tr}/{m}/".format(tr = train, m = mass, out_folder = out_folder)
        fname = massdir + "Run_{ch}_M{m}.sh".format(ch = channel, m = mass)
        script_all.write("source "+fname + "\n")
        script = open(fname, "write")
        script.write("#!/bin/bash\n")
        script.write("echo 'start channel {ch}'\n".format(ch = channel))
        script.write("pushd "+massdir+"\n")
        script.write("eval `scramv1 runtime -sh`\n")
        if channel == "ElEl_MuEl_MuMu":
          script.write("# If ElEL MuEl MuMu channel and needs datacard, create\n")
          script.write("if [ ! -f GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt ]; then\n".format(tr = train, m = mass, ch = channel))
          script.write("combineCards.py *.txt &> GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt\n".format(tr = train, m = mass, ch = channel))
          script.write("fi\n\n")
        script.write("# If workspace does not exist, create it once\n")
        script.write("if [ ! -f {out_prefix}_{ch}_M{m}_combine_workspace.root ]; then\n".format(ch = channel, m = mass, out_prefix = out_prefix))
        script.write("text2workspace.py GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt -m {m} -o {out_prefix}_{ch}_M{m}_combine_workspace.root -P HiggsAnalysis.CombinedLimit.DYEstimationCombineModelTao:DYDataDrivenEstimationModelInstance\n".format(ch = channel, m = mass, out_prefix = out_prefix, tr = train))
        script.write("fi\n\n")
        script.write("echo 'finished text2workspace, starting combine' \n")
        script.write("date \n\n")
        script.write("#Run limit\n\n")
        script.write("combine -M AsymptoticLimits -t -1 -m {m} -n {ch}_M{m} {out_prefix}_{ch}_M{m}_combine_workspace.root &> {out_prefix}_{ch}_M{m}.log\n".format(m = mass, ch = channel, out_prefix = out_prefix))
        script.write("popd\n")
        script.write("echo 'finish channel {ch}'\n".format(ch = channel))
        script.write("date \n\n")

        os.system("chmod 775 "+fname)

        condorname = out_folder+"/condor/Condor_{tr}_{ch}_M{m}.sh".format(tr = train, ch = channel, m = mass)
        script_condor = open(condorname, "write")
        script_condor.write("""universe                = vanilla
executable              = {fname}
arguments               = no
output                  = {out_folder}/condor/output/{out_prefix}_{tr}_{ch}_M{m}.$(ClusterId).$(ProcId).out
error                   = {out_folder}/condor/error/{out_prefix}_{tr}_{ch}_M{m}.$(ClusterId).$(ProcId).err
log                     = {out_folder}/condor/log/{out_prefix}_{tr}_{ch}_M{m}.$(ClusterId).log
request_memory          = 4000M
+JobFlavour             = "testmatch"
queue""".format(fname = fname, out_prefix = out_prefix, tr = train, ch = channel, m = mass, out_folder = out_folder))
        os.system("chmod 755 "+condorname)

        condor_all.write("condor_submit "+condorname + "\n")


  os.system("chmod 775 "+fname_all)
  os.system("chmod 775 "+cname_all)


if create_impacts_scripts:
  print "Starting to create impacts scripts"
  ################## Create scripts to run combine jobs
  out_folder = pwd+directory["datacard_output"]
  channellist = ["ElEl", "MuEl", "MuMu", "ElEl_MuEl_MuMu"]

  fname_all = out_folder+"/Run_all_impacts_scripts.sh"
  cname_all = out_folder+"/Run_all_impacts_condors.sh"

  script_all = open(fname_all, "write")
  script_all.write("#!/bin/bash\n")
  script_all.write("cd "+out_folder+"\n")
  script_all.write("eval `scramv1 runtime -sh`\n")

  condor_all = open(cname_all, "write")
  condor_all.write("#!/bin/bash\n")
  condor_all.write("cd "+out_folder+"/condor\n")
  condor_all.write("eval `scramv1 runtime -sh`\n")

  for train in trainlist:
    for mass in masslist:
      for channel in channellist:
        massdir = "{out_folder}/{tr}/{m}/".format(tr = train, m = mass, out_folder = out_folder)
        fname = massdir + "Run_{ch}_M{m}_impacts.sh".format(ch = channel, m = mass)
        script_all.write("source "+fname + "\n")
        script = open(fname, "write")
        script.write("#!/bin/bash\n")
        script.write("echo 'start channel {ch}'\n".format(ch = channel))
        script.write("pushd "+massdir+"\n")
        script.write("eval `scramv1 runtime -sh`\n")
        if channel == "ElEl_MuEl_MuMu":
          script.write("# If ElEL MuEl MuMu channel and needs datacard, create\n")
          script.write("if [ ! -f GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt ]; then\n".format(tr = train, m = mass, ch = channel))
          script.write("combineCards.py *.txt &> GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt\n".format(tr = train, m = mass, ch = channel))
          script.write("fi\n\n")
        script.write("# If workspace does not exist, create it once\n")
        script.write("if [ ! -f {out_prefix}_{ch}_M{m}_combine_workspace.root ]; then\n".format(ch = channel, m = mass, out_prefix = out_prefix))
        script.write("text2workspace.py GGToX0ToHHTo2B2L2Nu_{tr}_M{m}_{ch}_1_13TeV.txt -m {m} -o {out_prefix}_{ch}_M{m}_combine_workspace.root -P HiggsAnalysis.CombinedLimit.DYEstimationCombineModelTao:DYDataDrivenEstimationModelInstance\n".format(ch = channel, m = mass, out_prefix = out_prefix, tr = train))
        script.write("fi\n\n")


        scriptsuffix = "test"
        script.write("echo 'finished text2workspace, starting combine' \n")
        script.write("date \n\n")
        script.write("#Run impacts\n\n")
        script.write("combine --rMin 0 -m {m} -M MultiDimFit {out_prefix}_M{m}_{ch}_combine_workspace.root -t -1 --verbose 3 &> {out_prefix}_M{m}_{ch}_bestfit.log ## bestfit \n".format(m = mass, ch = channel, out_prefix = out_prefix))
        script.write("combine {out_prefix}_{ch}_M{m}_combine_workspace.root -n M{m}_{ch}_rscan_expectS1 -M MultiDimFit --algo grid --points 2000 --rMin -1 --rMax 100 -m {m} --autoRange 1 --squareDistPoiStep -t -1 --fastScan  --expectSignal=1 \n\n".format(m = mass, ch = channel, out_prefix = out_prefix))  ### r scan 
        script.write("combineTool.py -M Impacts --rMax 100 -m {m} -n {out_prefix}_M{m}_{ch}_impacts_S1 -d {out_prefix}_{ch}_M{m}_combine_workspace.root --doInitialFit --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> {out_prefix}_M{m}_{ch}_impacts_{suffix}_step1.log\n".format(m = mass, ch = channel,  suffix =scriptsuffix, out_prefix = out_prefix))
        script.write("combineTool.py -M Impacts --rMax 100 -m {m} -n {out_prefix}_M{m}_{ch}_impacts_S1 -d {out_prefix}_{ch}_M{m}_combine_workspace.root --doFits --parallel 4 --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> {out_prefix}_M{m}_{ch}_impacts_{suffix}_step2.log\n".format(m = mass, ch = channel,  suffix =scriptsuffix, out_prefix = out_prefix))
        script.write("combineTool.py -M Impacts --rMax 100 -m {m} -n {out_prefix}_M{m}_{ch}_impacts_S1 -d {out_prefix}_{ch}_M{m}_combine_workspace.root -o {out_prefix}_M{m}_{ch}_impacts_{suffix}.json --cminDefaultMinimizerStrategy 0 --robustFit 1 --X-rtd MINIMIZER_analytic --expectSignal=1 -t -1 &> {out_prefix}_M{m}_{ch}_impacts_{suffix}_step3.log\n".format(m = mass, ch = channel,  suffix =scriptsuffix, out_prefix = out_prefix))


        script.write("echo 'finish channnel {ch}'\n".format(ch = channel))


        script.write("\n")
        script.write("plotImpacts.py -i {out_prefix}_M{m}_{ch}_impacts_{suffix}.json -o {out_prefix}_M{m}_{ch}_impacts_{suffix} &> {out_prefix}_M{m}_{ch}_impacts_{suffix}_makeplot.log\n".format(m = mass, ch = channel,  suffix = scriptsuffix, out_prefix = out_prefix))
        script.write("\n")
        #script.write("cp {out_prefix}_M{m}_{ch}_impacts_{suffix}.pdf {plotdir} \n".format(m = mass, ch = channel,  suffix = scriptsuffix, plotdir = plotdir, out_prefix = out_prefix))


        script.write("popd\n")
        script.write("date \n\n")



        os.system("chmod 775 "+fname)

        condorname = out_folder+"/condor/Condor_{tr}_{ch}_M{m}_impacts.sh".format(tr = train, ch = channel, m = mass)
        script_condor = open(condorname, "write")
        script_condor.write("""universe                = vanilla
executable              = {fname}
arguments               = no
output                  = {out_folder}/condor/output/{out_prefix}_{tr}_{ch}_M{m}_impacts.$(ClusterId).$(ProcId).out
error                   = {out_folder}/condor/error/{out_prefix}_{tr}_{ch}_M{m}_impacts.$(ClusterId).$(ProcId).err
log                     = {out_folder}/condor/log/{out_prefix}_{tr}_{ch}_M{m}_impacts.$(ClusterId).log
request_memory          = 4000M
+JobFlavour             = "testmatch"
queue""".format(fname = fname, out_prefix = out_prefix, tr = train, ch = channel, m = mass, out_folder = out_folder))
        os.system("chmod 755 "+condorname)

        condor_all.write("condor_submit "+condorname + "\n")


  os.system("chmod 775 "+fname_all)
  os.system("chmod 775 "+cname_all)
