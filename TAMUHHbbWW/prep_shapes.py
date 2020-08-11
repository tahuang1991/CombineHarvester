import ROOT
import os
import sys

masslist =      [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
channellist =   ["ElEl", "MuEl", "MuMu"]
trainlist =     ["MTonly", "MTandMT2", "MTandMT2_MJJ"]
processlist =   ["TTbar","SingleTop","Drell_Yan","data_untagged","TTbar_untagged","SingleTop_untagged", "ttV","VV","Signal"]
statslist =     ["CMS_eff_b_heavy","CMS_eff_b_light","CMS_pu", "CMS_pdf", "CMS_eff_trigger","CMS_eff_e","CMS_eff_mu","CMS_iso_mu","QCDscale"]

indir = "raw_input/"

outdir = "shapes/"

out_prefix = "nnstep0p04"

indir = indir+"HHbbWW_20200807_NNoutput_MjjCR_NNcutstudy1D/"
out_folder = outdir+"HHbbWW_20200401_NNoutput_MjjCR_NNcutstudy1D/"
out_folder = outdir+"HHbbWW_20200807_NNoutput_MjjCR_NNcutstudy1D/"

for train in trainlist:
  os.system("mkdir -p {out_folder}/{tr}".format(out_folder = out_folder, tr = train))
  for mass in masslist:
    processnames = ["TT","sT","DY","data_untagged","TT_untagged","sT_untagged","ttV","VV", "RadionM{}".format(mass)]
    infile_name = indir+"Hhh_FinalBGYield_xsec1pb_NN_nnout_{tr}_nnstep0p04_nncut0p0_SignalM{m}.root".format(m = mass, tr = train)
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
