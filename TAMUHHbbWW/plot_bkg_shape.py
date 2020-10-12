import os
import ROOT

filename = "/afs/cern.ch/work/d/daebi/diHiggs/CMSSW_8_1_0/src/CombineHarvester/TAMUHHbbWW/shapes/2D_MjjCR_0p1_autoRebinTrue_threshold0_frac0p8/MTandMT2_MJJ/450/GGToX0ToHHTo2B2L2Nu_MTandMT2_MJJ_M450_MuMu_1_13TeV_input.root"

f = ROOT.TFile(filename)

event = f.Get("GGToX0ToHHTo2B2L2Nu_MTandMT2_MJJ_M450_MuMu_1_13TeV")
