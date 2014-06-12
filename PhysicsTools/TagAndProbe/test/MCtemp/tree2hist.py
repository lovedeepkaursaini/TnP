# Configuration file for tree2hists
# Created Jun 04, 2013.
try:
    ## the normal way to import from rootplot
    from rootplot.tree2hists import RootTree, Plot
except ImportError:
    ## special import for CMSSW installations of rootplot
    from PhysicsTools.PythonAnalysis.rootplot.tree2hists import RootTree, Plot
from array import array      # to allow making Float_t arrays for ROOT hists
from math import pi
from ROOT import TH1F, TH2F  # import other kinds of hists as neeeded

list_of_files = [
RootTree("GsfElectronToIdMedium/fitter_tree", 
fileName=
"/hdfs/store/user/kaur/mrgdTnPTreeElec_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Aug24.root"
#"/hdfs/store/user/kaur/mrgdTnPTreeElec_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraphSummer12_DR53X-PU_S10_START53_V7A-v1_Aug24.root"
, scale=1.0, cuts="probe_sc_abseta < 0.8 && ( probe_gsfEle_pt > 10. && probe_gsfEle_pt <15.) && (mass > 60. && mass < 120.)")]

#cuts="probe_sc_abseta < 0.8 && ( probe_gsfEle_pt > 10. && probe_gsfEle_pt <15.) && tag_sc_abseta < 0.8 && (tag_gsfEle_pt > 10. && tag_gsfEle_pt < 15.) && (mass > 60. && mass < 120.)"

output_filename = "Hists_TnP_powheg.root"

cut_for_all_files = ""

# All plots are made for each "cut set".
# A "cut set" is 3 things: folder name to store hists in, string to add to hist titles, and cuts for these hists.
# Let cut_sets = [] to make all plots.
cut_sets = [    ]

# Define histograms to plot
bins_et     = array("f", [15.0, 20.0, 30.0, 50.0, 80.0, 120.0]) # example custom bins
list_of_plots = [
    Plot("event_nPV"           , TH1F("event_nPV"           , "event_nPV", 20, 0., 40.0)),
    Plot("probe_dRTau"        ,TH1F("probe_dRTau_close"        , "probe_dRTau close",100,0.,5.0)),
    Plot("probe_dRTau"        ,TH1F("probe_dRTau_all"        , "probe_dRTau all",100,0.,110.)),
    #Plot("probe_dRTau"        ,TH1F("probe_dRTau"        , "probe_dRTau PASS",100,0.,110.),"probe_isWPMedium ==1"),
    #Plot("tag_dRTau"        ,TH1F("tag_dRTau"        , "tag_dRTau PASS",100,0.,5.), "probe_isWPMedium ==1"),
    #Plot("probe_dRTau"        ,TH1F("probe_dRTau F"        , "probe_dRTau FAIL",100,0.,5.), "probe_isWPMedium !=1"),
    #Plot("tag_dRTau"        ,TH1F("tag_dRTau F"        , "tag_dRTau FAIL",100,0.,5.), "probe_isWPMedium !=1"),

]


