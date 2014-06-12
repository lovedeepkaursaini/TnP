## This file contains all the necessary calls to the rootplot API to produce
## the same set of plots that were created from the command-line.

## You can use this file to intercept the objects and manipulate them before
## the figure is saved, making any custom changes that are not possible from
## the command-line.

## 'objects' is a python dictionary containing all the elements used in the
## plot, including 'hists', 'legend', etc.
##   ex: objects['hists'] returns a list of histograms

try:
  ## the normal way to import rootplot
  from rootplot import plot, plotmpl
except ImportError:
  ## special import for CMSSW installations of rootplot
  from PhysicsTools.PythonAnalysis.rootplot import plot, plotmpl

import os
os.chdir('..')  # return to the directory with the ROOT files

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hcounter', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hcounter.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'dump', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/dump.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0To0.8_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0.8To1.4442_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0To0.8_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0.8To1.4442_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.4442To1.566_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.566To2_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0To0.8_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_2To2.5_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.4442To1.566_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0To0.8_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0To0.8_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0To0.8_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.566To2_Fail.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_2To2.5_Pass.png')

canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0To0.8_Pass.png')
