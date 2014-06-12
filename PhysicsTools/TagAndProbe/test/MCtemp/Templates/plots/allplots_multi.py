## This file is the same as allplots.py, except that it uses multiprocessing
## to make better use of machines with multiple cores

try:
  ## the normal way to import rootplot
  from rootplot import plot, plotmpl
  from rootplot.core import report_progress
except ImportError:
  ## special import for CMSSW installations of rootplot
  from PhysicsTools.PythonAnalysis.rootplot import plot, plotmpl
  from PhysicsTools.PythonAnalysis.rootplot.core import report_progress
import ROOT
import multiprocessing as multi

import os
os.chdir('..')  # return to the directory with the ROOT files

calls = []

calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hcounter', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hcounter.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'dump', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/dump.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0To0.8_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0.8To1.4442_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0.8To1.4442_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_50To200_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_50To200_0To0.8_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0.8To1.4442_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0.8To1.4442_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.4442To1.566_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.4442To1.566_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_1.566To2_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_1.566To2_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0To0.8_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_2To2.5_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_2To2.5_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_1.4442To1.566_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_1.4442To1.566_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_0To0.8_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_0To0.8_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_10To15_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_10To15_0To0.8_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_0To0.8_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_40To50_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_40To50_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_30To40_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_30To40_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_20To30_1.566To2_Fail', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_20To30_1.566To2_Fail.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_2To2.5_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_2To2.5_Pass.png')
""")
calls.append("""
canvas, objects = plot('mcTruemadgraph_Medium.root', 'hMass_15To20_0To0.8_Pass', legend_entries='mcMatched-TauRemoved(dR<0.2)', draw='hist')
canvas.SaveAs('plots/hMass_15To20_0To0.8_Pass.png')
""")


queue = multi.JoinableQueue()
qglobals = multi.Manager().Namespace()
qglobals.nfinished = 0
qglobals.ntotal = len(calls)
for call in calls:
    queue.put(call)

def qfunc(queue, qglobals):
    from Queue import Empty
    while True:
        try: mycall = queue.get(timeout=5)
        except (Empty, IOError): break
        exec(mycall)
        ROOT.gROOT.GetListOfCanvases().Clear()
        qglobals.nfinished += 1
        report_progress(qglobals.nfinished, qglobals.ntotal, 
                        'plots', 'png')
        queue.task_done()

for i in range(24):
    p = multi.Process(target=qfunc, args=(queue, qglobals))
    p.daemon = True
    p.start()
queue.join()
report_progress(len(calls), len(calls), 'plots', 'png')
print ''
