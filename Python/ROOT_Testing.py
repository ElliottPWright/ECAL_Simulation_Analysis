import ROOT
import numpy as np


x = np.linspace(1, 100, num=100)
y = 2*x + 1 

graph = ROOT.TGraph(len(x), x, y)
graph.GetXaxis().SetTitle("X")
graph.GetYaxis().SetTitle("Y")
graph.SetMarkerStyle(20)
graph.SetMarkerColor(ROOT.kBlue)
graph.SetName("mygraph")
graph.Write()