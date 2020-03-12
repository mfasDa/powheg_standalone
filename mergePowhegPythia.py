#! /usr/bin/env python3
from ROOT import TFile, gDirectory
import sys

def readScaledHistos(inputfile) :
    result = [];
    reader = TFile.Open(inputfile, "READ")
    # get the number of events
    norm = reader.Get("hNEventsPos")
    for obj in gDirectory.GetListOfKeys():
        histo = obj.ReadObj()
        histname = histo.GetName()
        if not "hNEvent" in histname and not "hSumWeights" in histname:
            histo.Scale(norm.GetBinContent(1));
        histo.SetDirectory(None)
        result.append(histo);
    reader.Close()
    return result;

def mergePowhegPythia(outputfile, inputlist):
    resulthists = []
    for infile in inputlist:
        histos = readScaledHistos(infile)
        if not len(resulthists):
            resulthists = histos
        else:
            for h in histos:
                target = [x for x in resulthists if x.GetName() == h.GetName()][0]
                target.Add(h)
    
    #Scale to total number of events
    normhist = [x for x in resulthists if x.GetName() == "hNEventsPos"][0]
    spechists = [x for x in resulthists if not "hNEvents" in x.GetName() and not "hSumWeights" in x.GetName()]
    for sh in spechists:
        sh.Scale(1./normhist.GetBinContent(1))

    # Write to file
    writer = TFile.Open(outputfile, "RECREATE")
    for hist in histos:
        hist.Write()
    writer.Close()

if __name__ == "__main__":
    mergePowhegPythia(sys.argv[1], sys.argv[2:])