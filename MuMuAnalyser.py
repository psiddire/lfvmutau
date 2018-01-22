from MuMuTree import MuMuTree
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.WARNING)
import os
import ROOT
import math
import glob
import array
import mcCorrections
import baseSelections as selections
from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.PlotTools.decorators import memo_last
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, cos

@memo
def getVar(name, var):
    return name+var

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI

def deltaR(phi1, phi2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);

pucorrector = mcCorrections.make_puCorrector('singlem', None)

class MuMuAnalyser(MegaBase):
    tree = 'mm/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        logging.debug('MuMuAnalyser constructor')
        self.channel='MM'
        target = os.path.basename(os.environ['megatarget'])
        self.is_data = target.startswith('data_')
        print "data", self.is_data
#self.is_embedded = ('Embedded' in target)
        self.is_mc = not (self.is_data)
        self.is_DY = bool('JetsToLL' in target)
        self.is_W = bool('JetsToLNu' in target)

        super(MuMuAnalyser, self).__init__(tree, outfile, **kwargs)
        self.tree = MuMuTree(tree)
        self.out=outfile
        self.histograms = {}
        self.histo_locations = {}
        self.hfunc   = {}
        ROOT.TH1.AddDirectory(True)

        self.DYweight = {
            "Njet0M-10to50" : 11.79816076,#0.641343812,#28.41047665, 
            "Njet0M-50" : 0.760199819,#0.760199819,#2.416349945, 
            "Njet1M-10to50" : 0.232050075,#0.012614159,#0.28919893, 
            "Njet1M-50" : 0.245670161,#0.245091596,#0.644408977, 
            "Njet2M-10to50" : 0.208336374,#0.011325091,#0.37013562, 
            "Njet2M-50" : 0.249449334,#0.250170259,#0.628679429, 
            "Njet3M-50" : 0.248033247,#0.257234053,#0.565321637, 
            "Njet4M-50" : 0.20138031#0.210895346#0.380039915
            }

        self.Wweight = {
            0 : 11.15941847,#0.709390278,#19.67132617,
            1 : 2.989895189,#0.190063899,#8.291202067,
            2 : 0.920734871,#0.058529964,#2.80761307,
            3 : 0.302136585,#0.019206445,#0.751832259,
            4 : 0.302593339#0.01923548#0.657715041
            }

    def event_weight(self, row, sys_shifts, weight_to_use):

        mcweight = weight_to_use
        mcweight = mcweight*pucorrector(row.nTruePU)

        if self.is_DY:
            if row.numGenJets < 5:                                                              
                if row.numGenJets == 0 and row.m1_m2_Mass < 50:
                    mcweight = mcweight*self.DYweight["Njet0M-10to50"]#*dyweight
                elif row.numGenJets == 0 and row.m1_m2_Mass > 50:
                    mcweight = mcweight*self.DYweight["Njet0M-50"]#*dyweight
                elif row.numGenJets == 1 and row.m1_m2_Mass < 50:
                    mcweight = mcweight*self.DYweight["Njet1M-10to50"]#*dyweight
                elif row.numGenJets == 1 and row.m1_m2_Mass > 50:
                    mcweight = mcweight*self.DYweight["Njet1M-50"]#*dyweight
                elif row.numGenJets == 2 and row.m1_m2_Mass < 50:
                    mcweight = mcweight*self.DYweight["Njet2M-10to50"]#*dyweight
                elif row.numGenJets == 2 and row.m1_m2_Mass > 50:
                    mcweight = mcweight*self.DYweight["Njet2M-50"]#*dyweight
                elif row.numGenJets == 3 and row.m1_m2_Mass > 50:
                    mcweight = mcweight*self.DYweight["Njet3M-50"]#*dyweight
                elif row.numGenJets == 4 and row.m1_m2_Mass > 50:
                    mcweight = mcweight*self.DYweight["Njet4M-50"]#*dyweight
            else:                                            
                mcweight = mcweight*self.DYweight[0]  

        if self.is_W:
            if row.numGenJets < 5:
                mcweight = mcweight*self.Wweight[row.numGenJets]
            else:
                mcweight = mcweight*self.Wweight[0]

        if self.is_data:
            mcweight = 1.

        weights = {'': mcweight} 
        
        return weights

    def begin(self):
        folder = ['mm','Iso1','Iso0p5','Iso0p25','Iso0p15','Iso0p1']
        for f in folder:
            
       # f='mt' this is the name of the folder with histo in the output file
            self.book(f,"m1Pt", "m1 Pt", 200, 0, 200 )
            self.book(f,"m1Eta", "m1 Eta", 100, -10, 10 )
            self.book(f,"m1Phi", "m1 Phi", 100, -10, 10 )
            self.book(f,"m2Pt", "m2 Pt", 200, 0, 200 )
            self.book(f,"m2Eta", "m2 Eta", 100, -10, 10 )
            self.book(f,"m2Phi", "m2 Phi", 100, -10, 10 )
            self.book(f,"m1_m2_Mass", "m1 m2 Mass", 200, 0, 200 )
           

        for key, value in self.histograms.iteritems():
            location = os.path.dirname(key)
            name     = os.path.basename(key)
            if location not in self.histo_locations:
                self.histo_locations[location] = {name : value}
            else:
                self.histo_locations[location][name] = value
                
    def fill_histos(self, folder_str, row, weight, filter_label = ''):
        
        for attr, value in self.histo_locations[folder_str].iteritems():
            name = attr

            if value.InheritsFrom('TH2'):
                if attr in self.hfunc:
                    try:
                        result, out_weight = self.hfunc[attr](row, weight)
                    except Exception as e:
                        raise RuntimeError("Error running function %s. Error: \n\n %s" % (attr, str(e)))
                    r1, r2 = result
                    if out_weight is None:
                        value.Fill( r1, r2 ) #saves you when filling NTuples!
                    else:
                        value.Fill( r1, r2, out_weight )
                else:
                    attr1, attr2 = split(attr, '_vs_')
                    v1 = getattr(row,attr1)
                    v2 = getattr(row,attr2)
                    value.Fill( v2, v1, weight ) if weight is not None else value.Fill( v2, v1 )
            else:
                if attr in self.hfunc:
                    try:
                        result, out_weight = self.hfunc[attr](row, weight)
                    except Exception as e:
                        raise RuntimeError("Error running function %s. Error: \n\n %s" % (attr, str(e)))
                    if out_weight is None:
                        value.Fill( result ) #saves you when filling NTuples!
                    else:
                        value.Fill( result, out_weight )
                else:
                    #print attr, weight
                    value.Fill( getattr(row,attr), weight ) if weight is not None else value.Fill( getattr(row,attr) )


        return None

    
    def process(self):
        sys_shifts = []
        for row in self.tree:
            if int(row.m1Charge*row.m2Charge)==1:
                continue

            if row.singleIsoMu27Pass == 0:# and row.singleMu27Pass == 0:
                continue      

            if row.bjetCISVVeto30Medium!=0 : continue
            if not selections.muSelection(row, 'm1') or not selections.muSelection(row, 'm2'): continue
            if abs(91.19 - row.m1_m2_Mass) > 30 : continue
            if not row.m1PFIDLoose == 1 or not row.m2PFIDLoose == 1: continue

            dirnames = ['mm']
            weight_to_use = row.GenWeight #it should be changed when using corrections
            weight_map = self.event_weight(row, sys_shifts, weight_to_use)
            if weight_map[''] == 0 : continue
            weight_to_use = weight_map['']
            
            if row.m1RelPFIsoDBDefaultR04 < 1 and row.m2RelPFIsoDBDefaultR04 < 1:
                dirnames.append('Iso1')
            if row.m1RelPFIsoDBDefaultR04 < 0.5 and row.m2RelPFIsoDBDefaultR04 < 0.5:
                dirnames.append('Iso0p5')
            if row.m1RelPFIsoDBDefaultR04 < 0.25 and row.m2RelPFIsoDBDefaultR04 < 0.25:
                dirnames.append('Iso0p25')
            if row.m1RelPFIsoDBDefaultR04 < 0.15 and row.m2RelPFIsoDBDefaultR04 < 0.15:
                dirnames.append('Iso0p15')
                #if oldevent[0] == row.evt: print "eventID", row.evt, "m1Pt", row.m1Pt, "m2Pt", row.m2Pt, "row.m1Phi", row.m1Phi, "row.m2Phi", row.m2Phi, "row.m1Eta", row.m1Eta, "row.m2Eta", row.m2Eta, "deltaR", row.m1_m2_DR, "oldevent deltaR", oldevent[8], "m1_m2_Mass", row.m1_m2_Mass, "Mass", oldevent[7]
                #oldevent = (row.evt,row.m1Pt,row.m2Pt,row.m1Phi,row.m2Phi,row.m1Eta,row.m2Eta,row.m1_m2_Mass,row.m1_m2_DR)
            if row.m1RelPFIsoDBDefaultR04 < 0.1 and row.m2RelPFIsoDBDefaultR04 < 0.1:
                dirnames.append('Iso0p1')
            for dirname in dirnames:
                self.fill_histos(dirname, row, 1.0)
            

    def finish(self):
        self.write_histos()
