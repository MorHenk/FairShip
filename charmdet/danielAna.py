
#import yep
import ROOT,os,time,sys,operator,atexit
#ROOT.gROOT.ProcessLine('typedef std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::vector<MufluxSpectrometerHit*>>>> nestedList;')

from decorators import *
import __builtin__ as builtin
ROOT.gStyle.SetPalette(ROOT.kGreenPink)
PDG = ROOT.TDatabasePDG.Instance()
# -----Timer--------------------------------------------------------
timer = ROOT.TStopwatch()
# stop printing errors
ROOT.gErrorIgnoreLevel = ROOT.kBreak
from argparse import ArgumentParser
import shipunit as u
import rootUtils as ut
from array import array


parser = ArgumentParser()
parser.add_argument("-f", "--files", dest="listOfFiles", help="list of files comma separated", required=True)
parser.add_argument("-l", "--fileCatalog", dest="catalog", help="list of files in file", default=False)
parser.add_argument("-c", "--cmd", dest="command", help="command to execute", default="")
parser.add_argument("-d", "--Display", dest="withDisplay", help="detector display", default=True)
parser.add_argument("-e", "--eos", dest="onEOS", help="files on EOS", default=False)
parser.add_argument("-u", "--update", dest="updateFile", help="update file", default=False)
parser.add_argument("-i", "--input", dest="inputFile", help="input histo file", default='residuals.root')
parser.add_argument("-g", "--geofile", dest="geoFile", help="input geofile", default='')
###parser.add_argument("-s", "--smearing", dest="MCsmearing", help="additional MC smearing", default=MCsmearing)





options = parser.parse_args()
###MCsmearing = options.MCsmearing
fnames = []
if options.catalog:
    tmp = open(options.listOfFiles)
    for x in tmp.readlines():
        fname = x.replace('\n','')
        if fname.find("root")<0:continue
        f=ROOT.TFile.Open(fname)
        sTree = f.cbmsim
        if not sTree.GetBranch("FitTracks"): 
            print "does not contain FitTracks",fname
            f.Close()
            continue
        fnames.append(fname)
        fnames.append(x.replace('\n',''))
    tmp.close()
else:
    fnames = options.listOfFiles.split(',')
fname = fnames[0]
if options.updateFile:
    f=ROOT.TFile(fname,'update')
    sTree=f.Get('cbmsim')
    if not sTree: 
        print "Problem with updateFile",f
        exit(-1)
else:
    sTree = ROOT.TChain('cbmsim')
    for f in fnames: 
        print "add ",f
        if options.onEOS: sTree.Add(os.environ['EOSSHIP']+f)
        else:             sTree.Add(f)


rnames = []
for fname in fnames:
    rnames.append(fname[fname.rfind('/')+1:])
rname = rnames[0]
#sTree.SetMaxVirtualSize(300000)

from ShipGeoConfig import ConfigRegistry
if options.geoFile=="":
    ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/charm-geometry_config.py")
else:
    from ShipGeoConfig import ConfigRegistry
    from rootpyPickler import Unpickler
#load Shipgeo dictionary
    fgeo = ROOT.TFile.Open(options.geoFile)
    upkl    = Unpickler(options.geoFile)
    ShipGeo = upkl.load('ShipGeo')
builtin.ShipGeo = ShipGeo
import charmDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for creating VMC field
rtdb = run.GetRuntimeDb()
modules = charmDet_conf.configure(run,ShipGeo)
# -----Create geometry----------------------------------------------
run.Init()
sGeo = ROOT.gGeoManager
nav = sGeo.GetCurrentNavigator()
top = sGeo.GetTopVolume()
top.SetVisibility(0)
#if options.withDisplay:
#    try: top.Draw('ogl')
#    except: pass

saveGeofile = False
import saveBasicParameters
if saveGeofile:
    run.CreateGeometryFile("muflux_geofile.root")
# save ShipGeo dictionary in geofile
    saveBasicParameters.execute("muflux_geofile.root",ShipGeo)



    
def drawRawDrifttimes():
    ROOT.helloDaniel()
    sTree.Draw("Digi_MufluxSpectrometerHits.fdigi")


def localAna():
    f=ROOT.TFile.Open(options.listOfFiles)
    t=ROOT.TTreeReader("cbmsim",f)
    ROOT.dtAnaChain(t)
    #ROOT.FilterDTSpectrum(t)



def localAnaEvent(i):
    f=ROOT.TFile.Open(options.listOfFiles)
    t=ROOT.TTreeReader("cbmsim",f)
    ROOT.dtAnaChain(t,i)


def testEP():
    #m=ROOT.MufluxSpectrometer()
    #a=ROOT.TVector3()
    #b=ROOT.TVector3()
    #m.TubeEndPointsSurvey(10002001,a,b)
    #m.TubeEndPointsSurvey(10002002,a,b)
    #m.TubeEndPointsSurvey(10002003,a,b)
    m=ROOT.MufluxSpectrometerDTSurvey()
    m.Init()
    m.Test(11102001)


    
    
#localAna()
#testEP()
