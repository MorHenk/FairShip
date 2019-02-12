#include "MufluxReco.h"
#include <TROOT.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TString.h"
#include "TObject.h"
#include "RPCTrack.h"
#include "RKTrackRep.h"
#include "MufluxSpectrometerHit.h"

TVector3* parallelToZ = new TVector3(0., 0., 1.);
TVector3* NewPosition = new TVector3(0., 0., 0.);
// -----   Standard constructor   ------------------------------------------ 
MufluxReco::MufluxReco() {}
// -----   Standard constructor   ------------------------------------------ 
MufluxReco::MufluxReco(TTreeReader* t)
{
  xSHiP = t;
  MCdata = false;
  if (xSHiP->GetTree()->GetBranch("MCTrack")){MCdata=true;}
  std::cout << "MufluxReco initialized for "<<xSHiP->GetEntries(true) << " events "<<std::endl;
  xSHiP->ls();
  FitTracks = 0;
  TrackInfos = 0;
  RPCTrackY = 0;
  RPCTrackX = 0;
  TTree* fChain = xSHiP->GetTree();
  fChain->SetBranchAddress("FitTracks", &FitTracks, &b_FitTracks);
  fChain->SetBranchAddress("TrackInfos", &TrackInfos, &b_TrackInfos);
  fChain->SetBranchAddress("RPCTrackY", &RPCTrackY, &b_RPCTrackY);
  fChain->SetBranchAddress("RPCTrackX", &RPCTrackX, &b_RPCTrackX);
}

/*
MufluxReco::getKeys(std::unordered_map<int,MufluxSpectrometerHit*>> map){
 std::vector<int> keys;
 keys.reserve(map.size());
 std::vector<std::string> vals;
 vals.reserve(allHits.size());
 for(auto kv : allHits) {
    keys.push_back(kv.first);
    vals.push_back(kv.second);  
}*/

Double_t MufluxReco::extrapolateToPlane(genfit::Track* fT,Float_t z, TVector3& pos, TVector3& mom){
// etrapolate to a plane perpendicular to beam direction (z)
  Double_t rc = -1;
  auto fst = fT->getFitStatus();
  if (fst->isFitConverged()){
   if (z > cuts["firstDTStation_z"]-10 and z < cuts["lastDTStation_z"] + 10){
// find closest measurement
    Float_t mClose = 0;
    Float_t mZmin = 999999.;
    for (Int_t m=0;m<fT->getNumPointsWithMeasurement();m++) {
      auto st  = fT->getFittedState(m);
      auto Pos = st.getPos();
      if (TMath::Abs(z-Pos.z())<mZmin){
       mZmin = TMath::Abs(z-Pos.z());
       mClose = m;
      }
    }
    genfit::StateOnPlane fstate =  fT->getFittedState(mClose);
    NewPosition->SetXYZ(0., 0., z);
    Int_t pdgcode = -int(13*fstate.getCharge());
    genfit::RKTrackRep* rep      = new genfit::RKTrackRep( pdgcode );
    genfit::StateOnPlane* state   = new genfit::StateOnPlane(rep);
    auto Pos = fstate.getPos();
    auto Mom = fstate.getMom();
    rep->setPosMom(*state,Pos,Mom);
    rc = rep->extrapolateToPlane(*state, *NewPosition, *parallelToZ );
    pos = (state->getPos());
    mom = (state->getMom());
    delete rep;
    delete state;
   }else{
    genfit::StateOnPlane fstate = fT->getFittedState(0);
    if ( z > cuts["lastDTStation_z"]){
     Int_t nmeas = fT->getNumPointsWithMeasurement();
     Int_t M = TMath::Min(nmeas-1,30);
     fstate = fT->getFittedState(M);
    }
    pos = fstate.getPos();
    mom = fstate.getMom();
// use linear extrap 
    Float_t lam = (z-pos[2])/mom[2];
    pos[2]=z;
    pos[0]=pos[0]+lam*mom[0];
    pos[1]=pos[1]+lam*mom[1];
    rc = 0;
  }
 }
  return rc;
}


StringVecIntMap MufluxReco::countMeasurements(TrackInfo* trInfo){
 StringVecIntMap mStatistics = StringVecIntMap();
 for (Int_t n=0;n<trInfo->N();n++) {
    Int_t detID = trInfo->detId(n);
    MufluxSpectrometerHit* hit = new MufluxSpectrometerHit(detID,0);
    auto info = hit->StationInfo();
    delete hit;
    Int_t s=info[0]; Int_t v=info[1]; Int_t p=info[2]; Int_t l=info[3]; Int_t channelNr=info[5]; 
    if (trInfo->wL(n) <0.1 && trInfo->wR(n) <0.1){ continue;}
    if (v != 0){ 
       mStatistics["uv"].push_back(detID);
       if (v == 1) { mStatistics["u"].push_back(detID);}
       if (v == 2) { mStatistics["v"].push_back(detID);}
    }else{
     mStatistics["xAll"].push_back(detID);
     TString x = "x";x+=s;
     mStatistics[x.Data()].push_back(detID);
    }
    if (s > 2){  
       mStatistics["xDown"].push_back(detID);
    }else{
       mStatistics["xUp"].push_back(detID);
    }
 }
 return mStatistics;
}

void MufluxReco::trackKinematics(Float_t chi2UL, Int_t nMax){
 Int_t N = xSHiP->GetEntries(true);
 TTree* sTree = xSHiP->GetTree();
 if (nMax<0){nMax=N;}
 xSHiP->Restart();
 gROOT->cd();
 std::cout<< "fill trackKinematics: "<< N <<std::endl;
/* TTreeReader framework cannot deal with unsplit TClonesArray, too badFiles
 TTreeReaderValue <TClonesArray> FitTracks(*xSHiP, "FitTracks");
 TTreeReaderValue <TClonesArray> RPCTrackX(*xSHiP, "RPCTrackX");
 TTreeReaderValue <TClonesArray> RPCTrackY(*xSHiP, "RPCTrackY");
 TTreeReaderValue <TClonesArray> TrackInfos(*xSHiP, "TrackInfos");*/

 TH1D* h_chi2 =  (TH1D*)(gDirectory->GetList()->FindObject("chi2"));
 TH1D* h_Nmeasurements = (TH1D*)(gDirectory->GetList()->FindObject("Nmeasurements"));
 TH2D* h_ppt = (TH2D*)(gDirectory->GetList()->FindObject("p/pt"));
 TH2D* h_ppx = (TH2D*)(gDirectory->GetList()->FindObject("p/px"));
 TH2D* h_xy = (TH2D*)(gDirectory->GetList()->FindObject("xy"));
 TH2D* h_pxpy = (TH2D*)(gDirectory->GetList()->FindObject("pxpy"));

 TH1D* h_chi2mu =  (TH1D*)(gDirectory->GetList()->FindObject("chi2mu"));
 TH1D* h_Nmeasurementsmu = (TH1D*)(gDirectory->GetList()->FindObject("Nmeasurementsmu"));
 TH2D* h_pptmu = (TH2D*)(gDirectory->GetList()->FindObject("p/ptmu"));
 TH2D* h_ppxmu = (TH2D*)(gDirectory->GetList()->FindObject("p/pxmu"));
 TH2D* h_xymu = (TH2D*)(gDirectory->GetList()->FindObject("xymu"));
 TH2D* h_pxpymu = (TH2D*)(gDirectory->GetList()->FindObject("pxpymu"));

 TH2D* h_p1p2 = (TH2D*)(gDirectory->GetList()->FindObject("p1/p2"));
 TH2D* h_pt1pt2 = (TH2D*)(gDirectory->GetList()->FindObject("pt1/pt2"));

 Int_t nx = 0;
 while (nx<nMax){
   sTree->GetEvent(nx);
   nx+=1;
   Int_t Ntracks = FitTracks->GetEntries();
   for (Int_t k=0;k<Ntracks;k++) {
     genfit::Track* aTrack = (genfit::Track*)FitTracks->At(k);
     auto fitStatus   = aTrack->getFitStatus();
// track quality
     if (!fitStatus->isFitConverged()){continue;}
     TrackInfo* info = (TrackInfo*)TrackInfos->At(k);
     StringVecIntMap hitsPerStation = countMeasurements(info);
     if (hitsPerStation["x1"].size()<2){ continue;}
     if (hitsPerStation["x2"].size()<2){ continue;}
     if (hitsPerStation["x3"].size()<2){ continue;}
     if (hitsPerStation["x4"].size()<2){ continue;}
     auto chi2 = fitStatus->getChi2()/fitStatus->getNdf();
     auto fittedState = aTrack->getFittedState();
     Float_t P = fittedState.getMomMag();
     Float_t Px = fittedState.getMom().x();
     Float_t Py = fittedState.getMom().y();
     Float_t Pz = fittedState.getMom().z();
     h_chi2->Fill(chi2);
     h_Nmeasurements->Fill(fitStatus->getNdf());
     if (chi2 > chi2UL){ continue;}
     h_ppt->Fill(P,TMath::Sqrt(Px*Px+Py*Py));
     h_ppx->Fill(P,Px);
     auto pos = fittedState.getPos();
     h_xy->Fill(pos[0],pos[1]);
     h_pxpy->Fill(Px/Pz,Py/Pz);
// check for muon tag
     TVector3 posRPC; TVector3 momRPC;
     Double_t rc = MufluxReco::extrapolateToPlane(aTrack,cuts["zRPC1"], posRPC, momRPC);
     Float_t X =-1000;
     Float_t Y =-1000;
     Int_t Nx = RPCTrackX->GetEntries();
     Int_t Ny = RPCTrackY->GetEntries();
     if (Nx>0){
        RPCTrack *hit = (RPCTrack*)RPCTrackX->At(0);
        X = hit->m()*cuts["zRPC1"]+hit->b();}
     if (Ny>0){
        RPCTrack *hit = (RPCTrack*)RPCTrackY->At(0);
        Y = hit->m()*cuts["zRPC1"]+hit->b();}
      if (TMath::Abs(posRPC[0]-X)<5. && TMath::Abs(posRPC[1]-Y)<10.) { // within ~3sigma  X,Y from mutrack
        h_chi2mu->Fill(chi2);
        h_Nmeasurementsmu->Fill(fitStatus->getNdf());
        h_pptmu->Fill(P,TMath::Sqrt(Px*Px+Py*Py));
        h_ppxmu->Fill(P,Px);
        h_xymu->Fill(pos[0],pos[1]);
        h_pxpymu->Fill(Px/Pz,Py/Pz);
      }
//
     if (Ntracks==2 && k==0){
      genfit::Track* bTrack = (genfit::Track*)FitTracks->At(1);
      auto fitStatusb   = bTrack->getFitStatus();
      if (!fitStatusb->isFitConverged()){ continue;}
      Float_t chi2b = fitStatusb->getChi2()/float(fitStatusb->getNdf());
      if (chi2b > chi2UL){ continue;}
      auto fittedStateb = bTrack->getFittedState();
      Float_t Pb = fittedStateb.getMomMag();
      Float_t Pbx = fittedStateb.getMom().x();
      Float_t Pby = fittedStateb.getMom().y();
      Float_t Pbz = fittedStateb.getMom().z();
      h_p1p2->Fill(P,Pb);
      h_pt1pt2->Fill(TMath::Sqrt(Px*Px+Py*Py),TMath::Sqrt(Pbx*Pbx+Pby*Pby));
     }
   }
 }
}


void MufluxReco::sortHits(TClonesArray* hits, nestedList* l, Bool_t flag){
  nestedList spectrHitsSorted = *l;
 //spectrHitsSorted = {'_x':{1:{0:[],1:[],2:[],3:[]},2: {0:[],1:[],2:[],3:[]},3: {0:[],1:[],2:[],3:[]},4: {0:[],1:[],2:[],3:[]}},\
 //                    '_u':{1:{0:[],1:[],2:[],3:[]},2: {0:[],1:[],2:[],3:[]},3: {0:[],1:[],2:[],3:[]},4: {0:[],1:[],2:[],3:[]}},\
 //                    '_v':{1:{0:[],1:[],2:[],3:[]},2: {0:[],1:[],2:[],3:[]},3: {0:[],1:[],2:[],3:[]},4: {0:[],1:[],2:[],3:[]}}}
 for (Int_t k=0;k<hits->GetEntries();k++) {
   MufluxSpectrometerHit* hit = (MufluxSpectrometerHit*)hits->At(k);
   if ( !hit->isValid() && MCdata){continue;}
   if (flag && !MCdata){
    if (!hit->hasTimeOverThreshold() || !hit->hasDelay() || !hit->hasTrigger() ){ continue;} // no reliable TDC measuerement
  // remove noise hits
    auto result = std::find( noisyChannels.begin(), noisyChannels.end(), hit->GetDetectorID() );
    if ( result !=  noisyChannels.end()){ continue;}
   }
   if (hit->GetTimeOverThreshold() < cuts["tot"]) { continue;}
   auto info = hit->StationInfo();
   // info[0],info[1],info[2],info[3],viewDict[info[4]],info[5],info[6],info[7]
   // statnb,   vnb,    pnb,    lnb,    view,         channelID,tdcId,moduleId
   if (info[2] > 1 || info[3] >1){
    std::cout<< "sortHits: unphysical detctor ID "<<hit->GetDetectorID()<<std::endl;
    hit->Dump();
   }else{
    spectrHitsSorted[info[4]][info[0]][info[2]*2+info[3]].push_back(hit);
  }
 }
 *l = spectrHitsSorted;
}

void MufluxReco::fillHitMaps(Int_t nMax)
{
 Float_t deadThreshold = 1.E-4; // ~1% typical occupancy
 Int_t N = xSHiP->GetEntries(true);
 if (nMax<0){nMax=N;}
 xSHiP->Restart();
 gROOT->cd();
 std::cout<< "fillHitMaps: "<< N  <<std::endl;
 TTreeReaderArray <MufluxSpectrometerHit> Digi_MufluxSpectrometerHits(*xSHiP, "Digi_MufluxSpectrometerHits");
 TTreeReaderValue<FairEventHeader>* rvShipEventHeader= NULL;
 if (MCdata){ rvShipEventHeader = new TTreeReaderValue<FairEventHeader>(*xSHiP, "ShipEventHeader");}
 Int_t nx = 0;
 while (xSHiP->Next()){
  nx+=1;
  if (nx>nMax){break;}
  // std::cout<< "next event. #hits "<< Digi_MufluxSpectrometerHits.GetSize()  <<std::endl;
  for (Int_t k=0;k<Digi_MufluxSpectrometerHits.GetSize();k++) {
     MufluxSpectrometerHit* hit = &(Digi_MufluxSpectrometerHits[k]);
     auto info = hit->StationInfo();
     Int_t s=info[0]; Int_t v=info[1]; Int_t p=info[2]; Int_t l=info[3]; Int_t channelNr=info[5]; 
     Int_t tdcId=info[6]; Int_t moduleId=info[7];
     TString view = "_x"; 
     if (info[4]==1){view="_u";}
     if (info[4]==2){view="_v";}
     TString tot = "";
     if (!hit->hasTimeOverThreshold()){ tot = "_noToT";}
     TString histo; histo.Form("%d",1000*s+100*p+10*l);histo+=view;
     TH1D* h = (TH1D*)(gDirectory->GetList()->FindObject(histo));
     if (!h){
       std::cout<< "fillHitMaps: ERROR histo not known "<< histo  <<std::endl;
       continue;
     }
     h->Fill(channelNr);
     auto check = std::find( noisyChannels.begin(), noisyChannels.end(), hit->GetDetectorID() );
     if (check!=noisyChannels.end()){ continue;}
     Float_t t0 = 0;
     if (MCdata){ rvShipEventHeader->Get()->GetEventTime(); }
     TString TDChisto = "TDC"; TDChisto+=histo;TDChisto+=moduleId+tot;
     h = (TH1D*)(gDirectory->GetList()->FindObject(TDChisto));
     if (!h){
       std::cout<< "fillHitMaps: ERROR histo not known "<< TDChisto  <<std::endl; 
       continue;
     }
     h->Fill( hit->GetDigi()-t0);
     TString channel = "TDC";channel+=hit->GetDetectorID();
     h = (TH1D*)(gDirectory->GetList()->FindObject(channel));
     h->Fill( hit->GetDigi()-t0);
   }
  }
}
// -----   Destructor   ----------------------------------------------------
MufluxReco::~MufluxReco() { }
// -------------------------------------------------------------------------

 
ClassImp(MufluxReco) 
