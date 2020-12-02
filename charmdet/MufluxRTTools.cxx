#include "MufluxRTTools.h"
#include "TVector3.h"
#include "Track.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerHit.h"
#include <vector>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"



void helloMorten(){
	std::cout<<"Hello"<<std::endl;
}

Double_t MufluxRTTools::calculate2D(int det_id, TVector3 vtop, TVector3 vbot, TVector3 mom, TVector3 pos)
{

	Double_t dist = 0;
	TVector3 cross;
	TVector3 a = pos-vtop;
	cross = mom.Cross(vbot-vtop);
	dist = abs(a.Dot(mom)/cross.Mag());
	if (a.Dot(mom)/cross.Mag()<0){
		std::cout<<"UNTER NULL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	}


	return dist;
}

TH1D* MufluxRTTools::detRT(TH2D timedist, int pm)
/*
 version 1.3 change: 0.1 acceptance now is 0.4
 */
{
	int plusminusbin = pm;
	TH1D *hrtrel = new TH1D("hrtrel", "RT-Relation",30000,-490,2450);
	TH1D *hrtrelsigma = new TH1D("hrtrelsigma", "Sigma of RT-Relation",30000,-490,2450);

	std::cout<<"TH2D entries: "<< timedist.GetEntries()<< std::endl;
	std::cout<<"TH2D first above 0: "<< timedist.FindFirstBinAbove(0)<<" last above 0: "<< timedist.FindLastBinAbove(0) << std::endl;
	for(int i=timedist.FindFirstBinAbove(0);i<timedist.FindLastBinAbove(0);i++){
		std::cout<<"Bin Number: "<< i<<std::endl;
		std::cout<<"From Bin number: "<<timedist.FindFirstBinAbove(0)<< " to bin number: "<< timedist.FindLastBinAbove(0)<<std::endl;
		int plusminusbinhier = plusminusbin;
		TH1D *projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
		while(projection->GetEntries() < 100 && plusminusbinhier <= 100){
			plusminusbinhier++;
	        projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
	    }
		std::cout<<"Number of Entries: "<< projection->GetEntries()<<std::endl;
		if(projection->GetEntries()>20){
			std::cout<<"More than 20 Entries"<<std::endl;
			TF1 *projectionfit = new TF1("fit","gaus");
			std::cout<<"Slide declared"<<std::endl;
	        projection->Fit(projectionfit,"Q");
	        std::cout<<"projection calculated"<<std::endl;
            float mean =projectionfit->GetParameter(1);
            float sigma =projectionfit->GetParameter(2);
            std::cout<<"Radius : " << mean<< " Sigma: "<< sigma <<" Eintraege: "<< projection->GetEntries()<< " fit: "<<i<<std::endl;
            if(mean>0){ // && mean <18.15){
            	hrtrelsigma->SetBinContent(i,sigma);
	            hrtrel->SetBinContent(i,mean);
	            std::cout<<"Normal mean used; Mean: "<< mean<<std::endl;
 /*
               	if(abs(rtinit.GetBinContent(i)-mean)>0.4){
            		hrtrel->SetBinContent(i,hrtrel->GetBinContent(i-1));
            		hrtrelsigma->SetBinContent(i,hrtrelsigma->GetBinContent(i-1));
            		if(mean<1.815){std::cout<<"Mean: "<< mean<< " Eintraege: "<< i <<" prev value:"<< hrtrel->GetBinContent(i-1)<<std::endl;}
            	} else {
            	hrtrelsigma->SetBinContent(i,sigma);
	            hrtrel->SetBinContent(i,mean);
	            std::cout<<"Normal mean used; Mean: "<< mean<<std::endl;
            }

  */
            }
		}
}
	std::cout<<"Plusminus von "<<pm<<std::endl;
	return hrtrel;

}
/* version 1.1
{
	int plusminusbin = pm;
	TH1D *hrtrel = new TH1D("hrtrel", "RT-Relation",55000,-490,4900);
	TH1D *hrtrelsigma = new TH1D("hrtrelsigma", "Sigma of RT-Relation",55000,-490,4900);

	std::cout<<"TH2D entries: "<< timedist.GetEntries()<< std::endl;

	for(int i=timedist.FindFirstBinAbove(0);i<timedist.FindLastBinAbove(0);i++){
		//std::cout<<"Bin Number: "<< i<<std::endl;
		int plusminusbinhier = plusminusbin;
		TH1D *projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
		while(projection->GetEntries() < 100 && plusminusbinhier <= 100){
			plusminusbinhier++;
	        projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
	    }
		//std::cout<<"Number of Entries: "<< projection->GetEntries()<<std::endl;
		if(projection->GetEntries()>20){
			TF1 *projectionfit = new TF1("fit","gaus");
	        projection->Fit(projectionfit,"Q");
            float mean =projectionfit->GetParameter(1);
            float sigma =projectionfit->GetParameter(2);
            std::cout<<"Radius : " << mean<< " Sigma: "<< sigma <<" Eintraege: "<< projection->GetEntries()<< " fit: "<<i<<std::endl;
            if(mean>0){ // && mean <18.15){
            	if(hrtrel->GetBinContent(i-1)>mean || (mean> (hrtrel->GetBinContent(i-1)+0.1 ))){
            		hrtrel->SetBinContent(i,hrtrel->GetBinContent(i-1));
            		hrtrelsigma->SetBinContent(i,hrtrelsigma->GetBinContent(i-1));
            		if(mean<1.815){std::cout<<"Mean: "<< mean<< " Eintraege: "<< i <<" prev value:"<< hrtrel->GetBinContent(i-1)<<std::endl;}
            	} else if(((hrtrel->GetBinContent(i-1)!=hrtrel->GetBinContent(i-2)) && (mean > hrtrel->GetBinContent(i-1)*2 - hrtrel->GetBinContent(i-2))) && (hrtrel->GetEntries()>10)){
            		hrtrel->SetBinContent(i,hrtrel->GetBinContent(i-1));
            		hrtrelsigma->SetBinContent(i,hrtrelsigma->GetBinContent(i-1));
            		if(mean<1.815){std::cout<<"Mean: "<< mean<< " Eintraege: "<< i <<" prev value:"<< hrtrel->GetBinContent(i-1)<<std::endl;}
            	} else {
            	hrtrelsigma->SetBinContent(i,sigma);
	            hrtrel->SetBinContent(i,mean);
	            std::cout<<"ELSE!"<<std::endl;
            }
            }
		}
}
	std::cout<<"Plusminus von "<<pm<<std::endl;
	return hrtrel;
}
*/
/* Version 1.0
{
	int plusminusbin = pm;
	TH1D *hrtrel = new TH1D("hrtrel", "RT-Relation",55000,-490,4900);
	TH1D *hrtrelsigma = new TH1D("hrtrelsigma", "Sigma of RT-Relation",55000,-490,4900);

	for(int i=0;i<timedist.GetNbinsX();i++){
		int plusminusbinhier = plusminusbin;
		TH1D *projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
		while(projection->GetEntries() < 100 && plusminusbinhier <= 100){
			plusminusbinhier++;
	        projection = timedist.ProjectionY("nbin",TMath::Max(1,i-plusminusbinhier),TMath::Min(i+plusminusbinhier,timedist.GetNbinsX()-1));
	    }
		if(projection->GetEntries()>20){
			TF1 *projectionfit = new TF1("fit","gaus");
	        projection->Fit(projectionfit,"Q");
            float mean =projectionfit->GetParameter(1);
            float sigma =projectionfit->GetParameter(2);
            //std::cout<<"Radius : " << mean<< " Sigma: "<< sigma <<" Eintraege: "<< projection->GetEntries()<< " fit: "<<i<<std::endl;
            if(mean>0){ // && mean <18.15){
            	if(hrtrel->GetBinContent(i-1)>mean || mean> (hrtrel->GetBinContent(i-1)+0.1)){
            		hrtrel->SetBinContent(i,hrtrel->GetBinContent(i-1));
            		hrtrelsigma->SetBinContent(i,hrtrelsigma->GetBinContent(i-1));
            		if(mean<1.815){std::cout<<"Mean: "<< mean<< " Eintraege: "<< i <<" prev value:"<< hrtrel->GetBinContent(i-1)<<std::endl;}
            	} else{
            	hrtrelsigma->SetBinContent(i,sigma);
	            hrtrel->SetBinContent(i,mean);
	            std::cout<<"ELSE!"<<std::endl;
            }
            }
		}
}
	std::cout<<"Plusminus von "<<pm<<std::endl;
	return hrtrel;
}
*/
int MufluxRTTools::chnumgasorder(int det_id)
{
	if ((det_id/10000000) < 3 && (det_id > 0)){
			chnum= (det_id/10000000-1)*96 + ((det_id/1000000)%10) * 48 + ((det_id/100000)%10)*2 + ((det_id/10000)%10) + (det_id%100)*4;
	}else if(det_id/10000000 == 4){
		chnum =192 + ((det_id/100000)%10)*96 + ((det_id/10000)%10)*48 + det_id%100/12 + det_id%100%12*4;
	}else if(det_id/10000000 == 3){
		chnum =192*2 + (1-(det_id/100000)%10)*96 + (1-(det_id/10000)%10)*48 + det_id%100/12 + det_id%100%12*4;
	}else{
		chnum=-1;
	}


	return chnum;
}

