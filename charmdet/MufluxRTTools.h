#ifndef MUFLUXRTTOOLS_H
#define MUFLUXRTTOOLS_H

#include <vector>
#include <iostream>
#include "TVector3.h"
#include "MufluxSpectrometer.h"
#include "MufluxSpectrometerHit.h"
#include "Track.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"

void helloMorten();

class MufluxRTTools {
	public:
	Double_t calculate2D(int det_id, TVector3 vtop, TVector3 vbot, TVector3 mom, TVector3 pos);
	TH1D* detRT(TH2D timedist, int pm);
	int chnumgasorder(int det_id);

	private:
	Double_t distance;
	TH1D* hrtrel;
	TH1D* hrtrelsigma;
	int chnum;

};

#endif
