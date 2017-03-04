// -*- C++ -*-
//
// Package:    HCALCommissioning2017/HESourcing
// Class:      HESourcing
// 
/**\class HESourcing HESourcing.cc HCALCommissioning2017/HESourcing/plugins/HESourcing.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Burak Bilki
//         Created:  Fri, 03 Feb 2017 19:32:36 GMT
//
//


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "EventFilter/HcalRawToDigi/interface/AMC13Header.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUHTRData.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/// Per Event Header Structure
struct eventHeader
{
	uint32_t cdf0;
	uint32_t cdf1;
	uint32_t cdf2;
	uint32_t cdf3;
	uint32_t h0;
	uint32_t h1;
	uint32_t h2;
	uint32_t h3;
};

/// Structure for Source Position Data
struct xdaqSourcePositionDataFormat
{
	uint32_t cdfHeader[4];
	uint16_t n_doubles;
	uint16_t n_strings;
	uint16_t key_length;
	uint16_t string_value_length;
	char     start_of_data; // see below
};

/// Structure for Sourcing histos
struct sourcingHistos
{
	int nevt;
	int reel;
	string name;
	TH1D* h;
	TH2F* hh;
	TH2F* hhp;
};
sourcingHistos sh;

/// Structure for Calib
struct calibPar
{
	int ch;
	double ped;
	double spe;
};
calibPar cp;

/// Structure for ReelPos
struct reelPos
{
	int ch;
	int b1;
	int b2;
	int bm;
	double rpm;
};
reelPos rps;

struct edata
{
	Int_t reel;
	string name;
// 	vector <TH1D*> adc;
	vector <vector <Int_t>> *adc;
	vector <string> *cname;
// 	vector <Int_t> adc[64];
// 	Int_t adc[576][64];
// 	string cname[576];
};
edata ed;

class HESourcing : public edm::one::EDAnalyzer<edm::one::SharedResources> 
{
	public:
		explicit HESourcing(const edm::ParameterSet&);
		~HESourcing();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		TFile *_file;
		TTree *_tree;
		int histoFED;
		int driverFED;
		int EID;
		int numChannels;
		string outFileName;
		int runType;
		string digiCollection;
		
		edm::EDGetTokenT<FEDRawDataCollection> raw_token;
		std::string electronicsMapLabel_;
};
HESourcing::HESourcing(const edm::ParameterSet& iConfig)
{
	raw_token = consumes<FEDRawDataCollection>(edm::InputTag("source"));
	runType = iConfig.getParameter<int>("RunType");
	outFileName=iConfig.getUntrackedParameter<string>("OutFileName");
	histoFED = iConfig.getParameter<int>("histoFED");
	driverFED = iConfig.getParameter<int>("driverFED");
	_file = new TFile(outFileName.c_str(), "recreate");
	_tree = new TTree("Events", "Events");
	_tree->Branch("reel", &ed.reel);
	_tree->Branch("name", &ed.name);
	_tree->Branch("adc", &ed.adc);
// 	_tree->Branch("cname", &ed.cname);
// 	_tree->Branch("adc", &ed.adc,"adc[576][64]/I");
	_tree->Branch("cname", &ed.cname);
	EID=0;
}

HESourcing::~HESourcing()
{
	_file->cd();
	_file->Write();
	_file->Close();
}

void HESourcing::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	edm::Handle<FEDRawDataCollection> rawraw;  
	edm::ESHandle<HcalElectronicsMap> item;
	edm::ESHandle<HcalDbService> pSetup;

	iEvent.getByToken(raw_token, rawraw);          
	iSetup.get<HcalDbRecord>().get(pSetup);
	iSetup.get<HcalElectronicsMapRcd>().get(item);

// 	const HcalElectronicsMap* readoutMap = item.product();
	if(EID==0)
	{
		int fedid=0;
		for(int i1=50;i1<=90;i1++)
		{
			if((rawraw->FEDData(i1)).size()>0){fedid=i1;break;}
		}
		histoFED=fedid;
	}
	const FEDRawData& raw = rawraw->FEDData(histoFED);
	
	char hname[500];
	
	if(runType==2)
	{
		//the source position
		
		const FEDRawData& spraw = rawraw->FEDData(driverFED);
		
		const struct xdaqSourcePositionDataFormat* sp = (const struct xdaqSourcePositionDataFormat*)(spraw.data());
		map<string,double> sp_dblmap;
		map<string,string> sp_strmap;

		// List of doubles:
		const char   *keyptr = &sp->start_of_data;
		const double *valptr = (const double *)(&sp->start_of_data + sp->n_doubles*sp->key_length);
		for (int i=0; i<sp->n_doubles; i++)
		{
	// 		cout<<i<<" " << keyptr << " = " << *valptr << endl;
			if(i==4) ed.reel=*valptr;
			sp_dblmap[keyptr] = *valptr;
			keyptr += sp->key_length;
			valptr++;
		}
		// List of strings:
		keyptr = (const char *)valptr;
		const char *strptr = (keyptr + sp->n_strings*sp->key_length);
		for (int i=0; i<sp->n_strings; i++)
		{
	// 		cout << i<< " "<<keyptr << " = " << strptr << endl;
			if(i==1) ed.name=strptr;
			sp_strmap[keyptr] = strptr;
			keyptr += sp->key_length;
			strptr += sp->string_value_length;
		}
	}
	else{ed.reel=0;ed.name="";}
	
	//the histos
	
// 	cout<<"Event: "<<Nevt<<endl;
	const struct eventHeader* eh =(const struct eventHeader*)(raw.data());
	const uint32_t* pData = (const uint32_t*) raw.data(); 
	//Read event header
	uint32_t numHistos  = ((eh->h3)>>16)&0xFFFF;
// 	cout << "Number of Histograms: " << numHistos << endl;
	uint32_t numBins    = ((eh->h3)>>1)&0x0000FFFE; //includes overflow and header word
// 	cout << "Bins per Histogram: " << numBins << endl;
// 	bool sepCapIds = eh->h3&0x00000001;
// 	cout << "Separate CapIds: " << sepCapIds << endl;
	
	//Set histogram word pointer to first histogram    
	uint32_t crate   = 0;
	uint32_t slot    = 0;
	uint32_t fiber   = 0;
	uint32_t channel = 0;
// 	uint32_t cap     = 0;
	//Loop over data
	pData+=8;
// 	vector <TH1D*> TD;TH1D* td;
	vector <int> td;
// 	Int_t td[64]={0};
	for (unsigned int iHist = 0; iHist<numHistos; iHist++)
	{
// 		cout << "Histogram " << iHist <<" header: "<< *pData << endl;
		crate   = ((*pData)>>16)&0x00FF;
// 		cout << "Crate: " << crate << endl;
		slot    = ((*pData)>>12)&0x0000F;
// 		cout << "Slot: " << slot << endl;
		fiber   = (*pData>>7)&0x1F;
// 		cout << "Fiber: " << fiber << " "<< (*pData&0x00000F80)<<endl;
		channel = (*pData>>2)&0x1F;
		
		sprintf(hname,"ADC_%d_%d_%d_%d",crate,slot,fiber,channel);
// 		string tn="ADC_"+crate+"_"+slot+"_"+fiber+"_"+channel;
		string tn(hname);
// 		ed.cname[iHist]=tn;
// 		cout<<tn<<endl;
		ed.cname->push_back(tn);
// 		td=new TH1D(hname,hname,256,-0.5,255.5);
// 		TD.push_back(td);
		
// 		cout << "Det Id: " << ((HcalDetId)did) << endl;
// 		HcalUHTRhistogramDigiMutable digi = histoDigiCollection->addHistogram( (HcalDetId)did );
// 		double bs=0.;double bn=0.;
// 		double bc=0.;
		
// 		cout<<tn<<" ";
// 		
		for(unsigned int iBin = 0; iBin<numBins+1; iBin++)
		{
// 			TD[iHist]->Fill(iBin,pData[iBin+1]);
			td.push_back(pData[iBin+1]);
// 			ed.adc[iHist][iBin]=pData[iBin+1];
// 			cout<<pData[iBin+1]<<" ";
// 			td[iBin]=pData[iBin+1];
		}
// 		cout<<endl;
		
		if(iHist<(numHistos-1)) 
		{
			pData+=(numBins+2);
		}
		
// 		ed.adc.push_back(TD[iHist]);
		ed.adc->push_back(td);
		td.clear();
// 		ed.adc.push_back(td);
	}
	_tree->Fill();
	ed.adc->clear();
	ed.cname->clear();
// 	TD.clear();
	EID++;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HESourcing::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HESourcing::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HESourcing::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HESourcing);
