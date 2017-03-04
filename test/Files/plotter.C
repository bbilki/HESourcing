#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TClonesArray.h"

#include <sys/types.h>
#include <dirent.h>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TH2Poly.h>
#include <TPaveLabel.h>
#include <TPaveText.h>

using namespace std;using namespace ROOT::Math;

struct edata
{
	string *name;
	Int_t reel;
	vector <vector <Int_t>> *adc;
	vector <string> *cname;
// 	Int_t adc[576][64];
// 	string cname[576];
};
// edata ed;

struct adcinfo
{
	int mdepth;
	int box;
	int tileind;
	int side;
	double base;
	int bsize;
	vector <int> bs[2];
	vector <int> rr[2];
	vector <double> qm;
};
adcinfo ai;
vector <adcinfo> AI;

struct mapdata
{
	int crate;
	int slot;
	int fiber;
	int channel;
	int ieta;
	int iphi;
	int depth;
	int box;
	int tileind;
	string tilename;
};
mapdata mp;
vector <mapdata> MP;

struct histos
{
	string name;
	int layer;
	int type;
	int color;
	int mapind;
	int badcap;
	TH1D* ADCperReel[3];
	TH1D* ADCSpec;
	int startreel;
	int endreel;
	int firstevent1;
	int firstevent2;
	double adcmax1;
	int eadcmax1;
	double adcmax2;
	int eadcmax2;
	int radcmax1;
	int radcmax2;
	int nHgtTh;
};
histos H;
vector <histos> HH;

struct NTupledata
{
	int reel;
	int layer;
	int type;//0 short 1 long
	string *name;//wire name
	
	vector <int> *mapind;
	vector <int> *badcap;
	vector <vector <float>> *meanadc;
	vector <vector <float>> *meanadcerr;
	vector <vector <float>> *base;
	vector <vector <float>> *baseerr;
	vector <vector <int>> *N;
};
// NTupledata ND;
NTupledata ed;

struct summary
{
	int RunNo;
	int crate;
	int slot;
	int fiber;
	int channel;
	int ieta;
	int iphi;
	int depth;
	int layer;
	int SL;//0 short; 1 long
	float meanQ;
	float meanQerr;
	float sigmaQ;
	float sigmaQerr;
};
summary ss;vector <summary> SS;

int RunNo=0;
bool isPhaseIHEP17=false;

int readmap()
{
	ifstream infile;
	if(RunNo<287582){infile.open("semap_v1.txt");}
	else{infile.open("semap_v2.txt");}
	while(!infile.eof())
	{
		infile>>mp.crate>>mp.slot>>mp.fiber>>mp.channel>>mp.ieta>>mp.iphi>>mp.depth>>mp.box>>mp.tileind>>mp.tilename;
// 		cout<<mp.crate<<" "<<mp.slot<<" "<<mp.fiber<<" "<<mp.channel<<" "<<mp.ieta<<" "<<mp.iphi<<" "<<mp.depth<<" "<<mp.box<<" "<<mp.tileind<<" "<<mp.tilename<<endl;
		MP.push_back(mp);
	}
	infile.close();
}

vector<int> findLocation(string sample, char findIt)
{
    vector<int> characterLocations;
    for(int i =0; i < sample.size(); i++)
        if(sample[i] == findIt)
            characterLocations.push_back(i);

    return characterLocations;
}

Double_t fitf(Double_t *x,Double_t *par)
{
	Double_t fitval=0.;
	Double_t sigma=0.;
	Double_t tt=(x[0]-par[1]);
	if(x[0]<par[1])
	{
		sigma=par[2]+par[3]*pow(par[1]-x[0],par[4]);
	}
	else
	{
		sigma=par[2]+par[5]*pow(x[0]-par[1],par[6]);
	}
	fitval=par[0]*exp(-0.5*pow(tt/sigma,2));
	return fitval;
}

int analyze()
{
	char hname[900];
	char hname2[500];
	sprintf(hname,"Results_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	sprintf(hname,"../NTuples2/N2_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("reel",&ed.reel);
	tree->SetBranchAddress("layer",&ed.layer);
	tree->SetBranchAddress("type",&ed.type);
	tree->SetBranchAddress("name",&ed.name);
	tree->SetBranchAddress("mapind",&ed.mapind);
	tree->SetBranchAddress("badcap",&ed.badcap);
	tree->SetBranchAddress("meanadc",&ed.meanadc);
	tree->SetBranchAddress("meanadcerr",&ed.meanadcerr);
	tree->SetBranchAddress("base",&ed.base);
	tree->SetBranchAddress("baseerr",&ed.baseerr);
	tree->SetBranchAddress("N",&ed.N);
	
	TH1F* NormChi2=new TH1F("NormChi2","NormChi2",200,0.,10.);
	TH1F* NHgtTh=new TH1F("NHgtTh","NHgtTh",500,0.,500.);
	TF1* tf2=new TF1("tf2","[0]",2500.,8500.);
	string stype[2]={"S","L"};
	
	int hind=0;string hname1="";
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<10;i++)
	{
		tree->GetEntry(i);
// 		cout<<i<<" "<<ed.reel<<" "<<ed.layer<<" "<<ed.type<<" "<<*ed.name<<" "<<ed.mapind->size()<<" "<<ed.mapind->at(0)<<" "<<ed.meanadc->at(0)<<" "<<ed.meanadcerr->at(0)<<" "<<ed.base->at(0)<<" "<<ed.baseerr->at(0)<<endl;
		for(int i1=0;i1<ed.mapind->size();i1++)
		{
			sprintf(hname,"%d %d %d %d %s",MP[ed.mapind->at(i1)].ieta,MP[ed.mapind->at(i1)].iphi,MP[ed.mapind->at(i1)].depth,ed.layer,stype[ed.type].c_str());
			string hname1(hname);
			isPhaseIHEP17=false;if(RunNo>=287582 && MP[H.mapind].iphi>=63 && MP[H.mapind].iphi<=66 && MP[H.mapind].ieta>0) isPhaseIHEP17=true;
			hind=-1;
			for(int i2=0;i2<HH.size();i2++)
			{
				if(HH[i2].name==hname1){hind=i2;break;}
			}
			if(hind==-1)
			{
				H.name=hname1;
				H.layer=ed.layer;
				H.type=ed.type;
				H.color=1;
				H.startreel=0;
				H.endreel=0;
				H.adcmax1=0;
				H.eadcmax1=0;
				H.adcmax2=0;
				H.eadcmax2=0;
				H.radcmax1=0;
				H.radcmax2=0;
				H.firstevent1=0;
				H.firstevent2=0;
				H.mapind=ed.mapind->at(i1);
				H.badcap=ed.badcap->at(i1);
// 				sprintf(hname,"QperReel %s",hname1.c_str());
// 				sprintf(hname,"QperReel iPhi %d Layer %d %s", MP[ed.mapind->at(i1)].iphi, ed.layer, stype[ed.type].c_str());
				sprintf(hname,"%d %d %d %d %s", MP[ed.mapind->at(i1)].ieta,MP[ed.mapind->at(i1)].iphi,MP[ed.mapind->at(i1)].depth,ed.layer,stype[ed.type].c_str());
				H.ADCperReel[0]=new TH1D(hname,hname,300,2500.5,8500.5);
// 				H.ADCperReel[0]=new TH1D(hname,hname,120,2500.5,8500.5);
				H.ADCperReel[0]->Sumw2();H.ADCperReel[0]->SetBit(TH1::kIsAverage);
				H.ADCperReel[0]->GetXaxis()->SetTitle("Wire Position (mm)");H.ADCperReel[0]->GetXaxis()->CenterTitle();
				H.ADCperReel[0]->GetYaxis()->SetTitle("Mean Charge / 20 mm (fC)");H.ADCperReel[0]->GetYaxis()->CenterTitle();
				sprintf(hname,"QperReel %s norm",hname1.c_str());
				H.ADCperReel[1]=new TH1D(hname,hname,300,2500.5,8500.5);
// 				H.ADCperReel[1]=new TH1D(hname,hname,120,2500.5,8500.5);
				sprintf(hname,"QperReel %s N",hname1.c_str());
				H.ADCperReel[2]=new TH1D(hname,hname,300,2500.5,8500.5);
// 				H.ADCperReel[2]=new TH1D(hname,hname,120,2500.5,8500.5);
				sprintf(hname,"QSpectrum %s",hname1.c_str());
				
				if(isPhaseIHEP17)
				{
					H.ADCSpec=new TH1D(hname,hname,100,0.,20.);
					H.ADCSpec->GetXaxis()->SetTitle("Mean Charge (fC)");H.ADCSpec->GetXaxis()->CenterTitle();
					H.ADCSpec->GetYaxis()->SetTitle("Events / 0.2 fC");H.ADCSpec->GetYaxis()->CenterTitle();
				}
				else
				{
					H.ADCSpec=new TH1D(hname,hname,250,-0.05,0.2);
					H.ADCSpec->GetXaxis()->SetTitle("Mean Charge (fC)");H.ADCSpec->GetXaxis()->CenterTitle();
					H.ADCSpec->GetYaxis()->SetTitle("Events / 10^{-5} fC");H.ADCSpec->GetYaxis()->CenterTitle();
				}
				H.nHgtTh=0;
				HH.push_back(H);
				hind=HH.size()-1;
			}
			for(int i4=0;i4<4;i4++)
			{
				if(i4!=HH[hind].badcap)
				{
					HH[hind].ADCperReel[0]->Fill(ed.reel,ed.meanadc->at(i1)[i4]);
					HH[hind].ADCperReel[1]->Fill(ed.reel);
// 					HH[hind].ADCperReel[2]->Fill(ed.reel,ed.N->at(i1)[i4]);
					HH[hind].ADCperReel[2]->Fill(ed.reel,pow(ed.meanadc->at(i1)[i4],2.));
					if(ed.meanadc->at(i1)[i4]>0.015 && !isPhaseIHEP17) {HH[hind].nHgtTh++;}
					if(ed.meanadc->at(i1)[i4]>1. && isPhaseIHEP17) {HH[hind].nHgtTh++;}
				}
			}
		}
		if(i%1000==0) cout<<"Event :"<<i<<" / "<<tree->GetEntries()<<endl;
	}
	outroot->cd();
	double nchi2=0.;
	vector <float> errs;
	for(int i2=0;i2<HH.size();i2++)
	{
// 		HH[i2].ADCperReel[0]->Write();
// 		HH[i2].ADCperReel[1]->Write();
		errs.clear();
		for(int i3=1;i3<=HH[i2].ADCperReel[0]->GetNbinsX();i3++)
		{
			errs.push_back(HH[i2].ADCperReel[0]->GetBinError(i3));
		}
		HH[i2].ADCperReel[0]->Divide(HH[i2].ADCperReel[1]);
		HH[i2].ADCperReel[2]->Divide(HH[i2].ADCperReel[1]);
		for(int i3=1;i3<=HH[i2].ADCperReel[0]->GetNbinsX();i3++)
		{
// 			HH[i2].ADCperReel[0]->SetBinError(i3,HH[i2].ADCperReel[0]->GetBinError(i3)/sqrt(HH[i2].ADCperReel[1]->GetBinContent(i3)));
// 			HH[i2].ADCperReel[0]->SetBinError(i3,HH[i2].ADCperReel[0]->GetBinError(i3)/sqrt(HH[i2].ADCperReel[2]->GetBinContent(i3)));
			
			
			
// 			if(HH[i2].ADCperReel[2]->GetBinContent(i3)>0)
// 			{
// 				HH[i2].ADCperReel[0]->SetBinError(i3,sqrt(errs[i3-1]/HH[i2].ADCperReel[2]->GetBinContent(i3)));
// 			}
// 			else
// 			{
// 				HH[i2].ADCperReel[0]->SetBinError(i3,0.);
// 			}
			
			
			if(HH[i2].ADCperReel[1]->GetBinContent(i3)>0)
			{
				HH[i2].ADCperReel[0]->SetBinError(i3,sqrt((HH[i2].ADCperReel[2]->GetBinContent(i3)-pow(HH[i2].ADCperReel[0]->GetBinContent(i3),2.))/HH[i2].ADCperReel[1]->GetBinContent(i3)));
			}
			else
			{
				HH[i2].ADCperReel[0]->SetBinError(i3,0.);
			}
		}
		
		int nPgtTh=0;
		for(int i3=1;i3<=HH[i2].ADCperReel[0]->GetNbinsX();i3++)
		{
			if(HH[i2].ADCperReel[0]->GetBinContent(i3)>0.015 && !isPhaseIHEP17) nPgtTh++;
			if(HH[i2].ADCperReel[0]->GetBinContent(i3)>1. && isPhaseIHEP17) nPgtTh++;
		}
		NHgtTh->Fill(HH[i2].nHgtTh);
		if(nPgtTh<3)
		{
			HH[i2]=HH[HH.size()-1];
			HH.pop_back();
			i2--;
			continue;
		}
// 		HH[i2].ADCperReel[0]->Fit(tf2,"q","q",2500.,8500.);
// 		nchi2=tf2->GetChisquare()/tf2->GetNDF();
// 		NormChi2->Fill(nchi2);
// 		sprintf(hname,"%s %f",HH[i2].ADCperReel[0]->GetName(),nchi2);
// 		HH[i2].ADCperReel[0]->SetTitle(hname);
// 		HH[i2].ADCperReel[0]->Write();
	}
// 	int colors[20]={1,2,3,4,5,6,7,8,9,21,11,12,13,14,15,16,17,18,19,20};
	int colors[22]={1,2,3,4,6,7,8,9,12,22,32,42,30,40,41,43,44,45,46,47,38,49};
	vector <int> processedLayers;
	TCanvas* cc2=new TCanvas("cc2","cc2",600,600);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend* Legend;TPaveLabel* pl;TPaveText *pt;
	for(int i2=0;i2<HH.size();i2++)
	{
		bool found=false;
		for(int i1=0;i1<processedLayers.size();i1++)
		{
			if(processedLayers[i1]==HH[i2].layer){found=true;break;}
		}
		if(found) continue;
		
		HH[i2].ADCperReel[0]->Draw();
		HH[i2].ADCperReel[0]->GetYaxis()->SetTitleOffset(1.5);
		HH[i2].ADCperReel[0]->GetYaxis()->SetLabelSize(0.03);
		cc2->Update();
// 		sprintf(hname,"QperReel iPhi %d Layer %d %s", MP[HH[i2].mapind].iphi, HH[i2].layer, stype[HH[i2].type].c_str());
// 		HH[i2].ADCperReel[0]->SetTitle(hname);
		double lprofymax=HH[i2].ADCperReel[0]->GetBinContent(HH[i2].ADCperReel[0]->GetMaximumBin());
		int ic=1;
		int fb=HH[i2].ADCperReel[0]->FindFirstBinAbove(0.01);
		int lb=HH[i2].ADCperReel[0]->FindLastBinAbove(0.01);
		for(int i3=i2+1;i3<HH.size();i3++)
		{
			if(HH[i2].layer==HH[i3].layer)
			{
				HH[i3].color=colors[ic];
// 				cout<<HH[i3].name<<" "<<ic<<endl;
				HH[i3].ADCperReel[0]->SetLineColor(HH[i3].color);
				HH[i3].ADCperReel[0]->Draw("same");
				ic++;
				if(HH[i3].ADCperReel[0]->GetBinContent(HH[i3].ADCperReel[0]->GetMaximumBin())>lprofymax) lprofymax=HH[i3].ADCperReel[0]->GetBinContent(HH[i3].ADCperReel[0]->GetMaximumBin());
				
				if(!isPhaseIHEP17)
				{
					if(HH[i3].ADCperReel[0]->FindFirstBinAbove(0.01)<fb) fb=HH[i3].ADCperReel[0]->FindFirstBinAbove(0.01);
					if(HH[i3].ADCperReel[0]->FindLastBinAbove(0.01)>lb) lb=HH[i3].ADCperReel[0]->FindLastBinAbove(0.01);
				}
				else
				{
					if(HH[i3].ADCperReel[0]->FindFirstBinAbove(5.)<fb) fb=HH[i3].ADCperReel[0]->FindFirstBinAbove(5.);
					if(HH[i3].ADCperReel[0]->FindLastBinAbove(5.)>lb) lb=HH[i3].ADCperReel[0]->FindLastBinAbove(5.);
				}
// 				HH[i3].ADCperReel[0]->Write();
			}
		}
		lprofymax+=0.01;
		if(isPhaseIHEP17){lprofymax+=1.5;}
		HH[i2].ADCperReel[0]->GetXaxis()->SetRangeUser(HH[i2].ADCperReel[0]->GetBinCenter(fb-5),HH[i2].ADCperReel[0]->GetBinCenter(lb+5));
		HH[i2].ADCperReel[0]->GetYaxis()->SetRangeUser(0.005,lprofymax);
		processedLayers.push_back(HH[i2].layer);
		Legend=cc2->BuildLegend();
		sprintf(hname,"QperReel iPhi %d Layer %d %s", MP[HH[i2].mapind].iphi, HH[i2].layer, stype[HH[i2].type].c_str());
// 		cout<<hname<<endl;
// 		pl = new TPaveLabel(0.1, 0.9, 0.9, 0.9,hname);
// 		pl->SetFillColor(16);
// 		pl->SetTextFont(52);
// 		pl->Draw();
		
		pt = new TPaveText(0.15, 0.9, 0.85, 1.,"NB NDC");
		pt->AddText(hname);
		pt->SetLineColor(0);
		pt->SetFillColor(0);
		pt->SetTextFont(42);
		pt->Draw();
		
		cc2->Update();
		sprintf(hname,"ADCperReel_%d.png",HH[i2].layer);
		cc2->SaveAs(hname);
// 		HH[i2].ADCperReel[0]->Write();
	}
	
	outroot->cd();
	vector <int> maxreelset;
	int startbin=0;int endbin=0;int maxbin=0;
	for(int i2=0;i2<HH.size();i2++)
	{
		maxbin=HH[i2].ADCperReel[0]->GetMaximumBin();
		maxreelset.push_back(((int)HH[i2].ADCperReel[0]->GetBinCenter(maxbin)));
// 		for(int i3=maxbin;i3>=1;i3--)
// 		{
// 			if(HH[i2].ADCperReel[0]->GetBinContent(i3)<0.01){startbin=i3;break;}
// 		}
// 		for(int i3=maxbin;i3<=HH[i2].ADCperReel[0]->GetNbinsX();i3++)
// 		{
// 			if(HH[i2].ADCperReel[0]->GetBinContent(i3)<0.01){endbin=i3;break;}
// 		}
// 		HH[i2].startreel=HH[i2].ADCperReel[0]->GetBinCenter(startbin);
// 		HH[i2].endreel=HH[i2].ADCperReel[0]->GetBinCenter(endbin);
		
		if(!isPhaseIHEP17)
		{
			HH[i2].startreel=HH[i2].ADCperReel[0]->GetBinCenter(HH[i2].ADCperReel[0]->FindFirstBinAbove(0.01));
			HH[i2].endreel=HH[i2].ADCperReel[0]->GetBinCenter(HH[i2].ADCperReel[0]->FindLastBinAbove(0.01));
		}
		else
		{
			HH[i2].startreel=HH[i2].ADCperReel[0]->GetBinCenter(HH[i2].ADCperReel[0]->FindFirstBinAbove(0.1));
			HH[i2].endreel=HH[i2].ADCperReel[0]->GetBinCenter(HH[i2].ADCperReel[0]->FindLastBinAbove(0.1));
		}
// 		cout<<i2<<endl;
// 		HH[i2].ADCperReel[0]->Print("all");
// 		cout<<endl<<endl;
		HH[i2].ADCperReel[0]->Write();
	}
	
	cout<<endl<<"starting second round"<<endl;
	
	vector <int> allreelset;
	bool foundreel=false;
	vector <int> eventlist;
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		bool addtoallreel=false;
		for(int i1=0;i1<ed.mapind->size();i1++)
		{
			sprintf(hname,"%d %d %d %d %s",MP[ed.mapind->at(i1)].ieta,MP[ed.mapind->at(i1)].iphi,MP[ed.mapind->at(i1)].depth,ed.layer,stype[ed.type].c_str());
			string hname1(hname);
			hind=-1;
			for(int i2=0;i2<HH.size();i2++)
			{
				if(HH[i2].name==hname1){hind=i2;break;}
			}
			if(hind==-1) continue;
			if(ed.reel>=HH[hind].startreel && ed.reel<=HH[hind].endreel)
			{
				for(int i4=0;i4<4;i4++)
				{
					if(i4!=ed.badcap->at(i1))
					{
						HH[hind].ADCSpec->Fill(ed.meanadc->at(i1)[i4]);
						if(HH[hind].firstevent1==0) HH[hind].firstevent1=i;
						if(HH[hind].firstevent2==0 && (i-HH[hind].firstevent1)>100) HH[hind].firstevent2=i;
					}
				}
				addtoallreel=true;
			}
		}
		if(addtoallreel)
		{
			foundreel=false;
			for(int is1=0;is1<allreelset.size();is1++)
			{
				if(ed.reel==allreelset[is1]){foundreel=true;break;}
			}
			if(!foundreel){allreelset.push_back(ed.reel);}
		}
		if(i%1000==0) cout<<"Event :"<<i<<" / "<<tree->GetEntries()<<endl;
	}
	cout<<"Max reels : "<<maxreelset.size()<<" All reels : "<<allreelset.size()<<endl;
	//modify maxreel so that the value really corresponds to a reel value
	for(int i1=0;i1<maxreelset.size();i1++)
	{
		bool reelvalid=false;
		while(!reelvalid)
		{
			for(int i2=0;i2<allreelset.size();i2++)
			{
				if(maxreelset[i1]==allreelset[i2]){reelvalid=true;break;}
			}
			if(!reelvalid) maxreelset[i1]+=1;
		}
	}
	cout<<maxreelset.size()<<" reel values to plot "<<endl;
	
	cout<<endl<<"starting third round"<<endl;
	
// 	string pnames[19]={"29","28","27","26","25","24","23","22","21","20_1","20_2","19_1","19_2","18_1","18_2","17_1","17_2","16_1","16_2"};
// 	TH2Poly *p = new TH2Poly("HE","HE",-200,200,-200,200);
// 	TFile *inrootP=new TFile("HEMG.root");
// 	TMultiGraph *mg;
// 	inrootP->GetObject("HE",mg);
// 	TH2Poly *PH;
// // 	PH=(TH2Poly *)p->Clone();
// 	TCanvas* cc1=new TCanvas("cc1","cc1",600,600);
// 	gStyle->SetOptStat(0);
// 	gStyle->SetPalette(kRainBow);
// 	outroot->cd();
// 	for(int i=0;i<tree->GetEntries();i++)
// 	{
// 		tree->GetEntry(i);
// 		bool plotthis=false;
// 		for(int is1=0;is1<maxreelset.size();is1++)
// 		{
// 			if(ed.reel==maxreelset[is1]){plotthis=true;break;}
// 		}
// 		if(plotthis)
// 		{
// 			inrootP->cd();
// 			PH = new TH2Poly("HE","HE",-200,200,-200,200);
// 			inrootP->GetObject("HE",mg);
// 			for(int i1=0;i1<36;i1++)
// 			{
// 				for(int i2=0;i2<19;i2++)
// 				{
// 					sprintf(hname,"W%d_%s",(i1+1),pnames[i2].c_str());
// 					int bin = PH->AddBin(mg->GetListOfGraphs()->FindObject(hname));
// 				}
// 			}
// 			sprintf(hname,"Event=%d Reel=%d",i,ed.reel);
// 			PH->SetName(hname);PH->SetTitle(hname);
// 			for(int i1=0;i1<ed.mapind->size();i1++)
// 			{
// 				sprintf(hname,"%d %d %d %d %s",MP[ed.mapind->at(i1)].ieta,MP[ed.mapind->at(i1)].iphi,MP[ed.mapind->at(i1)].depth,ed.layer,stype[ed.type].c_str());
// 				string hname1(hname);
// 				hind=-1;
// 				for(int i2=0;i2<HH.size();i2++)
// 				{
// 					if(HH[i2].name==hname1){hind=i2;break;}
// 				}
// 				if(hind==-1) continue;
// 				sprintf(hname,"W%d_%s",MP[HH[hind].mapind].box,MP[HH[hind].mapind].tilename.c_str());
// 				for(int i4=0;i4<4;i4++)
// 				{
// 					if(i4!=ed.badcap->at(i1))
// 					{
// 						PH->Fill(hname,ed.meanadc->at(i1)[i4]);break;
// 					}
// 				}
// 			}
// 			PH->Draw("colz");
// 			PH->SetMinimum(-0.05);PH->SetMaximum(0.2);
// 	// 		mg->Draw("same");
// 			outroot->cd();
// 			PH->Write();
// 			sprintf(hname,"PH_%d.png",i+100000);
// 			cc1->SaveAs(hname);
// 			delete PH;
// 		}
// 		if(i%1000==0) cout<<"Event :"<<i<<" / "<<tree->GetEntries()<<endl;
// 	}
	
	cc2->cd();
	gStyle->SetOptTitle(1);
// 	gStyle->SetOptStat(1);
// 	gStyle->SetOptFit(1);
	outroot->cd();
	
	TF1* tf1=new TF1("tf1","gaus",-0.05,0.2);
// 	TF1 *func = new TF1("fit",fitf,0.,0.2,7);
	TF1 *func = new TF1("fit",fitf,0.,20.,7);
	func->SetParameters(30.,0.1,0.01,0.1,0.1,0.1,0.1);
// 	func->FixParameter(2,-0.0652);
// 	func->SetParameter(2,-0.052);
	
	
	double fitparams[7]={0.};double fnchi2=10000.;
	for(int i2=0;i2<HH.size();i2++)
	{
		HH[i2].ADCSpec->SetLineColor(1);
		
		func->SetParameters(30.,0.1,0.01,0.1,0.1,0.1,0.1);
		func->SetParLimits(1,0.,0.2);
		func->SetParLimits(2,0.,1.);
		func->SetParLimits(3,0.,1.);
		func->SetParLimits(5,0.,5.);
		HH[i2].ADCSpec->Fit(func,"q","q",HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindFirstBinAbove(0.)),HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindLastBinAbove(0.)));
// 		HH[i2].ADCSpec->Fit(func,"q","q",0.,0.2);
		if((func->GetChisquare()/func->GetNDF())<fnchi2)
		{
			fnchi2=(func->GetChisquare()/func->GetNDF());
			for(int i1=0;i1<7;i1++)
			{
				fitparams[i1]=func->GetParameter(i1);
			}
		}
	}
	double rchi2=0.;
	for(int i2=0;i2<HH.size();i2++)
	{
		HH[i2].ADCSpec->SetLineColor(1);
		
// 		func->SetParameters(30.,0.1,0.01,0.1,0.1,0.1,0.1);
		func->SetParameter(0,fitparams[0]);
		func->SetParLimits(0,0.,500.);
		func->SetParameter(1,fitparams[1]);
		func->SetParLimits(1,0.,0.2);
		func->SetParameter(2,fitparams[2]);
		func->SetParLimits(2,1e-3,1.);
// 		func->FixParameter(3,fitparams[3]);
		func->SetParameter(3,fitparams[3]);
// 		func->SetParLimits(3,0.,1.);
		func->SetParLimits(3,1e-4,10.);
		func->FixParameter(4,fitparams[4]);
// 		func->FixParameter(5,fitparams[5]);
		func->SetParameter(5,fitparams[5]);
// 		func->SetParLimits(5,0.,5.);
// 		func->SetParLimits(5,0.,10.);
		func->FixParameter(5,fitparams[5]);
		func->FixParameter(6,fitparams[6]);
		
// 		for(int i1=0;i1<7;i1++)
// 		{
// 			func->SetParameter(i1,fitparams[i1]);
// 		}
		
		HH[i2].ADCSpec->Fit(func,"q","q",HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindFirstBinAbove(0.)),HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindLastBinAbove(0.)));
		HH[i2].ADCSpec->Fit(func,"q","q",HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindFirstBinAbove(0.)),HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->FindLastBinAbove(0.)));
// 		HH[i2].ADCSpec->Fit(func,"q","q",0.,0.2);
// 		HH[i2].ADCSpec->Fit(func,"q","q",0.,0.2);
		
		rchi2=func->GetChisquare()/func->GetNDF();
		
// 		tf1->SetParameter(1,HH[i2].ADCSpec->GetBinCenter(HH[i2].ADCSpec->GetMaximumBin()));
// 		tf1->SetParameter(2,0.001);
// 		HH[i2].ADCSpec->Fit(tf1,"q","q",-0.05,0.2);
// 		HH[i2].ADCSpec->Fit(tf1,"q","q",tf1->GetParameter(1)-1.5*tf1->GetParameter(2),tf1->GetParameter(1)+1.5*tf1->GetParameter(2));
// // 		outfile<<MP[HH[i2].mapind].crate<<" "<<MP[HH[i2].mapind].slot<<" "<<MP[HH[i2].mapind].fiber<<" "<<MP[HH[i2].mapind].channel<<" "<<MP[HH[i2].mapind].ieta<<" "<<MP[HH[i2].mapind].iphi<<" "<<MP[HH[i2].mapind].depth<<" "<<HH[i2].layer<<" "<<HH[i2].type<<" "<<tf1->GetParameter(1)<<" "<<tf1->GetParError(1)<<" "<<tf1->GetParameter(2)<<" "<<tf1->GetParError(2)<<endl;
		func->SetNpx(1000);
		HH[i2].ADCSpec->Draw();
		sprintf(hname,"ADCSpectrum_%d.png",i2);
		cc2->SaveAs(hname);
		HH[i2].ADCSpec->Write();
		
		ss.RunNo=RunNo;
		ss.crate=MP[HH[i2].mapind].crate;
		ss.slot=MP[HH[i2].mapind].slot;
		ss.fiber=MP[HH[i2].mapind].fiber;
		ss.channel=MP[HH[i2].mapind].channel;
		ss.ieta=MP[HH[i2].mapind].ieta;
		ss.iphi=MP[HH[i2].mapind].iphi;
		ss.depth=MP[HH[i2].mapind].depth;
		ss.layer=HH[i2].layer;
		ss.SL=HH[i2].type;//0 short; 1 long
		
// 		ss.meanQ=tf1->GetParameter(1);
// 		ss.meanQerr=tf1->GetParError(1);
// 		ss.sigmaQ=tf1->GetParameter(2);
// 		ss.sigmaQerr=tf1->GetParError(2);
		
		ss.meanQ=func->GetParameter(1);
		ss.meanQerr=func->GetParError(1);
		ss.sigmaQ=func->GetParameter(2);
		ss.sigmaQerr=func->GetParError(2);
		
// 		int SSind=-1;
// 		for(int i1=0;i1<SS.size();i1++)
// 		{
// 			if(ss.RunNo==SS[i1].RunNo && ss.crate==SS[i1].crate && ss.slot==SS[i1].slot && ss.fiber==SS[i1].fiber && ss.channel==SS[i1].channel && ss.ieta==SS[i1].ieta && ss.iphi==SS[i1].iphi && ss.depth==SS[i1].depth && ss.layer==SS[i1].layer && ss.SL==SS[i1].SL)
// 			{SSind=i1;break;}
// 		}
// 		if(SSind==-1){SS.push_back(ss);}
// 		else{SS[SSind]=ss;}
		if(rchi2<5.)
		{
			SS.push_back(ss);
		}
	}
	ofstream outfile("MeanQperTile.txt");
// 	outfile<<"RunNo crate slot fiber channel ieta iphi depth layer SL(0)/(1) meanQ meanQerr sigmaQ sigmaQerr"<<endl;
	for(int i1=0;i1<SS.size();i1++)
	{
		outfile<<SS[i1].RunNo<<" "<<SS[i1].crate<<" "<<SS[i1].slot<<" "<<SS[i1].fiber<<" "<<SS[i1].channel<<" "<<SS[i1].ieta<<" "<<SS[i1].iphi<<" "<<SS[i1].depth<<" "<<SS[i1].layer<<" "<<SS[i1].SL<<" "<<SS[i1].meanQ<<" "<<SS[i1].meanQerr<<" "<<SS[i1].sigmaQ<<" "<<SS[i1].sigmaQerr<<endl;
	}
	outfile.close();
	
	sprintf(hname,"mv MeanQperTile.txt ../Plots/%d",RunNo);system(hname);
	
// 	inrootP->Close();
	outroot->cd();
	NormChi2->Write();
	NHgtTh->Write();
	outroot->Close();
	inroot->Close();
}

int main(int argc, char *argv[])
{
// 	int opt=atoi(argv[1]);
	RunNo=atoi(argv[1]);
	readmap();
	
// 	if(opt==1){analyze();}
	analyze();
}

















