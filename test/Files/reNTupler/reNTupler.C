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
edata ed;

struct capdata
{
	string name;
	int mapind;
	float base[4];
	float baseerr[4];
	int Nevt[4];
	int badcap;
	float max[4];
	TH1F* QperEvent[4];
	TH1F* QperEventBS[4];
};
capdata cd;
vector <capdata> CD;

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
NTupledata ND;

int RunNo=0;

double adc2fC_QIE10[256]={1.58, 4.73, 7.88, 11.0, 14.2, 17.3, 20.5, 23.6,26.8, 29.9, 33.1, 36.2, 39.4, 42.5, 45.7, 48.8,53.6, 60.1, 66.6, 73.0, 79.5, 86.0, 92.5, 98.9,105, 112, 118, 125, 131, 138, 144, 151,157, 164, 170, 177, 186, 199, 212, 225, 238, 251, 264, 277, 289, 302, 315, 328,341, 354, 367, 380, 393, 406, 418, 431,444, 464, 490, 516, 542, 568, 594, 620,569, 594, 619, 645, 670, 695, 720, 745,771, 796, 821, 846, 871, 897, 922, 947,960, 1010, 1060, 1120, 1170, 1220, 1270, 1320,1370, 1430, 1480, 1530, 1580, 1630, 1690, 1740,1790, 1840, 1890, 1940,2020, 2120, 2230, 2330,2430, 2540, 2640, 2740, 2850, 2950, 3050, 3150, 3260, 3360, 3460, 3570, 3670, 3770, 3880, 3980,4080, 4240, 4450, 4650, 4860, 5070, 5280, 5490,5080, 5280, 5480, 5680, 5880, 6080, 6280, 6480,6680, 6890, 7090, 7290, 7490, 7690, 7890, 8090,8400, 8810, 9220, 9630, 10000, 10400, 10900, 11300,11700, 12100, 12500, 12900, 13300, 13700, 14100, 14500,15000, 15400, 15800, 16200, 16800, 17600, 18400, 19300,20100, 20900, 21700, 22500, 23400, 24200, 25000, 25800,26600, 27500, 28300, 29100, 29900, 30700, 31600, 32400,33200, 34400, 36100, 37700, 39400, 41000, 42700, 44300,41100, 42700, 44300, 45900, 47600, 49200, 50800, 52500,54100, 55700, 57400, 59000, 60600, 62200, 63900, 65500,68000, 71300, 74700, 78000, 81400, 84700, 88000, 91400,94700, 98100, 101000, 105000, 108000,111000, 115000, 118000,121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000};

double adc2fC_QIE8=0.32;

int readmap()
{
	ifstream infile;
	if(RunNo<287582){infile.open("../semap_v1.txt");}
	else{infile.open("../semap_v2.txt");}
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

string between(string const &in, string const &before, string const &after)
{
	return in.substr(in.find(before)+before.size(), 2);
}

int reNTuple()
{
	char hname[500];
	sprintf(hname,"N2_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	TTree* outtree = new TTree("Events", "Events");
	outtree->Branch("reel", &ND.reel);
	outtree->Branch("layer", &ND.layer);
	outtree->Branch("type", &ND.type);
	outtree->Branch("name", &ND.name);
	outtree->Branch("mapind", &ND.mapind);
	outtree->Branch("badcap", &ND.badcap);
	outtree->Branch("meanadc", &ND.meanadc);
	outtree->Branch("meanadcerr", &ND.meanadcerr);
	outtree->Branch("base", &ND.base);
	outtree->Branch("baseerr", &ND.baseerr);
	outtree->Branch("N", &ND.N);
	
	sprintf(hname,"H2_%d.root",RunNo);
	TFile* outroot2=new TFile(hname,"recreate");
	TH1F* NormChi2=new TH1F("NormChi2","NormChi2",100,0.,20.);
	TH1F* MeanADCSpec=new TH1F("MeanQSpec","MeanQSpec",2000,0.,20.);
	
	sprintf(hname,"../../NTuples/N_%d.root",RunNo);
	TFile* inroot=new TFile(hname);
	TTree *tree = (TTree*)inroot->Get("Events");
	tree->SetBranchAddress("name",&ed.name);
	tree->SetBranchAddress("reel",&ed.reel);
	tree->SetBranchAddress("adc",&ed.adc);
	tree->SetBranchAddress("cname",&ed.cname);
	
	string name;vector <int> loc;
	int crate=0;int slot=0;int fiber=0;int channel=0;
	int MI=-1;
	
	TCanvas* cc1=new TCanvas("cc1","cc1",600,600);
	TF1* tf1=new TF1("tf1","gaus",0.,35.);
	TF1* tf2=new TF1("tf2","[0]",0.,30000.);
	
	tree->GetEntry(0);
	string hnamepre="";
	int capind=0;int Nevt=0;int cdind=0;
	for(int i1=0;i1<ed.cname->size();i1++)
	{
		name=ed.cname->at(i1);
		loc=findLocation(name,'_');
		crate=atoi(name.substr(loc[0]+1,loc[1]-loc[0]-1).c_str());
		slot=atoi(name.substr(loc[1]+1,loc[2]-loc[1]-1).c_str());
		fiber=atoi(name.substr(loc[2]+1,loc[2]-loc[2]-1).c_str());
		channel=atoi(name.substr(loc[3]+1,name.size()-loc[3]-1).c_str());
		loc.clear();
		MI=-1;
		for(int i2=0;i2<MP.size();i2++)
		{
			if(MP[i2].crate==crate && MP[i2].slot==slot && MP[i2].fiber==fiber && MP[i2].channel==channel)
			{
				MI=i2;break;
			}
		}
		if(MI==-1) continue;
		cdind=-1;
		for(int is1=0;is1<CD.size();is1++)
		{
			if(name==CD[is1].name){cdind=is1;break;}
		}
		if(cdind==-1)
		{
			cd.name=name;
			cd.mapind=MI;
			cd.badcap=-1;
			for(int i4=0;i4<4;i4++)
			{
				cd.base[i4]=0;
				cd.baseerr[i4]=0;
				cd.Nevt[i4]=0;
				cd.max[i4]=0;
				sprintf(hname,"QperEvent_%d_%d_%d_cap%d",MP[MI].ieta,MP[MI].iphi,MP[MI].depth,i4);
				cd.QperEvent[i4]=new TH1F(hname,hname,tree->GetEntries(),-0.5,tree->GetEntries()-0.5);
				cd.QperEvent[i4]->GetXaxis()->SetTitle("Event ID");cd.QperEvent[i4]->GetXaxis()->CenterTitle();
				cd.QperEvent[i4]->GetYaxis()->SetTitle("Mean Charge / Event (fC)");cd.QperEvent[i4]->GetYaxis()->CenterTitle();
				sprintf(hname,"QperEventBS_%d_%d_%d_cap%d",MP[MI].ieta,MP[MI].iphi,MP[MI].depth,i4);
				cd.QperEventBS[i4]=new TH1F(hname,hname,tree->GetEntries(),-0.5,tree->GetEntries()-0.5);
				cd.QperEventBS[i4]->GetXaxis()->SetTitle("Event ID");cd.QperEventBS[i4]->GetXaxis()->CenterTitle();
				cd.QperEventBS[i4]->GetYaxis()->SetTitle("Mean Charge / Event (fC)");cd.QperEventBS[i4]->GetYaxis()->CenterTitle();
			}
			CD.push_back(cd);
			cdind=CD.size()-1;
		}
		if(name!=hnamepre) {capind=0;}
		else capind++;
		Nevt=0;
		bool badcap=false;
		for(int i2=0;i2<ed.adc->at(i1).size();i2++)
		{
			Nevt+=ed.adc->at(i1)[i2];
			if(i2>35 && ed.adc->at(i1)[i2]>0) badcap=true;
		}
		CD[cdind].Nevt[capind]=Nevt;
		if(badcap) CD[cdind].badcap=capind;
		hnamepre=name;
	}
// 	cout<<CD.size()<<endl;
	TH1F* hh;
	for(int i=0;i<tree->GetEntries();i++)
// 	for(int i=0;i<100;i++)
	{
		tree->GetEntry(i);
		string hnamepre="";int capind=0;
		for(int i1=0;i1<ed.cname->size();i1++)
		{
			name=ed.cname->at(i1);
			cdind=-1;
			for(int i2=0;i2<CD.size();i2++)
			{
				if(CD[i2].name==name){cdind=i2;break;}
			}
			if(cdind==-1)continue;
			
			if(name!=hnamepre) capind=0;
			else capind++;
			MI=CD[cdind].mapind;
			hh=new TH1F("hh","hh",2000,0.,200.);
			if(RunNo>=287582 && MP[MI].iphi>=63 && MP[MI].iphi<=66 && MP[MI].ieta>0)
			{
				for(int iz2=0;iz2<35;iz2++)
				{
					hh->Fill(adc2fC_QIE10[iz2],ed.adc->at(i1)[iz2]);
				}
			}
			else
			{
				for(int iz2=0;iz2<35;iz2++)
				{
					hh->Fill(((double)iz2)*adc2fC_QIE8,ed.adc->at(i1)[iz2]);
				}
			}
			CD[cdind].QperEvent[capind]->SetBinContent(i+1,hh->GetMean());
			CD[cdind].QperEvent[capind]->SetBinError(i+1,hh->GetRMS()/sqrt(hh->Integral()));
			delete hh;
			hnamepre=name;
		}
		if(i%1000==0) cout<<"Event :"<<i<<" / "<<tree->GetEntries()<<endl;
	}
	
	cout<<endl<<"Starting second round"<<endl<<endl;
	
	double nchi2=0.;double base=0.;double baseerr=0.;int nmaxOK=0;
	
	outroot2->cd();
	for(int i2=0;i2<CD.size();i2++)
	{
		for(int i4=0;i4<4;i4++)
		{
			CD[i2].QperEvent[i4]->Fit(tf2,"q","q",0.,tree->GetEntries());
			nchi2=tf2->GetChisquare()/tf2->GetNDF();
			NormChi2->Fill(nchi2);
			
			CD[i2].QperEvent[i4]->Write();
			CD[i2].QperEvent[i4]->Fit(tf2,"q","q",0.,200.);
			base=tf2->GetParameter(0);
			baseerr=tf2->GetParError(0);
			CD[i2].base[i4]=base;
			CD[i2].baseerr[i4]=baseerr;
			CD[i2].max[i4]=CD[i2].QperEvent[i4]->GetBinContent(CD[i2].QperEvent[i4]->GetMaximumBin())-base;
		}
		nmaxOK=0;
		for(int i4=0;i4<4;i4++)
		{
			if(CD[i2].max[i4]>0.015) nmaxOK++;
		}
		if(nmaxOK>1)
		{
			for(int i4=0;i4<4;i4++)
			{
				for(int iz2=1;iz2<=CD[i2].QperEvent[i4]->GetNbinsX();iz2++)
				{
					CD[i2].QperEventBS[i4]->SetBinContent(iz2,CD[i2].QperEvent[i4]->GetBinContent(iz2)-CD[i2].base[i4]);
					CD[i2].QperEventBS[i4]->SetBinError(iz2,sqrt(pow(CD[i2].QperEvent[i4]->GetBinError(iz2),2.)+pow(CD[i2].baseerr[i4],2.)));
					MeanADCSpec->Fill(CD[i2].QperEventBS[i4]->GetBinContent(iz2));
				}
				CD[i2].QperEventBS[i4]->Write();
			}
		}
		else
		{
			CD[i2]=CD[CD.size()-1];
			CD.pop_back();
			i2--;
		}
	}
	
	cout<<endl<<"Starting final round"<<endl<<endl;
	
	vector <float> meanadcs;
	vector <float> meanadcerrs;
	vector <float> bases;
	vector <float> baseerrs;
	vector <int> Ns;
	for(int i=0;i<tree->GetEntries();i++)
	{
		tree->GetEntry(i);
		string hnamepre="";int capind=0;
		ND.mapind->clear();
		ND.badcap->clear();
		ND.meanadc->clear();
		ND.meanadcerr->clear();
		ND.base->clear();
		ND.baseerr->clear();
		ND.N->clear();
		
		ND.reel=ed.reel;
		ND.name=ed.name;
		name=*ed.name;
		
		ND.layer=atoi(between(name,"LAYER","_").c_str());
		if(name.back()=='S') ND.type=0;
		else ND.type=1;
		
		for(int i2=0;i2<CD.size();i2++)
		{
			meanadcs.clear();
			meanadcerrs.clear();
			bases.clear();
			baseerrs.clear();
			Ns.clear();
			for(int i4=0;i4<4;i4++)
			{
				meanadcs.push_back(CD[i2].QperEventBS[i4]->GetBinContent(i+1));
				meanadcerrs.push_back(CD[i2].QperEventBS[i4]->GetBinError(i+1));
				bases.push_back(CD[i2].base[i4]);
				baseerrs.push_back(CD[i2].baseerr[i4]);
				Ns.push_back(CD[i2].Nevt[i4]);
			}
			ND.mapind->push_back(CD[i2].mapind);
			ND.badcap->push_back(CD[i2].badcap);
			ND.meanadc->push_back(meanadcs);ND.meanadcerr->push_back(meanadcerrs);
			ND.base->push_back(bases);ND.baseerr->push_back(baseerrs);
			ND.N->push_back(Ns);
		}
		outtree->Fill();
		if(i%1000==0) cout<<"Event :"<<i<<" / "<<tree->GetEntries()<<endl;
	}
	outroot->Write();
	outroot->Close();
	
	outroot2->cd();
	
	NormChi2->GetXaxis()->SetTitle("#chi^{2} / ndf");NormChi2->GetXaxis()->CenterTitle();
	NormChi2->GetYaxis()->SetTitle("Events / 0.05");NormChi2->GetYaxis()->CenterTitle();
	NormChi2->Write();
	
	MeanADCSpec->GetXaxis()->SetTitle("Mean Charge (fC)");MeanADCSpec->GetXaxis()->CenterTitle();
	MeanADCSpec->GetYaxis()->SetTitle("Events / 10^{-3}");MeanADCSpec->GetYaxis()->CenterTitle();
	MeanADCSpec->Write();
	
	outroot2->Close();
	inroot->Close();
}

int main(int argc, char *argv[])
{
	RunNo=atoi(argv[1]);
	readmap();
	
	reNTuple();
}

















