#include <iostream>
#include <stdio.h>
#include <cmath>

using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TRandom.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TProfile.h>

Double_t xm=0;
Double_t ym=0;


void ToFdE2(){

  gROOT->Reset();

  Double_t TrigThreshFWC1[24] =
    {
      4000,3500,4000,4000,2500,4000,
      3200,3400,3500,4000,2600,2800,
      3500,2500,3500,3500,3000,3200,
      2800,3500,3300,2800,2500,2200
    };

  Double_t TrigThreshFWC2[24] =
    {
      2400,2600,2200,3000,3200,2200,
      2200,2200,1200,1600,1600,1500,
      2000,1700,1900,2200,1500,1600,
      1200,1700,1600,1500,1500,1500
    };
    
    Double_t merrx, merry;
    Double_t cutx, cuty; 
    Int_t nbinsx, nbinsy, nbinx0;
    Bool_t IsFWC1 = kTRUE;



  Double_t ToFFWC1Theta[6] =
    {
      2.088228e+01,2.092216e+01,2.112093e+01,2.149825e+01,2.212018e+01,2.297665e+01
    };
  
  Double_t ToFFWC2Theta[6] =
   {
     1.976750e+01,1.988218e+01,2.025691e+01,2.078644e+01,2.158632e+01,2.263857e+01
   };
  
  Double_t ToFFWCTheta[6];
  for(Int_t i = 0; i<6; i++){
    if(IsFWC1)  ToFFWCTheta[i] = ToFFWC1Theta[i];
   else ToFFWCTheta[i] = ToFFWC2Theta[i];
  }


  TFile*  f = new TFile("run41800dETI.root");
  TString name;
  TCanvas* c1[24];
  TH2D* hToFdE[6][24];
  TH2D* hdE;
  TProfile* pp[24][6];
  TProfile* ppLin[24];
  Int_t noPeakFlag=0;

  TProfile* ppfail[24][6];
  Int_t NEDC=0; //Not Enough Data Counter

  ym=0;
  xm=0;

  for(Int_t i=0; i<24; i++){
    c1[i] = new TCanvas(Form("c1_%02d",i+1), Form("c1_%02d",i+1),1350,950);
    c1[i]->Divide(2,3);
    Int_t FailCounter=0;
    for(Int_t j=0; j<6; j++){
      c1[i]->cd(j+1);

      name = Form("hToFFWCdEFWC_bin%02d_ifwc%02d",j+1,i+1);
      hToFdE[j][i] = (TH2D*)f->Get(name);


      hdE=hToFdE[j][i];

      //err x and y
      //So here I fit x coordinate of THE point
      //Then I calculate where I would like to cut the whole thing off
      TF1* g1 = new TF1("g1","gaus",15.,28.);
      Double_t maxt = ToFFWCTheta[j];
      hdE->ProjectionX("hproj1")->Fit("g1", "IQ", "", maxt-2.5, maxt+2.5);
      merrx = g1->GetParameter(2);
      cutx = maxt+1*merrx;
      

      //Same for y.
      TF1* g2 = new TF1("g2","gaus",0.,4000.);
      Int_t BinDown = hdE->GetXaxis()->FindBin(maxt-2.5);
      Int_t BinUp = hdE->GetXaxis()->FindBin(maxt+2.5);
      Double_t maxy = hdE->ProjectionY()->GetMaximumBin();
      Double_t maxdE = hdE->GetYaxis()->GetBinCenter(maxy);
      hdE->ProjectionY("hproj2",BinDown,BinUp)->Fit("g2", "IQ", "", maxdE-500., maxdE+500);
      merry = g2->GetParameter(2);
      cuty = maxdE+5*merry;
      maxy = g2->GetParameter(1);

      if(abs(g1->GetParameter(0))<200 || abs(g2->GetParameter(0)<200)){
      cout<<"i: "<<i+1<<" j: "<<j+1<<endl;
      cout<<"Amp: "<<g1->GetParameter(0) <<" "<< g2->GetParameter(0)<<endl;
      //cout<<"mean: "<<g1->GetParameter(1) <<" "<< g2->GetParameter(1)<<endl;
      //cout<<"std: "<<merrx<<" "<<merry<<endl;
      //cout<<"cut: "<<cutx<<" "<<cuty<<endl;
      cout<<"Number: "<<hdE->GetEntries()<<endl;
      cout<<"Ratio:"<<hdE->GetEntries()/g1->GetParameter(0) << " "<<hdE->GetEntries()/g2->GetParameter(0)<<endl<<endl;
      }

      //Check if there is a distinct peak.
      if(abs(g1->GetParameter(0))<100 || abs(g2->GetParameter(0)<100))
        noPeakFlag=1;
      else
        noPeakFlag=0;

      nbinsx = hdE->GetXaxis()->FindBin(38.);     
      nbinx0 = hdE->GetXaxis()->FindBin(maxt+cutx);

      //I set range of fitting in y
      Int_t up = hdE->GetYaxis()->GetNbins();
      nbinsy = up;            
      Int_t nentr = hdE->GetEntries();

      for(Int_t counter=0; counter<nbinsy; counter++){
      if(i<maxy) continue;
      if(hdE->ProjectionX("hproj",counter+1,counter+1)->GetEntries()<0.0005*nentr){
          nbinsy = counter+1; 
  		}
      }

    //  cout<<"x range: "<<maxt+4*merrx<<" - "<<38.<<endl;
    //  cout<<"Y range: "<<hdE->GetYaxis()->GetBinCenter(nbinsy)<<endl;
     
      ym = maxdE;
      xm = hdE->GetXaxis()->GetBinCenter(hdE->ProjectionX()->GetMaximumBin());
    //  cout<<" xm, ym: "<<xm<<" "<<ym<<endl;
    //  cout<<endl;

      //If there is a distinct peak, we are fitting function fixing it in the peak.
      if(noPeakFlag==0){
  	  TF1 *g0 = new TF1("g0","[0]*x*x + [1]*x + ([3] - [2]*[2]*[0] - [2]*[1])",0,40);
      g0->SetParameters(1,1,xm,ym);
      g0->FixParameter(2,xm);
      g0->FixParameter(3,ym);
      pp[i][j] = hdE->ProfileX();
      pp[i][j]->Fit("g0","IQ","",cutx, 38.);
      }
      //If not, we can fit pol2 or pol1
      else {
      //TF1 *g0 = new TF1("g0","pol2(0)",0,40);
      TF1 *g0 = new TF1("g0","pol1(0)",0,40);

      //g0->SetParameters(1,1,1);
      g0->SetParameters(1,1);

      pp[i][j] = hdE->ProfileX();
      pp[i][j]->Fit("g0","IQ","",maxt, 38.);

      }


      hdE->GetYaxis()->SetRangeUser(TrigThreshFWC2[i]-400, TrigThreshFWC2[i]+2000.);
      hdE->GetXaxis()->SetRangeUser(0., 50.);

        //If there is a peak, function is drawn in black
      	hdE->Draw("colz");
        if(noPeakFlag==0){
          pp[i][j]->GetFunction("g0")->Draw("same");
          FailCounter=0;
        }
        //If there is no peak, previous function is drawn in red (commented) on its canvas and on top of the current canvas.
        else if(noPeakFlag==1 && j!=0 && FailCounter<=j){
          FailCounter++;
          NEDC++;
          //cout<<"Achtung! "<<i+1<<" "<<j+1<<" Fail no.: "<<FailCounter<<" Total Fails: "<<NEDC<<endl;
          //pp[i][j-FailCounter]->GetFunction("g0")->SetLineColor(2);
          //pp[i][j-FailCounter]->GetFunction("g0")->Draw("same");
          pp[i][j]->GetFunction("g0")->Draw("same");
        }
        //If there is no peak and we are in the first bin, function is drawn in green.
        else if(noPeakFlag==1 && j==0){
          FailCounter++;
          NEDC++;
          //cout<<"Achtung! "<<i+1<<" "<<j+1<<" Fail no.: "<<FailCounter<<" Total Fails: "<<NEDC<<endl;
          pp[i][j]->GetFunction("g0")->SetLineColor(3);
          pp[i][j]->GetFunction("g0")->Draw("same");

        }
        //We shouldn't really get here.
        else{ 
          cout << "___________________SOMETHING WRONG !!!!!!!!!!!!!!"<<endl;
          NEDC++;}
        
      
      if(j==5){    
            //c1[i]->SaveAs(Form("canvas%d.pdf",i+1));
      }
    }
  }

  return;
}
