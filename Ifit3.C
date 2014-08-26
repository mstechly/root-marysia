//Basing on:
//   Example of a program to fit non-equidistant data points
//   =======================================================
//
//   The fitting function fcn is a simple chisquare function
//   The data consists of 5 data points (arrays x,y,z) + the errors in errorsz
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//Author: Rene Brun
//2nd Author: Maria Zurek

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
#include "TMinuit.h"
#include <TProfile.h>

TFile *file;
TH2D  *hdE;
Double_t merrx, merry;
Double_t cutx, cuty; 
Double_t xm, ym;
Int_t nbinsx, nbinsy, nbinx0;

Double_t TrigThreshFWC2[24] =
  {
    2400,2600,2200,3000,3200,2200,
    2200,2200,1200,1600,1600,1500,
    2000,1700,1900,2200,1500,1600,
    1200,1700,1600,1500,1500,1500
  };

//______________________________________________________________________________
Double_t func(float x ,Double_t *par){
  //par[2] = ym - xm*xm*par[0] - xm*par[1];
  Double_t value=par[0]*x*x + par[1]*x + (ym - xm*xm*par[0] - xm*par[1]);
 return value;
}
//______________________________________________________________________________
Double_t xc(float x, float y, Double_t *par){
  Bool_t complex;
  Double_t xout[3];
  Double_t coef[4];
 
  //  par[2] = ym - xm*xm*par[0] - xm*par[1];
  coef[3] = 2*par[0]*par[0];
  coef[2] = 3*par[0]*par[1];
  coef[1] = 1 + par[1]*par[1] + 2*par[0]*(ym - xm*xm*par[0] - xm*par[1]) - 2*par[0]*y;
  coef[0] = par[1]*(ym - xm*xm*par[0] - xm*par[1]) - par[1]*y - x;

  complex = TMath::RootsCubic(coef,xout[0],xout[1],xout[2]);
  if(complex){
    return xout[0];
  }
  else{
    return TMath::MinElement(3,xout);
  }
}
//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  Int_t i, j;

  //calculate chisquare
   Double_t chisq = 0;
   Double_t delta1, delta2;
   Double_t xp, yp, wp;
   Double_t x0;

   for(j=nbinx0; j<=nbinsx;j++){
     xp = hdE->GetXaxis()->GetBinCenter(j);
     for (i=1;i<=nbinsy; i++) {
       yp = hdE->GetYaxis()->GetBinCenter(i);
       //   if(xp < cutx && yp > cuty) continue;
       wp = hdE->GetBinContent(j,i);      
       //wp =1.;
       //if(xp < cutx){
       //wp*=0.1;
       //}
       x0 = xc(xp,yp,par);
       delta1  = (xp-x0)/merrx;
       delta2  = (yp-func(x0,par))/merry;      
       Double_t delta = (delta1*delta1+delta2*delta2);
       //if (delta>4) delta *= exp(9-delta);
       chisq += wp*delta;
     }
   }
   f = chisq;
}

//______________________________________________________________________________
void Ifit3()
{
  //Do I need it?
  //gROOT->Reset();

  file = new TFile("run41800dETIFWC1.root");
  Bool_t IsFWC1 = kTRUE;

  TMinuit *minuit = new TMinuit(3);
  minuit->SetFCN(fcn);
  minuit->SetPrintLevel(-1);
  Double_t arglist[10];
  Int_t ierflg;      
     
  //ToF peak position from MC
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

  //  Int_t nbin = 3;

  TH2D* hdEToF[6][24];
  TF1*  fFit[6][24];
  TProfile* pp[24][6];
  Int_t noPeakFlag=0;


 
  for(Int_t el = 1; el<25; el++){
    Double_t vstart[3]; 
    for(Int_t bin = 1; bin<7; bin++){
       
      hdEToF[bin-1][el-1] = (TH2D*)file->Get(Form("hToFFWCdEFWC_bin%02d_ifwc%02d",bin,el));
      hdE = hdEToF[bin-1][el-1];

      //err x and y
      TF1* g1 = new TF1("g1","gaus",15.,28.);  
      Double_t maxt = ToFFWCTheta[bin-1];
      hdE->ProjectionX("hproj1")->Fit("g1", "IQ", "", maxt-2.5, maxt+2.5);
      merrx = g1->GetParameter(2);
      cutx = maxt+1*merrx;
      
      TF1* g2 = new TF1("g2","gaus",0.,4000.);
      Int_t BinDown = hdE->GetXaxis()->FindBin(maxt-2.5);
      Int_t BinUp = hdE->GetXaxis()->FindBin(maxt+2.5);
      Double_t maxy = hdE->ProjectionY()->GetMaximumBin();
      Double_t maxdE = hdE->GetYaxis()->GetBinCenter(maxy);
      hdE->ProjectionY("hproj2",BinDown,BinUp)->Fit("g2", "IQ", "", maxdE-500., maxdE+500);
      merry = g2->GetParameter(2);
      cuty = maxdE+5*merry;
      maxy = g2->GetParameter(1);

      if(abs(g1->GetParameter(0))<200 || abs(g2->GetParameter(0)<200))
        {noPeakFlag=1;
          cout<<"TUTAJ TUTAJ ________________________________________________________"<<endl;
        }
      else
        noPeakFlag=0;
      //cout<<"err: "<<merrx<<" "<<merry<<endl;
      //cout<<"cut: "<<cutx<<" "<<cuty<<endl;
      //Ranges of fitting 
      nbinsx = hdE->GetXaxis()->FindBin(38.);     
      nbinx0 = hdE->GetXaxis()->FindBin(maxt+4*merrx);

      Int_t up = hdE->GetYaxis()->GetNbins();
      nbinsy = up;            
      Int_t nentr = hdE->GetEntries();
      for(Int_t i=0; i<nbinsy; i++){
			if(i<maxy) continue;
			if(hdE->ProjectionX("hproj",i+1,i+1)->GetEntries()<0.0005*nentr){
	  			nbinsy = i+1; 
	}
      }
      
      //cout<<"x range: "<<maxt+4*merrx<<" - "<<38.<<endl;
      //cout<<"Y range: "<<hdE->GetYaxis()->GetBinCenter(nbinsy)<<endl;
     
      //3Hen peak position
      ym = maxdE;
      xm =  hdE->GetXaxis()->GetBinCenter(hdE->ProjectionX()->GetMaximumBin());
      cout<<" xm, ym: "<<xm<<" "<<ym<<endl;
      /*
      //Fit to profile - start parameters
      TF1* g0 = new TF1("g0","pol2",15.,38.);
      TProfile* pp;
      pp = hdE->ProfileX();
      pp->Fit("g0","I","0",maxt-5*merrx,38.);   
      */

      if(noPeakFlag==0){
      TF1 *g0 = new TF1("g0","[0]*x*x + [1]*x + ([3] - [2]*[2]*[0] - [2]*[1])",0,40);
      g0->SetParameters(1,1,xm,ym);
      pp[el-1][bin-1] = hdE->ProfileX();
      pp[el-1][bin-1]->Fit("g0","IQ","",cutx, 38.);
      }
      //If not, we can fit pol2 or pol1
      else {
      TF1 *g0 = new TF1("g0","pol2(0)",0,40);
      //TF1 *g0 = new TF1("g0","pol1(0)",0,40);

      g0->SetParameters(1,1,1);
      //g0->SetParameters(1,1);
      pp[el-1][bin-1] = hdE->ProfileX();

      //ACHTUNG! Arbitralnie wybrane 20 !!!
      pp[el-1][bin-1]->Fit("g0","IQ","",20., 38.);

      }


      //initialize TMinuit with a maximum of 3 params
      ierflg = 0;
      arglist[0] = 1;
      minuit->mnexcm("SET ERR", arglist ,1,ierflg);
      
      // Set starting values and step sizes for parameters
      if(bin<5){
	vstart[0] = g0->GetParameter(2);
	vstart[1] = g0->GetParameter(1);
	//	vstart[2] = g0->GetParameter(0);
      }

      Double_t step[3] = {0.01 , 1. , 10.};
      minuit->mnparm(0, "a", vstart[0], step[0], 0,0,ierflg);
      minuit->mnparm(1, "b", vstart[1], step[1], 0,0,ierflg);
      // minuit->mnparm(2, "c", vstart[2], step[2], 0,0,ierflg);
      
      // Now ready for minimization step
      arglist[0] = 1500;
      arglist[1] = 1;
      minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      
      // Print results
      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
      //gMinuit->mnprin(3,amin);
      
      // get parameters
      
      double pval[2],perr[2],plo[2],phi[2];
      
      TString para0,para1,para2;
      
      int istat;
      
      minuit->mnpout(0,para0,pval[0],perr[0],plo[0],phi[0],istat);
      minuit->mnpout(1,para1,pval[1],perr[1],plo[1],phi[1],istat);
      //minuit->mnpout(2,para2,pval[2],perr[2],plo[2],phi[2],istat);
      
      fFit[bin-1][el-1] = new TF1(Form("fFit_bin%02d_el%02d",bin,el), "pol2", 0., 40.); 
      fFit[bin-1][el-1]->FixParameter(2, pval[0]);
      fFit[bin-1][el-1]->FixParameter(1, pval[1]);
      fFit[bin-1][el-1]->FixParameter(0, ym - pval[0]*xm*xm - pval[1]*xm);
    }
  }

  TCanvas* c[24];
  for(Int_t el=1; el<25; el++){
    c[el-1] = new TCanvas(Form("c_el%02d",el),Form("c_el%02d",el),1350,950);
    c[el-1]->Divide(2,3);
    for(Int_t bin = 1; bin<7; bin++){  
      c[el-1]->cd(bin);
      hdEToF[bin-1][el-1]->Draw("colz");
      fFit[bin-1][el-1]->Draw("same");
    }  
  }

}


