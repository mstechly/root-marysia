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
#include <TMinuit.h>
#include <TProfile.h>

Double_t xm, ym;
Bool_t noPeakFlag;

Double_t xc(float x, float y, Double_t *par){
  Bool_t complex;
  Double_t xout[3];
  Double_t coef[4];
  Double_t d[3]; 
  //  par[2] = ym - xm*xm*par[0] - xm*par[1];
  coef[3] = 2*par[0]*par[0];
  coef[2] = 3*par[0]*par[1];
  if(!noPeakFlag){
    coef[1] = 1 + par[1]*par[1] + 2*par[0]*(ym - xm*xm*par[0] - xm*par[1]) - 2*par[0]*y;
    coef[0] = par[1]*(ym - xm*xm*par[0] - xm*par[1]) - par[1]*y - x;
  }
  else if(par[0]!=0){
    coef[1] = 1 + par[1]*par[1] + 2*par[0]*par[2] - 2*par[0]*y;
    coef[0] = par[1]*par[2] - par[1]*y - x;
  }
  else{
    return (par[1]*y+x-par[1]*par[2])/(par[1]*par[1]);
  }

  complex = TMath::RootsCubic(coef,xout[0],xout[1],xout[2]);
  if(complex){
    return xout[0];
  }
  else{
    for(Int_t i=0; i<3; i++){
      d[i] = (x-xout[i])*(x-xout[i])+pow(y-par[0]*xout[i]*xout[i]-par[1]*xout[i]-(ym - xout[i]*xout[i]*par[0] - xout[i]*par[1]),2);
      //d[i] = (x-xout[i])*(x-xout[i])+pow(y-par[0]*xout[i]*xout[i]-par[1]*xout[i]-(par[2]),2);
      //      cout<<"distance: "<<d[i]<<endl;
    }
        //cout<<"distance min: "<<d[TMath::LocMin(3,d)]<<" xc: "<<xout[TMath::LocMin(3,d)]<<endl;   
    return xout[TMath::LocMin(3,d)];
 
  }
}

void run_profiles(){

//READING PARAMETERS FROM FILE
Double_t ParameterTable[6][24][3];
Int_t NBinsTable[6][24][4];
Bool_t noPeakTable[6][24];
Double_t xmTable[6][24];
Double_t ymTable[6][24];

ifstream parametersFile ( "parametersFile.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
string value;
Int_t FileRow=0;
Int_t FileCol=0;
Int_t FileEl=0;
Int_t FileBin=0;

string HelpTable[10];

	for(Int_t j=0; j<24; j++){
		for(Int_t i=0; i<6; i++){
			for(Int_t k=0; k<10; k++){
				getline ( parametersFile, value, ',' ); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
    			//cout << string( value, 0, value.length() )<<" "; // display value removing the first and the last character from it
    			HelpTable[k]=value;
    			//cout<<" "<<i<<" "<<j<<" "<<k<<" "<<HelpTable[k]<<endl;
			}
			xmTable[i][j]=atof(HelpTable[0].c_str());
			ymTable[i][j]=atof(HelpTable[1].c_str());
			noPeakTable[i][j]=atof(HelpTable[2].c_str());
			NBinsTable[i][j][0]=atof(HelpTable[3].c_str());
			NBinsTable[i][j][1]=atof(HelpTable[4].c_str());
			NBinsTable[i][j][2]=atof(HelpTable[5].c_str());
			NBinsTable[i][j][3]=atof(HelpTable[6].c_str());
			ParameterTable[i][j][0]=atof(HelpTable[7].c_str());
			ParameterTable[i][j][1]=atof(HelpTable[8].c_str());
			ParameterTable[i][j][2]=atof(HelpTable[9].c_str());
		}	
	}
    


//hdEToF z każdego runu


//MAIN HISTOGRAM CREATING PART

	file = new TFile("run41800dETIFWC1.root");
	TH2D* hdEToF[6][24];
  	TCanvas* c[24];
  	TH1D* DistHist[6][24];
	TF1*  fFit[6][24];



	Int_t StartX, EndX, StartY, EndY;
  	Double_t CurrentX, CurrentY, CurrentW;
  	Double_t CurrentA, CurrentB, CurrentC;
  	Double_t CurrentD, CurrentXC, CurrentYC;

  	for(Int_t el=1; el<25; el++){
    	c[el-1] = new TCanvas(Form("c_el%02d",el),Form("c_el%02d",el),1350,950);
    	c[el-1]->Divide(4,3);
    	for(Int_t bin = 1; bin<7; bin++){

    		hdEToF[bin-1][el-1] = (TH2D*)file->Get(Form("hToFFWCdEFWC_bin%02d_ifwc%02d",bin,el));
    		c[el-1]->cd(2*bin-1);
      		DistHist[bin-1][el-1]= new TH1D(Form("DistHist_%02d_%02d",el,bin),Form("DistHist_%02d_%02d",el,bin),150,-15,15);
     		fFit[bin-1][el-1] = new TF1(Form("fFit_bin%02d_el%02d",bin,el), "pol2", 0., 40.); 
     		fFit[bin-1][el-1]->FixParameter(0,ParameterTable[bin-1][el-1][2]);
     		fFit[bin-1][el-1]->FixParameter(1,ParameterTable[bin-1][el-1][1]);
     		fFit[bin-1][el-1]->FixParameter(2,ParameterTable[bin-1][el-1][0]);


      		StartX=NBinsTable[bin-1][el-1][0];
      		EndX=NBinsTable[bin-1][el-1][1];
      		StartY=NBinsTable[bin-1][el-1][2];
      		EndY=NBinsTable[bin-1][el-1][3];
      		noPeakFlag=noPeakTable[bin-1][el-1];
      		CurrentA=ParameterTable[bin-1][el-1][0];
      		CurrentB=ParameterTable[bin-1][el-1][1];
      		CurrentC=ParameterTable[bin-1][el-1][2];
      		xm=xmTable[bin-1][el-1];
      		ym=ymTable[bin-1][el-1];



    		for(Int_t i=StartX; i<=EndX; i++){
       			for(Int_t j=StartY; j<=EndY; j++){
            		CurrentX=hdEToF[bin-1][el-1]->GetXaxis()->GetBinCenter(i);
            		CurrentY=hdEToF[bin-1][el-1]->GetYaxis()->GetBinCenter(j);
            		CurrentW=hdEToF[bin-1][el-1]->GetBinContent(i,j);
            		CurrentXC=xc(CurrentX,CurrentY,ParameterTable[bin-1][el-1]);
            		CurrentYC=CurrentA*pow(CurrentXC,2.0)+CurrentB*CurrentXC+CurrentC;
            		CurrentD=sqrt(pow((CurrentXC-CurrentX),2.0)+pow((CurrentYC-CurrentY),2.0));
            		if(CurrentXC>CurrentX)
              			CurrentD=-CurrentD;
            		DistHist[bin-1][el-1]->Fill(CurrentD,CurrentW);
            
          		}
        	}
      		hdEToF[bin-1][el-1]->Draw("colz");
      		fFit[bin-1][el-1]->Draw("same");

      		c[el-1]->cd(2*bin);
      		DistHist[bin-1][el-1]->Draw("colz");
      		TF1* ProfileGauss = new TF1("ProfileGauss","gaus",-15.,15.);  
      		DistHist[bin-1][el-1]->Fit("ProfileGauss", "IQ", "", 15., 15.);
      		ProfileGauss->SetLineColor(2);
      		ProfileGauss->Draw("same");
      		cout<<"el: "<<el<<" bin: "<<bin;
      		cout<<" Amplitude: "<<ProfileGauss->GetParameter(0)<<" Mean: "<<ProfileGauss->GetParameter(1)<<" STD: "<<ProfileGauss->GetParameter(2)<<endl;

      		if(bin==6){    
        		//c[el-1]->SaveAs(Form("canvas%d.pdf",el));
      		}
		}
	}
	parametersFile.close();
}