//Create PMT reference photon# vs photon# of MPPC observation
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>

void saturation(){
  //Set parameters
  int n = 30; //data#
  double V[n];
  V[0] = 2.40; //initial Vpp of LED
   for (int i=1; i<n; i++){
    V[i] = V[i-1] + 0.01;
  }  
  int num_ch = 3; //channel# for each data
  double mean[num_ch][n];
  int MPPC_ch[2] = {0, 1}; //0: 25 um, 1: 50 um
  int PMT_ch = 2; //PMT channel
  int Npix[2] = {14400, 3600}; //pixel number
  double co[2] = {1.602176634e-7*4.0e6/4.5, 1.602176634e-7*4.0e6/1.26};  //coefficient for photon #
  double a[2], b[2];

  //Read files and obtain mean values
  for (int ich = 0; ich < num_ch; ich++){
    for (int i = 0; i < n; i++){
      //Load files and set histogram
      TFile* ifile = new TFile(Form("./LED%1.2lfVpp.root", V[i]), "read");
      TTree* LED_tree= (TTree*)ifile -> Get(Form("LED_tree_ch%d", ich));
      TH1D* hinteg = new TH1D(Form("hinteg_%dch", ich), "", 4010, -100, 40000);
      //Fill the data into histogram
      int N = LED_tree -> GetEntries("integ");
      double integ;
      LED_tree -> SetBranchAddress("integ", &integ);
      for (int clk = 0; clk < N; clk++){
        LED_tree -> GetEntry(clk);
	hinteg -> Fill(integ);
      }

      //Obtain mean value: MPPC->photon#, PMT->raw output
      //Using adjusted value of amplifier
      if (ich==0){ 
	if (V[i] < 2.44) {mean[ich][i] = (hinteg -> GetMean())/co[ich]/10.034/10.748;}
	else if (V[i] < 2.59) {mean[ich][i] = (hinteg -> GetMean())/co[ich]/10.034;}
	else {mean[ich][i] = (hinteg -> GetMean())/co[ich];}
      }
      else if (ich==1){ 
	if (V[i] < 2.44) {mean[ich][i] = (hinteg -> GetMean())/co[ich]/10.074/10.607;}
	else if (V[i] < 2.59) {mean[ich][i] = (hinteg -> GetMean())/co[ich]/10.074;}
	else {mean[ich][i] = (hinteg -> GetMean())/co[ich];}
      }
      else{
	if (V[i] < 2.59) mean[ich][i] = (hinteg -> GetMean())/10.393;
	else mean[ich][i] = hinteg -> GetMean();
      }
      //printf("%lf\n", mean[ich][i]);
    }
  }

  //Make reference photon# from PMT: When small light, MPPC photon# will be true
  TGraph*g[2];
  for (int i=0; i<2; i++){
    g[i] = new TGraph();
    for (int j = 0; j < 5; j++){
      g[i] -> SetPoint(j, mean[PMT_ch][j], mean[MPPC_ch[i]][j]);
    }
    //g1 -> Draw("ap");
    g[i] -> Fit("pol1"); //True MPPC photon = p0 + p1*(PMT observation) 
    b[i] = ((TF1*)(gROOT->FindObject("pol1"))) -> GetParameter(0);
    a[i] = ((TF1*)(gROOT->FindObject("pol1"))) -> GetParameter(1);
  }

  //Make Objective plot & Display linear function
  TGraph* gr[2];
  TGraph* g_data[2];
  TMultiGraph* part[2];
  TMultiGraph* sum = new TMultiGraph();

  for (int i = 0; i < 2; i++){
    gr[i] = new TGraph();
    g_data[i] = new TGraph();
    part[i] = new TMultiGraph();
    part[i] -> SetTitle(Form("%d um MPPC;",25*(i+1)));
    for (int j = 0; j < n; j++){
      gr[i] -> SetPoint(j, b[i] + a[i]*mean[PMT_ch][j], b[i] + a[i]*mean[PMT_ch][j]); //Set linear plot for guide
      g_data[i] -> SetPoint(j, b[i] + a[i]*mean[PMT_ch][j], mean[MPPC_ch[i]][j]); //Set target plot
    }
    g_data[i] -> SetMarkerStyle(20);
    if (i == 0){
      g_data[i] -> SetMarkerColor(kRed);
      gr[i] -> SetLineColor(kRed);
    }
    else{
      g_data[i] -> SetMarkerColor(kBlue);
      gr[i] -> SetLineColor(kBlue);
    }
    TCanvas* c = new TCanvas(Form("c%d",i),Form("c%d",i),10,10,900,600);
    g_data[i] -> Draw("ap");
    gr[i] -> Draw("c");
    part[i] -> Add(g_data[i],"p");
    part[i] -> Add(gr[i],"lp");
    part[i] -> Draw("ap");
    part[i] -> GetXaxis() -> SetTitle("N_{ref} [/us]");
    part[i] -> GetYaxis() -> SetTitle("N_{obs} [/us]");
    part[i] -> GetYaxis() -> SetTitleOffset(1.5);
    gPad -> Modified();
  }

  TCanvas* c = new TCanvas("c2","c2",10,10,900,600);
  sum -> SetTitle("Comparison;N_{ref} [/us];N_{obs} [/us]"); 
  for (int i = 0; i < 2; i++){
    sum -> Add(g_data[i],"p");
    sum -> Add(gr[i],"lp");
  }
  sum -> Draw("ap");
  gPad -> Modified();
  sum -> GetXaxis() -> SetLimits(0., 20000.);
  sum -> GetYaxis() -> SetTitleOffset(1.5);
  sum -> SetMaximum(20000.);

  //Make fitting function
  TF1* f[2];
  double k[2];
  for (int i = 0; i < 2; i++){
    f[i] = new TF1(Form("f%d", i), "x/(1. + [0]*x)");
    f[i] -> SetParameter(0,1.39*1e-5);
    if (i == 0) f[i] -> SetLineColor(kRed);
    else f[i] -> SetLineColor(kBlue);
    g_data[i] -> Fit(Form("f%d", i));
    k[i] = ((TF1*)(gROOT->FindObject(Form("f%d",i)))) -> GetParameter(0);
  }

  //Data/Fit plot
  TGraph* residual[2];
  for (int i = 0; i < 2; i++){
    residual[i] = new TGraph();
    for (int j = 0; j < n; j++){
      residual[i] -> SetPoint(j, b[i] + a[i]*mean[PMT_ch][j], mean[MPPC_ch[i]][j]/((b[i] + a[i]*mean[PMT_ch][j])/(1.+ (b[i] + a[i]*mean[PMT_ch][j])*k[i])));
    }
    residual[i] -> SetMarkerStyle(20);
    if (i == 0){
      residual[i] -> SetMarkerColor(kRed);
      residual[i] -> SetTitle("25 um MPPC;Nref[p.e./us];Nobs/Fit");
    }
    else{
      residual[i] -> SetMarkerColor(kBlue);
      residual[i] -> SetTitle("50 um MPPC;Nref[p.e./us];Nobs/Fit");
    }
    TCanvas* c = new TCanvas(Form("res_c%d",i),Form("res_c%d",i),10,10,900,600);
    residual[i] -> GetYaxis() -> SetTitleOffset(1.5);
    residual[i] -> Draw("ap");
    gPad -> Modified();
  }
  
  //Nobs/Nref plot
  TGraph* ratio[2];
  TMultiGraph* ratio_sum = new TMultiGraph();
  for (int i = 0; i < 2; i++){
    ratio[i] = new TGraph();
    for (int j = 0; j < n; j++){
      ratio[i] -> SetPoint(j, b[i] + a[i]*mean[PMT_ch][j], mean[MPPC_ch[i]][j]/ (b[i] + a[i]*mean[PMT_ch][j]));
    }
    ratio[i] -> SetMarkerStyle(20);
    if (i == 0){
      ratio[i] -> SetMarkerColor(kRed);
      ratio[i] -> SetTitle("25 um MPPC;Nref[p.e./us];Nobs/Nref");
    }
    else{
      ratio[i] -> SetMarkerColor(kBlue);
      ratio[i] -> SetTitle("50 um MPPC;Nref[p.e./us];Nobs/Nref");
    }
    TCanvas* c = new TCanvas(Form("ratio_c%d",i),Form("ratio_c%d",i),10,10,900,600);
    ratio[i] -> GetYaxis() -> SetTitleOffset(1.5);
    ratio[i] -> Draw("ap");
    gPad -> Modified();
  }
  TCanvas* c2 = new TCanvas("ratio_sum","ratio_sum",10,10,900,600);
  ratio_sum -> SetTitle("Comparison;N_{ref} [/us];N_{obs}/N_{ref}"); 
  for (int i = 0; i < 2; i++){
    ratio_sum -> Add(ratio[i],"p");
  }
  ratio_sum -> Draw("ap");
  gPad -> Modified();
  //  sum -> GetXaxis() -> SetLimits(0., 20000.);
  
  //Recovery time
  double t;
  for (int i = 0; i < 2; i++){
    t = k[i]*Npix[i]*1e3;
    cout<<"Recovery time = "<<t<<" ns"<<endl;
  }

}

