//Obtain x-talk, after pulse rate, dark count rtate, and 1 p.e. gain
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TText.h>
#include <TLine.h>

void dark_prop_nocorr(TString input_file_name){
  //parapeter
  double window_width = 800.; // Depend on setting of dark_analysis.cc
  double gain_fit_plus = 8.; //fitting parameter for 1 p.e. peak
  double gain_fit_minus = 8.; //fitting parameter for 1 p.e. peak

  //Read file
  TFile* ifile = new TFile(input_file_name, "read");
  const int num_ch = ifile -> GetNkeys();
  TTree* dark_tree[num_ch];
  for(int ich = 0; ich < num_ch; ich++){
    dark_tree[ich] = (TTree*)ifile -> Get(Form("dark_tree_ch%d", ich));
  }

  TH1D* signal[num_ch];
  TH1D* pedestal[num_ch];
  TH1D* hinteg[num_ch];
  TH1D* hist_gaus[num_ch];
  TH1D* sub_ped[num_ch];

  //Read files
  for (int ich = 0; ich < num_ch; ich++){
    hinteg[ich] = new TH1D(Form("hinteg_%dch", ich), Form("hinteg_%dch;charge [pC];Entry/bin", ich), 175, -100, 250);
    dark_tree[ich] -> Draw(Form("integ>>hinteg_%dch", ich), "", "");
    double a;
    bool flag;
    dark_tree[ich]->SetBranchAddress("integ",&a);
    dark_tree[ich]->SetBranchAddress("hit_flag",&flag);
    signal[ich] = new TH1D(Form("signal%d",ich), "MPPC output events;Charge [pC];Entry/bin",200,-100,300);
    sub_ped[ich] = new TH1D(Form("sub_ped%d",ich), "Difference b/w pulse_dis and its gaus fit;Charge [pC];Entry/bin",200,-100,300);
    pedestal[ich] = new TH1D(Form("pedestal%d",ich), "Pedestal events;Charge [pC];Entry/bin",100,-100,100);
    int N = dark_tree[ich]->GetEntries();

  //Fill the data into histogram
    for (int j=0;j<N;++j){
      dark_tree[ich]->GetEntry(j);
      if (flag==1){ 
	signal[ich]->Fill(a); //Over threshold event = MPPC output
	sub_ped[ich] -> Fill(a);
      }
      else pedestal[ich]->Fill(a);  //Under threshold event = pedestal
    }
    //sub_ped[ich] = (TH1D*)signal[ich]->Clone();
  }
 
  TCanvas* c1 = new TCanvas("c1", "c1", 20, 20, 1200, 800);
  c1 -> Divide(num_ch,1);
  TCanvas* c2 = new TCanvas("c2", "c2", 20, 20, 1200, 800);
  c2 -> Divide(num_ch,1);
  //  TCanvas* c3 = new TCanvas("c3", "c3", 20, 20, 1200, 800);
  //  c3 -> Divide(num_ch,1);
  // TCanvas* c4 = new TCanvas("c4", "c4", 20, 20, 1200, 800);
  //c4 -> Divide(num_ch,1);

  TF1* f = new TF1("f", "gaus", -100, 300);  
  for (int ich=0; ich<num_ch; ich++){
   //Make pedestal destribution
    c1 -> cd(ich+1);
    pedestal[ich] -> Draw();
    double ped_center = pedestal[ich]->GetXaxis()->GetXmin() + pedestal[ich]->GetBinWidth(0)*pedestal[ich]->GetMaximumBin(); //peak search for pedestal (or baseline)
    pedestal[ich] -> Fit("gaus", "Q", "", ped_center-20, ped_center+20);
    double offset = ((TF1*)(gROOT->FindObject("gaus"))) -> GetParameter(1);

   //Make signal distribution
    c2 -> cd(ich+1);
    signal[ich] -> Draw();
    double hist_center = signal[ich]->GetXaxis()->GetXmin() + signal[ich]->GetBinWidth(0)*signal[ich]->GetMaximumBin(); //peak search for 1 p.e. peak
    signal[ich] -> Fit("f","","",hist_center-gain_fit_minus, hist_center+gain_fit_plus);
    gStyle->SetOptStat(0);
    double scale = f -> GetParameter(0);
    double onepe = f -> GetParameter(1);
    double onepe_err =  f -> GetParError(1);
    double onepe_sigma = f -> GetParameter(2);
    double under_onepe = onepe - 1.*onepe_sigma;
    double upper_onepe = onepe + 3*onepe_sigma;

    /*
   //Make 1 p.e. gaus histogram
    c3 -> cd(ich+1);
    int N_gaus = (int)(scale*onepe_sigma*sqrt(2.*acos(-1)));
    hist_gaus[ich] = new TH1D(Form("gaus%d",ich), "Gaus;charge [pC];Entry/bin", 200, -100, 300);
    hist_gaus[ich] -> FillRandom("f", N_gaus/2);
    hist_gaus[ich] -> Draw("n");

   //Get potential pedestal
    c4 -> cd(ich+1);
    sub_ped[ich] -> Add(hist_gaus[ich],-1);
    sub_ped[ich] -> Draw("n");
    TLine *l1 = new TLine(under_onepe,0.1,under_onepe, 100); //0.5* 1 p.e. < dark pulse event
    l1 -> SetLineColor(kRed);
    l1 -> SetLineWidth(3);
    l1->Draw("n");
    */

    double N_ped = (double)(dark_tree[ich]->GetEntries("!hit_flag"));
    double measured_one = signal[ich]->Integral(signal[ich]->FindFixBin(under_onepe), signal[ich]->FindFixBin(upper_onepe)); //Count 1 p.e. dark events as 0.5 p.e. ~1.5 p.e.
    double measured_all = signal[ich]->Integral(signal[ich]->FindFixBin(under_onepe),300);
    double mean =  log((double)(dark_tree[ich]->GetEntries())/N_ped);
    double all_entries = (double) dark_tree[ich]->GetEntries();
    double g_onepe = (onepe - offset)/100/1.602176634e-7*1e-6;
    double g_onepe_err = (onepe_err - offset)/100/1.602176634e-7*1e-6;
    double geff = ((hinteg[ich]->GetMean())-offset)/1.602176634e-7/mean/1e2*1e-6;
    double geff_err = geff*sqrt(pow(1./sqrt(N_ped)/mean,2)+pow(hinteg[ich]->GetMeanError()/(hinteg[ich]->GetMean()-offset),2));

    cout<<"\n"<<endl;
    cout<<"offset: "<<offset<<" [pC]"<<endl;
    cout<<"1 p.e. gain: " << g_onepe<<" [e6]"<<endl;
    cout<<"1 p.e. gain error: "<< g_onepe_err<<" [e6]"<<endl;
    cout<<"g_eff gain: "<< geff <<"[e6]"<<endl;
    cout<<"g_eff error: "<< geff_err <<"[e6]"<<endl;
    cout<<"ca_rate: "<< (1 - g_onepe/geff)*100 <<" %"<<endl;
    cout<<"ca_rate_error: "<<(g_onepe/geff*sqrt(pow(g_onepe_err/g_onepe,2)+pow(geff_err/geff,2)))*100<<" %"<<endl;
    cout<<"Expected DCR: "<< mean/(window_width*1e-9)/1e3<<" [kHz]"<<endl;
    cout<<"DCR error: "<< 1./sqrt(N_ped)/(window_width*1e-9)/1e3<<" [kHz]"<<endl;
    
    /*
    cout<<"Measured dark count rate: "<<measured_all/(window_width*1e-9*(double)(dark_tree[ich]->GetEntries()))/1e3<<" [kHz]"<<endl;
    cout<<"Srd Dev for all window: "<<hinteg[ich]->GetStdDev()<<endl;
    cout<<"measured 1 p.e. event: "<< measured_one<<endl;
    cout<<"expected 1 p.e. event: "<< expected_one<<endl;
    cout<<"Measured all pulse: "<<measured_all<<endl;
    cout<<"Expected all pulse: "<<expected_all<<endl;
    cout<<"Mean pulse #: "<< mean <<endl;
    cout<<"AllEntries: " << dark_tree[ich] -> GetEntries() <<endl;
    cout<<"Pedestal entries: "<< dark_tree[ich]->GetEntries("!hit_flag") <<endl;
    cout<<"All entries: "<< all_entries<<endl;
    cout<<"1 p.e. charge: "<< onepe-offset <<" [pC]"<<endl;
    cout<<"Pulse mean: "<< hinteg[ich]->GetMean() - offset <<" [pC]" << endl;
    cout<<"Mean charge for all entries: "<< (hinteg[ich]->GetMean())-offset <<" [pC]" << endl;
    cout<<"Mean charge error: "<< hinteg[ich]->GetMeanError() << " [pC]" << endl;
    cout<<"Effective charge: "<< ((hinteg[ich]->GetMean())-offset)/mean <<" [pC]"<<endl; 
    
    TText* t1 = new TText(0,0,"");
    t1->SetTextColor(14);
    t1->DrawTextNDC(0.38,0.86,Form("Crosstalk & Afterpulse rate = %.2f %%",ct_ap_rate*100));
    TText* t2 = new TText(0,0,"");
    t2->SetTextColor(14);
    t2->DrawTextNDC(0.38,0.82,Form("1 p.e. gain = %1.2e", gain));
    */
  }
}
