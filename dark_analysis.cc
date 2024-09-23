// Convert rawdata into integrated histogram by Hybrid method
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>

void dark_analysis_mod(TString input_file_name, TString output_file_name){
  gStyle -> SetLabelSize(0.04, "xyz");
  gStyle -> SetTitleSize(0.045, "xyz");

  const double range = 2.00; //[V]
  const int depth = 1000; //sample # in 1 cycle
  const int division = 5; //window # in 1 cycle
  const int width = depth/division; //window sample #
  const double bl_fit_half_width = 10; //range of baseline-fit
  double m_th = 25.; //hit threshold from baseline

  //Read rawdata file
  TFile* ifile = new TFile(input_file_name, "read");
  TTree* rawwave = (TTree*)ifile -> Get("fadc_tree");
  const int num_ch = (rawwave -> GetNbranches()) -1;
  short wave[num_ch][depth];
  for(int ich = 0; ich < num_ch; ich++){
      rawwave -> SetBranchAddress(Form("fadc_ch%d", ich), wave[ich]);
  }

  //Create output file as dark.root -> integrated value & no signal or not
  TFile* ofile = new TFile(output_file_name, "recreate");
  TTree* dark_tree[num_ch];
  for(int ich = 0; ich < num_ch; ich++){
    dark_tree[ich] = new TTree(Form("dark_tree_ch%d", ich), Form("dark_tree_ch%d", ich));
  }
  double event;
  double integ;
  double bl;
  double sigma;
  bool hit_flag; //check there's any pilse in a window
  bool pulse_flag; //check pulse is divided by a window
  for(int ich = 0; ich < num_ch; ich++){
    dark_tree[ich] -> Branch("event", &event);
    dark_tree[ich] -> Branch("integ", &integ);
    dark_tree[ich] -> Branch("baseline", &bl);
    dark_tree[ich] -> Branch("bl_sigma", &sigma);
    dark_tree[ich] -> Branch("hit_flag", &hit_flag);
    dark_tree[ich] -> Branch("pulse_flag", &pulse_flag);
  }

  TH1D* hbase = new TH1D("hbase", "", (int)pow(2, 11), 0, (int)pow(2, 12));
  TH1D* hinteg[num_ch];
  for(int ich = 0; ich < num_ch; ich++){
    hinteg[ich] = new TH1D(Form("hinteg_%dch", ich), Form("hinteg_%dch;charge [pC];Entry/bin", ich), 500, -50, 250);
  }
  TCanvas* c = new TCanvas("c", "c", 20, 20, 1200, 800);
  c -> cd();
  for(int ev = 0; ev < rawwave->GetEntries(); ev++){
    if(ev%1000==0&&ev!=0)cout<<"Processing: evt" <<ev<<endl;
    rawwave -> GetEntry(ev);
    event = ev;
    for(int ich = 0; ich < num_ch; ich++){
      hbase -> Reset();
      for(int clk = 0; clk < depth; clk++){
        hbase -> Fill(wave[ich][clk]); //Make histogram for gaussian fitting
      }
      double mean = hbase -> GetBinLowEdge(hbase -> GetMaximumBin());
      hbase -> Fit("gaus", "QNLI", "", mean - bl_fit_half_width, mean + bl_fit_half_width); //Fit with gaussian in the width = 2bl_fit_half_width
      bl = ((TF1*)(gROOT->FindObject("gaus"))) -> GetParameter(1); //Define baseline as mean value 
      sigma = ((TF1*)(gROOT->FindObject("gaus"))) -> GetParameter(2);
      //      th = bl - 5 * sigma; //Define threshold for integration
      double th = bl - m_th;
      if(wave[ich][0]<th) pulse_flag = true;
      else pulse_flag = false;
      for(int idiv = 0; idiv < division; idiv++){
        integ = 0;
        hit_flag = false;
        for(int clk=idiv*width; clk<(idiv + 1)*width; clk++){ //Window method
          integ += (bl - wave[ich][clk]); //Integrate  voltage values
          if(wave[ich][clk]<th && pulse_flag==false) hit_flag = true; //If the voltage value is under threshold, hit_frag -> true: there's signal
	  if(clk==0) continue;
	  else if(wave[ich][clk-1]>th && wave[ich][clk]<th) pulse_flag=true;
	  else if(wave[ich][clk-1]<th && wave[ich][clk]>th) pulse_flag=false;
        }
        integ *= (1. * 0.488 / 50. * 4.); // --> mA x ns = pC: convert voltage into charge
        hinteg[ich] -> Fill(integ);
        dark_tree[ich] -> Fill();
	//if(ev%100==0&&ev!=0) dark_tree[ich] -> Show(); //debug 
      }
    }
  }
  cout << "Process completed!!" <<endl;

  //Output settings
  c->Divide(num_ch);
  for(int ich = 0; ich < num_ch; ich++){
    c -> cd(ich+1);
    hinteg[ich] -> Draw();
    gPad -> SetLeftMargin(0.15);
    gPad -> SetBottomMargin(0.15);
  }
  for(int ich = 0; ich < num_ch; ich++){
    dark_tree[ich] -> Write();
  }
  //ofile -> Close();
}
