//Convert rawdata into integrated histogram
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>

void LED_analysis(TString input_file_name, TString output_file_name){
  gStyle -> SetLabelSize(0.04, "xyz");
  gStyle -> SetTitleSize(0.045, "xyz");

  //parameters
  const double range = 2.00; //[V] for v1724
  const int depth = 1000; //sample # in 1 cycle
  const double bl_fit_half_width = 5;
  const int start = 400; //Start time of LED light
  const int end = 900; //End time of LED light

  //Read rawdata file
  TFile* ifile = new TFile(input_file_name, "read");
  TTree* rawwave = (TTree*)ifile -> Get("fadc_tree");
  const int num_ch = (rawwave -> GetNbranches()) -1;
  short wave[num_ch][depth];
  for(int ich = 0; ich < num_ch; ich++){
      rawwave -> SetBranchAddress(Form("fadc_ch%d", ich), wave[ich]);
  }

  //Create output file as LED.root -> integrated value & no signal or not
  TFile* ofile = new TFile(output_file_name, "recreate");
  TTree* LED_tree[num_ch];
  for(int ich = 0; ich < num_ch; ich++){
    LED_tree[ich] = new TTree(Form("LED_tree_ch%d", ich), Form("LED_tree_ch%d", ich));
  }
  double integ;
  for(int ich = 0; ich < num_ch; ich++){
    LED_tree[ich] -> Branch("integ", &integ);
  }

  //Create LED distribution
  TH1D* hbase = new TH1D("hbase", "", (int)pow(2, 12), 0, (int)pow(2, 12));
  TH1D* hinteg[num_ch];
  for(int ich = 0; ich < num_ch; ich++){
    hinteg[ich] = new TH1D(Form("hinteg_%dch", ich), Form("hinteg_%dch;charge [pC];Entry/bin", ich), 3000, 0, 30000);
  }
  TCanvas* c = new TCanvas("c", "c", 20, 20, 800, 600);
  c -> cd();
  for(int ev = 0; ev < rawwave->GetEntries(); ev++){
    if(ev%100==0&&ev!=0) cout<<"Processing: evt" <<ev<<endl;
    rawwave -> GetEntry(ev);
    for(int ich = 0; ich < num_ch; ich++){
      hbase -> Reset();
      for(int clk = 0; clk < depth; clk++){
        hbase -> Fill(wave[ich][clk]); //Make histogram for gaussian fitting
      }
      double mean = hbase -> GetBinLowEdge(hbase -> GetMaximumBin());
      hbase -> Fit("gaus", "QNLI", "", mean - bl_fit_half_width, mean + bl_fit_half_width); //Fit with gaussian in the width = 2bl_fit_half_width
      double bl = ((TF1*)(gROOT->FindObject("gaus"))) -> GetParameter(1); //Define baseline as gaussian mean value
      double sigma = ((TF1*)(gROOT->FindObject("gaus"))) -> GetParameter(2);
      //if(ev%100==0&&ev!=0) cout<<mean<<endl; 

      // Integration
      integ = 0;
      for(int clk = start; clk < end; clk++){
        integ += (bl - wave[ich][clk]); //Integrate voltage values
      } 
      integ *= (1. * 0.488 / 50. * 4.); // --> mA x ns = pC: convert voltage into charge
      //if(ev%100==0&&ev!=0) cout<<"integ="<<integ<<endl; //debug
      hinteg[ich] -> Fill(integ);
      LED_tree[ich] -> Fill();
      //if(ev%100==0&&ev!=0) LED_tree[ich] -> Show(); //debug
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
    LED_tree[ich] -> Write();
  }
  //ofile -> Close();
}
