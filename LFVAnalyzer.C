#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>


#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);


// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------

void LFVAnalyzer() {

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  SetPlotStyle();
  
  
  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  // ----------------------------------------------------------------------------------------------------------------

  TChain* tree = new TChain("IIHEAnalysis");
  tree->Add("STJets_13TeV_LFV_utemu_vector_Madgraph_outfile.root");
  
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
  

  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches
  // ----------------------------------------------------------------------------------------------------------------

  // reconstructed electrons
  vector<float>* el_pt;
  vector<float>* el_eta;
  vector<float>* el_phi;  
  vector<float>* el_charge;  
  vector<int>*   el_isGood;

  TBranch* b_el_pt;
  TBranch* b_el_eta;
  TBranch* b_el_phi;
  TBranch* b_el_charge;
  TBranch* b_el_isGood;

  el_pt = 0;
  el_eta = 0;
  el_phi = 0;
  el_charge = 0;
  el_isGood = 0;
        
  tree->SetBranchAddress("gsf_pt",        &el_pt,        &b_el_pt);
  tree->SetBranchAddress("gsf_eta",       &el_eta,       &b_el_eta);
  tree->SetBranchAddress("gsf_phi",       &el_phi,       &b_el_phi);
  tree->SetBranchAddress("gsf_charge",    &el_charge,    &b_el_charge);
  tree->SetBranchAddress("gsf_VID_cutBasedElectronID_Summer16_80X_V1_medium", &el_isGood, &b_el_isGood);

  // reconstructed muons
  vector<float>* mu_pt;
  vector<float>* mu_eta;
  vector<float>* mu_phi;
  vector<int>*   mu_charge;
  vector<int>*   mu_isGood;

  TBranch* b_mu_pt;
  TBranch* b_mu_eta;
  TBranch* b_mu_phi;
  TBranch* b_mu_charge;
  TBranch* b_mu_isGood;

  mu_pt = 0;
  mu_eta = 0;
  mu_phi = 0;
  mu_charge = 0;
  mu_isGood = 0;

  tree->SetBranchAddress("mu_gt_pt",        &mu_pt,        &b_mu_pt);
  tree->SetBranchAddress("mu_gt_eta",       &mu_eta,       &b_mu_eta);
  tree->SetBranchAddress("mu_gt_phi",       &mu_phi,       &b_mu_phi);
  tree->SetBranchAddress("mu_gt_charge",    &mu_charge,    &b_mu_charge);
  tree->SetBranchAddress("mu_isMediumMuon", &mu_isGood,    &b_mu_isGood);

  // reconstructed jets
  vector<float>* jet_pt;
  vector<float>* jet_eta;
  vector<float>* jet_phi;
  vector<float>* jet_mass;
  
  TBranch* b_jet_pt;
  TBranch* b_jet_eta;
  TBranch* b_jet_phi;
  TBranch* b_jet_mass;

  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_mass = 0;

  tree->SetBranchAddress("jet_pt",        &jet_pt,        &b_jet_pt);
  tree->SetBranchAddress("jet_eta",       &jet_eta,       &b_jet_eta);
  tree->SetBranchAddress("jet_phi",       &jet_phi,       &b_jet_phi);
  tree->SetBranchAddress("jet_mass",      &jet_mass,      &b_jet_mass);

  
  // ----------------------------------------------------------------------------------------------------------------
  // define histograms
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h_nlep = new TH1F("nlep", ";# leptons per event; Events", 5,  -0.5, 4.5);
  TH1F* h_lep1_pt = new TH1F("lep1_pt", ";Lepton p_{T} [GeV]; Leptons / 50 GeV", 20, 0, 600.0);
  TH1F* h_lep2_pt = new TH1F("lep2_pt", ";Lepton p_{T} [GeV]; Leptons / 50 GeV", 20, 0, 600.0);
  TH1F* h_lep3_pt = new TH1F("lep3_pt", ";Lepton p_{T} [GeV]; Leptons / 50 GeV", 20, 0, 600.0);
  

  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  // ----------------------------------------------------------------------------------------------------------------

  for (int i=0; i<nevt; i++) {

    tree->GetEntry(i,0);


    // ----------------------------------------------------------------------------------------------------------------
    // loop over reconstructed detector quantities
    // ----------------------------------------------------------------------------------------------------------------

    int nlep = 0;                                  // counter of number of leptons in each events
    vector<float>* lep_pt = new vector<float>;     // vector to store the pt of the leptons

    int n_el = (int)el_pt->size();   // number of electrons in the event
    int n_mu = (int)mu_pt->size();   // number of muons in the event
    int n_jet = (int)jet_pt->size(); // number of jets in the event
    
    // electrons
    for (int iel=0; iel<n_el; iel++) {

      // only consider electrons that fulfill some minimal criteria, here set as pt>15 GeV, |eta|<2.4 and pass a quality variable
      if (el_pt->at(iel)<15.0 || fabs(el_eta->at(iel))>2.4 || !el_isGood->at(iel)) continue;

      nlep++;
      lep_pt->push_back(el_pt->at(iel));
    }

    // muons
    for (int imu=0; imu<n_mu; imu++) {

      // similar as for electrons, require minimum quality 
      if (mu_pt->at(imu)<15.0 || fabs(mu_eta->at(imu))>2.4 || !mu_isGood->at(imu)) continue;

      nlep++;
      lep_pt->push_back(mu_pt->at(imu));
    }

    // jets
    for (int ij=0; ij<n_jet; ij++) {

      if (jet_pt->at(ij)<30.0 || fabs(jet_eta->at(ij)>2.4)) continue;

      // these are good jets
      
    }


    
    h_nlep->Fill(nlep);  // histogram to stores number of leptons per event

    // ==> events that have three leptons are what we want to look at
    if (nlep == 3) {

      sort(lep_pt->begin(),lep_pt->end()); // here I chose to sort the lepton pt vector ("sort" function sorts in increasing order)
      h_lep1_pt->Fill(lep_pt->at(2));      // leading lepton pt
      h_lep2_pt->Fill(lep_pt->at(1));      // second highest lepton pt in event
      h_lep3_pt->Fill(lep_pt->at(0));      // third lepton pt 

    }

    // ----------------------------------------------------------------------------------------------------------------
    // EXAMPLE OF GETTING THE MASS CORRESPONDING TO THREE PARTICLES
    // ----------------------------------------------------------------------------------------------------------------

    // see e.g. https://root.cern.ch/root/htmldoc/guides/users-guide/PhysicsVectors.html#tlorentzvector

    // initialize lorentzvectors
    TLorentzVector lv_el_1;
    TLorentzVector lv_mu_1;
    TLorentzVector lv_jet_1;
    TLorentzVector lv_elmujet_1;
    
    float el_mass = 0;
    float mu_mass = 0;

    float invmass_1 = 0;
    
    
    // set the pt, eta, phi, mass (here taking the 1st electron / muon / jet just as an example) 
    if (n_el > 0 && n_mu > 0 && n_jet > 0) {
      lv_el_1.SetPtEtaPhiM(el_pt->at(0),el_eta->at(0),el_phi->at(0),el_mass);
      lv_mu_1.SetPtEtaPhiM(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_mass);
      lv_jet_1.SetPtEtaPhiM(jet_pt->at(0),jet_eta->at(0),jet_phi->at(0),jet_mass->at(0));  // jet mass is different for each jet (not like muon/electron)

      lv_elmujet_1 = lv_el_1 + lv_mu_1 + lv_jet_1; // to get the combined electron-muon-jet object, we can add their Lorentz vectors
      invmass_1 = lv_elmujet_1.M();                // this is the invariant mass of the electron-muon-jet object

      //cout << "Invariant mass = " << invmass_1 << " GeV" << endl; 
    }
    
   
  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------


  // canvas for drawing the histograms 
  TCanvas c;

  // draw histogram with number of leptons
  h_nlep->Draw();
  c.SaveAs("nlep.pdf");

  // draw histograms with leading, 2nd highest, 3rd lepton pt in the same canvas
  // first set the y-axis range

  float max = h_lep1_pt->GetMaximum(); 
  if (h_lep2_pt->GetMaximum()>max) max = h_lep2_pt->GetMaximum();
  if (h_lep3_pt->GetMaximum()>max) max = h_lep3_pt->GetMaximum();
  h_lep1_pt->SetAxisRange(0,max*1.05,"Y");  

  h_lep1_pt->Draw();
  h_lep2_pt->SetLineColor(2);
  h_lep2_pt->Draw("same");
  h_lep3_pt->SetLineColor(4);
  h_lep3_pt->Draw("same");

  // add a legend to distinguish the different curves
  TLegend* l = new TLegend(0.7,0.74,0.9,0.89);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_lep1_pt,"1st lepton","l");
  l->AddEntry(h_lep2_pt,"2nd lepton","l");
  l->AddEntry(h_lep3_pt,"3rd lepton","l");
  l->SetTextFont(42);
  l->Draw();	

  c.SaveAs("lep_pt.pdf");


}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
