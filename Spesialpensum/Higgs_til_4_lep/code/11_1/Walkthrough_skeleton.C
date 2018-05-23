/*
 * Project:        Exercises 11.1-11.3
 * File:           Walkthrough_skeleton.C
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * A code skeleton for the searches part.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */


#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <iostream>
using namespace std;

//-------------------------------------------------------------------------------------------------------------------------
//-- full functions
TH1D * GetMassDistribution(int Itype = 1, double scalefactor = 1.00);
void MassPlot(int Irebin = 20);

//-- skeleton functions
void SideBandFit(int Irebin = 10);
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig);
TH1D * GenerateToyDataSet(TH1D *h_mass_template);
double ExpectedSignificance_ToyMC(double b = 1.0, double s = 10., double db = 0.0, double Ntoy = 1e6, int Iplot = 0);

//-- some fusefull functions
double IntegratePoissonFromRight(double mu, int N_obs);
double IntegrateFromRight(TH1D * h_X_bgr, double X_value);
vector<double> Get_Quantiles( TH1D* hist );
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045, 
	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );
//-------------------------------------------------------------------------------------------------------------------------



//========================================
// S O M E   F I N A L   F U N C T I O N S 
//========================================


//========================================================
TH1D * GetMassDistribution(int Itype, double scalefactor){
//========================================================
 //----------------------------------------------------------
 // Goal: return the histogram of the 4-lepton invariant mass
 //       for given type with an optional scale factor
 //
 //       Itype 1 = ZZ SM background
 //             2 = data
 //           125 = Higgs 125
 //           200 = Higgs 200
 //
 //      scalefactor: histograms will be scaled with this number
 //
 //  Note: Histograms have ~200 MeV bins, so need to rebin
 //---------------------------------------------------------

  //-- [1] Get histogram from the file
  TH1D *h_mass = 0;
  TDirectory* dir = gDirectory;   
  TFile *file = new TFile("Histograms_fake.root", "READ");
  dir->cd();

  //-- Higgs 125
  if(Itype == 125){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");     
  }
  //-- Higgs 200
  if(Itype == 200){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");     
  }
  //-- ZZ SM background
  if(Itype == 1){
    h_mass  = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");     
  }
  //-- data
  if(Itype == 2){
    h_mass  = (TH1D*) file->Get("h_m4l_data_fake")->Clone("h_mass");     
  }
  
  //-- [2] scale histograms
  int Nbins = h_mass->GetNbinsX();
  for (int i_bin = 1; i_bin < Nbins; i_bin++){
    double mu_bin = h_mass->GetBinContent(i_bin);
    h_mass -> SetBinContent( i_bin, scalefactor * mu_bin);
  }


  file->Close();
  //-- [3] return histogram   
  return h_mass;

  //===========================
} // end GetMassDistribution()
  //===========================




//========================
void MassPlot(int Irebin){
//========================
  // ------------------------------------------
  // Goal: produce SM+Higgs+data plot
  //       Note: rebinning is only for plotting
  // ------------------------------------------

  //------------------------------------
  //-- Standard stuff and prepare canvas
  //------------------------------------
  gROOT->Clear();
  gROOT->Delete();

  //-- Prepare canvas and plot histograms
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o Make cumulative histograms (for signal and background)
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125);
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  //-----------------------------------
  //-- [2] Plot histograms and make gif
  //--     o rebin histograms 
  //--     o prepare cumulative histogram
  //--     o make plot + opsmuk + gif
  //-----------------------------------

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);

  //-- Prepare cumulative histogram for signal + background 
  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset();
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
       printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
  }

  //-- prepare histograms and plot them on canvas
  double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
  double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
  h_sig_plus_bgr->SetFillColor(7); 
  h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
  h_sig_plus_bgr->SetAxisRange(0.,400.,"X");
  h_bgr->SetFillColor(2); 
  h_sig_plus_bgr->Draw("hist");  
  h_bgr->Draw("same");  
  h_bgr->Draw("axis same");  
  h_data->Draw("e same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-lepton invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
  leg1->Draw();

  //-- prepare gif
  canvas1->Print(Form("./MassPlot_rebin%d.gif",Irebin));

  return;

   //===============
 } // end MassPlot()
   //===============




//===============================================
// S O M E   S K E L E T O N    F U N C T I O N S 
//===============================================


//=============================================================
void Significance_Optimization(double Lumi_scalefactor = 1.00){
//=============================================================

  printf("\n Significance_Optimization()\n\n");

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o scale to correct luminosity
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
  h_sig  = GetMassDistribution(125);
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  
  // scale histograms to required luminosity
  TH1D *h_sig_plus_bgr = (TH1D*)h_bgr->Clone("h_sig_plus_bgr"); h_sig_plus_bgr->Reset();
  Int_t Nbins = h_bgr->GetNbinsX();
  for(int i_bin = 1; i_bin <= Nbins; i_bin++){
    h_sig->SetBinContent(i_bin, Lumi_scalefactor*h_sig->GetBinContent(i_bin));
    h_bgr->SetBinContent(i_bin, Lumi_scalefactor*h_bgr->GetBinContent(i_bin));
    h_sig_plus_bgr->SetBinContent(i_bin, h_sig->GetBinContent(i_bin)+h_bgr->GetBinContent(i_bin));
  }
  

  //-------------------------------------------------------
  //-- [2] Compute significance for various mass windows
  //--     o try various options for window (use histogram)
  //--     o compute expected and observed significance 
  //-------------------------------------------------------

  //-- Define histogram that defines mass windows to try
  TH1D *h_masswindow          = new TH1D("h_masswindow","",250,0.,25.);           // make a mass window - full width between 0 and 25 GeV
  TH1D *h_masswindow_expected = new TH1D("h_masswindow_expected","",250,0.,25.);  // histogram to hold results for expected
  TH1D *h_masswindow_observed = new TH1D("h_masswindow_observed","",250,0.,25.);  // histogram to hold results for observed

  // loop over various mass window options

  double mass_central = 125.;
  double masswindow_expected_optimal    = 0.00;
  double masswindow_observed_optimal    = 0.00;
  double significance_expected_optimal  = 0.00;
  double significance_observed_optimal  = 0.00;

  //---------------------------------
  //-- Loop over various mass windows
  //---------------------------------
  for (int i_bin = 1; i_bin<h_masswindow->GetNbinsX(); i_bin++ ){
  
    //-- get full width of mass window (bin center of our histogram) and the number of events in mass window for each event type
    double masswindow_fullwidth = h_masswindow->GetBinCenter(i_bin);    
    printf("   Trying as mass window: %5.2f GeV\n", masswindow_fullwidth);

    //-- [a] determine the number of events in the mass window for each event type
    //       Ndata_win, Nbgr_win and Nsig_win

    int Bin_low = h_data->FindBin(mass_central - 0.5*masswindow_fullwidth);
    int Bin_up = h_data->FindBin(mass_central + 0.5*masswindow_fullwidth);
    double Ndata_win = h_data->Integral(Bin_low,Bin_up);
    double Nbgr_win = h_bgr->Integral(Bin_low,Bin_up);
    double Nsig_win = h_sig->Integral(Bin_low,Bin_up);
    double Nsig_plus_bgr_win = Nsig_win + Nbgr_win;

    //-- compute EXPECTED significance and save in histogram
    double pvalue_expected       = IntegratePoissonFromRight(Nbgr_win,(int)Nsig_plus_bgr_win); // you need to do this yourself
    double significance_expected = 0.00;
    if( (1.000-pvalue_expected) < 1e-12 || pvalue_expected < 1e-12){
      printf(" undef expected pvalue -> put significance to zero\n");
    } else{
      significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    }
    // save optimal values
    if(significance_expected > significance_expected_optimal){
      significance_expected_optimal = significance_expected;
      masswindow_expected_optimal = masswindow_fullwidth;
    }
    h_masswindow_expected->SetBinContent(i_bin,significance_expected);

    //--  compute OBSERVED significance and save in histogram
    double pvalue_observed       = IntegratePoissonFromRight(Nbgr_win,(int)Ndata_win); // you need to do this yourself
    double significance_observed = 0.00;
    if( (1.000-pvalue_observed) < 1e-12 || pvalue_observed < 1e-12){
      printf("udef observed pvalue -> put sig to 0.00\n");
    } else {
      significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
    }
    // save optimal values
    if(significance_observed > significance_observed_optimal){
      significance_observed_optimal = significance_observed;
      masswindow_observed_optimal = masswindow_fullwidth;
    }
    h_masswindow_observed->SetBinContent(i_bin, significance_observed);


    //printing
    printf(" Window: (%5.2f) %4.1f < mh < %4.1f || b = %5.2f s = %5.2f | s+b = %5.2f d = %d \n",masswindow_fullwidth,mass_central - 0.5*masswindow_fullwidth, mass_central + 0.5*masswindow_fullwidth, Nbgr_win, Nsig_win,Nsig_plus_bgr_win, (int)Ndata_win);
    printf("  ||  siginf = %5.2f(%5.2f) sigma exp(obs) \n", significance_expected, significance_observed);

  } // end loop over width mass window

  //-- print optimum to the screen
  //
  printf("\n");
  printf("  --------------\n");
  printf("  Optimal mass window:\n");
  printf("  --------------\n");
  printf("Expected: optimal mass window (Lumi scale factor = %5.2f) = %5.2f GeV -> expected significance = %5.2f sigma \n",Lumi_scalefactor,masswindow_expected_optimal,significance_expected_optimal );
  printf("Observed: optimal mass window (Lumi scale factor = %5.2f) = %5.2f GeV -> observed significance = %5.2f sigma \n",Lumi_scalefactor,masswindow_observed_optimal, significance_observed_optimal );
  printf("\n");


  //----------------------------------
  //-- [3] Plot histogram and make gif
  //----------------------------------
  

  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 
 
  h_masswindow_expected->SetLineColor(1);
  h_masswindow_expected->SetLineWidth(2);
  h_masswindow_observed->SetLineColor(4);
  h_masswindow_observed->SetLineWidth(2);

  h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_masswindow_expected->Draw("l");
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    h_masswindow_observed->Draw("l same");
  }
  //-- axes
  AddText( 0.900, 0.035, "Mass window GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    
  AddText( 0.225, 0.825, Form("Luminosity scalefactor = %5.1f",Lumi_scalefactor),0.050, 0.,"left");            

  AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);                        
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);                        
  }
  //-- prepare gif
  canvas1->Print(Form("./Significance_Optimization_lumiscalefactor%d.gif",int(Lumi_scalefactor)));

  return;

  //================================
} // end Significance_Optimization()
  //================================






//===========================
void SideBandFit(int Irebin){
//===========================

  printf("\n SideBandFit()\n\n");

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr, *h_data;
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);
 
  //-- rebin histograms if necessary
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);
  printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.2f GeV\n", Irebin, h_data->GetBinWidth(1));


  // General variables
  double m4l_bin      = 0.;
  double Nobs_bin     = 0.;
  double mean_bgr_bin = 0.;
  int Nbins = h_data->GetNbinsX();

  double loglik = 0;

  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",10000,0.1,3.1); // you'll want to put more here I guess
  int Nbins_sf = h_scalefactor_bgr->GetNbinsX();
  double scalefactor_bgr  = 1.00; 
  for (int i_bin_sf = 1; i_bin_sf <= Nbins_sf; i_bin_sf++){
    loglik = 1e-10;
    scalefactor_bgr = h_scalefactor_bgr->GetBinCenter(i_bin_sf);


    //Likelihood
    for(int i_bin = 1; i_bin <=Nbins; i_bin++){
      m4l_bin = h_data->GetBinCenter(i_bin);
      Nobs_bin = h_data->GetBinContent(i_bin);
      if(m4l_bin >= 150. && m4l_bin<= 400. ){
        mean_bgr_bin = scalefactor_bgr * h_bgr->GetBinContent(i_bin);
        if(mean_bgr_bin > 0.){
          loglik += TMath::Log( TMath::Poisson( Nobs_bin,mean_bgr_bin));
        } else {
          printf("\n Warning: Mean expected bin content is zero for bin at mass %f\n", m4l_bin);
          printf("Likelihood will be corrupted\n");
        }
      }
    }

    h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik);
  }

  //---------------------------------------------
  // [2a] Loop 1: loop over scalefactors in alpha
  //---------------------------------------------
/*
  TH1D *h_scalefactor_bgr_rescaled
  for (int i_bin_sf = 1; i_bin_sf <=Nbins_sf; i_bin_sf++){

    //-- determine the scale factor for the background
    scalefactor_bgr = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
    printf(" Loop 1: I am now trying alpha = %5.2f\n",scalefactor_bgr);
  
    //-----------------------------------------------------------------------------------
    // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
    //-----------------------------------------------------------------------------------
    double loglik = 1e-10;
    for (int i_bin = 1; i_bin <=  h_data->GetNbinsX(); i_bin++){
      printf("        bin %d, m = %5.2f, Ndata = %5.2f, Nbgr = %5.2f\n",i_bin, h_bgr->GetBinCenter(i_bin),h_bgr->GetBinContent(i_bin),h_data->GetBinContent(i_bin));
      loglik += 1.; // you'll have to do this yourself;  	
    } // end loop over bins

    h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik);   
  } // end loop over scale factors for the background
  */

    //----------------------------------------------------
    //-- [3] Interpret the -2Log (Likelihood distribution)
    //----------------------------------------------------
    TH1D *h_scalefactor_bgr_rescaled  = (TH1D*) h_scalefactor_bgr->Clone("h_scalefactor_bgr_rescaled");  
    h_scalefactor_bgr_rescaled->Reset();
    double min2loglik = 0;
    double min2loglik_min = 1e90;
    double scalefactor_bgr_best = -999;   
    
    //-- Find minimum
    for (int i_bin_sf = 1; i_bin_sf <= Nbins_sf; i_bin_sf++){
      min2loglik = h_scalefactor_bgr->GetBinContent(i_bin_sf);
      if(min2loglik < min2loglik_min){
        min2loglik_min = min2loglik;
        scalefactor_bgr_best = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
      }
    }

    //--Rescale and find +/- 1 sigma errors
    double left_threshold = -1;
    double right_threshold = -1;
    for(int i_bin_sf = 1; i_bin_sf <= Nbins_sf; i_bin_sf++){
      h_scalefactor_bgr_rescaled->SetBinContent(i_bin_sf, h_scalefactor_bgr->GetBinContent(i_bin_sf) - min2loglik_min);
      if(left_threshold < 0. && h_scalefactor_bgr_rescaled->GetBinContent(i_bin_sf) < 1.000){
        left_threshold = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
      }
      if( left_threshold > 0. && right_threshold < 0. && h_scalefactor_bgr_rescaled->GetBinContent(i_bin_sf) > 1.000){
        right_threshold = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
      }
    }

    //-- print summary to screen
    double error_left = scalefactor_bgr_best - left_threshold;
    double error_right = right_threshold - scalefactor_bgr_best;

    printf("------------\n");
    printf("Result fit: \n");
    printf("------------\n");
    printf("Background scale factor from sideband fit = %5.2f -%5.2f +%5.2f\n",scalefactor_bgr_best, error_left, error_right );
    printf("\n");

  // Plot histgram and make gif
  TCanvas * canvas1 = new TCanvas("canvas1","Standard Canvas",600,400);
  canvas1->SetLeftMargin(0.175);
  canvas1->SetBottomMargin(0.125);
  canvas1->cd();
  h_scalefactor_bgr_rescaled->Draw();
  canvas1->Print("./SideBandFit.gif");

  
  double mass_central = 125.;
  double masswindow_fullwidth = 7.15;
  int Bin_low = h_data->FindBin(mass_central - 0.5*masswindow_fullwidth);
  int Bin_up = h_data->FindBin(mass_central + 0.5*masswindow_fullwidth);
  double Nbgr_unscaled = h_bgr->Integral(Bin_low,Bin_up);
  double Nbgr_scaled = scalefactor_bgr_best * Nbgr_unscaled;
  double Nbgr_scaled_lowererror = error_left * Nbgr_unscaled;
  double Nbgr_scaled_uppererror = error_right * Nbgr_unscaled;

  printf("------------------\n");
  printf("SM background in mass window: width = %5.2f\n", masswindow_fullwidth);
  printf("------------------\n");
  printf("Unscaled: Nbgr = %5.2f\n", Nbgr_unscaled);
  printf("Scaled: Nbgr = %5.2f -%5.2f +%5.2f\n", Nbgr_scaled, Nbgr_scaled_lowererror, Nbgr_scaled_uppererror);
  printf("\n");

  return;

  //==================
} // end SideBandFit()
  //==================



//=========================================================================================
double Get_TestStatistic(TH1D *h_mass_dataset = GenerateToyDataSet(GetMassDistribution(2,1.0)), TH1D *h_template_bgr = GetMassDistribution(1,1.0), TH1D *h_template_sig = GetMassDistribution(125,1.0)){
//=========================================================================================
  printf(" dummy = %d\n",h_mass_dataset->GetNbinsX() + h_template_bgr->GetNbinsX() + h_template_sig->GetNbinsX());
   
  int Nbins = h_mass_dataset->GetNbinsX();
  double massrange_fit_min = 100.;
  double massrange_fit_max = 400.;

  //-- do likelihood fit
  double loglik_bgr           = 0.;
  double loglik_sig_plus_bgr  = 0.;
  for (int i_bin = 1; i_bin <= Nbins; i_bin++){
    double mass_bin = h_mass_dataset->GetBinCenter(i_bin);
    double Nobs_bin = h_mass_dataset->GetBinContent(i_bin);
    double mu_bgr_bin = h_template_bgr->GetBinContent(i_bin);
    double mu_sig_bin = h_template_sig->GetBinContent(i_bin);
    if((mass_bin >= massrange_fit_min && mass_bin <= massrange_fit_max)){
      loglik_bgr +=           TMath::Log( TMath::Poisson(Nobs_bin, mu_bgr_bin + 0.0 * mu_sig_bin)); // likelihood for the mu=0 (b-only) scenario
      loglik_sig_plus_bgr +=  TMath::Log( TMath::Poisson(Nobs_bin, mu_bgr_bin + 1.0 * mu_sig_bin)); // likelihood for the mu=1 (s+b) scenario
    }
  } // end loop over bins
 
  //-- compute test statistic
  double test_statistic  =  -2. * (loglik_sig_plus_bgr - loglik_bgr);

  printf(" test_statistic = %5.2f\n", test_statistic);

  //-- return test_statistic
  return test_statistic;

} // end Get_TestStatistic()





//===============================================
TH1D * GenerateToyDataSet(TH1D *h_mass_template){
//===============================================
  //-------------------------------------------------
  // Goal: Generate Toy data set from input histogram
  // How:  dumb way -> draw random Poisson in each bin       
  //--------------------------------------------------
  TRandom3 *R = new TRandom3(0);

  //-- Create new histogram for the data-set
  TH1D *h_mass_toydataset = (TH1D*) h_mass_template->Clone("h_mass_toydataset"); h_mass_toydataset->Reset();
  int Nbins = h_mass_toydataset->GetNbinsX();

  //-- Loop over bins and draw Poisson number of event in each bin
  for(int i_bin = 1; i_bin <= Nbins; i_bin++){
    double mu_bin = h_mass_template->GetBinContent(i_bin);
    int Nevt_bin = R->Poisson(mu_bin);
    h_mass_toydataset->SetBinContent(i_bin, Nevt_bin);
  }

  //-- return histogram of toy data-set
  return h_mass_toydataset;
  
} // end GenerateToyDataSet()

/*
void Significance_LikelihoodRatio_ToyMC(int Ntoys = 100){
  TRandom3 *R = new TRandom3();

  // Get mass distribution for the background
  TH1D* h_bgr = GenerateToyDataSet(GetMassDistribution(1));

  // Histogram that contain significances around 125 and 00 GeV (and the maximum one)
  TH1D *h_significance_125 = new TH1D("h_significance_125", "significance 125", 1000., -6.0, 6.0);
  TH1D *h_significance_175 = new TH1D("h_significance_175", "significance 175", 1000., -6.0, 6.0);
  TH1D *h_significance_225 = new TH1D("h_significance_225", "significance 225", 1000., -6.0, 6.0);
  TH1D *h_significance_275 = new TH1D("h_significance_275", "significance 275", 1000., -6.0, 6.0);
  TH1D *h_significance_max = new TH1D("h_significance_max", "significance max", 1000., -6.0, 6.0);

  // 10 GeV mass window 
  double masswindow_fullwidth = 10.00;

  int Bin_low_125 = h_bgr->FindBin(125.0 - 0.5*masswindow_fullwidth);
  int Bin_up_125 = h_bgr->FindBin(125.0 + 0.5*masswindow_fullwidth);

  int Bin_low_175 = h_bgr->FindBin(175.0 - 0.5*masswindow_fullwidth);
  int Bin_up_175 = h_bgr->FindBin(175.0 + 0.5*masswindow_fullwidth);

  int Bin_low_225 = h_bgr->FindBin(225.0 - 0.5*masswindow_fullwidth);
  int Bin_up_225 = h_bgr->FindBin(225.0 + 0.5*masswindow_fullwidth);

  int Bin_low_275 = h_bgr->FindBin(275.0 - 0.5*masswindow_fullwidth);
  int Bin_up_275 = h_bgr->FindBin(275.0 + 0.5*masswindow_fullwidth);

  double mean_Nbgr_win_125 = h_bgr->Integral(Bin_low_125, Bin_up_125);
  double mean_Nbgr_win_175 = h_bgr->Integral(Bin_low_175, Bin_up_175);
  double mean_Nbgr_win_225 = h_bgr->Integral(Bin_low_225, Bin_up_225);
  double mean_Nbgr_win_275 = h_bgr->Integral(Bin_low_275, Bin_up_275);


  // Generate toy datasets

  // 125 GeV stuff
  double b_toy_125          = 0.;
  double significance_125 = 0.;
  double pvalue_125         = 0.;
  // 175 GeV stuff
  double b_toy_175          = 0.;
  double significance_175 = 0.;
  double pvalue_175         = 0.;
  // 225 GeV stuff
  double b_toy_225          = 0.;
  double significance_225 = 0.;
  double pvalue_225         = 0.;
  // 275 GeV stuff
  double b_toy_275          = 0.;
  double significance_275 = 0.;
  double pvalue_275         = 0.;

  // start loop over toys
  for(int i_toy = 0; i_toy < Ntoys; i_toy++){

    // toy at 125 GeV
    b_toy_125 = R->Poisson(mean_Nbgr_win_125);
    pvalue_125 = IntegratePoissonFromRight(mean_Nbgr_win_125, (int)b_toy_125);
    if( (1.000-pvalue_125) < 1e-12){ significance_125 = -5.;} else if ( pvalue_125 < 1e-12){significance_125 = 5.;}
    else { significance_125 = ROOT::Math::gaussian_quantile_c(pvalue_125,1);}
    h_significance_125->Fill(significance_125);

    // toy at 175 GeV
    b_toy_175 = R->Poisson(mean_Nbgr_win_175);
    pvalue_175 = IntegratePoissonFromRight(mean_Nbgr_win_175, (int)b_toy_175);
    if( (1.000-pvalue_175) < 1e-12){ significance_175 = -5.;} else if ( pvalue_175 < 1e-12){significance_175 = 5.;}
    else { significance_175 = ROOT::Math::gaussian_quantile_c(pvalue_175,1);}
    h_significance_175->Fill(significance_175);

    // toy at 225 GeV
    b_toy_225 = R->Poisson(mean_Nbgr_win_225);
    pvalue_225 = IntegratePoissonFromRight(mean_Nbgr_win_225, (int)b_toy_225);
    if( (1.000-pvalue_225) < 1e-12){ significance_225 = -5.;} else if ( pvalue_225 < 1e-12){significance_225 = 5.;}
    else { significance_225 = ROOT::Math::gaussian_quantile_c(pvalue_225,1);}
    h_significance_225->Fill(significance_225);

    // toy at 275 GeV
    b_toy_275 = R->Poisson(mean_Nbgr_win_275);
    pvalue_275 = IntegratePoissonFromRight(mean_Nbgr_win_275, (int)b_toy_275);
    if( (1.000-pvalue_275) < 1e-12){ significance_275 = -5.;} else if ( pvalue_275 < 1e-12){significance_275 = 5.;}
    else { significance_275 = ROOT::Math::gaussian_quantile_c(pvalue_275,1);}
    h_significance_275->Fill(significance_275);

    // Find Maximum significance of 125 and 200

    double sig_max1 = (significance_125 > significance_175) ? significance_125 : significance_175;
    double sig_max2 = (significance_225 > significance_275) ? significance_225 : significance_275;
    double significance_max = (sig_max1 > sig_max2) ? sig_max1 : sig_max2;

    h_significance_max->Fill(significance_max);

  }

  // Prepare plots
  h_significance_125->SetLineColor(1);
  h_significance_175->SetLineColor(2);
  h_significance_225->SetLineColor(3);
  h_significance_275->SetLineColor(4);
  h_significance_max->SetLineColor(8);

  h_significance_125->Draw();
  h_significance_175->Draw("same");
  h_significance_225->Draw("same");
  h_significance_275->Draw("same");
  h_significance_max->Draw("same");

  // compute what fraction of toys have a significance above 2 sigma

  printf("\n");
  printf(" Fraction above 2 sigma (125) = %6.4f\n", IntegrateFromRight(h_significance_125, 2.00));
  printf(" Fraction above 2 sigma (175) = %6.4f\n", IntegrateFromRight(h_significance_175, 2.00));
  printf(" Fraction above 2 sigma (225) = %6.4f\n", IntegrateFromRight(h_significance_225, 2.00));
  printf(" Fraction above 2 sigma (275) = %6.4f\n", IntegrateFromRight(h_significance_275, 2.00));
  printf(" Fraction above 2 sigma (max) = %6.4f\n", IntegrateFromRight(h_significance_max, 2.00));


  // Find translation

  TH1D *h_significance_tranlation = (TH1D*)h_significance_max->Clone("H_significance_tranlation"); h_significance_tranlation->Reset();
  int Nbins = h_significance_max->GetNbinsX();
  double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.);
  double new_2sigmapoint = -999.;
  for(int i_bin = 1; i_bin < Nbins; i_bin++){
    double gauss_nsigma_original = h_significance_max->GetBinCenter(i_bin);
    double pvalue_new = IntegrateFromRight(h_significance_max, gauss_nsigma_original);
    double gauss_nsigma_new = (pvalue_new < 1e-12 || pvalue_new > 1.-1e-12) ? 0. : ROOT::Math::gaussian_quantile_c(pvalue_new, 1);
    h_significance_tranlation->SetBinContent(i_bin, gauss_nsigma_new);

    if(new_2sigmapoint < 0. && pvalue_new < fraction_2sigma){
      new_2sigmapoint = gauss_nsigma_original;
    }
  }

  // option 2

  double pvalue_new_2 = IntegrateFromRight(h_significance_max, 2.00);
  double gauss_nsigma_new_2 = (pvalue_new_2 < 1e-12 || pvalue_new_2 > 1-1e-12) ? 0. : ROOT::Math::gaussian_quantile_c(pvalue_new_2, 1);


  TCanvas * canvas1 = new TCanvas( "canvas1", "Standard Canvas", 600, 400);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
  canvas1->cd();
  h_significance_tranlation->SetAxisRange(-2.,3.,"X");
  h_significance_tranlation->Draw("chist");
  AddText(0.900,0.035,"Local significance",0.060,0.,"right");
  AddText(0.040,0.900,"Global significance",0.060,90.,"right");

  canvas1->Print("Significance_LikelihoodRatio_ToyMC.gif");

  return;

}
*/

void Significance_LikelihoodRatio_ToyMC(int Ntoys = 100., int = 2.){
  for(int toy = 1; toy < Ntoys; toy++){
    TH1D* h_dataset = GenerateToyDataSet(GetMassDistribution(2));
    TH1D* h_bgr = GenerateToyDataSet(GetMassDistribution(1));
    TH1D* h_sig_plus_bgr = GenerateToyDataSet(GetMassDistribution(125));

    Get_TestStatistic(h_dataset, h_bgr, h_sig_plus_bgr);
  }
  // Use Test Statistic to get a distribution for b-only and s+b for 10000 experiments. Call TS, calling Generate Toy Dataset inside. Use templates. 
}

//============================================
double ExpectedSignificance_ToyMC(double b = 1.0, double s = 10., double db = 0.0, double Ntoy = 1e6, int Iplot = 0){
//============================================  



  // Standard stuff
  gROOT->Clear();
  gROOT->Delete();


  // Define histograms for number of events
  TH1D* h_Nbgr          = new TH1D("h_Nbgr",          "Number of background events        ", 500, -0.5, 499.5);
  TH1D* h_Nsig_plus_bgr = new TH1D("h_Nsig_plus_bgr", "Number of signal + background events", 500, -0.5, 499.5);

  // Toy experiments
  TRandom3 *R = new TRandom3();


  for(int i_toy = 0; i_toy < Ntoy; i_toy++){

    // Determine for this experiment the mean of the Poissons for b-only and s+b
    double mean_b = R->Gaus(b, db);
    double mean_splusb = R->Gaus(b, db) + s;

    // compute number of observed events in this experiment
    if(mean_b >= 0. && mean_splusb >= 0.){ // only throw Poisson if mean is positive
      double b_toy = R->Poisson(mean_b);
      double splusb_toy = R->Poisson(mean_splusb);

      h_Nbgr->Fill(b_toy);
      h_Nsig_plus_bgr->Fill(splusb_toy);
    }
    else{
      printf("the mean of b and s+b is less than zero, disregaring toy\n");
    }
    if ((i_toy%100000) == 0 && Iplot){ cout << "ExpectedSignificance(): Toy " << i_toy << " out of " << Ntoy << endl;}
  }


  // Compute p-values for median signal + background experiment

  double Nevts_median_splusb =s+b;
  double pvalue = IntegrateFromRight(h_Nbgr,Nevts_median_splusb); // Integrate #events for background from median
  double gauss_nsigma = ROOT::Math::gaussian_quantile_c(pvalue, 1);

  printf("p_value = %5.2e  ->   %5.2f sigma \n", pvalue, gauss_nsigma);

  if(Iplot == 1){
    TCanvas * canvas1 = new TCanvas("canvas1", "Standard Canvas", 600, 400);
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125);


    // Create separate histogram with excess
    TH1D* h_Nbgr_excess = (TH1D*)h_Nbgr->Clone("h_Nbgr_excess"); h_Nbgr_excess->Reset();
    for(int i_bin = 1; i_bin < h_Nbgr->GetNbinsX(); i_bin++){
      if(h_Nbgr->GetBinCenter(i_bin) >= (int)Nevts_median_splusb){
        h_Nbgr_excess->SetBinContent(i_bin,h_Nbgr->GetBinContent(i_bin));
      }
    }


    double Ymax_plot = 1.30 * h_Nbgr->GetBinContent(h_Nbgr->GetMaximumBin());
    h_Nbgr->SetAxisRange(0.,(s+b+15.),"X");
    h_Nbgr->SetAxisRange(0.,Ymax_plot,"Y");
    h_Nbgr->Draw("hist");
    h_Nbgr_excess->SetFillColor(2);
    h_Nbgr_excess->Draw("hist same");
    h_Nsig_plus_bgr->SetLineStyle(2);
    h_Nsig_plus_bgr->Draw("same");

    // Add axis labels

    AddText(0.900, 0.035, "Number of events", 0.060,0.,"right"); // X-axis
    AddText(0.040, 0.900, "Number of toy-experiments", 0.060, 90.,"right"); // Y-axis


    // Add Significance
    AddText(0.225, 0.825, "background-only", 0.050, 0., "left");
    AddText(0.225, 0.775, "experiments", 0.050, 0., "left");

    AddText(0.550, 0.825, "Signal+background", 0.050, 0.,"left");
    AddText(0.550, 0.775, "experiments", 0.050, 0.,"left");

    AddText(0.625, 0.675, Form("p-value = %7.2e", pvalue), 0.050, 0., "left");
    AddText(0.625, 0.625, Form("signif. = %5.2f sigma", gauss_nsigma), 0.050, 0.,"left");

    canvas1->Print("ExpectedSignificance_ToyMC.gif");


  }

  return pvalue;
}











//================================================
// S O M E   U S E F U L L   F U N C T I O N S 
//================================================



//====================================================
double IntegratePoissonFromRight(double mu, int N_obs){
//====================================================
// --------------------------------------------------------------------
// Compute p-value for case zero background uncertainty, i.e. 
//         just integrate Poisson from the right from N_obs to infinity
// --------------------------------------------------------------------

  double integral = 1.; 
  for(int i_obs = 0; i_obs < N_obs; i_obs++){
    integral -= TMath::Poisson(i_obs,mu);
  }
  
  return integral;
  
} // end IntegratePoissonFromRight()


//========================================================
double IntegrateFromRight(TH1D * h_X_bgr, double X_value){
//========================================================
// --------------------------------------------------------------------
// Compute p-value: integrate number of events from X_value to infinity
// --------------------------------------------------------------------

  //-- Integrate distributions
  int Nbins = h_X_bgr->GetNbinsX();
  int X_bin = h_X_bgr->FindBin(X_value); 

  //-- Compute integral from X-value to infinity
  double pvalue = h_X_bgr->Integral(X_bin,Nbins) / h_X_bgr->Integral();
  
  return pvalue;

} // end IntegrateFrom Right()




//=========================================
vector<double> Get_Quantiles( TH1D* hist ){
//=========================================
// Quantiles returns a vector<double> with 5 entries.
// Entries 0 and 4 are the values on the histogram x-axis
// so that 95% of the content lies between these values.
// Entries 1 and 3 bracket 68% in the same way.
// Entry 2 is the median of the histogram.

  //-- define quantiles
  double fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // 15.8655 %
  double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.); //  2.2750 %
  double probs[5] = {fraction_2sigma, fraction_1sigma, 0.50, 1.00-fraction_1sigma, 1.00-fraction_2sigma };

  //-- output of the quantiles
  double Xvalues[5];

  //-- extract quantiles
  hist->GetQuantiles( 5, Xvalues, probs );
  
  vector<double> Xvalues_output(5);
  for (int i=0; i<5; i++) 
    {
      Xvalues_output[i] = Xvalues[i];
    }

  return Xvalues_output;
} // end Get_Quantiles()





//=======================================================================================================================
void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
              Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
//=======================================================================================================================
{
  Int_t txt_align = 12;
  if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left 
  if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
  if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

  TLatex* t1 = new TLatex( txt_x, txt_y, txt);
  if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
  t1->SetTextSize(txt_size);
  t1->SetTextAlign(txt_align);
  t1->SetTextAngle(txt_angle);
  t1->SetTextColor(txt_color);
  t1->Draw();

} // end AddText()

