//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <numeric>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVirtualFFT.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TF1.h"

using namespace std;
using namespace std::chrono;

//const int numBins = 15;
const int numBins = 30;
const float driftVelE500 = 0.155; // Bern SingleCube and Module 0 Test
//const float driftVelE500 = 0.155*(186.0/201.0); // Module 0 HV Test
const float driftTimeMaxE500 = 186.0; // Bern SingleCube and Module 0 Test
const float driftTimeRangeE500 = 8.0; // Bern SingleCube and Module 0 Test
//const float driftTimeMaxE500 = 201.0; // Module 0 HV Test
//const float driftTimeRangeE500 = 5.0; // Module 0 HV Test
const float gain = 250.0*3.9;
const float engConv = 0.0000236;
//const float pedestal = 78.0; // 77 for CSU SingleCube?
//const float pedestal = 78.0;
const float pedestal = 0.0; // if pedestal already corrected in TrackMaker

const float dEdx = 2.65;
//const float LAr_density = 1.3692; // 91.5 K (Module 0 HV Test)
const float LAr_density = 1.3942; // 87.5 K (Module 0 Test)
//const float LAr_density = 1.3849; // 89 K
//const float LAr_density = 1.4095; // 85 K (Bern SingleCube)
const float ModBoxA = 0.930;
const float ModBoxB = 0.212;
const float A_Birks = 0.800;
const float k_Birks = 0.0486;

float FindVecMedian(vector<float> inputvec);

int main(int argc, char **argv)
{
  ///////////////////////////////////
  // Set Plot Formatting Options
  ///////////////////////////////////

  gErrorIgnoreLevel = kError;
  double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);
  gStyle->SetOptStat(0);

  ///////////////////////////////////
  // Get Input File Name
  ///////////////////////////////////

  char *inputfilename = (char*) "";
  float Efield = 500.0;
  if (argc < 2)
  {
    cout << endl << "No input file name specified!  Aborting." << endl << endl;
    return -1;
  }
  else
  {
    inputfilename = argv[1];
    if (argc >= 3)
    {
      Efield = (float) atof(argv[2]);
    }
  }
  
  ///////////////////////////////////
  // Set Correct Drift Velocity/Time
  ///////////////////////////////////

  TF1 driftVelFit("driftVelFit","pol5",-0.05,1.1);
  driftVelFit.SetParameter(0,0.0);
  driftVelFit.SetParameter(1,5.53416);
  driftVelFit.SetParameter(2,-6.53093);
  driftVelFit.SetParameter(3,3.20752);
  driftVelFit.SetParameter(4,0.389696);
  driftVelFit.SetParameter(5,-0.556184);

  float corrFactor = driftVelFit.Eval(Efield/1000.0)/driftVelFit.Eval(0.5);
  if(Efield > 500)
  {
    corrFactor *= 1.0-((141.5-140.325)/141.5)*((Efield-500.0)/(1000.0-500.0));
  }
  else if(Efield < 500)
  {
    corrFactor *= 1.0+((330.0-320.0)/320.0)*((500.0-Efield)/(500.0-200.0));
  }
  
  const float driftVel = driftVelE500*corrFactor;
  const float driftTimeMax = driftTimeMaxE500/corrFactor;
  const float driftTimeRange = driftTimeRangeE500/corrFactor;
  
  ///////////////////////////////////
  // Define Histograms
  ///////////////////////////////////

  TFile outfile("results.root","RECREATE");
  outfile.cd();

  TH1F *PitchHist = new TH1F("PitchHist","",50,0.0,5.0);
  
  TH2F *LifetimeHist2D = new TH2F("LifetimeHist2D","",numBins,0.0,driftTimeMax,50,0.0,160.0);
  TH1F *LifetimeHist2D_ProjX = (TH1F*) LifetimeHist2D->ProjectionX();
  TProfile *LifetimeHist2D_ProfileX;
  TH1F *LifetimeHist1D = new TH1F("LifetimeHist1D","",numBins,0.0,driftTimeMax);
    
  TH1F *dEdxHist = new TH1F("dEdxHist","",50,0.0,6.0);
  TH2F *dEdxHist2D = new TH2F("dEdxHist2D","",numBins,0.0,driftTimeMax,50,0.0,6.0);
  TH1F *dEdxHist2D_ProjX = (TH1F*) dEdxHist2D->ProjectionX();
    
  ///////////////////////////////////
  // Load Input Data
  ///////////////////////////////////

  TFile* inputfile = new TFile(inputfilename,"READ");
  
  TTreeReader readerTracks("trackTree", inputfile);
  TTreeReaderValue<int> eventNum(readerTracks, "eventNum");
  TTreeReaderValue<int> trackNum(readerTracks, "trackNum");
  TTreeReaderValue<float> minT_X(readerTracks, "minT_X");
  TTreeReaderValue<float> minT_Y(readerTracks, "minT_Y");
  TTreeReaderValue<float> minT_T(readerTracks, "minT_T");
  TTreeReaderValue<float> maxT_X(readerTracks, "maxT_X");
  TTreeReaderValue<float> maxT_Y(readerTracks, "maxT_Y");
  TTreeReaderValue<float> maxT_T(readerTracks, "maxT_T");
  TTreeReaderArray<float> trackHitX(readerTracks, "trackHitX");
  TTreeReaderArray<float> trackHitY(readerTracks, "trackHitY");
  TTreeReaderArray<float> trackHitT(readerTracks, "trackHitT");
  TTreeReaderArray<float> trackHitC(readerTracks, "trackHitC");

  ///////////////////////////////////
  // First Loop Over Data
  ///////////////////////////////////

  vector<float> dQdxVec[numBins];
  
  while (readerTracks.Next())
  {
    if((*maxT_T - *minT_T < 10.0*(driftTimeMax-driftTimeRange)) || (*maxT_T - *minT_T > 10.0*(driftTimeMax+driftTimeRange)))
    {
      continue;
    }

    float charge[numBins] = {0.0};
	
    int numHits = trackHitT.GetSize();
    for(int k = 0; k < numHits; k++)
    {
      int index = round(((((float) numBins)/driftTimeMax/10.0)*(trackHitT[k]-*minT_T))-0.5);
      if(index < 0)
      {
        index = 0;
      }
      else if(index > numBins-1)
      {
        index = numBins-1;
      }
      charge[index] += trackHitC[k]-pedestal;
    }

    float pitch = 0.1*(10.0*(dEdxHist2D_ProjX->GetBinCenter(2)-dEdxHist2D_ProjX->GetBinCenter(1))/(*maxT_T-*minT_T))*sqrt(pow(*maxT_X-*minT_X,2)+pow(*maxT_Y-*minT_Y,2)+pow(driftVel*(*maxT_T-*minT_T),2));
    PitchHist->Fill(pitch);
    
    for(int i = 1; i < numBins-1; i++)
    {
      LifetimeHist2D->Fill(LifetimeHist2D_ProjX->GetBinCenter(i+1),gain*charge[i]/pitch/1000.0);
      dQdxVec[i].push_back(gain*charge[i]/pitch/1000.0);
    }
  }

  LifetimeHist2D_ProfileX = LifetimeHist2D->ProfileX();
  for(int i = 0; i < numBins; i++)
  {
    float result = FindVecMedian(dQdxVec[i]);
    if(result > 0.0)
    {
      LifetimeHist1D->SetBinContent(i+1,result);
      LifetimeHist1D->SetBinError(i+1,sqrt(3.14159/2.0)*LifetimeHist2D_ProfileX->GetBinError(i+1));
    }
  }
  
  ///////////////////////////////////
  // Extract Electron Lifetime
  ///////////////////////////////////

  outfile.cd();
  TF1* lifetime_fit = new TF1("lifetime_fit","[0]*exp(-x/[1])");
  lifetime_fit->SetParameters(80.0,1000.0);
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",0,driftTimeMax);
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(5),LifetimeHist1D->GetBinLowEdge(13));
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(7),LifetimeHist1D->GetBinLowEdge(15));
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(6),LifetimeHist1D->GetBinLowEdge(14));
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(7),LifetimeHist1D->GetBinLowEdge(14)); // BEST SO FAR WITH NUMBINS=15
  TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(14),LifetimeHist1D->GetBinLowEdge(28)); // BEST SO FAR WITH NUMBINS=30 (USE FOR PAPER)
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(14),LifetimeHist1D->GetBinLowEdge(27)); // SYST 1
  //TFitResultPtr r = LifetimeHist1D->Fit("lifetime_fit","MQSE","",LifetimeHist1D->GetBinLowEdge(15),LifetimeHist1D->GetBinLowEdge(28)); // SYST 2

  float lifetime_val = r->Parameter(1)/1000.0;
  float lifetime_uncert_pos = fabs(r->UpperError(1)/1000);
  float lifetime_uncert_neg = fabs(r->LowerError(1)/1000);

  //cout << "Electron Lifetime:  " << lifetime_val << " + " << lifetime_uncert_pos << " - " << lifetime_uncert_neg << " ms" << endl;

  string temp = (string) inputfilename;
  temp = temp.substr(temp.find("_")+1);
  string year = temp.substr(0,temp.find("_")); 
  temp = temp.substr(temp.find("_")+1);
  string month = temp.substr(0,temp.find("_")); 
  temp = temp.substr(temp.find("_")+1);
  string day = temp.substr(0,temp.find("_")); 
  temp = temp.substr(temp.find("_")+1);
  string hour = temp.substr(0,temp.find("_")); 
  temp = temp.substr(temp.find("_")+1);
  string minute = temp.substr(0,temp.find("_")); 
  temp = temp.substr(temp.find("_")+1);
  string second = temp.substr(0,temp.find("_")); 

  cout << month << "-" << day << "-" << year << " " << hour << ":" << minute << ":" << second << "," << lifetime_val << "," << lifetime_uncert_pos << "," << lifetime_uncert_neg << "," << Efield << " V/cm" << endl;

  ///////////////////////////////////
  // Second Loop Over Data
  ///////////////////////////////////

  //float recombCorr = log(ModBoxA+(ModBoxB*dEdx)/(LAr_density*Efield/1000.0))/((ModBoxB*dEdx)/(LAr_density*Efield/1000.0)); // Modified Box Model
  float recombCorr = A_Birks/(1.0+(k_Birks*dEdx)/(LAr_density*Efield/1000.0)); // ICARUS Birks Model

  //cout << endl << recombCorr << endl;
  
  readerTracks.Restart();  
  while (readerTracks.Next())
  {
    if((*maxT_T - *minT_T < 10.0*(driftTimeMax-driftTimeRange)) || (*maxT_T - *minT_T > 10.0*(driftTimeMax+driftTimeRange)))
    {
      continue;
    }

    float charge[numBins] = {0.0};
	
    int numHits = trackHitT.GetSize();
    for(int k = 0; k < numHits; k++)
    {
      int index = round(((((float) numBins)/driftTimeMax/10.0)*(trackHitT[k]-*minT_T))-0.5);
      if(index < 0)
      {
        index = 0;
      }
      else if(index > numBins-1)
      {
        index = numBins-1;
      }
      charge[index] += trackHitC[k]-pedestal;
    }

    float pitch = 0.1*(10.0*(dEdxHist2D_ProjX->GetBinCenter(2)-dEdxHist2D_ProjX->GetBinCenter(1))/(*maxT_T-*minT_T))*sqrt(pow(*maxT_X-*minT_X,2)+pow(*maxT_Y-*minT_Y,2)+pow(driftVel*(*maxT_T-*minT_T),2));
    
    for(int i = 0; i < numBins; i++)
    {
      float lifetimeCorr = exp(dEdxHist2D_ProjX->GetBinCenter(i+1)/(1000.0*lifetime_val));
      dEdxHist2D->Fill(dEdxHist2D_ProjX->GetBinCenter(i+1),engConv*gain*lifetimeCorr*charge[i]/recombCorr/pitch);
      if((i > 0) && (i < numBins-1))
      {
        dEdxHist->Fill(engConv*gain*lifetimeCorr*charge[i]/recombCorr/pitch);
      }
    }
  }

  inputfile->Close();
  
  ///////////////////////////////////
  // Make Plots
  ///////////////////////////////////

  for(int i = 2; i <= LifetimeHist2D->GetNbinsX()-1; i++)
  {
    for(int j = 1; j <= LifetimeHist2D->GetNbinsY(); j++)
    {
      if(LifetimeHist2D->GetBinContent(i,j) <= 0.0)
      {
        LifetimeHist2D->SetBinContent(i,j,0.001);
      }
    }
  }

  for(int i = 1; i <= dEdxHist2D->GetNbinsX(); i++)
  {
    for(int j = 1; j <= dEdxHist2D->GetNbinsY(); j++)
    {
      if(dEdxHist2D->GetBinContent(i,j) <= 0.0)
      {
        dEdxHist2D->SetBinContent(i,j,0.001);
      }
    }
  }

  TCanvas c1;
  c1.cd();
  LifetimeHist2D->Draw("COLZ");
  LifetimeHist2D->SetTitle("");
  LifetimeHist2D->GetXaxis()->SetTitle("Drift Time [#mus]");
  LifetimeHist2D->GetXaxis()->SetTitleSize(0.045);
  LifetimeHist2D->GetXaxis()->SetTitleOffset(1.05);
  LifetimeHist2D->GetXaxis()->SetLabelSize(0.04);
  LifetimeHist2D->GetYaxis()->SetTitle("dQ/dx [ke#lower[-0.5]{-}/cm]");
  LifetimeHist2D->GetYaxis()->SetTitleSize(0.045);
  LifetimeHist2D->GetYaxis()->SetTitleOffset(1.12);
  LifetimeHist2D->GetYaxis()->SetLabelSize(0.04);
  c1.SaveAs("LifetimeHist2D.png");
  TPaletteAxis *palette_LifetimeHist2D = (TPaletteAxis*) LifetimeHist2D->GetListOfFunctions()->FindObject("palette");
  palette_LifetimeHist2D->SetX1NDC(0.885);
  palette_LifetimeHist2D->SetX2NDC(0.92);
  palette_LifetimeHist2D->SetY1NDC(0.10);
  palette_LifetimeHist2D->SetY2NDC(0.95);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.13);
  gPad->Modified();
  gPad->Update();
  c1.SaveAs("LifetimeHist2D.png");

  TCanvas c2;
  c2.cd();
  LifetimeHist1D->Draw();
  LifetimeHist1D->SetMarkerColor(kBlack);
  LifetimeHist1D->SetLineColor(kBlack);
  LifetimeHist1D->SetLineWidth(3.0);
  LifetimeHist1D->SetTitle("");
  LifetimeHist1D->GetXaxis()->SetTitle("Drift Time [#mus]");
  LifetimeHist1D->GetXaxis()->SetTitleSize(0.045);
  LifetimeHist1D->GetXaxis()->SetTitleOffset(1.05);
  LifetimeHist1D->GetXaxis()->SetLabelSize(0.04);
  LifetimeHist1D->GetYaxis()->SetTitle("Mean dQ/dx [ke#lower[-0.5]{-}/cm]");
  LifetimeHist1D->GetYaxis()->SetTitleSize(0.045);
  LifetimeHist1D->GetYaxis()->SetTitleOffset(1.4);
  LifetimeHist1D->GetYaxis()->SetLabelSize(0.04);
  LifetimeHist1D->GetYaxis()->SetRangeUser(0.0,160.0);
  TPaveText* text = new TPaveText(0.45,0.7,0.85,0.85,"nbNDC");
  text->AddText(Form("Elec. Lifetime:  %.2f^{+%.2f}_{-%.2f} ms",lifetime_val,lifetime_uncert_pos,lifetime_uncert_neg));
  text->SetFillColor(kWhite);
  text->SetTextSize(0.04);
  text->Draw("SAME");
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  c2.SaveAs("LifetimeHist1D.png");

  TCanvas c3;
  c3.cd();
  dEdxHist->Draw("HIST");
  dEdxHist->SetLineColor(kBlue);
  dEdxHist->SetLineWidth(3.0);
  dEdxHist->SetTitle("");
  dEdxHist->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  dEdxHist->GetXaxis()->SetTitleSize(0.045);
  dEdxHist->GetXaxis()->SetTitleOffset(1.05);
  dEdxHist->GetXaxis()->SetLabelSize(0.04);
  dEdxHist->GetYaxis()->SetTitle("# of Entries");
  dEdxHist->GetYaxis()->SetTitleSize(0.045);
  dEdxHist->GetYaxis()->SetTitleOffset(1.4);
  dEdxHist->GetYaxis()->SetLabelSize(0.04);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  c3.SaveAs("dEdxHist.png");

  TCanvas c4;
  c4.cd();
  dEdxHist2D->Draw("COLZ");
  dEdxHist2D->SetTitle("");
  dEdxHist2D->GetXaxis()->SetTitle("Drift Time [#mus]");
  dEdxHist2D->GetXaxis()->SetTitleSize(0.045);
  dEdxHist2D->GetXaxis()->SetTitleOffset(1.05);
  dEdxHist2D->GetXaxis()->SetLabelSize(0.04);
  dEdxHist2D->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  dEdxHist2D->GetYaxis()->SetTitleSize(0.045);
  dEdxHist2D->GetYaxis()->SetTitleOffset(1.12);
  dEdxHist2D->GetYaxis()->SetLabelSize(0.04);
  c4.SaveAs("dEdxHist2D.png");
  TPaletteAxis *palette_dEdxHist2D = (TPaletteAxis*) dEdxHist2D->GetListOfFunctions()->FindObject("palette");
  palette_dEdxHist2D->SetX1NDC(0.885);
  palette_dEdxHist2D->SetX2NDC(0.92);
  palette_dEdxHist2D->SetY1NDC(0.10);
  palette_dEdxHist2D->SetY2NDC(0.95);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.13);
  gPad->Modified();
  gPad->Update();
  c4.SaveAs("dEdxHist2D.png");
  
  ///////////////////////////////////
  // Write Output File
  ///////////////////////////////////
  
  outfile.cd();

  PitchHist->Write();
  LifetimeHist1D->Write();
  LifetimeHist2D->Write();
  dEdxHist->Write();
  dEdxHist2D->Write();
  
  outfile.Close();

  return 0;
}

float FindVecMedian(vector<float> inputvec)
{
  int size = inputvec.size();

  float result;
  if(size == 0)
  {
    result = 0.0;
  }
  else if(size == 1)
  {
    result = inputvec[0];
  }
  else
  {
    sort(inputvec.begin(), inputvec.end());
    if(size % 2 == 0)
    {
      result = (inputvec[size / 2 - 1] + inputvec[size / 2]) / 2.0;
    }
    else
    {
      result = inputvec[size / 2];
    }
  }

  return result;
}
