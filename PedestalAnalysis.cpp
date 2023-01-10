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

using namespace std;
using namespace std::chrono;

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

  Char_t *inputfilename = (Char_t*)"";
  if (argc < 2)
  {
    cout << endl << "No input file name specified!  Aborting." << endl << endl;
    return -1;
  }
  else
  {
    inputfilename = argv[1];
  }

  ///////////////////////////////////
  // Set Up Pedestal Mean/RMS Maps
  ///////////////////////////////////
  
  double pedMean[17][130][64] = {0.0};
  double pedRMS[17][130][64] = {0.0};
  double pedNum[17][130][64] = {0.0};

  ///////////////////////////////////
  // Load Input Data
  ///////////////////////////////////

  TFile* inputfile = new TFile(inputfilename,"READ");
  
  TTreeReader reader("tree", inputfile);
  TTreeReaderValue<unsigned char> adc_counts(reader, "adc_counts");
  TTreeReaderValue<unsigned char> tile_id(reader, "tile_id");
  TTreeReaderValue<unsigned char> chip_id(reader, "chip_id");
  TTreeReaderValue<unsigned char> channel_id(reader, "channel_id");
  TTreeReaderValue<unsigned long long> timestamp(reader, "timestamp");

  ///////////////////////////////////
  // Set Up Output File
  ///////////////////////////////////

  ofstream outfile;
  outfile.open("pedestals_multitile.dat"); 
  
  ///////////////////////////////////
  // Loop Over Data
  ///////////////////////////////////

  while(reader.Next())
  {
    if((*tile_id < 1) || (*tile_id > 16) || (*chip_id < 11) || (*chip_id > 110) || (*channel_id < 0) || (*channel_id > 63))
    {
      continue;
    }

    pedMean[*tile_id][*chip_id][*channel_id] += (double) *adc_counts;
    pedRMS[*tile_id][*chip_id][*channel_id] += (double) pow(*adc_counts,2.0);
    pedNum[*tile_id][*chip_id][*channel_id] += 1.0;
  }

  for(int i = 0; i < 17; i++)
  {
    for(int j = 0; j < 130; j++)
    {
      for(int k = 0; k < 64; k++)
      {
        if(pedNum[i][j][k] > 0)
	{
          pedMean[i][j][k] /= pedNum[i][j][k];
          pedRMS[i][j][k] = sqrt(pedRMS[i][j][k]/pedNum[i][j][k] - pow(pedMean[i][j][k],2.0));
	  outfile << i << " " << j << " " << k << " " << pedMean[i][j][k] << " " << pedRMS[i][j][k] << endl;
	}
      }
    }
  }

  ///////////////////////////////////
  // Close Output File
  ///////////////////////////////////

  outfile.close();

  return 0;
}
