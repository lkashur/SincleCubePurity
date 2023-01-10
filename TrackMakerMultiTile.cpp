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
#include "TStopwatch.h"

//special includes
#include "dbscan.h"

using namespace std;
using namespace std::chrono;

const float timeSF = 6.5;
//const float scanRadius = 10.0;
const float scanRadius = 25.0;
const int minClusterSize = 50;

//#const int dataChunkSize = 38400;
const int dataChunkSize = 20000;
//const int dataChunkSize = 5000;
const int maxEvents = 10000000;

int main(int argc, char **argv)
{
  TStopwatch timer;
  timer.Start();
  
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
  // Set Up Pedestal Map
  ///////////////////////////////////
  
  float pedMean[17][130][64] = {0.0};

  ifstream pedmapfile;
  pedmapfile.open("pedestals_multitile.dat");

  string string_tile;
  string string_chip;
  string string_channel;
  string string_pedMean;
  string string_pedRMS;

  int input_tile;
  int input_chip;
  int input_channel;
  float input_pedMean;
  float input_pedRMS;

  while(getline(pedmapfile,string_tile,' '))
  {
    getline(pedmapfile,string_chip,' ');
    getline(pedmapfile,string_channel,' ');
    getline(pedmapfile,string_pedMean,' ');
    getline(pedmapfile,string_pedRMS);
    
    input_tile = atoi(string_tile.c_str());
    input_chip = atoi(string_chip.c_str());
    input_channel = atoi(string_channel.c_str());
    input_pedMean = atof(string_pedMean.c_str());
    input_pedRMS = atof(string_pedRMS.c_str());

    pedMean[input_tile][input_chip][input_channel] = input_pedMean;
  }

  ///////////////////////////////////
  // Set Up Pixel Map
  ///////////////////////////////////
  
  float pixelXvals[17][130][64] = {0.0};
  float pixelYvals[17][130][64] = {0.0};
  float pixelZvals[17][130][64] = {0.0};

  ifstream pixelmapfile;
  pixelmapfile.open("channelmap_multitile.dat");

  string string_thetile;
  string string_thechip;
  string string_thechannel;
  string string_X;
  string string_Y;
  string string_Z;

  int input_thetile;
  int input_thechip;
  int input_thechannel;
  float input_X;
  float input_Y;
  float input_Z;

  while(getline(pixelmapfile,string_thetile,' '))
  {
    getline(pixelmapfile,string_thechip,' ');
    getline(pixelmapfile,string_thechannel,' ');
    getline(pixelmapfile,string_X,' ');
    getline(pixelmapfile,string_Y,' ');
    getline(pixelmapfile,string_Z);
    
    input_thetile = atoi(string_thetile.c_str());
    input_thechip = atoi(string_thechip.c_str());
    input_thechannel = atoi(string_thechannel.c_str());
    input_X = atof(string_X.c_str());
    input_Y = atof(string_Y.c_str());
    input_Z = atof(string_Z.c_str());

    pixelXvals[input_thetile][input_thechip][input_thechannel] = input_X;
    pixelYvals[input_thetile][input_thechip][input_thechannel] = input_Y;
    pixelZvals[input_thetile][input_thechip][input_thechannel] = input_Z;
  }

  ///////////////////////////////////
  // Set Up Output File
  ///////////////////////////////////

  TFile *outputFile = new TFile("analysis.root","RECREATE");
  outputFile->cd();

  int eventNum, trackNum, whichAnode;
  int unixT;
  float minT_X, minT_Y, minT_T, maxT_X, maxT_Y, maxT_T;
  float minY_X, minY_Y, minY_T, maxY_X, maxY_Y, maxY_T;
  vector<float> trackHitX, trackHitY, trackHitT, trackHitC;

  TTree *trackTree = new TTree("trackTree","");
  trackTree->Branch("eventNum",&eventNum);
  trackTree->Branch("trackNum",&trackNum);
  trackTree->Branch("whichAnode",&whichAnode);
  trackTree->Branch("unixT",&unixT);
  trackTree->Branch("minT_X",&minT_X);
  trackTree->Branch("minT_Y",&minT_Y);
  trackTree->Branch("minT_T",&minT_T);
  trackTree->Branch("maxT_X",&maxT_X);
  trackTree->Branch("maxT_Y",&maxT_Y);
  trackTree->Branch("maxT_T",&maxT_T);
  trackTree->Branch("minY_X",&minY_X);
  trackTree->Branch("minY_Y",&minY_Y);
  trackTree->Branch("minY_T",&minY_T);
  trackTree->Branch("maxY_X",&maxY_X);
  trackTree->Branch("maxY_Y",&maxY_Y);
  trackTree->Branch("maxY_T",&maxY_T);
  trackTree->Branch("trackHitX",&trackHitX);
  trackTree->Branch("trackHitY",&trackHitY);
  trackTree->Branch("trackHitT",&trackHitT);
  trackTree->Branch("trackHitC",&trackHitC);

  int trigNum;
  int trigT;
  int recentUnixTime;
  int ioGroup;

  TTree *triggerTree = new TTree("triggerTree","");
  triggerTree->Branch("trigNum",&trigNum);
  triggerTree->Branch("trigT",&trigT);
  triggerTree->Branch("unixT",&recentUnixTime);
  triggerTree->Branch("ioGroup",&ioGroup);
  
  ///////////////////////////////////
  // Load Input Data
  ///////////////////////////////////

  TFile* inputfile = new TFile(inputfilename,"READ");
  
  TTreeReader reader("tree", inputfile);
  TTreeReaderValue<unsigned char> adc_counts(reader, "adc_counts");
  TTreeReaderValue<unsigned char> tile_id(reader, "tile_id");
  TTreeReaderValue<unsigned char> chip_id(reader, "chip_id");
  TTreeReaderValue<unsigned char> channel_id(reader, "channel_id");
  TTreeReaderValue<unsigned char> io_group(reader, "io_group");
  TTreeReaderValue<unsigned char> packet_type(reader, "packet_type");
  TTreeReaderValue<unsigned long long> timestamp(reader, "timestamp");
  
  ///////////////////////////////////
  // Loop Over Data
  ///////////////////////////////////

  eventNum = 0;
  trackNum = 0;
  trigNum = 0;
  recentUnixTime = 0;
  int entryNum = 0;
  
  vector<vector<float>> hits;
  vector<vector<double>> hits_scaled;

  while(reader.Next())
  {
    if(eventNum >= maxEvents)
    {
      break;
    }
    else if((entryNum % dataChunkSize == 0) && (entryNum != 0))
    {
      auto dbscan = DBSCAN<std::vector<double>, double>();
      
      dbscan.Run(&hits_scaled,3,scanRadius,4);
      auto noise = dbscan.Noise;
      auto clusters = dbscan.Clusters;
      
      for(int i = 0; i < (int) clusters.size(); i++)
      {
        const int clusterSize = clusters.at(i).size();
        if(clusterSize > minClusterSize)
        {
          minT_T = 9999999999;
          maxT_T = -9999999999;
          minY_Y = 9999999999;
          maxY_Y = -9999999999;
      
          trackHitX.clear();
          trackHitY.clear();
          trackHitT.clear();
          trackHitC.clear();
      
          for(int j = 0; j < clusterSize; j++)
          {
            const int hitIndex = clusters.at(i).at(j);
      
            trackHitX.push_back(hits[hitIndex][0]);
            trackHitY.push_back(hits[hitIndex][1]);
            trackHitT.push_back(hits[hitIndex][2]);
            trackHitC.push_back(hits[hitIndex][3]);
      
            if(hits[hitIndex][2] < minT_T)
      	    {
              minT_X = hits[hitIndex][0];
              minT_Y = hits[hitIndex][1];
              minT_T = hits[hitIndex][2];
      	    }
            if(hits[hitIndex][2] > maxT_T)
      	    {
              maxT_X = hits[hitIndex][0];
              maxT_Y = hits[hitIndex][1];
              maxT_T = hits[hitIndex][2];
      	    }
            if(hits[hitIndex][1] < minY_Y)
      	    {
              minY_X = hits[hitIndex][0];
              minY_Y = hits[hitIndex][1];
              minY_T = hits[hitIndex][2];
      	    }
            if(hits[hitIndex][1] > maxY_Y)
      	    {
              maxY_X = hits[hitIndex][0];
              maxY_Y = hits[hitIndex][1];
              maxY_T = hits[hitIndex][2];
      	    }
          }

	  if(hits_scaled[clusters.at(i).at(0)][2] > 0.0)
	  {
            whichAnode = 1;
	  }
	  else
	  {
            whichAnode = -1;
	  }

	  unixT = (int) hits_scaled[clusters.at(i).at(0)][3];

	  if(unixT > 0)
	  {
            trackTree->Fill();    
            trackNum++;
	  }
        }
      }

      hits.clear();
      hits_scaled.clear();
      entryNum = 0;
      eventNum++;
    }
    else
    {
      if((*packet_type == 0) && ((*timestamp < 1000000) || (*timestamp > 11000000)))
      {
        continue;
      }
      else if((*tile_id < 1) || (*tile_id > 16) || (*chip_id < 11) || (*chip_id > 110) || (*channel_id < 0) || (*channel_id > 63) || (pedMean[*tile_id][*chip_id][*channel_id] == 0.0) || (*packet_type > 0))
      {
        if(*packet_type == 4)
	{
          recentUnixTime = (int) *timestamp;
	}
	else if((*packet_type == 7) && (recentUnixTime > 0))
	{
          trigT = (int) *timestamp;
	  ioGroup = (int) *io_group;

	  triggerTree->Fill();
	  trigNum++;
	}
	
	entryNum++;
        continue;
      }
      
      float Xval, Yval, Zval, Tval, Cval;
      Xval = pixelXvals[*tile_id][*chip_id][*channel_id];
      Yval = pixelYvals[*tile_id][*chip_id][*channel_id];
      Zval = pixelZvals[*tile_id][*chip_id][*channel_id];
      Tval = (float) *timestamp;
      Cval = ((float) *adc_counts) - pedMean[*tile_id][*chip_id][*channel_id];

      if((Xval == 0.0) && (Yval == 0.0) && (Zval == 0.0))
      {
	entryNum++;
        continue;
      }

      vector<float> hit;
      hit.push_back(Xval);
      hit.push_back(Yval);
      hit.push_back(Tval);
      hit.push_back(Cval);

      vector<double> hit_scaled;
      hit_scaled.push_back(Xval);
      hit_scaled.push_back(Yval);
      if(Zval > 0.0)
      {
        hit_scaled.push_back(Tval/timeSF);
      }
      else
      {
        hit_scaled.push_back(-1.0*Tval/timeSF);
      }
      hit_scaled.push_back(recentUnixTime);

      hits.push_back(hit);
      hits_scaled.push_back(hit_scaled);
      
      entryNum++;
    }
  }

  ///////////////////////////////////
  // Write Output File
  ///////////////////////////////////

  outputFile->cd();
  trackTree->Write();
  triggerTree->Write();
  outputFile->Close();

  timer.Stop();
  cout << "Track Clustering Time:  " << timer.CpuTime() << " sec." << endl;

  return 0;
}
