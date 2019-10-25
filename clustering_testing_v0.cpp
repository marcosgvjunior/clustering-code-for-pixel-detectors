//Institute of Physics - Federal University of Rio de Janeiro
//Interdisciplinary Academic Master's Degree in Applied Physics
//Student: Marcos Vieira
//October, 2019
//Working on ROOT 5.34/36, Windows 10 x64

#include <stdio.h>
#include <map>
#include <string>
#include <iterator>

// ROOT
#include "TH2D.h"
#include "TH1F.h"

#define MAX_NAME 100

void readBinMatrix( char* inputFile, int totalFrameNumber );

void readBinMatrix( char* inputFile, int totalFrameNumber )
{
  ofstream myfile;

  char *outfilename = new char[ MAX_NAME ];
  sprintf( outfilename,   "clustering_%s.dat", inputFile );
  myfile.open( outfilename, ios::out );

  myfile << Form( "frame\ti\tj\tlabel\tsize\tiTOT\tTOT\n" );

  // Open file
  FILE *file;
  file = fopen( inputFile, "rb" );

  // output file
  char *filename = new char[ MAX_NAME ];
  sprintf( filename,   "clustering_%s.root", inputFile );
  TFile *clusteringHistos = new TFile( filename, "RECREATE" );

  // histograms
  TH1F *multipixelhistogram   = new TH1F( "multipixelhistogram",  "multipixel Histogram",  10000, 0, 10000 );
  TH1F *singlepixelhistogram  = new TH1F( "singlepixelhistogram", "singlepixel Histogram", 10000, 0, 10000 );
  TH1F *allnpixelshistogram   = new TH1F( "allnpixelshistogram",  "All pixel Histogram",   10000, 0, 10000 );
  TH1F *clustersizehistogram  = new TH1F( "clustersizehistogram", "All pixel Histogram",     100, 0,   100 );
  TH1F *clusterperframehisto  = new TH1F( "clusterperframehisto", "All pixel Histogram",     400, 0,   400 );
  TH2D *pixelchargehist       = new TH2D( "pixelchargehist",      "nPixels x De",          10000, 0, 10000, 100, 0, 100 );
  TH2D *AllClustFrame         = new TH2D(  "AllClustFrame",       "Pixel Matrix",            256, 0,   255, 256, 0, 255 );

  // frame buffer
  unsigned short *frame;

  // allocate memory for one frame
  frame = ( unsigned short * )malloc( sizeof( unsigned short )*256*256 );

  // reading in sequency
  int col = 0, row = 0;

  // last pixel TOT/energy
  int lastpixel = 0, actualpixel = 0, clustersizekeep = 0, clusterchargekeep = 0;

  // map for pixel and cluster's label
  std::map<int, int> pixelLabel;

  // map for cluster and TOT
  std::map<int, int> clusterTOT;

  // map for index
  std::map<int, int> clusterIndexI;
  std::map<int, int> clusterIndexJ;

  // cluster size
  std::map<int, int> clusterSize;

  // clusters merged
  std::map<int, int> clustersMerged;
  std::map<int, int> pixelMerged;
  std::map<int, int> ipixelTOT;

  // total number of pixels
  int npixels    = 256*256;

  int wronglabel = 0, correctlabel = 0, wrongindex = 0;
  int uplast = 0, up = 0, upafter = 0, last = 0;

  for( int frameCounter = 0; frameCounter < totalFrameNumber; frameCounter++ )
  {
    //reading one frame into buffer
    fread( frame, sizeof( unsigned short )*256*256, 1, file );

    for( int n = 0; n < npixels; n++ )
    {
      col         = ( n%256 );
      row         = floor( n/256 );
      actualpixel = frame[n];

      if( actualpixel != 0 )
      {

        up      = pixelLabel[ col + 256 * ( row - 1 ) ];      uplast = pixelLabel[ col - 1 + 256 * ( row - 1 ) ];
        upafter = pixelLabel[ col + 1 + 256 * ( row - 1 ) ];  last   = pixelLabel[n-1];

        if( frameCounter == 0 ){ AllClustFrame -> Fill( row, col ); }

        if( lastpixel != 0 )
        {
          pixelLabel[n]               = last;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];

          frame[n] = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( uplast > 0 )
        {
          pixelLabel[n]               = uplast;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];

          frame[n] = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( up > 0  )
        {
          pixelLabel[n]               = up;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];

          frame[n] = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( upafter > 0 )
        {
          pixelLabel[n]               = upafter;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];

          frame[n] = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else
        {
          pixelLabel[n]               = ( col + 256 * row );
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];

          frame[n] = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        // to treat the case when more than one condition satisfies
        if( ( ( upafter > 0 ) && ( last > 0 ) )
         && ( upafter != last ) )
        {
          if( upafter > last ) {
            correctlabel = last; wronglabel = upafter; wrongindex = col + 1 + 256 * ( row - 1 );
          } else{ wronglabel = last; correctlabel = upafter; wrongindex = n-1; }

          clusterTOT[correctlabel]   += clusterTOT[wronglabel];
          clusterTOT.erase( wronglabel );

          clusterSize[correctlabel]  += clusterSize[wronglabel];
          clusterSize.erase( wronglabel );

          clustersMerged[correctlabel] = wronglabel;

          pixelMerged[wronglabel] = correctlabel;
          pixelLabel[wrongindex]  = correctlabel;
          pixelLabel[n] = correctlabel;
        }

        else if( ( ( upafter > 0 ) && ( uplast > 0 ) )
             &&  ( upafter != uplast ) )
        {
          if( upafter > uplast ) {
            correctlabel = uplast; wronglabel = upafter; wrongindex = col + 1 + 256 * ( row - 1 );
          } else{ wronglabel = uplast; correctlabel = upafter; wrongindex = col - 1 + 256 * ( row - 1 ); }

          clusterTOT[correctlabel]   += clusterTOT[wronglabel];
          clusterTOT.erase( wronglabel );

          clusterSize[correctlabel]  += clusterSize[wronglabel];
          clusterSize.erase( wronglabel );

          clustersMerged[correctlabel] = wronglabel;

          pixelMerged[wronglabel] = correctlabel;
          pixelLabel[wrongindex]  = correctlabel;
          pixelLabel[n] = correctlabel;
        }

        else if( ( ( upafter > 0 ) && ( up > 0 ) )
             &&  ( upafter != up ) )
        {
          if( upafter > up ) {
            correctlabel = up; wronglabel = upafter; wrongindex = col + 1 + 256 * ( row - 1 );
          } else{ wronglabel = up; correctlabel = upafter; wrongindex = col + 256 * ( row - 1 ); }

          clusterTOT[correctlabel]   += clusterTOT[wronglabel];
          clusterTOT.erase( wronglabel );

          clusterSize[correctlabel]  += clusterSize[wronglabel];
          clusterSize.erase( wronglabel );

          clustersMerged[correctlabel] = wronglabel;

          pixelMerged[wronglabel] = correctlabel;
          pixelLabel[wrongindex]  = correctlabel;
          pixelLabel[n] = correctlabel;
        }

        else if( ( ( upafter > 0 ) && ( up > 0 ) )
             &&  ( upafter == up ) && ( pixelMerged.count( up ) > 0 ))
        {
          pixelLabel[col + 256 * ( row - 1 ) ]    = pixelMerged.find( up ) -> second ;
          pixelLabel[col + 1 + 256 * ( row - 1 )] = pixelMerged.find( upafter ) -> second ;
        }
      }
      lastpixel  = actualpixel;
    }

    std::map<int,int>::iterator it3=clustersMerged.begin();
    while ( it3 != clustersMerged.end() )
    {
      if( clusterTOT.count( it3->second ) != 0 )
      {
        clusterTOT[it3->first]  += clusterTOT[it3->second];
        clusterTOT.erase( it3->second );

        clusterSize[it3->first] += clusterSize[it3->second];
        clusterSize.erase( it3->second );
      }
      ++it3;
    }

    if( frameCounter%1 == 0 ){ cout << "Frame number: " << frameCounter << "\tNumber of clusters: " << clusterTOT.size() << '\n'; }

    std::map<int,int>::iterator it1=clusterTOT.begin();
    std::map<int,int>::iterator it2=clusterSize.begin();
    while ( it1 != clusterTOT.end() )
    {
      clustersizekeep   =  it2->second;
      clusterchargekeep =  it1->second;

      if( clustersizekeep == 1 ){ singlepixelhistogram -> Fill( clusterchargekeep );
      } else{ multipixelhistogram -> Fill( clusterchargekeep ); }

      allnpixelshistogram   -> Fill( clusterchargekeep );
      clustersizehistogram  -> Fill( clustersizekeep );
      pixelchargehist       -> Fill( clusterchargekeep, clustersizekeep );

      ++it1;  ++it2;
    }
    clusterperframehisto  -> Fill( clusterTOT.size() );

    std::map<int,int>::iterator it6=pixelLabel.begin();
    while ( it6 != pixelLabel.end() )
    {
      if( pixelMerged.count( it6->second ) != 0 )
      {
        pixelLabel[it6->first] = pixelMerged.find( it6->second )->second;
      }
      ++it6;
    }

    std::map<int,int>::iterator it7=ipixelTOT.begin();
    while ( it7 != ipixelTOT.end() )
    {
      //myfile << Form( " i= %d\tj = %d\tlabel = %d\tsize = %d\tiTOT = %d\tTOT = %d \n", clusterIndexI[it7->first], clusterIndexJ[it7->first], pixelLabel[it7->first], clusterSize[pixelLabel[it7->first]], it7->second, clusterTOT[pixelLabel[it7->first]] );
      myfile << Form( "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", frameCounter, clusterIndexI[it7->first], clusterIndexJ[it7->first], pixelLabel[it7->first], clusterSize[pixelLabel[it7->first]], it7->second, clusterTOT[pixelLabel[it7->first]] );
      ++it7;
    }

    pixelLabel.clear();
    clusterTOT.clear();
    clusterSize.clear();
    clustersMerged.clear();
    pixelMerged.clear();
    clusterIndexI.clear();
    clusterIndexJ.clear();
    ipixelTOT.clear();
  }

  TProfile* pixelchargehistprofile = pixelchargehist -> ProfileX();

  clusteringHistos -> Write();
  //clusteringHistos -> Close();

  fclose( file );
  free( frame );
}
