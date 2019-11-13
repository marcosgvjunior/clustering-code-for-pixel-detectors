//Institute of Physics - Federal University of Rio de Janeiro
//Interdisciplinary Academic Master's Degree in Applied Physics
//Student: Marcos Vieira
//October, 2019
//Working on ROOT 5.34/36, Windows 10 x64

#include <stdio.h>
#include <map>
#include <string>
#include <iterator>
#include <time.h>

// ROOT
#include "TH2D.h"
#include "TH1F.h"

#define MAX_NAME 100

void readBinMatrix( char* inputFile, int totalFrameNumber );

void readBinMatrix( char* inputFile, int totalFrameNumber, int frameSetSize )
{
  time_t start = time(NULL);

  std::string str = inputFile;
  str = str.substr(0,str.length()-4);
  str = str.erase(0,3);

  // Open file
  FILE *file;
  file = fopen( inputFile, "rb" );

  // output file
  char *filename = new char[ MAX_NAME ];
  sprintf( filename, "clustering_%s.root", str );
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
  frame = ( unsigned short * )malloc( sizeof( unsigned short )*frameSetSize*256*256 );

  // reading in sequency
  int col = 0, row = 0;

  // last pixel TOT/energy
  int lastpixel = 0, dePixel = 0, clustersizekeep = 0, clusterchargekeep = 0;

  // map for pixel and cluster's label
  std::map<int, int> pixelLabel;

  // map for cluster and TOT
  std::map<int, int> clusterTOT;

  // cluster size
  std::map<int, int> clusterSize;

  // clusters merged
  std::map<int, int> clustersMerged;
  std::map<int, int> pixelMerged;

  // total number of pixels
  int npixels    = frameSetSize*256*256;

  int wronglabel = 0, correctlabel = 0, wrongindex = 0;
  int uplast     = 0, up           = 0, upafter    = 0, last = 0;

  for( int frameCounter = 0; frameCounter < totalFrameNumber; frameCounter++ )
  {
    //reading one frame into buffer
    fread( frame, sizeof( unsigned short )*frameSetSize*256*256, 1, file );

    for( int n = 0; n < npixels; n++ )
    {
      dePixel   = frame[n];

      if( dePixel != 0 )
      {
        col     = ( n%256 );
        row     = floor( n/256 );

        up      = pixelLabel[ col + 256 * ( row - 1 ) ];      uplast = pixelLabel[ col - 1 + 256 * ( row - 1 ) ];
        upafter = pixelLabel[ col + 1 + 256 * ( row - 1 ) ];  last   = pixelLabel[n-1];

        if( frameCounter == 0 ){ AllClustFrame -> Fill( row, col ); }

        if( (lastpixel != 0) && (n%256 != 0) )
        {
          pixelLabel[n]        = last;
          clusterTOT[last]    += dePixel;
          clusterSize[last]   += 1;
          frame[n]             = 0;
        }

        else if( (uplast > 0) && (n%256 != 0) && (row%256 != 0) )
        {
          pixelLabel[n]        = uplast;
          clusterTOT[uplast]  += dePixel;
          clusterSize[uplast] += 1;
          frame[n]             = 0;
        }

        else if( ( up > 0 ) && (row%256 != 0) )
        {
          pixelLabel[n]         = up;
          clusterTOT[up]       += dePixel;
          clusterSize[up]      += 1;
          frame[n]              = 0;
        }

        else if( (upafter > 0) && (n%256 != 255) && (row%256 != 0) )
        {
          pixelLabel[n]         = upafter;
          clusterTOT[upafter]  += dePixel;
          clusterSize[upafter] += 1;
          frame[n]              = 0;
        }

        else
        {
          pixelLabel[n]                     = ( col + 256 * row );
          clusterTOT[( col + 256 * row )]  += dePixel;
          clusterSize[( col + 256 * row )] += 1;
          frame[n]                          = 0;
        }

        // to treat the case when more than one condition satisfies
        if( ( ( upafter > 0 ) && ( last > 0 ) ) && ( upafter != last ) && (n%256 != 0) && (n%256 != 255) && (row%256 != 0) )
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

        else if( ( ( upafter > 0 ) && ( uplast > 0 ) ) && ( upafter != uplast ) && (n%256 != 0) && (n%256 != 255) && (row%256 != 0) )
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

        else if( ( ( upafter > 0 ) && ( up > 0 ) ) && ( upafter != up ) && (n%256 != 255) && (row%256 != 0) )
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

        else if( ( ( upafter > 0 ) && ( up > 0 ) ) && ( upafter == up ) && ( pixelMerged.count( up ) > 0 ) && (row%256 != 0))
        {
          pixelLabel[col + 256 * ( row - 1 ) ]    = pixelMerged.find( up )      -> second ;
          pixelLabel[col + 1 + 256 * ( row - 1 )] = pixelMerged.find( upafter ) -> second ;
        }
      }
      lastpixel  = dePixel;
    }

    std::map<int,int>::iterator it1 = clustersMerged.begin();
    while ( it1 != clustersMerged.end() )
    {
      if( clusterTOT.count( it1->second ) != 0 )
      {
        clusterTOT[it1->first]  += clusterTOT[it1->second];
        clusterTOT.erase( it1->second );

        clusterSize[it1->first] += clusterSize[it1->second];
        clusterSize.erase( it1->second );
      }
      ++it1;
    }

    if( frameCounter%1 == 0 ){ cout << "Frame number: " << frameCounter << "\tNumber of clusters: " << clusterTOT.size() << '\n'; }

    std::map<int,int>::iterator it2 = clusterTOT.begin();
    std::map<int,int>::iterator it3 = clusterSize.begin();
    while ( it2 != clusterTOT.end() )
    {
      clusterchargekeep =  it2->second;
      clustersizekeep   =  it3->second;

      if( clustersizekeep == 1 ){ singlepixelhistogram -> Fill( clusterchargekeep );
      } else{ multipixelhistogram -> Fill( clusterchargekeep ); }

      allnpixelshistogram   -> Fill( clusterchargekeep );
      clustersizehistogram  -> Fill( clustersizekeep );
      pixelchargehist       -> Fill( clusterchargekeep, clustersizekeep );

      ++it2;  ++it3;
    }
    clusterperframehisto  -> Fill( clusterTOT.size() );

    pixelLabel.clear();
    clusterTOT.clear();
    clusterSize.clear();
    clustersMerged.clear();
    pixelMerged.clear();
  }
  TProfile* pixelchargehistprofile = pixelchargehist -> ProfileX();

  clusteringHistos -> Write();
  clusteringHistos -> Close();

  fclose( file );
  free( frame );

  time_t end = time(NULL);
  cout << "Elapsed time: " << (double)(end-start) << " seconds\n";
}
