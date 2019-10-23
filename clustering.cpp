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

#define MAX_NAME 100

void readBinMatrix( char* inputFile, int totalFrameNumber );

void readBinMatrix( char* inputFile, int totalFrameNumber )
{
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
  //TH2D *AllClustFrame         = new TH2D(  "AllClustFrame",       "Pixel Matrix",            256, 0,   255, 256, 0, 255 );

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

  // cluster size
  std::map<int, int> clusterSize;

  // clusters merged
  std::map<int, int> clustersMerged;

  // total number of pixels
  int npixels = 256*256;

  for( int frameCounter = 0; frameCounter < totalFrameNumber; frameCounter++ )
  {
    //reading one frame into buffer
    fread( frame, sizeof( unsigned short )*256*256, 1, file );

    for( int n = 0; n < npixels; n++ )
    {
      col         = ( n%256 );
      row         = floor(n/256);
      actualpixel = frame[ n ];

      if( actualpixel != 0 )
      {
        //if( frameCounter == 2 ){ AllClustFrame -> Fill( row, col); }

        if( lastpixel != 0 )
        {
          pixelLabel[n]                        = pixelLabel[n-1];
          clusterTOT[pixelLabel[n]]           += frame[ n ];
          clusterSize[pixelLabel[n]]          += 1;
          frame[ n ] = 0;
        }

        else if( pixelLabel[ col + 256 * ( row - 1 ) ] > 0  )
        {
          pixelLabel[n]               = pixelLabel[col + 256 * ( row - 1 )];
          clusterTOT[pixelLabel[n]]  += frame[ n ];
          clusterSize[pixelLabel[n]] += 1;
          frame[ n ] = 0;
        }

        else if( pixelLabel[ col - 1 + 256 * ( row - 1 ) ] > 0 )
        {
          pixelLabel[n]               = pixelLabel[col - 1 + 256 * ( row - 1 )];
          clusterTOT[pixelLabel[n]]  += frame[ n ];
          clusterSize[pixelLabel[n]] += 1;
          frame[ n ] = 0;
        }

        else if( pixelLabel[ col + 1 + 256 * ( row - 1 ) ] > 0 )
        {
          pixelLabel[n]               = pixelLabel[col + 1 + 256 * ( row - 1 )];
          clusterTOT[pixelLabel[n]]  += frame[ n ];
          clusterSize[pixelLabel[n]] += 1;
          frame[ n ] = 0;
        }

        else
        {
          pixelLabel[n]               = ( col + 256 * row );
          clusterTOT[pixelLabel[n]]  += frame[ n ];
          clusterSize[pixelLabel[n]] += 1;
          frame[ n ] = 0;
        }

        // to treat the case when more than one condition satisfies
        if( ( ( pixelLabel[ col + 1 + 256 * ( row - 1 ) ] > 0 ) && ( lastpixel > 0 ) )
         && ( pixelLabel[col + 1 + 256 * ( row - 1 )] != pixelLabel[n-1] ) )
        {
          clusterTOT[pixelLabel[n-1]]   += clusterTOT[pixelLabel[col + 1 + 256 * ( row - 1 )]];
          clusterTOT.erase(pixelLabel[col + 1 + 256 * ( row - 1 )]);

          clusterSize[pixelLabel[n-1]]  += clusterSize[pixelLabel[col + 1 + 256 * ( row - 1 )]];
          clusterSize.erase(pixelLabel[col + 1 + 256 * ( row - 1 )]);

          clustersMerged[pixelLabel[n-1]] = pixelLabel[col + 1 + 256 * ( row - 1 )];
        }

        else if( ( ( pixelLabel[ col + 1 + 256 * ( row - 1 ) ] > 0 ) && ( pixelLabel[ col - 1 + 256 * ( row - 1 ) ] > 0 ) )
             &&  ( pixelLabel[col + 1 + 256 * ( row - 1 )] != pixelLabel[col - 1 + 256 * ( row - 1 )] ) )
        {
          clusterTOT[pixelLabel[col + 1 + 256 * ( row - 1 )]]  += clusterTOT[pixelLabel[col - 1 + 256 * ( row - 1 )]];
          clusterTOT.erase(pixelLabel[col - 1 + 256 * ( row - 1 )]);

          clusterSize[pixelLabel[col + 1 + 256 * ( row - 1 )]]  += clusterSize[pixelLabel[col - 1 + 256 * ( row - 1 )]];
          clusterSize.erase(pixelLabel[col - 1 + 256 * ( row - 1 )]);

          clustersMerged[pixelLabel[col + 1 + 256 * ( row - 1 )]] = pixelLabel[col - 1 + 256 * ( row - 1 )];
        }
      }
      lastpixel  = actualpixel;
    }

    std::map<int,int>::iterator it3=clustersMerged.begin();

    while (it3 != clustersMerged.end())
    {
      if(clusterTOT.count( it3->second ) != 0)
      {
        clusterTOT[it3->first]  += clusterTOT[it3->second];
        clusterTOT.erase(it3->second);

        clusterSize[it3->first] += clusterSize[it3->second];
        clusterSize.erase(it3->second);
      }
      ++it3;
    }

    if( frameCounter%100 == 0 ){ cout << "Frame number: " << frameCounter << "\tNumber of clusters: " << clusterTOT.size() << '\n'; }

    std::map<int,int>::iterator it1=clusterTOT.begin();
    std::map<int,int>::iterator it2=clusterSize.begin();

    while (it1 != clusterTOT.end())
    {
      clustersizekeep   =  it2->second;
      clusterchargekeep =  it1->second;

      if( clustersizekeep == 1 ){ singlepixelhistogram -> Fill( clusterchargekeep);
      } else{ multipixelhistogram -> Fill( clusterchargekeep ); }

      allnpixelshistogram   -> Fill( clusterchargekeep );
      clustersizehistogram  -> Fill( clustersizekeep );
      pixelchargehist       -> Fill( clusterchargekeep, clustersizekeep );

      ++it1;  ++it2;
    }
    clusterperframehisto  -> Fill( clusterTOT.size() );

    pixelLabel.clear();
    clusterTOT.clear();
    clusterSize.clear();
    clustersMerged.clear();
  }

  TProfile* pixelchargehistprofile = pixelchargehist -> ProfileX();

  clusteringHistos -> Write();
  //clusteringHistos -> Close();

  fclose( file );
  free( frame );
}

