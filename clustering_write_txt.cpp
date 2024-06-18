#include <stdio.h>
#include <map>
#include <string>
#include <iterator>
#include <time.h>

// ROOT
// #include "TH2D.h"
// #include "TH1F.h"

#define MAX_NAME 100

int mapiterator( std::map<int, int> &map1, std::map<int, int> &map2 );
void readBinMatrix( char* inputFile, int totalFrameNumber );

void readBinMatrix( char* inputFile, int totalFrameNumber )
{
  time_t start = time(NULL);

  ofstream myfile;

  std::string str = inputFile;
  str = str.substr(0,str.length()-4);
  str = str.erase(0,6);

  char *outfilename = new char[ MAX_NAME ];
  sprintf( outfilename,   "output/clustering_%s_write.dat", (str).c_str() );
  myfile.open( outfilename, ios::out );

  myfile << Form( "frame\ti\tj\tlabel\tsize\tiTOT\tTOT\n" );

  // Open file
  FILE *file;
  file = fopen( inputFile, "rb" );

  // output file
  char *filename = new char[ MAX_NAME ];
  sprintf( filename, "output/clustering_%s_write.root", (str).c_str() );
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
  int uplast     = 0, up           = 0, upafter    = 0, last        = 0;
  int correction = 1, outfirst     = 0, outsecond  = 0, outpixlabel = 0;

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

        if( (lastpixel != 0) && (n%256 != 0) )
        {
          pixelLabel[n]               = last;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];
          frame[n]                    = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( (uplast > 0) && (n%256 != 0) )
        {
          pixelLabel[n]               = uplast;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];
          frame[n]                    = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( up > 0  )
        {
          pixelLabel[n]               = up;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];
          frame[n]                    = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else if( (upafter > 0) && (n%256 != 255) )
        {
          pixelLabel[n]               = upafter;
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];
          frame[n]                    = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        else
        {
          pixelLabel[n]               = ( col + 256 * row );
          clusterTOT[pixelLabel[n]]  += frame[n];
          clusterSize[pixelLabel[n]] += 1;
          ipixelTOT[n]                = frame[n];
          frame[n]                    = 0;

          clusterIndexI[n] = row;
          clusterIndexJ[n] = col;
        }

        // to treat the case when more than one condition satisfies
        if( ( ( upafter > 0 ) && ( last > 0 ) ) && ( upafter != last ) && (n%256 != 0) && (n%256 != 255) )
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

        else if( ( ( upafter > 0 ) && ( uplast > 0 ) ) && ( upafter != uplast ) && (n%256 != 0) && (n%256 != 255) )
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

        else if( ( ( upafter > 0 ) && ( up > 0 ) ) && ( upafter != up ) && (n%256 != 255) )
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

        else if( ( ( upafter > 0 ) && ( up > 0 ) ) && ( upafter == up ) && ( pixelMerged.count( up ) > 0 ))
        {
          pixelLabel[col + 256 * ( row - 1 ) ]    = pixelMerged.find( up ) -> second ;
          pixelLabel[col + 1 + 256 * ( row - 1 )] = pixelMerged.find( upafter ) -> second ;
        }
      }
      lastpixel  = actualpixel;
    }

    std::map<int,int>::iterator it1=clustersMerged.begin();
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

    if( frameCounter%100 == 0 ){ cout << "Frame number: " << frameCounter << "\tNumber of clusters: " << clusterTOT.size() << '\n'; }

    std::map<int,int>::iterator it2=clusterTOT.begin();
    std::map<int,int>::iterator it3=clusterSize.begin();
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

    while( correction != 0 )
    {
      correction = mapiterator( pixelLabel, pixelMerged );
    }
    correction = 1;

    std::map<int,int>::iterator it4=ipixelTOT.begin();
    while ( it4 != ipixelTOT.end() )
    {
      outfirst = it4->first; outsecond = it4->second; outpixlabel = pixelLabel[outfirst];
      myfile << Form( "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", frameCounter, clusterIndexI[outfirst], clusterIndexJ[outfirst], pixelLabel[outfirst], clusterSize[outpixlabel], outsecond, clusterTOT[outpixlabel] );
      ++it4;
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

  time_t end = time(NULL);
  cout << "Elapsed time: " << (double)(end-start) << " seconds\n";
}

int mapiterator( std::map<int, int> &map1, std::map<int, int> &map2 )
{
  int endcondition = 0;

  std::map<int,int>::iterator it5=map1.begin();
  while ( it5 != map1.end() )
  {
    if( map2.count( it5->second ) != 0 )
    {
      map1[it5->first] = map2.find( it5->second )->second;
    }
    if( map2.count( map1[it5->first] ) != 0 ){ endcondition +=1; }
    ++it5;
  }

  return endcondition;
}
