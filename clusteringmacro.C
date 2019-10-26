// root -l -b -q clusteringmacro.C

void clusteringmacro(){
  //gROOT->ProcessLine(".L clustering.cpp");
  gROOT->ProcessLine(".L clustering_write_txt.cpp");
  //readBinMatrix("am241.tpx", 900);
  readBinMatrix("../cs137.tpx", 10);
}
