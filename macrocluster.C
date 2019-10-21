// root -l -b -q macrocluster.C

void macrocluster(){
  gROOT->ProcessLine(".L clustering.cpp");
  readBinMatrix("am241.tpx", 800);
}
