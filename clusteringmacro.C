// root -l -b -q clusteringmacro.C

void clusteringmacro(){

  gROOT->ProcessLine(".L clustering.cpp");
  gROOT->ProcessLine("readBinMatrix(\"input/cs137.tpx\", 1200);");

  // gROOT->ProcessLine(".L clustering_cal.cpp");
  // gROOT->ProcessLine("readBinMatrix(\"input/cs137.tpx\", 1200, true, 0.11, 2.47);");

  // gROOT->ProcessLine(".L clustering_write_txt.cpp");
  // gROOT->ProcessLine("readBinMatrix(\"input/cs137.tpx\", 1200);");

  // gROOT->ProcessLine(".L clustering_multiframes.cpp");
  // gROOT->ProcessLine("readBinMatrix(\"input/cs137.tpx\", 1200, 2);");
  
}