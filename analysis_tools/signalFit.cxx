void signalFit(string file, double lo = 0, double hi = 1E8)
{
  using namespace RooFit;

  TFile* dataFile = new TFile(file.c_str(), "READ");
  TTree* dataTree = (TTree*)dataFile->Get("save");

  RooWorkspace w("w", kTRUE);
  //w.factory("CBShape::sig(x[1.6,6],mean[3.097,3.0,3.2],sigma[0.285,0.2,0.35],alpha[1.151,-10,10],n[5,0,150])");
  w.factory("Gaussian::sig(mass[1.6,6],mean[3.134,3.,3.2],sigma[0.255,0.2,0.35])");
  //w.factory("CBShape::sig(x[1.6,6],mean[3.097,3.0,3.2],sigma[0.285,0.25,0.35],alpha[2.43],n[13.869])");
  w.factory("CEXPR::bkg('exp(a*mass)/(exp(b*mass+c)+d)',mass[1.6,6],a[-0.9,-5,0],b[-10,-50,0],c[1.,0,100],d[1.,0,100])");
  w.factory("SUM::model(nsig[0,25000]*sig,nbkg[0,120000]*bkg)");

  RooDataSet data("data", "data", dataTree, w::mass, "mass>1.6 && mass<6.");
  
  RooFitResult* res = w::model.fitTo(data, Save());
  res->Print();

  RooPlot* frame = w::mass.frame();
  data.plotOn(frame, Binning(44));
  w::model.plotOn(frame);
  w::model.plotOn(frame, Components(w::sig), LineColor(kRed));
  w::model.plotOn(frame, Components(w::bkg), LineColor(kBlue), LineStyle(kDashed));
  //w::model.paramOn(frame);

  frame->Draw();
}
