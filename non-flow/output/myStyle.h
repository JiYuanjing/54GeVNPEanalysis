#ifndef __myStyle_C__
#define __myStyle_C__

#include <string>
#include <iostream>
#include <sstream>

#include "TPaveText.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TText.h"

using namespace std;

void myStyle() {

  // define names for colours
  //  Int_t black  = 1;
  //Int_t red    = 2;
  //Int_t green  = 3;
  //Int_t blue   = 4;
  //Int_t yellow = 5; 
  //Int_t magenta= 6;
  //Int_t cyan   = 7;
  //Int_t purple = 9;
  

////////////////////////////////////////////////////////////////////


  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;
  // Line thickness
  Double_t lhcbWidth    = 2.00; 
  // Text size
  Double_t lhcbTSize    = 0.05; // 0.06  
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","My plots style");
  
  
  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.05); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);
  
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  //lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(0.05,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  //lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(2.,"Y");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  // add LHCb label
  TPaveText *lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);

  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  cout << "-------------------------" << endl;  
  cout << "Set LHCb Style - Feb 2012" << endl;
  cout << "-------------------------" << endl;  
  
}


TPaveText* printMyLHCb(string optLR="L", string optPrelim="Final", string optText="") {
//////////////////////////////////////////////////////////////////////////
// routine to print 'LHCb', 'LHCb Preliminary' on plots 
// options: optLR=L (top left) / R (top right) of plots
//          optPrelim= Final (LHCb), Prelim (LHCb Preliminary), Other
//          optText= text printed if 'Other' specified
////////////////////////////////////////////////////////////////////
  TPaveText* lhcbName = 0;
  const double rightmargin = gStyle->GetPadRightMargin();
  gStyle->SetPadTopMargin(0.05);
  const double topmargin = gStyle->GetPadTopMargin();
  
  if (optLR==string("R") && optPrelim=="Prelim"){    
    lhcbName = new TPaveText(0.70 - rightmargin,
			     0.80 - topmargin,
			     0.95 - rightmargin,
			     0.90 - topmargin,
			     "BRNDC");
  }
  else if (optLR==string("R")){    
    lhcbName = new TPaveText(0.70 - rightmargin,
			     0.85 - topmargin,
			     0.95 - rightmargin,
			     0.95 - topmargin,
			     "BRNDC");
  }
  else if (optLR==string("L")){
    lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			     0.87 - topmargin,
			     gStyle->GetPadLeftMargin() + 0.30,
			     0.95 - gStyle->GetPadTopMargin(),
			     "BRNDC");
  }
  else{
   cout << "printLHCb: option unknown" << optLR << endl;   
  }
  if (optPrelim==string("Final")){
    lhcbName->AddText("#splitline{LHCb}{#scale[1.0]{#sqrt{s} = 7 TeV data}}");
  }
  else if (optPrelim==string("Prelim")){
    lhcbName->AddText("#splitline{LHCb}{#scale[1.0]{Preliminary}}");  
  }
  else if (optPrelim==string("Other")){
    lhcbName->AddText(optText.c_str());
  }
  else{
    cout << "printLHCb: option unknown " << optPrelim << endl;   
  }

  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
  return lhcbName;
}
#endif
