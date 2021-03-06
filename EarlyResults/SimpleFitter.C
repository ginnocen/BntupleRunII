#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TLatex.h>

void SimpleFitter()
{
   
   TFile*output=new TFile("histoX.root","read");    
   TH1D*h=(TH1D*)output->Get("h");

   h->Sumw2();
   TF1 *f = new TF1("f","[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) +[3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5])+[6]+[7]*x+[8]*x*x+[9]*x*x*x");
   f->SetLineColor(4);
   f->SetParameters(21,3.68,0.0057,1.34714e+00,3.86,0.006,-2.21334e+04,1.01661e+04,-1.05364e+03,1.95474e+04);
   f->FixParameter(1,3.686);
   f->FixParameter(2,0.00357);
   f->FixParameter(4,3.8725);
   f->FixParameter(5,0.005);

   h->Fit("f","LL");
   h->Fit("f","");
   h->Fit("f","LL");
   
   h->Fit("f","LL","",3.65,3.94);
   h->Fit("f","LL","",3.65,3.94);

   f->ReleaseParameter(1);
   f->ReleaseParameter(2);
   f->ReleaseParameter(4);
   f->ReleaseParameter(5);
   h->Fit("f","LL","",3.65,3.94);
   h->SetXTitle("m(J/#psi#pi^{+}#pi^{-}) [GeV]");
   h->SetYTitle("Entries");
   h->SetStats(0);

   h->SetAxisRange(0,h->GetMaximum()*1.3	,"Y");
   TF1 *f2 = new TF1("f2","[6]+[7]*x+[8]*x*x+[9]*x*x*x");
   f2->SetParameter(6,f->GetParameter(6));
   f2->SetParameter(7,f->GetParameter(7));
   f2->SetParameter(8,f->GetParameter(8));
   f2->SetParameter(9,f->GetParameter(9));
   TF1 *f3 = new TF1("f3","[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])");
   f3->SetParameter(0,f->GetParameter(0));
   f3->SetParameter(1,f->GetParameter(1));
   f3->SetParameter(2,f->GetParameter(2));

   TF1 *f4 = new TF1("f4","[3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5])");
   f4->SetParameter(3,f->GetParameter(3));
   f4->SetParameter(4,f->GetParameter(4));
   f4->SetParameter(5,f->GetParameter(5));
      
   f->SetLineColor(4);
   f2->SetLineColor(4);
   f3->SetRange(3.65,3.94);
   f4->SetRange(3.65,3.94);
   f2->SetRange(3.65,3.94);
   f2->SetLineStyle(2);
   f3->SetLineStyle(2);
   f4->SetLineStyle(2);
   f2->Draw("same");
   f3->SetLineColor(2);
   f3->SetFillStyle(3004);
   f3->SetFillColor(2);
   f3->Draw("same");
   f4->SetLineColor(2);
   f4->SetFillStyle(3004);
   f4->SetFillColor(2);
   f4->Draw("same");

   TLatex *l = new TLatex(3.7,70./80*h->GetMaximum(),"#psi(2S)");
   l->Draw();
   TLatex *l2 = new TLatex(3.875,50./80*h->GetMaximum(),"X(3872)");
   l2->Draw();
   TLatex *l3 = new TLatex(3.812,70./80*h->GetMaximum(),"CMS Preliminary");
   l3->Draw();
   TLatex *l4 = new TLatex(3.78,60./80*h->GetMaximum(),"pp #sqrt{s_{NN}}=5.02 TeV");
   l4->Draw();
   
   cout<<"********************* SIGNIFICANCE *********************"<<endl;
   
   double meanX=f4->GetParameter(4);
   double sigmaX=f4->GetParameter(5);
   double signalX=f4->Integral(meanX-3*sigmaX,meanX+3*sigmaX)/h->GetBinWidth(0);
   double bkgX=f2->Integral(meanX-3*sigmaX,meanX+3*sigmaX)/h->GetBinWidth(0);
   double signalpsi2=f4->Integral(meanX-3*sigmaX,meanX+3*sigmaX)/h->GetBinWidth(0);
   double bkgpsi2=f2->Integral(meanX-3*sigmaX,meanX+3*sigmaX)/h->GetBinWidth(0);
   double significanceX=signalX/TMath::Sqrt(signalX+bkgX);   
   
   double meanpsi2=f3->GetParameter(1);
   double sigmapsi2=f3->GetParameter(2);
   double signalpsi2=f3->Integral(meanpsi2-3*sigmapsi2,meanpsi2+3*sigmapsi2)/h->GetBinWidth(0);
   double bkgpsi2=f2->Integral(meanpsi2-3*sigmapsi2,meanpsi2+3*sigmapsi2)/h->GetBinWidth(0);
   double signalpsi2=f3->Integral(meanpsi2-3*sigmapsi2,meanpsi2+3*sigmapsi2)/h->GetBinWidth(0);
   double bkgpsi2=f2->Integral(meanpsi2-3*sigmapsi2,meanpsi2+3*sigmapsi2)/h->GetBinWidth(0);
   double significancepsi2=signalpsi2/TMath::Sqrt(signalpsi2+bkgpsi2);


   cout<<"significance X ="<<significanceX<<endl;
   cout<<"significance psi(2S)="<<significancepsi2<<endl;
   
   cout<<"yield X="<<signalX<<endl;
   cout<<"yield 2S="<<signalpsi2<<endl;
   cout<<"yield X/yield 2S="<<signalX/signalpsi2<<endl;

      
}
