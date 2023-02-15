void invMgg_evMix(){
  TFile *file_in[4];
  

  file_in[0]=new TFile("./photons_day060_01.root");

  TFile *out1=new TFile("pp_gg_evMix.root","recreate");

   
 file_in[0]->cd();
 TH1F *hm_gg= (TH1F*)file_in[0]->Get("hinvM_gg_PT3");
 TH1F *hm_ggMix= (TH1F*)file_in[0]->Get("hinvM_ggMix");

 hm_gg->SetLineWidth(2);
 hm_gg->SetLineColor(1);
 hm_ggMix->SetLineWidth(2);
 hm_ggMix->SetLineColor(2);
 
 //normalization
 hm_ggMix->Scale(hm_gg->Integral(200,230)/hm_ggMix->Integral(200,230));
   
 hm_ggMix->Sumw2(kFALSE);
 hm_gg->Sumw2(kFALSE);

 TH1F *hm_gg_cl=(TH1F*) hm_gg->Clone("hm_gg_cl");
 hm_gg_cl->Add(hm_ggMix,-1);

 hm_gg_cl->SetLineWidth(2);
 hm_gg_cl->SetLineColor(4);
 hm_gg_cl->Sumw2(kFALSE);
 
 //*********************************** 
  TLegend *leg1a;
  leg1a = new TLegend(0.17,0.65,0.48,0.88);
  leg1a->AddEntry(hm_gg,"M_{#gamma#gamma}","l");
  leg1a->AddEntry(hm_ggMix,"M^{MIX}_{#gamma#gamma}","l");
  leg1a->AddEntry(hm_gg_cl,"M_{#gamma#gamma}-M^{MIX}_{#gamma#gamma}","l");

 
  TCanvas *c1=new TCanvas("mass_gg_pi0","mass_gg_pi0");
  c1->Divide(1,1);

  c1->cd(1);	
  hm_gg->Draw();
  hm_gg_cl->Draw("same");
  hm_ggMix->Draw("same");
  leg1a->Draw("same");


 //**************************************
 TH1F *hm_gg_w= (TH1F*)file_in[0]->Get("hMpippimpi0_PT3");
 TH1F *hm_gg_w_Mix= (TH1F*)file_in[0]->Get("hMpippimpi0_Mix");

 hm_gg_w->SetLineWidth(2);
 hm_gg_w->SetLineColor(1);
 hm_gg_w_Mix->SetLineWidth(2);
 hm_gg_w_Mix->SetLineColor(2);


 int a=500;
 int b=530;
 //cout<<hm_gg_w->Integral(a,b)/hm_gg_w_Mix->Integral(a,b)<<endl;
 hm_gg_w_Mix->Scale(hm_gg_w->Integral(a,b)/hm_gg_w_Mix->Integral(a,b));

 hm_gg_w_Mix->Sumw2(kFALSE);
 hm_gg_w->Sumw2(kFALSE);

 TH1F *hm_gg_w_cl=(TH1F*) hm_gg_w->Clone("hm_gg_w_cl");
 hm_gg_w_cl->Add(hm_gg_w_Mix,-1);

 hm_gg_w_cl->SetLineWidth(2);
 hm_gg_w_cl->SetLineColor(4);
 hm_gg_w_cl->Sumw2(kFALSE);


 TLegend *leg2a;
  leg2a = new TLegend(0.17,0.65,0.48,0.88);
  leg2a->AddEntry(hm_gg_w,"M_{#gamma#gamma}","l");
  leg2a->AddEntry(hm_gg_w_Mix,"M^{MIX}_{#gamma#gamma}","l");
  leg2a->AddEntry(hm_gg_w_cl,"M_{#gamma#gamma}-M^{MIX}_{#gamma#gamma}","l");

 
 TCanvas *c2=new TCanvas("mass_gg_w_omega","mass_gg_omega");
 c2->Divide(1,1);

 c2->cd(1);	
 hm_gg_w->Draw("same");
 hm_gg_w_Mix->Draw("same");
 hm_gg_w_cl->Draw("same");
 leg2a->Draw("same");

 //***************************************************	
 out1->cd();
 
 c1->Write();
 c2->Write();

 out1->Close();
	
}

