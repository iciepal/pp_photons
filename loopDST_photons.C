#include "Includes.h" 

using namespace std; 

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;


TLorentzVector TLV_p(float th, float ph, float mom, float mass){

    float px, py, pz, E;
    
    TLorentzVector lv;
    
    px = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Cos(ph*TMath::DegToRad());
    py = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Sin(ph*TMath::DegToRad());
    pz = mom*TMath::Cos(th*TMath::DegToRad());
    E=sqrt(mom*mom + mass*mass);
    lv.SetPxPyPzE(px,py,pz,E);

    return lv;
}



Int_t loopDST_photons() {

    TROOT trDSTAnalysis("trDSTAnalysis", "Simple DST photons analysis Macro");
    new HLoop(kTRUE); 

    TString input = "list_060.list"; 
    TString outFile = "output_photons.root";

    Long64_t nEventsDesired = -1;
    //Long64_t nEventsDesired = 100000;

    gLoop->addFilesList(input);



// read selected categories

    if (!gLoop->setInput("-*,+HStart2Hit,+HParticleCand,+HEmcNeutralCand,+HEmcCluster,+HEventHeader")){ 

        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } 

    //*********************************************  
    // Setting the cache size of the HLoop internal TChain to read data from the DSTs
    // Improves performance of the lustre storage by decreasing load on lustre META servers
    
    gLoop->getChain()->SetCacheSize(8l*1024l*1024l); // 8 MB
    gLoop->getChain()->AddBranchToCache("*", kTRUE);
    gLoop->getChain()->StopCacheLearningPhase();
    gLoop->printCategories(); 

    //************************************************
    // Get category pointers
    //*****************************************
    
  
    HCategory* fParticleCand = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if (!fParticleCand) { cout << "No catParticleCand!" << endl; }
   
    HCategory* fEmcNeutralCand = HCategoryManager::getCategory(catEmcNeutralCand, kTRUE, "catEmcNeutralCand");
    if (!fEmcNeutralCand) { cout << "No catEmcNeutralCand!" << endl; }
 
    HCategory* fEmcCluster = HCategoryManager::getCategory(catEmcCluster, kTRUE, "catEmcCluster");
    if (!fEmcCluster) { cout << "No catEmcCluster!" << endl; }
    
    HCategory * fStart2Hit = HCategoryManager::getCategory(catStart2Hit, kTRUE, "catStart2Hit");
    if (!fStart2Hit) { cout << "No catStart2Hit!" << endl; }

    //*********************************************************
    
    TH1F *htof_PT2=new TH1F("htof_PT2","htof_PT2",2200,-200,2000);
    TH1F *htof_PT3=new TH1F("htof_PT3","htof_PT3",2200,-200,2000);
    TH1F *hcl_size=new TH1F("hcl_size","hcl_size",20,0,20);
    TH1F *hbeta=new TH1F("hbeta","hbeta",100,0,2);
    TH1F *htracklength=new TH1F("htracklength","htracklength",2000,2000,4000);
    TH1F *hg_energy=new TH1F("hg_energy","hg_energy",2000,0,2000);

    TH2F* hmass_mom = new TH2F("hmass_mom","hmass_mom; p*q [MeV/c]; mass",1000,-2000,4000,1000,0,2000);
    TH1F *hmass=new TH1F("hmass","hmass",1000,-2000,2000);

    TH2F *hbeta_mom= new TH2F("hbeta_mom","hbeta_mom",1000,-4000.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_p= new TH2F("hbeta_mom_p","hbeta_mom_p",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pip= new TH2F("hbeta_mom_pip","hbeta_mom_pip",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pim= new TH2F("hbeta_mom_pim","hbeta_mom_pim",1000,-4000.,0.,700,0.,1.4);

    TH1F *hMM_pp=new TH1F("hMM_pp","hMM_pp",1400,0,1400);
    TH1F *hMM2_pp=new TH1F("hMM2_pp","hMM2_pp",900,-1,2);
    TH1F *hinvM_gg_PT3=new TH1F("hinvM_gg_PT3","hinvM_gg_PT3",700,0,700);
    TH1F *hinvM_ggMix=new TH1F("hinvM_ggMix","hinvM_ggMix",700,0,700);

    TH1F* hMpippimpi0_PT3=new TH1F("hMpippimpi0_PT3","hMpippimpi0_PT3; M_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",1000,0,2000);
    TH1F* hMpippimpi0_Mix =new TH1F("hMpippimpi0_Mix","hMpippimpi0_Mix; M_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",1000,0,2000);   

    TH1F *hinvM_gg_pi0_PT3=new TH1F("hinvM_gg_pi0_PT3","hinvM_gg_pi0_PT3",700,0,700);
    TH1F *hinvM_gg_omega_PT3=new TH1F("hinvM_gg_omega_PT3","hinvM_gg_omega_PT3",400,300,1500);
    TH1F *hinvM_gg_omega_pi0_PT3=new TH1F("hinvM_gg_omega_pi0_PT3","hinvM_gg_omega_pi0_PT3",700,0,700);

    //*****************************************************
    
    //PROT PID - mass vs mom

    TCutG *cutP = new TCutG("cutP",11);
    cutP->SetPoint(0,188.0792,1187.5);
    cutP->SetPoint(1,1559.183,1378.472);
    cutP->SetPoint(2,2203.665,1517.361);
    cutP->SetPoint(3,2677.548,1579.861);
    cutP->SetPoint(4,2886.057,1392.361);
    cutP->SetPoint(5,2544.861,1177.083);
    cutP->SetPoint(6,2052.022,888.8889);
    cutP->SetPoint(7,188.0792,593.75);
    cutP->SetPoint(8,86.98397,642.3611);
    cutP->SetPoint(9,200.7161,1190.972);
    cutP->SetPoint(10,188.0792,1187.5);

    //****************************************
    //PI+ PID - mass vs mom
    TCutG *cutPIP = new TCutG("cutPIP",10);
     cutPIP->SetPoint(0,42.75482,218.75);
     cutPIP->SetPoint(1,207.0345,520.8333);
     cutPIP->SetPoint(2,1489.68,791.6667);
     cutPIP->SetPoint(3,1913.016,760.4167);
     cutPIP->SetPoint(4,1495.998,364.5833);
     cutPIP->SetPoint(5,794.6504,27.77775);
     cutPIP->SetPoint(6,598.7784,-20.83336);
     cutPIP->SetPoint(7,30.11792,10.41664);
     cutPIP->SetPoint(8,49.07327,225.6944);
     cutPIP->SetPoint(9,42.75482,218.75);
    //****************************************
    //PI- PID - mass vs mom
     TCutG *cutPIM = new TCutG("cutPIM",11);
     cutPIM->SetPoint(0,-52.02194,208.3333);
     cutPIM->SetPoint(1,-816.5544,684.0278);
     cutPIM->SetPoint(2,-1543.176,836.8056);
     cutPIM->SetPoint(3,-1783.277,701.3889);
     cutPIM->SetPoint(4,-1505.265,395.8333);
     cutPIM->SetPoint(5,-1189.343,170.1389);
     cutPIM->SetPoint(6,-645.9562,3.472193);
     cutPIM->SetPoint(7,-64.65884,-6.944474);
     cutPIM->SetPoint(8,-52.02194,229.1666);
     cutPIM->SetPoint(9,-70.97729,215.2778);
     cutPIM->SetPoint(10,-52.02194,208.3333);
     //*************************************************
    const float oAngleCut=6.;
    vector<TLorentzVector> lv_neutr, lv_neutr1, lv_prot, lv_pip, lv_pim ;
    vector<TLorentzVector> lv_prot1, lv_prot2, lv_prot3, lv_prot4 ;
    vector<TLorentzVector> lv_gMix;

    //event Mixer settings
    HGenericEventMixerObj < TLorentzVector > eventmixer;
    eventmixer.setPIDs(1, 1, 7);
    eventmixer.setBuffSize(30);

    const double mp=938.27231;
    const double beam_energy = 4536. + mp;
    const double beam_momentum = sqrt(beam_energy*beam_energy-mp*mp);

    TLorentzVector proj(0,0,beam_momentum, beam_energy);
    TLorentzVector targ(0,0,0, mp);
    TLorentzVector beam(0,0,0,0);
    beam = proj + targ;

    //*************************************

    TStopwatch timer;
    timer.Reset();
    timer.Start();

      
    Long64_t nEvents = gLoop->getEntries();
    if (nEventsDesired < nEvents && nEventsDesired >= 0)
        nEvents = nEventsDesired;

    //****************************************
    // The global event loop which loops over all events in the DST files added to HLoop
    //********************************************

    for (Long64_t event = 0; event < nEvents; event++) {

      if (gLoop->nextEvent(event) <= 0) {
            cout << " Last event processed " << endl;
            break;
        }

	lv_neutr.clear();
	lv_prot1.clear();
	lv_prot2.clear();
	lv_prot3.clear();
	lv_prot4.clear();
	lv_prot.clear();
	lv_pip.clear();
	lv_pim.clear();
	lv_gMix.clear();


        HTool::printProgress(event, nEvents, 1, "Analyzed events: ");

	//Start - no iTOF
	HStart2Hit * fstart = nullptr;
	fstart = (HStart2Hit *) fStart2Hit->getObject(0);
	if (!fstart || fstart->getCorrFlag() == -1) continue;
	
	//selection of VertexZ  
        HEventHeader* event_header = gLoop->getEventHeader();
        if (event_header->getId() != 1)   continue;
	Int_t TBit   = (Int_t) event_header->getTBit();
        Double_t VertexZ = event_header->getVertexReco().getZ();
	if(VertexZ<-200 || VertexZ>-0) continue;
	
	//trigger bit info
 	int trigbit=-1;

	if((TBit&4096)==4096) trigbit=0;
	if((TBit&8192)==8192) trigbit=1;

	//*************************************************

	Int_t nNeutral_ev = fEmcNeutralCand->getEntries();

	for (int j = 0; j < nNeutral_ev; ++j)
	  {
	    HEmcNeutralCand* neutr_cand = HCategoryManager::getObject(neutr_cand, fEmcNeutralCand, j);

		Int_t ind=neutr_cand->getEmcClusterIndex();
		//Float_t dist  = neutr_cand->getDistanceToEmc();

		HEmcCluster *cl=nullptr;
		cl=HCategoryManager::getObject(cl, fEmcCluster, ind);
		Int_t cl_size = cl->getNCells();
		//Int_t sec = cl->getSector();
		//Float_t theta = cl->getTheta();
		//Float_t phi = cl->getPhi();

		Int_t cel = cl->getCell();		
		if(cel<33)continue;

		HGeomVector trackVec(cl->getXLab(), cl->getYLab(), cl->getZLab() - VertexZ);
		
		Double_t energy  = cl->getEnergy();
		Double_t tof  =cl->getTime();
		Double_t trackLength = trackVec.length();

		//Beta:
		Double_t beta = neutr_cand->getBeta();

		
		hcl_size->Fill(cl_size);
		hbeta->Fill(beta);
		
		if(trigbit==0)htof_PT2->Fill(tof);
		if(trigbit==1)htof_PT3->Fill(tof);
		htracklength->Fill(trackLength);
		hg_energy->Fill(energy);


		TLorentzVector lvg1;
		lvg1.SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);

		if (energy>150 && beta>0.8 && beta<1.2 && tof>0 && tof<100.){

		  lv_neutr.push_back(lvg1);	    
		  lv_gMix.push_back(lvg1);

		}

	  }//nNeutral_ev
	//*********************************************************
	if (fParticleCand)
	      {
		Int_t nPart_ev  = fParticleCand->getEntries();

		for (int j = 0; j < nPart_ev; ++j)
		  {
		    HParticleCand* fparticlecand = HCategoryManager::getObject(fparticlecand, fParticleCand, j);
		   
		    Float_t theta = fparticlecand->getTheta();
		    Float_t phi = fparticlecand->getPhi();
		    Float_t mom = fparticlecand->getMomentum();
		    Float_t beta = fparticlecand->getBeta();
		    //Float_t tof = fparticlecand->getTof();
		    Float_t charge=fparticlecand->getCharge(); 
		    Float_t mass= fparticlecand->getMass();


		    if(!fparticlecand->isFlagBit(kIsUsed))continue;


			hbeta_mom->Fill(charge*mom,beta);
			hmass_mom->Fill(charge*mom,mass);
			hmass->Fill(charge*mass);


			if(charge>0 &&  cutP->IsInside(charge*mom,mass) && beta < 1.){//PID protony
			  hbeta_mom_p->Fill(charge*mom,beta);

			  //part_pos.push_back(fparticlecand);
			  TLorentzVector lv1;
			  lv1=TLV_p(theta, phi, mom, 938.);

			  lv_prot.push_back(lv1);
			}			


			if(charge>0 &&  cutPIP->IsInside(charge*mom,mass) && beta < 1.0){//PID pi+
			  
			  TLorentzVector lv1a;
			  lv1a=TLV_p(theta, phi, mom, 140.);

			  lv_pip.push_back(lv1a);
			  hbeta_mom_pip->Fill(charge*mom,beta);
			  
			}			

			
			if(charge<0 && cutPIM->IsInside(charge*mom,mass) && beta < 1.0){//PID pi-
			  
			  TLorentzVector lv1b;
			  lv1b=TLV_p(theta, phi, mom, 140.);

			  lv_pim.push_back(lv1b);
			  hbeta_mom_pim->Fill(charge*mom,beta);

			}			
			
			
			//******************************************
					      
		  }//nPart_ev
    
	      }//fParticleCand

	//***********************************************
	//**************  gg  *****************
	//Inclusive:

	if(lv_neutr.size()>=2){
	    
	  for (long unsigned int ii=0;ii<lv_neutr.size();ii++){
	    for (long unsigned int jj=0;jj<lv_neutr.size();jj++){
	      if(ii==jj || ii>jj)continue;

	      float oAngle = lv_neutr[ii].Angle(lv_neutr[jj].Vect())*TMath::RadToDeg();
	      if (oAngle<oAngleCut) continue;
	      
	      TLorentzVector lvg2;
	      lvg2=lv_neutr[ii]+lv_neutr[jj];
	      float mass_gg=lvg2.M();
	      
	      //************** pi0 ************************
	      
	      if (trigbit==1) hinvM_gg_PT3->Fill(mass_gg);
	      
	      //**************** pip-pim-pi0  ********************
	      if(lv_pip.size() && lv_pim.size()){
		TLorentzVector lv3b;
		
		for (long unsigned int i=0;i<lv_pip.size();i++){
		  for (long unsigned int j=0;j<lv_pim.size();j++){
		    lv3b=lvg2+lv_pip[i]+lv_pim[j];
		    
		    if (trigbit==1) hMpippimpi0_PT3->Fill(lv3b.M());
		    
		  }
		}
	      }
	      //*****************************************
	    }
	  }
	}//lv_neutr.size()>=2

	//*****************************************
	//*************** MIX EVENTS  ******************	    

	int nNeutralCand=lv_gMix.size();
	
	eventmixer.nextEvent(nNeutralCand + 10*trigbit);
	eventmixer.addVector(lv_gMix, 1);
	vector < pair < TLorentzVector, TLorentzVector > >& pairsVec = eventmixer.getMixedVector();
	for (long unsigned int j = 0; j < pairsVec.size(); j++) {
	  pair < TLorentzVector, TLorentzVector >& pair = pairsVec[j];
	  TLorentzVector particleMix1 = pair.first;
	  TLorentzVector particleMix2 = pair.second;
	  
	  float oAngle = particleMix1.Angle(particleMix2.Vect())*TMath::RadToDeg();
	  if(oAngle< oAngleCut) continue;
	  
	  TLorentzVector lvMix;
	  lvMix=particleMix1+particleMix2;
	  
	  float mass_ggMix=lvMix.M();
	  
	  if(trigbit==1)hinvM_ggMix->Fill(mass_ggMix); //pi0
	  
	  //**********************************
	  
	  if(mass_ggMix>100. && mass_ggMix<200.){
	    
	    TLorentzVector lv4b;
	    if(lv_pip.size() && lv_pim.size()){
	      for (long unsigned int i=0;i<lv_pip.size();i++){
		for (long unsigned int j=0;j<lv_pim.size();j++){
		  lv4b=lvMix+lv_pip[i]+lv_pim[j];
		  if(trigbit==1)hMpippimpi0_Mix->Fill(lv4b.M());
		}
	      }
	    }
	    //***********************************		
	  }//if(mass_ggMix>100. && mass_ggMix<200.)
      
	}//End of event mixing

	
	//*******************************************************************
	//********* EXCLUSIVE ANALYSIS - pp selection ***********************
	
	if( lv_prot.size()>=2 ){
	  for (long unsigned int i=0;i<lv_prot.size();i++){
	    for (long unsigned int j=0;j<lv_prot.size();j++){
	      if(i==j || i>j)continue;
	      
	      
	      TLorentzVector lv_MMpp;
	      lv_MMpp=beam-lv_prot[i]-lv_prot[j];
	      
	      float mm=lv_MMpp.M();
	      float mm2=lv_MMpp.M2()*1e-6;
	     
	      if(trigbit==1){
		hMM_pp->Fill(mm);
		hMM2_pp->Fill(mm2);
	      }
	      

	      //protons from pi0  mass 
	      if(mm2>-0.1 && mm2<0.1){
		lv_prot1.push_back(lv_prot[i]);
		lv_prot2.push_back(lv_prot[j]);
	      }

	      //protons from omega  mass
	      if( mm2>0.5 && mm2<0.73){
		lv_prot3.push_back(lv_prot[i]);
		lv_prot4.push_back(lv_prot[j]);
	      }

	     
	    }
	  }
	}//lv_prot.size()>=2

	//*******************************************
	//*****************************************
	//exclusive: pi0 mass selected on MM2(pp) in Hades
	
	if(lv_neutr.size()>=2 && lv_prot1.size()==1 && lv_prot2.size()==1){
	  for (long unsigned int ii=0;ii<lv_neutr.size();ii++){
	    for (long unsigned int jj=0;jj<lv_neutr.size();jj++){
	      if(ii==jj || ii>jj)continue;
	      
	      float oAngle = lv_neutr[ii].Angle(lv_neutr[jj].Vect())*TMath::RadToDeg();
	      if (oAngle<oAngleCut) continue;
	      
	      TLorentzVector lvg2;
	      lvg2=lv_neutr[ii]+lv_neutr[jj];
	      float mass_gg=lvg2.M();
	      
	      if(trigbit==1) hinvM_gg_pi0_PT3->Fill(mass_gg);
	      
	    }
	  }
	}

	//**********************************************************
	//exclusive: omega selected in MM2(pp) in Hades
	
	if(lv_neutr.size()>=2 && lv_prot3.size()==1 && lv_prot4.size()==1 ){
	  for (long unsigned int ii=0;ii<lv_neutr.size();ii++){
	    for (long unsigned int jj=0;jj<lv_neutr.size();jj++){
	      if(ii==jj || ii>jj)continue;
	      
		  float oAngle = lv_neutr[ii].Angle(lv_neutr[jj].Vect())*TMath::RadToDeg();
		  if (oAngle<oAngleCut) continue;
		  
		  TLorentzVector lvg2;
		  lvg2=lv_neutr[ii]+lv_neutr[jj];
		  float mass_gg=lvg2.M();
		  if(trigbit==1) hinvM_gg_omega_pi0_PT3->Fill(mass_gg);
		  
		  for (long unsigned int i=0;i<lv_pip.size();i++){
		    for (long unsigned int j=0;j<lv_pim.size();j++){
		      TLorentzVector lv_om;
		      lv_om=lvg2+lv_pip[i]+lv_pim[j];
		      float mm_om=lv_om.M();
		      if(trigbit==1) hinvM_gg_omega_PT3->Fill(mm_om);
		    }
		  }
		  //*************************
		  	      		  
	    }
	  }
	}
	
	//***********************************************
	
    } // End of event loop

    
    
    timer.Stop();
    cout << "Finished DST processing" << endl;

   
    // Creating output file and storing results there
   
    TFile* out = new TFile(outFile.Data(), "RECREATE");
    out->cd();



    hbeta_mom->Write();
    hbeta_mom_p->Write();
    hbeta_mom_pip->Write();
    hbeta_mom_pim->Write();
    hmass_mom->Write();
    hmass->Write();

    htof_PT2->Write();
    htof_PT3->Write();
    hcl_size->Write();
    hbeta->Write();
    htracklength->Write();
    hg_energy->Write();

    hinvM_gg_PT3->Write();
    hinvM_ggMix->Write();

    //hMM_pp->Write();
    hMM2_pp->Write();
    hinvM_gg_pi0_PT3->Write();

    
    hMpippimpi0_PT3->Write();
    hMpippimpi0_Mix->Write();

    hinvM_gg_omega_PT3->Write();
    hinvM_gg_omega_pi0_PT3->Write();;

    
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
    return 0;
}
