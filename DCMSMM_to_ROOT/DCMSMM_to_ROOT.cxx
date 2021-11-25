/*
	DCMSMM_to_ROOT parser
	Vladislav Sandul, 2021
	vladislav.sandul@cern.ch

	To run the script, one have to modify filename_path and filename_body
	according to their needs, then type in terminal
		$ root "DCMSMM_to_ROOT.cxx(FILENUM)"
	where FILENUM is a number of the file and NEVENTS is a number of events in
	the input file.

	Output file will be save with the same name as input file, but with the
	.root format.
*/

void DCMSMM_to_ROOT(int filenum)
{
	TString filename_path = "./";
	string filename_body = Form("DCMSMM_BiBi_ss9.2GeV_mb_2k_%d", filenum);
	TString input_filename = Form("%s.r12", filename_body.c_str());
	TString output_filename = Form("%s.root", filename_body.c_str());

	auto fullpath = 	filename_path + input_filename;
	ifstream fInputFile(fullpath);
	cout << "Input file: " << fullpath << endl;
	if(!fInputFile){
		cout << "Cannot open the file: no such file. Exit." << endl;
		return;
	}

	//  ROOT file to write histograms or trees
	unique_ptr<TFile> out_file( TFile::Open(filename_path+output_filename, "RECREATE") );
	auto tree_final = make_unique<TTree>("FinalState", "FinalState");
	auto tree_proj_res = make_unique<TTree>("ProjectileResidual", "ProjectileResidual");
	auto tree_targ_res = make_unique<TTree>("TargetResidual", "TargetResidual");

		// ##### prepare tree_final
		const int MAX_N_TRACKS = 10000;
		int ntracks ;
		float b ;
		float bx ;
		float by ;
		int charge[MAX_N_TRACKS];
		int lepton_num[MAX_N_TRACKS];
		int strangeness[MAX_N_TRACKS];
		int baryon_num[MAX_N_TRACKS];
		int pdg[MAX_N_TRACKS];
		float px[MAX_N_TRACKS];
		float py[MAX_N_TRACKS];
		float pz[MAX_N_TRACKS];
		float pzlab[MAX_N_TRACKS];
		float m[MAX_N_TRACKS];

	tree_final -> Branch("nTracks",  &ntracks , "nTracks/I");
	tree_final -> Branch("Impact",   &b ,				"Impact/F");
	tree_final -> Branch("Impact_x", &bx ,			"Impact_x/F");
	tree_final -> Branch("Impact_y", &by ,			"Impact_y/F");
	tree_final -> Branch("Charge",   charge  ,	"Charge[nTracks]/I");
	tree_final -> Branch("PDGID",    pdg    ,		"PDGID[nTracks]/I");
	tree_final -> Branch("Px",       px      , 	"Px[nTracks]/F");
	tree_final -> Branch("Py",     	 py      ,	"Py[nTracks]/F");
	tree_final -> Branch("Pz",       pz      , 	"Pz[nTracks]/F");
	tree_final -> Branch("Pzlab",    pzlab    , "Pzlab[nTracks]/F");
	tree_final -> Branch("M",      	 m       , 	"M[nTracks]/F");


		// ##### prepare tree_proj_res
		const int MAX_N_FRAG = 1000;
		int proj_num_of_frags;
		float proj_atomic_num;
		float proj_charge;
		float proj_strangeness;
		float proj_excit_energy;
		float proj_px;
		float proj_py;
		float proj_pz;
		int charge_proj[MAX_N_FRAG];
		int lepton_num_proj[MAX_N_FRAG];
		int strangeness_proj[MAX_N_FRAG];
		int baryon_num_proj[MAX_N_FRAG];
		int pdg_proj[MAX_N_FRAG];
		float px_proj[MAX_N_FRAG];
		float py_proj[MAX_N_FRAG];
		float pz_proj[MAX_N_FRAG];
		float pzlab_proj[MAX_N_FRAG];
		float m_proj[MAX_N_FRAG];

	tree_proj_res -> Branch("nFrags",  					 &proj_num_of_frags,	"nFrags/I");
	tree_proj_res -> Branch("ProjAtomicNum",  	 &proj_atomic_num , 	"ProjAtomicNum/F");
	tree_proj_res -> Branch("ProjCharge",  			 &proj_charge ,				"ProjCharge/F");
	tree_proj_res -> Branch("ProjStrangeness",   &proj_strangeness ,	"ProjStrangeness/F");
	tree_proj_res -> Branch("ProjExcitEnergy",   &proj_excit_energy ,	"ProjExcitEnergy/F");
	tree_proj_res -> Branch("ProjPx",   					&proj_px ,					"ProjPx/F");
	tree_proj_res -> Branch("ProjPy",   					&proj_py ,					"ProjPy/F");
	tree_proj_res -> Branch("ProjPz",   					&proj_pz ,					"ProjPz/F");
	tree_proj_res -> Branch("Charge",  					 charge_proj  ,				"Charge[nFrags]/I");
	tree_proj_res -> Branch("PDGID",   					 pdg_proj    ,				"PDGID[nFrags]/I");
	tree_proj_res -> Branch("Px",       				 	px_proj      , 			"Px[nFrags]/F");
	tree_proj_res -> Branch("Py",     	   			  py_proj      ,			 "Py[nFrags]/F");
	tree_proj_res -> Branch("Pz",      					 pz_proj      , 			"Pz[nFrags]/F");
	tree_proj_res -> Branch("Pzlab",        			pzlab_proj     , 		"Pzlab[nFrags]/F");
	tree_proj_res -> Branch("M",      		 			 m_proj       , 			"M[nFrags]/F");

		// ##### prepare tree_proj_res
		int targ_num_of_frags;
		float targ_atomic_num;
		float targ_charge;
		float targ_strangeness;
		float targ_excit_energy;
		float targ_px;
		float targ_py;
		float targ_pz;
		int charge_targ[MAX_N_FRAG];
		int lepton_num_targ[MAX_N_FRAG];
		int strangeness_targ[MAX_N_FRAG];
		int baryon_num_targ[MAX_N_FRAG];
		int pdg_targ[MAX_N_FRAG];
		float px_targ[MAX_N_FRAG];
		float py_targ[MAX_N_FRAG];
		float pz_targ[MAX_N_FRAG];
		float pzlab_targ[MAX_N_FRAG];
		float m_targ[MAX_N_FRAG];

	tree_targ_res -> Branch("nFrags",   				&targ_num_of_frags ,	"nFrags/I");
	tree_targ_res -> Branch("TargAtomicNum",  	 &targ_atomic_num ,		"TargAtomicNum/F");
	tree_targ_res -> Branch("TargCharge",   		&targ_charge ,				"TargCharge/F");
	tree_targ_res -> Branch("TargStrangeness",   &targ_strangeness ,	"TargStrangeness/F");
	tree_targ_res -> Branch("TargExcitEnergy",   &targ_excit_energy ,	"TargExcitEnergy/F");
	tree_targ_res -> Branch("TargPx",   				&targ_px ,						"TargPx/F");
	tree_targ_res -> Branch("TargPy",   				&targ_py ,						"TargPy/F");
	tree_targ_res -> Branch("TargPz",   				&targ_pz ,						"TargPz/F");
	tree_targ_res -> Branch("Charge",  					 charge_targ  ,				"Charge[nFrags]/I");
	tree_targ_res -> Branch("PDGID",   					 pdg_targ    ,				"PDGID[nFrags]/I");
	tree_targ_res -> Branch("Px",       				 	px_targ      , 			"Px[nFrags]/F");
	tree_targ_res -> Branch("Py",     	    		 py_targ      ,				 "Py[nFrags]/F");
	tree_targ_res -> Branch("Pz",      					 pz_targ      , 			"Pz[nFrags]/F");
	tree_targ_res -> Branch("Pzlab",        		pzlab_targ      , 		"Pzlab[nFrags]/F");
	tree_targ_res -> Branch("M",      		  		m_targ       , 				"M[nFrags]/F");

	  // ---> Define event variables to be read from file - for tree_final
    int nevt_ = 0, ntracks_=0;
    float b_= 0., bx_ = 0., by_ = 0.;
    int chg_=0, lepnum_=0, stnum_=0, barnum_=0, pdgid_=0;
    float  px_=0., py_=0., pz_=0., pzlab_=0., m_=0.;

		// ---> Define event variables to be read from file - for tree_proj_res
		int nfrags_proj_ = 0;
		float proj_atomnum_ = 0., proj_chrg_ = 0.,  proj_st_ = 0.,
		 		proj_ex_enrg_ = 0., proj_px_ = 0., proj_py_ = 0., proj_pz_ =0.;
    int chg_proj_=0, lepnum_proj_=0, stnum_proj_=0, barnum_proj_=0, pdgid_proj_=0;
    float  px_proj_=0., py_proj_=0., pz_proj_=0., pzlab_proj_=0., m_proj_=0.;

		// ---> Define event variables to be read from file - for tree_targ_res
		int nfrags_targ_ = 0;
		float targ_atomnum_ = 0., targ_chrg_ = 0.,  targ_st_ = 0.,
		 		targ_ex_enrg_ = 0., targ_px_ = 0., targ_py_ = 0., targ_pz_ =0.;
		int chg_targ_=0, lepnum_targ_=0, stnum_targ_=0, barnum_targ_=0, pdgid_targ_=0;
		float  px_targ_=0., py_targ_=0., pz_targ_=0., pzlab_targ_=0., m_targ_=0.;

		// Canvas to save information about input file, colliding systems
		// and energy of collisions
  TCanvas *canv_info = new TCanvas("Info");
 	TPaveText *pt = new TPaveText(.05,.1,.95,.8);

	string s;
	string word;
	int n_event_counter = 0;
	// Reading from the input file line-by-line
	if (fInputFile ) {
		pt->AddText(input_filename);
		getline(fInputFile, s);
		pt->AddText(s.c_str());
		getline(fInputFile, s);
		pt->AddText(s.c_str());
   	pt->Draw();  // to add info in the output file

 		for (int i = 0; i < 16; i++){
 			getline(fInputFile, s); // to skip some initial lines in header
 		}
		getline(fInputFile, s);				// before going through the events
		stringstream ss_to_check(s); // we want to check up the correctness of file
		ss_to_check >> word;				// assuming that the 19th string of the .r12 input file
		if(word != "clusters"){			// must start with the word "clusters"
			cout << "ERROR! Wrong format of the input file. Please, check it up." << endl;
			return;
		}
		getline(fInputFile, s);

		// Event-by-event loop
 		while(getline(fInputFile, s)){
			stringstream ss(s);
			ss >> nevt_ >> ntracks_ >> b_ >> bx_ >> by_; //event info
			ntracks = ntracks_;
			b = b_;
			bx = bx_;
			by = by_;
			if (nevt_%1 == 0)
				cout << "processing  " << nevt_ <<  " event ..."<< "\r";
			cout.flush();
			n_event_counter++;
			ss.clear();
			ss.str(std::string());

			getline(fInputFile,s);
			ss.str(s);
			ss >> nfrags_proj_ >> proj_atomnum_ >> // Projectile residual nucleus info
					proj_chrg_ >>  proj_st_ >> proj_ex_enrg_ >>
					proj_px_ >> proj_py_ >> proj_pz_;
			proj_num_of_frags = nfrags_proj_;
			proj_atomic_num = proj_atomnum_;
			proj_charge = proj_chrg_;
			proj_strangeness = proj_st_;
			proj_excit_energy = proj_ex_enrg_;
			proj_px = proj_px_;
			proj_py = proj_py_;
			proj_pz = proj_pz_;
			ss.clear();
 			ss.str(std::string());
			for(int ifrag = 0; ifrag < nfrags_proj_; ifrag++){ // fragments loop
				getline(fInputFile,s);
				ss.str(s);
				ss >> chg_proj_ >> lepnum_proj_ >> stnum_proj_ >>
						barnum_proj_ >>	pdgid_proj_ >> px_proj_ >>
						py_proj_ >> pz_proj_ >> pzlab_proj_ >> m_proj_;

				charge_proj[ifrag] = chg_proj_;
			 	lepton_num_proj[ifrag] = lepnum_proj_;
			 	strangeness_proj[ifrag] = stnum_proj_;
			 	baryon_num_proj[ifrag] = barnum_proj_;
			 	pdg_proj[ifrag] = pdgid_proj_;
			  px_proj[ifrag] = px_proj_;
			  py_proj[ifrag] = py_proj_;
			 	pz_proj[ifrag] = pz_proj_;
			 	pzlab_proj[ifrag] = pzlab_proj_;
			 	m_proj[ifrag] = m_proj_;

				ss.clear();
				ss.str(std::string());
			}


			getline(fInputFile,s);
			ss.str(s);
			ss >> nfrags_targ_ >> targ_atomnum_ >> // Target residual nucleus info
					targ_chrg_ >>  targ_st_ >>	targ_ex_enrg_ >>
					targ_px_ >> targ_py_ >> targ_pz_;
			targ_num_of_frags = nfrags_targ_;
			targ_atomic_num = targ_atomnum_;
			targ_charge = targ_chrg_;
			targ_strangeness = targ_st_;
			targ_excit_energy = targ_ex_enrg_;
			targ_px = targ_px_;
			targ_py = targ_py_;
			targ_pz = targ_pz_;
			ss.clear();
 			ss.str(std::string());
			for(int ifrag = 0; ifrag < nfrags_targ_; ifrag++){ // fragments loop
				getline(fInputFile,s);
				ss.str(s);
				ss >> chg_targ_ >> lepnum_targ_ >> stnum_targ_ >>
				 		barnum_targ_ >>	pdgid_targ_ >> px_targ_ >>
						py_targ_ >> pz_targ_ >> pzlab_targ_ >> m_targ_;

				charge_targ[ifrag] = chg_targ_;
			 	lepton_num_targ[ifrag] = lepnum_targ_;
			 	strangeness_targ[ifrag] = stnum_targ_;
			 	baryon_num_targ[ifrag] = barnum_targ_;
			 	pdg_targ[ifrag] = pdgid_targ_;
			  px_targ[ifrag] = px_targ_;
			  py_targ[ifrag] = py_targ_;
			 	pz_targ[ifrag] = pz_targ_;
			 	pzlab_targ[ifrag] = pzlab_targ_;
			 	m_targ[ifrag] = m_targ_;

				ss.clear();
				ss.str(std::string());
			}

			// Characteristics of particles after cascade and light
  		// clusters after coalescence stages
			getline(fInputFile,s);
			for (int itrack = 0; itrack<ntracks_; itrack++){
				getline(fInputFile,s);
				ss.str(s);
				ss >> chg_ >> lepnum_ >> stnum_ >> barnum_ >>
				 		pdgid_ >> px_ >> py_ >> pz_ >> pzlab_ >> m_;

				charge[itrack] = chg_;
			 	lepton_num[itrack] = lepnum_;
			 	strangeness[itrack] = stnum_;
			 	baryon_num[itrack] = barnum_;
			 	pdg[itrack] = pdgid_;
			  px[itrack] = px_;
			  py[itrack] = py_;
			 	pz[itrack] = pz_;
			 	pzlab[itrack] = pzlab_;
			 	m[itrack] = m_;

				ss.clear();
				ss.str(std::string());
			}

			tree_final -> Fill();
			tree_proj_res -> Fill();
			tree_targ_res -> Fill();

		}

	}
	cout << "\nFinished!" << endl;

	pt->AddText(Form("%d events", n_event_counter));

	fInputFile.close();

	canv_info->Write();
	tree_final -> Write();
	tree_proj_res -> Write();
	tree_targ_res -> Write();

}
