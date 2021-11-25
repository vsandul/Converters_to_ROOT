/*
	DQGSM_to_ROOT parser
	Vladislav Sandul, 2021
	vladislav.sandul@cern.ch

	To run the script, one have to modify filename_path and filename_body
	according to their needs, then type in terminal
		$ root "DQGSM_to_ROOT.cxx(FILENUM)"
	where FILENUM is a number of the file .
	Output file will be save with the same name as input file, but with the
	.root format.
*/

void DQGSM_to_ROOT(int filenum)
{
	TString filename_path = "./";
	string filename_body = Form("DQGSM_AuAu_11_mb_2k_%d", filenum);
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
		float pzalab[MAX_N_TRACKS];
		float m[MAX_N_TRACKS];

	tree_final -> Branch("nTracks",  &ntracks , "nTracks/I");
	tree_final -> Branch("Impact",   &b ,				"Impact/F");
	tree_final -> Branch("Impact_x", &bx ,			"Impact_x/F");
	tree_final -> Branch("Impact_y", &by ,			"Impact_y/F");
	tree_final -> Branch("Charge",   charge  ,	"Charge[nTracks]/I");
	tree_final -> Branch("Baryon_num",   baryon_num  ,	"Baryon_num[nTracks]/I");
	tree_final -> Branch("Lepton_num",   lepton_num  ,	"Lepton_num[nTracks]/I");
	tree_final -> Branch("Strangeness",   strangeness  ,	"Strangeness[nTracks]/I");
	tree_final -> Branch("PDGID",    pdg    ,		"PDGID[nTracks]/I");
	tree_final -> Branch("Px",       px      , 	"Px[nTracks]/F");
	tree_final -> Branch("Py",     	 py      ,	"Py[nTracks]/F");
	tree_final -> Branch("Pz",       pz      , 	"Pz[nTracks]/F");
	tree_final -> Branch("Pzlab",    pzlab    , "Pzlab[nTracks]/F");
	tree_final -> Branch("Pzalab",    pzalab  , "Pzalab[nTracks]/F");
	tree_final -> Branch("M",      	 m       , 	"M[nTracks]/F");

	  // ---> Define event variables to be read from file - for tree_final
    int nevt_ = 0, ntracks_=0;
    float b_= 0., bx_ = 0., by_ = 0.;
    int chg_=0, lepnum_=0, stnum_=0, barnum_=0, pdgid_=0;
    float  px_=0., py_=0., pz_=0., pzlab_=0., pzalab_=0., m_=0.;

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
		//pt->AddText(Form("%d events", n_events_in_file));
		getline(fInputFile, s);
		pt->AddText(s.c_str());
		getline(fInputFile, s);
		pt->AddText(s.c_str());
   	pt->Draw();  // to add info in the output file


		getline(fInputFile, s);								// before going through the events
		stringstream ss_to_check(s); 					// we want to check up the correctness of file
		ss_to_check >> word >> word >> word;	// assuming that the 3rd string of the .r12 input file
		if(word != "ID=4"){										// must contain the word "ID=4" at the end.
			cout << "ERROR! Wrong format of the input file. Please, check it up." << endl;
			return;
		}
		getline(fInputFile, s);
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
				cout << "processing  " << nevt_ << " event ..."<< "\r";
			cout.flush();
			n_event_counter++;//= nevt_;
			ss.clear();
			ss.str(std::string());

			// Characteristics of particles after cascade and light
  		// clusters after coalescence stages
			for (int itrack = 0; itrack<ntracks_; itrack++){
				getline(fInputFile,s);
				ss.str(s);
				ss >> chg_ >> lepnum_ >> stnum_ >> barnum_ >>
				 		pdgid_ >> px_ >> py_ >> pz_ >> pzlab_ >> pzalab_ >> m_;

				charge[itrack] = chg_;
			 	lepton_num[itrack] = lepnum_;
			 	strangeness[itrack] = stnum_;
			 	baryon_num[itrack] = barnum_;
			 	pdg[itrack] = pdgid_;
			  px[itrack] = px_;
			  py[itrack] = py_;
			 	pz[itrack] = pz_;
			 	pzlab[itrack] = pzlab_;
				pzalab[itrack] = pzalab_;
			 	m[itrack] = m_;

				ss.clear();
				ss.str(std::string());
			}

			tree_final -> Fill();

		}

	}
	cout << "\nFinished!" << endl;

	pt->AddText(Form("%d events", n_event_counter));

	fInputFile.close();

	canv_info->Write();
	//tree_final -> Write();  // If in your output file you have no tree - uncomment this line!

}
