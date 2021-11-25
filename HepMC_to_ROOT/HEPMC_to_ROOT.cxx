/*
	HEPMC_to_ROOT parser
	Version 1.0 (only final particles info)
	Vladislav Sandul, 2021
	vladislav.sandul@cern.ch

	To run the script, one have to modify filename_path and filename_body
	according to their needs, then type in terminal
		$ root "HEPMC_to_ROOT.cxx(FILENUM)"
	where FILENUM is a number of the file.

	Output file will be save with the same name as input file, but with the
	.root format.
*/

void HEPMC_to_ROOT(int filenum)
{
	TString filename_path = "./";
	string filename_body = Form("crmc_epos199_300651974_A197Z79_p_%d", filenum);
	TString input_filename = Form("%s.hepmc", filename_body.c_str());
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
		//int charge[MAX_N_TRACKS];
		int pdg[MAX_N_TRACKS];
		float px[MAX_N_TRACKS];
		float py[MAX_N_TRACKS];
		float pz[MAX_N_TRACKS];
		float p0[MAX_N_TRACKS];
		float m[MAX_N_TRACKS];

	tree_final -> Branch("nTracks",  &ntracks , "nTracks/I");
	tree_final -> Branch("Impact",   &b ,				"Impact/F");
	//tree_final -> Branch("Charge",   charge  ,	"Charge[nTracks]/I");
	tree_final -> Branch("PDGID",    pdg    ,		"PDGID[nTracks]/I");
	tree_final -> Branch("Px",       px      , 	"Px[nTracks]/F");
	tree_final -> Branch("Py",     	 py      ,	"Py[nTracks]/F");
	tree_final -> Branch("Pz",       pz      , 	"Pz[nTracks]/F");
	tree_final -> Branch("P0",    	 p0    , 		"P0[nTracks]/F");
	tree_final -> Branch("M",      	 m       , 	"M[nTracks]/F");

	  // ---> Define event variables to be read from file - for tree_final
		char format_char = '0';

		// Event string variables
    int event_number_ = 0, mpi_=0, signal_process_id_ = 0, signal_process_vertex_ = 0,
				vertexes_size_ = 0, beam_particle1_ = 0, beam_particle2_ = 0,
				random_states_size_ = 0, event_weights_size_ = 0;
		float event_scale_ = 0., alphaQCD_ = 0., alphaQED_ = 0.;

		// Heavy-ion string variable
		int ncoll_hard_ = 0, npart_proj_ = 0., npart_targ_ = 0, ncoll_ = 0,
				spec_neutr_ = 0, spec_prot_ = 0, n_nwounded_col_ = 0, nwounded_n_col_ = 0;
		float impact_= 0., event_plane_angle_ = 0., eccentricity_ = 0, sigma_inel_NN_ = 0;

		// Vertex string variables
		int vertex_id_ = 0, id_ = 0, num_in_ = 0, num_out_ = 0, vertex_weights_size_ = 0;
		float pos_x_ = 0., pos_y_=0., pos_z_ = 0., pos_t_ = 0.;

		// Particle string variables
    int track_number_ = 0, pdgid_ = 0, status_ = 0, polar_theta_ = 0,
		 		polar_phi_ = 0, end_vertex_ = 0, flow_ = 0;
    float  px_=0., py_=0., pz_=0., p0_=0., m_=0.;

		// Canvas to save information about input file, colliding systems
		// and energy of collisions
  TCanvas *canv_info = new TCanvas("Info");
 	TPaveText *pt = new TPaveText(.05,.1,.95,.8);

	string s;
	string word;
	int n_event_counter = 0;
	// Reading from the input file line-by-line
	if (fInputFile ) {
		getline(fInputFile, s);
		pt->AddText(input_filename);
		getline(fInputFile, s);
		pt->AddText(s.c_str());
		getline(fInputFile, s);
		pt->AddText(s.c_str());
   	pt->Draw();  // to add info in the output file

		// before going through the events we want to check up the correctness of file
		// assuming that the 3rd string of the .hepmc input file  must have a certin view
		if(s != "HepMC::IO_GenEvent-START_EVENT_LISTING"){
			cout << "ERROR! Wrong format of the input file. Please, check it up." << endl;
			return;
		}

		// Event-by-event loop
 		while(getline(fInputFile, s)){
			if (s == "HepMC::IO_GenEvent-END_EVENT_LISTING"){
				break;
			}
			int ntracks_ = 0;
			// Read event info
			stringstream ss(s);
			ss >> format_char >> event_number_ >> mpi_ >> event_scale_ >> alphaQCD_ >>
						alphaQED_ >> signal_process_id_ >> signal_process_vertex_ >>
						vertexes_size_ >> beam_particle1_ >> beam_particle2_ >>
						random_states_size_ >> event_weights_size_;
			if (format_char != 'E'){
				cout << "Error! format_char != 'E'" << endl;
				return;
			}
			if (n_event_counter%1 == 0){
				cout << "processing  " << n_event_counter << " event ..."<< "\r";
				return;
			}
			cout.flush();
			n_event_counter++;//= nevt_;
			int final_tracks_counter = 0;
			ss.clear();
			ss.str(std::string());

			getline(fInputFile, s); // skip 'U' flaged string
			getline(fInputFile, s); // skip 'C' flaged string

			// Read Heavy-ion info
			getline(fInputFile, s);
			ss.str(s);
			ss >> format_char >> ncoll_hard_ >> npart_proj_ >> npart_targ_ >> ncoll_ >>
			 			spec_neutr_ >> spec_prot_ >> n_nwounded_col_ >> nwounded_n_col_ >>
						impact_ >> event_plane_angle_ >> eccentricity_ >> sigma_inel_NN_;
			if (format_char != 'H'){
				cout << "Error! format_char != 'H'" << endl;
				return;
			}
			b = impact_;
			ss.clear();
			ss.str(std::string());

			getline(fInputFile, s); // skip 'F' flaged string

			for (int ivert = 0; ivert < vertexes_size_; ivert++){
					// Read vertex
					getline(fInputFile, s);
					ss.str(s);
					ss >> format_char >> vertex_id_ >> id_ >> pos_x_ >> pos_y_ >> pos_z_ >>
					 			pos_t_ >> num_in_ >> num_out_ >> vertex_weights_size_;
					if (format_char != 'V'){
						cout << "Error! format_char != 'V'" << endl;
						return;
					}
					ss.clear();
					ss.str(std::string());

					for (int ipart = 0; ipart < num_in_+num_out_; ipart++){
						  // Read particle
							getline(fInputFile, s);
							ss.str(s);
							ss >> format_char >> track_number_ >> pdgid_ >> px_ >> py_ >> pz_ >>
							 			p0_ >> m_ >> status_ >> polar_theta_ >> polar_phi_ >>
										end_vertex_ >> flow_;
							if (format_char != 'P'){
								cout << "Error! format_char != 'P'" << endl;
								return;
							}
							if (status_ == 1){
																//cout << s << endl;
								final_tracks_counter++;
								pdg[final_tracks_counter] = pdgid_;
								px[final_tracks_counter] = px_;
								py[final_tracks_counter] = py_;
								pz[final_tracks_counter] = pz_;
								p0[final_tracks_counter] = p0_;
								m[final_tracks_counter] = m_;
							}
							ss.clear();
							ss.str(std::string());
					}
			}

			ntracks = final_tracks_counter;
			tree_final -> Fill();
			cout << ntracks << endl;

		}

	}
	cout << "\nFinished!" << endl;

	pt->AddText(Form("%d events", n_event_counter));

	fInputFile.close();

	canv_info->Write();
	tree_final -> Write();  // If in your output file you have no tree - uncomment this line!

}
