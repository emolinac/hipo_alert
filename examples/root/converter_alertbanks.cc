//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************
//* Example program for reading HIPO-4 Files..
//* Reads the file and converts the file to ROOT.
//*--
//* Author: G.Gavalian
//*

#include <cstdlib>
#include <iostream>
#include "reader.h"
#include "writer.h"
#include "utils.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "Compression.h"

#define maxhits        10000
#define maxpreclusters  1000
#define maxclusters     1000
#define maxmctracks      100
#define maxtracks        100
#define maxkftracks      100
#define maxahdcadc       500 // TODO: check this eventually!
#define maxrectracks    1000


void    convert(const char *hipoFile, bool do_mc);
double  benchmarkRoot(const char *file);
double  benchmarkHipo(const char *file);

int main(int argc, char **argv)
{
	std::cout << "reading file example program (HIPO) "
		  << __cplusplus << std::endl;

	char inputFile[512];
	// char outputFile[512];
	// char outputFileHipo[512];
	bool do_mc = false;

	if (argc > 2)
		do_mc = argv[2];

	if (argc > 1) {
		sprintf(inputFile, "%s", argv[1]);
		// sprintf(outputFile, "%s.root", argv[1]);
		// sprintf(outputFileHipo, "%s_writer.hipo", argv[1]);
		// sprintf(outputFile,"%s",argv[2]);
	} else {
		std::cout << " *** please provide a file name (and an optional flag to add MC)..." << std::endl;
		exit(1);
	}

	std::cout << inputFile << " " << do_mc << std::endl;
	convert(inputFile, do_mc);
}

void convert(const char *hipoFile, bool do_mc)
{
	int compression = 404;
	char outputFile[512];
	char outputFileHipo[512];

	sprintf(outputFile, "%s_w.root", hipoFile);
	sprintf(outputFileHipo, "%s_w.hipo", hipoFile);

	auto fileOutput = new TFile(outputFile, "RECREATE");
	fileOutput->SetCompressionSettings(compression);
	// fileOutput->SetCompressionLevel(ROOT::RCompressionSetting::ELevel::kDefaultLZ4);

	ROOT::TIOFeatures features;
	features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);

	auto treeOutput = new TTree("clas12", "");
	// if (iofeatures)
	//	treeOutput->SetIOFeatures(features);

	/* AHDC hits */
	int nhits;
	int hitid[maxhits];
	int hitlayer[maxhits];
	int hitsuperlayer[maxhits];
	int hitwire[maxhits];
	float hitdoca[maxhits];

	treeOutput->Branch("nhits", &nhits, "nhits/I", 512000);
	treeOutput->Branch("hitid", &hitid, "hitid[nhits]/I", 512000);
	treeOutput->Branch("hitlayer", &hitlayer, "hitlayer[nhits]/I", 512000);
	treeOutput->Branch("hitsuperlayer", &hitsuperlayer, "hitsuperlayer[nhits]/I", 512000);
	treeOutput->Branch("hitwire", &hitwire, "hitwire[nhits]/I", 512000);
	treeOutput->Branch("hitdoca", &hitdoca, "hitdoca[nhits]/F", 512000);

	/* AHDC preclusters */
	int npreclusters;
	float preclusterX[maxpreclusters];
	float preclusterY[maxpreclusters];

	treeOutput->Branch("npreclusters", &npreclusters, "npreclusters/I", 512000);
	treeOutput->Branch("preclusterX", &preclusterX, "preclusterX[npreclusters]/F", 512000);
	treeOutput->Branch("preclusterY", &preclusterY, "preclusterY[npreclusters]/F", 512000);

	/* AHDC clusters */
	int nclusters;
	float clusterX[maxclusters];
	float clusterY[maxclusters];
	float clusterZ[maxclusters];

	treeOutput->Branch("nclusters", &nclusters, "nclusters/I", 512000);
	treeOutput->Branch("clusterX", &clusterX, "clusterX[nclusters]/F", 512000);
	treeOutput->Branch("clusterY", &clusterY, "clusterY[nclusters]/F", 512000);
	treeOutput->Branch("clusterZ", &clusterZ, "clusterZ[nclusters]/F", 512000);

	// /* MC tracks */
	// int nmctracks;
	// int mctrackpid[maxmctracks];
	// float mctrackpx[maxmctracks];
	// float mctrackpy[maxmctracks];
	// float mctrackpz[maxmctracks];
	// float mctrackvx[maxmctracks];
	// float mctrackvy[maxmctracks];
	// float mctrackvz[maxmctracks];
	// float mctrackvt[maxmctracks];
	// if (do_mc) {
	// 	treeOutput->Branch("nmctracks", &nmctracks, "nmctracks/I", 512000);
	// 	treeOutput->Branch("mctrackpx", &mctrackpx, "mctrackpx[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackpy", &mctrackpy, "mctrackpy[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackpz", &mctrackpz, "mctrackpz[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackvx", &mctrackvx, "mctrackvx[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackvy", &mctrackvy, "mctrackvy[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackvz", &mctrackvz, "mctrackvz[nmctracks]/F", 512000);
	// 	treeOutput->Branch("mctrackvt", &mctrackvt, "mctrackvt[nmctracks]/F", 512000);
	// }

	/* AHDC tracks */
	int ntracks;
	float trackx[maxtracks];
	float tracky[maxtracks];
	float trackz[maxtracks];
	float trackpx[maxtracks];
	float trackpy[maxtracks];
	float trackpz[maxtracks];

	treeOutput->Branch("ahdc_ntracks", &ntracks, "ntracks/I",          512000); 
	treeOutput->Branch("ahdc_trackx",  &trackx,  "trackx[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_tracky",  &tracky,  "tracky[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_trackz",  &trackz,  "trackz[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_trackpx", &trackpx, "trackpx[ntracks]/F", 512000);
	treeOutput->Branch("ahdc_trackpy", &trackpy, "trackpy[ntracks]/F", 512000);
	treeOutput->Branch("ahdc_trackpz", &trackpz, "trackpz[ntracks]/F", 512000);

	/* AHDC KF tracks */
	int nkftracks;
	int kftrack_id[maxkftracks];
	float kftrackx[maxkftracks];
	float kftracky[maxkftracks];
	float kftrackz[maxkftracks];
	float kftrackpx[maxkftracks];
	float kftrackpy[maxkftracks];
	float kftrackpz[maxkftracks];
	float kftrackpath[maxkftracks];
	float kftrackdedx[maxkftracks];

	treeOutput->Branch("ahdc_nkftracks",   &nkftracks,   "nkftracks/I", 		  512000);
	treeOutput->Branch("ahdc_kftrack_id",  &kftrack_id,  "kftrack_id[nkftracks]/I",   512000);
	treeOutput->Branch("ahdc_kftrackx",    &kftrackx,    "kftrackx[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftracky",    &kftracky,    "kftracky[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftrackz",    &kftrackz,    "kftrackz[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftrackpx",   &kftrackpx,   "kftrackpx[nkftracks]/F",    512000);
	treeOutput->Branch("ahdc_kftrackpy",   &kftrackpy,   "kftrackpy[nkftracks]/F",    512000);
	treeOutput->Branch("ahdc_kftrackpz",   &kftrackpz,   "kftrackpz[nkftracks]/F",    512000);
	treeOutput->Branch("ahdc_kftrackpath", &kftrackpath, "kftrackpath[nkftracks]/F",  512000);
	treeOutput->Branch("ahdc_kftrackdedx", &kftrackdedx, "kftrackdedx[nkftracks]/F",  512000);

	/* AHDC ADC information */
	int   n_ahdcrows;
	int   AHDC_sector[maxahdcadc];
	int   AHDC_component[maxahdcadc];
	int   AHDC_order[maxahdcadc];
	int   AHDC_maxADC[maxahdcadc];
	float AHDC_leadingEdgeTime[maxahdcadc];
	float AHDC_timeOverThreshold[maxahdcadc];

	treeOutput->Branch("ahdc_nadcrows",          &n_ahdcrows,             "n_ahdcrows/I", 			       512000);
	treeOutput->Branch("ahdc_sector",            &AHDC_sector,            "AHDC_sector[n_ahdcrows]/I",             512000);
	treeOutput->Branch("ahdc_component",         &AHDC_component,         "AHDC_component[n_ahdcrows]/I",          512000);
	treeOutput->Branch("ahdc_order",             &AHDC_order,             "AHDC_order[n_ahdcrows]/I",              512000);
	treeOutput->Branch("ahdc_maxADC",            &AHDC_maxADC,            "AHDC_maxADC[n_ahdcrows]/I",             512000);
	treeOutput->Branch("ahdc_leadingEdgeTime",   &AHDC_leadingEdgeTime,   "AHDC_leadingEdgeTime[n_ahdcrows]/F",    512000);
	treeOutput->Branch("ahdc_timeOverThreshold", &AHDC_timeOverThreshold, "AHDC_timeOverThreshold[n_ahdcrows]/F",  512000);

	// Particle banks
	int   nrectracks;
	int   rec_track_pid[maxrectracks];
	float rec_track_x[maxrectracks];
	float rec_track_y[maxrectracks];
	float rec_track_z[maxrectracks];
	float rec_track_vt[maxrectracks];
	float rec_track_px[maxrectracks];
	float rec_track_py[maxrectracks];
	float rec_track_pz[maxrectracks];
	float rec_track_beta[maxrectracks];
	float rec_track_chi2pid[maxrectracks];
	int   rec_track_status[maxrectracks];
	int   rec_track_charge[maxrectracks];

	treeOutput->Branch("nrectracks",        &nrectracks,        "nrectracks/I",                    512000);
	treeOutput->Branch("rec_track_pid",     &rec_track_pid,     "rec_track_pid[nrectracks]/I",     512000);
	treeOutput->Branch("rec_track_x",       &rec_track_x,       "rec_track_x[nrectracks]/F",       512000);
	treeOutput->Branch("rec_track_y",       &rec_track_y,       "rec_track_y[nrectracks]/F",       512000);
	treeOutput->Branch("rec_track_z",       &rec_track_z,       "rec_track_z[nrectracks]/F",       512000);
	treeOutput->Branch("rec_track_vt",      &rec_track_vt,      "rec_track_vt[nrectracks]/F",      512000);
	treeOutput->Branch("rec_track_px",      &rec_track_px,      "rec_track_px[nrectracks]/F",      512000);
	treeOutput->Branch("rec_track_py",      &rec_track_py,      "rec_track_py[nrectracks]/F",      512000);
	treeOutput->Branch("rec_track_pz",      &rec_track_pz,      "rec_track_pz[nrectracks]/F",      512000);
	treeOutput->Branch("rec_track_beta",    &rec_track_beta,    "rec_track_beta[nrectracks]/F",    512000);
	treeOutput->Branch("rec_track_chi2pid", &rec_track_chi2pid, "rec_track_chi2pid[nrectracks]/F", 512000);
	treeOutput->Branch("rec_track_status",  &rec_track_status,  "rec_track_status[nrectracks]/I",  512000);
	treeOutput->Branch("rec_track_charge",  &rec_track_charge,  "rec_track_charge[nrectracks]/I",  512000);

	int nEv = 0;
	long runtime = 0;

	hipo::reader reader;
	reader.open(hipoFile);
	hipo::dictionary factory;
	reader.readDictionary(factory);
	factory.show();

	hipo::writer writer;
	writer.getDictionary().addSchema(factory.getSchema("REC::Particle")); // EFMC: Could this be an issue ????
	writer.open(outputFileHipo);

	hipo::bank hits(factory.getSchema("AHDC::hits"));
	hipo::bank preclusters(factory.getSchema("AHDC::preclusters"));
	hipo::bank clusters(factory.getSchema("AHDC::clusters"));
	// hipo::bank mctracks(factory.getSchema("MC::Particle"));
	hipo::bank tracks(factory.getSchema("AHDC::track"));
	hipo::bank kftracks(factory.getSchema("AHDC::kftrack"));
	hipo::bank ahdc_adc(factory.getSchema("AHDC::adc"));
	hipo::bank rec_tracks(factory.getSchema("REC::Particle"));

	hipo::event event;
	int counter = 0;

	hipo::benchmark writerBenchmark;
	hipo::benchmark readerBenchmark;
	hipo::benchmark transferBenchmark;
	hipo::benchmark restBenchmark;
	hipo::benchmark writerHipoBenchmark;



	while (reader.next() == true) {
		readerBenchmark.resume();
		reader.read(event);
		event.getStructure(hits);
		event.getStructure(preclusters);
		event.getStructure(clusters);
		// if (do_mc)
		// 	event.getStructure(mctracks);
		event.getStructure(tracks);
		event.getStructure(kftracks);
		event.getStructure(ahdc_adc);
		event.getStructure(rec_tracks);

		readerBenchmark.pause();

		nhits = hits.getRows();
		assert(nhits < maxhits);
		for (int i = 0; i < nhits; i++) {
			hitid[i] = hits.getInt("id", i);
			hitlayer[i] = hits.getInt("layer", i);
			hitsuperlayer[i] = hits.getInt("superlayer", i);
			hitwire[i] = hits.getInt("wire", i);
			hitdoca[i] = hits.getFloat("doca", i);
		}

		npreclusters = preclusters.getRows();
		assert(npreclusters < maxpreclusters);
		for (int i = 0; i < nhits; i++) {
			preclusterX[i] = preclusters.getFloat("x", i);
			preclusterY[i] = preclusters.getFloat("y", i);
		}
		nclusters = clusters.getRows();
		assert(nclusters < maxclusters);
		for (int i = 0; i < nhits; i++) {
			clusterX[i] = clusters.getFloat("x", i);
			clusterY[i] = clusters.getFloat("y", i);
			clusterZ[i] = clusters.getFloat("z", i);
		}

		// if (do_mc) {
		// 	nmctracks = mctracks.getRows();
		// 	assert(nmctracks < maxmctracks);
		// 	for (int i = 0; i < nhits; i++) {
		// 		mctrackpid[i] = mctracks.getInt("pid", i);
		// 		mctrackpx[i] = mctracks.getFloat("px", i);
		// 		mctrackpy[i] = mctracks.getFloat("py", i);
		// 		mctrackpz[i] = mctracks.getFloat("pz", i);
		// 		mctrackvx[i] = mctracks.getFloat("vx", i);
		// 		mctrackvy[i] = mctracks.getFloat("vy", i);
		// 		mctrackvz[i] = mctracks.getFloat("vz", i);
		// 		mctrackvt[i] = mctracks.getFloat("vt", i);
		// 	}
		// }

		ntracks = tracks.getRows();
		assert(ntracks < maxtracks);
		for (int i = 0; i < nhits; i++) {
			trackx[i] = tracks.getFloat("x", i);
			tracky[i] = tracks.getFloat("y", i);
			trackz[i] = tracks.getFloat("z", i);
			trackpx[i] = tracks.getFloat("px", i);
			trackpy[i] = tracks.getFloat("py", i);
			trackpz[i] = tracks.getFloat("pz", i);
		}

		nrectracks = rec_tracks.getRows();
		assert(nrectracks < maxrectracks);
		for (int i = 0; i < nrectracks; i++) {
			rec_track_pid[i]     = rec_tracks.getInt("pid", i);
			rec_track_x[i]       = rec_tracks.getFloat("vx", i);
			rec_track_y[i]       = rec_tracks.getFloat("vy", i);
			rec_track_z[i]       = rec_tracks.getFloat("vz", i);
			rec_track_vt[i]      = rec_tracks.getFloat("vt", i);
			rec_track_px[i]      = rec_tracks.getFloat("px", i);
			rec_track_py[i]      = rec_tracks.getFloat("py", i);
			rec_track_pz[i]      = rec_tracks.getFloat("pz", i);
			rec_track_beta[i]    = rec_tracks.getFloat("beta", i);
			rec_track_chi2pid[i] = rec_tracks.getFloat("chi2pid", i);
			rec_track_status[i]  = rec_tracks.getInt("status", i);
			rec_track_charge[i]  = rec_tracks.getInt("charge", i);
		}

		nkftracks = kftracks.getRows();
		assert(nkftracks < maxkftracks);
		for (int i = 0; i < nhits; i++) {
			kftrack_id[i] = kftracks.getInt("trackid", i);
			kftrackx[i]   = kftracks.getFloat("x", i);
			kftracky[i]   = kftracks.getFloat("y", i);
			kftrackz[i]   = kftracks.getFloat("z", i);
			kftrackpx[i]  = kftracks.getFloat("px", i);
			kftrackpy[i]  = kftracks.getFloat("py", i);
			kftrackpz[i]  = kftracks.getFloat("pz", i);
			kftrackpath[i]= kftracks.getFloat("path", i);
			kftrackdedx[i]= kftracks.getFloat("dEdx", i);
		}

		n_ahdcrows = ahdc_adc.getRows();
		assert(n_ahdcrows < maxahdcadc);
		for (int i = 0; i < n_ahdcrows; i++) {
			AHDC_sector[i]            = ahdc_adc.getInt("sector", i);
			AHDC_component[i]         = ahdc_adc.getInt("component", i);
			AHDC_order[i]             = ahdc_adc.getInt("order", i);
			AHDC_maxADC[i]            = ahdc_adc.getInt("integral", i);
			AHDC_leadingEdgeTime[i]   = ahdc_adc.getFloat("leadingEdgeTime", i);
			AHDC_timeOverThreshold[i] = ahdc_adc.getFloat("timeOverThreshold", i);
		}
                
                if (nhits > 0) {
			writerBenchmark.resume();
			treeOutput->Fill();
			writerBenchmark.pause();

			writerHipoBenchmark.resume();
			event.reset();
			event.addStructure(hits);
			event.addStructure(preclusters);
			event.addStructure(clusters);
			writer.addEvent(event);
			writerHipoBenchmark.pause();
		}
	}
	writer.close();
	fileOutput->Write();
	fileOutput->Close();

	printf("processed events = %d, root      (WRITE) : time = %10.2f sec , count = %d\n",
	       counter, writerBenchmark.getTimeSec(), writerBenchmark.getCounter());
	printf("processed events = %d, hipo      (WRITE) : time = %10.2f sec , count = %d\n",
	       counter, writerHipoBenchmark.getTimeSec(), writerHipoBenchmark.getCounter());
	printf("processed events = %d, benchmark (READ)  : time = %10.2f sec , count = %d\n",
	       counter, readerBenchmark.getTimeSec(), readerBenchmark.getCounter());
	printf("processed events = %d, benchmark (COPY)  : time = %10.2f sec , count = %d\n",
	       counter, transferBenchmark.getTimeSec(), transferBenchmark.getCounter());
	printf("processed events = %d, benchmark (REST)  : time = %10.2f sec , count = %d\n",
	       counter, restBenchmark.getTimeSec(), restBenchmark.getCounter());

	double total_time = writerBenchmark.getTimeSec() +
			   readerBenchmark.getTimeSec() +
			   transferBenchmark.getTimeSec() +
			   restBenchmark.getTimeSec();
	printf("\n total time = %10.2f\n", total_time);

	printf("[FIN] (C++ ) hipo write : %10.4f sec\n", writerHipoBenchmark.getTimeSec());
	printf("[FIN] (C++ ) root write : %10.4f sec\n", writerBenchmark.getTimeSec());
}

//### END OF GENERATED CODE
