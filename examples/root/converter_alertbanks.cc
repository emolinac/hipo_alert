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
#define maxtracks       1000
#define maxkftracks     1000
#define maxahdcadc      1000
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
	int   n_ahdcrecohits;
	int   ahdc_hitid[maxhits];
	int   ahdc_trackid[maxhits];
	int   ahdc_layer[maxhits];
	int   ahdc_wire[maxhits];
	int   ahdc_component[maxhits];
	float ahdc_time[maxhits];
	float ahdc_energy[maxhits];

	treeOutput->Branch("n_ahdcrecohits", &n_ahdcrecohits, "n_ahdcrecohits/I",                 512000);
	treeOutput->Branch("ahdc_hitid",     &ahdc_hitid,     "ahdc_hitid[n_ahdcrecohits]/I",     512000);
	treeOutput->Branch("ahdc_trackid",   &ahdc_trackid,   "ahdc_trackid[n_ahdcrecohits]/I",   512000);
	treeOutput->Branch("ahdc_layer",     &ahdc_layer,     "ahdc_layer[n_ahdcrecohits]/I",     512000);
	treeOutput->Branch("ahdc_wire",      &ahdc_wire,    "ahdc_wire[n_ahdcrecohits]/I",    512000);
	treeOutput->Branch("ahdc_time",      &ahdc_time,      "ahdc_time[n_ahdcrecohits]/F",      512000);
	
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

	/* AHDC tracks */
	int ntracks;
	float trackx[maxtracks];
	float tracky[maxtracks];
	float trackz[maxtracks];
	float trackpx[maxtracks];
	float trackpy[maxtracks];
	float trackpz[maxtracks];
	int   trackid[maxtracks];

	treeOutput->Branch("ahdc_ntracks", &ntracks, "ntracks/I",          512000); 
	treeOutput->Branch("ahdc_track_x",  &trackx,  "trackx[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_track_y",  &tracky,  "tracky[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_track_z",  &trackz,  "trackz[ntracks]/F",  512000);
	treeOutput->Branch("ahdc_track_px", &trackpx, "trackpx[ntracks]/F", 512000);
	treeOutput->Branch("ahdc_track_py", &trackpy, "trackpy[ntracks]/F", 512000);
	treeOutput->Branch("ahdc_track_pz", &trackpz, "trackpz[ntracks]/F", 512000);
	treeOutput->Branch("ahdc_track_id", &trackid, "trackid[ntracks]/I", 512000);

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
	int kftracksumadc[maxkftracks];

	treeOutput->Branch("ahdc_nkftracks",      &nkftracks,     "nkftracks/I", 	        512000);
	treeOutput->Branch("ahdc_kftrack_id",     &kftrack_id,    "kftrack_id[nkftracks]/I",    512000);
	treeOutput->Branch("ahdc_kftrack_x",      &kftrackx,      "kftrackx[nkftracks]/F",      512000);
	treeOutput->Branch("ahdc_kftrack_y",      &kftracky,      "kftracky[nkftracks]/F",      512000);
	treeOutput->Branch("ahdc_kftrack_z",      &kftrackz,      "kftrackz[nkftracks]/F",      512000);
	treeOutput->Branch("ahdc_kftrack_px",     &kftrackpx,     "kftrackpx[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftrack_py",     &kftrackpy,     "kftrackpy[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftrack_pz",     &kftrackpz,     "kftrackpz[nkftracks]/F",     512000);
	treeOutput->Branch("ahdc_kftrack_path",   &kftrackpath,   "kftrackpath[nkftracks]/F",   512000);
	treeOutput->Branch("ahdc_kftrack_dedx",   &kftrackdedx,   "kftrackdedx[nkftracks]/F",   512000);
	treeOutput->Branch("ahdc_kftrack_sumadc", &kftracksumadc, "kftracksumadc[nkftracks]/I", 512000);

	/* AHDC ADC information */
	int   n_ahdcadcrows;
	int   ahdc_adc_sector[maxahdcadc];
	int   ahdc_adc_layer[maxahdcadc];
	int   ahdc_adc_component[maxahdcadc];
	int   ahdc_adc_order[maxahdcadc];
	int   ahdc_adc_sumadc[maxahdcadc];
	float ahdc_adc_leadingEdgeTime[maxahdcadc];
	float ahdc_adc_timeOverThreshold[maxahdcadc];
	int   ahdc_adc_wfType[maxahdcadc];

	treeOutput->Branch("ahdc_adc_nadcrows",          &n_ahdcadcrows,              "n_ahdcadcrows/I",                 	      512000);
	treeOutput->Branch("ahdc_adc_sector",            &ahdc_adc_sector,            "ahdc_adc_sector[n_ahdcadcrows]/I",             512000);
	treeOutput->Branch("ahdc_adc_layer",             &ahdc_adc_layer,             "ahdc_adc_layer[n_ahdcadcrows]/I",             512000);
	treeOutput->Branch("ahdc_adc_component",         &ahdc_adc_component,         "ahdc_adc_component[n_ahdcadcrows]/I",          512000);
	treeOutput->Branch("ahdc_adc_order",             &ahdc_adc_order,             "ahdc_adc_order[n_ahdcadcrows]/I",              512000);
	treeOutput->Branch("ahdc_adc_sumadc",            &ahdc_adc_sumadc,            "ahdc_adc_sumadc[n_ahdcadcrows]/I",             512000);
	treeOutput->Branch("ahdc_adc_leadingEdgeTime",   &ahdc_adc_leadingEdgeTime,   "ahdc_adc_leadingEdgeTime[n_ahdcadcrows]/F",    512000);
	treeOutput->Branch("ahdc_adc_timeOverThreshold", &ahdc_adc_timeOverThreshold, "ahdc_adc_timeOverThreshold[n_ahdcadcrows]/F",  512000);
	treeOutput->Branch("ahdc_adc_wfType",            &ahdc_adc_wfType,            "ahdc_adc_wfType[n_ahdcadcrows]/I",             512000);

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

		n_ahdcrecohits = hits.getRows();
		assert(n_ahdcrecohits < maxhits);
		for (int i = 0; i < n_ahdcrecohits; i++) {
			ahdc_hitid[i]     = hits.getInt("id", i);
			ahdc_trackid[i]   = hits.getInt("trackid", i);
			ahdc_layer[i]     = hits.getInt("layer", i);
			ahdc_wire[i]      = hits.getInt("wire", i);
			ahdc_time[i]      = hits.getFloat("time", i);
		}

		npreclusters = preclusters.getRows();
		assert(npreclusters < maxpreclusters);
		for (int i = 0; i < npreclusters; i++) {
			preclusterX[i] = preclusters.getFloat("x", i);
			preclusterY[i] = preclusters.getFloat("y", i);
		}

		nclusters = clusters.getRows();
		assert(nclusters < maxclusters);
		for (int i = 0; i < nclusters; i++) {
			clusterX[i] = clusters.getFloat("x", i);
			clusterY[i] = clusters.getFloat("y", i);
			clusterZ[i] = clusters.getFloat("z", i);
		}

		ntracks = tracks.getRows();
		assert(ntracks < maxtracks);
		for (int i = 0; i < ntracks; i++) {
			trackx[i] = tracks.getFloat("x", i);
			tracky[i] = tracks.getFloat("y", i);
			trackz[i] = tracks.getFloat("z", i);
			trackpx[i] = tracks.getFloat("px", i);
			trackpy[i] = tracks.getFloat("py", i);
			trackpz[i] = tracks.getFloat("pz", i);
			trackid[i] = tracks.getInt("trackid", i);
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
		for (int i = 0; i < nkftracks; i++) {
			kftrack_id[i]    = kftracks.getInt("trackid", i);
			kftrackx[i]      = kftracks.getFloat("x", i);
			kftracky[i]      = kftracks.getFloat("y", i);
			kftrackz[i]      = kftracks.getFloat("z", i);
			kftrackpx[i]     = kftracks.getFloat("px", i);
			kftrackpy[i]     = kftracks.getFloat("py", i);
			kftrackpz[i]     = kftracks.getFloat("pz", i);
			kftrackpath[i]   = kftracks.getFloat("path", i);
			kftrackdedx[i]   = kftracks.getFloat("dEdx", i);
			kftracksumadc[i] = kftracks.getInt("sum_adc", i);
		}

		n_ahdcadcrows = ahdc_adc.getRows();
		assert(n_ahdcadcrows < maxahdcadc);
		for (int i = 0; i < n_ahdcadcrows; i++) {
			ahdc_adc_sector[i]            = ahdc_adc.getInt("sector", i);
			ahdc_adc_layer[i]             = ahdc_adc.getInt("layer", i);
			ahdc_adc_component[i]         = ahdc_adc.getInt("component", i);
			ahdc_adc_wfType[i]            = ahdc_adc.getInt("wfType", i);
			ahdc_adc_order[i]             = ahdc_adc.getInt("order", i);
			ahdc_adc_sumadc[i]            = ahdc_adc.getInt("integral", i);
			ahdc_adc_leadingEdgeTime[i]   = ahdc_adc.getFloat("leadingEdgeTime", i);
			ahdc_adc_timeOverThreshold[i] = ahdc_adc.getFloat("timeOverThreshold", i);
		}
                
                if (n_ahdcrecohits > 0) {
			writerBenchmark.resume();
			treeOutput->Fill();
			writerBenchmark.pause();

			writerHipoBenchmark.resume();
			event.reset();
			event.addStructure(hits);
			event.addStructure(preclusters);
			event.addStructure(clusters);
			event.addStructure(tracks);
			event.addStructure(kftracks);
			event.addStructure(ahdc_adc);
			event.addStructure(rec_tracks);
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
