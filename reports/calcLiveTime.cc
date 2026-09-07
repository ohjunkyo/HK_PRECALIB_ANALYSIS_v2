void calcLiveTime(TString rootfile)
{

	bool verbose = false;
	double NsPerTriggerTimeTag = 8.; // [ns]

	TChain * ch = new TChain("T");
	ch -> Add(rootfile);
	ch -> LoadTree(0);

	unsigned int EventNumber;
	ch -> SetBranchAddress("EventNumber", &EventNumber);

	unsigned int TriggerTimeTag;
	unsigned int TriggerTimeTag_prev = 0.;
	ch -> SetBranchAddress("TriggerTimeTag", &TriggerTimeTag);

	int NEntry = ch -> GetEntries();

	double LiveTime = 0.;
	for(int iEntry = 0; iEntry < NEntry; iEntry++)
	{

		if ((iEntry % 1000) == 0) {
			double percent = (iEntry * 100.0) / NEntry;
			std::cout << "Processing " << "... " << percent << "%\r" << std::flush;
		}

		ch -> GetEntry(iEntry);


		double TriggerTimeTagDiff = (double)TriggerTimeTag - (double)TriggerTimeTag_prev;
		if(TriggerTimeTagDiff < 0.) TriggerTimeTagDiff += (double)pow(2, 31);
		if(iEntry != 0) LiveTime += TriggerTimeTagDiff * NsPerTriggerTimeTag * 1e-9;
		if(verbose){
			if(iEntry % 1000 == 0)
				std::cout << iEntry << '\t' << EventNumber << '\t' << TriggerTimeTag << '\t' << TriggerTimeTagDiff << '\t' << LiveTime << std::endl;
		}

		TriggerTimeTag_prev = TriggerTimeTag;
	}

	double Rate = 0;
	if (LiveTime > 0) Rate = NEntry / LiveTime;

	std::cout << "=====[RESULT]=====" << std::endl;
	std::cout << " NEntry : " << NEntry << std::endl;
	std::cout << " LiveTime [s] : " << LiveTime << std::endl;
	std::cout << " Rate [Hz] : " << Rate << std::endl;
	std::cout << "==================" << std::endl;

	delete ch;

}
