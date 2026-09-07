// ./Base/analysisCode/Analysis.cpp

namespace Config {
	static const double ADC_to_mV = 0.1220703125;
	static const double Impedence = 50.0; // Ohm
	static const double TimePerSample = 2.0; // ns
}

double GetPedestal(unsigned int* wave_data, int ped_start, int ped_end, int nSamples) {
	if (!wave_data || ped_start < 0 || ped_end > nSamples || ped_start >= ped_end) return 0.0;

	double sum = 0.0;
	for (int i = ped_start; i < ped_end; ++i) {
		sum += static_cast<double>(wave_data[i]);
	}
	return sum / (ped_end - ped_start);
}

double GetCharge(unsigned int* wave_data, double pedestal, int sig_start, int sig_end, int nSamples) {
	if (!wave_data || sig_start >= sig_end || sig_end > nSamples) {
		std::cout << "Check for range! " << std::endl;
		return 0.0;
	}

	double integral = 0.0;
	for (int i = sig_start; i < sig_end; ++i) {
		integral += (pedestal - static_cast<double>(wave_data[i]));
	}

	// Charge [pC] = TimePerSample[ns] * (ADC_to_mV * integral) / Impedence[Ohm]
	double charge_pC = Config::TimePerSample * (Config::ADC_to_mV * integral) / Config::Impedence;
	return charge_pC;
}

std::pair<double, int> GetPeakInfo(unsigned int* wave_data, double pedestal, int sig_start, int sig_end, int nSamples) {
	if (!wave_data || sig_start >= sig_end || sig_end > nSamples) {
		std::cout << "Check for range! " << std::endl;
		return {-1e9, -1}; 
	}

	double max_amp = -1e9;
	int max_time = -1; 

	for (int i = sig_start; i < sig_end; ++i) {
		double signal = pedestal - static_cast<double>(wave_data[i]);
		if (signal > max_amp) {
			max_amp = signal;
			max_time = i; 		
		}
	}
	return std::make_pair(max_amp, max_time); 
}

int GetTimeAtThreshold(unsigned int* wave_data, double threshold, int nSamples, int searchStart = 0, int searchEnd = -1) {
    if (!wave_data) return -1;

    if (searchEnd == -1) searchEnd = nSamples;
    if (searchStart < 0) searchStart = 0;
    if (searchEnd > nSamples) searchEnd = nSamples;
    if (searchStart >= searchEnd) return -1;

    int validDurationRequired = 3; 
    int consecutiveBelow = 0;
    int firstCrossingTime = -1;

    for (int iSample = searchStart; iSample < searchEnd; ++iSample) {
        if (static_cast<double>(wave_data[iSample]) < threshold) {
            if (consecutiveBelow == 0) {
                firstCrossingTime = iSample;
            }
            consecutiveBelow++;
            
            if (consecutiveBelow >= validDurationRequired) {
                return firstCrossingTime;
            }
        } else {
            consecutiveBelow = 0;
            firstCrossingTime = -1;
        }
    }
    return -1;
}

std::vector<int> GetTimesBelowThreshold(unsigned int* wave_data, double threshold, int nSamples, int searchStart = 0, int searchEnd = -1) {
    std::vector<int> hitTimes;
    
    if (!wave_data) return hitTimes;

    if (searchEnd == -1) searchEnd = nSamples;
    if (searchStart < 0) searchStart = 0;
    if (searchEnd > nSamples) searchEnd = nSamples;
    if (searchStart >= searchEnd) return hitTimes;

    bool inPulse = false;
    double local_lowest_adc = 1e9;
    double local_highest_adc = -1e9;
    int last_peak_time = -1000;
    
    int min_distance_samples = 4;
    double prominence_adc = 1.0 / Config::ADC_to_mV;

    int consecutiveBelow = 0;
    int validDurationRequired = 3;
    int potentialStartTime = -1;
    bool pulseValidated = false;

    for (int iSample = searchStart; iSample < searchEnd; ++iSample) {
        double current_adc = static_cast<double>(wave_data[iSample]);

        if (current_adc < threshold) {
            if (!inPulse) {
                if (consecutiveBelow == 0) {
                    potentialStartTime = iSample;
                }
                consecutiveBelow++;

                if (consecutiveBelow >= validDurationRequired && !pulseValidated) {
                    inPulse = true;
                    pulseValidated = true;
                    hitTimes.push_back(potentialStartTime);
                    local_lowest_adc = current_adc;
                    local_highest_adc = current_adc;
                    last_peak_time = iSample;
                }
            } else {
                if (current_adc < local_lowest_adc) local_lowest_adc = current_adc;
                if (current_adc > local_highest_adc) local_highest_adc = current_adc;

                if ((local_highest_adc - current_adc) > prominence_adc && (iSample - last_peak_time) >= min_distance_samples) {
                    hitTimes.push_back(iSample);
                    local_lowest_adc = current_adc;
                    local_highest_adc = current_adc;
                    last_peak_time = iSample;
                }
            }
        } else {
            inPulse = false;
            consecutiveBelow = 0;
            pulseValidated = false;
        }
    }
    return hitTimes;
}
