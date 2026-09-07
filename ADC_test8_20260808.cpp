/* =============================================================================
 * V7 Optimized DAQ (ROOT & CSV Dual-Format Support Version)
 * Features:
 * - Perfectly matched the exact console output style of ADC_test6.
 * - Fixed the configuration failure bug caused by external trigger (-e) flag.
 * - Single continuous 16-bit memory block used for both ROOT and CSV parsing.
 * =============================================================================
 */

extern "C" {
#include "CAENDigitizer.h"
#include "CAENComm.h"
}
#include <iostream>
#include <fstream>
#include <bitset>
#include <ctime>
#include <boost/program_options.hpp>
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <chrono>   
#include <iomanip>  
#include <cstring>
#include <unistd.h> 

#include "/home/precalkor/Integrated_Control_SW/DAQ_Control_SW/config3.h"

using namespace std;
namespace po = boost::program_options;

uint32_t activechannels[8];
uint32_t trigger_value2[8] = {15000, 15030, 15000, 15000, 15000, 15000, 15000, 15000}; 
uint32_t offset_value2[8]  = {7050, 7050, 7050, 7050, 7050, 7050, 7050, 7050};

int getchannels(uint32_t* activechannels_pointer, uint32_t channels) {
    int nchannels = 0;
    for (int i = 0; i < 8; i++) {
        if ((channels >> i) & 1) {
            activechannels_pointer[nchannels++] = (uint32_t)i;
        }
    }
    return nchannels;
}

int main(int ac, char *av[]) {
    uint32_t recordlength, post_trigger, nevents_user, trigger_value, offset_value;
    string outputfile, channel_input, run_mode_input, shifter_input, expert_input, note_input;
    bool external_trigger = false, root = false;
    bool add_time_to_filename = false, trigger_risingedge = true;
    int trgOnOff = 0, offsetOnOff = 0;

    // Metadata Variables for RunInfo Tree
    char b_SN1[50], b_SN2[50], b_SN3[50], b_NOTE[300], b_RunMode[20], b_Shifter[20], b_Expert[20];
    // Cable direction (A~H) per channel — needed downstream to convert the
    // rotation-stage angle into the PMT incidence (Hamamatsu) angle.
    char b_Direction1[8], b_Direction2[8], b_Direction3[8];
    // [DEVELOP] double b_ElapsedTime = 0.0;   // wall-clock acquisition time [s], for rate-scan anchor validation
    int b_HV1, b_HV2, b_HV3, b_Laser, b_Wavelength, b_RawRotateAngle2, b_RawTiltAngle2, b_RawRotateAngle3, b_RawTiltAngle3;

    // Initialize Metadata from config3.h
    strncpy(b_SN1, SN1.c_str(), 49); strncpy(b_SN2, SN2.c_str(), 49); strncpy(b_SN3, SN3.c_str(), 49);
    strncpy(b_Direction1, direction1.c_str(), 7); b_Direction1[7] = '\0';
    strncpy(b_Direction2, direction2.c_str(), 7); b_Direction2[7] = '\0';
    strncpy(b_Direction3, direction3.c_str(), 7); b_Direction3[7] = '\0';

    int def_rot2 = 0, def_tilt2 = 0, def_rot3 = 0, def_tilt3 = 0;
    try {
        b_HV1 = stoi(HV1); b_HV2 = stoi(HV2); b_HV3 = stoi(HV3); b_Laser = stoi(Laser); b_Wavelength = stoi(Wavelength);
        def_rot2 = stoi(RotateAngle2); def_tilt2 = stoi(TiltAngle2);
        def_rot3 = stoi(RotateAngle3); def_tilt3 = stoi(TiltAngle3);
    } catch (...) { /* Handled warnings if needed */ }

    // Boost Command Line Options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("just-check-connection,j", "Only check if the digitizer connection works, then exit (no data taking). Used by the GUI's DAQ Status indicator.")
        ("external-trigger,e", "Use external trigger")
        ("run-mode,M", po::value<string>(&run_mode_input)->default_value("Laser"), "Run mode: Laser or Dark")
        ("falling-edge,f", "Trigger on falling instead of rising edge")
        ("timestamp,m", "Add timestamp to output file name")
        ("root,r", "Create ROOT file output")
        ("output-file,o", po::value<string>(&outputfile)->default_value("buffer"), "Output filename (without extension)")
        ("channels,c", po::value<string>(&channel_input)->default_value("00001111"), "Channel Mask (e.g., 00001111)")
        ("triggervalue,t", po::value<uint32_t>(&trigger_value)->default_value(0xFFF), "Default trigger value")
        ("offset,s", po::value<uint32_t>(&offset_value)->default_value(0x7FFF), "Default DC offset")
        ("recordlength,d", po::value<uint32_t>(&recordlength)->default_value(1024), "Record length in samples")
        ("posttrigger,p", po::value<uint32_t>(&post_trigger)->default_value(60), "Post trigger percentage")
        ("neventsuser,n", po::value<uint32_t>(&nevents_user)->default_value(0), "Number of events to take")
        ("EachCH-offset-OnOff,y", po::value<int>(&offsetOnOff)->default_value(0), "Use hardcoded offset values")
        ("EachCH-trggerval-OnOff,z", po::value<int>(&trgOnOff)->default_value(0), "Use hardcoded trigger thresholds")
        ("rot2", po::value<int>(&b_RawRotateAngle2)->default_value(def_rot2), "Overwrite Rotation Angle 2")
        ("tilt2", po::value<int>(&b_RawTiltAngle2)->default_value(def_tilt2), "Overwrite Tilt Angle 2")
        ("rot3", po::value<int>(&b_RawRotateAngle3)->default_value(def_rot3), "Overwrite Rotation Angle 3")
        ("tilt3", po::value<int>(&b_RawTiltAngle3)->default_value(def_tilt3), "Overwrite Tilt Angle 3")
        ("laser", po::value<int>(&b_Laser)->default_value(b_Laser), "Overwrite Laser mA")
        ("wavelength", po::value<int>(&b_Wavelength)->default_value(b_Wavelength), "Overwrite Wavelength")
        ("hv1", po::value<int>(&b_HV1)->default_value(b_HV1), "Overwrite HV1")
        ("hv2", po::value<int>(&b_HV2)->default_value(b_HV2), "Overwrite HV2")
        ("hv3", po::value<int>(&b_HV3)->default_value(b_HV3), "Overwrite HV3")
        ("shifter", po::value<string>(&shifter_input)->default_value(b_Shifter), "Overwrite Shifter")
        ("expert", po::value<string>(&expert_input)->default_value(b_Expert), "Overwrite Expert")
        ("note", po::value<string>(&note_input)->default_value(b_NOTE), "Overwrite NOTE");
    
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) { cout << desc << endl; return 0; }
    bool just_check = vm.count("just-check-connection");
    external_trigger = vm.count("external-trigger");
    trigger_risingedge = !vm.count("falling-edge");
    add_time_to_filename = vm.count("timestamp");
    root = vm.count("root");

    strncpy(b_RunMode, run_mode_input.c_str(), 19);
    strncpy(b_Shifter, shifter_input.c_str(), 19); b_Shifter[19] = '\0';
    strncpy(b_Expert, expert_input.c_str(), 19);   b_Expert[19] = '\0';
    strncpy(b_NOTE, note_input.c_str(), 299);      b_NOTE[299] = '\0';

    // Calculate channel mask and identify active channels
    uint32_t channels_mask = 0;
    for (int u = 0; u < 8; u++) {
        channels_mask = (channels_mask << 1) + (channel_input.at(u) - '0');
    }
    uint16_t nchannels = getchannels(&activechannels[0], uint32_t(bitset<8>(channel_input).to_ulong()));

    // Open Digitizer
    int handle;
    cout << "Connect to device..." << endl;
    CAEN_DGTZ_ErrorCode open_ret = CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_USB, 0, 0, 0, &handle);

    // "Already open" (-25) means the board is present but busy (e.g. a run is in
    // progress), which still counts as "connected" for a status check.
    bool connected = (open_ret == CAEN_DGTZ_Success ||
                      open_ret == CAEN_DGTZ_DigitizerAlreadyOpen);

    if (just_check) {
        // Connection-check mode (-j): report status and exit WITHOUT taking data.
        // The GUI's DAQ Status indicator relies on the exit code (0 = connected).
        if (connected) {
            cout << "[INFO] DAQ connection OK." << endl;
            // Only release the board if WE actually opened it (not when -25 busy).
            if (open_ret == CAEN_DGTZ_Success) CAEN_DGTZ_CloseDigitizer(handle);
            return 0;
        } else {
            cout << "Communication error (open code " << open_ret << ")." << endl;
            return 1;
        }
    }

    if (open_ret != CAEN_DGTZ_Success) {
        cout << "[ERROR] Device could not be opened." << endl;
        return 1;
    }

    // Get board info
    CAEN_DGTZ_BoardInfo_t BoardInfo;
    CAEN_DGTZ_GetInfo(handle, &BoardInfo);

    // Format Filename accurately 
    string filename = outputfile;
    string target_ext = root ? ".root" : ".csv";
    
    if (filename.size() > 5 && filename.substr(filename.size() - 5) == ".root") filename.erase(filename.size() - 5);
    if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".csv")  filename.erase(filename.size() - 4);

    if (add_time_to_filename) {
        time_t t_fn = time(0);
        tm* tm_fn = localtime(&t_fn);
        char time_str[50];
        sprintf(time_str, "%d%02d%02d%02d%02d%02d", 1900 + tm_fn->tm_year, tm_fn->tm_mon + 1, tm_fn->tm_mday, tm_fn->tm_hour, tm_fn->tm_min, tm_fn->tm_sec);
        filename += time_str;
    }
    filename += target_ext;

    TFile *MyFile = nullptr;
    TTree *tree = nullptr;
    TTree *infoTree = nullptr;
    ofstream stream;

    uint16_t* adc_buffer = new uint16_t[8 * recordlength];
    memset(adc_buffer, 0, 8 * recordlength * sizeof(uint16_t)); 

    uint32_t TriggerTimeTag = 0;
    int nevents = 0;

    if (root) {
        MyFile = new TFile(filename.c_str(), "RECREATE");

        infoTree = new TTree("RunInfo", "Metadata for this run");
        infoTree->Branch("SN1", b_SN1, "SN1/C"); infoTree->Branch("HV1", &b_HV1, "HV1/I");
        infoTree->Branch("Direction1", b_Direction1, "Direction1/C");
        infoTree->Branch("SN2", b_SN2, "SN2/C"); infoTree->Branch("HV2", &b_HV2, "HV2/I");
        infoTree->Branch("Direction2", b_Direction2, "Direction2/C");
        infoTree->Branch("RawRotateAngle2", &b_RawRotateAngle2, "RawRotateAngle2/I");
        infoTree->Branch("RawTiltAngle2", &b_RawTiltAngle2, "RawTiltAngle2/I");
        infoTree->Branch("SN3", b_SN3, "SN3/C"); infoTree->Branch("HV3", &b_HV3, "HV3/I");
        infoTree->Branch("Direction3", b_Direction3, "Direction3/C");
        infoTree->Branch("RawRotateAngle3", &b_RawRotateAngle3, "RawRotateAngle3/I");
        infoTree->Branch("RawTiltAngle3", &b_RawTiltAngle3, "RawTiltAngle3/I");
        infoTree->Branch("Wavelength", &b_Wavelength, "Wavelength/I");
        infoTree->Branch("Laser_mA", &b_Laser, "Laser_mA/I");
        infoTree->Branch("NOTE", b_NOTE, "NOTE/C");
        infoTree->Branch("RunMode", b_RunMode, "RunMode/C");
        infoTree->Branch("Expert", b_Expert, "Expert/C");
        infoTree->Branch("Shifter", b_Shifter, "Shifter/C");
        // [DEVELOP] infoTree->Branch("ElapsedTime", &b_ElapsedTime, "ElapsedTime/D");
        //
        // RecordLength, PostTrigger, OffsetValue<N> and TriggerValue<N> all used
        // to be branches on tree "T" (Digitizer Raw Data), rewritten unchanged
        // into EVERY event -- but each is set once from the CAEN_DGTZ_Set*()
        // calls before the acquisition loop and never changes during a run. At
        // 300k events/run x up to 8 channels x 92+ runs that is a real chunk of
        // Data/RAW for zero information gain, so they now live in RunInfo,
        // written exactly once. Only TriggerTimeTag genuinely varies per event.
        //
        // ChannelMask is new here: prod_ntp_v7.C used to infer the active channel
        // set from which OffsetValue<N> branches existed on "T", which is exactly
        // what kept those branches pinned there. Recording the real channel_input
        // string gives analysis an explicit source instead. prod_ntp_v7.C reads
        // RunInfo first and falls back to the old per-event branches, so files
        // written before this change still analyse identically.
        char b_ChannelMask[16];
        strncpy(b_ChannelMask, channel_input.c_str(), sizeof(b_ChannelMask) - 1);
        b_ChannelMask[sizeof(b_ChannelMask) - 1] = '\0';
        infoTree->Branch("ChannelMask", b_ChannelMask, "ChannelMask/C");
        infoTree->Branch("RecordLength", &recordlength, "RecordLength/i");
        infoTree->Branch("PostTrigger", &post_trigger, "PostTrigger/i");
        for (int i = 0; i < nchannels; i++) {
            int ch = activechannels[i];
            infoTree->Branch(Form("OffsetValue%d", ch), &offset_value2[ch], Form("OffsetValue%d/i", ch));
            infoTree->Branch(Form("TriggerValue%d", ch), &trigger_value2[ch], Form("TriggerValue%d/i", ch));
        }
        infoTree->Fill();   // restored: fill immediately at run start

        tree = new TTree("T", "Digitizer Raw Data");
        tree->Branch("EventNumber", &nevents, "EventNumber/I");
        tree->Branch("TriggerTimeTag", &TriggerTimeTag, "TriggerTimeTag/i");
        tree->Branch("ADC", adc_buffer, Form("ADC[8][%d]/s", recordlength));
    } else { 
        stream.open(filename.c_str());
        if (!stream.is_open()) {
            cout << "[ERROR] CSV File could not be opened." << endl; 
            delete[] adc_buffer;
            return 1;
        }
        
        stream << "Active channels:";
        for (int i = 0; i < nchannels; i++) stream << "CH" << activechannels[i] << ",";
        stream << "\nDevice info:\nModel: " << BoardInfo.Model << "\nChannel number: " << BoardInfo.Channels << "\nSerial number: " << BoardInfo.SerialNumber << "\nADC Bits: " << BoardInfo.ADC_NBits << "\nRecord length: " << recordlength << "\n";
        
        if (offsetOnOff == 1) {
            for (int k = 0; k < nchannels; k++) stream << "Offset value - CH" << activechannels[k] << ": " << offset_value2[activechannels[k]] << "\n";
        }
        if (trgOnOff == 1) {
            for (int k = 0; k < nchannels; k++) stream << "Trigger value - CH" << activechannels[k] << ": " << trigger_value2[activechannels[k]] << "\n";
        }
        stream << "Post trigger: " << post_trigger << "\n\n";
    }

    // Hardware Initialization (Safely processed without broken success counters)
    bool config_ok = true;
    if (CAEN_DGTZ_Reset(handle) != 0) config_ok = false;
    if (CAEN_DGTZ_SetRecordLength(handle, recordlength) != 0) config_ok = false;
    if (CAEN_DGTZ_SetChannelEnableMask(handle, uint32_t(bitset<8>(channel_input).to_ulong())) != 0) config_ok = false;
    if (CAEN_DGTZ_SetPostTriggerSize(handle, post_trigger) != 0) config_ok = false;
    if (CAEN_DGTZ_SetAcquisitionMode(handle, CAEN_DGTZ_SW_CONTROLLED) != 0) config_ok = false;

    CAEN_DGTZ_TriggerPolarity_t polarity;
    polarity = trigger_risingedge ? CAEN_DGTZ_TriggerOnRisingEdge : CAEN_DGTZ_TriggerOnFallingEdge;

    for (int k = 0; k < nchannels; k++) {
        uint32_t ch = activechannels[k];
        uint32_t ch_offset = (offsetOnOff == 1 ? offset_value2[ch] : 0x7FFF);
        uint32_t ch_trigger = (trgOnOff == 1 ? trigger_value2[ch] : 4095);

        if (CAEN_DGTZ_SetChannelDCOffset(handle, ch, ch_offset) != 0) config_ok = false;
        if (CAEN_DGTZ_SetTriggerPolarity(handle, ch, polarity) != 0) config_ok = false;
        if (CAEN_DGTZ_SetChannelTriggerThreshold(handle, ch, ch_trigger) != 0) config_ok = false;
    }

    if (external_trigger) {
        if (CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY) != 0) config_ok = false;
        if (CAEN_DGTZ_SetChannelSelfTrigger(handle, CAEN_DGTZ_TRGMODE_DISABLED, uint32_t(bitset<8>(channel_input).to_ulong())) != 0) config_ok = false;
    } else {
        if (CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED) != 0) config_ok = false;
        if (CAEN_DGTZ_SetChannelSelfTrigger(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY, uint32_t(bitset<8>(channel_input).to_ulong())) != 0) config_ok = false;
    }

    if (CAEN_DGTZ_SetSWTriggerMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY) != 0) config_ok = false;

    if (!config_ok) {
        cout << "Configuration not successful." << endl;   
        if (root && MyFile) MyFile->Close(); 
        else if (!root) stream.close();
        delete[] adc_buffer;
        return 1;
    }

    // Buffer Allocation
    char* buffer; buffer = NULL; uint32_t buffer_size, bsize;
    CAEN_DGTZ_MallocReadoutBuffer(handle, &buffer, &buffer_size);

    cout << "\n[INFO] Waiting 3 seconds for DAC stabilization..." << endl;
    sleep(3);

    // Start DAQ
    CAEN_DGTZ_SWStartAcquisition(handle);
    auto start_time = chrono::high_resolution_clock::now();
    CAEN_DGTZ_EventInfo_t eventinfo;
    void *event; event = NULL;

    // Output message matched exactly with test6 format
    cout << "[INFO] DAQ Started (Mode: " << b_RunMode << ", Format: 16-bit Optimized)" << endl;

    std::vector<double> ped_sum_1000(8, 0.0);
    int ped_calc_samples; ped_calc_samples = (recordlength > 50) ? 50 : recordlength;

    int read_errors = 0;
    while (1) {
        // Check the readout return code. On error, bsize may be left at its previous
        // (stale) value; falling through would re-decode the SAME buffer and write
        // duplicate/corrupt events. Skip the cycle instead, and abort if the board
        // keeps failing so a long automated scan can't hang forever on a dead link.
        CAEN_DGTZ_ErrorCode read_ret =
            CAEN_DGTZ_ReadData(handle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer, &bsize);
        if (read_ret != CAEN_DGTZ_Success) {
            if (++read_errors > 1000) {
                cout << "\n[ERROR] CAEN_DGTZ_ReadData failed repeatedly (code "
                     << read_ret << "). Aborting acquisition to avoid a hang." << endl;
                break;
            }
            usleep(1000);
            continue;
        }
        read_errors = 0;
        if (bsize == 0) { usleep(1000); continue; }

        uint32_t numEv;
        CAEN_DGTZ_GetNumEvents(handle, buffer, bsize, &numEv);
        if (numEv == 0) continue;

        for (uint32_t e = 0; e < numEv; e++) {
            char* eventptr;
            CAEN_DGTZ_GetEventInfo(handle, buffer, bsize, e, &eventinfo, &eventptr);
            CAEN_DGTZ_DecodeEvent(handle, eventptr, &event);
            CAEN_DGTZ_UINT16_EVENT_t *event16 = (CAEN_DGTZ_UINT16_EVENT_t *)event;
            TriggerTimeTag = eventinfo.TriggerTimeTag;

            for (int k = 0; k < nchannels; k++) {
                int ch = activechannels[k];
                for (int s = 0; s < (int)recordlength; s++) {
                    adc_buffer[ch * recordlength + s] = (uint16_t)event16->DataChannel[ch][s];
                }

                double event_ped = 0;
                for (int s = 0; s < ped_calc_samples; ++s) {
                    event_ped += event16->DataChannel[ch][s];
                }
                ped_sum_1000[ch] += (event_ped / (double)ped_calc_samples);
            }

            if (root) {
                tree->Fill();
            } else { 
                stream << "Event " << (nevents + 1) << "\n";
                stream << "Event info:\nEvent size: " << eventinfo.EventSize << " samples.\nChannel mask: " << eventinfo.ChannelMask << "\nEvent counter: " << eventinfo.EventCounter << "\nTime tag: " << TriggerTimeTag << "\n";

                for (int s = 0; s < (int)recordlength; s++) {
                    for (int k = 0; k < nchannels; k++) {
                        stream << adc_buffer[activechannels[k] * recordlength + s] << ",";
                    }
                    stream << "\n";
                }
                stream << "\n\n\n\n";
            } 

            nevents++;

            if (nevents % 1000 == 0) {
                cout << "\033[1;32m[Evt: " << setw(6) << nevents << "]\033[0m Pedestal -> ";
                for (int k = 0; k < nchannels; k++) {
                    int ch = activechannels[k];
                    double mean_val = ped_sum_1000[ch] / 1000.0;
                    cout << "Ch" << ch << ":" << fixed << setprecision(1) << mean_val << " | ";
                    ped_sum_1000[ch] = 0.0; 
                }
                cout << std::endl;
            }

            CAEN_DGTZ_FreeEvent(handle, &event);
            if (nevents_user != 0 && nevents >= (int)nevents_user) break;
        }
        if (nevents_user != 0 && nevents >= (int)nevents_user) break;
    }

    cout << endl; 

    CAEN_DGTZ_SWStopAcquisition(handle);

    if (root && MyFile) {
        MyFile->Write();
        MyFile->Close();
    } else if (!root) {
        stream.close();
    }

    CAEN_DGTZ_CloseDigitizer(handle);
    CAEN_DGTZ_FreeReadoutBuffer(&buffer);

    auto end_time = chrono::high_resolution_clock::now();
    // [DEVELOP] double elapsed_s = chrono::duration<double>(end_time - start_time).count();
    cout << "\n\n[DONE] DAQ Finished." << endl;
    cout << "  - Total Events : " << nevents << endl;
    cout << "  - Elapsed Time : " << chrono::duration<double>(end_time - start_time).count() << " s" << endl;

    delete[] adc_buffer;
    return 0;
}
