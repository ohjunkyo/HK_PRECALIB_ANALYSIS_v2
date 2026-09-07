# Set compiler
CXX = g++

CXXFLAGS = -std=c++11

CAEN_INCLUDE = -I/home/precalkor/ADC/include

#CAEN_LIBS_PATH = -L/home/precalkor/ADC/lib
CAEN_LIBS_PATH = -L/home/precalkor/ADC/lib -Wl,-rpath,/home/precalkor/ADC/lib

CAEN_LIBS = -lCAENComm -lCAENVME -lCAENDigitizer

OTHER_LIBS = -lboost_program_options

ROOT_FLAGS = $(shell root-config --cflags --glibs)

#TARGET = testrun # Execute program
#TARGET = execute_DAQ # Execute program
TARGET = execute_DAQ_v2 # Execute program


#SRC = ADC_test3.cpp # <== Edit this file /
#SRC = ADC_test6.cpp # <== Edit this file v.2026. 03. 17 by Junkyo
#SRC = ADC_test7.cpp # <== v.2026. 06. 10 by Junkyo for .csv file
SRC = ADC_test8_20260808.cpp # <== v.2026. 08. 08: TriggerValue<N> moved off the
# per-event tree "T" into "RunInfo" (it is set once before acquisition and never
# changes, so writing it into all 300k events was pure waste), and the channel
# mask is now saved explicitly as RunInfo/ChannelMask. RecordLength, PostTrigger
# and OffsetValue<N> deliberately STAY on "T": prod_ntp_v7.C reads the first two
# from there and infers the active channel list from OffsetValue<N>'s branch
# existence, so moving them would break analysis of both new and old files.


all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(CAEN_INCLUDE) $(CAEN_LIBS_PATH) $(CAEN_LIBS) $(OTHER_LIBS) $(ROOT_FLAGS)

clean:
	rm -f $(TARGET)
