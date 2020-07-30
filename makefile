CXX = g++
CXXFLAGS = -Wall -g -O0 -std=c++11
SRCDIR = srcs
INCS = -I$(SRCDIR)
LIBS = -lm
LDFLAGS = -fPIC -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO  -lMLP 
TARGET = preprocess
SOURCES = $(TARGET).cxx
SOURCES += $(wildcard $(SRCDIR)/*.cxx)
OBJS = $(SOURCES:.cxx=.o)

TRAIN = BDT_combine_training
SOURCES_train = $(TRAIN).cxx
OBJS_train = $(SOURCES_train:.cxx=.o)

DATA = BDT_combine_data
SOURCES_data = $(DATA).cxx
OBJS_data = $(SOURCES_data:.cxx=.o)


#INCS += -I$(shell root-config --incdir)
CXXFLAGS += $(shell root-config --cflags)  
LDFLAGS += $(shell root-config --libs)

default: $(TARGET) #clean0 --> .o files can be reused if the source codes have no updates 
train: $(TRAIN)
data: $(DATA)

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS)

$(OBJS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) $(INCS) $< -c -o $@

$(TRAIN):$(OBJS_train)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS_train)

$(OBJS_train): %.o : %.cxx
	$(CXX) $(CXXFLAGS) $(INCS) $< -c -o $@

$(DATA):$(OBJS_data)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS_data)

$(OBJS_data): %.o : %.cxx
	$(CXX) $(CXXFLAGS) $(INCS) $< -c -o $@

clean0:
	rm -f $(OBJS)
	rm -f $(OBJS_train)
	rm -f $(OBJS_data)

clean:
	rm -f $(TARGET)
	rm -f $(TRAIN)
	rm -f $(OBJS)
	rm -f $(OBJS_train)
	rm -f $(OBJS_data)
