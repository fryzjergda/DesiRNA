CXXFLAGS =	-O2 -pipe -march=native
LIBS = 

LDFLAGS =

OBJS =  SS_record_and_matrix.o SS_consensus_and_output.o find_SS_consensus.o functions.o

TARGET = find_SS_consensus

all:	find_SS_consensus

find_SS_consensus: $(OBJS)
	$(CXX) -o find_SS_consensus $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
