AR = ar
CXX = g++
LIBS = -lginac -lcln -lgmp -lmpfr -lmpc
#INCLUDE = -I/home/jgarcia/boost_1_55_0 -I. 
INCLUDE = -I.
#FLAGS = -O2 -w
FLAGS = -O0 -g -w
OBJS = read_input.o 
HEADERS = problem.hpp options.hpp ricpad.hpp read_input.hpp

all: ricpad

ricpad: $(HEADERS) ricpad.cpp $(OBJS)
	$(CXX) $(FLAGS) -o ricpad ricpad.cpp $(INCLUDE) $(OBJS) $(LIBS) 

libricpad.a: $(OBJS)
	$(AR) rvs libricpad.a $(OBJS) 

clean: 
	rm -rf *.o *.a test test1 test2 test3 

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) -o $@ -c $< $(LIBS) $(INCLUDE)

%.hpp: ;
