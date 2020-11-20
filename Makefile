AR = ar
CXX = g++
LIBS = -lginac -lcln -lgmp -lmpfr
INCLUDE = -I/home/jgarcia/boost_1_55_0 -I. 
#FLAGS = -O2 -w
FLAGS = -O0 -g
OBJS = read_input.o ricpad.o

all: test

test: libricpad.a test.o
	$(CXX) $(FLAGS) -o test test.o $(INCLUDE) $(OBJS) $(LIBS) 

libricpad.a: $(OBJS)
	$(AR) rvs libricpad.a $(OBJS) 

clean: 
	rm -rf *.o *.a test1

%.o: %.cpp ricpad.hpp
	$(CXX) $(FLAGS) -o $@ -c $< $(LIBS) $(INCLUDE)
