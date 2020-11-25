AR = ar
CXX = g++
LIBS = -lginac -lcln -lgmp -lmpfr -lmpc
INCLUDE = -I/home/jgarcia/boost_1_55_0 -I. 
#FLAGS = -O2 -w
FLAGS = -O0 -g -w
OBJS = read_input.o 
HEADERS = problem.hpp options.hpp

all: test

test3: test3.cpp
	$(CXX) -o test3 test3.cpp -lmpc -lmpfr -lgmp $(INCLUDE)

test2: $(HEADERS) test2.cpp $(OBJS)
	$(CXX) $(FLAGS) -o test2 test2.cpp $(INCLUDE) $(OBJS) $(LIBS) 

test: libricpad.a test.o
	$(CXX) $(FLAGS) -o test test.o $(INCLUDE) $(OBJS) $(LIBS) 

libricpad.a: $(OBJS)
	$(AR) rvs libricpad.a $(OBJS) 

clean: 
	rm -rf *.o *.a test1

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) -o $@ -c $< $(LIBS) $(INCLUDE)

%.hpp: ;
