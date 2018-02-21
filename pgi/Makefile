#Variables--------------------------------------------------------|
SRCS = main.cpp vec.cpp mden.cpp butter.cpp ts2means.cpp
OBJS = $(SRCS:.cpp=.o)
CC = g++				#C++ compiler
CFLAGS = -c -fopenmp			#Compiler flags
LFLAGS = -lsndfile -fopenmp -lfftw3	#linker flags
TARGET = MDEN
#-----------------------------------------------------------------|

#Build using g++--------------------------------------------------|
mden: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o mden libwavelib2.a

main.o: main.cpp 
	$(CC) $(CFLAGS) main.cpp

vec.o: vec.cpp vec.h
	$(CC) $(CFLAGS) vec.cpp

mden.o: mden.cpp mden.h
	$(CC) $(CFLAGS) mden.cpp

mbutter.o: butter.cpp
	$(CC) $(CFLAGS) butter.cpp

ts2means.o: ts2means.cpp ts2means.h
	$(CC) $(CFLAGS) ts2means.cpp
#-----------------------------------------------------------------|

#Build using OpenACC compiler and target-specific variable values-|
mdenACC: CC = pgc++
mdenACC: CFLAGS = -c -fast -ta=tesla -Minfo=accel
mdenACC: LFLAGS = -lsndfile -lfftw3 -ta=tesla

mdenACC: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o mdenACC libwavelib2.a
#-----------------------------------------------------------------|

#-----------------------------------------------------------------|

#Debug-Build using OpenACC cmplr and target-specific variable val-|
dbgACC: CC = pgc++
dbgACC: CFLAGS = -c -g -ta=tesla -Minfo=accel
dbgACC: LFLAGS = -lsndfile -lfftw3 -ta=tesla

dbgACC: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o dbgACC libwavelib2.a
#-----------------------------------------------------------------|


#Build OpenACC test file
testACC: mexa.cpp butter.cpp
	pgc++ -lsndfile -fast -ta=tesla -Minfo=accel  butter.cpp mexa.cpp -o testACC libwavelib2.a

#Build pgi forum mexa version 1 (just parallelizing)
pgimexa1: pgimexa.1.cpp
	pgc++ -lsndfile -fast -ta=tesla -Minfo=accel pgimexa.1.cpp -o pgimexa1 libwavelib2.a

#Build pgi forum mexa version 2 (using async)
pgimexa2: pgimexa.2.cpp butter.cpp ts2means.cpp ts2means.h vec.cpp
	pgc++ -lsndfile -fast -ta=tesla -Minfo=accel vec.cpp ts2means.cpp butter.cpp pgimexa.2.cpp -o pgimexa2 libwavelib2.a


clean: 
	rm -f *.o *.gch mden mdenACC dbgACC testACC pgimexa1 pgimexa2
