CC := gcc
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:= sand
#COMPARER := comparer

CFLAGS	:=	-g -std=c99 -Wno-deprecated-declarations -fopenmp
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include 
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS	:= -lm -lOpenCL -lGL -lGLU -lglut -fopenmp
endif

.phony: default clean

default: $(PROG)

$(PROG): sand.o display.o treatment.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS) 

#$(COMPARER): seq.o parallel.o comparer.o display.o parallel_task.o
#	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS) 


#sand.o: display.h seq.h parallel.h parallel_task.h
#comparer.o: display.h seq.h parallel.h parallel_task.h

clean:
	rm -rf *.o $(PROG) $(COMPARER)





