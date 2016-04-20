
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:= sand
COMPARER := comparer

CFLAGS	:=	-g -O3 -std=c99 -Wno-deprecated-declarations
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS	:= -lm -lOpenCL -lGL -lGLU -lglut
endif

.phony: default clean

default: $(PROG)

$(PROG): seq.o parallel.o sand.o display.o parallel_task.o
	$(CC) -o $@ $(LDFLAGS) -fopenmp $^ $(LDLIBS) 

$(COMPARER): seq.o parallel.o comparer.o display.o parallel_task.o
	$(CC) -o $@ $(LDFLAGS) -fopenmp $^ $(LDLIBS) 


sand.o: display.h seq.h parallel.h parallel_task.h
comparer.o: display.h seq.h parallel.h parallel_task.h

seq.o: seq.h

parallel.o: parallel.h

parallel_task.o: parallel_task.h

display.o: display.h

clean:
	rm -rf *.o $(PROG) $(PROG_SEQ) $(PROG_PARALLEL) 
