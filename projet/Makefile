
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:= sand

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

$(PROG): seq.o parallel.o sand.o display.o
	$(CC) -o $@ $(LDFLAGS) -fopenmp $^ $(LDLIBS) 

sand.o: display.h seq.h parallel.h

seq.o: seq.h

parallel.o: parallel.h

display.o: display.h

clean:
	rm -rf *.o $(PROG) $(PROG_SEQ) $(PROG_PARALLEL) 
