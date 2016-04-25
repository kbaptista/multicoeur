CC := gcc
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:= sand

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

$(PROG): sand.o display.o treatment_doubleline.o treatment_checkboard.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS) 

clean:
	rm -rf *.o $(PROG) $(COMPARER)





