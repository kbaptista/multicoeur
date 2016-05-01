CC := gcc
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:= sand

CFLAGS	:=	-g  -O3 -std=c99 -Wno-deprecated-declarations -fopenmp
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include 
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS	:= -lm -lOpenCL -lGL -lGLU -lglut -fopenmp
endif

.phony: default clean

default: $(PROG)

clean:
	rm -rf *.o $(PROG) $(COMPARER)

$(PROG): treatment.o display.o treatment.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS) 

$(EXEC): util.o display.o



