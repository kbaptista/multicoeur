
ARCH    := $(shell uname -s | tr a-z A-Z)

PROG	:=	sable

PROG_SEQ  := seq

CFLAGS	:=	-g -O3 -std=c99 -Wno-deprecated-declarations
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS		:= -lm -lOpenCL -lGL -lGLU -lglut
endif

.phony: default clean

default: $(PROG) $(PROG_SEQ)

$(PROG): main.o display.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

$(PROG_SEQ): seq.o display.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

main.o: display.h

display.o: display.h

clean:
	rm -rf *.o $(PROG) $(PROG_SEQ)
