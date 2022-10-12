# nautilus makefile

SHARED = \
nautilus-join \
nautilus-rebuild-bases \
nautilus-sequence \
nautilus-ss-find \
nautilus-target \
nautilus-tidy \
nautilus-tools \
nautilus-util \
nucleicacid_db \

LIBS = \
-l:libccp4c.so.8.0 \
-l:libclipper-ccp4.so.2 \
-l:libclipper-contrib.so.2 \
-l:libclipper-core.so.2 \
-l:libclipper-minimol.so.2 \
-l:libclipper-mmdb.so.2 \
-l:libfftw.so.2 \
-l:libmmdb2.so.0 \
-l:librfftw.so.2 \
-lm \
-lstdc++ \

FLAGS = \
-std=c++11 \
-O2 \
-Wall \
-Wno-sign-compare \
-fPIC \
-ftemplate-depth-50 \
-D_GLIBCXX_USE_CXX11_ABI=0 \

INCDIR = -I${CCP4}/include
LIBDIR = -L${CCP4}/lib
CFLAGS = ${FLAGS} ${INCDIR} -c
LFLAGS = ${FLAGS} ${LIBDIR} ${LIBS}

SRC_DIR = src
OBJS = $(SHARED:=.o)

TARGET_OBJS = $(addprefix ${SRC_DIR}/,${OBJS})

cnautilus: ${TARGET_OBJS} ${SRC_DIR}/cnautilus.o
	g++ ${TARGET_OBJS}  ${SRC_DIR}/cnautilus.o -o $@ ${LFLAGS}

cnautilus.o: cnautilus.cpp
	g++ ${CFLAGS} $<

%.o: %.cpp %.h
	g++ ${CFLAGS} $<

clean:
	rm *.o cnautilus
