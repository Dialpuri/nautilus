# nautilus makefile - Jordan Dialpuri Edit 12/10/22

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
-g

INCDIR = -I${CCP4}/include
LIBDIR = -L${CCP4}/lib
CFLAGS = ${FLAGS} ${INCDIR} -c
LFLAGS = ${FLAGS} ${LIBDIR} ${LIBS}

SRC_DIR = src
BIN_DIR = bin

OBJS = $(SHARED:=.o)

TARGET_OBJS = $(addprefix ${BIN_DIR}/,${OBJS})
TARGET_SRC = $(addprefix ${SRC_DIR}/,${OBJS})

BUILD_PRINT = @echo "\e[1;34mBuilding $<\e[0m"
COMPLETE_PRINT = @echo "\e[1;32mBuilding complete!\e[0m"

MKDIR_P = mkdir -p

#Source CCP4
IGNORE := $(shell bash -c "source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
include makeenv

$(shell mkdir -p $(BIN_DIR))

all: cnautilus

cnautilus: ${TARGET_SRC} cnautilus.o
	$(BUILD_PRINT)
	@g++ ${TARGET_OBJS} ${BIN_DIR}/cnautilus.o -o $@ ${LFLAGS}
	$(COMPLETE_PRINT)

cnautilus.o: ${SRC_DIR}/cnautilus.cpp
	$(BUILD_PRINT)
	@g++ ${CFLAGS} $< -o ${BIN_DIR}/$(@F)

%.o: %.cpp %.h
	$(BUILD_PRINT)
	@g++ ${CFLAGS} $< -o ${BIN_DIR}/$(@F)

clean:
	rm bin/*.o cnautilus
	rm makeenv
