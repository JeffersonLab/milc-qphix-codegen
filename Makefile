
DRIVER=dslash_xyzt 
milc=1
ifeq ($(milc),1)
DRIVER_CG=ks_dslash_xyzt
endif

mode=mic

mode:=$(strip $(mode))

CONFFILE=customMake.$(mode)
include $(CONFFILE)

CXXHOST  = icpc -O3 -g

ifeq ($(mode),mic)

ifeq ($(PRECISION),1)
override VECLEN=16
else
override VECLEN=8
endif
yesnolist += AVX512
endif

ifeq ($(mode),avx)
ifeq ($(PRECISION),1)
override VECLEN=8
else
override VECLEN=4
endif
#DEFS += -DNO_HW_MASKING
yesnolist += AVX2
endif

ifeq ($(mode),sse)
ifeq ($(PRECISION),1)
override VECLEN=4
override SOALEN=4
else
override VECLEN=2
override SOALEN=2
endif
#DEFS += -DNO_HW_MASKING
yesnolist += NO_MASKS
endif

ifeq ($(mode),scalar)
override VECLEN=1
override SOALEN=1
#DEFS += -DNO_HW_MASKING
DEFS += -DNO_MASKS
endif


# If streaming stores are enabled we
# should definitely disable prefetching 
# of output spinors
ifeq ($(ENABLE_STREAMING_STORES),1)
PREF_L1_SPINOR_OUT=0
PREF_L2_SPINOR_OUT=0
endif

ifeq ($(ENABLE_LOW_PRECISION),1)
USE_LP_SPINOR=1
USE_LP_GAUGE=1
USE_LP_CLOVER=1
endif


yesnolist += PREF_L1_SPINOR_IN 
yesnolist += PREF_L2_SPINOR_IN 
yesnolist += PREF_L1_SPINOR_OUT 
yesnolist += PREF_L2_SPINOR_OUT 
yesnolist += PREF_L1_GAUGE
yesnolist += PREF_L2_GAUGE
yesnolist += PREF_L1_CLOVER
yesnolist += PREF_L2_CLOVER
yesnolist += USE_LDUNPK
yesnolist += USE_PKST
yesnolist += USE_PACKED_GAUGES
yesnolist += USE_PACKED_CLOVER
yesnolist += USE_SHUFFLES
yesnolist += NO_GPREF_L1
yesnolist += NO_GPREF_L2
yesnolist += ENABLE_STREAMING_STORES

yesnolist += USE_LP_SPINOR
yesnolist += USE_LP_GAUGE
yesnolist += USE_LP_CLOVER
yesnolist += SERIAL_SPIN
yesnolist += TESTFLAG

deflist += SOALEN
deflist += VECLEN
deflist += PRECISION

DEFS += $(strip $(foreach var, $(yesnolist), $(if $(filter 1, $($(var))), -D$(var))))
DEFS += $(strip $(foreach var, $(deflist), $(if $($(var)), -D$(var)=$($(var)))))

SOURCES = data_types_xyzt.cc data_utils.cc complex.cc complex_str.cc mat_utils.cc mat_utils_str.cc dslash_common.cc inst_dp_vec8.cc inst_sp_vec16.cc inst_dp_vec4.cc inst_sp_vec8.cc inst_sp_vec4.cc inst_dp_vec2.cc inst_scalar.cc

SOURCES_CG = $(SOURCES) $(DRIVER_CG).cc

HEADERS_CG = address_types.h  data_types_xyzt.h  dslash_common.h  instructions.h data_utils.h complex.h mat_utils.h  Makefile $(CONFFILE)


all: codegen_cg

codegen_cg: $(SOURCES_CG) $(HEADERS_CG) 
	$(CXXHOST) $(DEFS) $(SOURCES_CG) -o ./codegen_cg


.PHONY: cgen mic avx avx2 avx512 sse scalar

cgen: mic avx avx2 avx512 sse scalar

mic:
	mkdir -p ./mic
	@make clean && make -f Makefile mode=mic PRECISION=2 && ./codegen_cg 
	@make clean && make -f Makefile mode=mic PRECISION=1 && ./codegen_cg 
#	@make clean && make -f Makefile mode=mic PRECISION=1 ENABLE_LOW_PRECISION=1 && ./codegen

avx512:
	mkdir -p ./avx512
	@make clean && make -f Makefile mode=mic AVX512=1 PRECISION=2 && ./codegen_cg
	@make clean && make -f Makefile mode=mic AVX512=1 PRECISION=1 && ./codegen_cg
#	@make clean && make -f Makefile mode=mic AVX512=1 PRECISION=1 ENABLE_LOW_PRECISION=1 && ./codegen


avx:
	mkdir -p ./avx
	@make clean && make -f Makefile mode=avx PRECISION=2 && ./codegen_cg 
	@make clean && make -f Makefile mode=avx PRECISION=1 && ./codegen_cg 

avx2:
	mkdir -p ./avx2
	@make clean && make -f Makefile mode=avx PRECISION=2 AVX2=1 && ./codegen_cg 
	@make clean && make -f Makefile mode=avx PRECISION=1 AVX2=1 && ./codegen_cg 
#	@make clean && make -f Makefile mode=avx PRECISION=1 AVX2=1 ENABLE_LOW_PRECISION=1 && ./codegen

sse:
	mkdir -p ./sse
	@make clean && make mode=sse PRECISION=2 && ./codegen_cg
	@make clean && make mode=sse PRECISION=1 && ./codegen_cg 

#	@make clean && make mode=sse PRECISION=2 && ./codegen
#	@make clean && make mode=sse PRECISION=1 && ./codegen

scalar:
	mkdir -p ./scalar
	@make clean && make mode=scalar PRECISION=2 && ./codegen_cg 
	@make clean && make mode=scalar PRECISION=1 && ./codegen_cg 

clean: 
	rm -rf *.o ./codegen*

cleanall: 
	rm -rf *.o ./codegen*
	rm -rf ./avx 
	rm -rf ./avx2
	rm -rf ./avx512
	rm -rf ./mic
	rm -rf ./sse
	rm -rf ./scalar
