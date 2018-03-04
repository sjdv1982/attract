SHELL := /bin/bash
# Find out the base directory

.PHONY: default
default: all

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

BASE_PATH = $(shell pwd)
BINDIR = $(BASE_PATH)/.
ATTRACTDIR = $(BASE_PATH)/.
GPU_DIR = $(BASE_PATH)/gpuATTRACT
GPU_DIR2 = ./gpuATTRACT

$(shell mkdir -p $(BINDIR))

.PHONY: help all gpuATTRACT emATTRACT mcATTRACT scATTRACT attract \
 clean_as clean_em clean_mc clean_sc clean_attract \
 cleanall_as cleanall_em cleanall_mc cleanall_sc cleanall_attract clean cleanall

em: as
	cd $(GPU_DIR)/emATTRACT/bin && $(MAKE) TARGET=$(TARGET)
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/emATTRACT/bin/emATTRACT

gpuATTRACTcheck: as
	cd $(GPU_DIR)/gpuATTRACTcheck/bin && $(MAKE) TARGET=$(TARGET)
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/gpuATTRACTcheck/bin/gpuATTRACTcheck

mc: as
	cd $(GPU_DIR)/mcATTRACT/bin && $(MAKE) TARGET=$(TARGET)
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/mcATTRACT/bin/mcATTRACT
	
sc: as
	cd $(GPU_DIR)/scATTRACT/bin && $(MAKE) TARGET=$(TARGET)
	cd $(BINDIR) && ln -sf $(GPU_DIR)/scATTRACT/bin/scATTRACT

as: check_cuda_version
	cd $(GPU_DIR)/AttractServer/lib && $(MAKE) TARGET=$(TARGET)
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/AttractServer/lib/libAttractServer.so

attract:
	cd $(ATTRACTDIR) && $(MAKE) -f _Makefile all
	

NVCC_RESULT = $(shell which nvcc 2> /dev/null)
NVCC_TEST = $(notdir $(NVCC_RESULT))

MIN_VERSION="6.5"
check_cuda_version:
	@echo "Makefile: Testing CUDA requirements"
	@if [ "$(NVCC_TEST)" == "nvcc" ]; then \
		echo -e >&2 "Makefile:\tOK: nvcc compiler found"; \
	else \
		echo -e >&2 "Makefile:\tERROR: nvcc compiler not found"; \
		false; \
	fi 
			
	@NVCC_VERSION=$$(nvcc --version | grep release | cut -d, -f2 | cut -d" " -f3); \
	if (( $$(echo "$$NVCC_VERSION >= $(MIN_VERSION)" | bc -l) )); then \
		echo -e >&2 "Makefile:\tOK: Compiler version $$NVCC_VERSION"; \
	else \
		echo -e >&2 "Makefile:\tERROR: Compiler version ($$NVCC_VERSION) too old. Upgrade to a version >= $(MIN_VERSION)"; \
		false; \
	fi

all: attract gpuATTRACT

gpuATTRACT: em mc sc gpuATTRACTcheck

clean: clean_gpu clean_attract
	
clean_as:
	@cd $(GPU_DIR)/AttractServer/lib && $(MAKE) clean
	
clean_em:
	@cd $(GPU_DIR)/emATTRACT/bin && $(MAKE) clean
	
clean_mc:
	@cd $(GPU_DIR)/mcATTRACT/bin && $(MAKE) clean
	
clean_sc:
	@cd $(GPU_DIR)/scATTRACT/bin && $(MAKE) clean

clean_attract:
	@cd $(ATTRACTDIR) && rm -f *.o

clean_gpu: clean_as clean_em clean_mc clean_sc
 
cleanall: cleanall_as cleanall_em cleanall_mc cleanall_sc cleanall_attract
	
cleanall_as:
	@cd $(GPU_DIR)/AttractServer/lib && make cleanall

	@rm -f $(BINDIR)/libAttractServer.so
	
cleanall_em:
	@cd $(GPU_DIR)/emATTRACT/bin && make cleanall

	@rm -f $(BINDIR)/emATTRACT
	
cleanall_mc:
	@cd $(GPU_DIR)/mcATTRACT/bin && make cleanall

	@rm -f $(BINDIR)/mcATTRACT
	
cleanall_sc:
	@cd $(GPU_DIR)/scATTRACT/bin && make cleanall

	@rm -f $(BINDIR)/scATTRACT

cleanall_attract:
	cd $(ATTRACTDIR) && make -f _Makefile clean

help:	
	@echo "Options:"
	@echo "make TARGET=DEBUG | RELEASE     choose target. (Default: RELEASE)"
	@echo
	@echo "e.g. make all TARGET=DEBUG -j4"
	@echo
	@echo "***Build targets:"
	@echo "make all            build all targets (attract and gpuATTRACT)"
	@echo "make attract        build original attract (CPU only)"
	@echo "make gpuATTRACT     build gpuATTRACT (libAttractServer.so, emATTRACT, mcATTRACT, scATTRACT)"
	@echo "make as             build libAttractServer.so"
	@echo "make em             build libAttractServer.so and emATTRACT."
	@echo "make mc             build libAttractServer.so and mcATTRACT."
	@echo "make sc             build libAttractServer.so and scATTRACT."
	@echo " "
	@echo "***Clean up: delete .o and .d files of ... " 
	@echo "make clean          all targets."
	@echo "make clean_gpuATTRACT" 
	@echo "                    gpuATTRACT"
	@echo "make clean_attract  attract"     
	@echo "make clean_as       libAttractServer.so"
	@echo "make clean_em       emATTRACT."
	@echo "make clean_mc       mcATTRACT."
	@echo "make clean_sc       scATTRACT."
	@echo " "
	@echo "***Clean up all: as Clean up + removal of executables of ... " 
	@echo "make cleanall       all targets."
	@echo "make cleanall_attract"
	@echo "                    attract"
	@echo "make cleanall_gpuATTACT"
	@echo "                    gpuATTRACT"
	@echo "make cleanall_as    libAttractServer.so"
	@echo "make cleanall_sc    scATTRACT."
	@echo "make cleanall_em    emATTRACT."
	@echo "make cleanall_mc    mcATTRACT."
	
	
