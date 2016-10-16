#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2016 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

#
# Makefile auxiliary file containing everything that we need for EFTCAMB
#

# get the eftcamb sources:
EFTCAMB_SOURCES_FILES := $(wildcard $(EFTCAMB_DIR)/*.f90) $(wildcard $(EFTCAMB_DIR)/*/*.f90)
EFTCAMB_SOURCES_FILES += $(wildcard $(EFTCAMB_DIR)/*.f) $(wildcard $(EFTCAMB_DIR)/*/*.f)

# replace the relevant suffix:
EFT_TEMP_1 := $(patsubst %.f90,%,$(EFTCAMB_SOURCES_FILES))
EFT_TEMP_1 := $(patsubst %.f,%,$(EFT_TEMP_1))
EFT_TEMP_2 := $(notdir $(EFT_TEMP_1))
EFT_FILES  := $(addsuffix .o, $(EFT_TEMP_2))

# files for eftcamb:
EFT_OBJ = $(addprefix $(OUTPUT_DIR)/, $(EFT_FILES))
EFT_SOJ = $(addprefix $(DLL_DIR)/, $(EFT_FILES))

# how to build EFTCAMB files:
$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/%.f90
	$(F90C) $(F90FLAGS) -c $< -o $(OUTPUT_DIR)/$*.o
$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/*/%.f90
	$(F90C) $(F90FLAGS) -c $< -o $(OUTPUT_DIR)/$*.o

$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/%.f90
	$(F90C) $(SF90FLAGS) -c $< -o $(DLL_DIR)/$*.o
$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/*/%.f90
	$(F90C) $(SF90FLAGS) -c $< -o $(DLL_DIR)/$*.o

$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/%.f
	$(F90C) $(F90FLAGS) -c $< -o $(OUTPUT_DIR)/$*.o
$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/*/%.f
	$(F90C) $(F90FLAGS) -c $< -o $(OUTPUT_DIR)/$*.o

$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/%.f
	$(F90C) $(SF90FLAGS) -c $< -o $(DLL_DIR)/$*.o
$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/*/%.f
	$(F90C) $(SF90FLAGS) -c $< -o $(DLL_DIR)/$*.o

# add the EFTCAMB files to the standard ones:
CAMBOBJ += $(EFT_OBJ)
CAMBSO  += $(EFT_SOJ)

camb: $(EFT_OBJ)

# EFTCAMB dependencies:
eftcamb_dep: $(EFTCAMB_SOURCES_FILES)
	@python $(EFTCAMB_DIR)/depend_gen.py -o $(EFTCAMB_DIR)

# include EFTCAMB dependencies:
include $(EFTCAMB_DIR)/eftcamb.dep

# additional EFTCAMB targets:

.PHONY: eftcamb
eftcamb: camb

.PHONY: eftcamb_doc
eftcamb_doc:
	@cd $(EFTCAMB_DIR)/../eftcamb_doc && doxygen EFTCAMB_doxyfile

# EFTCAMB APPS:

APPS_SOURCES := $(wildcard $(EFTCAMB_APPS)/*.f90)
APPS_SOURCES += $(wildcard $(EFTCAMB_APPS)/*.F90)

APPS_NAMES   := $(notdir $(patsubst %.f90,%,$(APPS_SOURCES)))
APPS_NAMES   := $(notdir $(patsubst %.F90,%,$(APPS_NAMES)))

APPS_BINS    := $(addsuffix .x, $(APPS_NAMES))
APPS_BINS    := $(addprefix $(CAMB_DIR)/, $(APPS_BINS))

.PHONY: eftcamb_apps

eftcamb_apps: directories $(CAMBOBJ) $(APPS_BINS) profile

$(CAMB_DIR)/%.x: directories $(CAMBOBJ) $(EFTCAMB_APPS)/%.f90 
	$(F90C) $(F90FLAGS) $(CAMBOBJ) $(EFTCAMB_APPS)/$*.f90 $(F90CRLINK) -o $(CAMB_DIR)/$*.x

$(CAMB_DIR)/%.x: directories $(CAMBOBJ) $(EFTCAMB_APPS)/%.F90 
	$(F90C) $(F90FLAGS) $(CAMBOBJ) $(EFTCAMB_APPS)/$*.F90 $(F90CRLINK) -o $(CAMB_DIR)/$*.x

# PROFILER: this requires a special target
ifeq ($(MAKECMDGOALS), profile)
F90FLAGS += -pg
endif

profile: directories $(CAMBOBJ) $(EFTCAMB_APPS)/benchmark.F90 
	$(F90C) $(F90FLAGS) $(CAMBOBJ) $(EFTCAMB_APPS)/benchmark.F90 $(F90CRLINK) -o $(CAMB_DIR)/profiler.x

clean_apps:
	@rm -f $(CAMB_DIR)/*.x
	
# add the apps targets to the main ones:
all: eftcamb_apps
clean: clean_apps

