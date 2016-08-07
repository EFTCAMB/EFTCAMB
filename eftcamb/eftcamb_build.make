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
EFTCAMB_SOURCES_FILES := $(wildcard $(EFTCAMB_DIR)/*.f90)
EFTCAMB_SOURCES_FILES +=$(wildcard $(EFTCAMB_DIR)/*.f)

# replace the relevant suffix:
EFT_TEMP_1 := $(patsubst %.f90,%,$(EFTCAMB_SOURCES_FILES))
EFT_TEMP_1 := $(patsubst %.f,%,$(EFT_TEMP_1))
EFT_TEMP_2 := $(notdir $(EFT_TEMP_1))
EFT_FILES  := $(addsuffix .o, $(EFT_TEMP_2))

# files for eftcamb:
EFT_OBJ = $(addprefix $(OUTPUT_DIR)/, $(EFT_FILES))
EFT_SOJ = $(addprefix $(DLL_DIR)/, $(EFT_FILES))

# include EFTCAMB dependencies:
include $(EFTCAMB_DIR)/eftcamb.dep

# how to build EFTCAMB files:
$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/%.f90
	$(F90C) $(F90FLAGS) -c $(EFTCAMB_DIR)/$*.f90 -o $(OUTPUT_DIR)/$*.o

$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/%.f90
	$(F90C) $(SF90FLAGS) -c $(EFTCAMB_DIR)/$*.f90 -o $(DLL_DIR)/$*.o

$(OUTPUT_DIR)/%.o: $(EFTCAMB_DIR)/%.f
	$(F90C) $(F90FLAGS) -c $(EFTCAMB_DIR)/$*.f -o $(OUTPUT_DIR)/$*.o

$(DLL_DIR)/%.o: $(EFTCAMB_DIR)/%.f
	$(F90C) $(SF90FLAGS) -c $(EFTCAMB_DIR)/$*.f -o $(DLL_DIR)/$*.o

# add the EFTCAMB files to the standard ones:
CAMBOBJ += $(EFT_OBJ)
CAMBSO  += $(EFT_SOJ)

# additional EFTCAMB targets:

.PHONY: eftcamb
eftcamb: camb

.PHONY: eftcamb_doc
eftcamb_doc:
	@cd $(EFTCAMB_DIR)/../eftcamb_doc && doxygen EFTCAMB_doxyfile
