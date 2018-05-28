#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013--2018, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Mingjian Wen


# load all basic KIM make configuration
KIM_API_BUILD_CONFIG = kim-api-v2-build-config
ifeq ($(shell $(KIM_API_BUILD_CONFIG) --version 2> /dev/null),)
  $(error $(KIM_API_BUILD_CONFIG) utility is not available.  Something is wrong with your KIM API package setup)
endif
include $(shell $(KIM_API_BUILD_CONFIG) --master-config)


# set model driver specific details
MODEL_DRIVER_NAME := RDP__MD_000000111111_000
MODEL_DRIVER_CREATE_FUNCTION_NAME := model_driver_create
MODEL_DRIVER_CREATE_FUNCTION_LANG := cpp

LOCALOBJ = RDP.o RDPImplementation.o helper.o

RDP.o: RDP.hpp RDPImplementation.hpp
RDPImplementation.o: RDPImplementation.hpp helper.h  \
                                 RDPImplementationComputeDispatch.cpp
helper.o: helper.h
RDPImplementationComputeDispatch.cpp: CreateDispatch.sh
	@./CreateDispatch.sh
	@printf "Creating... $@.\n"

LOCALCLEAN =  RDPImplementationComputeDispatch.cpp

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=
#LDLIBS   +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ModelDriver
