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
# Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    Ellad B. Tadmor
#    Valeriu Smirichinski
#
#
#

# load all basic KIM make configuration
include ../Makefile.KIM_Config

# set model driver specific details
MODEL_DRIVER_NAME := RDP_Kolmogorov_Crespi__MD_000000111111_000
MODEL_DRIVER_KIM_FILE_TEMPLATE := RDP_Kolmogorov_Crespi.kim.tpl
MODEL_DRIVER_INIT_FUNCTION_NAME := model_driver_init

LOCALOBJ = RDP_Kolmogorov_Crespi.o

LOCALCLEAN =

# APPEND to compiler option flag lists
#FFLAGS   +=
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.ModelDriver
