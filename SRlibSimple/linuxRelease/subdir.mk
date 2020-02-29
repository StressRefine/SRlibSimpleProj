################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SRbasis.cpp \
../SRbasisBrickWedge.cpp \
../SRconstraint.cpp \
../SRcoord.cpp \
../SRedge.cpp \
../SRelemBrickWedge.cpp \
../SRelement.cpp \
../SRface.cpp \
../SRfaceQuad.cpp \
../SRfile.cpp \
../SRforce.cpp \
../SRmachDep.cpp \
../SRmap.cpp \
../SRmapBrickWedge.cpp \
../SRmaterial.cpp \
../SRmath.cpp \
../SRmodel.cpp \
../SRnode.cpp \
../SRstring.cpp \
../SRutil.cpp 

OBJS += \
./SRbasis.o \
./SRbasisBrickWedge.o \
./SRconstraint.o \
./SRcoord.o \
./SRedge.o \
./SRelemBrickWedge.o \
./SRelement.o \
./SRface.o \
./SRfaceQuad.o \
./SRfile.o \
./SRforce.o \
./SRmachDep.o \
./SRmap.o \
./SRmapBrickWedge.o \
./SRmaterial.o \
./SRmath.o \
./SRmodel.o \
./SRnode.o \
./SRstring.o \
./SRutil.o 

CPP_DEPS += \
./SRbasis.d \
./SRbasisBrickWedge.d \
./SRconstraint.d \
./SRcoord.d \
./SRedge.d \
./SRelemBrickWedge.d \
./SRelement.d \
./SRface.d \
./SRfaceQuad.d \
./SRfile.d \
./SRforce.d \
./SRmachDep.d \
./SRmap.d \
./SRmapBrickWedge.d \
./SRmaterial.d \
./SRmath.d \
./SRmodel.d \
./SRnode.d \
./SRstring.d \
./SRutil.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


