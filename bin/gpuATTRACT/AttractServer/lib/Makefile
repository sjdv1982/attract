
.PHONY: default
default: all

# Find out the base directory
CURDIR = $(realpath $(PWD) )

OBJDIR = build
$(shell mkdir -p $(OBJDIR))

NAME = AttractServer
BINARY = lib$(NAME).so

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

SOURCES_CPP = $(shell find $(CURDIR)/../src/ -name "*.cpp")
SOURCES_CU = $(shell find $(CURDIR)/../src/ -name "*.cu")

#OBJECTS_CPP = $(SOURCES_CPP:.cpp=.o)
#OBJS_CPP = $(notdir $(OBJECTS_CPP))
OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))

#OBJECTS_CU = $(SOURCES_CU:.cu=.o)
OBJECTS_CU = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CU:.cu=.o)))

OBJECTS = $(OBJECTS_CPP) $(OBJECTS_CU) 

INCLUDES = -I../include

# AttractServer 
CXX = nvcc
ARCHFLAGS = -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50
ARCHFLAGS2 = -gencode arch=compute_30,code=compute_30 -gencode arch=compute_35,code=compute_35 -gencode arch=compute_50,code=compute_50
ifeq ($(TARGET), RELEASE)
	OFLAGS = -O3 -DNDEBUG -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED
else ifeq ($(TARGET), DEBUG) 
	OFLAGS = -O0 -G -g -Xcompiler -Wall -Xcompiler -Wextra
else
	OFLAGS =
endif	
CXXFLAGS =  $(OFLAGS) -Xcompiler -fPIC -std=c++11
 
LDFLAGS = --cudart static -shared --relocatable-device-code=true $(ARCHFLAGS) $(ARCHFLAGS2) -link
LIBS = -lnvToolsExt

# search directories for source files
VPATH = $(shell find ../src/ -type d)

$(BINARY): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC Linker'
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJECTS_CPP): $(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo 'Building file: "$<"'
	@echo 'Invoking: NVCC Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) --compile -x c++ -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '	
	
$(OBJECTS_CU): $(OBJDIR)/%.o: %.cu $(OBJDIR)/%.d
	@echo 'Building file: $<"'
	@echo 'Invoking: NVCC Compiler'
	$(CXX) $(CXXFLAGS) $(ARCHFLAGS) $(ARCHFLAGS2) $(INCLUDES) --compile --relocatable-device-code=true -x cu -o "$@" "$<"
	@echo 'Finished building: $<"'
	@echo ' '

depend: $(OBJECTS:.o=.d)
	@echo -n ''

$(OBJECTS_CPP:.o=.d): $(OBJDIR)/%.d: %.cpp
	@echo 'Building dependency file for: "$<"'
	@echo 'Invoking: NVCC Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -M -o "$(@:%.o=%.d)" "$<"
	@sed -i '1s,^,$(OBJDIR)/,' $(@:%.o=%.d)  
	@echo 'Finished building: "$<"'
	@echo ' '	
	
$(OBJECTS_CU:.o=.d): $(OBJDIR)/%.d: %.cu
	@echo 'Building dependency file for: $<"'
	@echo 'Invoking: NVCC Compiler'
	$(CXX) $(CXXFLAGS) $(ARCHFLAGS) $(INCLUDES) -M -o "$(@:%.o=%.d)" "$<"
	@sed -i '1s,^,$(OBJDIR)/,' $(@:%.o=%.d)  
	@echo 'Finished building: $<"'
	@echo ' '	

# do not include if we execute for cleanall or clean
ifneq (,$(filter cleanall clean,$(MAKECMDGOAL)))
	-include $(OBJECTS:.o=.d)
endif

	
.PHONY: all clean cleanall clean_deps

all: $(BINARY)

clean:
	rm -r $(OBJDIR)
	
cleanall: clean
	find . -type f -name "$(BINARY)" -delete

clean_deps:
	find $(OBJDIR) -type f -name "*.d" -delete

