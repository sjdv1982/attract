
.PHONY: default
default: all

# Find out the base directory
CURDIR = $(realpath $(PWD) )

OBJDIR = build
$(shell mkdir -p $(OBJDIR))

NAME = emATTRACT
BINARY = $(NAME)



# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

SOURCES_CPP = $(shell find $(CURDIR)/../src/ -name "*.cpp")
OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))

SOURCES_F = $(shell find $(CURDIR)/../src/ -name "*.f")
OBJECTS_F = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_F:.f=.o)))

OBJECTS = $(OBJECTS_CPP) $(OBJECTS_F)

ASDIR = ../../AttractServer
INCLUDES = -I../src -I$(CUDADIR)/include -I$(ASDIR)/include

# AttractServer 
CXX = g++
ifeq ($(TARGET), RELEASE)
	CXX_OFLAGS = -O3 -DNDEBUG
	FXX_OFLAGS = -O2 
else ifeq ($(TARGET), DEBUG) 
	CXX_OFLAGS = -O0 -g -Wall -Wextra
	FXX_OFLAGS = -O0 -g -Wall -Wextra
else
	CXX_OFLAGS =
	FXX_OFLAGS = 
endif	
CXXFLAGS =  $(CXX_OFLAGS) -std=c++11 -fmessage-length=0
 
LDFLAGS = -L$(ASDIR)/lib -L$(CUDADIR)/lib64
LIBS = -lboost_coroutine -lnvToolsExt -lcudart -lAttractServer -lboost_context -lgfortran

FXX = gfortran
FXXFLAGS = $(FXX_OFLAGS) -fcray-pointer

# search directories for source files
VPATH = $(shell find ../src/ -type d)

# pull in dependency info for *existing* .o files


$(BINARY): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJECTS_CPP): $(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '	

$(OBJDIR)/%.d: ;

$(OBJECTS_F): $(OBJDIR)/%.o: %.f
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC Fortran Compiler'
	$(FXX) $(FXXFLAGS) -c -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '		
	

-include $(OBJECTS:.o=.d)
	
.PHONY: all clean cleanall clean_deps

all: $(BINARY)

clean:
	rm -r $(OBJDIR)
	
cleanall: clean
	#rm $(BINARY)
	find . -type f -name "$(BINARY)" -delete

clean_deps:
	find $(OBJDIR) -type f -name "*.d" -delete



#EM_LDFLAGS = -L"." -L"${CUDADIR}/lib64"
#MC_LDFLAGS = -L"." -L"${CUDADIR}/lib64"