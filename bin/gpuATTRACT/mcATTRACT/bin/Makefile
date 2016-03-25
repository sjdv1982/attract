
.PHONY: default
default: all

# Find out the base directory
CURDIR = $(realpath $(PWD) )

OBJDIR = build
$(shell mkdir -p $(OBJDIR))

NAME = mcATTRACT
BINARY = $(NAME)



# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE

SOURCES_CPP = $(shell find $(CURDIR)/../src/ -name "*.cpp")

#OBJECTS_CPP = $(SOURCES_CPP:.cpp=.o)
#OBJS_CPP = $(notdir $(OBJECTS_CPP))
OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))

#OBJECTS_CU = $(SOURCES_CU:.cu=.o)

OBJECTS = $(OBJECTS_CPP)

ASDIR = ../../AttractServer
INCLUDES = -I../src -I$(CUDADIR)/include -I$(ASDIR)/include

# AttractServer 
CXX = g++
ifeq ($(TARGET), RELEASE)
	OFLAGS = -O3 -DNDEBUG
else ifeq ($(TARGET), DEBUG) 
	OFLAGS = -O0 -g -Wall -Wextra
else
	OFLAGS =
endif	
CXXFLAGS =  $(OFLAGS) -std=c++11 -fmessage-length=0
 
LDFLAGS = -L$(ASDIR)/lib -L$(CUDADIR)/lib64 
LIBS = -lnvToolsExt -lcudart -lAttractServer

# search directories for source files
VPATH = $(shell find ../src/ -type d)

# pull in dependency info for *existing* .o files


$(BINARY): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJECTS): $(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '	

$(OBJDIR)/%.d: ;	
	

-include $(OBJECTS:.o=.d)
	
.PHONY: all clean cleanall clean_deps

all: $(BINARY)

clean:
	rm -r $(OBJDIR)
	
cleanall: clean
	find . -type f -name "$(BINARY)" -delete

clean_deps:
	find $(OBJDIR) -type f -name "*.d" -delete



#EM_LDFLAGS = -L"." -L"${CUDADIR}/lib64"
#MC_LDFLAGS = -L"." -L"${CUDADIR}/lib64"
