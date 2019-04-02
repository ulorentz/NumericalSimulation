#for multithreaded compiling                                                       
OS := $(shell uname)
ifeq ($(OS), Linux)                    
    export MAKEFLAGS="-j $(grep -c ^processor /proc/cpuinfo)" \
else($(OS), Darwin)                                    
    export MAKEFLAGS="-j $(sysctl-n hw.ncpu)" \
else                                                        
    #Windows?  not worth...
endif 

TOPTARGETS := all clean

SUBDIRS := $(wildcard */.)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)
