CXXFLAGS := $(CXXFLAGS) -std=c++11  -O3 -Wall 
#it allows to put object files in a subfolder
obj=$(patsubst %.o,$(OBJPATH)%.o, random.o main.o Ising1D.o)

		#source path
VPATH=./source/
		#object files path
OBJPATH=./obj/
		#executable name
EXENAME=Ising1D

all: $(obj)
ifeq ($(wildcard outputs*),) #search for outputs path, create it if there isn't   
    @mkdir -p outputs                                                           
endif  
	$(CXX) $(OBJPATH)*.o -o $(EXENAME) $(CXXFLAGS)

$(OBJPATH)%.o: %.cpp
ifeq ($(wildcard $(OBJPATH)*),) #search for obj path, create it if there isn't
	@mkdir -p $(OBJPATH)
endif
	$(CXX) -c $< -o $@ $(CXXFLAGS) 

.PHONY: clean
clean:
	rm $(OBJPATH)*.o $(EXENAME)
