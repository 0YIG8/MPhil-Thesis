# Credits: Dr James Fergusson (Advanced Research Computing lecture slides)

BUILDS = compile build clean clean-VTK

.PHONY := $(BUILDS)

.DEFAULT_GOAL := compile

compile :
	@echo "Compiling Project"
	cmake --build build
	@echo "Compilation complete!"
	@echo "Run ./bin/main_Cavity cavity.cfg OR ./bin/main_MRF MRF.cfg to run the simulation."

build :
	@echo "Building CMake Project"
	cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -S . -B build
	@echo "CMake build complete!"

clean :
	@echo "Cleaning CMake Project"
	rm -rf build bin
	@echo "CMake clean complete!"
	@echo "Run 'make build' to build project."

clean-VTK :
	@if [ "$(test_case)" = "" ]; then \
		echo "Please specify a test case output folder to clean (i.e., make clean-VTK test_case=mixerSRF)"; \
	elif [ -d "$(test_case)/VTK" ]; then \
		echo "Cleaning $(test_case) output folder!"; \
		rm $(test_case)/VTK/*.vt*; \
		echo "$(test_case) VTK output folder cleared!"; \
	else \
		echo "Invalid output folder!"; \
	fi