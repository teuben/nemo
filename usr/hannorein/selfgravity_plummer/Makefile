export OPENGL=0# Set this to 1 to enable OpenGL
export SERVER=1# Set this to 1 to enable the visualization web server
include ../rebound/src/Makefile.defs

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
#CCPROBLEM=$(CC) -I../rebound/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) $< -L$(NEMOLIB) -lrebound $(LIB) -o rebound_plummer
CCPROBLEM=$(CC) -I$(NEMOINC) -I$(NEMOLIB) -I../rebound/src/ $(OPT) $(PREDEF) $< -L$(NEMOLIB) -lnemo -lrebound $(LIB) -o rebound_plummer

.PHONY: all librebound
all: rebound_plummer.c
	@echo "Compiling $< ..."
	$(CCPROBLEM)
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

clean:
	@echo No clean here yet

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_console.html -o rebound_console.html
