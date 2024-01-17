export OPENGL=0# Set this to 1 to enable OpenGL
export SERVER=1# Set this to 1 to enable the visualization web server
include ../rebound/src/Makefile.defs

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
CCPROBLEM0=$(CC) -I$(NEMOINC) -I$(NEMOLIB) $(OPT) $(PREDEF) $< -L$(NEMOLIB) -lnemo -lrebound $(LIB)

.PHONY: all librebound

all:   reb_viewer reb2s

reb_viewer: reb_viewer.c
	@echo "Compiling $< ..."
	$(CCPROBLEM0) -o reb_viewer

reb_integrate: reb_integrate.c
	@echo "Compiling $< ..."
	$(CCPROBLEM0) -o reb_integrate


reb2s: reb2s.c
	@echo "Compiling $< ..."
	$(CCPROBLEM0) -o reb2s

s2reb: s2reb.c
	@echo "Compiling $< ..."
	$(CCPROBLEM0) -o s2reb

clean:
	@echo No clean here yet

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_console.html -o rebound_console.html
