%% This file is to be included by latex in nemo.tex
%
% chapter: Directory Structure
\myfile{dirs.inc}
\mylabel{a:dirs}

In this Appendix we show the directory tree structure as  you should
find it under NEMO. 
Directory names preceded with a (*) are being fully exported by the
export procedures. 
The ``{\tt src}'' tree contains all of the official 
source code, and is expanded a bit on following pages. 
The ``{\tt usr}'' (or ``{\tt contrib}'') tree is what used 
to be the {\tt \$NEMO/src} tree in the first
releases of NEMO (Version 1.x). 
Some of the code in {\tt \$NEMO/usr} is now
outdated, and superseded by more recent versions
somewhere in {\tt \$NEMO/src} tree.
Any 'personal' (most public domain also) code 
should now be placed in the appropriate {\tt \$NEMO/usr} tree before
it can be accepted (moved over)
into the {\tt \$NEMO/src} tree. We advice programmers
to stick as closely as possible to the install procedures
in the {\tt \$NEMO/src} tree, if at all possible.
See the Makefiles and the discussion
on these in the Programmers Guide (see \ref{ss:makefile}).


\index{nemo, directory structure}
\newpage
\begin{table}

\footnotesize\begin{verbatim}
    adm/                        local administrativia
        export/                 - exported tar files
        import/                 - imported tar files from other sites
    bin/                        executables ($NEMOBIN)
  * bugs/                       ... the known ones (optional)
  * csh/                        example shell scripts
  * data/                       NEMO data files
  * demo/                       ... demo scripts
    etc/                        various /etc type files
  * inc/                        standard include dir for NEMO
        fortran/                -- fortran I/O routines
        multicode/              -- josh' multicode
        snapshot/               -- snapshot
    lib/                        object libraries ($NEMOLIB)
    local/                      things local to your site
  * man/                        manual pages
        man1/                   -- programs
        man3/                   -- library routines
        man5/                   -- file structure
        man8/                   -- miscellaneous/maintenance
        doc/                    -- inline .doc for miriad/nemotool
  * news/                       news on new/modified software
    obj/                        dynamic object loader files ($NEMOOBJ)
        bodytrans/              
        potential/
  * src/                        The whole source tree: (see also next page)
        kernel/                  Basic general stuff (i/o, utils)
        nbody/                   Nbody (snapshot)
        image/                   Images
        orbit/                   Orbit calulations, potentials
        tools/                   Various borrowed tools
        scripts/                 Various (maintenance) scripts
        tutorial/                Tutorial material - toy programs
        ...                     ==See next page for expanded view==
 (*)usr/                        Public donation programs (old NEMO/src)
        josh/			
        pjt/
        piet/
        ...
  * text/                       alternate doc, mostly (la)tex stuff
        manuals/                latex documents (e.g. this manual)
    tmp/                        temporary storage
\end{verbatim}\normalsize

\caption{Overview of the {\tt \$NEMO/} tree}
\end{table}

\newpage

\begin{table}
\footnotesize\begin{verbatim}
  * src/                  THE SOURCE TREE
        kernel/            All of these needed for full NEMO implementation
             cores/         independent but otherwise essential utilities for IO
             io/            getparam & filestruct (user interface, database i/o)
             misc/          miscellaneous library utilities
             yapp/          YAPP plotting device drivers
             loadobj/       Dynamic object loader
             tab/           Table manipulation stuff
             fortran/       Fortran to C interface
        nbody/             Nbody
             cores/         essential utilities for Nbody
             io/            I/O get_snap etc.
             init/          create models
             evolve/        Evolve (integrate) models
                hackcode/    Treecode
                   export/    --export version - no bells and whistles
                   source1/   Barnes and Hut (official c version)
                   treecode/  Lars Hernquist fortran version of treecode
                multicode/   Simon White`s multipole expansion (Barnes version)
                aarseth/     Sverre Aarseth
                   nbody0/
                   nbody1/
                potcode/     Non-selfconsistent particle integrator
                silly/       Silly demo version of N^2 code
             trans/         Transformation utilities for snapshots
             reduc/         Analysis utilities
        image/             2- and 3D image analysis
             cores/
             io/            
             fits/
             misc/
        orbit/             Orbit analysis
             io/
             cores/
             potential/      Potentials
                data/         repository of potential source code files
        tools/             Tools, independent and mostly (PD) borroware
             movie/          movie for sun raster files (PD)
             movietool/      movietool for sun raster files (PD)
             ds/             image display program (Sebok)
        tutor/             Example programs and Makefile to toy with
\end{verbatim}\normalsize

\caption{Some pieces of the {\tt \$NEMO/src/} tree (this branch is officially
exported)}
\end{table}

\newpage

\begin{table}
\footnotesize\begin{verbatim}
 (*)usr/                        User contributed software tree
        aarseth/                --Sverre Aarseth 
        josh/                   --Josh Barnes 
            clib/
            hackcode/
                export/
                source1/
            multicode/
            nbody/
                bodytrans/
                snapshot/
                tools/
        lars/                   --Lars Hernquist 
            treecode/
        makino/                 --Jun Makino
            nbody/
        mbellon/                --Mark Bellon
            gravsim/            
        nemo/                   --NEMO
            demo/
            maint/
            news/
            util/
                shar/
                f2c/
            yapp/
        piet/                   --Piet Hut
            clib/
            newton/
        pjt/                    --Peter Teuben
            clib/
            image/
            nbody/
            orbit/
        pswisnov/               --Peter Wisnovksi 
            ccddisplay/
            snapsmooth/
        unemo/                  -- Micro Nemo (old version)
            lib/
            nbody/
\end{verbatim}\normalsize
\caption{Some pieces of the {\tt \$NEMO/usr/} tree (this branch is officially
not exported)}
\end{table}

\index{Hernquist L.}\index{Barnes J.}\index{Bellon M.}\index{Aarseth S.}
\index{Teuben P.} \index{Wisnovski P.} \index{Hut P.}
\index{Newton I.} \index{Makino J.}
