%%% This file is to be included by latex in nemo.tex
%
% Import/Export: unofficial document

\myfile{system.inc}

{\it For sake of generality, the root directory of the package is designated
with the {\tt MIR} environment variable in this documentation, but may be
anything you like ({\tt GIPSY, NEMO}).}

\section{Tightly coupled machines: HOSTTYPE}

The classical solution to the problem of maintaining and expanding a
software package on a single computer by a group of programmers has been
approached in several ways.  In addition to commercially available
source control systems [e.g.  SCCS and RCS on UNIX, similar things
available under VMS, NOS(/BE)/(VE) etc.], there are various home
brewn schemes that can be cooked up depending on local needs and
resources.

This situation gets a bit more complicated when the package is also
to be maintained and expanded accoss a variety of loosely linked 
computers with a variety of operating systems. This is exactly the
situation the MIRIAD (henceforth referred to as MIR) package
finds itself under. 

Currently the package is actively being used, maintained and developed
on a Cray2 (UNICOS), SUN 3's and 4's (SUN OS 4.x), various VAX machines
(VMS 4.x, 5.x), a Multiflow (BSD 4.3), and to a lesser extent an Alliant
FX/8(?) and a Convex C-2.  {\it any comments Bobs?}.  Moreover, development
is now active on three campuses: Berkeley, Illinois and Maryland. 

We have decided to maintain a central repository, where the master copy
of all sources should reside, and cloning off is the preferred site from.
Currently this is CASTOR, the VMS engine at Illinois, although this will
become a bottleneck now that MIR is increasingly being installed on UNIX
operating systems. Every cluster of machines that can be seen as tightly
connected [i.e. directly accessible through an NFS type mechanism - 
perhaps 'rcp' should be seen as loosely connected, as {\tt ftp} is] can
be maintained and developed in a similar way as was done previously,
with a small modification.

The situation for a cluster of tightly connected machines, even with
different CPU's, has now successfully been used for over a year within
NEMO and more recently, MIR.  In this case all machines share a common
disk, using NFS or the like, which allows the source to be shared by all
machines.  Machine dependant object code lives in files in subdirectories
of directories where you normally would find such files. 

An example of a directory tree will hopefully clarify this 
[in unix notation]:

\begin{verbatim}
MIR/src/                The source tree:
        prog/             where programs live
        subs/             where subroutines live
    inc/                common include files
    adm/                Local administrativia 
        export/          export files to other sites
        import/          import files from other sites
    bin/                System dependant Executables
        sun3/            for sun3
        sun4/            for sun4
        mf/              for Multiflow
    lib/                System dependant Object Code (lib's, etc's)
        sun3/            for sun3
        sun4/            for sun4
        mf/              for mf
\end{verbatim}

One can see that, depending on which machine one works, the machine
dependant files live in separate subdirectories with a common name,
which we refer to as the {\tt HOSTYTPE} [The UNIX tcsh shell already has a
predefined shell variable {\tt HOSTTYPE}]. The process of 'seeing' the right
object files is made transparent to the user by the login script.

A user (programmers too) commonly adds a package to his environment
when he logs into the system by executing some kind of startup script. 
Under VMS these are often called {\tt LOGIN.COM}, under Unix
{\tt .cshrc},
{\tt .login} or {\tt .profile}, depending on which shell one uses.  In MIR we
call it {\tt MIRRC} (for historic reasons capitalized under UNIX).  For UNIX
one normally also adds the name of the directory where the binaries live
in ones search path by adding to the {\tt PATH} environment variable.  For VMS
life is a bit more tedious though. 

The Unix startup script for the csh-shell normally contains the
following lines (For Unix freaks: note the backquotes below):

\begin{verbatim}
    setenv HOSTTYPE `/usr/local/bin/hosttype`
    setenv MIR /user/mirth
    source $MIR/MIRRC
    ...
    set path=( ... $MIRBIN ... )
    ...    
\end{verbatim}

The first line, where the {\tt HOSTTYPE} is determined, is only needed if you
do not use the tcsh-shell, but will have to be supplied as a shell
script hosttype in /usr/local/bin.  If you would only use SUN
workstations, the SUN program {\tt arch} could be used 
instead of {\t hosttype}, {\it i.e.}:
\begin{verbatim}
    setenv HOSTTYPE `arch`
\end{verbatim}
The solution with a 'hosttype' program seems a bit more general,
especially since the newest shells have that notion too. 
Somewhere in the {\tt MIRRC} script there will be a line saying
\begin{verbatim}
    setenv MIRBIN $MIR/bin/$HOSTTYPE
\end{verbatim}
etc.

A last word on system dependant include files: we use a preprocessor for
Fortran-code (called ratty in MIR), and the normal 'cpp' preprocessor
for C-code.  The '-Dname' switch can be used to keep at least the source
code common for all machines by using a replacement script for the 'cc'
C compiler in the MIRBIN area.  This requires the UNIX system areas to
be after the MIRBIN area in the search path.  This solution has been
effectively used in NEMO, to hide compile-time differences between e.g. 
SUN3 and SUN4 machines. 

{\it
    The question of system dependant code has not been solved
    satisfactory, since we have subdirectories with files like
    src/subs/misc/unix/unicos/errmsg.c and
    src/subs/misc/vms/errmsg.c and
    src/subs/misc/errmsg.for

    Not very elegant for maintaince and installation purposes. 
    We shall have to come back to this at some later point, or merge
    code and use the preprocessor to isolate code.

    At some point it is perhaps unavoidable to have system dependant
    code, e.g. the MIR menu drivers for SUN and VMS VT100 screens
    live in separate unix/sun and vms diectories.
}

\section{Loosely coupled machines: import and export}

\subsection{Example (UNIX)}

\subsubsection{Export}

Although run by a nightly daemon, one can run it manually by issuing:
\begin{verbatim}
    % cd $MIR/adm
    % export.csh $MIR teuben bobs wakker wright
\end{verbatim}
It will send a tarfile with name 'DDMMYY.HOSTNAME' to all sites listed
in the {\tt \$MIR/adm/Netrc} file using ftp.

\subsubsection{Import}

Every morning the site manager gets mail what has been new in all the
sites of the project since their last update.  The new update files can
be found in {\tt \$MIR/adm/import} and are of the form {\t DDMMMYY.HOSTTYPE}.
Each
file is a tar file.  The site manager then runs a script 'import' to
merge the files within that tar file into the package tree:
\begin{verbatim}
    % cd $MIR/adm
    % import.csh import/30mar90.saturn $MIR "submit -m" 1
\end{verbatim}
The MIR program {\tt "submit -m"} can also be replaced by unix program 
{\tt mv} if
you don't mind overwriting old versions. 

For a site manger who has been away for a long time, and find zillions
of little tar files the following combination needs to be done:
\begin{verbatim}
    % cd $MIR/adm/import            # goto import
    % foreach file (*)              # 
      tar xf $file                  # untar each tar file ..
      rm $file                      # and remove them
      end
    % tar cf ../new.tar *           # combine all stuff into new tar file
    % cd ..                         # go back
    % import.csh new.tar $MIR "submit -m" 1
\end{verbatim}


\subsection{Implementation (UNIX)}

For a group a loosely coupled machines [i.e.  no shared disks, but e.g. 
programs like {\tt ftp} and {\tt rcp} can be used to directly communicate (COPY)
files] we are now using an experimental scheme of nightly daemons
which send updates around to a variety of sites.  For N sites which
independantly maintain AND develop code, there will then be (at most) N*(N-1)
update packages going around the network. 

This does however assume that any of two sites do not modify the SAME
code.  Other means of communication are needed to avoid this clash.  The
classical method of 'checking out source code' (cf.  RCS or SCCS) seems
a bit cumbersume, especially in our work environment, but is in
principle also possible.  We have listed a 'Person Responsible' in each
file, which will have to be contacted in case of major revisions.  There
are still a number of problems with this scheme, but experience will
learn if this is sufficient for our purposes. 

It is installed in the form of an import and export script. The
'export' script is normally run by an automatic daemon (we call it
'Nightly' which in turn calls export).  The 'import' script currently
has to be run manually, and will mergen an update package from another
site with the current version of the package. 

All activity connected with import and export occurs in a subdirectory
'adm', which should never be part of the import/export procedure itself. 

The export script probes a selected number of files/directories within
the package if they were modified since the last probe, and creates a tar
file out of them.  We will use a tar file, because it not only
simplifies maintaining the hierarchical file structure, but also retains
the time stamp on the files, when sending the update accross the
network. The name of the tar file if of the form {\tt DDMMMYY.HOSTNAME},
to uniquely identify their creation date and origin.

Files/Directories to be found in {\tt \$MIR/adm}:

\begin{verbatim}
    Nightly      shell script, normally run by cron, to do nightly work
    crontab      example crontab entry, see your local system how to use this
    export.csh   the actual export script, normally run by Nightly!
    import.csh   script to import an export file from another site
    Netrc        list of sites where to export to; used by 'export' script
    import/      Directory where imported tar files are stored
    export/      Directory where exported tar files are stored
\end{verbatim}

Export usage:
-------------

  "export root-name [mail-user(s)]"

e.g.

  % export.csh /saturn/mirth  teuben wakker bobs


Import usage:
-------------

  "import tar-file root-name import-program mode"

e.g.

  % import.csh nemo.tar /saturn/nemo "submit -m" 1

The 'export' script assumes the following:

  - in the root directory of the package there is a makefile, which if
    called as 'make new' produces a list of files which are new or
    modified since the last update. The time of the last update is found
    from the timestamp of a file 'Last\_ntar' in the 'adm' subdirectory.
  - There is a file 'Netrc' in the 'adm' subdirectory which contains per
    valid line 5 entries: the hostname (as provided by the system),
    internet name/number ftp needs, login name, password
    and subdirectory where the export tar files are to be placed.
    Example of such a line:

        saturn 128.8.251.43 mirth zonnig adm/import

    The 'ftp' program, together with a temporary '~/.netrc' will be
    used for automatic transfer of programs. If the implementor so
    wishes, this part of the script can be modified to use any other
    scheme to transport files to the other sites. A viable option could
    be rcp, if available, since it is much easier to use. 

  - It will create, within the export directory, a file of the name
    'DDMMMYY.HOSTNAME' if any files were found to be newer than the 
    check timefile 'Last\_ntar'. [later we may use a checkfile different
    for each site]

The 'import' script assumes the following:
  
  - all files to be checked are in a tar file somewhere

  - the principle directory structure in the tar file and the
    package are the same.

  - bla bla

\subsection{Implementation (VMS)}

[a more detailed descriptiption will follow a some later day]

\subsubsection{Example (Export)}

The current VMS implementation uses an export script named 'UPDATE.
It calls a number of other small files, some of them are site dependant:

\begin{verbatim}
    update.com              master script
    update_com.saturn       site dependant call, called by update
    hierarchy.saturn        translation list, needed by update_com.saturn
    last_update.saturn      timecheck file for it's last export
    update.lock             lockfile to prevent multiple runs

    rcp_update.com          copies all files from a list to a site
    
    note.com                ??
    oneoff.com              ??

    csh_update.com          builds 
    dcl_update.com  
\end{verbatim}

\subsubsection{Example (import)}

Keep track of what came in and run either 'CADD' or 'CLNK' on them to
add modules to library or add tasks to executable area respectively. 


\section{Improvements to be made to IMPORT/EXPORT}

- more symmetry between the UNIX and VMS version, to easy maintenance. 

- each export site can have their own update period. This will have
to go into the Netrc file. A different name for this file is also
more appropriate, since it ties in with the internet FTP program
right now. Perhaps and 'rcp' implementation will bring new light
The VMS implementation has a different check file for each
export site.

- save logfiles of what import is importing - this could perhaps
    be used by a script which then installs the modules into the 
    system. Within MIRTH the new subs files can be added to
    the mirth library in MIRLIB, whereas all prog files can be added as
    a new task to the MIRBIN area.


The whole timestamp business is clumsy, since a 'cp', as opposed to
'mv', will screw up the whole scheme.  A file which was updated at site
B will be sent to site A.  If site A then copies ('cp') the file into
the main directory, it would find it to be new the next time an export
was run, and would have been sent back to site B! Hence the files need
to be 'copied' with 'mv', i.e.  moved, from the exploded tar file to the
directory tree, of course using a scheme where a numbered backup file
will be created.  This could be done using 'submit' This process is
manual however, and hence 'export' will use the technique used in
'submit' to stick files back into the tree. 

To add a 'Nightly' script, one has to study the local implementation
of cron/crontabs. An example of a crontab entry can be found in the
'adm' subdirectory of your package.

makefiles screw up in the scheme where the same file is in two
dir's, e.g. to hide system dependance. When one is modified,
and we execute make from the other, we get screwed !!! Long live
manual labor.

