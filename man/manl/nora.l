.TH NORA 1 "7 June 1994"
.SH NAME
nora -\ analysis tool for N-body MODELS 
.SH SYNTAX
.B nora
.SH DESCRIPTION
.PP
.I nora
reads lines of text from the standard input 
and uses any standard
.IR pgplot (3)
device to display and analyze Nbody models.
.SH COMMANDS
.nf
                !               macro
                ?               menu            quit
                show            list            device
                data            getmodel
                datar8          getmodelr8
                cmcent          medcent         potcent
                setview         levels
                bodsrange       bodsreloc
                space           addspace
                vel             addvel          velRT
                dens            sbr             Rm
                innertia
                rotc            addrotc
                vdist           addvdist
                prof4           multiprof4
                contour         gray            points
                outFITS         outBDF
                edhead          outmodel
 
.fi
.SH AUTHOR
Marc Balcells
.SH "SEE ALSO"
pgplot (3)
