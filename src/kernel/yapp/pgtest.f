        program curtest
c
        character ch*1
        integer pgcurse
c
        call pgask(.FALSE.)
        call pgbegin(0,'/ps',1,1)
        call pgenv(0.,10.,0.,20.,0,1)
        call pglabel("x","y","Cursor test")
c
c
        k = pgcurse(x,y,ch)
        write(*,*) 'pgcurse: k,x,y,ch=',k,x,y,ch
c
        call pgend
        end
