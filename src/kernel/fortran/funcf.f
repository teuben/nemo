        DOUBLE PRECISION FUNCTION FUNCF(X)
        DOUBLE PRECISION X, FUNCC
c	WRITE (*,*) 'FuncF', X
        IF (X.LT.1.0d0) THEN
            FUNCF = SQRT(X)
        ELSE
            FUNCF = FUNCC(1.0d0/X)
        ENDIF 
        END
