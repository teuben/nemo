
c
c standard colour indices
c
	integer BLACK, RED, GREEN, YELLOW, BLUE, MAGENT, 
     +		CYAN, WHITE
	parameter(BLACK = 0, RED = 1, GREEN = 2, YELLOW = 3)
	parameter(BLUE = 4, MAGENT = 5, CYAN = 6, WHITE = 7)

c
c polygon modes
c
	integer PYM_PO, PYM_LI, PYM_FI, PYM_HO
	parameter(PYM_PO = 0, PYM_LI = 0, PYM_FI = 1)
	parameter(PYM_HO = 1)
c
c Misc defines
c
	integer FLAT, SMOOTH, GDXPMA, GDYPMA
	parameter(FLAT = 0, SMOOTH = 1, GDXPMA = 1, GDYPMA = 2)
