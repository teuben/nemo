	interface to subroutine C_bgnpoint[c]
	end

	interface to subroutine C_bgnline[c]
	end

	interface to subroutine C_bgnclosedline[c]
	end

	interface to subroutine C_bgnpolygon[c]
	end

	interface to subroutine C_endpoint[c]
	end

	interface to subroutine C_endline[c]
	end

	interface to subroutine C_endclosedline[c]
	end

	interface to subroutine C_endpolygon[c]
	end

	interface to subroutine C_endtmesh[c]
	end

	interface to subroutine C_endqstrip[c]
	end

	interface to subroutine C_bgntmesh[c]
	end

	interface to subroutine C_bgnqstrip[c]
	end

	interface to subroutine C_swaptmesh[c]
	end

	subroutine bgntmesh
	call C_bgntmesh
	end

	subroutine bgntme
	call C_bgntmesh
	end

	subroutine endtmesh
	call C_endtmesh
	end

	subroutine endtme
	call C_endtmesh
	end

	subroutine swaptmesh
	call C_swaptmesh
	end

	subroutine swaptm
	call C_swaptmesh
	end

	subroutine bgnqstrip
	call C_bgnqstrip
	end

	subroutine bgnqst
	call C_bgnqstrip
	end

	subroutine endqstrip
	call C_endqstrip
	end

	subroutine endqst
	call C_endqstrip
	end

	subroutine bgnpoint
	call C_bgnpoint
	end

	subroutine bgnline
	call C_bgnline
	end

	subroutine bgnclosedline
	call C_bgnclosedline
	end

	subroutine bgnpolygon
	call C_bgnpolygon
	end

	subroutine endpoint
	call C_endpoint
	end

	subroutine endline
	call C_endline
	end

	subroutine endclosedline
	call C_endclosedline
	end

	subroutine endpolygon
	call C_endpolygon
	end
	subroutine bgnpoi
	call C_bgnpoint
	end

	subroutine bgnlin
	call C_bgnline
	end

	subroutine bgnclo
	call C_bgnclosedline
	end

	subroutine bgnpol
	call C_bgnpolygon
	end

	subroutine endpoi
	call C_endpoint
	end

	subroutine endlin
	call C_endline
	end

	subroutine endclo
	call C_endclosedline
	end

	subroutine endpol
	call C_endpolygon
	end
