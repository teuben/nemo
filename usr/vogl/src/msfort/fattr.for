	interface to subroutine C_pushattributes[c]
	end

	interface to subroutine C_popattributes[c]
	end

	subroutine pushattributes
	call C_pushattributes
	end

	subroutine popattributes
	call C_popattributes
	end

	subroutine pushat
	call C_pushattributes
	end

	subroutine popatt
	call C_popattributes
	end
