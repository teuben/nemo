		TITLE	EGAPADDR - Return video buffer address of a pixel.
		NAME	EGAPADDR
		PAGE	55,132


		COMMENT	$

	Name:		EGAPADDR

	Function:	Determine buffer address of pixel in native EGA and VGA modes:
			320 x 200 16 colour
			640 x 200 16 colour
			640 x 350 16 colour
			640 x 350 monochrome (4 colour)
			640 x 480 2 colour
			640 x 480 16 colour

	Caller:		AX = y co-ord
			BX = x co-ord

	Returns:
			AH = bit mask
			CL = number of bits to shift left
			BX = byte offset in buffer
			ES = video buffer segment

								$

BytesPerLine	EQU		80

		extrn	__buffer_segment:word

VEGA_TEXT	SEGMENT	byte public 'CODE'
		ASSUME	cs:VEGA_TEXT

		PUBLIC	egapaddr

egapaddr	PROC	far

	mov		cl,bl
	push		dx

	mov		dx,BytesPerLine
	mul		dx

	pop		dx
	shr		bx,1
	shr		bx,1
	shr		bx,1
	add		bx,ax

	mov		ax,__buffer_segment
	mov		es,ax

	and		cl,7
	xor		cl,7
	mov		ah,1
	ret

egapaddr	endp

VEGA_TEXT	ends

		end


