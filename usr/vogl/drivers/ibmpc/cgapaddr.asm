
	TITLE	CGAPADDR - Return video buffer address of a pixel.
	NAME	CGAPADDR
	PAGE	55,132


	COMMENT	$

	Name:		CGAPADDR

	Function:	Determine buffer address of pixel in 640x200 2 color mode

	Caller:			AX = y co-ord (0-199)
				BX = x co-ord (0-639)

	Returns:
				AH = bit mask
				BX = byte offset in buffer
				CL = number of bits to shift left
				ES = video buffer segment

		$

		extrn	__buffer_segment:word
		extrn	__buffer_offset:word

		public	__cga_set_buffer
		public	cgapaddr

CGAA_TEXT	SEGMENT	byte public 'CODE'
		ASSUME  CS: CGAA_TEXT

__cga_set_buffer	proc	far

	push	bp
	mov	bp, sp

	les	di, [bp + 6]
	mov	__buffer_segment, es
	mov	__buffer_offset, di

	mov	sp, bp
	pop	bp
	ret
__cga_set_buffer	endp

cgapaddr	PROC	far

	mov		cl,bl		; CL = low-order byte of x
	xchg		ah,al		; AX = 100h * y
	shr		bx,1		; BX = x/2
	shr		ax,1		; AL = 80h*(y&1)
	add		bh,al		; BX = x/2 + 8000h*(y&1)
	xor		al,al		; AX = 100h*(y/2)
	add		bx,ax		; BX = x/2 + 8000h*(y&1) + 100h*(y/2)
	shr		ax,1
	shr		ax,1		; AX = 40h*(y/2)
	add		bx,ax		; BX = x/2 + 8000h*(y&1) + 140h*(y/2)
	shr		bx,1
	shr		bx,1		; BX = x/8 + 2000h*(y&1) + 50h*(y/2)

	add		bx, __buffer_offset

	mov		ax, __buffer_segment
	mov		es,ax		; ES:BX = byte address of pixel

	and		cl,7		; CL = x &  
	xor		cl,7		; CL = number of bits to shift left
	mov		ah,1		; AH unshifted bit mask
	ret

cgapaddr	endp

CGAA_TEXT	ends

		end
