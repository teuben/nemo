	COMMENT	$

		Handles Double Buffering on the CGA.
	$

		extrn	__buffer_segment:word
		extrn	__buffer_offset:word

		public	_cga_swapbuf, _cga_frontbuf

CGAA_TEXT	SEGMENT	byte	public	'CODE'
 		ASSUME	cs:CGAA_TEXT

	COMMENT	$
		Copies the RAM buffer to display buffer.
	$

_cga_swapbuf	proc	far
	push		bp
	mov		bp, sp
	sub		sp, 8

	push		es
	push		ds
	push		si
	push		di

	mov		ax, 0b800h	
	mov		es, ax		; ES:DI -> video buffer
	xor		di, di

	mov		ax, __buffer_segment
	mov		ds, ax
	mov		si, __buffer_offset	; DS:SI -> buffer in RAM

	;jmp		bonk


	; Wait for verticle blanking interval
	mov		dx, 3DAh	; CGA status port
L01:
	mov		cx, 22		; HRT timeout val
L02:	
	in		al, dx
	test		al, 8
	jnz		L02
	test		al, 1
	jnz		L02

	cli
L03:
	in		al, dx
	test		al, 1
	loopnz		L03

	sti

	jz		L01
	
	; Blank the display
	mov		dl, 0D8h
	mov		al, ((1 SHL 4) OR (1 SHL 5))
	out		dx, al

	; Move the data
	mov		cx, 2000h		; Size of buffer.
bonk:
	cld
	rep		movsw 

	; Reenable the Display

	or		al, (1 SHL 5)
	out		dx, al

	pop		di
	pop		si
	pop		ds
	pop		es
	
	mov		sp, bp
	pop		bp
	ret
_cga_swapbuf		endp


	COMMENT	$
		Write into the screen buffer.
	$

_cga_frontbuf	proc	far
	push		bp
	mov		bp, sp
	
	mov		__buffer_segment, 0b800h
	mov		__buffer_offset, 0

	mov		sp, bp
	pop		bp
	ret
_cga_frontbuf		endp

CGAA_TEXT		ends

			end
