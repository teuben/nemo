	TITLE	HGCLINE - Fast line drawing routine.
	NAME	HGCLINE
	PAGE	55,132


	COMMENT	$

	Name:		HGCLINE

	Function:	Draw a line in 720x348 Hercules graphics mode

	Caller:		Microsoft C:

			void 	hgcline(x1, y1, x2, y2, n)

			int	x1, y1, x2, y2;		/* pixel co-ords */
			int	n;			/* pixel value */


		$


; Stack frame addressing - LARGE CODE MODEL

ARGx1		EQU	word ptr [bp+6]
ARGy1		EQU	word ptr [bp+8]
ARGx2		EQU	word ptr [bp+10]
ARGy2		EQU	word ptr [bp+12]
ARGn		EQU	byte ptr [bp+14]

VARleafincr	EQU	word ptr [bp-6]
VARincr1	EQU	word ptr [bp-8]
VARincr2	EQU	word ptr [bp-10]
VARroutine	EQU	word ptr [bp-12]

ByteOffsetShift	EQU	3	; Used to convert pixels to byte offset

DGROUP		GROUP	HERC_DATA

HERC_TEXT	SEGMENT	byte public 'CODE'
		ASSUME	cs:HERC_TEXT,ds:DGROUP

		EXTRN	hgcpaddr:far
		PUBLIC	_hgcline

_hgcline	PROC	far

		push	bp	; Set up stack frame
		mov	bp,sp
		sub	sp,12	; Stack space for local variables
		push	si
		push	di

		mov	si,2000h	; Increment for video buffer interleave
		mov	di,90-8000h	; Increment for last to first interleave

		mov	cx,ARGx2
		sub	cx,ARGx1	; CX = x2 - x1
		jz	VertLine	; jump if vertical line

; force x1 < x2

		jns	L01		; jump if x2 > x1

		neg	cx		; CX = x1 - x2

		mov	bx,ARGx2	; exchange x1 and x2
		xchg	bx,ARGx1
		mov	ARGx2,bx

		mov	bx,ARGy2	; exchange y1 and y2
		xchg	bx,ARGy1
		mov	ARGy2,bx

; calc dy = abs(y2 - y1)

L01:		mov	bx,ARGy2
		sub	bx,ARGy1	; BX = y2 - y1
;;;;;;;;;;;;;;	jz	HorizLine
		jnz	L02
		jmp	HorizLine	; jump if Horizontal
L02:
		jns	L03

		neg	bx		; BX = y1 - y2
		neg	si		; negate increment for buffer interleave
		neg	di

; select appropriate routine for slope of line

L03:		mov	VARleafincr,di	; save increment for buffer interleave
		mov	VARroutine,offset LoSlopeLine
		cmp	bx,cx
		jle	L04		; jump if dy <= dx (slopr <= 1)
		mov	VARroutine,offset HiSlopeLine
		xchg	bx,cx		; exchange dy and dx

; calc initial decision variable and increments

L04:		shl	bx,1		; BX = 2 * dy
		mov	VARincr1,bx	; incr1 = 2 * dy
		sub	bx,cx
		mov	di,bx		; DI = d = 2*dy - dx
		sub	bx,cx
		mov	VARincr2,bx	; incr2 = 2*(dy - dx)

; calc first pixel address

		push	cx		; preserve the register
		mov	ax,ARGy1	; ax = y
		mov	bx,ARGx1	; bx = x
		call	hgcpaddr	; ah = bit mask
					; es:bx -> buffer
					; CL = # bits to shift left

		mov	al,ARGn		; AL = unshifted pixel value
		shl	ax,cl		; AH = bitmask in proper position
					; AL = pixel value in proper position

		mov	dx,ax		; DH = bitmask
					; DL = pixel value
		not	dh		; DH = inverse bitmask

		pop	cx
		inc	cx		; CX = # of pixels to draw
					
		jmp	VARroutine	; Jump to the appropiate routine for slope

; routine for verticle lines

VertLine:	mov	ax,ARGy1	; AX = y1
		mov	bx,ARGy2	; BX = y2
		mov	cx,bx		; 
		sub	cx,ax		; CX = dy
		jge	L31		; jump if dy >= 0

		neg	cx		; Force dy >= 0
		mov	ax,bx		; AX = y2

L31:		inc	cx		; CX = # of pixels to draw
		mov	bx,ARGx1	; BX = x
		push	cx
		call	hgcpaddr	; AH = bit mask
					; ES:BX -> video buffer
					; CL = # bits to shift left
		mov	al,ARGn		; AL = pixel value
		shl	ax,cl		; AH = bit mask in proper position
					; AL =pixel value in proper position
		not	ah		; AH = inverse bit mask
		pop	cx

; draw the line

		test	al,al
		jz	L34		; jump if pixel value is zero

L32:		or	es:[bx],al	; set pixel values in buffer

		add	bx,si		; increment to next portion of interleave
		jns	L33
		add	bx,di		; increment to first portion of interleave
L33:		loop	L32
		jmp	short L36

L34:		and	es:[bx],ah	; reset pixel values in buffer
		add	bx,si		; incremant to next portion of interleave
		jns	L35
		add	bx,di		; increment to first portion of interleave
L35:		loop	L34

L36:		jmp	Lexit

; routine for horizontal line (slope=0)

HorizLine:	mov	ax,ARGy1
		mov	bx,ARGx1
		call	hgcpaddr	; AH = bit mask
					; ES:BX -> video buffer
					; CL = # bits to shift left
		mov	di,bx		; ES:DI -> buffer
		mov	dh,ah
		not	dh		; DH = unshifted bit mask for left byte

		mov	dl,0FFh		; DL = unshifted bit mask for right byte

		shl	dh,cl		; DH = reverse bit mask for first byte
		not	dh		; DH = bit mask for first byte

		mov	cx,ARGx2
		and	cl,7
		xor	cl,7		; CL = # bits to shift left
		shl	dl,cl		; DL = bit mask for last byte

; determine byte offset of first and last pixel in line

		mov	ax,ARGx2	; AX = x2
		mov	bx,ARGx1	; BX = x1

		mov	cl,ByteOffsetShift	; # of bits to shift to convert
						; pixels to bytes

		shr	ax,cl		; AX = byte offset of x2
		shr	bx,cl		; BX = byte offset of x1
		mov	cx,ax
		sub	cx,bx		; CX = (# bytes in line) - 1

; propagate pixel value throughout one byte

		mov	bx,offset DGROUP:PropagatedPixel
		mov	al,ARGn		; AL = pixel value
		xlat			; AL = propagated pixel value

; set pixels in leftmost byte of line

		or	dh,dh
		js	L43		; jump if byte-aligned (x1 is leftmost
					; pixel in byte)
		or	cx,cx
		jnz	L42		; jump if more than one byte in line

		and	dl,dh		; bit mask for the line
		jmp	short L44

L42:		mov	ah,al
		and	ah,dh		; AH = masked pixel bits
		not	dh		; DH = reverse bit mask for 1st byte
		and	es:[di],dh	; zero masked pixels in buffer
		or	es:[di],ah	; update masked pixels in buffer
		inc	di
		dec	cx

; draw remainder of the line

L43:		rep	stosb		; update all pixels in the line

; set pixels in rightmost byte of line

L44:		and	al,dl		; AL = masked pixels for last byte
		not	dl
		and	es:[di],dl	; zero masked pixels in buffer
		or	es:[di],al	; update masked pixels in buffer

		jmp	Lexit


; routine for dy <= dx (slope <= 1)	; ES:BX -> video buffer
					; CX = # pixels to draw
					; DH = inverse bit mask
					; DL = pixel value in proper position
					; SI = buffer interleave increment
					; DI = decision variable

LoSlopeLine:

L10:		mov	ah,es:[bx]	; AH = byte from video buffer

L11:		and	ah,dh		; zero pixel value at current bit offset
		or	ah,dl		; set pixel value in byte

		ror	dl,1		; rotate pixel value
		ror	dh,1		; rotate bit mask
		jnc	L14		; jump if bit make rotated to
					; leftmost pixel position

; bit mask not shifted out

		or	di,di		; test sign of d
		jns	L12		; jump if d >= 0

		add	di,VARincr1	; d = d + incr1
		loop	L11

		mov	es:[bx],ah	; store remaining pixels in buffer
		jmp	short Lexit

L12:		add	di,VARincr2	; d = d + incr2
		mov	es:[bx],ah	; update buffer

		add	bx,si		; increment y
		jns	L13		; jump if not in last interleave

		add	bx,VARleafincr	; increment into next interleave

L13:		loop	L10
		jmp	short Lexit

; bit mask shifted out

L14:		mov	es:[bx],ah	; update buffer
		inc	bx		; BX = offset of next byte

		or	di,di		; test sign of d
		jns	L15		; jump if non-negative

		add	di,VARincr1	; d = d + incr1
		loop	L10
		jmp	short Lexit

L15:		add	di,VARincr2	; d = d + incr2

		add	bx,si		; increment y
		jns	L16		; jump if not in last interleave

		add	bx,VARleafincr	; increment into next interleave
		
L16:		loop	L10		; loop until all pixels are set
		jmp	short Lexit


; routine for dy > dx (slope > 1)	; ES:BX -> video buffer
					; CX = # pixels to draw
					; DH = inverse bit mask
					; DL = pixel value in proper position
					; SI = buffer interleave increment
					; DI = decision variable

HiSlopeLine:

L21:		and	es:[bx],dh	; zero pixel value in video buffer
		or	es:[bx],dl	; set pixel value in byte

		add	bx,si		; increment y
		jns	L22		; jump if not in last interleave

		add	bx,VARleafincr	; increment into next interleave

L22:		or	di,di
		jns	L23		; jump if d >= 0

		add	di,VARincr1	; d = d + incr1
		loop	L21
		jmp	short Lexit

L23:		add	di,VARincr2	; d = d + incr2
		ror	dl,1		; rotate pixel value
		ror	dh,1		; rotate bit mask
		cmc			; cf set if bit mask not rotated to
					; leftmost pixel position
		adc	bx,0		; BX = offset of next byte

		loop	L21


Lexit:
		pop	di		; restore registers and return
		pop	si
		mov	sp,bp
		pop	bp
		ret

_hgcline	ENDP

HERC_TEXT	ENDS


HERC_DATA	SEGMENT	word public 'DATA'
PropagatedPixel	DB	00000000b	;0
		DB	11111111b	;1
HERC_DATA	ENDS

		end
