	TITLE	CGALINE - Fast line drawing routine.
	NAME	CGALINE
	PAGE	55,132

	COMMENT	$

	Name:		CGALINE

	Function:	Draw a line in 640x200 2-color mode

	Caller:		Microsoft C:

			void 	cgaline(x1, y1, x2, y2, n);

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

DGROUP		GROUP	CGAA_DATA

		EXTRN	cgapaddr:far

CGAA_TEXT	SEGMENT	byte public 'CODE'
		ASSUME	cs:CGAA_TEXT,ds:DGROUP

		PUBLIC	_cgaline

_cgaline	PROC	far

		push	bp	; Set up stack frame
		mov	bp,sp
		sub	sp,12	; Stack space for local variables
		push	si
		push	di

		mov	si,2000h	; Increment for video buffer interleave
		mov	di,80-2000h	; Increment for last to first interleave
; check for vertical line

		mov	cx,ARGx2
		sub	cx,ARGx1
		jz	VertLine

; force x1 < x2

		jns	L01	; jump if x2 > x1

		neg	cx	; CX = x1 - x2

		mov	bx,ARGx2; exchange x1 and x2
		xchg	bx,ARGx1
		mov	ARGx2,bx

		mov	bx,ARGy2; exchange y1 and y2
		xchg	bx,ARGy1
		mov	ARGy2,bx

; calc dy = abs(y2 - y1)

L01:
		mov	bx,ARGy2
		sub	bx,ARGy1; BX = y2 - y1
		jnz	L02
		jmp	HorizLine
L02:
		jns	L03

		neg	bx	; BX = y1 - y2
		neg	si	; negate increment for buffer interleave
		neg	di
		xchg	si,di

; select appropriate routine for slope of line

L03:
		mov	VARleafincr,di; save increment for buffer interleave
		mov	VARroutine,offset LoSlopeLine
		cmp	bx,cx
		jle	L04	; jump if dy <= dx (slopr <= 1)
		mov	VARroutine,offset HiSlopeLine
		xchg	bx,cx	; exchange dy and dx

; calc initial decision variable and increments

L04:
		shl	bx,1		; BX = 2 * dy
		mov	VARincr1,bx	; incr1 = 2 * dy
		sub	bx,cx
		mov	di,bx		; DI = d = 2*dy - dx
		sub	bx,cx
		mov	VARincr2,bx	; incr2 = 2*(dy - dx)

; calc first pixel address

		push	cx			; preserve the register
		mov	ax,ARGy1	; ax = y
		mov	bx,ARGx1	; bx = x
		call	cgapaddr	; ah = bit mask
					; es:bx -> buffer
					; CL = # bits to shift left
		mov	al,ARGn		; AL = unshifted pixel value
		shl	ax,cl		; AH = bitmask in proper position
					; AL = pixel value in proper position
		mov	dx,ax		; DH = bitmask
					; DL = pixel value

		not	dh		; DH = inverse bitmask
		pop	cx
		inc	cx
		test	bx,2000h	; set zero flag if BX in 1st interleave
		jz	L05
		xchg	si,VARleafincr	; exchange increament 
					; values if 1st pixel lies in 1st interleave
		
L05:		jmp	VARroutine	; Jump to the appropiate routine for slope

; routine for verticle lines

VertLine:
		mov	ax,ARGy1
		mov	bx,ARGy2
		mov	cx,bx
		sub	cx,ax
		jge	L31
		neg	cx
		mov	ax,bx

L31:
		inc	cx
		mov	bx,ARGx1
		push	cx
		call	cgapaddr

		mov	al,ARGn
		shl	ax,cl
		not	ah
		pop	cx

		test	bx,si
		jz	L32

		xchg	si,di

L32:
		test	al,al
		jz	L34
L33:
		or	es:[bx],al
		add	bx,si
		xchg	si,di
		loop	L33
		jmp	short	L35
L34:
		and	es:[bx],ah
		add	bx,si
		xchg	si,di
		loop	L34

L35:
		jmp	Lexit

; routine for horizontal line

HorizLine:
		mov	ax,ARGy1
		mov	bx,ARGx1
		call	cgapaddr
		mov	di,bx
		mov	dh,ah
		not	dh
		mov	dl,0FFh

		shl	dh,cl
		not	dh

		mov	cx,ARGx2
		and	cl,7
		xor	cl,7
		shl	dl,cl

; determine byte offset of first and last pixel in line

		mov	ax,ARGx2
		mov	bx,ARGx1
		mov	cl,ByteOffsetShift

		shr	ax,cl
		shr	bx,cl
		mov	cx,ax
		sub	cx,bx

; propagate pixel value throughout one byte

		mov	bx,offset DGROUP:PropagatedPixel
		mov	al,ARGn
		xlat

; set pixels in leftmost byte of line

		or	dh,dh
		js	L43
		or	cx,cx
		jnz	L42
		and	dl,dh
		jmp	short L44

L42:
		mov	ah,al
		and	ah,dh
		not	dh
		and	es:[di],dh
		or	es:[di],ah
		inc	di
		dec	cx

; draw remainder of the line

L43:
		rep	stosb

; set pixels in rightmost byte of line

L44:
		and	al,dl
		not	dl
		and	es:[di],dl
		or	es:[di],al
		jmp	Lexit



; routine for dy >= dx (slope <= 1)

LoSlopeLine:

L10:
		mov	ah,es:[bx]

L11:
		and	ah,dh
		or	ah,dl
		ror	dl,1
		ror	dh,1
		jnc	L14

; bit mask not shifted out

		or	di,di
		jns	L12
		add	di,VARincr1
		loop	L11

		mov	es:[bx],ah
		jmp	short Lexit

L12:
		add	di,VARincr2
		mov	es:[bx],ah
		add	bx,si
		xchg	si,VARleafincr
		loop	L10
		jmp	short Lexit

; bit mask shifted out

L14:
		mov	es:[bx],ah
		inc	bx

		or	di,di
		jns	L15
		add	di,VARincr1
		loop	L10
		jmp	short Lexit

L15:
		add	di,VARincr2
		add	bx,si
		xchg	si,VARleafincr
		
		loop	L10
		jmp	short Lexit



; routine for dy > dx (slope > 1)

HiSlopeLine:

L21:
		and	es:[bx],dh
		or	es:[bx],dl

		add	bx,si
		xchg	si,VARleafincr

L22:
		or	di,di
		jns	L23

		add	di,VARincr1
		loop	L21

		jmp	short Lexit

L23:
		add	di,VARincr2
		ror	dl,1
		ror	dh,1
		cmc

		adc	bx,0

		loop	L21


Lexit:
		pop	di
		pop	si
		mov	sp,bp
		pop	bp
		ret

_cgaline	ENDP

CGAA_TEXT	ENDS


CGAA_DATA			SEGMENT	word public 'DATA'
PropagatedPixel	DB	00000000b	;0
		DB	11111111b	;1
CGAA_DATA			ENDS

		END
