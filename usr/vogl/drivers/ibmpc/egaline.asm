	TITLE	EGALINE - Fast line drawing routine.
	NAME	EGALINE
	PAGE	55,132


	COMMENT	$

	Name:		EGALINE

	Function:	Draw a line in the following EGA/VGA graphics modes:
			200 line 16 colour modes
			350 line modes
			640 x 480 16 colour modes

	Caller:		Microsoft C:

			void 	egaline(x1, y1, x2, y2, n);
				int	x1, y1, x2, y2;	/* pixel co-ords */
				int	n;		/* pixel value */
		$


; Stack frame addressing - LARGE CODE MODEL

ARGx1		EQU	word ptr [bp+6]
ARGy1		EQU	word ptr [bp+8]
ARGx2		EQU	word ptr [bp+10]
ARGy2		EQU	word ptr [bp+12]
ARGn		EQU	byte ptr [bp+14]

VARvertincr	EQU	word ptr [bp-6]
VARincr1	EQU	word ptr [bp-8]
VARincr2	EQU	word ptr [bp-10]
VARroutine	EQU	word ptr [bp-12]

ByteOffsetShift	EQU	3
BytesPerLine	EQU	80
RMWbits		EQU	0


VEGA_TEXT	SEGMENT	byte public 'CODE'
		ASSUME	cs:VEGA_TEXT

		EXTRN	egapaddr:far
		PUBLIC	_egaline

_egaline	PROC	far

		push	bp		; Set up stack frame
		mov	bp,sp
		sub	sp,10
		push	si
		push	di

; config graphics controller

		mov	dx,3CEH
		mov	ah,ARGn
		xor	al,al
		out	dx,ax

		mov	ax,0F01H
		out	dx,ax

		mov	ah,RMWbits
		mov	al,3
		out	dx,ax

; check for vertical line

		mov	si,BytesPerLine
		mov	cx,ARGx2
		sub	cx,ARGx1
		jz	VertLine

; force x1 < x2

		jns	L01

		neg	cx

		mov	bx,ARGx2
		xchg	bx,ARGx1
		mov	ARGx2,bx

		mov	bx,ARGy2
		xchg	bx,ARGy1
		mov	ARGy2,bx

; calc dy = abs(y2 - y1)

L01:
		mov	bx,ARGy2
		sub	bx,ARGy1
		jz	HorizLine
		jns	L03

		neg	bx
		neg	si

; select appropriate routine for slope of line

L03:
		mov	VARvertincr,si
		mov	VARroutine,offset LoSlopeLine
		cmp	bx,cx
		jle	L04
		mov	VARroutine,offset HiSlopeLine
		xchg	bx,cx

; calc initial decision variable and increments

L04:
		shl	bx,1
		mov	VARincr1,bx
		sub	bx,cx
		mov	si,bx
		sub	bx,cx
		mov	VARincr2,bx

; calc first pixel address

		push	cx
		mov	ax,ARGy1
		mov	bx,ARGx1
		call	egapaddr
		mov	di,bx
		shl	ah,cl
		mov	bl,ah
		mov	al,8
		pop	cx
		inc	cx
		jmp	VARroutine

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
		call	egapaddr

; set up graphics controller

		shl	ah,cl
		mov	al,8
		out	dx,ax
		pop	cx

; draw the line

L32:
		or	es:[bx],al
		add	bx,si
		loop	L32
		jmp	Lexit

; routine for horizontal line

HorizLine:
		push	ds

		mov	ax,ARGy1
		mov	bx,ARGx1
		call	egapaddr
		mov	di,bx
		mov	dh,ah
		not	dh
		shl	dh,cl
		not	dh

		mov	cx,ARGx2
		and	cl,7
		xor	cl,7
		mov	dl,0FFH
		shl	dl,cl

; determine byte offset of first and last pixel in line

		mov	ax,ARGx2
		mov	bx,ARGx1
		mov	cl,ByteOffsetShift

		shr	ax,cl
		shr	bx,cl
		mov	cx,ax
		sub	cx,bx

; get graphics controller port address into DX

		mov	bx,dx
		mov	dx,3CEH
		mov	al,8

; make video buffer addressable through DS:SI

		push	es
		pop	ds
		mov	si,di

; set pixels in leftmost byte of line

		or	bh,bh
		js	L43
		or	cx,cx
		jnz	L42
		and	bl,bh
		jmp	short L44

L42:
		mov	ah,bh
		out	dx,ax
		movsb
		dec	cx

; draw remainder of the line

L43:
		mov	ah,0FFH
		out	dx,ax
		rep	movsb

; set pixels in rightmost byte of line

L44:
		mov	ah,bl
		out	dx,ax
		movsb
		pop	ds
		jmp	short Lexit



; routine for dy >= dx (slope <= 1)

LoSlopeLine:

L10:
		mov	ah,bl

L11:
		or	ah,bl
		ror	bl,1
		jc	L14

; bit mask not shifted out

		or	si,si
		jns	L12
		add	si,VARincr1
		loop	L11

		out	dx,ax
		or	es:[di],al
		jmp	short Lexit

L12:
		add	si,VARincr2
		out	dx,ax
		or	es:[di],al
		add	di,VARvertincr
		loop	L10
		jmp	short Lexit

; bit mask shifted out

L14:
		out	dx,ax
		or	es:[di],al
		inc	di
		or	si,si
		jns	L15
		add	si,VARincr1
		loop	L10
		jmp	short Lexit

L15:
		add	si,VARincr2
		add	di,VARvertincr
		loop	L10
		jmp	short Lexit



; routine for dy > dx (slope > 1)

HiSlopeLine:
		mov	bx,VARvertincr

L21:
		out	dx,ax
		or	es:[di],al
		add	di,bx

L22:
		or	si,si
		jns	L23

		add	si,VARincr1
		loop	L21
		jmp	short Lexit

L23:
		add	si,VARincr2
		ror	ah,1
		adc	di,0
		loop	L21

; restore default graphics controoler state and return to caller

Lexit:

COMMENT $
		xor	ax,ax
		out	dx,ax
		inc	ax
		out	dx,ax
		mov	ax,3
		out	dx,ax
		mov	ax,0FF08H
		out	dx,ax
		$

		pop	di
		pop	si
		mov	sp,bp
		pop	bp
		ret

_egaline	endp

VEGA_TEXT	ends

		end

