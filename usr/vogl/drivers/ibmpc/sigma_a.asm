	TITLE	SIG_LINE - Line drawing routine.
	NAME	SIG_LINE
	PAGE	55,132


	COMMENT	$

	Name:		SIG_LINE

	Function:	Draw a line in 640x400 SIGMA graphics mode.

	Caller:		Microsoft C:

			void 	sigline(x1, y1, x2, y2, n);
				int	x1, y1, x2, y2;	/* pixel co-ords */
				int	n;		/* pixel value */
		$



SIGMA_DATA	SEGMENT	WORD	PUBLIC	'DATA'
	color	dw	?

row_table	label	word			; define row table
	line_num = 0
	rept	100
		dw	line_num * 80
		dw	2000h + line_num * 80
		dw	4000h + line_num * 80
		dw	6000h + line_num * 80
		line_num = line_num + 1		; next row
	endm

	negerr	dw	?
	poserr	dw	?
SIGMA_DATA	ENDS

extrn	__cur_color:word
extrn	_chartab:word
extrn	__fontsize:word

argbase	equ	2			; 2 for large model, 0 for small
ARGx1	equ	[bp + argbase + 4]
ARGy1	equ	[bp + argbase + 6]
ARGx2	equ	[bp + argbase + 8]
ARGy2	equ	[bp + argbase + 10]
ARGcol	equ	[bp + argbase + 12]

inc_di	equ	47h
inc_cx	equ	41h
dec_di	equ	4fh
dec_cx	equ	49h
do_nop	equ	90h

public	_sig_line
public	_sigmachar
public	_sigmaclear
public	_sigma_set_colors

DGROUP	GROUP	SIGMA_DATA

SIGMA_TEXT	SEGMENT	byte public 'CODE'
		ASSUME	cs:SIGMA_TEXT, ds:DGROUP


_sig_line	proc	far
	push	bp
	mov	bp, sp
	push	di
	push	si
	push	bx
	push	cx
	push	dx

	mov	ax, ARGcol
	mov	color, ax
	mov	ax, 0b800h			; display mem segment
	mov	es, ax				; display mem segment
	
	mov	bx, inc_di * 256 + inc_cx	; ystep + 1 in bh
						; xstep + 1 in bl
; | x2 - x1 |
	mov	cx, ARGx2			; x2
	sub	cx, ARGx1			; x1
	jge	pos_delta_x			; jump if not negative
	mov	bl, dec_cx			; cx = delta_x
	neg	cx

pos_delta_x:
; | y2 - y1 |
;	mov	bh, inc_di
	mov	dx, ARGy2			; y2
	sub	dx, ARGy1			; y1
	jge	pos_delta_y
	mov	bh, dec_di			; ystep - 1
	neg	dx				; | y2 - y1 |

pos_delta_y:					; dx = delta_y
; self modifying code
	mov	di, offset cs:positive
	mov	word ptr cs:[di], bx		; self mod

; find long axis
	cmp	cx, dx				; delta_x, delta_y
	jge	big_delta_x
	mov	bl, do_nop			; if vert, no inc/dec so nop
	xchg	cx, dx				; swap delta_x and delta_y
	jmp	short got_long_axis

big_delta_x:
	mov	bh, do_nop			; if horiz, no inc/dec so nop

got_long_axis:
; self modifying code
	mov	di, offset cs:negative		; location to modify
	mov	word ptr cs:[di], bx		; make mod.
	shl	dx, 1				; (delta_y * 2)
	mov	negerr, dx
	sub	dx, cx				; (delta_y * 2) - delta_x
	mov	bx, dx				; save initial value
	sub	dx, cx				; (delta_y * 2) - (delta_x * 2)
	mov	poserr, dx
	mov	dx, 02DEh			; plane select port

; adjust count
	mov	si, cx				; si = length of long axis
	cmp	si, 0
	jne	not_0
	inc	si
not_0:
	mov	di, ARGy1			; y
	mov	cx, ARGx1			; x
	push	bp
	mov	bp, negerr

main_loop:
	push	di				; save y
	push	cx				; save x
	mov	ax, cx				; x
	shl	di, 1				; word table index
	mov	di, row_table[di]
	shr	ax, 1
	shr	ax, 1
	shr	ax, 1
	add	di, ax
	and	cl, 7
	mov	al, 10000000b			; mask
	ror	al ,cl
	
	mov	cl, al				; bitORmask
	mov	ch, al				
	not	ch				; bitANDmask

	xor	al, al
	out	dx, al				; plane 0
	
	mov	ah, byte ptr color
	and	ah, 1
	jz	np1a
	or	es:[di], cl
	jmp	short np1b
np1a:
	and	es:[di], ch
np1b:
	mov	al, 1
	out	dx, al				; plane 1
	
	mov	ah, byte ptr color
	and	ah, 2
	jz	np2a
	or	es:[di], cl
	jmp	short np2b
np2a:
	and	es:[di], ch
np2b:
	mov	al, 2
	out	dx, al				; plane 2
	
	mov	ah, byte ptr color
	and	ah, 4
	jz	np3a
	or	es:[di], cl
	jmp	short np3b
np3a:
	and	es:[di], ch
np3b:
	mov	al, 3
	out	dx, al				; plane 3
	
	mov	ah, byte ptr color
	and	ah, 8
	jz	np4a
	or	es:[di], cl
	jmp	short np4b
np4a:
	and	es:[di], ch
np4b:

	pop	cx
	pop	di

	cmp	bx, 0
	jge	positive

negative:
	inc	cx				; update x
						; inc cx / dec cx / nop
						; update y
	inc	di				; inc di / dec di / nop
	
	add	bx, bp				; negerr

dec_long_axis:
	dec	si				; dec count
	jge	main_loop
	jmp	short end_line			; end of line

positive:
	inc	cx				; update x
						; inc cx / dec cx
	inc	di				; update y
						; inc di / dec di
	add	bx, poserr			; (dx)
	dec	si				; dec count
	jz	short end_line
	jmp	main_loop			; not end of line

end_line:
	pop	bp
	pop	dx
	pop	cx
	pop	bx
	pop	si
	pop	di
	mov	sp, bp
	pop	bp
	ret
_sig_line	endp

;
;			void
;			sigmachar(c, x, y, col)
;				int	c;		/* character code */
;				int	x, y;		/* upper left pixel */
;				int	col		/* colour	*/
;

ARGC		equ	word ptr [bp + argbase + 4]	; stack frame addressing
ARGX		equ	word ptr [bp + argbase + 6]
ARGY		equ	word ptr [bp + argbase + 8]
ARGcol		equ	         [bp + argbase + 10]

_sigmachar	proc	far

	push	bp		; preserve caller registers
	mov	bp,sp
	push	di
	push	si
	push	ds

	mov	ax, ARGcol
	mov 	color, ax

	mov	ax, 0b800h			; display mem segment
	mov	es, ax				; display mem segment

; calculate first pixel address

	mov	di,ARGY		; di = y
	mov	cx,ARGX		; aX = x
	mov	ax, cx
	sub	di,__fontsize	; Make lower left the origin by
	shl	di, 1		; Index into row table
	mov	di, row_table[di]
	shr	ax, 1
	shr	ax, 1
	shr	ax, 1
	add	di, ax
	and	cl, 7

; set up character definition table addressing

	mov	ch,byte ptr __fontsize
	cmp	ch,16
	je	largel

	mov	ax,0f000h
	mov	ds,ax
	mov	ax,0fa6eh
	mov	si,ax
	
	jmp	cont1
largel:
	mov	ax, SEG _chartab
	mov	ds, ax
	mov	si, OFFSET _chartab
cont1:
	mov	ax,ARGC
	mul	ch
	add	si,ax
		
; mask and set pixels in the video buffer


L20:
	xor	ax,ax
	lodsb			; AX = bit pattern for next pixel row

	ror	ax, cl		; rotate pixels into position
	mov	dx, ax		; bitORmask
	mov	bx, ax
	not	bx		; bitANDmask

	xor	al, al
	push	dx
	mov	dx, 02DEh
	out	dx, al		; plane 0
	pop	dx
	
	mov	ah, byte ptr ARGcol
	and	ah, 1
	jz	np11a
	or	es:[di], dx
	jmp	short np11b
np11a:
	and	es:[di], bx
np11b:
	mov	al, 1
	push	dx
	mov	dx, 02DEh
	out	dx, al		; plane 1
	pop	dx
	
	mov	ah, byte ptr ARGcol
	and	ah, 2
	jz	np22a
	or	es:[di], dx
	jmp	short np22b
np22a:
	and	es:[di], bx
np22b:
	mov	al, 2
	push	dx
	mov	dx, 02DEh
	out	dx, al		; plane 2
	pop	dx
	
	mov	ah, byte ptr ARGcol
	and	ah, 4
	jz	np33a
	or	es:[di], dx
	jmp	short np33b
np33a:
	and	es:[di], bx
np33b:
	mov	al, 3
	push	dx
	mov	dx, 02DEh
	out	dx, al		; plane 3
	pop	dx
	
	mov	ah, byte ptr ARGcol
	and	ah, 8
	jz	np44a
	or	es:[di], dx
	jmp	short np44b
np44a:
	and	es:[di], bx
np44b:

	add	di, 2000h
	jns	L22
	add	di, 80-8000h

L22:	dec	ch
	jnz	L20
Lexit:
	pop	ds		; restore caller registers and return
	pop	si
	pop	di
	mov	sp,bp
	pop	bp
	ret

_sigmachar	endp

_sigmaclear	proc	far
	push	bp
	mov	bp, sp
	push	di
	push	si
	push	cx
	push	dx

	mov	ax, 0b800h			; display mem segment
	mov	es, ax				; display mem segment
	xor	di, di
        mov 	cx, 04000h
	mov	dx, 02DEh

	mov	al, 0				; plane 0
	out	dx, al
	mov	ah, byte ptr __cur_color
	and	ah, 1
	jz	n1
	mov	ax, 0FFFFh
	jmp	short	n2
n1:
	xor	ax, ax
n2:

        cld
        rep	stosw                     
	xor	di, di
        mov 	cx, 04000h

	mov	al, 1				; plane 1
	out	dx, al
	mov	ah, byte ptr __cur_color
	and	ah, 2
	jz	n3
	mov	ax, 0FFFFh
	jmp	short	n4
n3:
	xor	ax, ax
n4:
        rep stosw                     
	xor	di, di
        mov 	cx, 04000h

	mov	al, 2				; plane 2
	out	dx, al
	mov	ah, byte ptr __cur_color
	and	ah, 4
	jz	n5
	mov	ax, 0FFFFh
	jmp	short	n6
n5:
	xor	ax, ax
n6:
        rep stosw                     
	xor	di, di
        mov 	cx, 04000h

	mov	al, 3				; plane 2
	out	dx, al
	mov	ah, byte ptr __cur_color
	and	ah, 8
	jz	n7
	mov	ax, 0FFFFh
	jmp	short	n8
n7:
	xor	ax, ax
n8:
        rep stosw                     

        pop di
        pop si
        pop cx
        pop dx
	mov sp, bp
	pop bp
        ret
_sigmaclear endp

_sigma_set_colors	proc	far
         push    bp
         mov     bp, sp
         les     dx,[bp+argbase+4]
         mov    al,2
         mov    ah,16;

         int     10H
         pop     bp   
         ret
_sigma_set_colors	endp

SIGMA_TEXT	ends

		end

