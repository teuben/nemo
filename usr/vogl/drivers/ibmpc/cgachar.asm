;
; Name:		_cgachar
;
; Function:	Display a character in 640*200 2-colour mode
;
; Call:		Microsoft C:
;
;		void
;		_cgachar(c, x, y, fgd, bkgd)
;		int	c;			/* character code */
;		int	x, y;			/* upper left pixel */
;		int	fgd, bkgd;		/* fore and background pixel values */
;
;

argc		equ	word ptr [bp+6]	; stack frame addressing
argx		equ	word ptr [bp+8]
argy		equ	word ptr [bp+10]
argfgd		equ	byte ptr [bp+12]
argbkgd		equ	byte ptr [bp+14]

varmask		equ	[bp-8]
vartoggle	equ	[bp-10]

		extrn	_chartab:word
		extrn	__fontsize:word
		extrn	cgapaddr:far

CGAA_TEXT	segment byte public 'CODE'
		assume	cs:CGAA_TEXT

		public	_cgachar
_cgachar	proc	far

		push	bp		; preserve call registers
		mov	bp,sp
		sub	sp,8		; stack space for local variables
		push	di
		push	si
		push	ds

; set up foreground pixel toggle mask

		mov	ah,argfgd	; AH = 0 or 1 (foreground pixel value)
		ror	ah,1		; high order bit of AH = 0 or 1
		cwd			; propogate high-order bit through DX
		not	dx
		mov	vartoggle,dx

; calculate first pixel address

		mov	ax,argy		; ax = y
		mov	bx,argx		; bx = x
		sub	ax,__fontsize
		call	cgapaddr	; ES:BX -> buffer
					; CL = # bits to shift left

		xor	cl,7		; CL = # bits to rotate right

		mov	ax,0FF00h
		ror	ax,cl		; AX = bit mask in proper position
		mov	varmask,ax

; set up video buffer addressing

		mov	dx,2000h	; increment for video buffer interleave
		mov	di,80-2000h	; increment from last to first interleave
		test	bx,2000h	; set zero flag if BX is 1st interleave
		jz	L01

		xchg	di,dx		; exchange increment values if 1st pixel
					; lies in 1st interleave

; set up character definition table addressing

L01:		mov	ch,byte ptr __fontsize
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
		mov	ax,argc
		mul	ch
		add	si,ax
		
		test	cl,cl		; test # bits to rotate
		jnz	L20		; jump if character is not byte aligned

; routine for byte-aligned characters

		mov	ah,vartoggle	; AH = foreground toggle mask
		xchg	ch,cl		; CX = points

L10:		lodsb			; AL = bit patterb for next pixel row
		xor	al,ah		; toggle pixels if foreground = 0
		mov	es:[bx],al	; store pixels in buffer

		add	bx,dx		; BX = next row in buffer
		xchg	di,dx		; swap buffer increments
		loop	L10
		jmp	short Lexit

; routine for non byte-aligned characters

L20:		mov	ax,varmask
		and	es:[bx],ax	; mask character pixels in buffer

		xor	ah,ah
		lodsb			; AX = bit pattern for next pixel row
		xor	al,vartoggle	; toggle pixels if foreground = 0

		ror	ax,cl		; rotate pixels into position
		or	es:[bx],ax	; store pixels in buffer

		add	bx,dx		; BX = next row in buffer
		xchg	di,dx		; swap buffer increments
		dec	ch
		jnz	L20

Lexit:		pop	ds		; restore caller registers and return
		pop	si
		pop	di
		mov	sp,bp
		pop	bp
		ret

_cgachar	endp

CGAA_TEXT	ends

		end
