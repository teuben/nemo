
;
; Name:		_hgcchar
;
; Function: Display a character in Hercules 720*348 monochrome graphics mode
;
; Call:		Microsoft C:
;
;			void
;			hgcchar(c, x, y, fgd, bkgd)
;				int	c;		/* character code */
;				int	x, y;		/* upper left pixel */
;				int	fgd, bkgd;	/* fore and background
;							  pixel values */
;

ARGC		equ	word ptr [bp+6]	; stack frame addressing
ARGX		equ	word ptr [bp+8]
ARGY		equ	word ptr [bp+10]
ARGFGD		equ	byte ptr [bp+12]
ARGBKGD		equ	byte ptr [bp+14]

VARmask		equ		 [bp-8]
VARtoggle	equ		 [bp-10]
;VAR9bits	equ	byte ptr [bp-12]

		extrn	_chartab:word
		extrn	__fontsize:word

HERC_TEXT	segment byte public 'CODE'
		assume	cs:HERC_TEXT

		extrn	hgcpaddr:far

		public	_hgcchar
_hgcchar	proc	far

		push	bp		; preserve caller registers
		mov	bp,sp
		sub	sp,10		; stack space for local variables
		push	di
		push	si
		push	ds

; calculate first pixel address

		mov	ax,ARGY		; AX = y
		mov	bx,ARGX		; BX = x
		sub	ax,__fontsize	; Make lower left the origin by
					; shifting up by the char size
		call	hgcpaddr	; ES:BX -> buffer
					; CL = # bits to shift left

		xor	cl,7		; CL = # bits to rotate right

; set up 8 bit mask

		mov	ax,0FF00h	; AX = 8-bit mask

    		ror	ax,cl		; AX = bit mask in proper position
		mov	VARmask,ax

; set up foreground pixel toggle mask

		mov	ah,ARGFGD	; AH = 0 or 1 (foreground pixel value)
		ror	ah,1		; high-order bit of AH = 0 or 1
		cwd			; propagate high-order bit through DX
		not	dx		; dx = 0 if foreground = 1
					; or = FFFh if foreground = 0
		mov	ax,VARmask
		not	ax
		and	dx,ax		; zero unused bits of toggle mask in DX
		mov	VARtoggle,dx

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

L20:		mov	ax,VARmask
		and	es:[bx],ax	; mask character pixels in buffer

		xor	ah,ah
		lodsb			; AX = bit pattern for next pixel row

    		ror	ax,cl		; rotate pixels into position
		xor	ax,VARtoggle	; toggle pixels if foreground = 0
		or	es:[bx],ax	; store pixels in buffer
		add	bx,2000h	; increment to next portion of interleave
		jns	L22
		add	bx,90-8000h	; increment to 1st portion of interleave

L22:		dec	ch
		jnz	L20

Lexit:		pop	ds		; restore caller registers and return
		pop	si
		pop	di
		mov	sp,bp
		pop	bp
		ret

_hgcchar	endp

HERC_TEXT	ENDS

		END
