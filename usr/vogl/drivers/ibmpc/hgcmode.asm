
	TITLE	HGCMODE	- Sets Hercules Graphics mode (on page 1)
	NAME	HGCMODE
	PAGE	55,132

	COMMENT	$
	
	Establish Various Hercules modes

	Caller:	Microsoft C

		void
		hgc_config()	/* Set Hercules to "full" config */

		void
		hgc_gmode()	/* Set Graphics mode (on page 1) */

		void
		hgc_tmode()	/* Set Text mode */

		void
		_hgc_clear()	/* Clear the graphics screen */

		void
		hgc_backbuf()	/* Displays page 1, Sets page 0 for writing */

		void
		hgc_frontbuf()	/* Displays page 0, Sets page 1 for writing */

		void
		hgc_swapbuf()	/* Swaps page 0 and page 1 on display */
		
	$
		
config_switch		equ 		03bfh
TextSeg			equ		0b000h
DisplaySeg0		equ		0b000h
DisplaySeg1		equ		0b800h
page0			equ		00000000b
page1			equ		10000000b
dmc_port		equ		03b8h
index_6845		equ		03b4h
disp_stat_port		equ		03bah
text			equ		00000000b
graphics		equ		00000010b
off_screen		equ		00000000b
on_screen		equ		00001000b

public _hgc_config, _hgc_clear, _hgc_gmode, _hgc_tmode
PUBLIC	_hgc_backbuf, _hgc_frontbuf, _hgc_swapbuf

herc_data segment	word public	'data'
ttable		db	61h,50h,52h,0fh,19h,06h,19h,19h,02h,0dh,0bh,0ch
gtable		db	35h,2dh,2eh,07h,5bh,02h,57h,57h,02h,03h,00h,00h
herc_data ends

	extrn	__buffer_segment:word
	extrn	__cur_color:word

herc_text	segment byte public 'code'

	dgroup	group	herc_data
	assume cs:herc_text, ds:dgroup

wait_v_retrace	macro
	local		in_v_retrace, not_in_v_retrace
	mov		dx,disp_stat_port
in_v_retrace:
	in		al,dx
	shl		al,1
	jnc		in_v_retrace
not_in_v_retrace:
	in		al,dx
	shl		al,1
	jc		not_in_v_retrace
	endm
	

vertical_retrace	proc far
	wait_v_retrace
	ret
vertical_retrace	endp

_hgc_config proc far
	push	bp
	mov		bp,sp
	mov		al,3
	mov		dx,config_switch
	out		dx,al
	mov		sp,bp
	pop		bp
	ret
_hgc_config	endp
	
_hgc_clear proc far
	push		bp
	mov		bp,sp
	push		di
	push		cx
	pushf
	mov		cx,__buffer_segment
	mov		es,cx
	xor		di,di
	cmp		__cur_color, 0
	jne		white
	xor		ax,ax
	jmp		next
white:
	mov		ax,0ffffh
next:
	mov		cx,16*1024
	cld
	rep		stosw
	popf
	pop		cx
	pop		di
	mov		ax, 14
	mov		sp,bp
	pop		bp
	ret	
_hgc_clear		endp

setmod	proc		far
	pushf
	mov		dx,dmc_port	
	out		dx,al
	mov		dx,index_6845
	mov		cx,12
	xor		ah,ah
	cld
next_6845:
	mov		al,ah
	out		dx,al
	inc		dx
	lodsb
	out		dx,al
	inc		ah
	dec 		dx
	loop		next_6845
	popf
	ret
setmod	endp

_hgc_tmode	proc far
	push		bp
	mov		bp,sp
	sub		sp, 10
	push		di	
	push		si
	push		bx
	push		cx
	push		dx
	lea		si,ttable
	call		setmod
	call		vertical_retrace
	xor		ax,ax
	mov		al,page0 + on_screen
	add		al,text
	mov		dx,dmc_port
	out		dx,al
	mov		__buffer_segment, displayseg0
	pop		dx
	pop		cx
	pop		bx
	pop		si
	pop		di
	mov		sp,bp
	pop		bp
	ret
_hgc_tmode	endp

_hgc_gmode	proc far
	push		bp
	mov		bp,sp
	sub		sp,10
	push		di	
	push		si
	push		bx
	push		cx
	push		dx
	xor		ax,ax
	mov		al,graphics
	lea		si,gtable
	call		setmod
	mov		al,page1 + on_screen
	add		al,graphics
	push		ax
 	call		vertical_retrace
	pop		ax
	mov		dx,dmc_port
	out		dx,al
	mov		__buffer_segment, displayseg1
	pop		dx	
	pop		cx	
	pop		bx	
	pop		si
	pop		di
	mov		sp,bp
	pop		bp
	ret
_hgc_gmode	endp

	comment	$
		displays page0 (internal routine)
	$

hgc_page0	proc	far
	push		bp
	mov		bp, sp

	push		dx
	call		vertical_retrace
	mov		dx, dmc_port
	mov		ax, page0
	add		al, graphics + on_screen
	out		dx, al
	pop		dx
	mov		sp, bp
	pop		bp
	ret
hgc_page0		endp

	comment	$
		displays page1 (internal routine).
	$

hgc_page1	proc	far
	push		bp
	mov		bp, sp

	push		dx
	call		vertical_retrace
	mov		dx, dmc_port
	mov		ax, page1
	add		al, graphics + on_screen
	out		dx, al
	pop		dx
	mov		sp, bp
	pop		bp
	ret
hgc_page1		endp
	
	
	comment	$
		displays page1, sets page 0 for writing.
	$

_hgc_backbuf	proc	far
	push		bp
	mov		bp, sp

	call		hgc_page1
	mov		__buffer_segment, displayseg0

	mov		ax, 1	; return value
	mov		sp, bp
	pop		bp
	ret
_hgc_backbuf		endp	

	comment	$
		displays page1, sets page 1 for writing.
		(This ensures that written page is the display page).
	$

_hgc_frontbuf	proc	far
	push		bp
	mov		bp, sp

	call		hgc_page1
	mov		__buffer_segment, displayseg1

	mov		sp, bp
	pop		bp
	ret
_hgc_frontbuf		endp	


	comment	$
		display current buffer (page)
		set drawing in other buffer (page)
	$

_hgc_swapbuf	proc	far
	push		bp
	mov		bp, sp

	cmp		__buffer_segment, displayseg0
	je		disp_page0
	call		hgc_page1
	mov		__buffer_segment, displayseg0
	jmp		go

disp_page0:
	call		hgc_page0
	mov		__buffer_segment, displayseg1

go:
	mov		sp, bp
	pop		bp
	ret
_hgc_swapbuf		endp	

herc_text ends

	end
