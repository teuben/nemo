	TITLE	VEGABUF	- Handles VGA/EGA double buffering / page swapping
	NAME	VEGABUF
	PAGE	55,132

	COMMENT	$
	
	Caller:	Microsoft C

		vega_backbuf()	/* Displays page 1, Sets page 0 for writing */

		void
		vega_frontbuf()	/* Displays page 0, Sets page 1 for writing */

		void
		vega_swapbuf()	/* Swaps page 0 and page 1 on display */
		
	$
		
displayseg0		equ		0A000h
displayseg1		equ		displayseg0 + 0800h

PUBLIC	_vega_backbuf, _vega_frontbuf, _vega_swapbuf

	extrn	__buffer_segment:word
	extrn	__buffer_offset:word
	extrn	__cur_mode:word

vega_text	segment byte public 'code'

	assume cs:vega_text

	comment	$
		displays page0, sets page 1 for writing.
	$

_vega_backbuf	proc	far
	push		bp
	mov		bp, sp

	mov		dx, 3dah	; status reg
vt1:	
	in		al, dx
	test		al, 8		; test for vert retrace
	jz		vt1

	mov		ax, 0500h	; Show page 0, write to page 1
        int		010h

	mov		__buffer_segment, displayseg1
	cmp		__cur_mode, 18	; Is this VGA mode?
	jne		not_vga
	add		__buffer_segment, 0800h
not_vga:

	mov		ax, 1	; return value
	pop		bp
	ret
_vega_backbuf		endp	

	comment	$
		displays page0, sets page 0 for writing.
		Makes sure that written page is the display page.
	$

_vega_frontbuf	proc	far
	push		bp
	mov		bp, sp

	mov		dx, 3dah	; status reg
vt2:	
	in		al, dx
	test		al, 8		; test for vert retrace
	jz		vt2

	mov		ax, 0500h	; show page 0, write to page 0
        int		010h

	mov		__buffer_segment, displayseg0

	pop		bp
	ret
_vega_frontbuf		endp	


	comment	$
		display current buffer (__buffer_segment)
		set drawing in other buffer.
	$

_vega_swapbuf	proc	far
	push		bp
	mov		bp, sp

	mov		dx, 3dah	; status reg

	cmp		__buffer_segment, displayseg0
	je		disp_page0
;  Written page is page 1, displayed page is page 0
;  so display page 1 and write to page 0.

vt3:	
	in		al, dx
	test		al, 8		; test for vert retrace
	jz		vt3

        mov		ax, 0501h	; Show page 1
        int		010h

	mov		__buffer_segment, displayseg0
	jmp		go

disp_page0:

vt4:	
	in		al, dx
	test		al, 8		; test for vert retrace
	jz		vt4

	mov		ax, 0500h
        int		010h

	mov		__buffer_segment, displayseg1
	cmp		__cur_mode, 18
	jne		go
	add		__buffer_segment, 0800h
go:
	pop		bp
	ret
_vega_swapbuf		endp	

vega_text ends

	end
