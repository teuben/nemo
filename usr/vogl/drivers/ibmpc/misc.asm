extrn	__buffer_segment:word
extrn	__buffer_offset:word
extrn	__cur_mode:word
extrn	__cur_color:word

MISC_TEXT      SEGMENT byte public 'CODE'
	ASSUME  CS: MISC_TEXT

public	_vega_clear, _cga_clear

_vega_clear proc far
	push bp
	mov bp, sp
        push es
        push ax
        push cx
        push di
	mov dx, 3ceh
	mov al, 0
	out dx,al
	mov dx, 3cfh
	mov al, byte ptr __cur_color
	out dx,al
	mov dx, 3ceh
	mov al, 8
	out dx,al
	mov dx, 3cfh
	mov al, 0ffh
	out dx, al
        mov ax, __buffer_segment
        mov es, ax
        xor di, di
	mov cx, 14000		; How many words on an EGA 640x350 mode.
	cmp __cur_mode, 18	; Is it a VGA ?
	jne not_vga
	mov cx, 19200		; How many words on a VGA 640x480 mode.
not_vga:
        mov ax, 0
        cld
        rep stosw                     
        pop di
        pop cx
        pop ax
        pop es
	mov sp, bp
	pop bp
        ret
_vega_clear endp

_cga_clear proc far
	push bp
	mov bp, sp
        push es
        push ax
        push cx
        push di
	mov ax, __buffer_segment
        mov es, ax
	mov di, __buffer_offset
        mov cx, 2000h
	cmp __cur_color, 0
	jne  white
        xor ax, ax		; Zero the screen
	jmp next
white:
	mov ax,0ffffh		; White the screen
next:
        cld
        rep stosw                     
        pop di
        pop cx
        pop ax
        pop es
	mov sp, bp
	pop bp
        ret
_cga_clear endp

MISC_TEXT	ENDS
END
