
MOUSE_TEXT	SEGMENT  BYTE PUBLIC 'CODE'
MOUSE_TEXT	ENDS
MOUSE_DATA	SEGMENT  WORD PUBLIC 'DATA'
MOUSE_DATA	ENDS
CONST	SEGMENT  WORD PUBLIC 'CONST'
CONST	ENDS
_BSS	SEGMENT  WORD PUBLIC 'BSS'
_BSS	ENDS
DGROUP	GROUP	CONST,	_BSS,	MOUSE_DATA
	ASSUME  CS: MOUSE_TEXT, DS: DGROUP, SS: DGROUP, ES: DGROUP

MOUSE_TEXT      SEGMENT
public	_ismouse, _showmouse, _unshowmouse, _readmouse, _set_loc

;
; The Microsoft mouse manual says that you should bung a 5 at location
; 40h:49h (Ie. The spot that says what mode we're in) for the mouse to work
; on a HERCULES Card.
;
_set_loc   proc	far
	push	bp
	mov	bp, sp
        push	si
        push	di
	push    bx

        mov 	bx,40h
        mov 	es,bx
        mov 	bx,49h
        mov 	si,bx
	mov	ah,es:[si]
        mov     ah, byte ptr [bp + 6]
        mov 	byte ptr es:[si],ah            ; Bung something there
 
	mov	ax, [bp + 6]
	pop 	bx
        pop 	di
        pop 	si
        mov     sp, bp
        pop     bp

        ret
_set_loc   endp

;
; Check if the mouse driver is present, if it is, then set
; the x and y limits for it.
;
_ismouse        proc far
        push    bp
        push    es
        mov     bp,sp

        mov     ax,03533h       ;Get int 33h by calling int 21
        int     21
        mov     ax,es           ;Check segment, offset of int 33
        or      ax,bx           ;vector. If 0 or pointing to IRET
        jnz     yesdrv          ;then driver not installed.
        mov     bl, es:[bx]
        cmp     bl,0cfh
        jne     yesdrv

        mov     ax,0            ; return 0 (false) if no driver
        jmp     pissoff

yesdrv: 
        mov     ax,0            ; Initialize mouse
        int     33h
        cmp     ax,0            ; Is mouse installed 
        jz      pissoff         ; pissoff if not installed

	mov	ax,7		;Function 7, set x limits
	mov	cx,0		
	mov	dx,[bp+6]
	int	33h
	mov	ax,8		;Function 8, set y limits
	mov	cx,0			
	mov	dx,[bp+8]	
	int	33h
        mov     ax,1            ; return 1 (true) if driver installed

pissoff:                        ; pissoff back to caller.
        mov     sp,bp
        pop     es
        pop     bp
        ret     

_ismouse        endp

;
; Turns the mouse cursor on
;
_showmouse	proc far
        push    bp
        mov     bp,sp
	mov	ax,1
	int 	33h
	mov	sp,bp
	pop	bp
	ret
_showmouse	endp

;
; Turns the mouse cursor off
;
_unshowmouse	proc far
        push    bp
        mov     bp,sp
	mov	ax,2
	int 	33h
	mov	sp,bp
	pop	bp
	ret
_unshowmouse	endp

;
;	readmouse(xmouse, ymouse)
;
_readmouse      proc far
        push    bp
        mov     bp,sp
	push	di
	push	bx
	push	cx
	push	dx
;
;	Get mouse status.
;
	mov     ax,3
	int     33h
	
    	mov	ax,bx				 ; return buttons 

;  	 mov     di,[bp+4]            ;x  <= this is for small mem model
;        mov     word ptr [di],cx
;        mov     di,[bp+6]            ;y
;        mov     word ptr [di],dx

	les	di,[bp+6]		; <= this is for large mem model
	mov	word ptr es:[di],cx		;x
	les	di,[bp+10]
	mov	word ptr es:[di],dx		;y

	pop	dx
	pop	cx
	pop	bx
	pop	di
        mov     sp,bp
        pop     bp
        ret     

_readmouse      endp

MOUSE_TEXT	ENDS
END
