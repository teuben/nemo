	TITLE	EGA/VGA Fast character drawing routines.
	NAME	EGACHAR
	PAGE	55,132

	COMMENT	$

	Name:	egaschar, egalchar

	Function: Draws characters from the BIOS character set
		  to the video buffer in EGA and VGA graphics modes.

		egaschar()	- 8 x 8 character cell
		egalchar()	- 8 x 14 character cell

	Caller:		Microsoft C

			void egaschar(x, y, c, n);
			void egalchar(x, y, c, n);

			int	x, y;	/* co-ords of upper left pixel of char */
			int	c;	/* the character */
			int	n;	/* its colour */


	Notes:
		The code in these routines has been optimized for speed,
		not readibility. As much as possible is put into registers.
		Replication of large blocks of code has been used in some
		places to avoid calling subroutines. 
		Taken from J.D. MacDonald's Fastega package.

		$

EGA_DATA	segment  word public 'DATA'
first8  	dw		0
first14		dw		0
EGA_DATA	ends


_CONST		segment  word public 'CONST'
fourt		dw		14
_CONST		ends


_BSS		segment  word public 'BSS'
es8		dw		?
bp8		dw		?
es14		dw		?
bp14		dw		?      
_BSS		ends

argbase		equ	2

DGROUP		GROUP	_CONST, _BSS, EGA_DATA
		assume	cs:VEGA_TEXT, ds: DGROUP, ss: DGROUP, es: DGROUP

		extrn	__buffer_segment:word

VEGA_TEXT	SEGMENT	BYTE PUBLIC	'CODE'

		public	_egaschar
_egaschar	proc	far

		push    bp
		mov     bp,sp
		push    di
		push    si
		cmp     first8,0        ;get address of bios table of characters
		jne     sec8
		inc     first8          ;done only first time this routine is called
		mov     ah, 011h
		mov     al, 030h
		mov     bh, 03
		push    bp              ;(gasp)Bios call destroys BP
		int     010h
		mov     ax,es
		mov     es8,ax
		mov     bp8,bp
		pop     bp
sec8:
		xor     ax,ax
		mov     al,byte ptr [bp+argbase+8]   ;this is the character
					;look up location of ascii characters
		mov     cl,3
		shl     ax,cl            ;in bios 
		add     ax,bp8             
		mov     si,ax              
		mov     ax,word ptr [bp+argbase+4]      ;this is x
		mov     dx,ax
		and     dx,7             ;break into bit and byte
		mov     cl,3
		sar     ax,cl
		mov     [bp+argbase+8],ax  ;byte #
		mov     [bp+argbase+4],dx  ;bit in byte     
		mov     ax, __buffer_segment       ; segment of display buffer
		mov     es,ax               
		push    ds
		mov     ax,es8          ; segment for bios reference
		mov     ds,ax           ; used to 
		mov     bx,8            ;find character
		jmp     commonl

_egaschar	endp


		public	_egalchar
_egalchar	proc	far
		push    bp
		mov     bp,sp
		push    di
		push    si

		cmp     first14,0       ;works just like slettr : done again for speed
		jne     sec14
		inc     first14
		mov     ah, 011h
		mov     al, 030h
		mov     bh, 02
		push    bp
		int     010h
		mov     ax,es
		mov     es14,ax
		mov     bp14,bp
		pop     bp
sec14:

		xor     ax,ax
		mov     al,byte ptr [bp+argbase+8]    
		mul     fourt        
		add     ax,bp14           
		mov     si,ax              
		mov     ax,word ptr [bp+argbase+4]  
		mov     dx,ax
		and     dx,7
		mov     cl,3
		sar     ax,cl
		mov     [bp+argbase+8],ax
		mov     [bp+argbase+4],dx              

		mov     ax, __buffer_segment     
		mov     es,ax               
		push    ds
		mov     ax,es14       
		mov     ds,ax      
		mov     bx,14       


commonl:
			; symbol, slettr, and llettr the same after this
			; compute location of byte to change

		mov     ax,word ptr [bp+argbase+6]         ;this is y
		mov     dx,ax
		shl     dx,1
		shl     dx,1
		add     ax,dx
		mov     cl,4
		shl     ax,cl
		add     ax,[bp+argbase+8]  ;this is horiz byte number
		mov     di,ax

		mov     dx,3ceh         ; change SET/RESET register
		mov     al,00h          ; to contain the color to write.
		out     dx,al
		mov     dx,3cfh
		mov     ax,word ptr [bp+argbase+10]     ; <-- color
		out     dx,al

		mov     dx,3ceh         ;point controller to bit mask register
		mov     al,08h
		out     dx,al
		mov     dx,3cfh
		cmp     word ptr [bp+argbase+8],0ffffh  ;special case:partially off
							;left edge of screen
		je      loop2m1
		cmp     word ptr [bp+argbase+8],79      ;Partially off right edge
		je      loop2p1
		cmp     word ptr [bp+argbase+8],80  ;totally off screen horizontally
		jae     commf
loop2:

		xor     ax,ax
		mov     ah,ds:[si+0]            ;get pattern of current row
		mov     cx,[bp+argbase+4]
		shr     ax,cl                   ;shr to get two bytes
		out     dx,al    
		inc     byte ptr es:[di+1]      ;actual screen write
		mov     al,ah
		out     dx,al
		inc     byte ptr es:[di]        ;write second byte
		add     di,80                   ;next row
		inc     si
		dec     bx
		jg      loop2
		jmp     commf
loop2m1:

		xor     ax,ax                   ;leftmost byte is off screen
		mov     ah,ds:[si+0]            ;get current row pattern
		mov     cx,[bp+argbase+4]
		shr     ax,cl                   ;shr to get two bytes
		out     dx,al    
		inc     byte ptr es:[di+1]      ;actual screen write
		add     di,80
		inc     si
		dec     bx
		jg      loop2m1
		jmp     commf
loop2p1:

		xor     ax,ax                  ;rightmost byte is off screen
		mov     ah,ds:[si+0]           ;get current row pattern
		mov     cx,[bp+argbase+4]
		shr     ax,cl                  ;shr to get two bytes
		mov     al,ah
		out     dx,al
		inc     byte ptr es:[di]
		add     di,80
		inc     si
		dec     bx
		jg      loop2p1

commf:
		pop     ds
		pop     si
		pop     di
		pop     bp
		ret

_egalchar	endp

	public	_zsetup
_zsetup	proc	far
        push    dx
        push    ax

        mov     dx,3c4h         ; allow writing to all planes
        mov     al,02h          ; can be changed by setmask
        out     dx,al
        mov     dx,3c5h
        mov     al,0fh
        out     dx,al

        mov     dx,3ceh         ; change ENABLE SET/RESET register
        mov     al,01h          ; to use SET/RESET feature
        out     dx,al
        mov     dx,3cfh
        mov     al,0fh
        out     dx,al
        pop     ax
        pop     dx
        ret
_zsetup	endp

	public	_setmask
_setmask 	proc	far
;  set mask register     setmask(mask)

        push    bp
        mov     bp,sp

        mov     cx,word ptr [bp+argbase+4]
        mov     dx,03c4h
        mov     al,2
        out     dx,al
        mov     dx,03c5h
        mov     al,cl
        out     dx,al
        pop     bp
        ret	

_setmask	endp


VEGA_TEXT	ENDS

		end

