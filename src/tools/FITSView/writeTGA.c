/* write a TRUEVISION TGA file, format from Truevision TGA V2.0 */
/* source image is a SunView Pixrect */

#include "ql.h"
#include <stdio.h>

/* Definitions for Image Types */
#define NO_IMAGE_DATA   0x00
#define TGA_MapRGB      0x01
#define TGA_RawRGB      0x02
#define TGA_RawMono     0x03
#define TGA_MapEnCode   0x09
#define TGA_RawEnCode   0xa

#define TGA_MonoEnCode  0xb
/* TGA File header definition */
        struct TGA_ImageHeader {
                unsigned char IDLength;  /* Length of Identifier String, max=255 */
                unsigned char CoMapType; /* 0 = No ColorMap included in image    */
                                         /* 1 = ColorMap included in image       */
                unsigned char ImgType;   /* Image Type (0,1,2,3,9,10)            */
                                         /* 0 No image data
                                          * 1 Uncompressed Color-Mapped Image
                                          * 2 Uncompressed True-Color Image
                                          * 3 Uncompressed Black&White Image
                                          * 9 Run-Length encoded CMI
                                          * 10 R-L Encoded TCI
                                          * 11 R-L Encoded B&WI
                                          */
/* for some reason making these shorts messes things up */
/*              unsigned short Index;  /* Index of first color map entry */
/*              unsigned short Length; /* Length of color map (number of entries)*/
                unsigned char Index_h;  /* Index of first color map entry */
                unsigned char Index_l;  /* Index of first color map entry */
                unsigned char Length_l; /* Length of color map (number of entries)*/                unsigned char Length_h; /* Length of color map (number of entries)*/                unsigned char CoSize;   /* bits/color map entry */
                unsigned short X_org;   /* X Origin of Image     */
                unsigned short Y_org;   /* Y Origin of Image     */
                unsigned short Width;   /* Width in pixels       */
                unsigned short Height;  /* Height in pixels      */
                unsigned char PixDepth; /* Bits/Pixel (8,16,24,32) */
/* the following 8 bits make up the Image Descriptor */
                unsigned AttBits : 4;   /* Number of Attribute Bits per pixel */
                unsigned OrgBit  : 2; /* Origin Bits:
                                        00 bottom left
                                        01 bottom right
                                        10 top left
                                        11 top right
                                       */
                unsigned Rsrvd   : 2; /* Reserved bits */
                };

/* TGA file footer definition as per 1989 standard */
struct TGA_FileFooter {
	unsigned long ExtAreaOffset;
	unsigned long DevAreaOffset;
	unsigned char Signature[16];
	unsigned char period; 
	unsigned char terminator;
} ;
struct TGA_FileFooter tgaf = {0, 0, "TRUEVISION-XFILE",'.',0x00};
struct TGA_ImageHeader  *tga;

writeTGA() 
{
	char TGA_ImageIDField[255];
	FILE *fout; 
	Pixrect *pr;
	register unsigned char *start, *image;
	int k, j, i, linebytes, nb;
	short xs,ys, swap;
/* must swap bytes for motorola -> intel, TGA uses intel byte order 
/* swab_(from,to,nbytes) nbytes even */
/* C library swab() is returning -2 for some reason */

	if(!is24Bit) {
		fprintf(stderr,"writeTGA  only works on 24bit displays\n");
		return(-1);
	}

/* Input the Targa file header */
        if ((tga=(struct TGA_ImageHeader *)
                 malloc(sizeof(struct TGA_ImageHeader))) == NULL)
        {
                fprintf(stderr,"Can't allocate TGA memory\n");
                return(-1);
        }

/* assign values to image header*/
	if(is24Bit) {
		pr = pic24.pr;
		j = sizeof(pic24.label);
		tga->IDLength = sizeof(pic24.label);
	}
        else {
		pr = pic[0].pr;
		j = sizeof(pic[0].label);
		tga->IDLength = j; 
	}
	tga->CoMapType = 0;
	tga->ImgType = TGA_RawRGB;

	tga->Index_l  = 0x00;
	tga->Index_h  = 0x00;
	tga->Length_l = 0x00;
	tga->Length_h = 0x00;
	tga->CoSize = 0;
/* no need to swap these if they are zero */
	tga->X_org  = (short) 0; 
	tga->Y_org  = (short) 0;
	xs = (short) pr->pr_size.x;
	nb=2;
	swab_(&xs,&swap,&nb);
	tga->Width=swap;
	ys = (short) pr->pr_size.y;
	swab_(&ys,&swap,&nb);
	tga->Height=swap;
	xs = (short) pr->pr_size.x;
	ys = (short) pr->pr_size.y;
/*fprintf(stderr,"xs = %d\n",xs);
fprintf(stderr,"ys = %d\n",ys);
fprintf(stderr,"W %d H %d\n",tga->Width, tga->Height);*/
/* "24-bit" pixrects are 32-bit, but we will throw away the top byte */
/*	tga->PixDepth = (is24Bit ? 24 : 8); */
	tga->PixDepth = 0x18;
	tga->AttBits=0;
	tga->OrgBit=0; 
	tga->Rsrvd=0;
/* don't swap character strings */
	if(is24Bit) 
		strcpy(TGA_ImageIDField,pic24.label);
        else {
		strcpy(TGA_ImageIDField,pic[0].label);
	}

	if( (fout = fopen("tga.out","w")) == NULL) {
		fprintf(stderr,"error opening file ./tga.out\n");
		return(-1);
	}
/* write out header */
	fwrite(tga,1,18,fout); 
	fwrite(TGA_ImageIDField,tga->IDLength,1,fout);

/* write out image, brute force */
	start = (unsigned char *) mpr_d(pr)->md_image;
        linebytes = mpr_d(pr)->md_linebytes;
/* i counts rows, j counts pixels/row, k counts bytes */
/* image is a pointer to a char[4]                    */
/* linebytes is bytes/line, NOT pixels/line           */
        for(i=ys-1;i>-1;i--) {
                for(j=0,k=0;j<xs;k+=4,j++) {
                        image=start+k+i*linebytes;
/*fprintf(stdout,"%d %d %d %d\n",image[0],image[1],image[2], image[3]);*/
/* NOTE: For some reason, GAD swapped the bytes on our
   first picture, (ie BGR instead of RGB) even though 
   this is definitely the right TGA format. TO quickly get
   around this problem, write out the bytes in the reverse
   order and figure out what GAD is doing wrong later.
	mwp 7/23/90

			putc(image[3],fout);
			putc(image[2],fout);
			putc(image[1],fout);
*/
			putc(image[1],fout);
			putc(image[2],fout);
			putc(image[3],fout);
		}	
	}
/* write out footer */
	fwrite(tgaf,1,26,fout);
	free(tga);
	fclose(fout);
	return(0);
}
swab_(src, dst, nbyt)
register char *src, *dst;
int *nbyt;
{
	register int nb = *nbyt >> 1;
	register int c;
	while(--nb >= 0) {
		c = *src++;
		*dst++ = *src++;
		*dst++ = c;
	}
}
getbits(x,p,n) /* get n bits from postion p ; 0 at right */
unsigned x, p ,n;
{
	return((x>> (p+1-n)) & ~(~0 << n));
}
