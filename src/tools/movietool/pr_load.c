/*
 * This file is part of
 * Movietool: display a succession of rasterfile frames in real-time
 *
 * Copyright 1989, 1990 by Ole H. Nielsen
 *
 * Author: Ole H. Nielsen
 *         Lab of Applied Physics, Bygn. 307
 *         Technical University of Denmark, DK-2800 Lyngby
 *         E-mail: ohnielse@ltf.dth.dk
 */

#include <stdio.h>
#include <pixrect/pixrect_hs.h>
 
int pr_load_header (input, rh)
FILE *input;			/* Input file */
struct rasterfile *rh;		/* Rasterfile header */
{
 
	/* Read and check file header   */
 
	if (fread(rh, sizeof(*rh), 1, input) != 1)
		return (PIX_ERR);
	if (rh->ras_magic != RAS_MAGIC)
		return (PIX_ERR);
	return (0);
}
 
int pr_load_colormap (input, rh, colormap)
FILE *input;				/* Input file */
struct rasterfile *rh;			/* Rasterfile header */
colormap_t *colormap;			/* Optional colourmap */
{
	int len;
 
	/* Load (or skip) colormap	*/
 
	if (colormap == NULL)
		fseek(input, rh->ras_maplength, 1);	/* skip map */
	else {
		colormap->type	= rh->ras_maptype;
		if (rh->ras_maptype == RMT_EQUAL_RGB) {
			/* Pixrect Reference Manual p. 62, 66:
			 * If colormap is non-NULL, colormap is read in.
			 * pr_load_colormap will allocate space for colormap
			 * if it does not match that of the file. */
			len = rh->ras_maplength;
			if (colormap->length < len / 3) {
				colormap->length = len / 3;
				if ((colormap->map[0] = (u_char *)malloc (len))
					== NULL)
					return (PIX_ERR);
				colormap->map[1] = colormap->map[0] + len / 3;
				colormap->map[2] = colormap->map[1] + len / 3;
			}
			fread(colormap->map[0], colormap->length, 1, input);
			fread(colormap->map[1], colormap->length, 1, input);
			fread(colormap->map[2], colormap->length, 1, input);
		} else
			fseek(input, rh->ras_maplength, 1); /*skip map*/
	}
	return (0);
}
 
/* The byte-encoded data format is as follows (Pixrect Ref. Manual p.59):
 *	0x80 is used as a flag to introduce run-length encoded data
 *	<0x80> <rept> <value>	represents rept+1 occurrences of <value>
 *	<0x80> <0>	 	represents a single occurrence of <0x80>
 *	<0x80> <1> <0x80> 	represents <0x80> <0x80>
 *	All other values	represent unencoded data
 * It is only worth encoding repetitions of 3 or more, except if the
 * value is <0x80> which must always be encoded.
 * Note that memory-pixrect lines must be aligned on 32-bit boundaries,
 * and that a repeated value (above) may cross line boundaries.
 */
 
#define FLAG 0x80
#define PIXRECT_NULL (struct pixrect *)NULL

struct pixrect *pr_load_image (input, rh, colormap)
register FILE *input;			/* Input file */
struct rasterfile *rh;			/* Rasterfile header */
colormap_t *colormap;			/* Optional colourmap */
{
	register char	*image;		/* Pointer to image data */
	register int	ch;		/* Input data character */
	register short	w, bw, h, count, padding, pending=0;
	struct pixrect	*pr;		/* Output pixrect */
	struct mpr_data	*data;		/* Pixrect data */
 
	/* Create pixrect and load data */
	/* We must take care because of the 32-bit line padding
	 * in a memory pixrect (See Pixrect Reference Manual p. 53) */
 
	if ((pr = mem_create (rh->ras_width,
		rh->ras_height, rh->ras_depth)) == NULL)
		return (PIXRECT_NULL);
	data = mpr_d(pr);	/* Nice for debugging */
	image = (char *) data->md_image;
	/* Byte width of a line.
	 * Note that rasterfiles contain short (16-bit) words of image data.
	 * (see /usr/include/rasterfile.h) */
	bw = ((pr->pr_depth * pr->pr_size.x + 15) / 16) * 2;
	h = (short)pr->pr_size.y;
 
	if (rh->ras_type == RT_STANDARD) {
		while (h--) {
			if (fread(image, bw, 1, input) != 1)
				return (PIXRECT_NULL);
			image += data->md_linebytes;
		}
	} else if (rh->ras_type == RT_BYTE_ENCODED) {
		padding = data->md_linebytes - bw;
		w = 0;
		while (h--) {
			w += bw;
			while (w > 0) {
				w--;
				if ((ch = getc(input)) == FLAG) {
					/* Count 0 means FLAG */
					if ((count = getc(input)) != 0) {
						ch = getc(input);
						/* count+1 of ch */
						w -= count;
						if (w < 0) {
							pending = -w;
							count -= pending;
						}
						while (count--)
							*image++ = ch;
					}
				}
				*image++ = ch;
			}
			image += padding;	/* Apply padding */
			if (pending > 0) {
				if (pending > bw) {
					count = bw;
					pending -= bw;
				} else {
					count = pending;
					pending = 0;
				}
				while (count--)
					*image++ = ch;
			}
		}
		if (ch == EOF)
			return (PIXRECT_NULL);		/* Check end-of-file */
	}

	return (pr);
}

/* Read the RT_BYTE_ENCODED image data without decoding */

struct pixrect *pr_load_encoded_image (input, rh)
FILE *input;
struct rasterfile *rh;
{
	struct pixrect *pr;

	if (rh->ras_type != RT_BYTE_ENCODED)
		return (PIXRECT_NULL);

	pr = (struct pixrect *) malloc(sizeof(struct pixrect));
	pr->pr_ops	= (struct pixrectops *)NULL;
	pr->pr_size.x	= rh->ras_width;
	pr->pr_size.y	= rh->ras_height;
	pr->pr_depth	= rh->ras_depth;
	pr->pr_data	= (caddr_t) malloc(rh->ras_length);
	/* Stuff the RT_BYTE_ENCODED image data into the pixrect data area */
	if (fread(pr->pr_data, rh->ras_length, 1, input) != 1) {
		free (pr->pr_data);
		free (pr);
		return (PIXRECT_NULL);
	}
	return (pr);
}

struct pixrect *pr_load(input, colormap)
FILE *input;				/* Input file */
colormap_t *colormap;			/* Optional colourmap */
{
	struct rasterfile rh;		/* Rasterfile header */
 
	/* Read and check file header   */
 
	if (pr_load_header (input, &rh) == PIX_ERR)
		return (PIXRECT_NULL);
 
	/* Load (or skip) colormap	*/
 
	if (pr_load_colormap (input, &rh, colormap) == PIX_ERR)
		return (PIXRECT_NULL);
 
	/* Create pixrect and load data */
 
	return (pr_load_image (input, &rh, colormap));
}
