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
#include <memory.h>
#include <pixrect/pixrect_hs.h>

/* The byte-encoded data format is as follows (Pixrect Ref. Manual p.59):
 *	0x80 is used as a flag to introduce run-length encoded data
 *	<0x80> <rept> <value>	represents rept+1 occurrences of <value>
 *	<0x80> <0>	 	represents a single occurrence of <0x80>
 *	<0x80> <1> <0x80> 	represents <0x80> <0x80>
 *	All other values	represent unencoded data
 * It is only worth encoding repetitions of 3 or more, except if the
 * value is <0x80> which must always be encoded.
 * Note that memory-pixrect lines must be aligned on 32-bit boundaries,
 * and that a repeated value (above) may cross one or more line boundaries.
 */

#define FLAG 0x80

int pr_decode (dpr, xoff, yoff, spr)
struct pixrect *dpr, *spr;
int xoff, yoff;		/* Offsets */
{
	register u_char *buf;		/* Input data buffer */
	register u_char *image;		/* Pointer to image data */
	register u_char ch;
	register short w, count;
	register short bw, h, padding, pending = 0;
	short line;
	u_char *lastbyte;
	struct mpr_data *data;
 
	if (spr->pr_depth != dpr->pr_depth)
		return (PIX_ERR);

	/* Input buffer of RT_BYTE_ENCODED data (My private definition) */
	buf = (u_char *) spr->pr_data;
	data = mpr_d(dpr);
	image = mprd8_addr (data, xoff, yoff, dpr->pr_depth);
	lastbyte = mprd8_addr(data, dpr->pr_size.x, dpr->pr_size.y,
		dpr->pr_depth);
	/* Byte width of a line.
	 * Note that rasterfiles contain short (16-bit) words of image data.
	 * (see /usr/include/rasterfile.h) */
	bw = ((spr->pr_depth * spr->pr_size.x + 15) / 16) * 2;
	h = (short)spr->pr_size.y;
 
	padding = data->md_linebytes - bw;
	w = 0;

	/* Big loop over rasterfile lines */

	for (line=0; line<h; line++) {
		if (pending > 0) {
			if (pending > bw) {
				count = bw;
				pending -= bw;
			} else {
				count = pending;
				pending = 0;
			}
			(void) memset (image, *(buf-1), count);
			image += count;
		}
		w += bw;

		while (w > 0) {
			w--;
			if (*buf == FLAG) {
				/* Count 0 means FLAG */
				if ((count = *(++buf)) != 0) {
					buf++;
					/* count+1 of ch */
					w -= count;
					if (w < 0) {
						pending = - w;
						count -= pending;
					}
					count++;
					(void) memset (image, *buf, count);
					image += count;
				} else
					*image++ = FLAG;
			} else
				*image++ = *buf;
			buf++;
		}
		image += padding;		/* Apply padding */
		if (image >= lastbyte)		/* Clip at end of dpr */
			break;
	}
	return (0);
}

int pr_decode_and_zoom (dpr, xoff, yoff, spr, zoom)
struct pixrect *dpr, *spr;
int xoff, yoff;		/* Offsets */
int zoom;		/* Zoom factor */
{
	register u_char *buf;		/* Input data buffer */
	register u_char *image;		/* Pointer to image data */
	register u_char ch;
	register short w, count;
	register short bw, h, padding, pending = 0;
	short line;
	u_char *lastbyte, *line_start;
	struct mpr_data *data;
 
	if (zoom < 2 || spr->pr_depth != dpr->pr_depth)
		return (PIX_ERR);

	/* Input buffer of RT_BYTE_ENCODED data (My private definition) */
	buf = (u_char *) spr->pr_data;
	data = mpr_d(dpr);
	image = mprd8_addr (data, xoff, yoff, dpr->pr_depth);
	lastbyte = mprd8_addr(data, dpr->pr_size.x, dpr->pr_size.y,
		dpr->pr_depth);
	/* Byte width of a line.
	 * Note that rasterfiles contain short (16-bit) words of image data.
	 * (see /usr/include/rasterfile.h) */
	bw = zoom * ((spr->pr_depth * spr->pr_size.x + 15) / 16) * 2;
	h = (short)spr->pr_size.y;
 
	padding = data->md_linebytes - bw;
	w = 0;

	/* Big loop over rasterfile lines */

	for (line=0; line<h; line++) {
		line_start = image;
		if (pending > 0) {
			if (pending > bw) {
				count = bw;
				pending -= bw;
			} else {
				count = pending;
				pending = 0;
			}
			(void) memset (image, *(buf-1), count);
			image += count;
		}
		w += bw;

		while (w > 0) {
			w -= zoom;
			if (*buf == FLAG) {
				/* Count 0 means FLAG */
				if ((count = *(++buf)) != 0) {
					buf++;
					/* count+1 of ch */
					count *= zoom;
					w -= count;
					if (w < 0) {
						pending = - w;
						count -= pending;
					}
					count += zoom;
					(void) memset (image, *buf, count);
					image += count;
				} else {
					(void) memset (image, FLAG, zoom);
					image += zoom;
				}
			} else {
				(void) memset (image, *buf, zoom);
				image += zoom;
			}
			buf++;
		}
		image += padding;		/* Apply padding */
		if (image >= lastbyte)		/* Clip at end of dpr */
			break;

		/* Zoom in the y-direction */
		for (count = 1; count < zoom; count++) {
			memcpy (image, line_start, bw);
			image += data->md_linebytes;
			if (image >= lastbyte) break;
		}
		if (image >= lastbyte) break;
	}
	return (0);
}
