
This directory contains the public domain Hershey font data, distributed on
the USENET network. 

The Hershey Fonts were originally created by Dr. A. V. Hershey while
working at the U. S. National Bureau of Standards.

The format of the Font data in this distribution
was originally created by:
			James Hurt
			Cognition, Inc.
			900 Technology Park Drive
			Billerica, MA 01821
			(mit-eddie!ci-dandelion!hurt)

The data is stored in the file hersh.oc, and the file hersh.or.  Each glyph
in the font file is stored as a record.  For example, glyph 2 has the entry:

    2 16MWOMOV ROMSMUNUPSQ ROQSQURUUSVOV

The first 5 characters are the glyph number (2), and the following 3 are the
number of data pairs (16).  The remaining characters are (x,y) pairs, each
character representing an x (or y) motion.  The characters are encoded
relative to "R" which is zero (Q = -1, S = 1, etc).  The first pair is a
left and right offset for proportional spacing.  The second pair is a move.
The program takes a space in the x coordinate to mean that the next pen
motion is a move.  Thus the line above encodes:

<glyph 2> <#pairs 16> 
<spacing -5, 5>	(MW)
(skip)  (-3,-5)	(OM)
(draw)  (-3, 4)	(OV)
<skip>  (-3,-5)	(OM)
<draw>  ( 1,-5)	(SM)
<draw>  ( 3,-4)	(UN)
<draw>  ( 3,-2)	(UP)
<draw>  ( 1,-1)	(SQ)
<skip>  (-3,-1)	(OQ)
<draw>  ( 1,-1)	(SQ)
<draw>  ( 3, 0)	(UR)
<draw>  ( 3, 3)	(UU)
<draw>  ( 1, 4)	(SV)
<draw>  (-3, 4)	(OV)

Note: some records contain new-lines whch should be ignored.

The Occidental Characters.

The file hersh.oc files contain approximately 19 different fonts in the
A-Z alphabet plus greek and cyrillic, along with hundreds of special symbols,
described generically below.

Symbols:
        Mathematical (227-229,232,727-779,732,737-740,1227-1270,2227-2270,
                        1294-1412,2294-2295,2401-2412)
        Daggers (for footnotes, etc) (1276-1279, 2276-2279)
        Astronomical (1281-1293,2281-2293)
        Astrological (2301-2312)
        Musical (2317-2382)
        Typesetting (ffl,fl,fi sorts of things) (miscellaneous places)
        Miscellaneous (mostly in 741-909, but also elsewhere):
                - Playing card suits
                - Meteorology
                - Graphics (lines, curves)
                - Electrical
                - Geometric (shapes)
                - Cartographic
                - Naval
                - Agricultural
                - Highways
                - Etc...

The Oriental Characters.

	The file hersh.or contains over 1000 Kanji, or Chinese characters,
and the Japanese Hiragana and Katakana alphabets.

