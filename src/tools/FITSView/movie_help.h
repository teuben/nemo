" \n\
    Right Mouse Button:  click in image frame to run movie\n\
    Left Mouse Button:   click in image frame to select previous image \n\
    Middle Mouse Button: click in image frame to select next image \n\
\n\
    Frame Speed:    Controls movie speed. Slider indicates number\n\
                    of tenths of seconds between each frame \n\
                    (e.g., a slider value of 5, will cause a pause \n\
                    of 0.5 sec between successive frames). Due to\n\
                    system loads, your actual mileage may vary.\n\
\n\
    Movie Cycles:   Controls number of movie cycles. The movie will be \n\
		    shown the number of times indicated by the slider.\n\
\n\
    Min/Max Action: Set whether changes in Min and Max will affect all\n\
                    images or just the current image. This is just a \n\
	            toggle with no direct effect on the images--use Max\n\
	            to reload images after changing this setting.\n\
\n\
\n\
    Min:    set minimum display value for image\n\
    Max:    set maximum display value for image and reload image with\n\
            new min and max values.\n\
\n\
    Gamma:  set display correction factor for non-linearity\n\
            of monitor. [Default = .37]\n\
            Correction to greyscale ramp is of the form:\n\
\n\
               ramp[i] = 255*{(i/255)^(gamma)}   (i = 0 to 255)\n\
\n\
            where a perfect linear monitor would need no correction and\n\
            have a greyscale ramp:\n\
\n\
               ramp[i] = i; (i = 0 to 255)\n\
\n\
\n\
    Grey-Scale/Pseudo-Color/3-Color: \n\
            Select either a grey-scale, pseudocolor or 3-color movie.\n\
\n\
            Choosing Pseudo-color allows the user to input the\n\
             number of color contours to divide the 256 brightness levels\n\
             of the display into. \n\
            Choosing 3-color will cause the (n-1)th frame to be displayed in\n\
             blue, the nth frame in green and the (n+1)th in red, where n is the\n\
             \"current\" \n\ frame.  On 24-bit displays, this is a full-color\n\
             representation; on 8-bit displays a dithering scheme is used to\n\
             mimic 24-bits of color (see documentation for 'rgb')."
