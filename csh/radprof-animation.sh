#! /bin/bash 
#
#   create animated gif from a set of PNG's
#   Warning: seems to run out of memory after about 290 images in full size,
#   for 512x512 size all 400 in my example worked fine.


file=p10ks.out
times=$(snaptime $file)
fuzz=0.01
xmin=-3
xmax=1
ymin=-6
ymax=6

yapp() {
    num=$(printf %04d  $1)
    if test $yapp = "xs"; then
        echo $1/$yapp
    elif test $yapp = "_ps"; then
        echo fig$num.ps
    else
        echo fig$num.$yapp/$yapp
    fi
}

export PGPLOT_GIF_WIDTH=512
export PGPLOT_GIF_HEIGHT=512
yapp=gif
#yapp=png
n=0
for t in $times; do
    snaptrim $file - times=$t timefuzz=$fuzz |\
	snapdens - - |\
	snapprint  - r,aux |\
	tabmath - - 'log(%1),log(%2)' all |\
	tabplot - xmin=$xmin xmax=$xmax ymin=$ymin ymax=$ymax headline=$t yapp=$(yapp $n) 
    #snaptrim $file - times=$t timefuzz=$fuzz | hackdens - - write_at_phi=t | snapprint  - r,phi | tabmath - - 'log(%1),log(%2)' all | tabplot -
    n=$(expr $n + 1)
    o=$(yapp $n)
    echo yapp=$o n=$n
done


echo 'Possible commands to make an animation:'
echo 'convert -delay 20 -loop 0 fig0*png test.gif'
echo 'ffmpeg -r 60 -f image2 -s 1920x1080 -i fig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4'
