rem balrk
rem  Usage is "mkfnts destination", eg. "mkfnts \lib\hershey"
rem  (Make sure that the directory exists first though)
rem

h2v ..\data\hersh.oc
h2v ..\data\hersh.or ..\fonts\japan.hmp japanese
copy astrolog %1
copy cursive %1
copy cyrillic %1
copy futura.l %1
copy futura.m %1
copy gothic.eng %1
copy gothic.ger %1
copy gothic.ita %1
copy greek %1
copy markers %1
copy math.low %1
copy math.upp %1
copy meteorol %1
copy music %1
copy script %1
copy symbolic %1
copy times.g %1
copy times.i %1
copy times.ib %1
copy times.r %1
copy times.rb %1
copy japanese %1
del astrolog
del cursive
del cyrillic
del futura.l
del futura.m
del gothic.eng
del gothic.ger
del gothic.ita
del greek
del markers
del math.low
del math.upp
del meteorol
del music
del script
del symbolic
del times.g
del times.i
del times.ib
del times.r
del times.rb
del japanese
