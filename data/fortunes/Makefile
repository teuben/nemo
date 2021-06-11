
#   Maintain the fortunes database

FORTUNES = /usr/share/games/fortunes

nemo.dat:    nemo
	strfile -c% nemo nemo.dat

#   if you have root permission, you could do this too:
#   this way "fortune nemo"  would work
#   else
#            "fortune $NEMO/data/fortunes"
#   is needed.

root:
	sudo cp nemo nemo.dat $(FORTUNES)
