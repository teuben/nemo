#	
# June 1995: the (Solaris?) windex file is slightly more sensitive
#            to the one-liner requirement of the .SH NAME entry
#	     Try:     grep '\\(em' * | awk -F: '{print $1}'
#            to find man entries that won't make it into the windex
#            and solaris-man won't subsequently find.
help:
	@echo "print	print all manual pages since file .lastprint:"
	@echo "            `ls -l .lastprint`"
	@echo "newman	show all new manual pages"
	@echo "catman	force new 'whatis' or 'windex' file"
	@echo "fmtman	force new man/catX files"	
	@echo "clean    clean out catX files"

print:	
	@find . -type f -newer .lastprint -name '*.[1358]' -exec manlaser '{}' \;
	@touch .lastprint

newman:
	@find . -type f -newer .lastprint -name '*.[1358]' -exec ls -l '{}' \;

bad:
	@echo Manual pages that do not make it into the catman index 
	grep '\\(em' man*/*

clean:
	rm -f cat?/*

CATMANFLAGS = 
SECTIONS=13568l

# install: easier target for updates
install:
	@make fmtman
	mkdir doc

# make new whatis
whatis: .whatis
windex: .windex
#
.whatis: catman
.windex: catman

catman:
	@echo Making new whatis file, takes few seconds 
	@echo "Use 'make fmtman' to also make new catX files"
	catman $(CATMANFLAGS) -w -M $(NEMO)/man -T $(NEMO)/man/tmac.an \
		$(SECTIONS)

# make new cat files
fmtman:
	@echo Making new cat files, please be patient...
	catman $(CATMANFLAGS) -M $(NEMO)/man -T $(NEMO)/man/tmac.an \
		$(SECTIONS)

# show what the new ones are:
new:
	@catman $(CATMANFLAGS) -p -M $(NEMO)/man -T $(NEMO)/man/tmac.an \
		$(SECTIONS)

fortune: nemo-fortunes.dat
	fortune nemo-fortunes

nemo-fortunes.dat: nemo-fortunes
	strfile nemo-fortunes nemo-fortunes.dat

#    Reminder how to make a PS file from a man page

MAN = tsf.1
ps: $(MAN)
	groff -man -Tps $(MAN) > man.ps

#    Peter's update for the website
chara:
	./make_html chara.astro.umd.edu:/chara/teuben/nemo-www/man_html/

#    updates the docs tree for https://teuben.github.io/nemo/
docs:
	./make_html ../docs/man_html/
