#

URL = apus:/n/apus1/nemo/local/www/mass99

FILES = index.html plot1.gif plot2.gif plot3.gif plot4.gif


help:
	@echo scp will copy to your website $(URL)

scp:	.last_scp

.last_scp:	$(FILES)
	scp $? $(URL)
	touch .last_scp
