

default:
	Rscript  -e 'library(devtools); build(); document(); '


check:
	R CMD check .


pdf:
	R CMD Rd2pdf .

install:
	cd .. && R CMD INSTALL BF

package:
	cd .. && R CMD build BF

