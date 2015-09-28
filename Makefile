
# build package documentation
doc:
	/unsup/R-3.2.1/bin/R -e 'roxygen2::roxygenize()'

clean:
	rm -f *~
	rm -f */*~
	rm -f .*~	
	rm -f .Rhistory
	rm -f inst/rscripts/*~
	rm -f vignettes/.Rhistory
	rm -f vignettes/.RData
	rm -f vignettes/vignette.aux
	rm -f vignettes/vignette.log
	rm -f vignettes/vignette.out
	rm -f vignettes/vignette.toc
	rm -fr vignettes/auto

# knit the vignettes
inst/doc/%.pdf:vignettes/%.Rnw
	cd vignettes;/unsup/R-3.2.1/bin/R CMD Sweave --engine=knitr::knitr --pdf $(<F);cp -f $(<F) ../inst/doc;mv -f $(<F:.Rnw=.pdf) ../inst/doc ; rm -fr figure ; rm -f vignette.aux ; rm -f vignette.log ; rm -f vignette.out; rm -f vignette.toc ; cd ..
