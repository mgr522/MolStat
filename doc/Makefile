manual.pdf: documentation.config
	doxygen documentation.config
	cd manual/latex; make
	mv manual/latex/refman.pdf manual.pdf
