all: manual.pdf fitter-v1.0.tgz

manual.pdf: documentation.config
	doxygen documentation.config
	cd manual/latex; make
	cp manual/latex/refman.pdf fitter/doc/manual.pdf

fitter-v1.0.tgz:
	COPYFILE_DISABLE=true tar -czf fitter-v1.0.tgz --exclude="intro.dox" fitter/
