if CAN_BUILD_DOC

doc_DATA = userman.pdf fullref.pdf

fullref.pdf: FORCE
	$(DOXYGEN) fullref.doxygen
	cd latex; make
	mv -f latex/refman.pdf fullref.pdf
	rm -rf latex

userman.pdf: 1-intro.dox 2-math-background.dox 3-run-molstat.dox 4-electron-transport.dox
	$(DOXYGEN) userman.doxygen
	cd latex; make
	mv -f latex/refman.pdf userman.pdf
	rm -rf latex

FORCE:
	-rm -f fullref.pdf userman.pdf

CLEANFILES = warn-fullref.log warn-userman.log

endif
