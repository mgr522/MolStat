if BUILD_DOC

doc_DATA = userman.pdf fullref.pdf

.PHONY: FORCE

userman.pdf: FORCE
	$(DOXYGEN) userman.doxygen
	cd latex; make
	mv -f latex/refman.pdf userman.pdf
	rm -rf latex

fullref.pdf: FORCE
	$(DOXYGEN) fullref.doxygen
	cd latex; make
	mv -f latex/refman.pdf fullref.pdf
	rm -rf latex

FORCE:
	-rm -f fullref.pdf userman.pdf

CLEANFILES = warn-fullref.log warn-userman.log

endif
