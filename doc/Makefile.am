if CAN_BUILD_DOC
.PHONY: pdf
pdf:
	$(DOXYGEN) fullref.doxygen
	sed -i.bak -e 's/latex_count=5/latex_count=8/g' latex/Makefile
	rm -f latex/Makefile.bak
	cd latex; make
	mv -f latex/refman.pdf fullref.pdf
	$(DOXYGEN) userman.doxygen
	sed -i.bak -e 's/latex_count=5/latex_count=8/g' latex/Makefile
	rm -f latex/Makefile.bak
	cd latex; make
	mv -f latex/refman.pdf userman.pdf

.PHONY: clean
clean::
	/bin/rm -rf latex

veryclean: clean

distclean: clean

install::
	-mkdir -p $(installroot)$(prefix)/doc
	-cp -r latex $(installroot)$(prefix)/doc

install_devel::

#else
#pdf:
#	echo "doxygen is not available; cannot build the documentation"
endif
