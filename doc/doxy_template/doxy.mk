# this is the default doxy_template dir. We only used another dir for fmatvec
# since for fmatvec it may not installed yet
fmatvec_share_dir ?= $(prefix)/share/fmatvec

doxytempl.all: doxyfile
	rm -rf html
	$(doxygen) doxyfile
	cp $(fmatvec_share_dir)/doxy_template/doxy-boot.js html/.
	$(fmatvec_share_dir)/checkHtml/checkMathJax.sh html
doxytempl.clean:
	rm -rf html
doxytempl.install:
	mkdir -p $(datadir)/doc/$(PACKAGE)
	cp -r html/* $(datadir)/doc/$(PACKAGE)/.
doxytempl.uninstall:
	rm -rf $(datadir)/doc/$(PACKAGE)
