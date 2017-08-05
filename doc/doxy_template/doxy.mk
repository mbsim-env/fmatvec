# this is the default doxy_template dir. We only used another dir for fmatvec
# since for fmatvec it may not installed yet
doxy_template_dir ?= $(prefix)/share/fmatvec/doxy_template

doxytempl.all: doxyfile
	rm -rf html
	$(doxygen) doxyfile
	cp $(doxy_template_dir)/doxy-boot.js html/.
doxytempl.clean:
	rm -rf html
doxytempl.install:
	mkdir -p $(datadir)/doc/$(PACKAGE)
	cp -r html/* $(datadir)/doc/$(PACKAGE)/.
doxytempl.uninstall:
	rm -rf $(datadir)/doc/$(PACKAGE)
