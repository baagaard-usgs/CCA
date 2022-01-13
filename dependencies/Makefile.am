# GNU Automake config

DEPS_URL=https://g-c662a6.a78b8.36fe.data.globus.org/ucvm/v21_10/lib
PROJ_VER=5.0.0
EUCLID_VER=3.1.4


all: installed_etree installed_proj

# ------------------------------------------------------------------------------
# ETREE/EUCLID
# ------------------------------------------------------------------------------
etree:
if INSTALL_ETREE
EUCLID_URL=https://github.com/baagaard-usgs/euclid-etree/releases/download/v$(EUCLID_VER)
EUCLID_TARBALL=euclid-$(EUCLID_VER).tar.gz
	if [ ! -f "$(EUCLID_TARBALL)" ]; then $(CURL) -L -O $(EUCLID_URL)/$(EUCLID_TARBALL); fi
	$(TAR) -xf $(EUCLID_TARBALL)
	$(MKDIR_P) euclid-build
	cd euclid-build && \
		../euclid-$(EUCLID_VER)/configure \
			--prefix=$(prefix) \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make -C libsrc && make install -C libsrc && \
		touch $(abs_builddir)/installed_$@
else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_etree:
	$(MAKE) $(AM_MAKEFLAGS) etree

# ------------------------------------------------------------------------------
# PROJ
# ------------------------------------------------------------------------------
proj:
if INSTALL_PROJ
PROJ_TARBALL=proj-$(PROJ_VER).tar.gz
	if [ ! -f "$(PROJ_TARBALL)" ]; then $(CURL) -O $(DEPS_URL)/$(PROJ_TARBALL); fi
	$(TAR) -xf $(PROJ_TARBALL)
	$(MKDIR_P) proj-build
	cd proj-build && \
		../proj-$(PROJ_VER)/configure --with-jni=no \
			--prefix=$(prefix) \
			CC=$(CC) CPPFLAGS="$(CPPFLAGS)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"; \
		make && make install && \
		touch $(abs_builddir)/installed_$@
else
	@echo "Using existing $@ installation."
	touch $(abs_builddir)/installed_$@
endif

installed_proj:
	$(MAKE) $(AM_MAKEFLAGS) proj


# End of file
