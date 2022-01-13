# GNU Automake file

MODEL = cca
MODEL_URL=https://g-c662a6.a78b8.36fe.data.globus.org/ucvm/v21_10/model
MODEL_TARBALL=cca.tar.gz
MODEL_DIR=$(datadir)/model/$(MODEL)

dist_noinst_DATA = config
dist_noinst_SCRIPTS = make_data_files.py

install-data-local: installed_data

installed_data:
	if [ ! -f "$(MODEL_TARBALL)" ]; then $(CURL) -O $(MODEL_URL)/$(MODEL_TARBALL); fi
	$(MKDIR_P) $(MODEL_DIR) && \
		cd $(MODEL_DIR) && \
		$(TAR) --strip 3 -xf $(abs_builddir)/$(MODEL_TARBALL) CCA/data/config CCA/data/i06/density.dat CCA/data/i06/vp.dat CCA/data/i06/vs.dat
	touch $(abs_builddir)/installed_data


# End of file
