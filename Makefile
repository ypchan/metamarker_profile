PREFIX ?= $(PWD)
BINDIR ?= $(PREFIX)/bin
REF_DIR ?= $(PREFIX)/refs

.PHONY: install uninstall check

install:
	install -d "$(BINDIR)"
	install -d "$(REF_DIR)"
	sed 's|^DEFAULT_REF_DIR=.*|DEFAULT_REF_DIR="$${META_MARKER_COUNT_REF_DIR:-$(REF_DIR)}"|' \
		meta_marker_count.sh > "$(BINDIR)/meta_marker_count"
	chmod 0755 "$(BINDIR)/meta_marker_count"
	install -m 0755 scripts/build_ref_files.py "$(BINDIR)/meta_marker_build_refs"

uninstall:
	rm -f "$(BINDIR)/meta_marker_count" "$(BINDIR)/meta_marker_build_refs"

check:
	bash -n meta_marker_count.sh
	python3 -m py_compile scripts/build_ref_files.py
	Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", full.names=TRUE); invisible(lapply(files, parse))'