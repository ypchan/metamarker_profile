PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin

.PHONY: install uninstall check

install:
	install -d "$(BINDIR)"
	install -m 0755 meta_marker_count.sh "$(BINDIR)/meta_marker_count"
	install -m 0755 scripts/build_ref_files.py "$(BINDIR)/meta_marker_build_refs"

uninstall:
	rm -f "$(BINDIR)/meta_marker_count" "$(BINDIR)/meta_marker_build_refs"

check:
	bash -n meta_marker_count.sh
	python3 -m py_compile scripts/build_ref_files.py
	Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", full.names=TRUE); invisible(lapply(files, parse))'
