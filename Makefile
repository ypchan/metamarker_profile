PREFIX ?= $(CURDIR)
BINDIR ?= $(HOME)/bin
REF_DIR ?= $(PREFIX)/refs
REF_CONFIG ?= $(PREFIX)/.meta_marker_count_ref_dir

SRC := meta_marker_count.sh
REF_BUILDER := scripts/build_ref_files.py

.PHONY: all install uninstall check print-paths configure-ref

all: install

install: configure-ref
	mkdir -p "$(BINDIR)"
	mkdir -p "$(REF_DIR)"
	ln -sfn "$(abspath $(SRC))" "$(BINDIR)/meta_marker_count"
	ln -sfn "$(abspath $(REF_BUILDER))" "$(BINDIR)/meta_marker_build_refs"
	@echo "[OK] Linked meta_marker_count:"
	@echo "  $(BINDIR)/meta_marker_count -> $(abspath $(SRC))"
	@echo "[OK] Linked meta_marker_build_refs:"
	@echo "  $(BINDIR)/meta_marker_build_refs -> $(abspath $(REF_BUILDER))"
	@echo "[OK] Reference directory:"
	@echo "  $(REF_DIR)"

configure-ref:
	mkdir -p "$(REF_DIR)"
	printf '%s\n' "$(abspath $(REF_DIR))" > "$(REF_CONFIG)"
	chmod 755 "$(SRC)"
	chmod 755 "$(REF_BUILDER)"

uninstall:
	rm -f "$(BINDIR)/meta_marker_count"
	rm -f "$(BINDIR)/meta_marker_build_refs"

check:
	bash -n "$(SRC)"
	python3 -m py_compile "$(REF_BUILDER)"
	@if command -v Rscript >/dev/null 2>&1; then \
		Rscript -e 'files <- list.files("scripts", pattern="[.]R$$", full.names=TRUE); invisible(lapply(files, parse))'; \
	else \
		echo "[WARN] Rscript not found; skip R script syntax check."; \
	fi

print-paths:
	@echo "PREFIX  = $(PREFIX)"
	@echo "BINDIR  = $(BINDIR)"
	@echo "REF_DIR = $(REF_DIR)"
	@echo "REF_CONFIG = $(REF_CONFIG)"
	@echo "SRC     = $(abspath $(SRC))"
