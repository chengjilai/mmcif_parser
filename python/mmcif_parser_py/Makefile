# Simple automation wrappers for maturin workflows

PYTHON ?= python
VENV ?= ../../.venv
MATURIN ?= $(VENV)/bin/maturin

.PHONY: develop build clean distclean

$(MATURIN):
	@echo "maturin executable not found at $(MATURIN)"
	@echo "Call 'python3 -m venv $(VENV) && source $(VENV)/bin/activate && pip install maturin' first."
	@exit 1

develop: $(MATURIN)
	. $(VENV)/bin/activate && maturin develop --release

build: $(MATURIN)
	. $(VENV)/bin/activate && maturin build --release --sdist --strip

clean:
	rm -rf target
	rm -rf *.egg-info
	rm -rf build

# remove local venv-installed package as well
.ONESHELL:
distclean: clean
	. $(VENV)/bin/activate && pip uninstall -y mmcif_parser || true
