SHELL := /bin/bash
CONDA_RUN ?= conda run -n cryoem-book
REFERENCE_DIR ?= /home/phonchi/cryo-em_general/References
CRYOEM_OUTPUT_DIR ?= /tmp/cryoem-primer-output
JUPYTER_DATA_DIR ?= /tmp/cryoem-primer-jupyter

.PHONY: ingest sync test qa-fast build-full linkcheck clean

ingest:
	$(CONDA_RUN) python scripts/reference_pipeline.py ingest --source-root "$(REFERENCE_DIR)"
	$(CONDA_RUN) python scripts/reference_pipeline.py validate --source-root "$(REFERENCE_DIR)"

sync:
	JUPYTER_DATA_DIR="$(JUPYTER_DATA_DIR)" $(CONDA_RUN) jupytext --sync book/*.py

test:
	$(CONDA_RUN) pytest -q

qa-fast: sync
	$(CONDA_RUN) pytest -q
	CRYOEM_NUM_IMAGES=100 CRYOEM_OUTPUT_DIR="$(CRYOEM_OUTPUT_DIR)/fast" \
		$(CONDA_RUN) jupyter-book build book --all -W --keep-going
	$(CONDA_RUN) python scripts/finalize_html.py book/_build/html
	$(CONDA_RUN) python scripts/check_site.py book/_build/html
	git diff --exit-code -- book/*.ipynb

build-full: sync
	$(CONDA_RUN) pytest -q
	CRYOEM_NUM_IMAGES=5000 CRYOEM_OUTPUT_DIR="$(CRYOEM_OUTPUT_DIR)/full" \
		$(CONDA_RUN) jupyter-book build book --all -W --keep-going
	$(CONDA_RUN) python scripts/finalize_html.py book/_build/html
	$(CONDA_RUN) python scripts/check_site.py book/_build/html

linkcheck:
	$(CONDA_RUN) jupyter-book build book --builder linkcheck --all --keep-going

clean:
	$(CONDA_RUN) jupyter-book clean book --all
