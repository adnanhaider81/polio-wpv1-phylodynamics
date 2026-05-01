R ?= Rscript

.PHONY: install install-optional synthetic demo demo-fast clean

install:
	$(R) scripts/install_packages.R

install-optional:
	$(R) scripts/install_packages.R --optional

synthetic:
	$(R) scripts/make_synthetic_data.R

demo:
	$(R) run_pipeline.R --root . --mode distance --nsim 5 --seed 1

demo-fast:
	$(R) run_pipeline.R --root . --mode distance --nsim 3 --seed 1 --skip_figures

clean:
	rm -rf data/processed/* results/figures/* results/tables/* results/trees/* logs/* tmp/*
