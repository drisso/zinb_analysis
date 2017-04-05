NCELLS=10000
RATIO=5
OFFS=5
DS=Allen
#RATIO=1 5 50
#OFFS=0 2 5 
#DS=Zeisel
NCORES=8

targets := $(foreach r,$(OFFS),$(patsubst %,%-offs$r,$(RATIO)))

all: $(targets)
    
$(targets):
	@echo "#!/bin/bash" > $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#$$ -cwd" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#$$ -V" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#$$ -j y" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#$$ -S /bin/bash" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "#" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	@echo "mpirun -n 1 Rscript --vanilla fitZinb10000.R $(DS) $(NCELLS) $@" >> $(DS)_nc$(NCELLS)_ratio$@.sh
	qsub -pe orte $(NCORES) $(DS)_nc$(NCELLS)_ratio$@.sh
