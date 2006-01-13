DEV_ROOT=.

-include make/Makefile

make/Makefile : make/Makefile.in configure
	@if [ "x$(GNUMAKE)" == "x" ]; then \
		@echo ;\
		@echo 'Building OOPSE requires the use of GNU Make!'; \
		@echo ;\
		@exit 1;\
	fi
	@echo
	@echo 'Please run (or re-run) configure'
	@echo
	@exit 1


