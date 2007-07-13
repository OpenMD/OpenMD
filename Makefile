DEV_ROOT=.

RootFiles =  \
	AUTHORS \
	LICENSE \
	README \
	ChangeLog

-include make/Makefile

make/Makefile : make/Makefile.in configure
	@echo
	@echo 'Please run (or re-run) configure'
	@echo
	@exit 1


