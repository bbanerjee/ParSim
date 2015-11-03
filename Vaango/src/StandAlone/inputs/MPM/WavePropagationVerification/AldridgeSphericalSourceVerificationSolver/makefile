RELEASE_FILES = explod_for.f  explod_for.in  heaviside.in  makefile  README.txt  INFO.txt  EXPLODreleaseInstructions.txt aldridge_sphere_source.pdf explod_for.in-example explod_for.out-example heaviside.in-example



run: explod_for.exe explod_for.in heaviside.in
	touch explod_for.out
	rm -f fort.1 fort.2 fort.3
	ln -s explod_for.in    fort.1
	ln -s heaviside.in     fort.2
	ln -s explod_for.out   fort.3
	./explod_for.exe
	rm -f fort.1 fort.2 fort.3



explod_for.exe: explod_for.f
	f77 explod_for.f -o explod_for.exe



wipe:
	rm -f fort.[123] a.out core.* explod_for.exe



release: $(RELEASE_FILES)
	@rm -rf EXPLODrelease
	@rm -f EXPLODrelease.tar
	@mkdir EXPLODrelease
	@cp $(RELEASE_FILES) EXPLODrelease
	@unix2dos EXPLODrelease/README.txt
	@unix2dos EXPLODrelease/INFO.txt
	@unix2dos EXPLODrelease/EXPLODreleaseInstructions.txt
	@tar -cvf EXPLODrelease.tar EXPLODrelease
	@mv EXPLODrelease.tar EXPLODrelease
	@echo ''
	@echo ''
	@echo ''
	@echo '#############################################################'
	@echo 'The release files have been written to directory EXPLODrelease'
	@echo 'That directory contains a file called EXPLODrelease.tar, which'
	@echo 'can be sent to colleagues who wish to use this code.'
	@echo 'Simply send them EXPLODrelease.tar'
	@echo 'along with the following instructions:'
	@cat EXPLODreleaseInstructions.txt
