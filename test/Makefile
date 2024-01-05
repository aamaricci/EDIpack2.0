TSTDIR=src
TSTBIN=bin

all: bin
	make -C $(TSTDIR)

bin:
	mkdir $(TSTBIN)
####### BATH_TYPE FOLDERS
	mkdir $(TSTBIN)/NORMAL
#	mkdir $(TSTBIN)/HYBRID
	mkdir $(TSTBIN)/REPLICA
	mkdir $(TSTBIN)/GENERAL
####### ED_MODE FOLDERS
	mkdir $(TSTBIN)/NORMAL/NORMAL
	mkdir $(TSTBIN)/NORMAL/SUPERC
#	mkdir $(TSTBIN)/HYBRID/NORMAL
#	mkdir $(TSTBIN)/HYBRID/SUPERC
	mkdir $(TSTBIN)/REPLICA/NORMAL
	mkdir $(TSTBIN)/REPLICA/SUPERC
	mkdir $(TSTBIN)/GENERAL/NORMAL
	mkdir $(TSTBIN)/GENERAL/SUPERC
####### INSERT FOLDER TO BE TESTED IN LIST_DIR
	echo NORMAL/NORMAL >> $(TSTBIN)/list_dir
	echo NORMAL/SUPERC >> $(TSTBIN)/list_dir
#	 echo HYBRID/NORMAL  >> $(TSTBIN)/list_dir
#	 echo HYBRID/SUPERC  >> $(TSTBIN)/list_dir
	echo REPLICA/NORMAL >> $(TSTBIN)/list_dir
	echo REPLICA/SUPERC >> $(TSTBIN)/list_dir
	echo GENERAL/NORMAL >> $(TSTBIN)/list_dir
	echo GENERAL/SUPERC >> $(TSTBIN)/list_dir

test:
	./test.sh

clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
