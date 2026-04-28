SUBDIRS = protein_digestion elemental_mass fragmentation calcisotope motif_search

all:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir || exit 1; done

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done
