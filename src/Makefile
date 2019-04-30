scpc : main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o write_species.o
	cc -lm -o scpc main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o write_species.o

main.o : main.c
	cc -c main.c

stack.o : stack.c
	cc -c stack.c

newick.o : newick.c
	cc -c newick.c

write_species.o : write_species.c
	cc -c write_species.c

correct_tip_number.o : correct_tip_number.c
	cc -c correct_tip_number.c

parse_fasta2bin.o : parse_fasta2bin.c
	cc -c parse_fasta2bin.c

opra_matrix.o : opra_matrix.c
	cc -c opra_matrix.c

dist_dna.o : dist_dna.c
	cc -c dist_dna.c

statistic.o : statistic.c
	cc -c statistic.c

write_nwk.o : write_nwk.c
	cc -c write_nwk.c

root.o : root.c
	cc -c root.c

install:
	cp scpc ../bin

clean:
	rm main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o write_species.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o scpc 
