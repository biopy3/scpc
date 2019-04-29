scpc : main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o fixing_tree.o
	cc -o $@ $^

main.o : main.c
	cc -c $^

stack.o : stack.c
	cc -c $^

newick.o : newick.c
	cc -c $^

fixing_tree.o : fixing_tree.c
	cc -c $^

correct_tip_number.o : correct_tip_number.c
	cc -c $^

parse_fasta2bin.o : parse_fasta2bin.c
	cc -c $^

opra_matrix.o : opra_matrix.c
	cc -c $^

dist_dna.o : dist_dna.c
	cc -c $^

statistic.o : statistic.c
	cc -c $^

write_nwk.o : write_nwk.c
	cc -c $^

root.o : root.c
	cc -c $^

install:
	cp scpc ../bin

clean:
	rm main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o fixing_tree.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o scpc 