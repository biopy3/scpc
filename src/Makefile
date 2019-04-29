scpc : main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o fixing_tree.o
	gcc -o $@ $^

main.o : main.c
	gcc -c $^

stack.o : stack.c
	gcc -c $^

newick.o : newick.c
	gcc -c $^

fixing_tree.o : fixing_tree.c
	gcc -c $^

correct_tip_number.o : correct_tip_number.c
	gcc -c $^

parse_fasta2bin.o : parse_fasta2bin.c
	gcc -c $^

opra_matrix.o : opra_matrix.c
	gcc -c $^

dist_dna.o : dist_dna.c
	gcc -c $^

statistic.o : statistic.c
	gcc -c $^

write_nwk.o : write_nwk.c
	gcc -c $^

root.o : root.c
	gcc -c $^

install:
	cp scpc ../bin

clean:
	rm main.o stack.o newick.o parse_fasta2bin.o opra_matrix.o fixing_tree.o correct_tip_number.o dist_dna.o statistic.o write_nwk.o root.o scpc 