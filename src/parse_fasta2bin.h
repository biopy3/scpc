int chec_file(char *path,int num_seqs, int len_seq[], int len_name[], int *know_bases);
unsigned char *init_seq(int num_seqs, int len_seq, int len_seq_name[]);
void init_seq_name(char *seq_name[], int num_seqs, int len_seq_name[]);
void free_sfas(char *Seqs);
void free_seq_name(char *seq_name[], int num_seqs);
int raw2DNAbin(unsigned char * Seqs,char *seq_name[], char * input_path, int len_seq);