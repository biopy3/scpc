#include "scpc.h"
#include "parse_fasta2bin.h"

extern int cterminals;

#define SSS(lable) \
    tempc = fgetc(f); \
    lable:if(tempc == EOF) \
    { \
        if(feof(f)) \
        { \
            ; \
        } \
        else { \
            printf("\n Something went wrong, please try again.\n"); \
            return -1; \
        } \
    }else if ( tempc == '\r' || tempc == '\n' || tempc == '\t' || tempc == ' ' ){ \
        tempc = fgetc(f); \
        goto lable; \
    }else if(tempc == '>') \
    { \
        count_seqs_num++; \
        goto read_and_encode_seq; \
    }else{ \
        Seqs[count_seqs_num*len_seq + count_seq_len++] = tab_trans[tempc]; \
        tempc = fgetc(f); \
        goto lable; \
    }

#define AAA(lable) \
    tempc = fgetc(f); \
    if( cterminals < num_seqs)\
    { printf("The number of sequences great than the number of terminals.\n"); \
     return EXIT_FAILURE;}\
    lable:if(tempc == EOF) \
    { \
        if(feof(f)) \
        { \
            if(len_seq[num_seqs-1] == len_seq[0]) \
            { \
            }else{ \
                printf("The sequence is not same length, please check!\n"); \
                return -1; \
            } \
        } \
        else { \
            printf("\n Something went wrong, please try again.\n"); \
            return -1; \
        } \
    }else if ( tempc == '\r' || tempc == '\n' || tempc == '\t' || tempc == ' ' ){ \
        tempc = fgetc(f); \
        goto lable; \
    }else if(tempc == '>') \
    { \
        if(len_seq[num_seqs-1] == len_seq[0]) \
        { \
        }else{ \
            printf("The sequence is not same length, please check!\n"); \
            return -1; \
        } \
        goto chec_seqs; \
    }else{ \
        if( (tab_trans[tempc]) != '\0') \
        { \
			int s = KnownBase(tab_trans[tempc]); \
			if(  KnownBase(tab_trans[tempc]) ) \
			{ \
				(*know_bases)++; \
			} \
            len_seq[num_seqs-1]++; \
            tempc = fgetc(f); \
            goto lable; \
        } \
        else{ \
            printf("The FASTA file is't correct.Only the IUAPC code (plus '-' and '?') is used, please check!\n"); \
			return -1; \
        } \
    }

// The initial code defining and initialising the translation table:
//
//	for (i = 0; i < 122; i++) tab_trans[i] = 0x00;
//
//	tab_trans[65] = 0x88; /* A */
//	tab_trans[71] = 0x48; /* G */
//	tab_trans[67] = 0x28; /* C */
// 	tab_trans[84] = 0x18; /* T */
// 	tab_trans[82] = 0xc0; /* R */
// 	tab_trans[77] = 0xa0; /* M */
// 	tab_trans[87] = 0x90; /* W */
// 	tab_trans[83] = 0x60; /* S */
// 	tab_trans[75] = 0x50; /* K */
// 	tab_trans[89] = 0x30; /* Y */
// 	tab_trans[86] = 0xe0; /* V */
// 	tab_trans[72] = 0xb0; /* H */
// 	tab_trans[68] = 0xd0; /* D */
//  tab_trans[66] = 0x70; /* B */
// 	tab_trans[78] = 0xf0; /* N */
//
//	tab_trans[97] = 0x88; /* a */
//	tab_trans[103] = 0x48; /* g */
//	tab_trans[99] = 0x28; /* c */
// 	tab_trans[116] = 0x18; /* t */
// 	tab_trans[114] = 0xc0; /* r */
// 	tab_trans[109] = 0xa0; /* m */
// 	tab_trans[119] = 0x90; /* w */
// 	tab_trans[115] = 0x60; /* s */
// 	tab_trans[107] = 0x50; /* k */
// 	tab_trans[121] = 0x30; /* y */
// 	tab_trans[118] = 0xe0; /* v */
// 	tab_trans[104] = 0xb0; /* h */
// 	tab_trans[100] = 0xd0; /* d */
//  	tab_trans[98] = 0x70; /* b */
// 	tab_trans[110] = 0xf0; /* n */
//
//  	tab_trans[45] = 0x04; /* - */
//  	tab_trans[63] = 0x02; /* ? */

static const unsigned char tab_trans[] = {
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 0-9 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 10-19 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 20-29 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 30-39 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, /* 40-49 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 50-59 */
	0x00, 0x00, 0x00, 0x02, 0x00, 0x88, 0x70, 0x28, 0xd0, 0x00, /* 60-69 */
	0x00, 0x48, 0xb0, 0x00, 0x00, 0x50, 0x00, 0xa0, 0xf0, 0x00, /* 70-79 */
	0x00, 0x00, 0xc0, 0x60, 0x18, 0x00, 0xe0, 0x90, 0x00, 0x30, /* 80-89 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x88, 0x70, 0x28, /* 90-99 */
	0xd0, 0x00, 0x00, 0x48, 0xb0, 0x00, 0x00, 0x50, 0x00, 0xa0, /* 100-109 */
	0xf0, 0x00, 0x00, 0x00, 0xc0, 0x60, 0x18, 0x00, 0xe0, 0x90, /* 110-119 */
	0x00, 0x30, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 120-129 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 130-139 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 140-149 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 150-159 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 160-169 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 170-179 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 180-189 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 190-199 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 200-209 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 210-219 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 220-229 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 230-239 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 240-249 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }; /* 250-255 */

static const unsigned char hook = 0x3e;
static const unsigned char lineFeed = 0x0a;
/* static const unsigned char space = 0x20; */

/* translation table DNAbin -> CHAR */
static const unsigned char tab_trans_rev[] = {
	0x00, 0x00, 0x3f, 0x00, 0x2d, 0x00, 0x00, 0x00, 0x00, 0x00, /* 0-9 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 10-19 */
	0x00, 0x00, 0x00, 0x00, 0x54, 0x00, 0x00, 0x00, 0x00, 0x00, /* 20-29 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 30-39 */
	0x43, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x59, 0x00, /* 40-49 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 50-59 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 60-69 */
	0x00, 0x00, 0x47, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 70-79 */
	0x4b, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 80-89 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x53, 0x00, 0x00, 0x00, /* 90-99 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 100-109 */
	0x00, 0x00, 0x42, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 110-119 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 120-129 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x41, 0x00, 0x00, 0x00, /* 130-139 */
	0x00, 0x00, 0x00, 0x00, 0x57, 0x00, 0x00, 0x00, 0x00, 0x00, /* 140-149 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 150-159 */
	0x4d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 160-169 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x48, 0x00, 0x00, 0x00, /* 170-179 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 180-189 */
	0x00, 0x00, 0x52, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 190-199 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x44, 0x00, /* 200-209 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 210-219 */
	0x00, 0x00, 0x00, 0x00, 0x56, 0x00, 0x00, 0x00, 0x00, 0x00, /* 220-229 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 230-239 */
	0x4e, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 240-249 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }; /* 250-255 */

int chec_file(char *path,int cterminals, int len_seq[], int len_name[], int *know_bases)
{
    int i, num_seqs = 0;
    *know_bases = 0;
    for( i=0; i<cterminals; i++)
    {
        len_name[i] = 0;
        len_seq[i] = 0;
    }
  
    FILE * f = fopen(path, "rt");
	if ( f == NULL )
	{
		printf("Can't open the input fasta file, may be no permission or something, please check then try again!\n");
		exit(1);
	}
    
    int tempc = fgetc(f);
	chec_seqs:if(tempc == EOF)
	{
		if(feof(f))
		{
			;
		}else{
			printf("The file handle is break or something, please check then try again!\n");
			return -1;
		}
	}else if(tempc == '>'){
		tempc = fgetc(f);
		(num_seqs)++;
//
		reac_seq_name_end:if(tempc == EOF)
		{
			if(feof(f))
			{
				return 0;
			}
			else{
				printf("The file handle is break or something, please check then try again!\n");
				return -1;
			}
		}else if(tempc == '\r'){
			if ((tempc = fgetc(f)) == '\n')
			{
				AAA(_01)
			}
			else{
				fseek(f, -1L, SEEK_CUR);
				AAA(_02)
			}
		}else if(tempc == '\n')
		{
			AAA(_03)
		}else{
			len_name[num_seqs-1]++;
			tempc = fgetc(f);
			goto reac_seq_name_end;
		}
	}else{
        tempc = fgetc(f);
        goto chec_seqs;
    }

	fclose(f);
  if( cterminals != num_seqs)\
  { 
    printf("The number of sequences isn't equal the number of terminals.\n");
    return EXIT_FAILURE;
  }
  return 0;
}

unsigned char *init_seq(int num_seqs, int len_seq, int len_seq_name[])
{
	unsigned char * Seqs = (unsigned char *)malloc(num_seqs*len_seq*sizeof(unsigned char));
	if(!Seqs)
	{
		printf("Can't malloc the memory for Seqs!\n");
		exit(1); 
	}
	printf("Init!\n");
	return Seqs;
}

void init_seq_name(char *seq_name[], int num_seqs, int len_seq_name[])
{
	int i;
    for( i=0; i<num_seqs; i++)
    {
        seq_name[i] = (char *)malloc((len_seq_name[i] + 1)*sizeof(char));
        if(!seq_name[i])
        {
            printf("\nCan't malloc the memory for seq_name[%d]\n", i);
            i--;
            for(; 0 <= i; i--)
            {
                free(seq_name[i]);
            }
            exit( 1);
        }else{
            seq_name[i][len_seq_name[i]] = '\0';
        }
    }
}

void free_sfas(char *Seqs)
{
	free(Seqs);
}

void free_seq_name(char *seq_name[], int num_seqs)
{
	int i;
    for( i=0; i<num_seqs; i++)
    {
        free(seq_name[i]);
    }
}
int raw2DNAbin(unsigned char * Seqs,char *seq_name[], char * input_path, int len_seq)
{
	FILE *f = fopen(input_path, "rt");
	int tempc, count_seqs_num = 0, count_seq_len = 0, count_seq_name = 0;
	if (!f)
	{
		printf("The file handle is break or something, please check then try again!\n");
		return -1;
	}

	tempc = fgetc(f);

	read_and_encode_seq:if(tempc == EOF)
	{
		if(feof(f))
		{
			;
		}
		else{
			printf("\n Something went wrong, please try again.\n");
			return -1;
		}
	}else if(tempc == '>')
	{
		count_seq_len = 0;
		count_seq_name = 0;
		tempc = fgetc(f);
		reac_seq_name_end:if(tempc == EOF)
		{
			if(feof(f))
			{
				;
			}
			else{
				printf("\n Something went wrong, please try again.\n");
				return -1;
			}
		}else if(tempc == '\r')
		{
			if ((tempc = fgetc(f)) == '\n')
			{
				SSS(read_seq2end01)
			}else{
				fseek(f, -1L, SEEK_CUR);
				SSS(read_seq2end02)
			}
		}else if(tempc == '\n')
		{
			SSS(read_seq2end03)
		}
		else{
            seq_name[count_seqs_num][count_seq_name++] = tempc;
			tempc = fgetc(f);
			goto reac_seq_name_end;
		}
	}else{
		tempc = fgetc(f);
		goto read_and_encode_seq;
	}
	
	fclose(f);
	return 0;
}