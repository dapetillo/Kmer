//#define PY_SSIZE_T_CLEAN
//#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash.h"

void kmers(char *seq, int k, ffp *hash_table[])
{
    char *subword;
    subword = malloc(k + 1);
    int seq_len = strlen(seq);
    
    ffp *ffp_array = malloc((seq_len - k + 1) * sizeof(ffp));
    for (int i = 0; i < seq_len - k + 1; i++)
    {
        int j = 0;
        while (j < k)
        {
            subword[j] = seq[i + j];
            j++;
        }
        strncpy(ffp_array[i].name, subword, k + 1);
        ffp_array[i].name[k] = '\0';
        ffp_array[i].count = 1;
        hash_table_insert_increment(&ffp_array[i], hash_table);


    }
    
    free(subword);
    free(ffp_array);
}


int main()
{
    ffp *hash_table[TABLE_SIZE];
    init_hash_table(hash_table);
    char seq[] = "ATTTCGATGATTT";
    kmers(seq, 4, hash_table);
    printf("AFTER FUNCTION\n");
    print_table(hash_table);
    printf("END\n");

    return 0;
}