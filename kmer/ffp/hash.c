#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define KORDER 20

typedef struct
{
    char name[KORDER];
    int count;
} ffp;


/* it is a perfect hash function. Provided that
 * table_size is large enough (at least 4**length),
 * collisions will not happen.
 */
unsigned int hash(char *name, int table_size)
{
    int length = strlen(name);

    unsigned int hash_value = 0;
    for (int i = 0; i < length; i++)
    {
        switch (name[i])
        {
            case 'A':
                hash_value += 0 * pow(4, i);
                break;
            case 'T':
                hash_value += 1 * pow(4, i);
                break;
            case 'C':
                hash_value += 2 * pow(4, i);
                break;
            case 'G':
                hash_value += 3 * pow(4, i);
                break;
        }
    }
    hash_value %= table_size;
    return hash_value;
}

void init_hash_table(ffp *hash_table[], int table_size)
{
    for (int i = 0; i < table_size; i++)
    {
        hash_table[i] = NULL;
    }
}


void print_table(ffp *hash_table[], int table_size)
{
    for (int i = 0; i < table_size; i++)
    {
        if (hash_table[i] == NULL)
        {
            printf("\t%i\t---\n", i);
        } else {
            printf("\t%i\t%s\t%d\n", i, hash_table[i]->name, hash_table[i]->count);
        }
    }
}


bool hash_table_insert_increment(ffp *p, ffp *hash_table[], int table_size)
{
    if (p == NULL)
    {
        return false;
    }
    int index = hash(p->name, table_size);
    if (hash_table[index] != NULL)
    {
        hash_table[index]->count += 1;
        return true;
    } else {
        hash_table[index] = p;
    }
    
    return true;
}