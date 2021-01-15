#ifndef _HASH_H_
#define _HASH_H_
#include <stdbool.h>

#define KORDER 20

typedef struct
{
    char name[KORDER];
    int count;
} ffp;

//extern ffp *hash_table[TABLE_SIZE];

int get_size(char *name);

unsigned int hash(char *name, int table_size);

void init_hash_table(ffp *hash_table[], int table_size);

void print_table(ffp *hash_table[], int table_size);

bool hash_table_insert_increment(ffp *p, ffp *hash_table[], int table_size);

#endif