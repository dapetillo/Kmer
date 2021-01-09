//#define PY_SSIZE_T_CLEAN
//#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define KORDER 10
#define TABLE_SIZE 256 //1048576 //4**KORDER

typedef struct
{
    char name[KORDER];
    int count;
} ffp;

//ffp *hash_table[TABLE_SIZE];


int get_size(char *name)
{
    int length = strlen(name);

    int table_size = 1;
    for (int i = 0; i < length; i++)
    {
        table_size *= length;
    }
    
    return table_size;
}

/* it is a perfect hash function. Provided that
 * TABLE_SIZE is large enough (at least 4**KORDER),
 * collisions will not happen.
 */
unsigned int hash(char *name)
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
    hash_value %= TABLE_SIZE;
    return hash_value;
}

void init_hash_table(ffp *hash_table[])
{
    for (int i = 0; i < TABLE_SIZE; i++)
    {
        hash_table[i] = NULL;
    }
}

void print_table(ffp *hash_table[])
{
    for (int i = 0; i < TABLE_SIZE; i++)
    {
        if (hash_table[i] == NULL)
        {
            printf("\t%i\t---\n", i);
        } else {
            printf("\t%i\t%s\t%d\n", i, hash_table[i]->name, hash_table[i]->count);
        }
    }
}

bool hash_table_insert_increment(ffp *p, ffp *hash_table[])
{
    if (p == NULL)
    {
        return false;
    }
    int index = hash(p->name);
    if (hash_table[index] != NULL)
    {
        printf("Increment %s\n", p->name);
        hash_table[index]->count += 1;
        return true;
    } else {
        printf("Add %s with index %d\n", p->name, index);
        hash_table[index] = p;
        //hash_table[index]->count = 1;
        //printf("From hash table is %d\n", hash_table[index]->count);
    }
    
    return true;
}


/*
static PyObject *InitHashTable(PyObject *self, PyObject *args)
{
    init_hash_table();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *InsertIncrement(PyObject *self, PyObject *args)
{
    char *word = NULL;

    if (!PyArg_ParseTuple(args, "s", &word))
    {
        return NULL;
    }

    ffp kmer = {.name = word, .count = 1};
    hash_table_insert_increment(&kmer);

    return Py_BuildValue("{s:i}", hash_table);


}
*/

/*
int main()
{
    init_hash_table();

    ffp boh = {.name = "AT", .count = 1};
    ffp ete = {.name = "TA", .count = 1};
    ffp tte = {.name = "GC", .count = 1};
    ffp aga = {.name = "TA", .count = 1};
    ffp taa = {.name = "TT", .count = 1};

    hash_table_insert_increment(&boh);
    hash_table_insert_increment(&ete);
    hash_table_insert_increment(&tte);
    hash_table_insert_increment(&aga);
    hash_table_insert_increment(&tte);
    hash_table_insert_increment(&taa);
    print_table();
    
    printf("ATTC => %u\n", hash("ATTC"));
    printf("AGTG => %u\n", hash("AGTG"));
    printf("ACTT => %u\n", hash("ACTT"));
    printf("GTAG => %u\n", hash("GTAG"));
    printf("AAAA => %u\n", hash("AAAA"));
    printf("GGTA => %u\n", hash("GGTA"));
    printf("CCTA => %u\n", hash("CCTA"));
    printf("TCTC => %u\n", hash("TCTC"));

    return 0;
}
*/