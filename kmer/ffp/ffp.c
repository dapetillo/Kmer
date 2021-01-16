#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash.h"

ffp *kmers(char *seq, int k, ffp *hash_table[], int table_size)
{
    char *subword;
    subword = malloc(k + 1);
    int seq_len = strlen(seq);
    
    ffp *ffp_array = malloc((seq_len - k + 1) * sizeof(*ffp_array));
    //*ffp_array = malloc((seq_len - k + 1) * sizeof(**ffp_array));
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
        hash_table_insert_increment(&ffp_array[i], hash_table, table_size);
    }
    free(subword);

    return ffp_array;
}


static PyObject *freq_feat_profile(PyObject *self, PyObject *args)
{
    char *seq = NULL;
    int k = 1;

    if (!PyArg_ParseTuple(args, "si", &seq, &k))
    {
        return NULL;
    }

    int table_size = 1;
    for (int i = 0; i < k; i++)
    {
        table_size *= 4;
    }

    ffp *hash_table[table_size];
    ffp *ffp_array;
    init_hash_table(hash_table, table_size);
    ffp_array = kmers(seq, k, hash_table, table_size);

    PyObject *dict = PyDict_New();
    for (int i = 0; i < table_size; i++)
    {
        if (hash_table[i] == NULL)
        {
            continue;
        }
        else
        {
            PyDict_SetItem(dict, PyUnicode_FromString(hash_table[i]->name), PyLong_FromLong(hash_table[i]->count));
        }
    }

    return dict;
}

static PyMethodDef FreqFeatMethods[] = {
    {"ffp", freq_feat_profile, METH_VARARGS, "Compute frequency feature profile"},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef ffpmodule = {
    PyModuleDef_HEAD_INIT,
    "ffp",
    NULL,
    -1,
    FreqFeatMethods
};

PyMODINIT_FUNC PyInit_ffp(void)
{
    return PyModule_Create(&ffpmodule);
}