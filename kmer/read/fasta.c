#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void ReadFasta(char *path, char **ref, char **seq)
{
    FILE *fasta = NULL;

    char c;
    int i = 0;
    int j = 0;

    *seq = malloc(2 * sizeof(**seq));
    *ref = malloc(2 * sizeof(**ref));
    if (seq == NULL)
    {
        perror("Error");
        exit(1);
    }

    if (ref == NULL)
    {
        perror("Error");
        exit(1);
    }

    fasta = fopen(path, "r");
    if (fasta == NULL)
    {
        perror("Error");
        exit(1);
    }

    while ((c = fgetc(fasta)) != EOF)
    {
        if (c == '>')
        {
            while (c != '\n' && c != '\r')
            {

                char *tmp = realloc(*ref, (j + 2) * sizeof(*ref));
                if (tmp)
                {
                    *ref = tmp;
                }
                else
                {
                    printf("Reallocation failed!");
                    exit(1);
                }

                (*ref)[j] = c;
                c = fgetc(fasta);
                j++;
            }
            (*ref)[j] = '\0';
        }


        if (c != '\n' && c != '\r')
        {

            char *tmp = realloc(*seq, (i + 1) * sizeof(*seq));
            if (tmp)
            {
                *seq = tmp;
            }
            else
            {
                printf("Reallocation failed!");
                exit(1);
            }

            (*seq)[i] = c;
            i++;
        }
    }
    (*seq)[i] = '\0';
}


static PyObject *read_fasta(PyObject *self, PyObject *args)
{

    char *path, *ref, *seq = NULL;

    if (!PyArg_ParseTuple(args, "s", &path))
    {
        return NULL;
    }

    ReadFasta(path, &ref, &seq);
    return Py_BuildValue("(ss)", ref, seq);
}

static PyMethodDef FastaMethods[] = {
    {"ReadFasta", read_fasta, METH_VARARGS, "Reads FASTA files"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fastamodule = {
    PyModuleDef_HEAD_INIT,
    "fasta",
    NULL,
    -1,
    FastaMethods
};

PyMODINIT_FUNC PyInit_fasta(void)
{
    return PyModule_Create(&fastamodule);
}