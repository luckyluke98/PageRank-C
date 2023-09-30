#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/sysinfo.h>
#include <string.h>

struct dictionary_entry {
    int key;
    int * value;
    int dim;
    int size;
};

struct dictionary {
    struct dictionary_entry * dictionary;
    int dim;
    int * occ;
};

struct dictionary new_dict(int dim){
    struct dictionary dict;

    dict.dictionary = (struct dictionary_entry *) calloc(dim, sizeof(struct dictionary_entry));

    if (dict.dictionary == NULL) {
        printf("Error malloc\n");
        exit(1);
    }

    dict.dim = dim;

    for (int i = 0; i < dim; i++) {
        dict.dictionary[i].value = (int *) calloc(1, sizeof(int));

        if (dict.dictionary[i].value == NULL) {
            printf("Error malloc\n");
            exit(1);
        }

        dict.dictionary[i].size = 0;
        dict.dictionary[i].dim = 1;
    }

    return dict;
}

void expand_dict(struct dictionary * dict) {
    dict->dictionary = (struct dictionary_entry *) realloc(dict->dictionary, sizeof(struct dictionary_entry) * dict->dim * 2);

    if (dict->dictionary == NULL) {
        printf("Error Realloc\n");
        exit(1);
    }

    for (int i = dict->dim; i < dict->dim * 2; i++) {

        dict->dictionary[i].value = (int *) calloc(1, sizeof(int));

        if (dict->dictionary[i].value == NULL) {
            printf("Error malloc\n");
            exit(1);
        }

        dict->dictionary[i].size = 0;
        dict->dictionary[i].dim = 1;
    }
    dict->dim = dict->dim * 2;
}

void add_to_dict(struct dictionary * dict, int key, int value) {

    if (key >= dict->dim)
        expand_dict(dict);

    if (dict->dictionary[key].dim == dict->dictionary[key].size) {
        dict->dictionary[key].value = (int *) realloc(dict->dictionary[key].value, sizeof(int) * dict->dictionary[key].dim * 2);
        
        if (dict->dictionary[key].value == NULL) {
            printf("Error Realloc\n");
            exit(1);
        }

        dict->dictionary[key].dim = dict->dictionary[key].dim * 2;
    }

    dict->dictionary[key].key = key;
    dict->dictionary[key].value[dict->dictionary[key].size] = value;
    dict->dictionary[key].size = dict->dictionary[key].size + 1;
}

FILE * skip_comment(char * filename, int * nodes, int * edges) {
    FILE * fp;
    char c, str[100];

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Cannot open the file");
        exit(1);
    }

    c = getc(fp);

    while (c == '#') {
        fgets(str, 100 - 1, fp);
        sscanf(str, "%*s %d %*s %d", nodes, edges);
        c = getc(fp);
    }
    ungetc(c, fp);

    return fp;
}

int check(char * filename) {
    FILE * fp;
    int nodes, edges;

    fp = skip_comment(filename, &nodes, &edges);

    int fromnode, tonode;
    int * occ = (int *) calloc(nodes, sizeof(int));

    while (!feof(fp)) {
        fscanf(fp, "%d%d", &fromnode, &tonode);
        occ[fromnode] = 1;
        occ[tonode] = 1;
    }

    if (occ[0] == 0)
        return 1;

    return 0;

}

void order_dataset(char * filename, char * out_filename, int * nodes, int * edges) {
    FILE * fp, * fp_out;
    
    int fromnode = 0, tonode = 0;
    int old_fromnode = -1, old_tonode = -1;

    fp = skip_comment(filename, nodes, edges);

    // Like "Single pass in memory"
    struct dictionary dict = new_dict(*nodes);

    fscanf(fp, "%d %d", &fromnode, &tonode);

    while (!feof(fp)) {
        add_to_dict(&dict, tonode, fromnode);
        fscanf(fp, "%d %d", &fromnode, &tonode);
    }

    if ((fp_out = fopen(out_filename, "w")) == NULL) {
        printf("Unable to create file.\n");
        exit(1);
    }
    
    int start_with_one = check(filename);

    for (int i = 0; i < dict.dim; i++) {
        for (int j = 0; j < dict.dictionary[i].size; j++) {
            if (start_with_one)
                fprintf(fp_out, "%d\t%d\n", dict.dictionary[i].key - 1, dict.dictionary[i].value[j] - 1);
            else
                fprintf(fp_out, "%d\t%d\n", dict.dictionary[i].key, dict.dictionary[i].value[j]);
        }
    }
    fclose(fp_out);

    for (int i = 0; i < dict.dim; i++) {
        free(dict.dictionary[i].value);
    }
    free(dict.dictionary);
}

void page_rank(char * filename, int nodes, int edges) {
    clock_t start, end;
    double total_time;
    start = clock();

    FILE *fp;
    int fromnode, tonode;

    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Cannot open the file");
        exit(1);
    }

    // CSR
    float *val = (float *) calloc(edges, sizeof(float));
    int *col_ind = (int *) calloc(edges, sizeof(int));
    int *row_ptr = (int *) calloc(nodes + 1, sizeof(int));
    int * out_link = (int *) calloc(nodes, sizeof(int));

    if (val == NULL || col_ind == NULL || row_ptr == NULL || out_link == NULL) {
        printf("Error malloc\n");
        exit(1);
    }

    int j, i = 0;
    int elrow = 0, curel = 0, row = 0;

    row_ptr[0] = 0;

    while (!feof(fp)) {
        fscanf(fp, "%d%d", &fromnode, &tonode);

        if (fromnode > row) {
            curel = curel + elrow;
            for (int k = row + 1; k <= fromnode; k++) {
                row_ptr[k] = curel;
            }
            elrow = 0;
            row = fromnode;
        }
        val[i] = 1.0;
        col_ind[i] = tonode;
        elrow++;
        i++;

    }
    row_ptr[row + 1] = curel + elrow - 1;
    row++;

    while (row < nodes) {
        row_ptr[row + 1] = row_ptr[row];
        row++;
    }

    for (int i = 0; i < edges; i++) {
        out_link[col_ind[i]] = out_link[col_ind[i]] + 1;
    }

    for (int i = 0; i < edges; i++) {
        val[i] = val[i] / out_link[col_ind[i]];
    }

    float p[nodes];
    for (i = 0; i < nodes; i++)
        p[i] = 1.0 / nodes;

    int looping = 1;
    int k = 0;
    float p_new[nodes];
    float d = 0.85; 
    int num_elem_row = 0;
    float dandling_nodes = 0;

    while (looping) {
        
        for (i = 0; i < nodes; i++)
            p_new[i] = 0.0;

        num_elem_row = 0;
        dandling_nodes = 0;
        int pt = 0;

        for (i = 0; i < nodes; i++) {
            num_elem_row = row_ptr[i + 1] - row_ptr[i];
            for (j = 0; j < num_elem_row; j++) {
                p_new[i] = p_new[i] + (val[pt] * p[col_ind[pt]]);
                pt++;
            }
        }

        for (int j = 0; j < nodes; j++)
            if (out_link[j] == 0)
                dandling_nodes = dandling_nodes + p[j] / nodes;

        for (int i = 0; i < nodes; i++)
            p_new[i] = p_new[i] + dandling_nodes;
        
        for (i = 0; i < nodes; i++)
            p_new[i] = d * p_new[i] + (1.0 - d) / nodes;
        
        float epsilon = 0.0;
        for (i = 0; i < nodes; i++)
            epsilon = epsilon + fabs(p_new[i] - p[i]);
       
        if (epsilon < 0.000001)
            looping = 0;

        for (i = 0; i < nodes; i++)
            p[i] = p_new[i];

        k = k + 1;
    }

    end = clock();
    total_time = (double)(end - start) / CLOCKS_PER_SEC;

    printf("\nNumber of iteration: %d \n\n", k);

    float sum = 0;
    for (i = 0; i < nodes; i++) {
        sum = sum + p[i];
    }

    printf("Sum: %f \nTime spent: %f seconds.\n", sum, total_time);
    
}

int main() {

    
    char filename[] = "web-Stanford.txt";
    char filename_out[] = "web-BerkStan-out.txt";
    char * ordered_dataset;
    int nodes, edges;

    order_dataset(filename, filename_out, &nodes, &edges);

    page_rank(filename_out, nodes, edges);

    printf("%d %d\n", nodes, edges);
    
    return 1;

}
