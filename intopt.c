#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPSILON 1e-6
#define ROW_SIZE (5)
#define COL_SIZE (2)

int glob = 0;

struct simplex_t {
int m ; /* Constraints. */
int n ; /* Decision variables. */
int *var; /* 0..n  1 are nonbasic. */
double **a; /* A. */
double *b; /* b. */
double *x; /* x. */
double *c; /* c. */
double y ; /* y. */
};

typedef struct simplex_t simplex_t;

int init (simplex_t *s, int m,  int n, double **a, double *b, double *c, double *x, double y, int *var)
{
    int i,k;
    //*s = (simplex_t) {s->m=m,s->n=n,s->var=var,s->a=a,s->b=b,s->x=x,s->c=c,s->y=y}; // assign each attribute
    *s=(simplex_t){m,n,var,a,b,x,c,y};
    if (s->var == NULL){
        s->var =  calloc(m+n+1, sizeof(int));
        for (i = 0; i < m+n; i = i + 1){
            s->var[i] = i;
        }
    }
    for (k = 0, i = 1; i < m; i = i + 1){
        if (s->b[i] < s->b[k]){
        k=i;
        }
    }
    
    return k;
}

int select_nonbasic (simplex_t *s){
    int i;
    for (i = 0; i < s->n; i = i + 1){
        if (s->c[i] > EPSILON) {
            return i;
        }
    }
    return -1;
}

void pivot (simplex_t *s,int row,int col) //changed to pointer so it affects the actual object
{
    glob += 1;
    double **a = s->a;
    double *b = s->b;
    double *c = s->c;
    int m = s->m;
    int n = s->n;
    int i,j,t;
    t = s->var[col];
    s->var[col] = s->var[n+row];
    s->var[n+row] = t;
    s->y = s->y + c[col] * b[row] / a[row][col];
    for (i = 0; i < n; i = i + 1){
        if (i != col){
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }
    c[col] = - c[col] / a[row][col];
    for (i = 0; i < m; i = i + 1){
        if (i != row){
            b[i] = b[i] - a[i][col] * b[row] / a[row][col];
        }
    }
    for (i = 0; i < m; i = i + 1){
        if (i != row){
            for (j = 0; j < n; j = j + 1){
                if (j != col){
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }
    for (i = 0; i < m; i = i + 1){
        if (i != row){
            a[i][col] = -a[i][col] / a[row][col];
        }
    }
    for (i = 0; i < n; i = i + 1){
        if (i != col){
            a[row][i] = a[row][i] / a[row][col];
        }
    }
    b[row] = b[row] / a[row][col];
    a[row][col] = 1 / a[row][col];
}

// void prepare (simplex_t *s, int k){
//     printf("prepare \n");
//     int m = s->m;
//     int n = s->n;
//     int i;
//     // make room for xm+n at s.var[n] by moving s.var[n..n+m-1] one
//     // step to the right.
//     for (i = m + n; i > n; i = i - 1){
//         s->var[i] = s->var[i-1];
//     }
//     s->var[n] = m + n;
//     // add xm+n to each constraint
//     n = n + 1;
//     for (i = 0; i < m; i = i + 1){
//         s->a[i][n-1]= -1;
//     }
//     s->x = calloc(m+n,sizeof(double));
//     s->c = calloc(n,sizeof(double));
//     s->c[n-1] = -1;
//     s->n = n;
//     pivot(s, k, n-1); //sent in the address of it so it modifies the actual s
// }

int initial(simplex_t *s, int m,  int n, double **a, double *b, double *c, double *x, double y, int *var){
    int i,j,k;
    double w;
    k = init(s, m, n, a, b, c, x, y, var);
    if (b[k] >= 0) return 1; // feasible
}

double xsimplex (int m,  int n, double **a, double *b, double *c, double *x, double y, int *var, int h){
    simplex_t s;
    int i,row,col;
    if (!initial(&s, m, n, a, b,c, x, y, var)){ //&s
        free(s.var);
        return NAN; // not a number
    }
    while ((col = select_nonbasic(&s)) >= 0) {
        row = -1;
        for (i = 0; i < m; i = i + 1){
            if (a[i][col] > EPSILON && (row < 0 || (b[i] / a[i][col] < b[row] / a[row][col]))){
            row = i;
            }
        }
        if (row < 0){
            free(s.var); //maybe ->?
            return INFINITY; // unbounded
        }
        pivot (&s,row, col);
    }
    if (h==0){
        for (i = 0; i < n; i = i + 1){
            if (s.var[i] < n){
                x[s.var[i]] = 0;
            }
        }
        for (i = 0; i < m; i = i + 1){
            if (s.var[n+i] < n){
                x[s.var[n+i]] = s.b[i];
            }
        }
        free(s.var);
    }
    else{
        for (i = 0; i < n; i = i + 1){
            x[i] = 0;
        }
        for (i = n; i < n+m; i = i + 1){
            x[i] = s.b[i-n];
        }
    }
    return s.y;
    }

double simplex (int m,  int n, double **a, double *b, double *c, double *x, double y){
    return xsimplex (m,n,a,b,c,x,y,NULL,0);
}

double** make_matrix(int m, int n)
{
    double**    a;
    int         i;
    a = calloc(m, sizeof(double*));
    for (i = 0; i < m; i += 1)
    a[i] = calloc(n, sizeof(double));
    return a;
}

int main(int agrc, char** argv)
{
    int m;
    int n;

    scanf("%d %d", &m, &n);

    double*     c = calloc(n, sizeof(double));
    double*     b = calloc(m, sizeof(double));
    double**    a = make_matrix(m,n);

    //scan in all the c coefficients for the max/min solution
    for (int i=0 ; i<n ; i+=1){
        scanf("%lf", &c[i]);
    }
    //scan the A matrix coefficients for the system
    for (int i=0 ; i<m; i+=1){
        for (int j=0 ; j<n; j+=1){
            scanf("%lf", &a[i][j]);
        }
    }
    //scan in the right hand side of the system (b)
    for (int i=0 ; i<n; i+=1){
        scanf("%lf", &b[i]);
    }
    printf("max z = ");
    for (int i=0 ; i<m ; i+=1){
        printf("%10.3lfx%d", c[i],i);
        if (i<n-1) printf(" + ");
    }
    printf("\n");

    for (int i=0 ; i<m ; i+=1){
        for (int j=0 ; j<n ; j+=1){
            printf("%10.3lfx%d", a[i][j], j);
            if (j<n-1) printf(" + ");
        }
        printf("<=");
        printf("%10.3lf\n", b[i]);
    }

    printf("hejsan \n");
    double *x = calloc(n+1,sizeof(double));

    double y = simplex (m,  n, a, b, c, x, 0);//work

    printf("%lf \n", y);

    for(int i=0; i<m; i+=1)
    {
        free(a[i]);
    }
    free(a);
    free(b);
    free(c);
    return 0;
}