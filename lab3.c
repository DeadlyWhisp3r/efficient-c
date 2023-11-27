#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPSILON 1e-6
#define ROW_SIZE (5)
#define COL_SIZE (2)
int local_array[10];

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
        // s->var =  malloc((m+n+1)* sizeof(int));
        //remove i (= 0) for assignment 6
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
        for (i = 0; i < 11; i += 1) local_array[i] = i;
    for (i = 0; i < s->n; i = i + 1){
        if (s->c[i] > EPSILON) {
            return i;
        }
    }
    return -1;
}

void pivot (simplex_t *s,int row,int col) //changed to pointer so it affects the actual object
{
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

void prepare (simplex_t *s, int k){
    printf("prepare \n");
    int m = s->m;
    int n = s->n;
    int i;
    // make room for xm+n at s.var[n] by moving s.var[n..n+m-1] one
    // step to the right.
    for (i = m + n; i > n; i = i - 1){
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m + n;
    // add xm+n to each constraint
    n = n + 1;
    for (i = 0; i < m; i = i + 1){
        s->a[i][n-1]= -1;
    }
    s->x = calloc(m+n,sizeof(double));
    s->c = calloc(n,sizeof(double));
    s->c[n-1] = -1;
    s->n = n;
    pivot(s, k, n-1); //sent in the address of it so it modifies the actual s
}
double xsimplex (int m,  int n, double **a, double *b, double *c, double *x, double y, int *var, int h); //forward declaration

int initial(simplex_t *s, int m,  int n, double **a, double *b, double *c, double *x, double y, int *var){
    int i,j,k;
    double w;
    k = init(s, m, n, a, b, c, x, y, var);
    if (b[k] >= 0){
        return 1; // feasible
    }
    prepare(s,k);
    n = s->n;
    s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var,1);
    for (i = 0; i < m+n; i = i + 1) {
        if (s->var[i] == m+n-1){
            if (abs(s->x[i]) > EPSILON) {
                free(s->x);
                free(s->c);
                return 0; // infeasible
            }
            else{
                break; // This i will be used on the next page.
            }
        }
    }
    if (i >= n){
        // xn+m is basic. ﬁnd good nonbasic.
        for (j = k = 0; k < n; k = k + 1){
            if (abs(s->a[i-n][k]) > abs(s->a[i-n][j])){
                j=k;
            }
        }
        pivot(s,i-n,j);
        i=j;
    }
    if (i < n-1){
        // xn+m is nonbasic and not last. swap columns i and n-1
        k = s->var[i]; 
        s->var[i] = s->var[n-1];
        s->var[n-1] = k;
        for (k = 0; k < m; k = k + 1){
            w = s->a[k][n-1];
            s->a[k][n-1] = s->a[k][i];
            s->a[k][i] = w;
        }
    }
    else{}
        // xn+m is nonbasic and last. forget it.
    free(s->c);
    s->c = c;
    s->y = y;
    for (k = n-1; k < n+m-1; k = k + 1){
        s->var[k] = s->var[k+1];
    }
    n = s->n = s->n - 1;
    double *t = calloc(n,sizeof(double));
    for (k = 0; k < n; k = k + 1) {
        for (j = 0; j < n; j = j + 1){
            if (k == s->var[j]){
                // xk is nonbasic. add ck
                t[j] = t[j] + s->c[k];
                goto next_k;
            }
        }
        // xk is basic.
        for (j = 0; j < m; j = j + 1){
            if (s->var[n+j] == k){
            // xk is at row j
                break;
            }
        }
        s->y = s->y + s->c[k] * s->b[j];
        for (i = 0; i < n; i = i + 1){
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }
    next_k:;
    }
    for (i = 0; i < n; i = i + 1){
        s->c[i] = t[i];
    }
    free(t);
    free(s->x);
    return 1;
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
    double**    a = make_matrix(m,n+1);

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
    for (int i=0 ; i<m; i+=1){ //changed to m since that is the number of constraints
        scanf("%lf", &b[i]);
    }
    printf("max z = ");
    for (int i=0 ; i<n ; i+=1){
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
        printf("%10.3lf\n", b[i]); //here the m was already chosen
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
    free(x);
    return 0;
}