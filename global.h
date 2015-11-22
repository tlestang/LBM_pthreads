#define IDX(i,j,l) ( 9*( (j) +Dy*(i) ) + (l) )
#define idx(i,j) (j+Dy*i)

#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })

extern const double w[9];
extern const int c[9][2];

