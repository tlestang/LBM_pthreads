#ifdef FLIPPED
#define IDX(i,j,k) ((j) + Dy*( k + 9*(i) ) )
#else
#define IDX(i,j,k) ((k) + 9*( (j)+ Dy*(i) ) )
#endif

#define idx(i,j) (j+Dy*i)

#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })

extern const double w[9];
extern const int c[9][2];

