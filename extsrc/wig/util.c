#include "common.h"

#ifdef KR_headers
integer pow_ii(ap, bp) integer *ap, *bp;
#else
integer pow_ii(integer *ap, integer *bp)
#endif
{
        integer pow, x, n;

        x = *ap;
        n = *bp;

        if (n <= 0) {
                if (n == 0 || x == 1)
                        return 1;
                if (x != -1)
                        return x == 0 ? 1/x : 0;
                n = -n;
                }
        for(pow = 1; ; )
                {
                if(n & 01)
                        pow *= x;
                if(n >>= 1)
                        x *= x;
                else
                        break;
                }
        return(pow);
        }

