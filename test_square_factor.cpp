// Standalone test of the decomposition stage. Compile the production source
// into this TU so no duplicate implementation of the tested code exists.
#define main conway_fast_embedded_main
#include "conway_fast.cpp"
#undef main
#include <cassert>
#include <random>

static std::pair<U64,U64> reference(U64 n){
    U64 rem=n,s=1;
    for(U64 p=2;p<=rem/p;++p){
        unsigned e=0;
        while(rem%p==0){rem/=p;++e;}
        for(unsigned j=0;j<e/2;++j)s*=p;
    }
    return {n/(s*s),s};
}
static bool squarefree(U64 n){
    for(U64 p=2;p*p<=n;++p)if(n%(p*p)==0)return false;
    return true;
}
int main(){
    for(U64 d=2;d<=100000;++d){
        auto actual=discover_square_factor(d);
        auto expected=reference(d);
        if(actual!=expected){
            std::cerr<<"failed exhaustive n="<<d<<"\n";return 1;
        }
    }
    std::mt19937_64 gen(0xC0DE);
    // All below 10^12; exact squarefree-part discovery must be complete.
    for(unsigned i=0;i<5000;++i){
        U64 d=2+gen()%2500;
        if(!squarefree(d))continue;
        U64 s=1+gen()%10000;
        U64 n=d*s*s;
        if(n>1000000000000ULL)continue;
        auto got=discover_square_factor(n);
        if(got.first!=d || got.second!=s){
            std::cerr<<"failed structured n="<<n<<" d="<<d<<" s="<<s<<" got="<<got.first<<","<<got.second<<"\n";
            return 1;
        }
    }
    // Above 10^12 the bounded method may return only a partial factor:
    // the decomposition MUST still be exact, and must not mislabel the result.
    for(unsigned i=0;i<1000;++i){
        U64 n=1000000000001ULL + gen()%100000000000000ULL;
        auto [d,s]=discover_square_factor(n);
        if((__uint128_t)d*s*s!=n){
            std::cerr<<"bad partial factor n="<<n<<"\n";return 1;
        }
    }
    // Check the fast and unreduced paths compute exactly the same minimal Pell
    // unit on a representative range containing many squareful inputs.
    for(U64 d=2; d<=2000;++d){
        if(isqrt64(d)*isqrt64(d)==d)continue;
        auto [dd,scale]=discover_square_factor(d);
        Pell baseline=direct(zz(d));
        if(scale>1){
            auto reduced=reduced_solve(zz(d),zz(scale));
            Unit a=unit_power({reduced.base.X,reduced.base.Y},reduced.exponent,reduced.d);
            mpz_divexact_ui(a.v.get_mpz_t(),a.v.get_mpz_t(),(unsigned long)scale);
            if(a.u!=baseline.X || a.v!=baseline.Y){
                std::cerr<<"Pell solution discrepancy at D="<<d<<"\n";return 1;
            }
        }
    }
    std::cout<<"Square-factor tests passed (100k exhaustive plus random and structured).\n";
}
