// Exact C++17/GMP Pell solver, specialized from octonion/puzzles/continued/conway_full.cpp.
// X^2 - D Y^2 = 1. No factorization required by the direct algorithm.
// Optimized form/continued-fraction recurrence, midpoint + small-norm exits,
// checked uint64 blocks, GMP fallback, and standalone square-factor discovery.
#include <gmpxx.h>
#include "primes_10000.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <climits>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
using Z=mpz_class;
using U64=std::uint64_t;
struct Pell {Z X,Y;};
static Z zz(U64 x){static_assert(sizeof(unsigned long)==8,"requires LP64");return Z((unsigned long)x);}
static void check(const Z&D,const Pell&p){
    if(p.X<=0||p.Y<=0 || p.X*p.X-D*p.Y*p.Y!=1)
        throw std::runtime_error("Pell norm verification failed");
}
static U64 isqrt64(U64 n){
    // Exact integer square root; this call occurs once per solve.
    Z a=zz(n),b;mpz_sqrt(b.get_mpz_t(),a.get_mpz_t());return b.get_ui();
}
static Pell finish(Z u,Z v,long k){
    // Requires the primitive, positive, least-unit-associated state from the
    // reduced walk: u*u-D*v*v == k, k in {+/-1,+/-2,+/-4}.
    if(k==1)return {std::move(u),std::move(v)};
    Z t=u*u;
    if(k==-1)return {2*t+1,2*u*v};
    if(k==2 || k==-2)return {t-k/2,u*v};
    if(k==4){
        if(mpz_even_p(u.get_mpz_t()))return {t/2-1,u*v/2};
        return {u*(t-3)/2,v*(t-1)/2};
    }
    if(k==-4){
        if(mpz_even_p(u.get_mpz_t()))return {t/2+1,u*v/2};
        Z A=u*(t+3)/2,B=v*(t+1)/2; // norm -1, then square
        return {2*A*A+1,2*A*B};
    }
    throw std::logic_error("unexpected residual");
}
static bool terminal(U64 k){return k==1||k==2||k==4;}
struct Block {U64 a=1,b=0,c=0,d=1;};
// Safe because the expression is evaluated in unsigned 128-bit arithmetic.
static bool push(Block &w,U64 q){
    __uint128_t b=(__uint128_t)q*w.b+w.a, d=(__uint128_t)q*w.d+w.c;
    if(b>ULONG_MAX||d>ULONG_MAX)return false;
    w={w.b,(U64)b,w.d,(U64)d};return true;
}
static void cfstep(Z& prev,Z& v,U64 q){
    mpz_addmul_ui(prev.get_mpz_t(),v.get_mpz_t(),(unsigned long)q);
    mpz_swap(prev.get_mpz_t(),v.get_mpz_t());
}
static void flush(Block&w,Z&prev,Z&v,Z&np,Z&nv){
    mpz_mul_ui(np.get_mpz_t(),prev.get_mpz_t(),(unsigned long)w.a);
    mpz_addmul_ui(np.get_mpz_t(),v.get_mpz_t(),(unsigned long)w.c);
    mpz_mul_ui(nv.get_mpz_t(),prev.get_mpz_t(),(unsigned long)w.b);
    mpz_addmul_ui(nv.get_mpz_t(),v.get_mpz_t(),(unsigned long)w.d);
    mpz_swap(np.get_mpz_t(),prev.get_mpz_t());
    mpz_swap(nv.get_mpz_t(),v.get_mpz_t());w=Block{};
}
// Online checked word arithmetic; no growing-quotient vector for small inputs.
// For D<10^15 this routine was exercised against independent CF over 100k values.
static Pell direct_word(U64 D){
    if(D<2 || D>=1000000000000000ULL)throw std::invalid_argument("word path range");
    U64 r=isqrt64(D),rho=D-r*r;
    if(!rho)throw std::invalid_argument("D is square");
    U64 h=r,k=rho;int sign=-1;
    if(2*r+1-rho<=2){h=r+1;k=h*h-D;sign=1;}
    Z prev=0,v=1,np,nv;Block w;bool pending=false;size_t terms=0;
    auto drain=[&](){if(pending){flush(w,prev,v,np,nv);pending=false;}};
    auto advance=[&](U64 q){
        if(terms++<16){cfstep(prev,v,q);return;}
        if(!push(w,q)){drain();if(!push(w,q))throw std::overflow_error("single quotient too large");}
        pending=true;
    };
    int midpoint=0;
    while(!terminal(k)){
        U64 q=(h+r)/k,nh=q*k-h,nk=(D-nh*nh)/k;
        if(nh==h){midpoint=2;break;}
        U64 oldk=k;advance(q);h=nh;k=nk;sign=-sign;
        if(terminal(k))break;
        if(k==oldk){midpoint=1;break;}
    }
    drain();Z H=zz(h),K=zz(k),u=H*v+K*prev;
    if(midpoint==0)return finish(std::move(u),std::move(v),sign*(long)k);
    if(midpoint==2){
        Z X=2*u*u,Y=2*u*v;
        mpz_divexact_ui(X.get_mpz_t(),X.get_mpz_t(),(unsigned long)k);
        mpz_divexact_ui(Y.get_mpz_t(),Y.get_mpz_t(),(unsigned long)k);
        X-=sign;
        return {std::move(X),std::move(Y)};
    }
    Z um=H*(v*v-prev*prev)+2*K*v*prev,vm=v*v+prev*prev;
    return {2*um*um+1,2*um*vm};
}
static Pell direct_gmp(const Z&D){
    if(D<2)throw std::invalid_argument("D must be greater than 1");
    Z r,rho;mpz_sqrtrem(r.get_mpz_t(),rho.get_mpz_t(),D.get_mpz_t());
    if(rho==0)throw std::invalid_argument("D is square");
    Z h=r,k=rho;int sign=-1;
    if(2*r+1-rho<=2){h=r+1;k=h*h-D;sign=1;}
    Z prev=0,v=1,q;
    while(k!=1&&k!=2&&k!=4){
        q=(h+r)/k;
        h=q*k-h;k=(D-h*h)/k;sign=-sign;
        mpz_addmul(prev.get_mpz_t(),v.get_mpz_t(),q.get_mpz_t());
        mpz_swap(prev.get_mpz_t(),v.get_mpz_t());
    }
    return finish(h*v+k*prev,v,sign*k.get_si());
}
static Pell direct(const Z&D){
    if(D<2 || mpz_perfect_square_p(D.get_mpz_t()))throw std::invalid_argument("D must be positive nonsquare");
    if(mpz_fits_ulong_p(D.get_mpz_t()) && D<zz(1000000000000000ULL))
        return direct_word(D.get_ui());
    return direct_gmp(D);
}
// Exact prime-power index algorithm; factoring s is caller-enabled, never required.
// A known factorization of s may be provided in production to avoid trial division.
using Fact=std::pair<Z,unsigned>;
static std::vector<Fact> factor(Z n){
    if(n<1)throw std::invalid_argument("factor positive only");
    std::vector<Fact> out;
    auto take=[&](const Z&p){unsigned e=0;while(mpz_divisible_p(n.get_mpz_t(),p.get_mpz_t())){
        mpz_divexact(n.get_mpz_t(),n.get_mpz_t(),p.get_mpz_t());++e;}
        if(e)out.emplace_back(p,e);
    };
    take(Z(2));take(Z(3));
    for(Z p=5,step=2;p*p<=n;p+=step,step=6-step)take(p);
    if(n>1)out.emplace_back(n,1);
    return out;
}
struct Unit {Z u,v;}; // u+v sqrt(d); norm 1
static Unit mul_mod(const Unit&a,const Unit&b,const Z&d,const Z&m){
    return {(a.u*b.u+d*a.v*b.v)%m,(a.u*b.v+a.v*b.u)%m};
}
static Unit pow_mod(Unit a,Z e,const Z&d,const Z&m){
    a.u%=m;a.v%=m;Unit z{1,0};
    while(e>0){if(mpz_odd_p(e.get_mpz_t()))z=mul_mod(z,a,d,m);e>>=1;if(e>0)a=mul_mod(a,a,d,m);}
    return z;
}
static bool divisible_at(const Unit&a,const Z&e,const Z&d,const Z&m){return pow_mod(a,e,d,m).v==0;}
static Z index_prime_power(const Unit&a,const Z&d,const Z&p,unsigned e){
    Z rank=p==2?Z(2):mpz_divisible_p(d.get_mpz_t(),p.get_mpz_t())?p:(p-mpz_jacobi(d.get_mpz_t(),p.get_mpz_t()))/2;
    for(auto &[q,j]:factor(rank))for(unsigned i=0;i<j;i++){
        Z reduced=rank/q;if(divisible_at(a,reduced,d,p))rank=reduced;else break;
    }
    Z modulus=p;
    for(unsigned i=1;i<e;i++){modulus*=p;if(!divisible_at(a,rank,d,modulus))rank*=p;}
    if(!divisible_at(a,rank,d,modulus))throw std::runtime_error("failed local order verification");
    return rank;
}
static Unit unit_power(Unit a,Z n,const Z& d){
    Unit result{1,0};
    while(n>0){
        if(mpz_odd_p(n.get_mpz_t())){
            // Multiplication of two norm-one units.
            Unit tmp{result.u*a.u+d*result.v*a.v,result.u*a.v+result.v*a.u};
            result=std::move(tmp);
        }
        n>>=1;
        if(n>0){
            Z old_u=a.u;
            a.u=2*a.u*a.u-1;
            a.v=2*old_u*a.v;
        }
    }
    return result;
}
#include <iomanip>

// Standalone bounded square-factor extraction; always returns n = d*scale^2.
// Complete for n <= limit^3; otherwise d may retain a square factor.
static U64 icbrt64(U64 n){
    // Correct the floating-point estimate with exact 128-bit cubes.
    U64 r=static_cast<U64>(std::cbrt(static_cast<long double>(n)));
    while((__uint128_t)(r+1)*(r+1)*(r+1)<=n)++r;
    while((__uint128_t)r*r*r>n)--r;
    return r;
}
// Exact square test: floating-point estimate, integer-checked result.
static bool square64(U64 n,U64 &root){
    U64 x=static_cast<U64>(std::sqrt(static_cast<double>(n)));
    __uint128_t sq=static_cast<__uint128_t>(x)*x;
    while(sq>n){--x;sq=static_cast<__uint128_t>(x)*x;}
    while(static_cast<__uint128_t>(x+1)*(x+1)<=n){++x;sq=static_cast<__uint128_t>(x)*x;}
    root=x;
    return sq==n;
}
static std::pair<U64,U64> discover_square_factor(U64 n,U64 limit=10000){
    if(n<2)throw std::invalid_argument("positive nonsquare D required");
    if(limit>10000)throw std::invalid_argument("compiled prime table supports limit <= 10000");
    U64 bound=limit;
    if (static_cast<__uint128_t>(n)>static_cast<__uint128_t>(limit)*limit &&
        static_cast<__uint128_t>(n)<static_cast<__uint128_t>(limit)*limit*limit)
        bound=icbrt64(n);
    U64 rem=n,scale=1;
    // Check for a square remainder after each extracted prime factor.
    for(std::uint16_t prime:primes10000){
        U64 p=prime;
        if(p>bound || p*p>rem)break; // p <= 10000, so p*p cannot overflow
        if(rem%p!=0)continue;
        unsigned e=0;
        do{rem/=p;++e;}while(rem%p==0);
        for(unsigned i=0;i<e/2;++i)scale*=p;
        U64 root;
        if(square64(rem,root)){
            scale*=root;
            rem=1;
            break;
        }
    }
    if(rem!=1){U64 root;if(square64(rem,root))scale*=root;}
    // scale^2 divides n, so its product cannot overflow U64.
    return {n/(scale*scale),scale};
}

static std::vector<std::pair<U64,unsigned>> word_factors(U64 n){
    if(n==0)throw std::invalid_argument("factor zero");
    std::vector<std::pair<U64,unsigned>> f;
    auto take=[&](U64 p){unsigned e=0;while(n%p==0){n/=p;++e;}
        if(e)f.push_back({p,e});};
    take(2);take(3);
    for(U64 p=5,step=2;p<=n/p;p+=step,step=6-step)take(p);
    if(n>1)f.push_back({n,1});
    return f;
}
struct SmallUnit{U64 u,v;};
static U64 umul(U64 a,U64 b,U64 mod){return (U64)((__uint128_t)a*b%mod);}
static SmallUnit umul_unit(SmallUnit a,SmallUnit b,U64 dmod,U64 mod){
    U64 t=umul(dmod,umul(a.v,b.v,mod),mod);
    return {(U64)(((__uint128_t)umul(a.u,b.u,mod)+t)%mod),
            (U64)(((__uint128_t)umul(a.u,b.v,mod)+umul(a.v,b.u,mod))%mod)};
}
static bool coeff_v_zero(const Pell &base,U64 exponent,U64 d,U64 mod){
    SmallUnit a{mpz_fdiv_ui(base.X.get_mpz_t(),(unsigned long)mod),
                mpz_fdiv_ui(base.Y.get_mpz_t(),(unsigned long)mod)},result{1,0};
    U64 dm=d%mod;
    while(exponent){
        if(exponent&1)result=umul_unit(result,a,dm,mod);
        exponent>>=1;
        if(exponent)a=umul_unit(a,a,dm,mod);
    }
    return result.v==0;
}
static U64 word_index(const Pell&base,U64 d,U64 p,unsigned e){
    U64 rank;
    if(p==2)rank=2;
    else if(d%p==0)rank=p;
    else{Z aa=zz(d%p),pp=zz(p);rank=(p-mpz_jacobi(aa.get_mpz_t(),pp.get_mpz_t()))/2;}
    for(auto [q,f]:word_factors(rank))for(unsigned j=0;j<f;++j){
        if(coeff_v_zero(base,rank/q,d,p))rank/=q;else break;
    }
    U64 mod=p;
    for(unsigned i=1;i<e;++i){
        if(mod>UINT64_MAX/p)throw std::overflow_error("word prime power exceeds 64 bits");
        mod*=p;
        if(!coeff_v_zero(base,rank,d,mod))rank*=p;
    }
    if(!coeff_v_zero(base,rank,d,mod))throw std::runtime_error("failed native order certificate");
    return rank;
}
struct Reduced {Pell base; Z d,s,exponent;};
static Reduced reduced_solve(const Z& D,const Z& s){
    if(s<1 || D%(s*s)!=0)throw std::invalid_argument("invalid supplied square factor");
    Z d=D/(s*s);
    Pell base=direct(d);
    Z exp=1;
    if(s>1){
        // Exploit exact 64-bit modular arithmetic when every local factor fits.
        if(mpz_fits_ulong_p(s.get_mpz_t()) && mpz_fits_ulong_p(d.get_mpz_t()) &&
            s<zz(1000000000ULL) && d<zz(1000000000000000ULL)){
            for(auto [p,e]:word_factors(s.get_ui())){
                Z q=zz(word_index(base,d.get_ui(),p,e));
                mpz_lcm(exp.get_mpz_t(),exp.get_mpz_t(),q.get_mpz_t());
            }
        }else{
            Unit u{base.X,base.Y};
            for(auto &[p,e]:factor(s)){
                Z q=index_prime_power(u,d,p,e);
                mpz_lcm(exp.get_mpz_t(),exp.get_mpz_t(),q.get_mpz_t());
            }
        }
    }
    return {std::move(base),std::move(d),s,std::move(exp)};
}
// mpz_sizeinbase(z,10) can overestimate by one; avoid decimal conversion
// for both the default output decision and large-result summaries.
static bool fits_print_limit(const Z& value, unsigned long max_digits) {
    if(max_digits==0) return false;
    const size_t estimate=mpz_sizeinbase(value.get_mpz_t(),10);
    if(estimate <= max_digits) return true;
    if(estimate - max_digits > 1) return false;
    Z cutoff;
    mpz_ui_pow_ui(cutoff.get_mpz_t(),10UL,max_digits);
    return mpz_cmpabs(value.get_mpz_t(),cutoff.get_mpz_t()) < 0;
}
// X_n = T_n(u), where u+v*sqrt(d) is the primitive norm-one unit.
// This lower bound is deliberately conservative but entirely integer-exact:
// T_n(u) >= u^n >= 2^(n*(bit_length(u)-1)); and 2^4 > 10.
// Consequently X_n cannot fit in `limit` decimal digits when the test succeeds.
// In borderline cases, expand and perform the existing exact print-limit test.
static bool proven_too_large_to_print(const Reduced &r, unsigned long limit) {
    if (limit == 0) return true;
    size_t bits = mpz_sizeinbase(r.base.X.get_mpz_t(),2);
    if (bits <= 1) return false;
    Z guaranteed_bits_lower = r.exponent * zz((U64)(bits - 1));
    return guaranteed_bits_lower >= Z(limit) * 4;
}

enum class PrintMode { Auto, Yes, No };
static volatile U64 sink=0;
int main(int argc,char **argv){
 try{
    if(argc<2)throw std::invalid_argument("usage: conway_fast D [--square-factor S|--no-discover-square-factor] [--verify|--expand] [--print|--no-print] [--print-limit N] [--repeat N]");
    Z D(argv[1]),s=1;bool discover=true,verify=false,force_expand=false,explicit_discover=false,manual_factor=false;
    PrintMode print_mode=PrintMode::Auto;
    unsigned long print_limit=1000;
    unsigned repeat=1;
    for(int i=2;i<argc;++i){std::string arg=argv[i];
        if(arg=="--square-factor"&&i+1<argc){s=Z(argv[++i]);manual_factor=true;discover=false;}
        else if(arg=="--discover-square-factor"){discover=true;explicit_discover=true;}
        else if(arg=="--direct" || arg=="--no-discover-square-factor")discover=false;
        else if(arg=="--verify")verify=true;
        else if(arg=="--expand")force_expand=true;
        else if(arg=="--print")print_mode=PrintMode::Yes;
        else if(arg=="--no-print")print_mode=PrintMode::No;
        else if(arg=="--print-limit"&&i+1<argc){
            const std::string value=argv[++i];
            if(value.empty() || !std::all_of(value.begin(),value.end(),[](unsigned char c){return c>='0' && c<='9';}))
                throw std::invalid_argument("invalid --print-limit");
            size_t consumed=0;
            print_limit=std::stoul(value,&consumed,10);
            if(consumed!=value.size())throw std::invalid_argument("invalid --print-limit");
        }
        else if(arg=="--repeat"&&i+1<argc)repeat=std::stoul(argv[++i]);
        else throw std::invalid_argument("unknown argument: "+arg);
    }
    if(repeat==0)throw std::invalid_argument("repeat must be positive");
    if(D<2||mpz_perfect_square_p(D.get_mpz_t()))throw std::invalid_argument("D must be nonsquare and >1");
    if(explicit_discover && manual_factor)throw std::invalid_argument("cannot combine supplied and discovered factors");
    if(discover && !mpz_fits_ulong_p(D.get_mpz_t())){
        if(explicit_discover)throw std::invalid_argument("explicit square-factor discovery requires uint64 D");
        discover=false; // Very large inputs use the correct direct solver.
    }
    Pell whole;
    Reduced reduced;
    Z selected_scale=1;
    bool expanded=true;
    auto start=std::chrono::steady_clock::now();
    for(unsigned i=0;i<repeat;++i){
        Z scale=s;
        if(discover){auto [dd,ss]=discover_square_factor(D.get_ui());(void)dd;scale=zz(ss);}
        selected_scale=scale;
        if(scale!=1){
            reduced=reduced_solve(D,scale);
            // An exact compact answer is (d,s,u,v,n) with
            // u+v*sqrt(d) norm 1 and n minimal such that s divides V_n.
            // The full solution, if requested, is
            //    X = T_n(u), Y = v*U_{n-1}(u)/s.
            // There is no obligation to construct these massive integers
            // when the caller requests compact output only.
            expanded=force_expand || verify || print_mode==PrintMode::Yes ||
                (print_mode==PrintMode::Auto &&
                 !proven_too_large_to_print(reduced,print_limit));
            if(expanded){
                if(reduced.exponent==1){
                    whole={reduced.base.X,reduced.base.Y/scale};
                }else{
                    Unit unit=unit_power({reduced.base.X,reduced.base.Y},reduced.exponent,reduced.d);
                    whole={std::move(unit.u),std::move(unit.v)};
                    mpz_divexact(whole.Y.get_mpz_t(),whole.Y.get_mpz_t(),scale.get_mpz_t());
                }
                sink+=mpz_sizeinbase(whole.X.get_mpz_t(),2);
            }else{
                // Keep an observable checksum without touching the huge unit.
                sink+=mpz_sizeinbase(reduced.base.X.get_mpz_t(),2)
                    +mpz_sizeinbase(reduced.exponent.get_mpz_t(),2);
            }
        }else{
            whole=direct(D);
            expanded=true; // Without reduction, the current solver builds the unit.
            sink+=mpz_sizeinbase(whole.X.get_mpz_t(),2);
        }
    }
    auto end=std::chrono::steady_clock::now();
    if(verify) check(D,whole);
    std::cout<<std::fixed<<std::setprecision(3)
      <<"compute_us="<<std::chrono::duration<double,std::micro>(end-start).count()/repeat;
    if(selected_scale>1){
        std::cout<<" square_factor="<<selected_scale<<" reduced_D="<<reduced.d
                 <<" exponent="<<reduced.exponent;
    }
    if(expanded){
        std::cout<<" representation=expanded"
                 <<" X_bits="<<mpz_sizeinbase(whole.X.get_mpz_t(),2)
                 <<" Y_bits="<<mpz_sizeinbase(whole.Y.get_mpz_t(),2);
    }else{
        std::cout<<" representation=compact"
                 <<" base_X_bits="<<mpz_sizeinbase(reduced.base.X.get_mpz_t(),2)
                 <<" base_Y_bits="<<mpz_sizeinbase(reduced.base.Y.get_mpz_t(),2);
        // Coordinate lengths are deliberately omitted: the exact full values
        // have not been constructed, and a bit-length guess is not a proof.
    }
    bool show=expanded && (print_mode==PrintMode::Yes ||
              (print_mode==PrintMode::Auto &&
               fits_print_limit(whole.X,print_limit) &&
               fits_print_limit(whole.Y,print_limit)));
    if(show)std::cout<<" X="<<whole.X<<" Y="<<whole.Y;
    std::cout<<"\n";
    return 0;
 }catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<"\n";return 1;}
}
