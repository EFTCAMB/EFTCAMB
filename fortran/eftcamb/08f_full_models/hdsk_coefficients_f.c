#include <math.h> 

/** functions to compute coefficients of perturbation equations*/
void horndeski_perturbeq_coefficients(
    double a,
    double H,
    double dH,
    double ddH,
    double dphi,
    double ddphi,
    double dddphi,
    double rptot,
    double drptot,
    double G2,
    double G2p,
    double G2x,
    double G2pp,
    double G2px,
    double G2xx,
    double G2ppp,
    double G2ppx,
    double G2pxx,
    double G2xxx,
    double G3p,
    double G3x,
    double G3pp,
    double G3px,
    double G3xx,
    double G3ppp,
    double G3ppx,
    double G3pxx,
    double G3xxx,
    double G4,
    double G4p,
    double G4x,
    double G4pp,
    double G4px,
    double G4xx,
    double G4ppp,
    double G4ppx,
    double G4pxx,
    double G4xxx,
    double G4pppp,
    double G4pppx,
    double G4ppxx,
    double G4pxxx,
    double G4xxxx,
    double G5,
    double G5p,
    double G5x,
    double G5pp,
    double G5px,
    double G5xx,
    double G5ppp,
    double G5ppx,
    double G5pxx,
    double G5xxx,
    double G5pppp,
    double G5pppx,
    double G5ppxx,
    double G5pxxx,
    double G5xxxx,
    double *spi1,
    double *spi2,
    double *spi3,
    double *spih,
    double *spie,
    double *spim,
    double *s00,
    double *s00k,
    double *s00p,
    double *s0i,
    double *s0ip,
    double *sii,
    double *siik,
    double *siip,
    double *siipp,
    double *sij,
    double *sijdot)
{  
   // Milestone spi1
   *spi1 = (-12*pow(a,11)*pow(G4p,2)*H + ddphi*pow(dphi,9)*(3*G5xx*G5xxx - 2*G5x*G5xxxx)*pow(H,4) + 
   a*pow(dphi,8)*pow(H,3)*(-2*ddphi*(6*G4xxxx*G5x - 3*G5pxxx*G5x - 6*G4xxx*G5xx + 3*G5pxx*G5xx - 6*G4xx*G5xxx + 3*G5px*G5xxx + 2*G4x*G5xxxx - G5p*G5xxxx) + 
      pow(dphi,2)*(-3*G5xx*G5xxx + 2*G5x*G5xxxx)*pow(H,2)) + pow(a,2)*pow(dphi,7)*pow(H,2)*
    (dphi*H*(6*dH*(pow(G5xx,2) - G5x*G5xxx) + dphi*(12*G4xxxx*G5x - 8*G5pxxx*G5x - 12*G4xxx*G5xx + 9*G5pxx*G5xx - 12*G4xx*G5xxx + 6*G5px*G5xxx + 4*G4x*G5xxxx - 2*G5p*G5xxxx)*H) + 
      3*ddphi*(16*G4xx*G4xxx - 8*G4x*G4xxxx + 4*G4xxxx*G5p - 8*G4xxx*G5px - 8*G4xx*G5pxx + 4*G5px*G5pxx + 4*G4x*G5pxxx - 2*G5p*G5pxxx - 2*G3xxx*G5x + 4*G4pxxx*G5x + G3xx*G5xx - 2*G4pxx*G5xx + G3x*G5xxx - 
         2*G4px*G5xxx + 7*pow(G5xx,2)*pow(H,2) - 5*G5x*G5xxx*pow(H,2))) + pow(a,6)*pow(dphi,3)*
    (6*dH*dphi*H*(-16*G4x*G4xx + 8*G4*G4xxx + 12*G4x*G5px - 2*G5p*G5px - 4*G4*G5pxx + G3x*G5x - 2*G4p*G5xx + 3*pow(G5x,2)*pow(H,2)) + 
      2*ddphi*(3*pow(G3x,2) + 2*G2xxx*G4 - 2*G3pxx*G4 - 3*G3xx*G4p - 15*G3x*G4px + 18*pow(G4px,2) + 6*G4p*G4pxx - 6*G2xx*G4x + 8*G3px*G4x + 3*G2xx*G5p - 4*G3px*G5p + 12*G4x*G4xx*pow(H,2) + 
         96*G4*G4xxx*pow(H,2) - 66*G4xx*G5p*pow(H,2) + 36*G5p*G5px*pow(H,2) - 54*G4*G5pxx*pow(H,2) + 12*G3x*G5x*pow(H,2) - 27*G4px*G5x*pow(H,2) - 21*G4p*G5xx*pow(H,2) + 
         21*pow(G5x,2)*pow(H,4)) - pow(dphi,2)*(6*G3x*G4ppx - 12*G4ppx*G4px + 4*G2pxx*G4x - 4*G3ppx*G4x - 2*G2pxx*G5p + 2*G3ppx*G5p + 12*G3xxx*G4*pow(H,2) - 48*G4*G4pxxx*pow(H,2) - 
         48*G3xx*G4x*pow(H,2) + 168*G4pxx*G4x*pow(H,2) + 126*G3x*G4xx*pow(H,2) - 372*G4px*G4xx*pow(H,2) - 24*G4p*G4xxx*pow(H,2) + 18*G3xx*G5p*pow(H,2) - 48*G4pxx*G5p*pow(H,2) + 
         72*G4xx*G5pp*pow(H,2) - 36*G4x*G5ppx*pow(H,2) + 6*G5p*G5ppx*pow(H,2) + 12*G4*G5ppxx*pow(H,2) - 81*G3x*G5px*pow(H,2) + 222*G4px*G5px*pow(H,2) - 36*G5pp*G5px*pow(H,2) + 
         18*G4p*G5pxx*pow(H,2) - 6*G2xx*G5x*pow(H,2) - 18*G4ppx*G5x*pow(H,2) + 3*G2x*G5xx*pow(H,2) - 6*G3p*G5xx*pow(H,2) + 12*G4pp*G5xx*pow(H,2) + 198*G4xx*G5x*pow(H,4) - 
         123*G5px*G5x*pow(H,4) + 18*G4x*G5xx*pow(H,4) - 72*G5p*G5xx*pow(H,4) + 48*G4*G5xxx*pow(H,4) + G3px*(-3*G3x + 6*G4px + 11*G5x*pow(H,2)))) + 
   pow(a,7)*pow(dphi,2)*(12*dH*dphi*(G3xx*G4 - 2*G4*G4pxx - G3x*G4x + 4*G4px*G4x - 2*G4p*G4xx - G4px*G5p + G4p*G5px + G4x*G5x*pow(H,2) - 4*G5p*G5x*pow(H,2) + 7*G4*G5xx*pow(H,2)) + 
      12*ddphi*H*(5*G3xx*G4 - 12*G4*G4pxx + G3x*G4x - 2*G4px*G4x - 8*G4p*G4xx - 2*G3x*G5p + 5*G4px*G5p + 5*G4p*G5px + 7*G4x*G5x*pow(H,2) - 8*G5p*G5x*pow(H,2) + 9*G4*G5xx*pow(H,2)) - 
      pow(dphi,2)*H*(15*pow(G3x,2) + 4*G2xxx*G4 - 16*G3pxx*G4 - 6*G3xx*G4p + 24*G4*G4ppxx + 120*pow(G4px,2) + 36*G4p*G4pxx - 12*G2xx*G4x + 28*G3px*G4x - 48*G4ppx*G4x + 12*G2x*G4xx - 24*G3p*G4xx + 
         48*G4pp*G4xx + 6*G2xx*G5p - 8*G3px*G5p + 12*G4ppx*G5p - 36*G4px*G5pp - 12*G4p*G5ppx - 6*G2x*G5px + 12*G3p*G5px - 24*G4pp*G5px + 2*G2px*G5x - 4*G3pp*G5x + 168*G4x*G4xx*pow(H,2) + 
         192*G4*G4xxx*pow(H,2) - 276*G4xx*G5p*pow(H,2) - 96*G4x*G5px*pow(H,2) + 180*G5p*G5px*pow(H,2) - 136*G4*G5pxx*pow(H,2) - 186*G4px*G5x*pow(H,2) + 42*G5pp*G5x*pow(H,2) - 
         54*G4p*G5xx*pow(H,2) + 69*pow(G5x,2)*pow(H,4) + G3x*(-90*G4px + 18*G5pp + 60*G5x*pow(H,2)))) + 
   pow(a,3)*pow(dphi,6)*H*(2*ddphi*(-6*G3xxx*G4x + 12*G4pxxx*G4x + 6*G3xx*G4xx - 12*G4pxx*G4xx + 6*G3x*G4xxx - 12*G4px*G4xxx + 3*G3xxx*G5p - 6*G4pxxx*G5p - 3*G3xx*G5px + 6*G4pxx*G5px - 3*G3x*G5pxx + 
         6*G4px*G5pxx - G2xxx*G5x + G3pxx*G5x - 30*G4xxx*G5x*pow(H,2) + 18*G5pxx*G5x*pow(H,2) + 66*G4xx*G5xx*pow(H,2) - 36*G5px*G5xx*pow(H,2) - 18*G4x*G5xxx*pow(H,2) + 6*G5p*G5xxx*pow(H,2) + 
         2*G4*G5xxxx*pow(H,2)) + dphi*H*(-6*dH*(4*G4xxx*G5x - 2*G5pxx*G5x - 6*G4xx*G5xx + 3*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx) + 
         dphi*H*(-48*G4xx*G4xxx + 24*G4x*G4xxxx - 12*G4xxxx*G5p + 24*G4xxx*G5px + 36*G4xx*G5pxx - 18*G5px*G5pxx - 16*G4x*G5pxxx + 8*G5p*G5pxxx + 6*G3xxx*G5x - 24*G4pxxx*G5x + 6*G5ppxx*G5x - 3*G3xx*G5xx + 
            18*G4pxx*G5xx - 6*G5ppx*G5xx - 3*G3x*G5xxx + 6*G4px*G5xxx - 24*pow(G5xx,2)*pow(H,2) + 15*G5x*G5xxx*pow(H,2)))) + 
   pow(a,4)*pow(dphi,5)*(ddphi*(-6*G3xx*G4px + 12*G4px*G4pxx - 4*G2xxx*G4x + 4*G3pxx*G4x + 2*G2xxx*G5p - 2*G3pxx*G5p + 192*pow(G4xx,2)*pow(H,2) - 144*G4x*G4xxx*pow(H,2) + 24*G4*G4xxxx*pow(H,2) + 
         48*G4xxx*G5p*pow(H,2) - 216*G4xx*G5px*pow(H,2) + 60*pow(G5px,2)*pow(H,2) + 84*G4x*G5pxx*pow(H,2) - 30*G5p*G5pxx*pow(H,2) - 12*G4*G5pxxx*pow(H,2) - 21*G3xx*G5x*pow(H,2) + 
         54*G4pxx*G5x*pow(H,2) - 60*G4px*G5xx*pow(H,2) - 6*G4p*G5xxx*pow(H,2) + 27*G5x*G5xx*pow(H,4) + 3*G3x*(G3xx - 2*G4pxx + 9*G5xx*pow(H,2))) + 
      dphi*H*(6*dH*(8*pow(G4xx,2) - 8*G4x*G4xxx + 4*G4xxx*G5p - 8*G4xx*G5px + 2*pow(G5px,2) + 4*G4x*G5pxx - 2*G5p*G5pxx - G3xx*G5x + 2*G4pxx*G5x + G3x*G5xx - 2*G4px*G5xx - G5x*G5xx*pow(H,2)) + 
         dphi*H*(12*G3xxx*G4x - 48*G4pxxx*G4x - 12*G3xx*G4xx + 72*G4pxx*G4xx - 12*G3x*G4xxx + 24*G4px*G4xxx - 6*G3xxx*G5p + 24*G4pxxx*G5p - 24*G4xx*G5ppx + 12*G4x*G5ppxx - 6*G5p*G5ppxx + 6*G3xx*G5px - 
            36*G4pxx*G5px + 12*G5ppx*G5px + 9*G3x*G5pxx - 18*G4px*G5pxx + 2*G2xxx*G5x - 8*G3pxx*G5x + 12*G4ppxx*G5x + 3*G3px*G5xx - 6*G4ppx*G5xx + 60*G4xxx*G5x*pow(H,2) - 41*G5pxx*G5x*pow(H,2) - 
            162*G4xx*G5xx*pow(H,2) + 99*G5px*G5xx*pow(H,2) + 36*G4x*G5xxx*pow(H,2) - 12*G5p*G5xxx*pow(H,2) - 4*G4*G5xxxx*pow(H,2)))) + 
   4*pow(a,9)*(6*ddphi*H*(G3x*G4 - 3*G4*G4px - G4p*G4x + G4p*G5p + G4*G5x*pow(H,2)) + 6*dH*dphi*(G3x*G4 - 3*G4*G4px - G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) + 
      pow(dphi,2)*H*(-3*G2xx*G4 + 10*G3px*G4 + 3*(3*G3x*G4p - 11*G4p*G4px - G2x*G4x + 2*G3p*G4x - 4*G4pp*G4x + G2x*G5p - 2*G3p*G5p + 4*G4pp*G5p + 3*G4p*G5pp - 10*pow(G4x,2)*pow(H,2) + 
            20*G4x*G5p*pow(H,2) - 10*pow(G5p,2)*pow(H,2) + 6*G4p*G5x*pow(H,2) - 2*G4*(3*G4ppx + (9*G4xx - 7*G5px)*pow(H,2))))) + 
   2*pow(a,10)*(24*dH*G4*(G4x - G5p)*H + dphi*(2*G2px*G4 - 4*G3pp*G4 + 3*G2x*G4p - 6*(G3p*G4p - 2*G4p*G4pp + 7*G4p*(-G4x + G5p)*pow(H,2) + 2*G4*pow(H,2)*(G3x - 4*G4px + G5pp + G5x*pow(H,2))))) - 
   pow(a,8)*dphi*(12*dH*dphi*(-16*G4*G4xx + 2*G4x*G5p - 2*pow(G5p,2) + 10*G4*G5px + 3*G4p*G5x)*H - 
      4*ddphi*(3*G2xx*G4 - 4*G3px*G4 - 3*G3x*G4p + 9*G4p*G4px + 12*pow(G4x,2)*pow(H,2) + 54*G4*G4xx*pow(H,2) - 24*G4x*G5p*pow(H,2) + 12*pow(G5p,2)*pow(H,2) - 36*G4*G5px*pow(H,2) - 
         9*G4p*G5x*pow(H,2)) + pow(dphi,2)*(3*G2x*(G3x - 2*G4px + 3*G5x*pow(H,2)) + 
         2*(-2*G2pxx*G4 + 2*G3ppx*G4 + 3*G3px*G4p + 6*G3x*G4pp - 6*G4p*G4ppx - 12*G4pp*G4px + 2*G2px*G4x - 4*G3pp*G4x - G2px*G5p + 2*G3pp*G5p + 30*G3xx*G4*pow(H,2) - 120*G4*G4pxx*pow(H,2) + 
            33*G3x*G4x*pow(H,2) - 102*G4px*G4x*pow(H,2) - 78*G4p*G4xx*pow(H,2) - 39*G3x*G5p*pow(H,2) + 126*G4px*G5p*pow(H,2) + 24*G4x*G5pp*pow(H,2) - 30*G5p*G5pp*pow(H,2) + 
            30*G4*G5ppx*pow(H,2) + 57*G4p*G5px*pow(H,2) + 18*G4pp*G5x*pow(H,2) + 87*G4x*G5x*pow(H,4) - 93*G5p*G5x*pow(H,4) + 54*G4*G5xx*pow(H,4) - 3*G3p*(G3x - 2*G4px + 3*G5x*pow(H,2))))) + 
   pow(a,5)*pow(dphi,4)*(-6*dH*dphi*(2*G3xx*G4x - 4*G4pxx*G4x - 2*G3x*G4xx + 4*G4px*G4xx - G3xx*G5p + 2*G4pxx*G5p + G3x*G5px - 2*G4px*G5px - 2*G4xx*G5x*pow(H,2) - G5px*G5x*pow(H,2) + 
         8*G4x*G5xx*pow(H,2) - G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) - H*
       (pow(dphi,2)*(36*G4px*G4pxx - 4*G2xxx*G4x + 16*G3pxx*G4x - 24*G4ppxx*G4x - 12*G3px*G4xx + 24*G4ppx*G4xx + 2*G2xxx*G5p - 8*G3pxx*G5p + 12*G4ppxx*G5p - 12*G4px*G5ppx + 6*G3px*G5px - 12*G4ppx*G5px + 
            2*G2pxx*G5x - 2*G3ppx*G5x + 264*pow(G4xx,2)*pow(H,2) - 144*G4x*G4xxx*pow(H,2) + 24*G4*G4xxxx*pow(H,2) + 48*G4xxx*G5p*pow(H,2) - 336*G4xx*G5px*pow(H,2) + 
            102*pow(G5px,2)*pow(H,2) + 100*G4x*G5pxx*pow(H,2) - 32*G5p*G5pxx*pow(H,2) - 16*G4*G5pxxx*pow(H,2) + 66*G4pxx*G5x*pow(H,2) - 12*G5ppx*G5x*pow(H,2) - 108*G4px*G5xx*pow(H,2) + 
            18*G5pp*G5xx*pow(H,2) - 6*G4p*G5xxx*pow(H,2) + 45*G5x*G5xx*pow(H,4) - 3*G3xx*(2*G4px + 7*G5x*pow(H,2)) + 3*G3x*(G3xx - 6*G4pxx + 2*G5ppx + 13*G5xx*pow(H,2))) + 
         2*ddphi*(-6*G3xxx*G4 + 24*G3xx*G4x - 60*G4pxx*G4x - 36*G3x*G4xx + 84*G4px*G4xx + 12*G4p*G4xxx - 9*G3xx*G5p + 24*G4pxx*G5p + 21*G3x*G5px - 48*G4px*G5px - 6*G4p*G5pxx + 3*G2xx*G5x - 4*G3px*G5x - 
            54*G4xx*G5x*pow(H,2) + 27*G5px*G5x*pow(H,2) + 6*G4x*G5xx*pow(H,2) + 21*G5p*G5xx*pow(H,2) + 12*G4*(G4pxxx - 2*G5xxx*pow(H,2))))))/
 (2.*pow(a,3)*H*(-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
     4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
     2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
        18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
     2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
        14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
        3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
      (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
        2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*
      (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
        18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))));

   // Milestone spi2
   *spi2 = (-4*pow(a,9)*pow(G4p,2) + 2*ddphi*pow(dphi,7)*(pow(G5xx,2) - G5x*G5xxx)*pow(H,3) + 
   a*(2*ddphi*pow(dphi,6)*(-4*G4xxx*G5x + 2*G5pxx*G5x + 6*G4xx*G5xx - 3*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2) + pow(dphi,8)*(-5*pow(G5xx,2) + 4*G5x*G5xxx)*pow(H,4)) + 
   2*pow(a,2)*pow(dphi,5)*H*(ddphi*(8*pow(G4xx,2) - 8*G4x*G4xxx + 4*G4xxx*G5p - 8*G4xx*G5px + 2*pow(G5px,2) + 4*G4x*G5pxx - 2*G5p*G5pxx - G3xx*G5x + 2*G4pxx*G5x + G3x*G5xx - 2*G4px*G5xx) + 
      dphi*H*(-(dH*G5x*G5xx) + 2*dphi*(5*G4xxx*G5x - 3*G5pxx*G5x - 9*G4xx*G5xx + 5*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*H)) + 
   16*pow(a,8)*(dH*G4*(G4x - G5p) - dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))) - 
   4*pow(a,7)*(-4*dH*dphi*G4*G5x*H - 2*ddphi*(G3x*G4 - 3*G4*G4px - G4p*G4x + G4p*G5p + G4*G5x*pow(H,2)) + 
      pow(dphi,2)*(G2xx*G4 - 2*G3px*G4 - 2*G3x*G4p + 2*G4*G4ppx + 6*G4p*G4px - G4p*G5pp + 12*pow(G4x,2)*pow(H,2) + 26*G4*G4xx*pow(H,2) - 24*G4x*G5p*pow(H,2) + 12*pow(G5p,2)*pow(H,2) - 
         18*G4*G5px*pow(H,2) - 7*G4p*G5x*pow(H,2))) + 2*pow(a,4)*pow(dphi,3)*
    (dH*dphi*(-8*G4x*G4xx + 4*G4xx*G5p + 4*G4x*G5px - 2*G5p*G5px + G3x*G5x - 2*G4px*G5x - pow(G5x,2)*pow(H,2)) + 
      2*ddphi*H*(-4*G4x*G4xx + 4*G4*G4xxx - 2*G4xx*G5p + 4*G4x*G5px - 2*G4*G5pxx + G4px*G5x - G4p*G5xx + 2*pow(G5x,2)*pow(H,2)) + 
      pow(dphi,2)*H*(-14*G3x*G4xx + 36*G4px*G4xx + G3xx*(8*G4x - 4*G5p) + 12*G4pxx*(-2*G4x + G5p) - 4*G4xx*G5pp + 4*G4x*G5ppx - 2*G5p*G5ppx + 8*G3x*G5px - 20*G4px*G5px + 2*G5pp*G5px + G2xx*G5x - 
         2*G3px*G5x + 2*G4ppx*G5x - 20*G4xx*G5x*pow(H,2) + 8*G5px*G5x*pow(H,2) + 4*G4x*G5xx*pow(H,2) + 6*G5p*G5xx*pow(H,2) - 4*G4*G5xxx*pow(H,2))) - 
   pow(a,5)*pow(dphi,2)*(8*dH*dphi*(2*G4x*G5x - G5p*G5x - G4*G5xx)*H + 4*ddphi*
       (-(G3xx*G4) + 2*G4*G4pxx + G3x*G4x - 4*G4px*G4x + 2*G4p*G4xx + G4px*G5p - G4p*G5px - 5*G4x*G5x*pow(H,2) + 6*G5p*G5x*pow(H,2) - 5*G4*G5xx*pow(H,2)) + 
      pow(dphi,2)*(3*pow(G3x,2) + 20*pow(G4px,2) - 4*G2xx*G4x + 8*G3px*G4x - 8*G4ppx*G4x + 2*G2xx*G5p - 4*G3px*G5p + 4*G4ppx*G5p + 8*G4x*G4xx*pow(H,2) + 40*G4*G4xxx*pow(H,2) - 
         60*G4xx*G5p*pow(H,2) + 8*G4x*G5px*pow(H,2) + 28*G5p*G5px*pow(H,2) - 24*G4*G5pxx*pow(H,2) + 6*G5pp*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 25*pow(G5x,2)*pow(H,4) + 
         2*G3x*(-8*G4px + G5pp + 6*G5x*pow(H,2)) - 4*G4px*(G5pp + 7*G5x*pow(H,2)))) - 
   2*pow(a,3)*pow(dphi,4)*(ddphi*(2*G3xx*G4x - 4*G4pxx*G4x - 2*G3x*G4xx + 4*G4px*G4xx - G3xx*G5p + 2*G4pxx*G5p + G3x*G5px - 2*G4px*G5px - 2*G4xx*G5x*pow(H,2) - G5px*G5x*pow(H,2) + 
         4*G4x*G5xx*pow(H,2) + G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + dphi*H*
       (dH*(4*G4x*G5xx - 2*G5p*G5xx) + dphi*H*(32*pow(G4xx,2) - 20*G4x*G4xxx + 10*G4xxx*G5p - 36*G4xx*G5px + 10*pow(G5px,2) + 12*G4x*G5pxx - 6*G5p*G5pxx - 4*G3xx*G5x + 12*G4pxx*G5x - 2*G5ppx*G5x + 
            4*G3x*G5xx - 10*G4px*G5xx + G5pp*G5xx + 3*G5x*G5xx*pow(H,2)))) - 4*pow(a,6)*dphi*
    (dH*dphi*(4*pow(G4x,2) - 4*G4*G4xx - 6*G4x*G5p + 2*pow(G5p,2) + 2*G4*G5px + G4p*G5x) + 
      2*H*(ddphi*(-2*pow(G4x,2) - 6*G4*G4xx + 4*G4x*G5p - 2*pow(G5p,2) + 4*G4*G5px + G4p*G5x) + 
         pow(dphi,2)*(2*G3xx*G4 + G3x*G4x - 2*G4px*G4x - 5*G4p*G4xx - 2*G3x*G5p + 5*G4px*G5p + G4x*G5pp - G5p*G5pp + 3*G4p*G5px + 8*G4x*G5x*pow(H,2) - 9*G5p*G5x*pow(H,2) + 
            G4*(-6*G4pxx + G5ppx + 5*G5xx*pow(H,2))))))/
 (a*(-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
     4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
     2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
        18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
     2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
        14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
        3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
      (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
        2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*
      (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
        18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))));

   // Milestone spi3
   *spi3 = (a*pow(dphi,16)*(-4*G5pxx*G5pxxx*pow(G5x,2) - 6*pow(G5pxx,2)*G5x*G5xx + 36*G5px*G5pxxx*G5x*G5xx - 36*G4xxx*G5px*pow(G5xx,2) - 12*G4x*G5pxxx*pow(G5xx,2) + 6*G5p*G5pxxx*pow(G5xx,2) - 
   36*G4pxxx*G5x*pow(G5xx,2) + 12*G5ppxx*G5x*pow(G5xx,2) + 36*G4pxx*pow(G5xx,3) - 9*G5ppx*pow(G5xx,3) + 4*G5ppxx*pow(G5x,2)*G5xxx + 36*pow(G5px,2)*G5xx*G5xxx + 24*G4pxx*G5x*G5xx*G5xxx - 
   18*G5ppx*G5x*G5xx*G5xxx - 12*G4px*pow(G5xx,2)*G5xxx + 6*G5pp*pow(G5xx,2)*G5xxx - 60*G4xx*G5xx*(G5pxxx*G5x + G5px*G5xxx) + 
   6*G5pxx*(6*G4xxx*G5x*G5xx + 21*G4xx*pow(G5xx,2) - 10*G5px*pow(G5xx,2) + 6*G4xx*G5x*G5xxx - 4*G5px*G5x*G5xxx + 2*G4x*G5xx*G5xxx - G5p*G5xx*G5xxx))*pow(H,8) + 
3*pow(dphi,17)*G5xx*(-2*G5pxxx*G5x*G5xx + 3*G5pxx*pow(G5xx,2) + 2*G5pxx*G5x*G5xxx - 2*G5px*G5xx*G5xxx)*pow(H,9) + 
pow(a,2)*pow(dphi,15)*pow(H,7)*(48*G5px*pow(G5pxx,2)*G5x - 48*pow(G5px,2)*G5pxxx*G5x - 16*G4x*G5pxx*G5pxxx*G5x + 8*G5p*G5pxx*G5pxxx*G5x - 24*G4pxxx*G5pxx*pow(G5x,2) - 
   24*G4pxx*G5pxxx*pow(G5x,2) + 12*G5ppx*G5pxxx*pow(G5x,2) + 96*pow(G5px,2)*G5pxx*G5xx - 12*G4x*pow(G5pxx,2)*G5xx + 6*G5p*pow(G5pxx,2)*G5xx + 72*G4x*G5px*G5pxxx*G5xx - 36*G5p*G5px*G5pxxx*G5xx + 
   216*G4pxxx*G5px*G5x*G5xx - 84*G5ppxx*G5px*G5x*G5xx + 18*G3xx*G5pxx*G5x*G5xx + 12*G4pxx*G5pxx*G5x*G5xx - 6*G5ppx*G5pxx*G5x*G5xx - 24*G3x*G5pxxx*G5x*G5xx + 72*G4px*G5pxxx*G5x*G5xx - 
   12*G5pp*G5pxxx*G5x*G5xx - 72*G4pxxx*G4x*pow(G5xx,2) + 36*G4pxxx*G5p*pow(G5xx,2) + 24*G4x*G5ppxx*pow(G5xx,2) - 12*G5p*G5ppxx*pow(G5xx,2) - 18*G3xx*G5px*pow(G5xx,2) - 
   288*G4pxx*G5px*pow(G5xx,2) + 108*G5ppx*G5px*pow(G5xx,2) + 45*G3x*G5pxx*pow(G5xx,2) - 102*G4px*G5pxx*pow(G5xx,2) + 6*G5pp*G5pxx*pow(G5xx,2) - 18*G3pxx*G5x*pow(G5xx,2) + 
   18*G5pppx*G5x*pow(G5xx,2) + 9*G3px*pow(G5xx,3) + 18*G4ppx*pow(G5xx,3) - 18*G5ppp*pow(G5xx,3) + 
   12*G4xxx*(-12*G5px*G5pxx*G5x + 2*G5ppxx*pow(G5x,2) + 18*pow(G5px,2)*G5xx + 3*G5xx*(2*G4x*G5pxx - G5p*G5pxx + 4*G4pxx*G5x - 3*G5ppx*G5x - 2*G4px*G5xx + G5pp*G5xx)) - 48*pow(G5px,3)*G5xxx - 
   48*G4x*G5px*G5pxx*G5xxx + 24*G5p*G5px*G5pxx*G5xxx + 16*G4x*G5ppxx*G5x*G5xxx - 8*G5p*G5ppxx*G5x*G5xxx - 96*G4pxx*G5px*G5x*G5xxx + 60*G5ppx*G5px*G5x*G5xxx + 18*G3x*G5pxx*G5x*G5xxx - 
   60*G4px*G5pxx*G5x*G5xxx + 12*G5pp*G5pxx*G5x*G5xxx + 24*G4ppxx*pow(G5x,2)*G5xxx - 12*G5pppx*pow(G5x,2)*G5xxx + 48*G4pxx*G4x*G5xx*G5xxx - 24*G4pxx*G5p*G5xx*G5xxx - 36*G4x*G5ppx*G5xx*G5xxx + 
   18*G5p*G5ppx*G5xx*G5xxx - 24*G3x*G5px*G5xx*G5xxx + 144*G4px*G5px*G5xx*G5xxx - 48*G5pp*G5px*G5xx*G5xxx + 6*G3px*G5x*G5xx*G5xxx - 36*G4ppx*G5x*G5xx*G5xxx + 12*G5ppp*G5x*G5xx*G5xxx - 
   144*pow(G4xx,2)*(G5pxxx*G5x - 4*G5pxx*G5xx + G5px*G5xxx) + 12*G4xx*(-5*pow(G5pxx,2)*G5x + 14*G5px*G5pxxx*G5x - 10*G4x*G5pxxx*G5xx + 5*G5p*G5pxxx*G5xx - 30*G4pxxx*G5x*G5xx + 11*G5ppxx*G5x*G5xx + 
      42*G4pxx*pow(G5xx,2) - 12*G5ppx*pow(G5xx,2) + 6*G4xxx*(3*G5pxx*G5x - 5*G5px*G5xx) + 14*pow(G5px,2)*G5xxx + 12*G4pxx*G5x*G5xxx - 8*G5ppx*G5x*G5xxx - 10*G4px*G5xx*G5xxx + 5*G5pp*G5xx*G5xxx + 
      G5pxx*(-43*G5px*G5xx + 6*G4x*G5xxx - 3*G5p*G5xxx)) - 48*G5pxxx*pow(G5x,2)*G5xx*pow(H,2) + 99*G5pxx*G5x*pow(G5xx,2)*pow(H,2) - 15*G5px*pow(G5xx,3)*pow(H,2) + 
   30*G5pxx*pow(G5x,2)*G5xxx*pow(H,2) - 30*G5px*G5x*G5xx*G5xxx*pow(H,2)) + 
pow(a,3)*pow(dphi,14)*pow(H,6)*(864*pow(G4xx,3)*G5pxx - 24*pow(G5px,3)*G5pxx + 96*G4x*G5px*pow(G5pxx,2) - 48*G5p*G5px*pow(G5pxx,2) - 96*G4x*pow(G5px,2)*G5pxxx + 
   48*G5p*pow(G5px,2)*G5pxxx - 16*pow(G4x,2)*G5pxx*G5pxxx + 16*G4x*G5p*G5pxx*G5pxxx - 4*pow(G5p,2)*G5pxx*G5pxxx - 288*G4pxxx*pow(G5px,2)*G5x + 120*G5ppxx*pow(G5px,2)*G5x - 
   96*G4pxxx*G4x*G5pxx*G5x + 48*G4pxxx*G5p*G5pxx*G5x - 72*G3xx*G5px*G5pxx*G5x + 192*G4pxx*G5px*G5pxx*G5x - 60*G5ppx*G5px*G5pxx*G5x - 42*G3x*pow(G5pxx,2)*G5x + 156*G4px*pow(G5pxx,2)*G5x - 
   36*G5pp*pow(G5pxx,2)*G5x - 96*G4pxx*G4x*G5pxxx*G5x + 48*G4pxx*G5p*G5pxxx*G5x + 48*G4x*G5ppx*G5pxxx*G5x - 24*G5p*G5ppx*G5pxxx*G5x + 60*G3x*G5px*G5pxxx*G5x - 168*G4px*G5px*G5pxxx*G5x + 
   24*G5pp*G5px*G5pxxx*G5x - 144*G4pxx*G4pxxx*pow(G5x,2) + 72*G4pxxx*G5ppx*pow(G5x,2) + 12*G3xx*G5ppxx*pow(G5x,2) + 48*G4pxx*G5ppxx*pow(G5x,2) - 36*G5ppx*G5ppxx*pow(G5x,2) - 
   12*G3pxx*G5pxx*pow(G5x,2) - 48*G4ppxx*G5pxx*pow(G5x,2) + 36*G5pppx*G5pxx*pow(G5x,2) - 12*G3px*G5pxxx*pow(G5x,2) + 24*G4ppx*G5pxxx*pow(G5x,2) + 432*G4pxxx*G4x*G5px*G5xx - 
   216*G4pxxx*G5p*G5px*G5xx - 168*G4x*G5ppxx*G5px*G5xx + 84*G5p*G5ppxx*G5px*G5xx + 108*G3xx*pow(G5px,2)*G5xx + 648*G4pxx*pow(G5px,2)*G5xx - 324*G5ppx*pow(G5px,2)*G5xx + 36*G3xx*G4x*G5pxx*G5xx + 
   24*G4pxx*G4x*G5pxx*G5xx - 18*G3xx*G5p*G5pxx*G5xx - 12*G4pxx*G5p*G5pxx*G5xx - 12*G4x*G5ppx*G5pxx*G5xx + 6*G5p*G5ppx*G5pxx*G5xx - 156*G3x*G5px*G5pxx*G5xx + 216*G4px*G5px*G5pxx*G5xx + 
   48*G5pp*G5px*G5pxx*G5xx - 48*G3x*G4x*G5pxxx*G5xx + 144*G4px*G4x*G5pxxx*G5xx + 24*G3x*G5p*G5pxxx*G5xx - 72*G4px*G5p*G5pxxx*G5xx - 24*G4x*G5pp*G5pxxx*G5xx + 12*G5p*G5pp*G5pxxx*G5xx + 
   72*G3xx*G4pxx*G5x*G5xx + 144*pow(G4pxx,2)*G5x*G5xx - 144*G3x*G4pxxx*G5x*G5xx + 432*G4px*G4pxxx*G5x*G5xx - 72*G4pxxx*G5pp*G5x*G5xx - 54*G3xx*G5ppx*G5x*G5xx - 180*G4pxx*G5ppx*G5x*G5xx + 
   72*pow(G5ppx,2)*G5x*G5xx + 60*G3x*G5ppxx*G5x*G5xx - 192*G4px*G5ppxx*G5x*G5xx + 36*G5pp*G5ppxx*G5x*G5xx + 108*G3pxx*G5px*G5x*G5xx - 72*G4ppxx*G5px*G5x*G5xx - 72*G5pppx*G5px*G5x*G5xx + 
   6*G2xx*G5pxx*G5x*G5xx + 24*G3px*G5pxx*G5x*G5xx + 12*G4ppx*G5pxx*G5x*G5xx - 36*G5ppp*G5pxx*G5x*G5xx - 6*G2x*G5pxxx*G5x*G5xx + 12*G3p*G5pxxx*G5x*G5xx - 12*G4pp*G5pxxx*G5x*G5xx - 
   36*G3xx*G4px*pow(G5xx,2) + 180*G3x*G4pxx*pow(G5xx,2) - 504*G4px*G4pxx*pow(G5xx,2) - 36*G3pxx*G4x*pow(G5xx,2) + 18*G3pxx*G5p*pow(G5xx,2) + 18*G3xx*G5pp*pow(G5xx,2) + 
   72*G4pxx*G5pp*pow(G5xx,2) + 36*G4x*G5pppx*pow(G5xx,2) - 18*G5p*G5pppx*pow(G5xx,2) - 63*G3x*G5ppx*pow(G5xx,2) + 234*G4px*G5ppx*pow(G5xx,2) - 54*G5pp*G5ppx*pow(G5xx,2) - 
   6*G2xx*G5px*pow(G5xx,2) - 84*G3px*G5px*pow(G5xx,2) - 36*G4ppx*G5px*pow(G5xx,2) + 108*G5ppp*G5px*pow(G5xx,2) + 9*G2x*G5pxx*pow(G5xx,2) - 18*G3p*G5pxx*pow(G5xx,2) + 
   18*G4pp*G5pxx*pow(G5xx,2) - 6*G2pxx*G5x*pow(G5xx,2) - 12*G3ppx*G5x*pow(G5xx,2) + 36*G4pppx*G5x*pow(G5xx,2) + 9*G3pp*pow(G5xx,3) - 18*G4ppp*pow(G5xx,3) + 16*pow(G4x,2)*G5ppxx*G5xxx - 
   16*G4x*G5p*G5ppxx*G5xxx + 4*pow(G5p,2)*G5ppxx*G5xxx - 192*G4pxx*G4x*G5px*G5xxx + 96*G4pxx*G5p*G5px*G5xxx + 120*G4x*G5ppx*G5px*G5xxx - 60*G5p*G5ppx*G5px*G5xxx + 60*G3x*pow(G5px,2)*G5xxx - 
   264*G4px*pow(G5px,2)*G5xxx + 72*G5pp*pow(G5px,2)*G5xxx + 36*G3x*G4x*G5pxx*G5xxx - 120*G4px*G4x*G5pxx*G5xxx - 18*G3x*G5p*G5pxx*G5xxx + 60*G4px*G5p*G5pxx*G5xxx + 24*G4x*G5pp*G5pxx*G5xxx - 
   12*G5p*G5pp*G5pxx*G5xxx + 72*G3x*G4pxx*G5x*G5xxx - 240*G4px*G4pxx*G5x*G5xxx + 96*G4ppxx*G4x*G5x*G5xxx - 48*G4ppxx*G5p*G5x*G5xxx + 48*G4pxx*G5pp*G5x*G5xxx - 48*G4x*G5pppx*G5x*G5xxx + 
   24*G5p*G5pppx*G5x*G5xxx - 42*G3x*G5ppx*G5x*G5xxx + 132*G4px*G5ppx*G5x*G5xxx - 24*G5pp*G5ppx*G5x*G5xxx - 24*G3px*G5px*G5x*G5xxx + 96*G4ppx*G5px*G5x*G5xxx - 24*G5ppp*G5px*G5x*G5xxx + 
   6*G2x*G5pxx*G5x*G5xxx - 12*G3p*G5pxx*G5x*G5xxx + 12*G4pp*G5pxx*G5x*G5xxx + 12*G3ppx*pow(G5x,2)*G5xxx - 24*G4pppx*pow(G5x,2)*G5xxx - 48*G3x*G4px*G5xx*G5xxx + 144*pow(G4px,2)*G5xx*G5xxx + 
   12*G3px*G4x*G5xx*G5xxx - 72*G4ppx*G4x*G5xx*G5xxx - 6*G3px*G5p*G5xx*G5xxx + 36*G4ppx*G5p*G5xx*G5xxx + 24*G3x*G5pp*G5xx*G5xxx - 96*G4px*G5pp*G5xx*G5xxx + 12*pow(G5pp,2)*G5xx*G5xxx + 
   24*G4x*G5ppp*G5xx*G5xxx - 12*G5p*G5ppp*G5xx*G5xxx - 6*G2x*G5px*G5xx*G5xxx + 12*G3p*G5px*G5xx*G5xxx - 12*G4pp*G5px*G5xx*G5xxx - 6*G3pp*G5x*G5xx*G5xxx + 12*G4ppp*G5x*G5xx*G5xxx - 
   48*pow(G4xx,2)*(18*G4xxx*G5px + 23*G5px*G5pxx + 6*G4x*G5pxxx - 3*G5p*G5pxxx + 18*G4pxxx*G5x - 7*G5ppxx*G5x - 48*G4pxx*G5xx + 15*G5ppx*G5xx + 6*G4px*G5xxx - 3*G5pp*G5xxx) - 
   82*pow(G5pxx,2)*pow(G5x,2)*pow(H,2) + 136*G5px*G5pxxx*pow(G5x,2)*pow(H,2) - 432*G5px*G5pxx*G5x*G5xx*pow(H,2) - 180*G4x*G5pxxx*G5x*G5xx*pow(H,2) + 120*G5p*G5pxxx*G5x*G5xx*pow(H,2) - 
   288*G4pxxx*pow(G5x,2)*G5xx*pow(H,2) + 136*G5ppxx*pow(G5x,2)*G5xx*pow(H,2) + 132*pow(G5px,2)*pow(G5xx,2)*pow(H,2) + 162*G4x*G5pxx*pow(G5xx,2)*pow(H,2) - 
   144*G5p*G5pxx*pow(G5xx,2)*pow(H,2) + 12*G4*G5pxxx*pow(G5xx,2)*pow(H,2) + 384*G4pxx*G5x*pow(G5xx,2)*pow(H,2) - 147*G5ppx*G5x*pow(G5xx,2)*pow(H,2) - 30*G4px*pow(G5xx,3)*pow(H,2) - 
   3*G5pp*pow(G5xx,3)*pow(H,2) + 60*pow(G5px,2)*G5x*G5xxx*pow(H,2) + 120*G4x*G5pxx*G5x*G5xxx*pow(H,2) - 78*G5p*G5pxx*G5x*G5xxx*pow(H,2) + 144*G4pxx*pow(G5x,2)*G5xxx*pow(H,2) - 
   82*G5ppx*pow(G5x,2)*G5xxx*pow(H,2) - 48*G4x*G5px*G5xx*G5xxx*pow(H,2) + 54*G5p*G5px*G5xx*G5xxx*pow(H,2) - 12*G4*G5pxx*G5xx*G5xxx*pow(H,2) - 84*G4px*G5x*G5xx*G5xxx*pow(H,2) + 
   30*G5pp*G5x*G5xx*G5xxx*pow(H,2) + 12*G4p*pow(G5xx,2)*G5xxx*pow(H,2) - 
   12*G4xxx*(24*pow(G5px,3) + 4*G5p*G5ppxx*G5x - 9*G3x*G5pxx*G5x + 30*G4px*G5pxx*G5x - 6*G5pp*G5pxx*G5x - 12*G4ppxx*pow(G5x,2) + 6*G5pppx*pow(G5x,2) + 12*G4pxx*G5p*G5xx - 9*G5p*G5ppx*G5xx - 
      3*G3px*G5x*G5xx + 18*G4ppx*G5x*G5xx - 6*G5ppp*G5x*G5xx - 2*G4x*(4*G5ppxx*G5x + 12*G4pxx*G5xx - 9*G5ppx*G5xx) - 15*G5pxx*pow(G5x,2)*pow(H,2) + 
      3*G5px*(8*G4x*G5pxx - 4*G5p*G5pxx + 16*G4pxx*G5x - 10*G5ppx*G5x + 4*G3x*G5xx - 24*G4px*G5xx + 8*G5pp*G5xx + 5*G5x*G5xx*pow(H,2))) + 
   6*G4xx*(64*pow(G5px,2)*G5pxx - 20*G4x*pow(G5pxx,2) + 10*G5p*pow(G5pxx,2) + 18*G3xx*G5pxx*G5x - 28*G4pxx*G5pxx*G5x + 8*G5ppx*G5pxx*G5x - 18*G3x*G5pxxx*G5x + 52*G4px*G5pxxx*G5x - 
      8*G5pp*G5pxxx*G5x - 120*G4pxxx*G4x*G5xx + 60*G4pxxx*G5p*G5xx + 44*G4x*G5ppxx*G5xx - 22*G5p*G5ppxx*G5xx + 66*G3x*G5pxx*G5xx - 136*G4px*G5pxx*G5xx + 2*G5pp*G5pxx*G5xx - 30*G3pxx*G5x*G5xx + 
      12*G4ppxx*G5x*G5xx + 24*G5pppx*G5x*G5xx + 21*G3px*pow(G5xx,2) + 30*G4ppx*pow(G5xx,2) - 36*G5ppp*pow(G5xx,2) + 
      12*G4xxx*(14*pow(G5px,2) + 6*G4x*G5pxx - 3*G5p*G5pxx + 12*G4pxx*G5x - 8*G5ppx*G5x - 10*G4px*G5xx + 5*G5pp*G5xx) + 48*G4pxx*G4x*G5xxx - 24*G4pxx*G5p*G5xxx - 32*G4x*G5ppx*G5xxx + 
      16*G5p*G5ppx*G5xxx + 6*G3px*G5x*G5xxx - 28*G4ppx*G5x*G5xxx + 8*G5ppp*G5x*G5xxx - 42*G5pxxx*pow(G5x,2)*pow(H,2) + 158*G5pxx*G5x*G5xx*pow(H,2) + 
      G5px*(56*G4x*G5pxxx - 28*G5p*G5pxxx + 168*G4pxxx*G5x - 68*G5ppxx*G5x - 30*G3xx*G5xx - 420*G4pxx*G5xx + 168*G5ppx*G5xx - 18*G3x*G5xxx + 108*G4px*G5xxx - 36*G5pp*G5xxx - 37*pow(G5xx,2)*pow(H,2) - 
         20*G5x*G5xxx*pow(H,2)))) + pow(a,4)*pow(dphi,13)*pow(H,5)*(-1728*G4px*pow(G4xx,2)*G4xxx + 864*pow(G4xx,2)*G4xxx*G5pp - 1152*pow(G4xx,3)*G5ppx - 1152*G4x*G4xx*G4xxx*G5ppx + 
   576*G4xx*G4xxx*G5p*G5ppx + 672*G4x*pow(G4xx,2)*G5ppxx + 96*pow(G4x,2)*G4xxx*G5ppxx - 336*pow(G4xx,2)*G5p*G5ppxx - 96*G4x*G4xxx*G5p*G5ppxx + 24*G4xxx*pow(G5p,2)*G5ppxx - 
   432*G3xx*pow(G4xx,2)*G5px - 648*G3x*G4xx*G4xxx*G5px + 3888*G4px*G4xx*G4xxx*G5px - 1296*G4xx*G4xxx*G5pp*G5px + 2304*pow(G4xx,2)*G5ppx*G5px + 720*G4x*G4xxx*G5ppx*G5px - 360*G4xxx*G5p*G5ppx*G5px - 
   816*G4x*G4xx*G5ppxx*G5px + 408*G4xx*G5p*G5ppxx*G5px + 504*G3xx*G4xx*pow(G5px,2) + 360*G3x*G4xxx*pow(G5px,2) - 1584*G4px*G4xxx*pow(G5px,2) + 432*G4xxx*G5pp*pow(G5px,2) - 
   1440*G4xx*G5ppx*pow(G5px,2) + 240*G4x*G5ppxx*pow(G5px,2) - 120*G5p*G5ppxx*pow(G5px,2) - 144*G3xx*pow(G5px,3) + 288*G5ppx*pow(G5px,3) + 216*G3xx*G4x*G4xx*G5pxx + 864*G3x*pow(G4xx,2)*G5pxx - 
   1632*G4px*pow(G4xx,2)*G5pxx + 216*G3x*G4x*G4xxx*G5pxx - 720*G4px*G4x*G4xxx*G5pxx - 108*G3xx*G4xx*G5p*G5pxx - 108*G3x*G4xxx*G5p*G5pxx + 360*G4px*G4xxx*G5p*G5pxx - 48*pow(G4xx,2)*G5pp*G5pxx + 
   144*G4x*G4xxx*G5pp*G5pxx - 72*G4xxx*G5p*G5pp*G5pxx + 96*G4x*G4xx*G5ppx*G5pxx - 48*G4xx*G5p*G5ppx*G5pxx - 144*G3xx*G4x*G5px*G5pxx - 660*G3x*G4xx*G5px*G5pxx + 792*G4px*G4xx*G5px*G5pxx + 
   72*G3xx*G5p*G5px*G5pxx + 264*G4xx*G5pp*G5px*G5pxx - 120*G4x*G5ppx*G5px*G5pxx + 60*G5p*G5ppx*G5px*G5pxx + 96*G3x*pow(G5px,2)*G5pxx + 48*G4px*pow(G5px,2)*G5pxx - 120*G5pp*pow(G5px,2)*G5pxx - 
   84*G3x*G4x*pow(G5pxx,2) + 312*G4px*G4x*pow(G5pxx,2) + 42*G3x*G5p*pow(G5pxx,2) - 156*G4px*G5p*pow(G5pxx,2) - 72*G4x*G5pp*pow(G5pxx,2) + 36*G5p*G5pp*pow(G5pxx,2) - 216*G3x*G4x*G4xx*G5pxxx + 
   624*G4px*G4x*G4xx*G5pxxx + 108*G3x*G4xx*G5p*G5pxxx - 312*G4px*G4xx*G5p*G5pxxx - 96*G4x*G4xx*G5pp*G5pxxx + 48*G4xx*G5p*G5pp*G5pxxx + 48*pow(G4x,2)*G5ppx*G5pxxx - 48*G4x*G5p*G5ppx*G5pxxx + 
   12*pow(G5p,2)*G5ppx*G5pxxx + 120*G3x*G4x*G5px*G5pxxx - 336*G4px*G4x*G5px*G5pxxx - 60*G3x*G5p*G5px*G5pxxx + 168*G4px*G5p*G5px*G5pxxx + 48*G4x*G5pp*G5px*G5pxxx - 24*G5p*G5pp*G5px*G5pxxx - 
   432*G3pxx*pow(G4xx,2)*G5x + 288*G4ppxx*pow(G4xx,2)*G5x + 576*G4ppxx*G4x*G4xxx*G5x + 216*G3px*G4xx*G4xxx*G5x - 1008*G4ppx*G4xx*G4xxx*G5x - 288*G4ppxx*G4xxx*G5p*G5x + 288*G4xx*G4xxx*G5ppp*G5x + 
   288*pow(G4xx,2)*G5pppx*G5x - 288*G4x*G4xxx*G5pppx*G5x + 144*G4xxx*G5p*G5pppx*G5x - 288*G3xx*G4xx*G5ppx*G5x - 252*G3x*G4xxx*G5ppx*G5x + 792*G4px*G4xxx*G5ppx*G5x - 144*G4xxx*G5pp*G5ppx*G5x + 
   288*G4xx*pow(G5ppx,2)*G5x + 48*G3xx*G4x*G5ppxx*G5x + 276*G3x*G4xx*G5ppxx*G5x - 840*G4px*G4xx*G5ppxx*G5x - 24*G3xx*G5p*G5ppxx*G5x + 144*G4xx*G5pp*G5ppxx*G5x - 144*G4x*G5ppx*G5ppxx*G5x + 
   72*G5p*G5ppx*G5ppxx*G5x + 504*G3pxx*G4xx*G5px*G5x - 432*G4ppxx*G4xx*G5px*G5x - 144*G3px*G4xxx*G5px*G5x + 576*G4ppx*G4xxx*G5px*G5x - 144*G4xxx*G5ppp*G5px*G5x - 288*G4xx*G5pppx*G5px*G5x + 
   180*G3xx*G5ppx*G5px*G5x - 144*pow(G5ppx,2)*G5px*G5x - 156*G3x*G5ppxx*G5px*G5x + 456*G4px*G5ppxx*G5px*G5x - 72*G5pp*G5ppxx*G5px*G5x - 144*G3pxx*pow(G5px,2)*G5x + 144*G4ppxx*pow(G5px,2)*G5x + 
   72*G5pppx*pow(G5px,2)*G5x + 54*G3x*G3xx*G5pxx*G5x - 180*G3xx*G4px*G5pxx*G5x - 48*G3pxx*G4x*G5pxx*G5x - 192*G4ppxx*G4x*G5pxx*G5x + 36*G2xx*G4xx*G5pxx*G5x + 48*G3px*G4xx*G5pxx*G5x + 
   120*G4ppx*G4xx*G5pxx*G5x + 36*G2x*G4xxx*G5pxx*G5x - 72*G3p*G4xxx*G5pxx*G5x + 72*G4pp*G4xxx*G5pxx*G5x + 24*G3pxx*G5p*G5pxx*G5x + 96*G4ppxx*G5p*G5pxx*G5x + 36*G3xx*G5pp*G5pxx*G5x - 
   144*G4xx*G5ppp*G5pxx*G5x + 144*G4x*G5pppx*G5pxx*G5x - 72*G5p*G5pppx*G5pxx*G5x + 66*G3x*G5ppx*G5pxx*G5x - 276*G4px*G5ppx*G5pxx*G5x + 72*G5pp*G5ppx*G5pxx*G5x - 24*G2xx*G5px*G5pxx*G5x - 
   96*G4ppx*G5px*G5pxx*G5x + 72*G5ppp*G5px*G5pxx*G5x - 18*G2x*pow(G5pxx,2)*G5x + 36*G3p*pow(G5pxx,2)*G5x - 36*G4pp*pow(G5pxx,2)*G5x - 18*pow(G3x,2)*G5pxxx*G5x + 96*G3x*G4px*G5pxxx*G5x - 
   120*pow(G4px,2)*G5pxxx*G5x - 48*G3px*G4x*G5pxxx*G5x + 96*G4ppx*G4x*G5pxxx*G5x - 24*G2x*G4xx*G5pxxx*G5x + 48*G3p*G4xx*G5pxxx*G5x - 48*G4pp*G4xx*G5pxxx*G5x + 24*G3px*G5p*G5pxxx*G5x - 
   48*G4ppx*G5p*G5pxxx*G5x - 12*G3x*G5pp*G5pxxx*G5x + 24*G4px*G5pp*G5pxxx*G5x + 12*G2x*G5px*G5pxxx*G5x - 24*G3p*G5px*G5pxxx*G5x + 24*G4pp*G5px*G5pxxx*G5x + 72*G3xx*G4ppxx*pow(G5x,2) + 
   72*G3ppx*G4xxx*pow(G5x,2) - 144*G4pppx*G4xxx*pow(G5x,2) - 36*G3xx*G5pppx*pow(G5x,2) + 36*G3pxx*G5ppx*pow(G5x,2) - 72*G4ppxx*G5ppx*pow(G5x,2) + 4*G2xx*G5ppxx*pow(G5x,2) + 
   32*G3px*G5ppxx*pow(G5x,2) - 72*G4ppx*G5ppxx*pow(G5x,2) - 4*G2pxx*G5pxx*pow(G5x,2) - 32*G3ppx*G5pxx*pow(G5x,2) + 72*G4pppx*G5pxx*pow(G5x,2) - 4*G2px*G5pxxx*pow(G5x,2) + 
   4*G3pp*G5pxxx*pow(G5x,2) - 360*G3xx*G4px*G4xx*G5xx - 360*G3pxx*G4x*G4xx*G5xx + 144*G4ppxx*G4x*G4xx*G5xx + 576*G3px*pow(G4xx,2)*G5xx + 576*G4ppx*pow(G4xx,2)*G5xx - 288*G3x*G4px*G4xxx*G5xx + 
   864*pow(G4px,2)*G4xxx*G5xx + 72*G3px*G4x*G4xxx*G5xx - 432*G4ppx*G4x*G4xxx*G5xx + 180*G3pxx*G4xx*G5p*G5xx - 72*G4ppxx*G4xx*G5p*G5xx - 36*G3px*G4xxx*G5p*G5xx + 216*G4ppx*G4xxx*G5p*G5xx + 
   180*G3xx*G4xx*G5pp*G5xx + 144*G3x*G4xxx*G5pp*G5xx - 576*G4px*G4xxx*G5pp*G5xx + 72*G4xxx*pow(G5pp,2)*G5xx - 864*pow(G4xx,2)*G5ppp*G5xx + 144*G4x*G4xxx*G5ppp*G5xx - 72*G4xxx*G5p*G5ppp*G5xx + 
   288*G4x*G4xx*G5pppx*G5xx - 144*G4xx*G5p*G5pppx*G5xx - 108*G3xx*G4x*G5ppx*G5xx - 576*G3x*G4xx*G5ppx*G5xx + 2016*G4px*G4xx*G5ppx*G5xx + 54*G3xx*G5p*G5ppx*G5xx - 432*G4xx*G5pp*G5ppx*G5xx + 
   144*G4x*pow(G5ppx,2)*G5xx - 72*G5p*pow(G5ppx,2)*G5xx + 120*G3x*G4x*G5ppxx*G5xx - 384*G4px*G4x*G5ppxx*G5xx - 60*G3x*G5p*G5ppxx*G5xx + 192*G4px*G5p*G5ppxx*G5xx + 72*G4x*G5pp*G5ppxx*G5xx - 
   36*G5p*G5pp*G5ppxx*G5xx - 72*G3x*G3xx*G5px*G5xx + 432*G3xx*G4px*G5px*G5xx + 216*G3pxx*G4x*G5px*G5xx - 144*G4ppxx*G4x*G5px*G5xx - 60*G2xx*G4xx*G5px*G5xx - 732*G3px*G4xx*G5px*G5xx - 
   144*G4ppx*G4xx*G5px*G5xx - 36*G2x*G4xxx*G5px*G5xx + 72*G3p*G4xxx*G5px*G5xx - 72*G4pp*G4xxx*G5px*G5xx - 108*G3pxx*G5p*G5px*G5xx + 72*G4ppxx*G5p*G5px*G5xx - 144*G3xx*G5pp*G5px*G5xx + 
   864*G4xx*G5ppp*G5px*G5xx - 144*G4x*G5pppx*G5px*G5xx + 72*G5p*G5pppx*G5px*G5xx + 360*G3x*G5ppx*G5px*G5xx - 1152*G4px*G5ppx*G5px*G5xx + 216*G5pp*G5ppx*G5px*G5xx + 36*G2xx*pow(G5px,2)*G5xx + 
   216*G3px*pow(G5px,2)*G5xx - 72*G4ppx*pow(G5px,2)*G5xx - 216*G5ppp*pow(G5px,2)*G5xx + 63*pow(G3x,2)*G5pxx*G5xx - 204*G3x*G4px*G5pxx*G5xx + 12*pow(G4px,2)*G5pxx*G5xx + 12*G2xx*G4x*G5pxx*G5xx + 
   48*G3px*G4x*G5pxx*G5xx + 24*G4ppx*G4x*G5pxx*G5xx + 72*G2x*G4xx*G5pxx*G5xx - 144*G3p*G4xx*G5pxx*G5xx + 144*G4pp*G4xx*G5pxx*G5xx - 6*G2xx*G5p*G5pxx*G5xx - 24*G3px*G5p*G5pxx*G5xx - 
   12*G4ppx*G5p*G5pxx*G5xx - 24*G3x*G5pp*G5pxx*G5xx + 192*G4px*G5pp*G5pxx*G5xx - 36*pow(G5pp,2)*G5pxx*G5xx - 72*G4x*G5ppp*G5pxx*G5xx + 36*G5p*G5ppp*G5pxx*G5xx - 18*G2x*G5px*G5pxx*G5xx + 
   36*G3p*G5px*G5pxx*G5xx - 36*G4pp*G5px*G5pxx*G5xx - 12*G2x*G4x*G5pxxx*G5xx + 24*G3p*G4x*G5pxxx*G5xx - 24*G4pp*G4x*G5pxxx*G5xx + 6*G2x*G5p*G5pxxx*G5xx - 12*G3p*G5p*G5pxxx*G5xx + 
   12*G4pp*G5p*G5pxxx*G5xx - 72*G3pxx*G3x*G5x*G5xx + 18*G3px*G3xx*G5x*G5xx - 108*G3xx*G4ppx*G5x*G5xx + 72*G3x*G4ppxx*G5x*G5xx + 216*G3pxx*G4px*G5x*G5xx - 288*G4ppxx*G4px*G5x*G5xx - 
   60*G2pxx*G4xx*G5x*G5xx - 84*G3ppx*G4xx*G5x*G5xx + 288*G4pppx*G4xx*G5x*G5xx - 36*G3pp*G4xxx*G5x*G5xx + 72*G4ppp*G4xxx*G5x*G5xx - 36*G3pxx*G5pp*G5x*G5xx + 72*G4ppxx*G5pp*G5x*G5xx + 
   36*G3xx*G5ppp*G5x*G5xx + 36*G3x*G5pppx*G5x*G5xx - 72*G4px*G5pppx*G5x*G5xx - 18*G2xx*G5ppx*G5x*G5xx - 90*G3px*G5ppx*G5x*G5xx + 216*G4ppx*G5ppx*G5x*G5xx + 18*G2x*G5ppxx*G5x*G5xx - 
   36*G3p*G5ppxx*G5x*G5xx + 36*G4pp*G5ppxx*G5x*G5xx + 36*G2pxx*G5px*G5x*G5xx + 36*G3ppx*G5px*G5x*G5xx - 144*G4pppx*G5px*G5x*G5xx + 12*G2px*G5pxx*G5x*G5xx + 6*G3pp*G5pxx*G5x*G5xx - 
   36*G4ppp*G5pxx*G5x*G5xx + 45*G3px*G3x*pow(G5xx,2) + 18*G3x*G4ppx*pow(G5xx,2) - 12*G2xx*G4px*pow(G5xx,2) - 150*G3px*G4px*pow(G5xx,2) + 108*G4ppx*G4px*pow(G5xx,2) - 
   12*G2pxx*G4x*pow(G5xx,2) - 24*G3ppx*G4x*pow(G5xx,2) + 72*G4pppx*G4x*pow(G5xx,2) + 108*G3pp*G4xx*pow(G5xx,2) - 216*G4ppp*G4xx*pow(G5xx,2) + 6*G2pxx*G5p*pow(G5xx,2) + 
   12*G3ppx*G5p*pow(G5xx,2) - 36*G4pppx*G5p*pow(G5xx,2) + 6*G2xx*G5pp*pow(G5xx,2) + 30*G3px*G5pp*pow(G5xx,2) - 72*G4ppx*G5pp*pow(G5xx,2) - 54*G3x*G5ppp*pow(G5xx,2) + 
   108*G4px*G5ppp*pow(G5xx,2) - 18*G2x*G5ppx*pow(G5xx,2) + 36*G3p*G5ppx*pow(G5xx,2) - 36*G4pp*G5ppx*pow(G5xx,2) - 6*G2px*G5px*pow(G5xx,2) - 48*G3pp*G5px*pow(G5xx,2) + 
   108*G4ppp*G5px*pow(G5xx,2) - 6*G2ppx*G5x*pow(G5xx,2) + 6*G3ppp*G5x*pow(G5xx,2) + 144*pow(G4pxx,2)*(2*G4xx*G5x + 2*G4x*G5xx - G5p*G5xx) + 96*G4ppxx*pow(G4x,2)*G5xxx - 
   216*G3x*G4px*G4xx*G5xxx + 624*pow(G4px,2)*G4xx*G5xxx + 72*G3px*G4x*G4xx*G5xxx - 336*G4ppx*G4x*G4xx*G5xxx - 96*G4ppxx*G4x*G5p*G5xxx - 36*G3px*G4xx*G5p*G5xxx + 168*G4ppx*G4xx*G5p*G5xxx + 
   24*G4ppxx*pow(G5p,2)*G5xxx + 108*G3x*G4xx*G5pp*G5xxx - 408*G4px*G4xx*G5pp*G5xxx + 48*G4xx*pow(G5pp,2)*G5xxx + 96*G4x*G4xx*G5ppp*G5xxx - 48*G4xx*G5p*G5ppp*G5xxx - 48*pow(G4x,2)*G5pppx*G5xxx + 
   48*G4x*G5p*G5pppx*G5xxx - 12*pow(G5p,2)*G5pppx*G5xxx - 84*G3x*G4x*G5ppx*G5xxx + 264*G4px*G4x*G5ppx*G5xxx + 42*G3x*G5p*G5ppx*G5xxx - 132*G4px*G5p*G5ppx*G5xxx - 48*G4x*G5pp*G5ppx*G5xxx + 
   24*G5p*G5pp*G5ppx*G5xxx - 18*pow(G3x,2)*G5px*G5xxx + 216*G3x*G4px*G5px*G5xxx - 456*pow(G4px,2)*G5px*G5xxx - 48*G3px*G4x*G5px*G5xxx + 192*G4ppx*G4x*G5px*G5xxx - 24*G2x*G4xx*G5px*G5xxx + 
   48*G3p*G4xx*G5px*G5xxx - 48*G4pp*G4xx*G5px*G5xxx + 24*G3px*G5p*G5px*G5xxx - 96*G4ppx*G5p*G5px*G5xxx - 72*G3x*G5pp*G5px*G5xxx + 240*G4px*G5pp*G5px*G5xxx - 24*pow(G5pp,2)*G5px*G5xxx - 
   48*G4x*G5ppp*G5px*G5xxx + 24*G5p*G5ppp*G5px*G5xxx + 12*G2x*pow(G5px,2)*G5xxx - 24*G3p*pow(G5px,2)*G5xxx + 24*G4pp*pow(G5px,2)*G5xxx + 12*G2x*G4x*G5pxx*G5xxx - 24*G3p*G4x*G5pxx*G5xxx + 
   24*G4pp*G4x*G5pxx*G5xxx - 6*G2x*G5p*G5pxx*G5xxx + 12*G3p*G5p*G5pxx*G5xxx - 12*G4pp*G5p*G5pxx*G5xxx + 18*G3px*G3x*G5x*G5xxx - 60*G3x*G4ppx*G5x*G5xxx - 60*G3px*G4px*G5x*G5xxx + 
   168*G4ppx*G4px*G5x*G5xxx + 48*G3ppx*G4x*G5x*G5xxx - 96*G4pppx*G4x*G5x*G5xxx - 24*G3pp*G4xx*G5x*G5xxx + 48*G4ppp*G4xx*G5x*G5xxx - 24*G3ppx*G5p*G5x*G5xxx + 48*G4pppx*G5p*G5x*G5xxx + 
   12*G3px*G5pp*G5x*G5xxx - 24*G4ppx*G5pp*G5x*G5xxx + 12*G3x*G5ppp*G5x*G5xxx - 24*G4px*G5ppp*G5x*G5xxx - 12*G2x*G5ppx*G5x*G5xxx + 24*G3p*G5ppx*G5x*G5xxx - 24*G4pp*G5ppx*G5x*G5xxx + 
   12*G3pp*G5px*G5x*G5xxx - 24*G4ppp*G5px*G5x*G5xxx + 4*G2ppx*pow(G5x,2)*G5xxx - 4*G3ppp*pow(G5x,2)*G5xxx - 12*G2x*G4px*G5xx*G5xxx + 24*G3p*G4px*G5xx*G5xxx - 24*G4pp*G4px*G5xx*G5xxx - 
   12*G3pp*G4x*G5xx*G5xxx + 24*G4ppp*G4x*G5xx*G5xxx + 6*G3pp*G5p*G5xx*G5xxx - 12*G4ppp*G5p*G5xx*G5xxx + 6*G2x*G5pp*G5xx*G5xxx - 12*G3p*G5pp*G5xx*G5xxx + 12*G4pp*G5pp*G5xx*G5xxx - 
   720*G4xx*G4xxx*G5px*G5x*pow(H,2) + 360*G4xxx*pow(G5px,2)*G5x*pow(H,2) + 2160*pow(G4xx,2)*G5pxx*G5x*pow(H,2) + 720*G4x*G4xxx*G5pxx*G5x*pow(H,2) - 468*G4xxx*G5p*G5pxx*G5x*pow(H,2) - 
   2004*G4xx*G5px*G5pxx*G5x*pow(H,2) + 480*pow(G5px,2)*G5pxx*G5x*pow(H,2) - 352*G4x*pow(G5pxx,2)*G5x*pow(H,2) + 206*G5p*pow(G5pxx,2)*G5x*pow(H,2) - 936*G4x*G4xx*G5pxxx*G5x*pow(H,2) + 
   612*G4xx*G5p*G5pxxx*G5x*pow(H,2) + 496*G4x*G5px*G5pxxx*G5x*pow(H,2) - 332*G5p*G5px*G5pxxx*G5x*pow(H,2) + 16*G4*G5pxx*G5pxxx*G5x*pow(H,2) - 492*G4xxx*G5ppx*pow(G5x,2)*pow(H,2) + 
   708*G4xx*G5ppxx*pow(G5x,2)*pow(H,2) - 396*G5ppxx*G5px*pow(G5x,2)*pow(H,2) + 90*G3xx*G5pxx*pow(G5x,2)*pow(H,2) + 210*G5ppx*G5pxx*pow(G5x,2)*pow(H,2) - 
   108*G3x*G5pxxx*pow(G5x,2)*pow(H,2) + 288*G4px*G5pxxx*pow(G5x,2)*pow(H,2) - 24*G5pp*G5pxxx*pow(G5x,2)*pow(H,2) - 1008*pow(G4xx,2)*G5px*G5xx*pow(H,2) - 
   288*G4x*G4xxx*G5px*G5xx*pow(H,2) + 324*G4xxx*G5p*G5px*G5xx*pow(H,2) + 1188*G4xx*pow(G5px,2)*G5xx*pow(H,2) - 360*pow(G5px,3)*G5xx*pow(H,2) + 1536*G4x*G4xx*G5pxx*G5xx*pow(H,2) - 
   72*G4*G4xxx*G5pxx*G5xx*pow(H,2) - 1344*G4xx*G5p*G5pxx*G5xx*pow(H,2) - 660*G4x*G5px*G5pxx*G5xx*pow(H,2) + 588*G5p*G5px*G5pxx*G5xx*pow(H,2) + 12*G4*pow(G5pxx,2)*G5xx*pow(H,2) - 
   168*pow(G4x,2)*G5pxxx*G5xx*pow(H,2) + 120*G4*G4xx*G5pxxx*G5xx*pow(H,2) + 228*G4x*G5p*G5pxxx*G5xx*pow(H,2) - 72*pow(G5p,2)*G5pxxx*G5xx*pow(H,2) - 72*G4*G5px*G5pxxx*G5xx*pow(H,2) - 
   504*G4px*G4xxx*G5x*G5xx*pow(H,2) + 180*G4xxx*G5pp*G5x*G5xx*pow(H,2) - 1452*G4xx*G5ppx*G5x*G5xx*pow(H,2) + 532*G4x*G5ppxx*G5x*G5xx*pow(H,2) - 332*G5p*G5ppxx*G5x*G5xx*pow(H,2) - 
   90*G3xx*G5px*G5x*G5xx*pow(H,2) + 834*G5ppx*G5px*G5x*G5xx*pow(H,2) + 372*G3x*G5pxx*G5x*G5xx*pow(H,2) - 768*G4px*G5pxx*G5x*G5xx*pow(H,2) - 30*G5pp*G5pxx*G5x*G5xx*pow(H,2) + 
   24*G4p*G5pxxx*G5x*G5xx*pow(H,2) - 144*G3pxx*pow(G5x,2)*G5xx*pow(H,2) + 240*G4ppxx*pow(G5x,2)*G5xx*pow(H,2) + 24*G5pppx*pow(G5x,2)*G5xx*pow(H,2) - 480*G4px*G4xx*pow(G5xx,2)*pow(H,2) + 
   72*G4p*G4xxx*pow(G5xx,2)*pow(H,2) - 12*G4xx*G5pp*pow(G5xx,2)*pow(H,2) - 276*G4x*G5ppx*pow(G5xx,2)*pow(H,2) + 210*G5p*G5ppx*pow(G5xx,2)*pow(H,2) - 
   24*G4*G5ppxx*pow(G5xx,2)*pow(H,2) - 105*G3x*G5px*pow(G5xx,2)*pow(H,2) + 606*G4px*G5px*pow(G5xx,2)*pow(H,2) - 18*G5pp*G5px*pow(G5xx,2)*pow(H,2) - 78*G4p*G5pxx*pow(G5xx,2)*pow(H,2) + 
   69*G3px*G5x*pow(G5xx,2)*pow(H,2) + 54*G4ppx*G5x*pow(G5xx,2)*pow(H,2) - 60*G5ppp*G5x*pow(G5xx,2)*pow(H,2) - 54*G4pp*pow(G5xx,3)*pow(H,2) - 168*G4x*G4xx*G5px*G5xxx*pow(H,2) + 
   228*G4xx*G5p*G5px*G5xxx*pow(H,2) + 72*G4x*pow(G5px,2)*G5xxx*pow(H,2) - 120*G5p*pow(G5px,2)*G5xxx*pow(H,2) + 120*pow(G4x,2)*G5pxx*G5xxx*pow(H,2) - 72*G4*G4xx*G5pxx*G5xxx*pow(H,2) - 
   156*G4x*G5p*G5pxx*G5xxx*pow(H,2) + 48*pow(G5p,2)*G5pxx*G5xxx*pow(H,2) + 48*G4*G5px*G5pxx*G5xxx*pow(H,2) - 360*G4px*G4xx*G5x*G5xxx*pow(H,2) + 108*G4xx*G5pp*G5x*G5xxx*pow(H,2) - 
   316*G4x*G5ppx*G5x*G5xxx*pow(H,2) + 206*G5p*G5ppx*G5x*G5xxx*pow(H,2) - 16*G4*G5ppxx*G5x*G5xxx*pow(H,2) - 30*G3x*G5px*G5x*G5xxx*pow(H,2) + 252*G4px*G5px*G5x*G5xxx*pow(H,2) - 
   48*G5pp*G5px*G5x*G5xxx*pow(H,2) - 12*G4p*G5pxx*G5x*G5xxx*pow(H,2) + 54*G3px*pow(G5x,2)*G5xxx*pow(H,2) - 180*G4ppx*pow(G5x,2)*G5xxx*pow(H,2) + 24*G5ppp*pow(G5x,2)*G5xxx*pow(H,2) - 
   144*G4px*G4x*G5xx*G5xxx*pow(H,2) + 120*G4p*G4xx*G5xx*G5xxx*pow(H,2) + 132*G4px*G5p*G5xx*G5xxx*pow(H,2) + 48*G4x*G5pp*G5xx*G5xxx*pow(H,2) - 54*G5p*G5pp*G5xx*G5xxx*pow(H,2) + 
   36*G4*G5ppx*G5xx*G5xxx*pow(H,2) - 48*G4p*G5px*G5xx*G5xxx*pow(H,2) + 12*G4pp*G5x*G5xx*G5xxx*pow(H,2) - 90*G5pxxx*pow(G5x,3)*pow(H,4) + 261*G5pxx*pow(G5x,2)*G5xx*pow(H,4) + 
   15*G5px*G5x*pow(G5xx,2)*pow(H,4) - 12*G4pxxx*(-72*pow(G4xx,2)*G5p + 84*G4xx*G5p*G5px - 24*G5p*pow(G5px,2) + 8*pow(G4x,2)*G5pxx + 2*pow(G5p,2)*G5pxx - 24*G4pxx*G5p*G5x + 12*G5p*G5ppx*G5x - 
      30*G3x*G5px*G5x + 84*G4px*G5px*G5x - 12*G5pp*G5px*G5x + 6*G3px*pow(G5x,2) - 12*G4ppx*pow(G5x,2) - 12*G3x*G5p*G5xx + 36*G4px*G5p*G5xx - 6*G5p*G5pp*G5xx + 3*G2x*G5x*G5xx - 6*G3p*G5x*G5xx + 
      6*G4pp*G5x*G5xx - 68*G5px*pow(G5x,2)*pow(H,2) - 60*G5p*G5x*G5xx*pow(H,2) - 6*G4*pow(G5xx,2)*pow(H,2) + 6*G4xx*G5x*(9*G3x - 26*G4px + 4*G5pp + 21*G5x*pow(H,2)) + 
      2*G4x*(72*pow(G4xx,2) - 84*G4xx*G5px + 24*pow(G5px,2) - 4*G5p*G5pxx + 24*G4pxx*G5x - 12*G5ppx*G5x + 12*G3x*G5xx - 36*G4px*G5xx + 6*G5pp*G5xx + 45*G5x*G5xx*pow(H,2))) + 
   12*G4pxx*(288*pow(G4xx,3) - 456*pow(G4xx,2)*G5px + 48*G4xxx*G5p*G5px - 36*pow(G5px,3) - 16*G5p*G5px*G5pxx - 8*pow(G4x,2)*G5pxxx - 2*pow(G5p,2)*G5pxxx + 36*G3x*G4xxx*G5x - 
      120*G4px*G4xxx*G5x + 24*G4xxx*G5pp*G5x - 8*G5p*G5ppxx*G5x - 24*G3xx*G5px*G5x + 18*G5ppx*G5px*G5x - 17*G3x*G5pxx*G5x + 70*G4px*G5pxx*G5x - 18*G5pp*G5pxx*G5x - 6*G3pxx*pow(G5x,2) + 
      6*G5pppx*pow(G5x,2) - 6*G3xx*G5p*G5xx + 15*G5p*G5ppx*G5xx - 66*G3x*G5px*G5xx + 156*G4px*G5px*G5xx - 12*G5pp*G5px*G5xx + 2*G2xx*G5x*G5xx + 13*G3px*G5x*G5xx - 18*G4ppx*G5x*G5xx - 6*G5ppp*G5x*G5xx + 
      3*G2x*pow(G5xx,2) - 6*G3p*pow(G5xx,2) + 6*G4pp*pow(G5xx,2) - 6*G3x*G5p*G5xxx + 20*G4px*G5p*G5xxx - 4*G5p*G5pp*G5xxx + 2*G2x*G5x*G5xxx - 4*G3p*G5x*G5xxx + 4*G4pp*G5x*G5xxx + 
      72*G4xxx*pow(G5x,2)*pow(H,2) - 43*G5pxx*pow(G5x,2)*pow(H,2) - 159*G5px*G5x*G5xx*pow(H,2) - 47*G5p*pow(G5xx,2)*pow(H,2) - 30*G5p*G5x*G5xxx*pow(H,2) - 4*G4*G5xx*G5xxx*pow(H,2) + 
      2*G4xx*(72*G4x*G4xxx - 36*G4xxx*G5p + 114*pow(G5px,2) - 14*G4x*G5pxx + 7*G5p*G5pxx + 18*G3xx*G5x - 24*G5ppx*G5x + 66*G3x*G5xx - 174*G4px*G5xx + 21*G5pp*G5xx + 156*G5x*G5xx*pow(H,2)) + 
      2*G4x*(-48*G4xxx*G5px + 16*G5px*G5pxx + 4*G5p*G5pxxx + 8*G5ppxx*G5x + 6*G3xx*G5xx - 15*G5ppx*G5xx + 6*G3x*G5xxx - 20*G4px*G5xxx + 4*G5pp*G5xxx + 26*pow(G5xx,2)*pow(H,2) + 
         24*G5x*G5xxx*pow(H,2)))) + 48*pow(a,16)*dphi*H*(4*G2px*G4*pow(G4p,2) - 8*G3pp*G4*pow(G4p,2) + 4*G2x*pow(G4p,3) - 8*G3p*pow(G4p,3) - 5*G2x*G4*G4p*G4pp + 10*G3p*G4*G4p*G4pp + 
   9*pow(G4p,3)*G4pp + 24*G3px*pow(G4,2)*G4p*pow(H,2) + 18*G3x*G4*pow(G4p,2)*pow(H,2) - 24*G3x*pow(G4,2)*G4pp*pow(H,2) - 72*pow(G4,2)*G4p*G4ppx*pow(H,2) - 
   78*G4*pow(G4p,2)*G4px*pow(H,2) + 72*pow(G4,2)*G4pp*G4px*pow(H,2) - 6*G2x*G4*G4p*G4x*pow(H,2) + 12*G3p*G4*G4p*G4x*pow(H,2) + 18*pow(G4p,3)*G4x*pow(H,2) - 30*G4*G4p*G4pp*G4x*pow(H,2) + 
   6*G2x*G4*G4p*G5p*pow(H,2) - 12*G3p*G4*G4p*G5p*pow(H,2) - 18*pow(G4p,3)*G5p*pow(H,2) + 30*G4*G4p*G4pp*G5p*pow(H,2) + 24*G4*pow(G4p,2)*G5pp*pow(H,2) - 36*G4*G4p*pow(G4x,2)*pow(H,4) + 
   72*G4*G4p*G4x*G5p*pow(H,4) - 36*G4*G4p*pow(G5p,2)*pow(H,4) + 24*pow(G4,2)*G4p*G5px*pow(H,4) + 18*G4*pow(G4p,2)*G5x*pow(H,4) - 24*pow(G4,2)*G4pp*G5x*pow(H,4) - 
   2*G2pp*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
   G2p*(2*G3px*pow(G4,2) + 3*G3x*G4*G4p - 6*pow(G4,2)*G4ppx - 13*G4*G4p*G4px - G2x*G4*G4x + 2*G3p*G4*G4x - 5*pow(G4p,2)*G4x - 4*G4*G4pp*G4x + G2x*G4*G5p - 2*G3p*G4*G5p + 5*pow(G4p,2)*G5p + 
      4*G4*G4pp*G5p + 4*G4*G4p*G5pp - 6*G4*pow(G4x,2)*pow(H,2) + 12*G4*G4x*G5p*pow(H,2) - 6*G4*pow(G5p,2)*pow(H,2) + 2*pow(G4,2)*G5px*pow(H,2) + 3*G4*G4p*G5x*pow(H,2)) - 
   3*G3px*G4*G4p*rptot - 3*G3x*pow(G4p,2)*rptot + 3*G3x*G4*G4pp*rptot + 9*G4*G4p*G4ppx*rptot - G2x*G4*G4px*rptot + 2*G3p*G4*G4px*rptot + 12*pow(G4p,2)*G4px*rptot - 9*G4*G4pp*G4px*rptot + 
   G2px*G4*G4x*rptot - 2*G3pp*G4*G4x*rptot + G2x*G4p*G4x*rptot - 2*G3p*G4p*G4x*rptot + 6*G4p*G4pp*G4x*rptot - G2px*G4*G5p*rptot + 2*G3pp*G4*G5p*rptot - G2x*G4p*G5p*rptot + 2*G3p*G4p*G5p*rptot - 
   6*G4p*G4pp*G5p*rptot + G2x*G4*G5pp*rptot - 2*G3p*G4*G5pp*rptot - 3*pow(G4p,2)*G5pp*rptot + 6*G4p*pow(G4x,2)*pow(H,2)*rptot - 12*G4p*G4x*G5p*pow(H,2)*rptot + 
   6*G4p*pow(G5p,2)*pow(H,2)*rptot - 3*G4*G4p*G5px*pow(H,2)*rptot - 3*pow(G4p,2)*G5x*pow(H,2)*rptot + 3*G4*G4pp*G5x*pow(H,2)*rptot) + 
pow(a,6)*pow(dphi,11)*pow(H,3)*(1728*G4px*pow(G4pxx,2)*G4x - 1440*pow(G4px,2)*G4pxxx*G4x - 288*G3pxx*G4pxx*pow(G4x,2) - 288*G3px*G4pxxx*pow(G4x,2) + 576*G4ppx*G4pxxx*pow(G4x,2) + 
   4032*pow(G4px,2)*G4pxx*G4xx + 1872*G3pxx*G4px*G4x*G4xx - 2592*G4ppxx*G4px*G4x*G4xx + 288*G2xx*G4pxx*G4x*G4xx + 1008*G3px*G4pxx*G4x*G4xx - 1440*G4ppx*G4pxx*G4x*G4xx - 288*G2x*G4pxxx*G4x*G4xx + 
   576*G3p*G4pxxx*G4x*G4xx - 576*G4pp*G4pxxx*G4x*G4xx - 288*G2xx*G4px*pow(G4xx,2) - 2592*G3px*G4px*pow(G4xx,2) + 2304*G4ppx*G4px*pow(G4xx,2) + 576*G2x*G4pxx*pow(G4xx,2) - 
   1152*G3p*G4pxx*pow(G4xx,2) + 1152*G4pp*G4pxx*pow(G4xx,2) - 288*G2pxx*G4x*pow(G4xx,2) - 288*G3ppx*G4x*pow(G4xx,2) + 1152*G4pppx*G4x*pow(G4xx,2) + 576*G3pp*pow(G4xx,3) - 
   1152*G4ppp*pow(G4xx,3) - 1440*pow(G4px,3)*G4xxx - 720*G3px*G4px*G4x*G4xxx + 2016*G4ppx*G4px*G4x*G4xxx + 288*G2x*G4pxx*G4x*G4xxx - 576*G3p*G4pxx*G4x*G4xxx + 576*G4pp*G4pxx*G4x*G4xxx + 
   288*G3ppx*pow(G4x,2)*G4xxx - 576*G4pppx*pow(G4x,2)*G4xxx - 288*G2x*G4px*G4xx*G4xxx + 576*G3p*G4px*G4xx*G4xxx - 576*G4pp*G4px*G4xx*G4xxx - 288*G3pp*G4x*G4xx*G4xxx + 576*G4ppp*G4x*G4xx*G4xxx - 
   864*G4px*pow(G4pxx,2)*G5p + 720*pow(G4px,2)*G4pxxx*G5p + 288*G3pxx*G4pxx*G4x*G5p + 288*G3px*G4pxxx*G4x*G5p - 576*G4ppx*G4pxxx*G4x*G5p - 936*G3pxx*G4px*G4xx*G5p + 1296*G4ppxx*G4px*G4xx*G5p - 
   144*G2xx*G4pxx*G4xx*G5p - 504*G3px*G4pxx*G4xx*G5p + 720*G4ppx*G4pxx*G4xx*G5p + 144*G2x*G4pxxx*G4xx*G5p - 288*G3p*G4pxxx*G4xx*G5p + 288*G4pp*G4pxxx*G4xx*G5p + 144*G2pxx*pow(G4xx,2)*G5p + 
   144*G3ppx*pow(G4xx,2)*G5p - 576*G4pppx*pow(G4xx,2)*G5p + 360*G3px*G4px*G4xxx*G5p - 1008*G4ppx*G4px*G4xxx*G5p - 144*G2x*G4pxx*G4xxx*G5p + 288*G3p*G4pxx*G4xxx*G5p - 288*G4pp*G4pxx*G4xxx*G5p - 
   288*G3ppx*G4x*G4xxx*G5p + 576*G4pppx*G4x*G4xxx*G5p + 144*G3pp*G4xx*G4xxx*G5p - 288*G4ppp*G4xx*G4xxx*G5p - 72*G3pxx*G4pxx*pow(G5p,2) - 72*G3px*G4pxxx*pow(G5p,2) + 144*G4ppx*G4pxxx*pow(G5p,2) + 
   72*G3ppx*G4xxx*pow(G5p,2) - 144*G4pppx*G4xxx*pow(G5p,2) - 576*pow(G4pxx,2)*G4x*G5pp + 288*G4px*G4pxxx*G4x*G5pp + 720*G4px*G4pxx*G4xx*G5pp - 288*G3pxx*G4x*G4xx*G5pp + 576*G4ppxx*G4x*G4xx*G5pp + 
   144*G2xx*pow(G4xx,2)*G5pp + 432*G3px*pow(G4xx,2)*G5pp - 1152*G4ppx*pow(G4xx,2)*G5pp + 1008*pow(G4px,2)*G4xxx*G5pp + 144*G3px*G4x*G4xxx*G5pp - 288*G4ppx*G4x*G4xxx*G5pp + 
   144*G2x*G4xx*G4xxx*G5pp - 288*G3p*G4xx*G4xxx*G5pp + 288*G4pp*G4xx*G4xxx*G5pp + 288*pow(G4pxx,2)*G5p*G5pp - 144*G4px*G4pxxx*G5p*G5pp + 144*G3pxx*G4xx*G5p*G5pp - 288*G4ppxx*G4xx*G5p*G5pp - 
   72*G3px*G4xxx*G5p*G5pp + 144*G4ppx*G4xxx*G5p*G5pp - 288*G4pxx*G4xx*pow(G5pp,2) - 144*G4px*G4xxx*pow(G5pp,2) - 576*G4pxx*G4x*G4xx*G5ppp + 1728*G4px*pow(G4xx,2)*G5ppp - 288*G4px*G4x*G4xxx*G5ppp + 
   288*G4pxx*G4xx*G5p*G5ppp + 144*G4px*G4xxx*G5p*G5ppp + 288*G4pxx*pow(G4x,2)*G5pppx - 576*G4px*G4x*G4xx*G5pppx - 288*G4pxx*G4x*G5p*G5pppx + 288*G4px*G4xx*G5p*G5pppx + 72*G4pxx*pow(G5p,2)*G5pppx - 
   432*G4px*G4pxx*G4x*G5ppx + 144*G3pxx*pow(G4x,2)*G5ppx - 288*G4ppxx*pow(G4x,2)*G5ppx - 3456*pow(G4px,2)*G4xx*G5ppx - 192*G2xx*G4x*G4xx*G5ppx - 672*G3px*G4x*G4xx*G5ppx + 
   1728*G4ppx*G4x*G4xx*G5ppx - 288*G2x*pow(G4xx,2)*G5ppx + 576*G3p*pow(G4xx,2)*G5ppx - 576*G4pp*pow(G4xx,2)*G5ppx - 144*G2x*G4x*G4xxx*G5ppx + 288*G3p*G4x*G4xxx*G5ppx - 288*G4pp*G4x*G4xxx*G5ppx + 
   216*G4px*G4pxx*G5p*G5ppx - 144*G3pxx*G4x*G5p*G5ppx + 288*G4ppxx*G4x*G5p*G5ppx + 96*G2xx*G4xx*G5p*G5ppx + 336*G3px*G4xx*G5p*G5ppx - 864*G4ppx*G4xx*G5p*G5ppx + 72*G2x*G4xxx*G5p*G5ppx - 
   144*G3p*G4xxx*G5p*G5ppx + 144*G4pp*G4xxx*G5p*G5ppx + 36*G3pxx*pow(G5p,2)*G5ppx - 72*G4ppxx*pow(G5p,2)*G5ppx + 288*G4pxx*G4x*G5pp*G5ppx + 864*G4px*G4xx*G5pp*G5ppx - 144*G4pxx*G5p*G5pp*G5ppx - 
   288*G4px*G4x*pow(G5ppx,2) + 144*G4px*G5p*pow(G5ppx,2) + 672*pow(G4px,2)*G4x*G5ppxx + 16*G2xx*pow(G4x,2)*G5ppxx + 128*G3px*pow(G4x,2)*G5ppxx - 288*G4ppx*pow(G4x,2)*G5ppxx + 
   144*G2x*G4x*G4xx*G5ppxx - 288*G3p*G4x*G4xx*G5ppxx + 288*G4pp*G4x*G4xx*G5ppxx - 336*pow(G4px,2)*G5p*G5ppxx - 16*G2xx*G4x*G5p*G5ppxx - 128*G3px*G4x*G5p*G5ppxx + 288*G4ppx*G4x*G5p*G5ppxx - 
   72*G2x*G4xx*G5p*G5ppxx + 144*G3p*G4xx*G5p*G5ppxx - 144*G4pp*G4xx*G5p*G5ppxx + 4*G2xx*pow(G5p,2)*G5ppxx + 32*G3px*pow(G5p,2)*G5ppxx - 72*G4ppx*pow(G5p,2)*G5ppxx - 144*G4px*G4x*G5pp*G5ppxx + 
   72*G4px*G5p*G5pp*G5ppxx - 1440*pow(G4px,2)*G4pxx*G5px - 1008*G3pxx*G4px*G4x*G5px + 1440*G4ppxx*G4px*G4x*G5px - 192*G2xx*G4pxx*G4x*G5px - 384*G3px*G4pxx*G4x*G5px + 576*G4ppx*G4pxx*G4x*G5px + 
   144*G2x*G4pxxx*G4x*G5px - 288*G3p*G4pxxx*G4x*G5px + 288*G4pp*G4pxxx*G4x*G5px + 648*G2xx*G4px*G4xx*G5px + 2664*G3px*G4px*G4xx*G5px - 3168*G4ppx*G4px*G4xx*G5px - 432*G2x*G4pxx*G4xx*G5px + 
   864*G3p*G4pxx*G4xx*G5px - 864*G4pp*G4pxx*G4xx*G5px + 336*G2pxx*G4x*G4xx*G5px + 240*G3ppx*G4x*G4xx*G5px - 1152*G4pppx*G4x*G4xx*G5px - 96*G2px*pow(G4xx,2)*G5px - 768*G3pp*pow(G4xx,2)*G5px + 
   1728*G4ppp*pow(G4xx,2)*G5px + 216*G2x*G4px*G4xxx*G5px - 432*G3p*G4px*G4xxx*G5px + 432*G4pp*G4px*G4xxx*G5px + 144*G3pp*G4x*G4xxx*G5px - 288*G4ppp*G4x*G4xxx*G5px + 504*G3pxx*G4px*G5p*G5px - 
   720*G4ppxx*G4px*G5p*G5px + 96*G2xx*G4pxx*G5p*G5px + 192*G3px*G4pxx*G5p*G5px - 288*G4ppx*G4pxx*G5p*G5px - 72*G2x*G4pxxx*G5p*G5px + 144*G3p*G4pxxx*G5p*G5px - 144*G4pp*G4pxxx*G5p*G5px - 
   168*G2pxx*G4xx*G5p*G5px - 120*G3ppx*G4xx*G5p*G5px + 576*G4pppx*G4xx*G5p*G5px - 72*G3pp*G4xxx*G5p*G5px + 144*G4ppp*G4xxx*G5p*G5px - 576*G4px*G4pxx*G5pp*G5px + 144*G3pxx*G4x*G5pp*G5px - 
   288*G4ppxx*G4x*G5pp*G5px - 216*G2xx*G4xx*G5pp*G5px - 360*G3px*G4xx*G5pp*G5px + 1152*G4ppx*G4xx*G5pp*G5px - 72*G2x*G4xxx*G5pp*G5px + 144*G3p*G4xxx*G5pp*G5px - 144*G4pp*G4xxx*G5pp*G5px - 
   72*G3pxx*G5p*G5pp*G5px + 144*G4ppxx*G5p*G5pp*G5px + 144*G4pxx*pow(G5pp,2)*G5px + 288*G4pxx*G4x*G5ppp*G5px - 1728*G4px*G4xx*G5ppp*G5px - 144*G4pxx*G5p*G5ppp*G5px + 288*G4px*G4x*G5pppx*G5px - 
   144*G4px*G5p*G5pppx*G5px + 1872*pow(G4px,2)*G5ppx*G5px + 120*G2xx*G4x*G5ppx*G5px + 312*G3px*G4x*G5ppx*G5px - 864*G4ppx*G4x*G5ppx*G5px + 288*G2x*G4xx*G5ppx*G5px - 576*G3p*G4xx*G5ppx*G5px + 
   576*G4pp*G4xx*G5ppx*G5px - 60*G2xx*G5p*G5ppx*G5px - 156*G3px*G5p*G5ppx*G5px + 432*G4ppx*G5p*G5ppx*G5px - 432*G4px*G5pp*G5ppx*G5px - 72*G2x*G4x*G5ppxx*G5px + 144*G3p*G4x*G5ppxx*G5px - 
   144*G4pp*G4x*G5ppxx*G5px + 36*G2x*G5p*G5ppxx*G5px - 72*G3p*G5p*G5ppxx*G5px + 72*G4pp*G5p*G5ppxx*G5px - 264*G2xx*G4px*pow(G5px,2) - 672*G3px*G4px*pow(G5px,2) + 1008*G4ppx*G4px*pow(G5px,2) + 
   72*G2x*G4pxx*pow(G5px,2) - 144*G3p*G4pxx*pow(G5px,2) + 144*G4pp*G4pxx*pow(G5px,2) - 96*G2pxx*G4x*pow(G5px,2) - 48*G3ppx*G4x*pow(G5px,2) + 288*G4pppx*G4x*pow(G5px,2) + 
   96*G2px*G4xx*pow(G5px,2) + 336*G3pp*G4xx*pow(G5px,2) - 864*G4ppp*G4xx*pow(G5px,2) + 48*G2pxx*G5p*pow(G5px,2) + 24*G3ppx*G5p*pow(G5px,2) - 144*G4pppx*G5p*pow(G5px,2) + 
   72*G2xx*G5pp*pow(G5px,2) + 72*G3px*G5pp*pow(G5px,2) - 288*G4ppx*G5pp*pow(G5px,2) + 432*G4px*G5ppp*pow(G5px,2) - 72*G2x*G5ppx*pow(G5px,2) + 144*G3p*G5ppx*pow(G5px,2) - 
   144*G4pp*G5ppx*pow(G5px,2) - 24*G2px*pow(G5px,3) - 48*G3pp*pow(G5px,3) + 144*G4ppp*pow(G5px,3) + 27*pow(G3x,3)*G5pxx + 312*pow(G4px,3)*G5pxx - 120*G2xx*G4px*G4x*G5pxx + 
   288*G3px*G4px*G4x*G5pxx - 624*G4ppx*G4px*G4x*G5pxx - 216*G2x*G4pxx*G4x*G5pxx + 432*G3p*G4pxx*G4x*G5pxx - 432*G4pp*G4pxx*G4x*G5pxx - 16*G2pxx*pow(G4x,2)*G5pxx - 128*G3ppx*pow(G4x,2)*G5pxx + 
   288*G4pppx*pow(G4x,2)*G5pxx + 96*G2px*G4x*G4xx*G5pxx + 48*G3pp*G4x*G4xx*G5pxx - 288*G4ppp*G4x*G4xx*G5pxx + 60*G2xx*G4px*G5p*G5pxx - 144*G3px*G4px*G5p*G5pxx + 312*G4ppx*G4px*G5p*G5pxx + 
   108*G2x*G4pxx*G5p*G5pxx - 216*G3p*G4pxx*G5p*G5pxx + 216*G4pp*G4pxx*G5p*G5pxx + 16*G2pxx*G4x*G5p*G5pxx + 128*G3ppx*G4x*G5p*G5pxx - 288*G4pppx*G4x*G5p*G5pxx - 48*G2px*G4xx*G5p*G5pxx - 
   24*G3pp*G4xx*G5p*G5pxx + 144*G4ppp*G4xx*G5p*G5pxx - 4*G2pxx*pow(G5p,2)*G5pxx - 32*G3ppx*pow(G5p,2)*G5pxx + 72*G4pppx*pow(G5p,2)*G5pxx - 408*pow(G4px,2)*G5pp*G5pxx + 24*G2xx*G4x*G5pp*G5pxx - 
   96*G3px*G4x*G5pp*G5pxx + 144*G4ppx*G4x*G5pp*G5pxx - 72*G2x*G4xx*G5pp*G5pxx + 144*G3p*G4xx*G5pp*G5pxx - 144*G4pp*G4xx*G5pp*G5pxx - 12*G2xx*G5p*G5pp*G5pxx + 48*G3px*G5p*G5pp*G5pxx - 
   72*G4ppx*G5p*G5pp*G5pxx + 72*G4px*pow(G5pp,2)*G5pxx + 144*G4px*G4x*G5ppp*G5pxx - 72*G4px*G5p*G5ppp*G5pxx + 72*G2x*G4x*G5ppx*G5pxx - 144*G3p*G4x*G5ppx*G5pxx + 144*G4pp*G4x*G5ppx*G5pxx - 
   36*G2x*G5p*G5ppx*G5pxx + 72*G3p*G5p*G5ppx*G5pxx - 72*G4pp*G5p*G5ppx*G5pxx - 36*G2x*G4px*G5px*G5pxx + 72*G3p*G4px*G5px*G5pxx - 72*G4pp*G4px*G5px*G5pxx - 48*G2px*G4x*G5px*G5pxx - 
   24*G3pp*G4x*G5px*G5pxx + 144*G4ppp*G4x*G5px*G5pxx + 24*G2px*G5p*G5px*G5pxx + 12*G3pp*G5p*G5px*G5pxx - 72*G4ppp*G5p*G5px*G5pxx + 36*G2x*G5pp*G5px*G5pxx - 72*G3p*G5pp*G5px*G5pxx + 
   72*G4pp*G5pp*G5px*G5pxx + 24*G2x*G4px*G4x*G5pxxx - 48*G3p*G4px*G4x*G5pxxx + 48*G4pp*G4px*G4x*G5pxxx - 16*G2px*pow(G4x,2)*G5pxxx + 16*G3pp*pow(G4x,2)*G5pxxx - 12*G2x*G4px*G5p*G5pxxx + 
   24*G3p*G4px*G5p*G5pxxx - 24*G4pp*G4px*G5p*G5pxxx + 16*G2px*G4x*G5p*G5pxxx - 16*G3pp*G4x*G5p*G5pxxx - 4*G2px*pow(G5p,2)*G5pxxx + 4*G3pp*pow(G5p,2)*G5pxxx - 360*G3pxx*pow(G4px,2)*G5x + 
   576*G4ppxx*pow(G4px,2)*G5x - 240*G2xx*G4px*G4pxx*G5x + 168*G3px*G4px*G4pxx*G5x - 144*G4ppx*G4px*G4pxx*G5x - 144*G2x*pow(G4pxx,2)*G5x + 288*G3p*pow(G4pxx,2)*G5x - 288*G4pp*pow(G4pxx,2)*G5x + 
   72*G2x*G4px*G4pxxx*G5x - 144*G3p*G4px*G4pxxx*G5x + 144*G4pp*G4px*G4pxxx*G5x - 144*G3px*G3pxx*G4x*G5x + 288*G3pxx*G4ppx*G4x*G5x + 96*G2xx*G4ppxx*G4x*G5x + 192*G3px*G4ppxx*G4x*G5x - 
   576*G4ppx*G4ppxx*G4x*G5x - 96*G2pxx*G4pxx*G4x*G5x - 192*G3ppx*G4pxx*G4x*G5x + 576*G4pppx*G4pxx*G4x*G5x - 96*G2px*G4pxxx*G4x*G5x + 96*G3pp*G4pxxx*G4x*G5x + 36*G2xx*G3px*G4xx*G5x + 
   108*pow(G3px,2)*G4xx*G5x - 72*G2x*G3pxx*G4xx*G5x + 144*G3p*G3pxx*G4xx*G5x - 144*G3pxx*G4pp*G4xx*G5x - 168*G2xx*G4ppx*G4xx*G5x - 408*G3px*G4ppx*G4xx*G5x + 576*pow(G4ppx,2)*G4xx*G5x + 
   144*G2x*G4ppxx*G4xx*G5x - 288*G3p*G4ppxx*G4xx*G5x + 288*G4pp*G4ppxx*G4xx*G5x + 312*G2pxx*G4px*G4xx*G5x - 24*G3ppx*G4px*G4xx*G5x - 576*G4pppx*G4px*G4xx*G5x + 192*G2px*G4pxx*G4xx*G5x - 
   48*G3pp*G4pxx*G4xx*G5x - 288*G4ppp*G4pxx*G4xx*G5x - 96*G2ppx*pow(G4xx,2)*G5x + 96*G3ppp*pow(G4xx,2)*G5x + 36*G2x*G3px*G4xxx*G5x - 72*G3p*G3px*G4xxx*G5x + 72*G3px*G4pp*G4xxx*G5x - 
   72*G2x*G4ppx*G4xxx*G5x + 144*G3p*G4ppx*G4xxx*G5x - 144*G4pp*G4ppx*G4xxx*G5x + 72*G3pp*G4px*G4xxx*G5x - 144*G4ppp*G4px*G4xxx*G5x + 96*G2ppx*G4x*G4xxx*G5x - 96*G3ppp*G4x*G4xxx*G5x + 
   72*G3px*G3pxx*G5p*G5x - 144*G3pxx*G4ppx*G5p*G5x - 48*G2xx*G4ppxx*G5p*G5x - 96*G3px*G4ppxx*G5p*G5x + 288*G4ppx*G4ppxx*G5p*G5x + 48*G2pxx*G4pxx*G5p*G5x + 96*G3ppx*G4pxx*G5p*G5x - 
   288*G4pppx*G4pxx*G5p*G5x + 48*G2px*G4pxxx*G5p*G5x - 48*G3pp*G4pxxx*G5p*G5x - 48*G2ppx*G4xxx*G5p*G5x + 48*G3ppp*G4xxx*G5p*G5x + 72*G3pxx*G4px*G5pp*G5x - 144*G4ppxx*G4px*G5pp*G5x + 
   48*G2xx*G4pxx*G5pp*G5x - 120*G3px*G4pxx*G5pp*G5x + 144*G4ppx*G4pxx*G5pp*G5x - 48*G2pxx*G4xx*G5pp*G5x + 48*G3ppx*G4xx*G5pp*G5x + 144*G4px*G4pxx*G5ppp*G5x + 48*G2xx*G4xx*G5ppp*G5x - 
   48*G3px*G4xx*G5ppp*G5x + 72*pow(G4px,2)*G5pppx*G5x - 48*G2xx*G4x*G5pppx*G5x + 48*G3px*G4x*G5pppx*G5x + 24*G2xx*G5p*G5pppx*G5x - 24*G3px*G5p*G5pppx*G5x + 132*G2xx*G4px*G5ppx*G5x + 
   84*G3px*G4px*G5ppx*G5x - 432*G4ppx*G4px*G5ppx*G5x + 72*G2x*G4pxx*G5ppx*G5x - 144*G3p*G4pxx*G5ppx*G5x + 144*G4pp*G4pxx*G5ppx*G5x + 48*G2pxx*G4x*G5ppx*G5x - 48*G3ppx*G4x*G5ppx*G5x - 
   96*G2px*G4xx*G5ppx*G5x + 96*G3pp*G4xx*G5ppx*G5x - 24*G2pxx*G5p*G5ppx*G5x + 24*G3ppx*G5p*G5ppx*G5x - 24*G2xx*G5pp*G5ppx*G5x + 24*G3px*G5pp*G5ppx*G5x - 36*G2x*G4px*G5ppxx*G5x + 72*G3p*G4px*G5ppxx*G5x - 
   72*G4pp*G4px*G5ppxx*G5x + 48*G2px*G4x*G5ppxx*G5x - 48*G3pp*G4x*G5ppxx*G5x - 24*G2px*G5p*G5ppxx*G5x + 24*G3pp*G5p*G5ppxx*G5x - 24*G2xx*G3px*G5px*G5x - 48*pow(G3px,2)*G5px*G5x + 
   36*G2x*G3pxx*G5px*G5x - 72*G3p*G3pxx*G5px*G5x + 72*G3pxx*G4pp*G5px*G5x + 96*G2xx*G4ppx*G5px*G5x + 192*G3px*G4ppx*G5px*G5x - 288*pow(G4ppx,2)*G5px*G5x - 72*G2x*G4ppxx*G5px*G5x + 
   144*G3p*G4ppxx*G5px*G5x - 144*G4pp*G4ppxx*G5px*G5x - 168*G2pxx*G4px*G5px*G5x + 24*G3ppx*G4px*G5px*G5x + 288*G4pppx*G4px*G5px*G5x - 96*G2px*G4pxx*G5px*G5x + 24*G3pp*G4pxx*G5px*G5x + 
   144*G4ppp*G4pxx*G5px*G5x + 96*G2ppx*G4xx*G5px*G5x - 96*G3ppp*G4xx*G5px*G5x + 24*G2pxx*G5pp*G5px*G5x - 24*G3ppx*G5pp*G5px*G5x - 24*G2xx*G5ppp*G5px*G5x + 24*G3px*G5ppp*G5px*G5x + 
   48*G2px*G5ppx*G5px*G5x - 48*G3pp*G5ppx*G5px*G5x - 24*G2ppx*pow(G5px,2)*G5x + 24*G3ppp*pow(G5px,2)*G5x + 6*G2x*G2xx*G5pxx*G5x - 12*G2xx*G3p*G5pxx*G5x - 24*G2x*G3px*G5pxx*G5x + 
   48*G3p*G3px*G5pxx*G5x + 12*G2xx*G4pp*G5pxx*G5x - 48*G3px*G4pp*G5pxx*G5x + 36*G2x*G4ppx*G5pxx*G5x - 72*G3p*G4ppx*G5pxx*G5x + 72*G4pp*G4ppx*G5pxx*G5x - 24*G2px*G4px*G5pxx*G5x - 12*G3pp*G4px*G5pxx*G5x + 
   72*G4ppp*G4px*G5pxx*G5x - 48*G2ppx*G4x*G5pxx*G5x + 48*G3ppp*G4x*G5pxx*G5x + 24*G2ppx*G5p*G5pxx*G5x - 24*G3ppp*G5p*G5pxx*G5x + 12*G2xx*G3ppx*pow(G5x,2) - 12*G2pxx*G3px*pow(G5x,2) - 
   12*G2px*G3pxx*pow(G5x,2) + 12*G3pp*G3pxx*pow(G5x,2) - 24*G2xx*G4pppx*pow(G5x,2) + 24*G3px*G4pppx*pow(G5x,2) + 24*G2pxx*G4ppx*pow(G5x,2) - 24*G3ppx*G4ppx*pow(G5x,2) + 
   24*G2px*G4ppxx*pow(G5x,2) - 24*G3pp*G4ppxx*pow(G5x,2) - 24*G2ppx*G4pxx*pow(G5x,2) + 24*G3ppp*G4pxx*pow(G5x,2) + 144*G2xx*pow(G4px,2)*G5xx + 396*G3px*pow(G4px,2)*G5xx - 
   648*G4ppx*pow(G4px,2)*G5xx - 72*G2x*G4px*G4pxx*G5xx + 144*G3p*G4px*G4pxx*G5xx - 144*G4pp*G4px*G4pxx*G5xx + 12*G2xx*G3px*G4x*G5xx + 60*pow(G3px,2)*G4x*G5xx - 36*G2x*G3pxx*G4x*G5xx + 
   72*G3p*G3pxx*G4x*G5xx - 72*G3pxx*G4pp*G4x*G5xx - 72*G2xx*G4ppx*G4x*G5xx - 216*G3px*G4ppx*G4x*G5xx + 288*pow(G4ppx,2)*G4x*G5xx + 72*G2x*G4ppxx*G4x*G5xx - 144*G3p*G4ppxx*G4x*G5xx + 
   144*G4pp*G4ppxx*G4x*G5xx + 144*G2pxx*G4px*G4x*G5xx - 288*G4pppx*G4px*G4x*G5xx + 96*G2px*G4pxx*G4x*G5xx - 24*G3pp*G4pxx*G4x*G5xx - 144*G4ppp*G4pxx*G4x*G5xx + 72*G2x*G3px*G4xx*G5xx - 
   144*G3p*G3px*G4xx*G5xx + 144*G3px*G4pp*G4xx*G5xx - 144*G2x*G4ppx*G4xx*G5xx + 288*G3p*G4ppx*G4xx*G5xx - 288*G4pp*G4ppx*G4xx*G5xx - 96*G2px*G4px*G4xx*G5xx - 336*G3pp*G4px*G4xx*G5xx + 
   864*G4ppp*G4px*G4xx*G5xx - 96*G2ppx*G4x*G4xx*G5xx + 96*G3ppp*G4x*G4xx*G5xx - 6*G2xx*G3px*G5p*G5xx - 30*pow(G3px,2)*G5p*G5xx + 18*G2x*G3pxx*G5p*G5xx - 36*G3p*G3pxx*G5p*G5xx + 36*G3pxx*G4pp*G5p*G5xx + 
   36*G2xx*G4ppx*G5p*G5xx + 108*G3px*G4ppx*G5p*G5xx - 144*pow(G4ppx,2)*G5p*G5xx - 36*G2x*G4ppxx*G5p*G5xx + 72*G3p*G4ppxx*G5p*G5xx - 72*G4pp*G4ppxx*G5p*G5xx - 72*G2pxx*G4px*G5p*G5xx + 
   144*G4pppx*G4px*G5p*G5xx - 48*G2px*G4pxx*G5p*G5xx + 12*G3pp*G4pxx*G5p*G5xx + 72*G4ppp*G4pxx*G5p*G5xx + 48*G2ppx*G4xx*G5p*G5xx - 48*G3ppp*G4xx*G5p*G5xx - 96*G2xx*G4px*G5pp*G5xx - 
   48*G3px*G4px*G5pp*G5xx + 288*G4ppx*G4px*G5pp*G5xx - 36*G2x*G4pxx*G5pp*G5xx + 72*G3p*G4pxx*G5pp*G5xx - 72*G4pp*G4pxx*G5pp*G5xx - 24*G2pxx*G4x*G5pp*G5xx + 24*G3ppx*G4x*G5pp*G5xx + 
   48*G2px*G4xx*G5pp*G5xx - 48*G3pp*G4xx*G5pp*G5xx + 12*G2pxx*G5p*G5pp*G5xx - 12*G3ppx*G5p*G5pp*G5xx + 12*G2xx*pow(G5pp,2)*G5xx - 12*G3px*pow(G5pp,2)*G5xx - 216*pow(G4px,2)*G5ppp*G5xx + 
   24*G2xx*G4x*G5ppp*G5xx - 24*G3px*G4x*G5ppp*G5xx - 12*G2xx*G5p*G5ppp*G5xx + 12*G3px*G5p*G5ppp*G5xx + 72*G2x*G4px*G5ppx*G5xx - 144*G3p*G4px*G5ppx*G5xx + 144*G4pp*G4px*G5ppx*G5xx - 
   48*G2px*G4x*G5ppx*G5xx + 48*G3pp*G4x*G5ppx*G5xx + 24*G2px*G5p*G5ppx*G5xx - 24*G3pp*G5p*G5ppx*G5xx - 6*G2x*G2xx*G5px*G5xx + 12*G2xx*G3p*G5px*G5xx - 30*G2x*G3px*G5px*G5xx + 60*G3p*G3px*G5px*G5xx - 
   12*G2xx*G4pp*G5px*G5xx - 60*G3px*G4pp*G5px*G5xx + 72*G2x*G4ppx*G5px*G5xx - 144*G3p*G4ppx*G5px*G5xx + 144*G4pp*G4ppx*G5px*G5xx + 72*G2px*G4px*G5px*G5xx + 144*G3pp*G4px*G5px*G5xx - 
   432*G4ppp*G4px*G5px*G5xx + 48*G2ppx*G4x*G5px*G5xx - 48*G3ppp*G4x*G5px*G5xx - 24*G2ppx*G5p*G5px*G5xx + 24*G3ppp*G5p*G5px*G5xx - 24*G2px*G5pp*G5px*G5xx + 24*G3pp*G5pp*G5px*G5xx - 6*G2pxx*G2x*G5x*G5xx + 
   12*G2pxx*G3p*G5x*G5xx - 6*G2xx*G3pp*G5x*G5xx + 6*G2x*G3ppx*G5x*G5xx - 12*G3p*G3ppx*G5x*G5xx + 12*G2px*G3px*G5x*G5xx - 6*G3pp*G3px*G5x*G5xx - 12*G2pxx*G4pp*G5x*G5xx + 12*G3ppx*G4pp*G5x*G5xx + 
   12*G2xx*G4ppp*G5x*G5xx - 12*G3px*G4ppp*G5x*G5xx - 24*G2px*G4ppx*G5x*G5xx + 24*G3pp*G4ppx*G5x*G5xx + 24*G2ppx*G4px*G5x*G5xx - 24*G3ppp*G4px*G5x*G5xx + 24*G2x*pow(G4px,2)*G5xxx - 
   48*G3p*pow(G4px,2)*G5xxx + 48*G4pp*pow(G4px,2)*G5xxx + 12*G2x*G3px*G4x*G5xxx - 24*G3p*G3px*G4x*G5xxx + 24*G3px*G4pp*G4x*G5xxx - 24*G2x*G4ppx*G4x*G5xxx + 48*G3p*G4ppx*G4x*G5xxx - 
   48*G4pp*G4ppx*G4x*G5xxx + 24*G3pp*G4px*G4x*G5xxx - 48*G4ppp*G4px*G4x*G5xxx + 16*G2ppx*pow(G4x,2)*G5xxx - 16*G3ppp*pow(G4x,2)*G5xxx - 6*G2x*G3px*G5p*G5xxx + 12*G3p*G3px*G5p*G5xxx - 
   12*G3px*G4pp*G5p*G5xxx + 12*G2x*G4ppx*G5p*G5xxx - 24*G3p*G4ppx*G5p*G5xxx + 24*G4pp*G4ppx*G5p*G5xxx - 12*G3pp*G4px*G5p*G5xxx + 24*G4ppp*G4px*G5p*G5xxx - 16*G2ppx*G4x*G5p*G5xxx + 
   16*G3ppp*G4x*G5p*G5xxx + 4*G2ppx*pow(G5p,2)*G5xxx - 4*G3ppp*pow(G5p,2)*G5xxx - 12*G2x*G4px*G5pp*G5xxx + 24*G3p*G4px*G5pp*G5xxx - 24*G4pp*G4px*G5pp*G5xxx - 
   5184*G4pxxx*pow(G4x,2)*G4xx*pow(H,2) + 1728*G4*G4pxxx*pow(G4xx,2)*pow(H,2) + 13824*G4pxx*G4x*pow(G4xx,2)*pow(H,2) - 3456*G4px*pow(G4xx,3)*pow(H,2) + 
   3456*G4pxx*pow(G4x,2)*G4xxx*pow(H,2) - 1728*G4*G4pxx*G4xx*G4xxx*pow(H,2) - 3456*G4px*G4x*G4xx*G4xxx*pow(H,2) + 1728*G4p*pow(G4xx,2)*G4xxx*pow(H,2) + 6912*G4pxxx*G4x*G4xx*G5p*pow(H,2) - 
   12096*G4pxx*pow(G4xx,2)*G5p*pow(H,2) - 4320*G4pxx*G4x*G4xxx*G5p*pow(H,2) + 3456*G4px*G4xx*G4xxx*G5p*pow(H,2) - 2160*G4pxxx*G4xx*pow(G5p,2)*pow(H,2) + 
   1296*G4pxx*G4xxx*pow(G5p,2)*pow(H,2) + 864*G4x*G4xx*G4xxx*G5pp*pow(H,2) - 1296*G4xx*G4xxx*G5p*G5pp*pow(H,2) - 5856*G4x*pow(G4xx,2)*G5ppx*pow(H,2) - 
   1824*pow(G4x,2)*G4xxx*G5ppx*pow(H,2) + 1152*G4*G4xx*G4xxx*G5ppx*pow(H,2) + 4656*pow(G4xx,2)*G5p*G5ppx*pow(H,2) + 2400*G4x*G4xxx*G5p*G5ppx*pow(H,2) - 
   744*G4xxx*pow(G5p,2)*G5ppx*pow(H,2) + 2592*pow(G4x,2)*G4xx*G5ppxx*pow(H,2) - 672*G4*pow(G4xx,2)*G5ppxx*pow(H,2) - 192*G4*G4x*G4xxx*G5ppxx*pow(H,2) - 3264*G4x*G4xx*G5p*G5ppxx*pow(H,2) + 
   96*G4*G4xxx*G5p*G5ppxx*pow(H,2) + 984*G4xx*pow(G5p,2)*G5ppxx*pow(H,2) + 2688*G4pxxx*pow(G4x,2)*G5px*pow(H,2) - 2016*G4*G4pxxx*G4xx*G5px*pow(H,2) - 13488*G4pxx*G4x*G4xx*G5px*pow(H,2) + 
   11136*G4px*pow(G4xx,2)*G5px*pow(H,2) + 1152*G4*G4pxx*G4xxx*G5px*pow(H,2) + 1728*G4px*G4x*G4xxx*G5px*pow(H,2) - 1296*G4p*G4xx*G4xxx*G5px*pow(H,2) - 3696*G4pxxx*G4x*G5p*G5px*pow(H,2) + 
   12216*G4pxx*G4xx*G5p*G5px*pow(H,2) - 2808*G4px*G4xxx*G5p*G5px*pow(H,2) + 1176*G4pxxx*pow(G5p,2)*G5px*pow(H,2) - 240*pow(G4xx,2)*G5pp*G5px*pow(H,2) - 144*G4x*G4xxx*G5pp*G5px*pow(H,2) + 
   720*G4xxx*G5p*G5pp*G5px*pow(H,2) + 6192*G4x*G4xx*G5ppx*G5px*pow(H,2) - 720*G4*G4xxx*G5ppx*G5px*pow(H,2) - 5400*G4xx*G5p*G5ppx*G5px*pow(H,2) - 1392*pow(G4x,2)*G5ppxx*G5px*pow(H,2) + 
   816*G4*G4xx*G5ppxx*G5px*pow(H,2) + 1800*G4x*G5p*G5ppxx*G5px*pow(H,2) - 552*pow(G5p,2)*G5ppxx*G5px*pow(H,2) + 576*G4*G4pxxx*pow(G5px,2)*pow(H,2) + 3504*G4pxx*G4x*pow(G5px,2)*pow(H,2) - 
   9720*G4px*G4xx*pow(G5px,2)*pow(H,2) + 144*G4p*G4xxx*pow(G5px,2)*pow(H,2) - 3120*G4pxx*G5p*pow(G5px,2)*pow(H,2) + 456*G4xx*G5pp*pow(G5px,2)*pow(H,2) - 
   1704*G4x*G5ppx*pow(G5px,2)*pow(H,2) + 1572*G5p*G5ppx*pow(G5px,2)*pow(H,2) - 240*G4*G5ppxx*pow(G5px,2)*pow(H,2) + 2544*G4px*pow(G5px,3)*pow(H,2) - 168*G5pp*pow(G5px,3)*pow(H,2) + 
   192*G4*G4pxxx*G4x*G5pxx*pow(H,2) - 2544*G4pxx*pow(G4x,2)*G5pxx*pow(H,2) + 336*G4*G4pxx*G4xx*G5pxx*pow(H,2) - 4704*G4px*G4x*G4xx*G5pxx*pow(H,2) - 1824*G4p*pow(G4xx,2)*G5pxx*pow(H,2) + 
   720*G4*G4px*G4xxx*G5pxx*pow(H,2) - 144*G4p*G4x*G4xxx*G5pxx*pow(H,2) - 96*G4*G4pxxx*G5p*G5pxx*pow(H,2) + 2712*G4pxx*G4x*G5p*G5pxx*pow(H,2) + 3984*G4px*G4xx*G5p*G5pxx*pow(H,2) + 
   72*G4p*G4xxx*G5p*G5pxx*pow(H,2) - 720*G4pxx*pow(G5p,2)*G5pxx*pow(H,2) - 432*G4x*G4xx*G5pp*G5pxx*pow(H,2) - 144*G4*G4xxx*G5pp*G5pxx*pow(H,2) + 264*G4xx*G5p*G5pp*G5pxx*pow(H,2) + 
   1008*pow(G4x,2)*G5ppx*G5pxx*pow(H,2) - 96*G4*G4xx*G5ppx*G5pxx*pow(H,2) - 1056*G4x*G5p*G5ppx*G5pxx*pow(H,2) + 276*pow(G5p,2)*G5ppx*G5pxx*pow(H,2) - 384*G4*G4pxx*G5px*G5pxx*pow(H,2) + 
   2760*G4px*G4x*G5px*G5pxx*pow(H,2) + 1848*G4p*G4xx*G5px*G5pxx*pow(H,2) - 1776*G4px*G5p*G5px*G5pxx*pow(H,2) - 132*G5p*G5pp*G5px*G5pxx*pow(H,2) + 120*G4*G5ppx*G5px*G5pxx*pow(H,2) - 
   432*G4p*pow(G5px,2)*G5pxx*pow(H,2) - 312*G4*G4px*pow(G5pxx,2)*pow(H,2) + 24*G4p*G4x*pow(G5pxx,2)*pow(H,2) - 12*G4p*G5p*pow(G5pxx,2)*pow(H,2) + 72*G4*G5pp*pow(G5pxx,2)*pow(H,2) + 
   192*G4*G4pxx*G4x*G5pxxx*pow(H,2) + 912*G4px*pow(G4x,2)*G5pxxx*pow(H,2) - 624*G4*G4px*G4xx*G5pxxx*pow(H,2) + 240*G4p*G4x*G4xx*G5pxxx*pow(H,2) - 96*G4*G4pxx*G5p*G5pxxx*pow(H,2) - 
   1224*G4px*G4x*G5p*G5pxxx*pow(H,2) - 120*G4p*G4xx*G5p*G5pxxx*pow(H,2) + 384*G4px*pow(G5p,2)*G5pxxx*pow(H,2) - 48*pow(G4x,2)*G5pp*G5pxxx*pow(H,2) + 96*G4*G4xx*G5pp*G5pxxx*pow(H,2) + 
   96*G4x*G5p*G5pp*G5pxxx*pow(H,2) - 36*pow(G5p,2)*G5pp*G5pxxx*pow(H,2) - 96*G4*G4x*G5ppx*G5pxxx*pow(H,2) + 48*G4*G5p*G5ppx*G5pxxx*pow(H,2) + 336*G4*G4px*G5px*G5pxxx*pow(H,2) - 
   144*G4p*G4x*G5px*G5pxxx*pow(H,2) + 72*G4p*G5p*G5px*G5pxxx*pow(H,2) - 48*G4*G5pp*G5px*G5pxxx*pow(H,2) + 576*G4*G4pxx*G4pxxx*G5x*pow(H,2) - 2880*pow(G4pxx,2)*G4x*G5x*pow(H,2) + 
   6192*G4px*G4pxxx*G4x*G5x*pow(H,2) - 15696*G4px*G4pxx*G4xx*G5x*pow(H,2) + 720*G4p*G4pxxx*G4xx*G5x*pow(H,2) - 2808*G3pxx*G4x*G4xx*G5x*pow(H,2) + 5040*G4ppxx*G4x*G4xx*G5x*pow(H,2) + 
   1728*G3px*pow(G4xx,2)*G5x*pow(H,2) + 288*G4ppx*pow(G4xx,2)*G5x*pow(H,2) - 576*G4*G4ppxx*G4xxx*G5x*pow(H,2) + 1584*pow(G4px,2)*G4xxx*G5x*pow(H,2) - 288*G4p*G4pxx*G4xxx*G5x*pow(H,2) + 
   1296*G3px*G4x*G4xxx*G5x*pow(H,2) - 4032*G4ppx*G4x*G4xxx*G5x*pow(H,2) + 144*G4pp*G4xx*G4xxx*G5x*pow(H,2) + 1296*pow(G4pxx,2)*G5p*G5x*pow(H,2) - 4032*G4px*G4pxxx*G5p*G5x*pow(H,2) + 
   1836*G3pxx*G4xx*G5p*G5x*pow(H,2) - 2808*G4ppxx*G4xx*G5p*G5x*pow(H,2) - 756*G3px*G4xxx*G5p*G5x*pow(H,2) + 2520*G4ppx*G4xxx*G5p*G5x*pow(H,2) - 432*G4pxxx*G4x*G5pp*G5x*pow(H,2) + 
   216*G4pxx*G4xx*G5pp*G5x*pow(H,2) - 216*G4px*G4xxx*G5pp*G5x*pow(H,2) + 360*G4pxxx*G5p*G5pp*G5x*pow(H,2) - 72*G4xxx*pow(G5pp,2)*G5x*pow(H,2) - 1152*pow(G4xx,2)*G5ppp*G5x*pow(H,2) + 
   432*G4x*G4xxx*G5ppp*G5x*pow(H,2) - 360*G4xxx*G5p*G5ppp*G5x*pow(H,2) + 288*G4x*G4xx*G5pppx*G5x*pow(H,2) + 288*G4*G4xxx*G5pppx*G5x*pow(H,2) - 432*G4xx*G5p*G5pppx*G5x*pow(H,2) - 
   288*G4*G4pxxx*G5ppx*G5x*pow(H,2) + 1608*G4pxx*G4x*G5ppx*G5x*pow(H,2) + 7608*G4px*G4xx*G5ppx*G5x*pow(H,2) + 216*G4p*G4xxx*G5ppx*G5x*pow(H,2) - 516*G4pxx*G5p*G5ppx*G5x*pow(H,2) - 
   288*G4xx*G5pp*G5ppx*G5x*pow(H,2) - 144*G5p*pow(G5ppx,2)*G5x*pow(H,2) - 192*G4*G4pxx*G5ppxx*G5x*pow(H,2) - 3144*G4px*G4x*G5ppxx*G5x*pow(H,2) - 264*G4p*G4xx*G5ppxx*G5x*pow(H,2) + 
   1992*G4px*G5p*G5ppxx*G5x*pow(H,2) + 216*G4x*G5pp*G5ppxx*G5x*pow(H,2) - 180*G5p*G5pp*G5ppxx*G5x*pow(H,2) + 144*G4*G5ppx*G5ppxx*G5x*pow(H,2) + 8520*G4px*G4pxx*G5px*G5x*pow(H,2) - 
   432*G4p*G4pxxx*G5px*G5x*pow(H,2) + 1488*G3pxx*G4x*G5px*G5x*pow(H,2) - 2976*G4ppxx*G4x*G5px*G5x*pow(H,2) - 120*G2xx*G4xx*G5px*G5x*pow(H,2) - 1836*G3px*G4xx*G5px*G5x*pow(H,2) + 
   312*G4ppx*G4xx*G5px*G5x*pow(H,2) - 996*G3pxx*G5p*G5px*G5x*pow(H,2) + 1704*G4ppxx*G5p*G5px*G5x*pow(H,2) - 288*G4pxx*G5pp*G5px*G5x*pow(H,2) + 1008*G4xx*G5ppp*G5px*G5x*pow(H,2) + 
   144*G5p*G5pppx*G5px*G5x*pow(H,2) - 4356*G4px*G5ppx*G5px*G5x*pow(H,2) + 216*G5pp*G5ppx*G5px*G5x*pow(H,2) + 168*G4p*G5ppxx*G5px*G5x*pow(H,2) + 60*G2xx*pow(G5px,2)*G5x*pow(H,2) + 
   504*G3px*pow(G5px,2)*G5x*pow(H,2) - 264*G4ppx*pow(G5px,2)*G5x*pow(H,2) - 216*G5ppp*pow(G5px,2)*G5x*pow(H,2) + 48*G3pxx*G4*G5pxx*G5x*pow(H,2) + 192*G4*G4ppxx*G5pxx*G5x*pow(H,2) + 
   1044*pow(G4px,2)*G5pxx*G5x*pow(H,2) - 24*G4p*G4pxx*G5pxx*G5x*pow(H,2) + 120*G2xx*G4x*G5pxx*G5x*pow(H,2) - 720*G3px*G4x*G5pxx*G5x*pow(H,2) + 2016*G4ppx*G4x*G5pxx*G5x*pow(H,2) + 
   372*G2x*G4xx*G5pxx*G5x*pow(H,2) - 744*G3p*G4xx*G5pxx*G5x*pow(H,2) + 216*G4pp*G4xx*G5pxx*G5x*pow(H,2) - 78*G2xx*G5p*G5pxx*G5x*pow(H,2) + 336*G3px*G5p*G5pxx*G5x*pow(H,2) - 
   1068*G4ppx*G5p*G5pxx*G5x*pow(H,2) - 60*G4px*G5pp*G5pxx*G5x*pow(H,2) + 36*pow(G5pp,2)*G5pxx*G5x*pow(H,2) - 216*G4x*G5ppp*G5pxx*G5x*pow(H,2) + 180*G5p*G5ppp*G5pxx*G5x*pow(H,2) - 
   144*G4*G5pppx*G5pxx*G5x*pow(H,2) + 12*G4p*G5ppx*G5pxx*G5x*pow(H,2) - 210*G2x*G5px*G5pxx*G5x*pow(H,2) + 420*G3p*G5px*G5pxx*G5x*pow(H,2) - 180*G4pp*G5px*G5pxx*G5x*pow(H,2) + 
   48*G3px*G4*G5pxxx*G5x*pow(H,2) - 96*G4*G4ppx*G5pxxx*G5x*pow(H,2) - 144*G4p*G4px*G5pxxx*G5x*pow(H,2) - 108*G2x*G4x*G5pxxx*G5x*pow(H,2) + 216*G3p*G4x*G5pxxx*G5x*pow(H,2) - 
   120*G4pp*G4x*G5pxxx*G5x*pow(H,2) + 66*G2x*G5p*G5pxxx*G5x*pow(H,2) - 132*G3p*G5p*G5pxxx*G5x*pow(H,2) + 84*G4pp*G5p*G5pxxx*G5x*pow(H,2) + 24*G4p*G5pp*G5pxxx*G5x*pow(H,2) + 
   864*G3pxx*G4px*pow(G5x,2)*pow(H,2) - 1728*G4ppxx*G4px*pow(G5x,2)*pow(H,2) + 144*G2xx*G4pxx*pow(G5x,2)*pow(H,2) - 252*G3px*G4pxx*pow(G5x,2)*pow(H,2) + 
   792*G4ppx*G4pxx*pow(G5x,2)*pow(H,2) - 180*G2x*G4pxxx*pow(G5x,2)*pow(H,2) + 360*G3p*G4pxxx*pow(G5x,2)*pow(H,2) - 216*G4pp*G4pxxx*pow(G5x,2)*pow(H,2) - 
   252*G2pxx*G4xx*pow(G5x,2)*pow(H,2) + 108*G3ppx*G4xx*pow(G5x,2)*pow(H,2) + 288*G4pppx*G4xx*pow(G5x,2)*pow(H,2) + 72*G2px*G4xxx*pow(G5x,2)*pow(H,2) - 
   252*G3pp*G4xxx*pow(G5x,2)*pow(H,2) + 216*G4ppp*G4xxx*pow(G5x,2)*pow(H,2) - 72*G3pxx*G5pp*pow(G5x,2)*pow(H,2) + 144*G4ppxx*G5pp*pow(G5x,2)*pow(H,2) - 
   144*G4pxx*G5ppp*pow(G5x,2)*pow(H,2) - 82*G2xx*G5ppx*pow(G5x,2)*pow(H,2) + 10*G3px*G5ppx*pow(G5x,2)*pow(H,2) + 72*G4ppx*G5ppx*pow(G5x,2)*pow(H,2) + 
   94*G2x*G5ppxx*pow(G5x,2)*pow(H,2) - 188*G3p*G5ppxx*pow(G5x,2)*pow(H,2) + 108*G4pp*G5ppxx*pow(G5x,2)*pow(H,2) + 136*G2pxx*G5px*pow(G5x,2)*pow(H,2) - 
   100*G3ppx*G5px*pow(G5x,2)*pow(H,2) - 72*G4pppx*G5px*pow(G5x,2)*pow(H,2) - 32*G2px*G5pxx*pow(G5x,2)*pow(H,2) + 126*G3pp*G5pxx*pow(G5x,2)*pow(H,2) - 
   108*G4ppp*G5pxx*pow(G5x,2)*pow(H,2) + 4*G2p*G5pxxx*pow(G5x,2)*pow(H,2) - 288*G4*pow(G4pxx,2)*G5xx*pow(H,2) - 864*G4*G4px*G4pxxx*G5xx*pow(H,2) - 5664*G4px*G4pxx*G4x*G5xx*pow(H,2) + 
   288*G4p*G4pxxx*G4x*G5xx*pow(H,2) - 504*G3pxx*pow(G4x,2)*G5xx*pow(H,2) + 1104*G4ppxx*pow(G4x,2)*G5xx*pow(H,2) + 360*G3pxx*G4*G4xx*G5xx*pow(H,2) - 144*G4*G4ppxx*G4xx*G5xx*pow(H,2) + 
   6240*pow(G4px,2)*G4xx*G5xx*pow(H,2) - 2160*G4p*G4pxx*G4xx*G5xx*pow(H,2) + 1080*G3px*G4x*G4xx*G5xx*pow(H,2) - 480*G4ppx*G4x*G4xx*G5xx*pow(H,2) - 2880*G4pp*pow(G4xx,2)*G5xx*pow(H,2) - 
   72*G3px*G4*G4xxx*G5xx*pow(H,2) + 432*G4*G4ppx*G4xxx*G5xx*pow(H,2) - 576*G4p*G4px*G4xxx*G5xx*pow(H,2) + 144*G4pp*G4x*G4xxx*G5xx*pow(H,2) + 4920*G4px*G4pxx*G5p*G5xx*pow(H,2) - 
   144*G4p*G4pxxx*G5p*G5xx*pow(H,2) + 684*G3pxx*G4x*G5p*G5xx*pow(H,2) - 1176*G4ppxx*G4x*G5p*G5xx*pow(H,2) - 1116*G3px*G4xx*G5p*G5xx*pow(H,2) - 336*G4ppx*G4xx*G5p*G5xx*pow(H,2) - 
   72*G4pp*G4xxx*G5p*G5xx*pow(H,2) - 216*G3pxx*pow(G5p,2)*G5xx*pow(H,2) + 312*G4ppxx*pow(G5p,2)*G5xx*pow(H,2) + 144*G4*G4pxxx*G5pp*G5xx*pow(H,2) - 192*G4pxx*G4x*G5pp*G5xx*pow(H,2) - 
   1344*G4px*G4xx*G5pp*G5xx*pow(H,2) - 156*G4pxx*G5p*G5pp*G5xx*pow(H,2) - 96*G4xx*pow(G5pp,2)*G5xx*pow(H,2) - 192*G4x*G4xx*G5ppp*G5xx*pow(H,2) - 144*G4*G4xxx*G5ppp*G5xx*pow(H,2) + 
   960*G4xx*G5p*G5ppp*G5xx*pow(H,2) - 48*pow(G4x,2)*G5pppx*G5xx*pow(H,2) - 288*G4*G4xx*G5pppx*G5xx*pow(H,2) - 96*G4x*G5p*G5pppx*G5xx*pow(H,2) + 60*pow(G5p,2)*G5pppx*G5xx*pow(H,2) + 
   360*G4*G4pxx*G5ppx*G5xx*pow(H,2) + 2880*G4px*G4x*G5ppx*G5xx*pow(H,2) + 288*G4p*G4xx*G5ppx*G5xx*pow(H,2) - 2448*G4px*G5p*G5ppx*G5xx*pow(H,2) + 24*G4x*G5pp*G5ppx*G5xx*pow(H,2) + 
   204*G5p*G5pp*G5ppx*G5xx*pow(H,2) - 144*G4*pow(G5ppx,2)*G5xx*pow(H,2) + 384*G4*G4px*G5ppxx*G5xx*pow(H,2) - 96*G4p*G4x*G5ppxx*G5xx*pow(H,2) + 48*G4p*G5p*G5ppxx*G5xx*pow(H,2) - 
   72*G4*G5pp*G5ppxx*G5xx*pow(H,2) - 216*G3pxx*G4*G5px*G5xx*pow(H,2) + 144*G4*G4ppxx*G5px*G5xx*pow(H,2) - 4932*pow(G4px,2)*G5px*G5xx*pow(H,2) + 1296*G4p*G4pxx*G5px*G5xx*pow(H,2) - 
   48*G2xx*G4x*G5px*G5xx*pow(H,2) - 492*G3px*G4x*G5px*G5xx*pow(H,2) + 312*G4ppx*G4x*G5px*G5xx*pow(H,2) - 300*G2x*G4xx*G5px*G5xx*pow(H,2) + 600*G3p*G4xx*G5px*G5xx*pow(H,2) + 
   2832*G4pp*G4xx*G5px*G5xx*pow(H,2) + 54*G2xx*G5p*G5px*G5xx*pow(H,2) + 612*G3px*G5p*G5px*G5xx*pow(H,2) - 84*G4ppx*G5p*G5px*G5xx*pow(H,2) + 1056*G4px*G5pp*G5px*G5xx*pow(H,2) + 
   12*pow(G5pp,2)*G5px*G5xx*pow(H,2) + 24*G4x*G5ppp*G5px*G5xx*pow(H,2) - 444*G5p*G5ppp*G5px*G5xx*pow(H,2) + 144*G4*G5pppx*G5px*G5xx*pow(H,2) - 288*G4p*G5ppx*G5px*G5xx*pow(H,2) + 
   174*G2x*pow(G5px,2)*G5xx*pow(H,2) - 348*G3p*pow(G5px,2)*G5xx*pow(H,2) - 660*G4pp*pow(G5px,2)*G5xx*pow(H,2) - 12*G2xx*G4*G5pxx*G5xx*pow(H,2) - 48*G3px*G4*G5pxx*G5xx*pow(H,2) - 
   24*G4*G4ppx*G5pxx*G5xx*pow(H,2) + 792*G4p*G4px*G5pxx*G5xx*pow(H,2) + 156*G2x*G4x*G5pxx*G5xx*pow(H,2) - 312*G3p*G4x*G5pxx*G5xx*pow(H,2) + 24*G4pp*G4x*G5pxx*G5xx*pow(H,2) - 
   114*G2x*G5p*G5pxx*G5xx*pow(H,2) + 228*G3p*G5p*G5pxx*G5xx*pow(H,2) - 84*G4pp*G5p*G5pxx*G5xx*pow(H,2) - 96*G4p*G5pp*G5pxx*G5xx*pow(H,2) + 72*G4*G5ppp*G5pxx*G5xx*pow(H,2) + 
   12*G2x*G4*G5pxxx*G5xx*pow(H,2) - 24*G3p*G4*G5pxxx*G5xx*pow(H,2) + 24*G4*G4pp*G5pxxx*G5xx*pow(H,2) + 72*G3pxx*G4p*G5x*G5xx*pow(H,2) - 84*G2xx*G4px*G5x*G5xx*pow(H,2) - 
   912*G3px*G4px*G5x*G5xx*pow(H,2) + 1032*G4ppx*G4px*G5x*G5xx*pow(H,2) + 408*G2x*G4pxx*G5x*G5xx*pow(H,2) - 816*G3p*G4pxx*G5x*G5xx*pow(H,2) + 264*G4pp*G4pxx*G5x*G5xx*pow(H,2) - 
   180*G2pxx*G4x*G5x*G5xx*pow(H,2) + 156*G3ppx*G4x*G5x*G5xx*pow(H,2) + 48*G4pppx*G4x*G5x*G5xx*pow(H,2) - 204*G2px*G4xx*G5x*G5xx*pow(H,2) + 744*G3pp*G4xx*G5x*G5xx*pow(H,2) - 
   672*G4ppp*G4xx*G5x*G5xx*pow(H,2) + 36*G2p*G4xxx*G5x*G5xx*pow(H,2) + 120*G2pxx*G5p*G5x*G5xx*pow(H,2) - 36*G3ppx*G5p*G5x*G5xx*pow(H,2) - 168*G4pppx*G5p*G5x*G5xx*pow(H,2) + 
   30*G2xx*G5pp*G5x*G5xx*pow(H,2) + 18*G3px*G5pp*G5x*G5xx*pow(H,2) - 24*G4ppx*G5pp*G5x*G5xx*pow(H,2) + 192*G4px*G5ppp*G5x*G5xx*pow(H,2) - 72*G4p*G5pppx*G5x*G5xx*pow(H,2) - 
   192*G2x*G5ppx*G5x*G5xx*pow(H,2) + 384*G3p*G5ppx*G5x*G5xx*pow(H,2) - 132*G4pp*G5ppx*G5x*G5xx*pow(H,2) + 108*G2px*G5px*G5x*G5xx*pow(H,2) - 366*G3pp*G5px*G5x*G5xx*pow(H,2) + 
   300*G4ppp*G5px*G5x*G5xx*pow(H,2) - 30*G2p*G5pxx*G5x*G5xx*pow(H,2) - 8*G2ppx*pow(G5x,2)*G5xx*pow(H,2) + 8*G3ppp*pow(G5x,2)*G5xx*pow(H,2) + 12*G2pxx*G4*pow(G5xx,2)*pow(H,2) + 
   24*G3ppx*G4*pow(G5xx,2)*pow(H,2) + 12*G2xx*G4p*pow(G5xx,2)*pow(H,2) - 30*G3px*G4p*pow(G5xx,2)*pow(H,2) - 72*G4*G4pppx*pow(G5xx,2)*pow(H,2) - 180*G4p*G4ppx*pow(G5xx,2)*pow(H,2) - 
   96*G2x*G4px*pow(G5xx,2)*pow(H,2) + 192*G3p*G4px*pow(G5xx,2)*pow(H,2) + 372*G4pp*G4px*pow(G5xx,2)*pow(H,2) - 48*G2px*G4x*pow(G5xx,2)*pow(H,2) + 120*G3pp*G4x*pow(G5xx,2)*pow(H,2) - 
   48*G4ppp*G4x*pow(G5xx,2)*pow(H,2) - 108*G2p*G4xx*pow(G5xx,2)*pow(H,2) + 24*G2px*G5p*pow(G5xx,2)*pow(H,2) - 114*G3pp*G5p*pow(G5xx,2)*pow(H,2) + 132*G4ppp*G5p*pow(G5xx,2)*pow(H,2) + 
   30*G2x*G5pp*pow(G5xx,2)*pow(H,2) - 60*G3p*G5pp*pow(G5xx,2)*pow(H,2) - 24*G4pp*G5pp*pow(G5xx,2)*pow(H,2) + 108*G4p*G5ppp*pow(G5xx,2)*pow(H,2) + 60*G2p*G5px*pow(G5xx,2)*pow(H,2) + 
   6*G2pp*G5x*pow(G5xx,2)*pow(H,2) + 480*G4*G4px*G4pxx*G5xxx*pow(H,2) - 192*G4*G4ppxx*G4x*G5xxx*pow(H,2) + 288*pow(G4px,2)*G4x*G5xxx*pow(H,2) - 96*G4p*G4pxx*G4x*G5xxx*pow(H,2) + 
   216*G3px*pow(G4x,2)*G5xxx*pow(H,2) - 624*G4ppx*pow(G4x,2)*G5xxx*pow(H,2) - 72*G3px*G4*G4xx*G5xxx*pow(H,2) + 336*G4*G4ppx*G4xx*G5xxx*pow(H,2) - 384*G4p*G4px*G4xx*G5xxx*pow(H,2) + 
   48*G4pp*G4x*G4xx*G5xxx*pow(H,2) + 96*G4*G4ppxx*G5p*G5xxx*pow(H,2) - 456*pow(G4px,2)*G5p*G5xxx*pow(H,2) + 48*G4p*G4pxx*G5p*G5xxx*pow(H,2) - 252*G3px*G4x*G5p*G5xxx*pow(H,2) + 
   792*G4ppx*G4x*G5p*G5xxx*pow(H,2) - 24*G4pp*G4xx*G5p*G5xxx*pow(H,2) + 72*G3px*pow(G5p,2)*G5xxx*pow(H,2) - 240*G4ppx*pow(G5p,2)*G5xxx*pow(H,2) - 96*G4*G4pxx*G5pp*G5xxx*pow(H,2) + 
   96*G4px*G4x*G5pp*G5xxx*pow(H,2) - 24*G4p*G4xx*G5pp*G5xxx*pow(H,2) + 156*G4px*G5p*G5pp*G5xxx*pow(H,2) - 48*G4x*pow(G5pp,2)*G5xxx*pow(H,2) + 48*pow(G4x,2)*G5ppp*G5xxx*pow(H,2) - 
   96*G4*G4xx*G5ppp*G5xxx*pow(H,2) - 96*G4x*G5p*G5ppp*G5xxx*pow(H,2) + 36*pow(G5p,2)*G5ppp*G5xxx*pow(H,2) + 96*G4*G4x*G5pppx*G5xxx*pow(H,2) - 48*G4*G5p*G5pppx*G5xxx*pow(H,2) - 
   264*G4*G4px*G5ppx*G5xxx*pow(H,2) + 72*G4p*G4x*G5ppx*G5xxx*pow(H,2) - 36*G4p*G5p*G5ppx*G5xxx*pow(H,2) + 48*G4*G5pp*G5ppx*G5xxx*pow(H,2) + 48*G3px*G4*G5px*G5xxx*pow(H,2) - 
   192*G4*G4ppx*G5px*G5xxx*pow(H,2) + 48*G4p*G4px*G5px*G5xxx*pow(H,2) + 12*G2x*G4x*G5px*G5xxx*pow(H,2) - 24*G3p*G4x*G5px*G5xxx*pow(H,2) + 24*G4pp*G4x*G5px*G5xxx*pow(H,2) + 
   6*G2x*G5p*G5px*G5xxx*pow(H,2) - 12*G3p*G5p*G5px*G5xxx*pow(H,2) + 12*G4pp*G5p*G5px*G5xxx*pow(H,2) + 48*G4p*G5pp*G5px*G5xxx*pow(H,2) + 48*G4*G5ppp*G5px*G5xxx*pow(H,2) - 
   12*G2x*G4*G5pxx*G5xxx*pow(H,2) + 24*G3p*G4*G5pxx*G5xxx*pow(H,2) - 24*G4*G4pp*G5pxx*G5xxx*pow(H,2) - 48*G3ppx*G4*G5x*G5xxx*pow(H,2) - 12*G3px*G4p*G5x*G5xxx*pow(H,2) + 
   96*G4*G4pppx*G5x*G5xxx*pow(H,2) + 72*G4p*G4ppx*G5x*G5xxx*pow(H,2) - 12*G2x*G4px*G5x*G5xxx*pow(H,2) + 24*G3p*G4px*G5x*G5xxx*pow(H,2) + 48*G4pp*G4px*G5x*G5xxx*pow(H,2) + 
   48*G2px*G4x*G5x*G5xxx*pow(H,2) - 156*G3pp*G4x*G5x*G5xxx*pow(H,2) + 120*G4ppp*G4x*G5x*G5xxx*pow(H,2) + 24*G2p*G4xx*G5x*G5xxx*pow(H,2) - 24*G2px*G5p*G5x*G5xxx*pow(H,2) + 
   90*G3pp*G5p*G5x*G5xxx*pow(H,2) - 84*G4ppp*G5p*G5x*G5xxx*pow(H,2) - 6*G2x*G5pp*G5x*G5xxx*pow(H,2) + 12*G3p*G5pp*G5x*G5xxx*pow(H,2) - 36*G4pp*G5pp*G5x*G5xxx*pow(H,2) - 
   24*G4p*G5ppp*G5x*G5xxx*pow(H,2) - 12*G2p*G5px*G5x*G5xxx*pow(H,2) - 4*G2pp*pow(G5x,2)*G5xxx*pow(H,2) + 12*G3pp*G4*G5xx*G5xxx*pow(H,2) + 12*G2x*G4p*G5xx*G5xxx*pow(H,2) - 
   24*G3p*G4p*G5xx*G5xxx*pow(H,2) + 24*G4p*G4pp*G5xx*G5xxx*pow(H,2) - 24*G4*G4ppp*G5xx*G5xxx*pow(H,2) + 12*G2p*G4x*G5xx*G5xxx*pow(H,2) - 6*G2p*G5p*G5xx*G5xxx*pow(H,2) + 
   1440*pow(G4xx,2)*G5px*G5x*pow(H,4) + 360*G4x*G4xxx*G5px*G5x*pow(H,4) + 180*G4xxx*G5p*G5px*G5x*pow(H,4) - 1440*G4xx*pow(G5px,2)*G5x*pow(H,4) + 360*pow(G5px,3)*G5x*pow(H,4) + 
   2880*G4x*G4xx*G5pxx*G5x*pow(H,4) - 504*G4*G4xxx*G5pxx*G5x*pow(H,4) - 3276*G4xx*G5p*G5pxx*G5x*pow(H,4) - 1748*G4x*G5px*G5pxx*G5x*pow(H,4) + 1660*G5p*G5px*G5pxx*G5x*pow(H,4) + 
   292*G4*pow(G5pxx,2)*G5x*pow(H,4) - 936*pow(G4x,2)*G5pxxx*G5x*pow(H,4) + 648*G4*G4xx*G5pxxx*G5x*pow(H,4) + 1404*G4x*G5p*G5pxxx*G5x*pow(H,4) - 504*pow(G5p,2)*G5pxxx*G5x*pow(H,4) - 
   328*G4*G5px*G5pxxx*G5x*pow(H,4) - 3024*G4pxxx*G4x*pow(G5x,2)*pow(H,4) + 4536*G4pxx*G4xx*pow(G5x,2)*pow(H,4) - 432*G4px*G4xxx*pow(G5x,2)*pow(H,4) + 
   2268*G4pxxx*G5p*pow(G5x,2)*pow(H,4) - 216*G4xxx*G5pp*pow(G5x,2)*pow(H,4) - 1212*G4xx*G5ppx*pow(G5x,2)*pow(H,4) + 1380*G4x*G5ppxx*pow(G5x,2)*pow(H,4) - 
   1008*G5p*G5ppxx*pow(G5x,2)*pow(H,4) - 2460*G4pxx*G5px*pow(G5x,2)*pow(H,4) + 702*G5ppx*G5px*pow(G5x,2)*pow(H,4) - 438*G4px*G5pxx*pow(G5x,2)*pow(H,4) + 
   120*G4p*G5pxxx*pow(G5x,2)*pow(H,4) - 270*G3pxx*pow(G5x,3)*pow(H,4) + 360*G4ppxx*pow(G5x,3)*pow(H,4) + 90*G5pppx*pow(G5x,3)*pow(H,4) + 768*G4x*G4xx*G5px*G5xx*pow(H,4) - 
   72*G4*G4xxx*G5px*G5xx*pow(H,4) + 84*G4xx*G5p*G5px*G5xx*pow(H,4) - 384*G4x*pow(G5px,2)*G5xx*pow(H,4) - 78*G5p*pow(G5px,2)*G5xx*pow(H,4) + 744*pow(G4x,2)*G5pxx*G5xx*pow(H,4) - 
   384*G4*G4xx*G5pxx*G5xx*pow(H,4) - 1404*G4x*G5p*G5pxx*G5xx*pow(H,4) + 660*pow(G5p,2)*G5pxx*G5xx*pow(H,4) + 144*G4*G5px*G5pxx*G5xx*pow(H,4) + 216*G4*G4x*G5pxxx*G5xx*pow(H,4) - 
   168*G4*G5p*G5pxxx*G5xx*pow(H,4) + 720*G4*G4pxxx*G5x*G5xx*pow(H,4) + 4176*G4pxx*G4x*G5x*G5xx*pow(H,4) - 1152*G4px*G4xx*G5x*G5xx*pow(H,4) + 648*G4p*G4xxx*G5x*G5xx*pow(H,4) - 
   3744*G4pxx*G5p*G5x*G5xx*pow(H,4) - 1836*G4xx*G5pp*G5x*G5xx*pow(H,4) - 1444*G4x*G5ppx*G5x*G5xx*pow(H,4) + 1286*G5p*G5ppx*G5x*G5xx*pow(H,4) - 400*G4*G5ppxx*G5x*G5xx*pow(H,4) + 
   96*G4px*G5px*G5x*G5xx*pow(H,4) + 1170*G5pp*G5px*G5x*G5xx*pow(H,4) - 636*G4p*G5pxx*G5x*G5xx*pow(H,4) + 261*G3px*pow(G5x,2)*G5xx*pow(H,4) + 90*G4ppx*pow(G5x,2)*G5xx*pow(H,4) - 
   138*G5ppp*pow(G5x,2)*G5xx*pow(H,4) - 120*G4*G4pxx*pow(G5xx,2)*pow(H,4) - 432*G4px*G4x*pow(G5xx,2)*pow(H,4) + 984*G4p*G4xx*pow(G5xx,2)*pow(H,4) + 348*G4px*G5p*pow(G5xx,2)*pow(H,4) - 
   204*G4x*G5pp*pow(G5xx,2)*pow(H,4) + 162*G5p*G5pp*pow(G5xx,2)*pow(H,4) + 132*G4*G5ppx*pow(G5xx,2)*pow(H,4) - 546*G4p*G5px*pow(G5xx,2)*pow(H,4) - 402*G4pp*G5x*pow(G5xx,2)*pow(H,4) + 
   120*pow(G4x,2)*G5px*G5xxx*pow(H,4) - 120*G4*G4xx*G5px*G5xxx*pow(H,4) - 36*G4x*G5p*G5px*G5xxx*pow(H,4) - 48*pow(G5p,2)*G5px*G5xxx*pow(H,4) + 96*G4*pow(G5px,2)*G5xxx*pow(H,4) - 
   168*G4*G4x*G5pxx*G5xxx*pow(H,4) + 120*G4*G5p*G5pxx*G5xxx*pow(H,4) - 432*G4*G4pxx*G5x*G5xxx*pow(H,4) - 144*G4px*G4x*G5x*G5xxx*pow(H,4) + 504*G4p*G4xx*G5x*G5xxx*pow(H,4) + 
   252*G4px*G5p*G5x*G5xxx*pow(H,4) - 216*G4x*G5pp*G5x*G5xxx*pow(H,4) + 54*G5p*G5pp*G5x*G5xxx*pow(H,4) + 220*G4*G5ppx*G5x*G5xxx*pow(H,4) - 228*G4p*G5px*G5x*G5xxx*pow(H,4) - 
   12*G4pp*pow(G5x,2)*G5xxx*pow(H,4) + 24*G4*G4px*G5xx*G5xxx*pow(H,4) + 192*G4p*G4x*G5xx*G5xxx*pow(H,4) - 156*G4p*G5p*G5xx*G5xxx*pow(H,4) + 12*G4*G5pp*G5xx*G5xxx*pow(H,4) - 
   45*G5pxx*pow(G5x,3)*pow(H,6) + 585*G5px*pow(G5x,2)*G5xx*pow(H,6) - 3*pow(G3x,2)*
    (72*G4pxxx*G4x - 360*G4pxx*G4xx + 72*G4px*G4xxx - 36*G4pxxx*G5p - 36*G4xxx*G5pp + 144*G4xx*G5ppx - 32*G4x*G5ppxx + 16*G5p*G5ppxx + 18*G3xx*G5px + 168*G4pxx*G5px - 84*G5ppx*G5px + 34*G4px*G5pxx + 
      10*G5pp*G5pxx + 18*G3pxx*G5x - 24*G4ppxx*G5x - 6*G5pppx*G5x - 21*G3px*G5xx + 6*G4ppx*G5xx + 18*G5ppp*G5xx - 87*G5pxx*G5x*pow(H,2) + 51*G5px*G5xx*pow(H,2)) + 
   6*G3xx*(36*G3px*G4x*G4xx - 168*G4ppx*G4x*G4xx - 18*G3px*G4xx*G5p + 84*G4ppx*G4xx*G5p + 12*G4ppxx*pow(-2*G4x + G5p,2) + 48*G4pxx*G4x*G5pp - 24*G4pxx*G5p*G5pp + 24*G4xx*pow(G5pp,2) + 
      48*G4x*G4xx*G5ppp - 24*G4xx*G5p*G5ppp - 24*pow(G4x,2)*G5pppx + 24*G4x*G5p*G5pppx - 6*pow(G5p,2)*G5pppx - 24*G4x*G5pp*G5ppx + 12*G5p*G5pp*G5ppx + 12*pow(G4px,2)*(26*G4xx - 19*G5px) - 
      24*G3px*G4x*G5px + 96*G4ppx*G4x*G5px - 12*G2x*G4xx*G5px + 24*G3p*G4xx*G5px - 24*G4pp*G4xx*G5px + 12*G3px*G5p*G5px - 48*G4ppx*G5p*G5px - 12*pow(G5pp,2)*G5px - 24*G4x*G5ppp*G5px + 
      12*G5p*G5ppp*G5px + 6*G2x*pow(G5px,2) - 12*G3p*pow(G5px,2) + 12*G4pp*pow(G5px,2) + 6*G2x*G4x*G5pxx - 12*G3p*G4x*G5pxx + 12*G4pp*G4x*G5pxx - 3*G2x*G5p*G5pxx + 6*G3p*G5p*G5pxx - 
      6*G4pp*G5p*G5pxx + 12*G2x*G4pxx*G5x - 24*G3p*G4pxx*G5x + 24*G4pp*G4pxx*G5x + 24*G3ppx*G4x*G5x - 48*G4pppx*G4x*G5x - 12*G3pp*G4xx*G5x + 24*G4ppp*G4xx*G5x - 12*G3ppx*G5p*G5x + 24*G4pppx*G5p*G5x + 
      6*G3px*G5pp*G5x - 12*G4ppx*G5pp*G5x - 6*G2x*G5ppx*G5x + 12*G3p*G5ppx*G5x - 12*G4pp*G5ppx*G5x + 6*G3pp*G5px*G5x - 12*G4ppp*G5px*G5x + 2*G2ppx*pow(G5x,2) - 2*G3ppp*pow(G5x,2) - 6*G3pp*G4x*G5xx + 
      12*G4ppp*G4x*G5xx + 3*G3pp*G5p*G5xx - 6*G4ppp*G5p*G5xx + 3*G2x*G5pp*G5xx - 6*G3p*G5pp*G5xx + 6*G4pp*G5pp*G5xx - 84*G4x*G4xx*G5px*pow(H,2) + 114*G4xx*G5p*G5px*pow(H,2) + 
      36*G4x*pow(G5px,2)*pow(H,2) - 60*G5p*pow(G5px,2)*pow(H,2) + 60*pow(G4x,2)*G5pxx*pow(H,2) - 36*G4*G4xx*G5pxx*pow(H,2) - 78*G4x*G5p*G5pxx*pow(H,2) + 24*pow(G5p,2)*G5pxx*pow(H,2) + 
      24*G4*G5px*G5pxx*pow(H,2) + 288*G4pxx*G4x*G5x*pow(H,2) - 180*G4pxx*G5p*G5x*pow(H,2) + 54*G4xx*G5pp*G5x*pow(H,2) - 158*G4x*G5ppx*G5x*pow(H,2) + 103*G5p*G5ppx*G5x*pow(H,2) - 
      8*G4*G5ppxx*G5x*pow(H,2) - 24*G5pp*G5px*G5x*pow(H,2) - 6*G4p*G5pxx*G5x*pow(H,2) + 27*G3px*pow(G5x,2)*pow(H,2) - 90*G4ppx*pow(G5x,2)*pow(H,2) + 12*G5ppp*pow(G5x,2)*pow(H,2) - 
      24*G4*G4pxx*G5xx*pow(H,2) + 60*G4p*G4xx*G5xx*pow(H,2) + 24*G4x*G5pp*G5xx*pow(H,2) - 27*G5p*G5pp*G5xx*pow(H,2) + 18*G4*G5ppx*G5xx*pow(H,2) - 24*G4p*G5px*G5xx*pow(H,2) + 
      6*G4pp*G5x*G5xx*pow(H,2) - 6*G4px*(40*G4pxx*G4x - 20*G4pxx*G5p + 34*G4xx*G5pp - 22*G4x*G5ppx + 11*G5p*G5ppx - 20*G5pp*G5px + 5*G3px*G5x - 14*G4ppx*G5x + 2*G5ppp*G5x + G2x*G5xx - 2*G3p*G5xx + 
         2*G4pp*G5xx + 30*G4xx*G5x*pow(H,2) - 21*G5px*G5x*pow(H,2) + 12*G4x*G5xx*pow(H,2) - 11*G5p*G5xx*pow(H,2))) + 
   3*G3x*(384*G4px*G4pxxx*G4x - 216*G3pxx*G4x*G4xx + 240*G4ppxx*G4x*G4xx + 288*G3px*pow(G4xx,2) + 384*pow(G4px,2)*G4xxx + 72*G3px*G4x*G4xxx - 240*G4ppx*G4x*G4xxx - 192*G4px*G4pxxx*G5p + 
      108*G3pxx*G4xx*G5p - 120*G4ppxx*G4xx*G5p - 36*G3px*G4xxx*G5p + 120*G4ppx*G4xxx*G5p + 48*pow(G4pxx,2)*(-2*G4x + G5p) - 48*G4pxxx*G4x*G5pp - 240*G4px*G4xxx*G5pp + 24*G4pxxx*G5p*G5pp + 
      24*G4xxx*pow(G5pp,2) - 288*pow(G4xx,2)*G5ppp + 48*G4x*G4xxx*G5ppp - 24*G4xxx*G5p*G5ppp + 96*G4x*G4xx*G5pppx - 48*G4xx*G5p*G5pppx + 864*G4px*G4xx*G5ppx - 144*G4xx*G5pp*G5ppx + 
      48*G4x*pow(G5ppx,2) - 24*G5p*pow(G5ppx,2) - 176*G4px*G4x*G5ppxx + 88*G4px*G5p*G5ppxx + 24*G4x*G5pp*G5ppxx - 12*G5p*G5pp*G5ppxx + 120*G3pxx*G4x*G5px - 144*G4ppxx*G4x*G5px - 36*G2xx*G4xx*G5px - 
      324*G3px*G4xx*G5px + 144*G4ppx*G4xx*G5px - 12*G2x*G4xxx*G5px + 24*G3p*G4xxx*G5px - 24*G4pp*G4xxx*G5px - 60*G3pxx*G5p*G5px + 72*G4ppxx*G5p*G5px + 288*G4xx*G5ppp*G5px - 48*G4x*G5pppx*G5px + 
      24*G5p*G5pppx*G5px - 480*G4px*G5ppx*G5px + 72*G5pp*G5ppx*G5px + 20*G2xx*pow(G5px,2) + 88*G3px*pow(G5px,2) - 72*G4ppx*pow(G5px,2) - 72*G5ppp*pow(G5px,2) - 20*pow(G4px,2)*G5pxx + 
      12*G2xx*G4x*G5pxx - 16*G3px*G4x*G5pxx + 56*G4ppx*G4x*G5pxx + 24*G2x*G4xx*G5pxx - 48*G3p*G4xx*G5pxx + 48*G4pp*G4xx*G5pxx - 6*G2xx*G5p*G5pxx + 8*G3px*G5p*G5pxx - 28*G4ppx*G5p*G5pxx + 
      88*G4px*G5pp*G5pxx - 12*pow(G5pp,2)*G5pxx - 24*G4x*G5ppp*G5pxx + 12*G5p*G5ppp*G5pxx - 6*G2x*G5px*G5pxx + 12*G3p*G5px*G5pxx - 12*G4pp*G5px*G5pxx - 4*G2x*G4x*G5pxxx + 8*G3p*G4x*G5pxxx - 
      8*G4pp*G4x*G5pxxx + 2*G2x*G5p*G5pxxx - 4*G3p*G5p*G5pxxx + 4*G4pp*G5p*G5pxxx + 96*G3pxx*G4px*G5x - 144*G4ppxx*G4px*G5x - 12*G2x*G4pxxx*G5x + 24*G3p*G4pxxx*G5x - 24*G4pp*G4pxxx*G5x - 
      36*G2pxx*G4xx*G5x - 12*G3ppx*G4xx*G5x + 96*G4pppx*G4xx*G5x - 12*G3pp*G4xxx*G5x + 24*G4ppp*G4xxx*G5x - 12*G3pxx*G5pp*G5x + 24*G4ppxx*G5pp*G5x - 24*G4px*G5pppx*G5x - 14*G2xx*G5ppx*G5x - 
      22*G3px*G5ppx*G5x + 72*G4ppx*G5ppx*G5x + 6*G2x*G5ppxx*G5x - 12*G3p*G5ppxx*G5x + 12*G4pp*G5ppxx*G5x + 20*G2pxx*G5px*G5x + 4*G3ppx*G5px*G5x - 48*G4pppx*G5px*G5x + 4*G2px*G5pxx*G5x + 
      2*G3pp*G5pxx*G5x - 12*G4ppp*G5pxx*G5x - 16*G2xx*G4px*G5xx - 116*G3px*G4px*G5xx + 120*G4ppx*G4px*G5xx - 16*G2pxx*G4x*G5xx - 8*G3ppx*G4x*G5xx + 48*G4pppx*G4x*G5xx + 72*G3pp*G4xx*G5xx - 
      144*G4ppp*G4xx*G5xx + 8*G2pxx*G5p*G5xx + 4*G3ppx*G5p*G5xx - 24*G4pppx*G5p*G5xx + 8*G2xx*G5pp*G5xx + 16*G3px*G5pp*G5xx - 48*G4ppx*G5pp*G5xx + 72*G4px*G5ppp*G5xx - 12*G2x*G5ppx*G5xx + 
      24*G3p*G5ppx*G5xx - 24*G4pp*G5ppx*G5xx - 4*G2px*G5px*G5xx - 32*G3pp*G5px*G5xx + 72*G4ppp*G5px*G5xx - 4*G2ppx*G5x*G5xx + 4*G3ppp*G5x*G5xx - 4*G2x*G4px*G5xxx + 8*G3p*G4px*G5xxx - 8*G4pp*G4px*G5xxx - 
      4*G3pp*G4x*G5xxx + 8*G4ppp*G4x*G5xxx + 2*G3pp*G5p*G5xxx - 4*G4ppp*G5p*G5xxx + 2*G2x*G5pp*G5xxx - 4*G3p*G5pp*G5xxx + 4*G4pp*G5pp*G5xxx - 576*pow(G4xx,2)*G5px*pow(H,2) - 
      48*G4x*G4xxx*G5px*pow(H,2) + 132*G4xxx*G5p*G5px*pow(H,2) + 676*G4xx*pow(G5px,2)*pow(H,2) - 200*pow(G5px,3)*pow(H,2) + 840*G4x*G4xx*G5pxx*pow(H,2) - 72*G4*G4xxx*G5pxx*pow(H,2) - 
      708*G4xx*G5p*G5pxx*pow(H,2) - 428*G4x*G5px*G5pxx*pow(H,2) + 324*G5p*G5px*G5pxx*pow(H,2) + 28*G4*pow(G5pxx,2)*pow(H,2) - 120*pow(G4x,2)*G5pxxx*pow(H,2) + 72*G4*G4xx*G5pxxx*pow(H,2) + 
      156*G4x*G5p*G5pxxx*pow(H,2) - 48*pow(G5p,2)*G5pxxx*pow(H,2) - 40*G4*G5px*G5pxxx*pow(H,2) - 792*G4pxxx*G4x*G5x*pow(H,2) - 216*G4px*G4xxx*G5x*pow(H,2) + 504*G4pxxx*G5p*G5x*pow(H,2) + 
      36*G4xxx*G5pp*G5x*pow(H,2) - 884*G4xx*G5ppx*G5x*pow(H,2) + 388*G4x*G5ppxx*G5x*pow(H,2) - 240*G5p*G5ppxx*G5x*pow(H,2) + 498*G5ppx*G5px*G5x*pow(H,2) - 348*G4px*G5pxx*G5x*pow(H,2) - 
      2*G5pp*G5pxx*G5x*pow(H,2) + 16*G4p*G5pxxx*G5x*pow(H,2) - 108*G3pxx*pow(G5x,2)*pow(H,2) + 192*G4ppxx*pow(G5x,2)*pow(H,2) + 12*G5pppx*pow(G5x,2)*pow(H,2) + 
      96*G4*G4pxxx*G5xx*pow(H,2) - 696*G4px*G4xx*G5xx*pow(H,2) + 96*G4p*G4xxx*G5xx*pow(H,2) + 84*G4xx*G5pp*G5xx*pow(H,2) - 356*G4x*G5ppx*G5xx*pow(H,2) + 274*G5p*G5ppx*G5xx*pow(H,2) - 
      40*G4*G5ppxx*G5xx*pow(H,2) + 692*G4px*G5px*G5xx*pow(H,2) - 64*G5pp*G5px*G5xx*pow(H,2) - 100*G4p*G5pxx*G5xx*pow(H,2) + 108*G3px*G5x*G5xx*pow(H,2) - 56*G4ppx*G5x*G5xx*pow(H,2) - 
      44*G5ppp*G5x*G5xx*pow(H,2) - 66*G4pp*pow(G5xx,2)*pow(H,2) - 48*G4px*G4x*G5xxx*pow(H,2) + 72*G4p*G4xx*G5xxx*pow(H,2) + 60*G4px*G5p*G5xxx*pow(H,2) - 18*G5p*G5pp*G5xxx*pow(H,2) + 
      28*G4*G5ppx*G5xxx*pow(H,2) - 24*G4p*G5px*G5xxx*pow(H,2) - 4*G4pp*G5x*G5xxx*pow(H,2) + 87*G5pxx*pow(G5x,2)*pow(H,4) + 80*G5px*G5x*G5xx*pow(H,4) + 
      6*G3xx*(24*G4pxx*G4x - 36*G4px*G4xx - 12*G4pxx*G5p + 18*G4xx*G5pp - 14*G4x*G5ppx + 7*G5p*G5ppx + 36*G4px*G5px - 12*G5pp*G5px + 3*G3px*G5x - 10*G4ppx*G5x + 2*G5ppp*G5x - 5*G5px*G5x*pow(H,2)) - 
      4*G4pxx*(6*G4x*G5ppx - 3*G5p*G5ppx + 12*G4px*(33*G4xx - 14*G5px) - 6*G2xx*G5x - 3*G3px*G5x + 6*G4ppx*G5x + 6*G5ppp*G5x - 6*G2x*G5xx + 12*G3p*G5xx - 12*G4pp*G5xx + 283*G5px*G5x*pow(H,2) - 
         204*G4x*G5xx*pow(H,2) + 168*G5p*G5xx*pow(H,2) + 12*G4*G5xxx*pow(H,2) - 18*G4xx*(G5pp + 30*G5x*pow(H,2)))) + 36*G4xxx*G5pxx*G5x*pow(H,2)*rptot - 18*pow(G5pxx,2)*G5x*pow(H,2)*rptot - 
   24*G4xx*G5pxxx*G5x*pow(H,2)*rptot + 12*G5px*G5pxxx*G5x*pow(H,2)*rptot - 36*G4xxx*G5px*G5xx*pow(H,2)*rptot + 72*G4xx*G5pxx*G5xx*pow(H,2)*rptot - 18*G5px*G5pxx*G5xx*pow(H,2)*rptot - 
   12*G4x*G5pxxx*G5xx*pow(H,2)*rptot + 6*G5p*G5pxxx*G5xx*pow(H,2)*rptot - 36*G4pxxx*G5x*G5xx*pow(H,2)*rptot + 18*G5ppxx*G5x*G5xx*pow(H,2)*rptot + 36*G4pxx*pow(G5xx,2)*pow(H,2)*rptot - 
   18*G5ppx*pow(G5xx,2)*pow(H,2)*rptot - 24*G4xx*G5px*G5xxx*pow(H,2)*rptot + 12*pow(G5px,2)*G5xxx*pow(H,2)*rptot + 12*G4x*G5pxx*G5xxx*pow(H,2)*rptot - 6*G5p*G5pxx*G5xxx*pow(H,2)*rptot + 
   24*G4pxx*G5x*G5xxx*pow(H,2)*rptot - 12*G5ppx*G5x*G5xxx*pow(H,2)*rptot - 12*G4px*G5xx*G5xxx*pow(H,2)*rptot + 6*G5pp*G5xx*G5xxx*pow(H,2)*rptot) + 
pow(a,5)*pow(dphi,12)*pow(H,4)*(3744*G4px*G4pxxx*G4x*G4xx - 864*G3xx*G4px*pow(G4xx,2) - 864*G3pxx*G4x*pow(G4xx,2) + 576*G4ppxx*G4x*pow(G4xx,2) + 864*G3px*pow(G4xx,3) + 
   576*G4ppx*pow(G4xx,3) + 576*G4ppxx*pow(G4x,2)*G4xxx + 3744*pow(G4px,2)*G4xx*G4xxx + 432*G3px*G4x*G4xx*G4xxx - 2016*G4ppx*G4x*G4xx*G4xxx - 1872*G4px*G4pxxx*G4xx*G5p + 
   432*G3pxx*pow(G4xx,2)*G5p - 288*G4ppxx*pow(G4xx,2)*G5p - 576*G4ppxx*G4x*G4xxx*G5p - 216*G3px*G4xx*G4xxx*G5p + 1008*G4ppx*G4xx*G4xxx*G5p + 144*G4ppxx*G4xxx*pow(G5p,2) - 576*G4pxxx*G4x*G4xx*G5pp + 
   432*G3xx*pow(G4xx,2)*G5pp - 2448*G4px*G4xx*G4xxx*G5pp + 288*G4pxxx*G4xx*G5p*G5pp + 288*G4xx*G4xxx*pow(G5pp,2) - 1152*pow(G4xx,3)*G5ppp + 576*G4x*G4xx*G4xxx*G5ppp - 288*G4xx*G4xxx*G5p*G5ppp + 
   576*G4x*pow(G4xx,2)*G5pppx - 288*pow(G4x,2)*G4xxx*G5pppx - 288*pow(G4xx,2)*G5p*G5pppx + 288*G4x*G4xxx*G5p*G5pppx - 72*G4xxx*pow(G5p,2)*G5pppx + 288*G4pxxx*pow(G4x,2)*G5ppx - 
   576*G3xx*G4x*G4xx*G5ppx + 4320*G4px*pow(G4xx,2)*G5ppx + 1584*G4px*G4x*G4xxx*G5ppx - 288*G4pxxx*G4x*G5p*G5ppx + 288*G3xx*G4xx*G5p*G5ppx - 792*G4px*G4xxx*G5p*G5ppx + 72*G4pxxx*pow(G5p,2)*G5ppx - 
   864*pow(G4xx,2)*G5pp*G5ppx - 288*G4x*G4xxx*G5pp*G5ppx + 144*G4xxx*G5p*G5pp*G5ppx + 576*G4x*G4xx*pow(G5ppx,2) - 288*G4xx*G5p*pow(G5ppx,2) + 48*G3xx*pow(G4x,2)*G5ppxx - 
   1680*G4px*G4x*G4xx*G5ppxx - 48*G3xx*G4x*G5p*G5ppxx + 840*G4px*G4xx*G5p*G5ppxx + 12*G3xx*pow(G5p,2)*G5ppxx + 288*G4x*G4xx*G5pp*G5ppxx - 144*G4xx*G5p*G5pp*G5ppxx - 144*pow(G4x,2)*G5ppx*G5ppxx + 
   144*G4x*G5p*G5ppx*G5ppxx - 36*pow(G5p,2)*G5ppx*G5ppxx - 2016*G4px*G4pxxx*G4x*G5px + 1944*G3xx*G4px*G4xx*G5px + 1008*G3pxx*G4x*G4xx*G5px - 864*G4ppxx*G4x*G4xx*G5px - 144*G2xx*pow(G4xx,2)*G5px - 
   1584*G3px*pow(G4xx,2)*G5px - 2736*pow(G4px,2)*G4xxx*G5px - 288*G3px*G4x*G4xxx*G5px + 1152*G4ppx*G4x*G4xxx*G5px - 144*G2x*G4xx*G4xxx*G5px + 288*G3p*G4xx*G4xxx*G5px - 288*G4pp*G4xx*G4xxx*G5px + 
   1008*G4px*G4pxxx*G5p*G5px - 504*G3pxx*G4xx*G5p*G5px + 432*G4ppxx*G4xx*G5p*G5px + 144*G3px*G4xxx*G5p*G5px - 576*G4ppx*G4xxx*G5p*G5px + 288*G4pxxx*G4x*G5pp*G5px - 648*G3xx*G4xx*G5pp*G5px + 
   1440*G4px*G4xxx*G5pp*G5px - 144*G4pxxx*G5p*G5pp*G5px - 144*G4xxx*pow(G5pp,2)*G5px + 1728*pow(G4xx,2)*G5ppp*G5px - 288*G4x*G4xxx*G5ppp*G5px + 144*G4xxx*G5p*G5ppp*G5px - 576*G4x*G4xx*G5pppx*G5px + 
   288*G4xx*G5p*G5pppx*G5px + 360*G3xx*G4x*G5ppx*G5px - 4896*G4px*G4xx*G5ppx*G5px - 180*G3xx*G5p*G5ppx*G5px + 864*G4xx*G5pp*G5ppx*G5px - 288*G4x*pow(G5ppx,2)*G5px + 144*G5p*pow(G5ppx,2)*G5px + 
   912*G4px*G4x*G5ppxx*G5px - 456*G4px*G5p*G5ppxx*G5px - 144*G4x*G5pp*G5ppxx*G5px + 72*G5p*G5pp*G5ppxx*G5px - 792*G3xx*G4px*pow(G5px,2) - 288*G3pxx*G4x*pow(G5px,2) + 288*G4ppxx*G4x*pow(G5px,2) + 
   168*G2xx*G4xx*pow(G5px,2) + 912*G3px*G4xx*pow(G5px,2) - 432*G4ppx*G4xx*pow(G5px,2) + 72*G2x*G4xxx*pow(G5px,2) - 144*G3p*G4xxx*pow(G5px,2) + 144*G4pp*G4xxx*pow(G5px,2) + 
   144*G3pxx*G5p*pow(G5px,2) - 144*G4ppxx*G5p*pow(G5px,2) + 216*G3xx*G5pp*pow(G5px,2) - 864*G4xx*G5ppp*pow(G5px,2) + 144*G4x*G5pppx*pow(G5px,2) - 72*G5p*G5pppx*pow(G5px,2) + 
   1368*G4px*G5ppx*pow(G5px,2) - 216*G5pp*G5ppx*pow(G5px,2) - 48*G2xx*pow(G5px,3) - 168*G3px*pow(G5px,3) + 144*G4ppx*pow(G5px,3) + 144*G5ppp*pow(G5px,3) - 360*G3xx*G4px*G4x*G5pxx - 
   48*G3pxx*pow(G4x,2)*G5pxx - 192*G4ppxx*pow(G4x,2)*G5pxx - 24*pow(G4px,2)*G4xx*G5pxx + 72*G2xx*G4x*G4xx*G5pxx + 96*G3px*G4x*G4xx*G5pxx + 240*G4ppx*G4x*G4xx*G5pxx + 144*G2x*pow(G4xx,2)*G5pxx - 
   288*G3p*pow(G4xx,2)*G5pxx + 288*G4pp*pow(G4xx,2)*G5pxx + 72*G2x*G4x*G4xxx*G5pxx - 144*G3p*G4x*G4xxx*G5pxx + 144*G4pp*G4x*G4xxx*G5pxx + 180*G3xx*G4px*G5p*G5pxx + 48*G3pxx*G4x*G5p*G5pxx + 
   192*G4ppxx*G4x*G5p*G5pxx - 36*G2xx*G4xx*G5p*G5pxx - 48*G3px*G4xx*G5p*G5pxx - 120*G4ppx*G4xx*G5p*G5pxx - 36*G2x*G4xxx*G5p*G5pxx + 72*G3p*G4xxx*G5p*G5pxx - 72*G4pp*G4xxx*G5p*G5pxx - 
   12*G3pxx*pow(G5p,2)*G5pxx - 48*G4ppxx*pow(G5p,2)*G5pxx + 72*G3xx*G4x*G5pp*G5pxx + 840*G4px*G4xx*G5pp*G5pxx - 36*G3xx*G5p*G5pp*G5pxx - 144*G4xx*pow(G5pp,2)*G5pxx - 288*G4x*G4xx*G5ppp*G5pxx + 
   144*G4xx*G5p*G5ppp*G5pxx + 144*pow(G4x,2)*G5pppx*G5pxx - 144*G4x*G5p*G5pppx*G5pxx + 36*pow(G5p,2)*G5pppx*G5pxx - 552*G4px*G4x*G5ppx*G5pxx + 276*G4px*G5p*G5ppx*G5pxx + 144*G4x*G5pp*G5ppx*G5pxx - 
   72*G5p*G5pp*G5ppx*G5pxx + 384*pow(G4px,2)*G5px*G5pxx - 48*G2xx*G4x*G5px*G5pxx - 192*G4ppx*G4x*G5px*G5pxx - 72*G2x*G4xx*G5px*G5pxx + 144*G3p*G4xx*G5px*G5pxx - 144*G4pp*G4xx*G5px*G5pxx + 
   24*G2xx*G5p*G5px*G5pxx + 96*G4ppx*G5p*G5px*G5pxx - 528*G4px*G5pp*G5px*G5pxx + 72*pow(G5pp,2)*G5px*G5pxx + 144*G4x*G5ppp*G5px*G5pxx - 72*G5p*G5ppp*G5px*G5pxx - 36*G2x*G4x*pow(G5pxx,2) + 
   72*G3p*G4x*pow(G5pxx,2) - 72*G4pp*G4x*pow(G5pxx,2) + 18*G2x*G5p*pow(G5pxx,2) - 36*G3p*G5p*pow(G5pxx,2) + 36*G4pp*G5p*pow(G5pxx,2) - 240*pow(G4px,2)*G4x*G5pxxx - 
   48*G3px*pow(G4x,2)*G5pxxx + 96*G4ppx*pow(G4x,2)*G5pxxx - 48*G2x*G4x*G4xx*G5pxxx + 96*G3p*G4x*G4xx*G5pxxx - 96*G4pp*G4x*G4xx*G5pxxx + 120*pow(G4px,2)*G5p*G5pxxx + 48*G3px*G4x*G5p*G5pxxx - 
   96*G4ppx*G4x*G5p*G5pxxx + 24*G2x*G4xx*G5p*G5pxxx - 48*G3p*G4xx*G5p*G5pxxx + 48*G4pp*G4xx*G5p*G5pxxx - 12*G3px*pow(G5p,2)*G5pxxx + 24*G4ppx*pow(G5p,2)*G5pxxx + 48*G4px*G4x*G5pp*G5pxxx - 
   24*G4px*G5p*G5pp*G5pxxx + 24*G2x*G4x*G5px*G5pxxx - 48*G3p*G4x*G5px*G5pxxx + 48*G4pp*G4x*G5px*G5pxxx - 12*G2x*G5p*G5px*G5pxxx + 24*G3p*G5p*G5px*G5pxxx - 24*G4pp*G5p*G5px*G5pxxx - 
   720*pow(G4px,2)*G4pxxx*G5x + 288*G3xx*G4ppxx*G4x*G5x - 288*G3px*G4pxxx*G4x*G5x + 576*G4ppx*G4pxxx*G4x*G5x + 108*G3px*G3xx*G4xx*G5x - 504*G3xx*G4ppx*G4xx*G5x + 936*G3pxx*G4px*G4xx*G5x - 
   1296*G4ppxx*G4px*G4xx*G5x - 144*G2x*G4pxxx*G4xx*G5x + 288*G3p*G4pxxx*G4xx*G5x - 288*G4pp*G4pxxx*G4xx*G5x - 144*G2pxx*pow(G4xx,2)*G5x - 144*G3ppx*pow(G4xx,2)*G5x + 576*G4pppx*pow(G4xx,2)*G5x - 
   360*G3px*G4px*G4xxx*G5x + 1008*G4ppx*G4px*G4xxx*G5x + 288*G3ppx*G4x*G4xxx*G5x - 576*G4pppx*G4x*G4xxx*G5x - 144*G3pp*G4xx*G4xxx*G5x + 288*G4ppp*G4xx*G4xxx*G5x - 144*G3xx*G4ppxx*G5p*G5x + 
   144*G3px*G4pxxx*G5p*G5x - 288*G4ppx*G4pxxx*G5p*G5x - 144*G3ppx*G4xxx*G5p*G5x + 288*G4pppx*G4xxx*G5p*G5x + 144*G4px*G4pxxx*G5pp*G5x - 144*G3pxx*G4xx*G5pp*G5x + 288*G4ppxx*G4xx*G5pp*G5x + 
   72*G3px*G4xxx*G5pp*G5x - 144*G4ppx*G4xxx*G5pp*G5x + 144*G3xx*G4xx*G5ppp*G5x - 144*G4px*G4xxx*G5ppp*G5x - 144*G3xx*G4x*G5pppx*G5x - 288*G4px*G4xx*G5pppx*G5x + 72*G3xx*G5p*G5pppx*G5x + 
   396*G3xx*G4px*G5ppx*G5x + 144*G3pxx*G4x*G5ppx*G5x - 288*G4ppxx*G4x*G5ppx*G5x - 96*G2xx*G4xx*G5ppx*G5x - 336*G3px*G4xx*G5ppx*G5x + 864*G4ppx*G4xx*G5ppx*G5x - 72*G2x*G4xxx*G5ppx*G5x + 
   144*G3p*G4xxx*G5ppx*G5x - 144*G4pp*G4xxx*G5ppx*G5x - 72*G3pxx*G5p*G5ppx*G5x + 144*G4ppxx*G5p*G5ppx*G5x - 72*G3xx*G5pp*G5ppx*G5x - 144*G4px*pow(G5ppx,2)*G5x + 336*pow(G4px,2)*G5ppxx*G5x + 
   16*G2xx*G4x*G5ppxx*G5x + 128*G3px*G4x*G5ppxx*G5x - 288*G4ppx*G4x*G5ppxx*G5x + 72*G2x*G4xx*G5ppxx*G5x - 144*G3p*G4xx*G5ppxx*G5x + 144*G4pp*G4xx*G5ppxx*G5x - 8*G2xx*G5p*G5ppxx*G5x - 
   64*G3px*G5p*G5ppxx*G5x + 144*G4ppx*G5p*G5ppxx*G5x - 72*G4px*G5pp*G5ppxx*G5x - 72*G3px*G3xx*G5px*G5x + 288*G3xx*G4ppx*G5px*G5x - 504*G3pxx*G4px*G5px*G5x + 720*G4ppxx*G4px*G5px*G5x + 
   72*G2x*G4pxxx*G5px*G5x - 144*G3p*G4pxxx*G5px*G5x + 144*G4pp*G4pxxx*G5px*G5x + 168*G2pxx*G4xx*G5px*G5x + 120*G3ppx*G4xx*G5px*G5x - 576*G4pppx*G4xx*G5px*G5x + 72*G3pp*G4xxx*G5px*G5x - 
   144*G4ppp*G4xxx*G5px*G5x + 72*G3pxx*G5pp*G5px*G5x - 144*G4ppxx*G5pp*G5px*G5x - 72*G3xx*G5ppp*G5px*G5x + 144*G4px*G5pppx*G5px*G5x + 60*G2xx*G5ppx*G5px*G5x + 156*G3px*G5ppx*G5px*G5x - 
   432*G4ppx*G5ppx*G5px*G5x - 36*G2x*G5ppxx*G5px*G5x + 72*G3p*G5ppxx*G5px*G5x - 72*G4pp*G5ppxx*G5px*G5x - 48*G2pxx*pow(G5px,2)*G5x - 24*G3ppx*pow(G5px,2)*G5x + 144*G4pppx*pow(G5px,2)*G5x + 
   18*G2x*G3xx*G5pxx*G5x - 36*G3p*G3xx*G5pxx*G5x + 36*G3xx*G4pp*G5pxx*G5x - 60*G2xx*G4px*G5pxx*G5x + 144*G3px*G4px*G5pxx*G5x - 312*G4ppx*G4px*G5pxx*G5x - 16*G2pxx*G4x*G5pxx*G5x - 
   128*G3ppx*G4x*G5pxx*G5x + 288*G4pppx*G4x*G5pxx*G5x + 48*G2px*G4xx*G5pxx*G5x + 24*G3pp*G4xx*G5pxx*G5x - 144*G4ppp*G4xx*G5pxx*G5x + 8*G2pxx*G5p*G5pxx*G5x + 64*G3ppx*G5p*G5pxx*G5x - 
   144*G4pppx*G5p*G5pxx*G5x + 12*G2xx*G5pp*G5pxx*G5x - 48*G3px*G5pp*G5pxx*G5x + 72*G4ppx*G5pp*G5pxx*G5x + 72*G4px*G5ppp*G5pxx*G5x + 36*G2x*G5ppx*G5pxx*G5x - 72*G3p*G5ppx*G5pxx*G5x + 
   72*G4pp*G5ppx*G5pxx*G5x - 24*G2px*G5px*G5pxx*G5x - 12*G3pp*G5px*G5pxx*G5x + 72*G4ppp*G5px*G5pxx*G5x + 12*G2x*G4px*G5pxxx*G5x - 24*G3p*G4px*G5pxxx*G5x + 24*G4pp*G4px*G5pxxx*G5x - 
   16*G2px*G4x*G5pxxx*G5x + 16*G3pp*G4x*G5pxxx*G5x + 8*G2px*G5p*G5pxxx*G5x - 8*G3pp*G5p*G5pxxx*G5x - 36*G3px*G3pxx*pow(G5x,2) + 36*G3ppx*G3xx*pow(G5x,2) - 72*G3xx*G4pppx*pow(G5x,2) + 
   72*G3pxx*G4ppx*pow(G5x,2) + 24*G2xx*G4ppxx*pow(G5x,2) + 48*G3px*G4ppxx*pow(G5x,2) - 144*G4ppx*G4ppxx*pow(G5x,2) - 24*G2px*G4pxxx*pow(G5x,2) + 24*G3pp*G4pxxx*pow(G5x,2) + 
   24*G2ppx*G4xxx*pow(G5x,2) - 24*G3ppp*G4xxx*pow(G5x,2) - 12*G2xx*G5pppx*pow(G5x,2) + 12*G3px*G5pppx*pow(G5x,2) + 12*G2pxx*G5ppx*pow(G5x,2) - 12*G3ppx*G5ppx*pow(G5x,2) + 
   12*G2px*G5ppxx*pow(G5x,2) - 12*G3pp*G5ppxx*pow(G5x,2) - 12*G2ppx*G5pxx*pow(G5x,2) + 12*G3ppp*G5pxx*pow(G5x,2) + 432*G3xx*pow(G4px,2)*G5xx + 36*G3px*G3xx*G4x*G5xx - 216*G3xx*G4ppx*G4x*G5xx + 
   432*G3pxx*G4px*G4x*G5xx - 576*G4ppxx*G4px*G4x*G5xx - 72*G2x*G4pxxx*G4x*G5xx + 144*G3p*G4pxxx*G4x*G5xx - 144*G4pp*G4pxxx*G4x*G5xx - 120*G2xx*G4px*G4xx*G5xx - 1248*G3px*G4px*G4xx*G5xx + 
   1008*G4ppx*G4px*G4xx*G5xx - 120*G2pxx*G4x*G4xx*G5xx - 168*G3ppx*G4x*G4xx*G5xx + 576*G4pppx*G4x*G4xx*G5xx + 432*G3pp*pow(G4xx,2)*G5xx - 864*G4ppp*pow(G4xx,2)*G5xx - 72*G2x*G4px*G4xxx*G5xx + 
   144*G3p*G4px*G4xxx*G5xx - 144*G4pp*G4px*G4xxx*G5xx - 72*G3pp*G4x*G4xxx*G5xx + 144*G4ppp*G4x*G4xxx*G5xx - 18*G3px*G3xx*G5p*G5xx + 108*G3xx*G4ppx*G5p*G5xx - 216*G3pxx*G4px*G5p*G5xx + 
   288*G4ppxx*G4px*G5p*G5xx + 36*G2x*G4pxxx*G5p*G5xx - 72*G3p*G4pxxx*G5p*G5xx + 72*G4pp*G4pxxx*G5p*G5xx + 60*G2pxx*G4xx*G5p*G5xx + 84*G3ppx*G4xx*G5p*G5xx - 288*G4pppx*G4xx*G5p*G5xx + 
   36*G3pp*G4xxx*G5p*G5xx - 72*G4ppp*G4xxx*G5p*G5xx - 288*G3xx*G4px*G5pp*G5xx - 72*G3pxx*G4x*G5pp*G5xx + 144*G4ppxx*G4x*G5pp*G5xx + 60*G2xx*G4xx*G5pp*G5xx + 228*G3px*G4xx*G5pp*G5xx - 
   576*G4ppx*G4xx*G5pp*G5xx + 36*G2x*G4xxx*G5pp*G5xx - 72*G3p*G4xxx*G5pp*G5xx + 72*G4pp*G4xxx*G5pp*G5xx + 36*G3pxx*G5p*G5pp*G5xx - 72*G4ppxx*G5p*G5pp*G5xx + 36*G3xx*pow(G5pp,2)*G5xx + 
   72*G3xx*G4x*G5ppp*G5xx + 864*G4px*G4xx*G5ppp*G5xx - 36*G3xx*G5p*G5ppp*G5xx - 144*G4px*G4x*G5pppx*G5xx + 72*G4px*G5p*G5pppx*G5xx - 828*pow(G4px,2)*G5ppx*G5xx - 36*G2xx*G4x*G5ppx*G5xx - 
   180*G3px*G4x*G5ppx*G5xx + 432*G4ppx*G4x*G5ppx*G5xx - 144*G2x*G4xx*G5ppx*G5xx + 288*G3p*G4xx*G5ppx*G5xx - 288*G4pp*G4xx*G5ppx*G5xx + 18*G2xx*G5p*G5ppx*G5xx + 90*G3px*G5p*G5ppx*G5xx - 
   216*G4ppx*G5p*G5ppx*G5xx + 216*G4px*G5pp*G5ppx*G5xx + 36*G2x*G4x*G5ppxx*G5xx - 72*G3p*G4x*G5ppxx*G5xx + 72*G4pp*G4x*G5ppxx*G5xx - 18*G2x*G5p*G5ppxx*G5xx + 36*G3p*G5p*G5ppxx*G5xx - 
   36*G4pp*G5p*G5ppxx*G5xx - 18*G2x*G3xx*G5px*G5xx + 36*G3p*G3xx*G5px*G5xx - 36*G3xx*G4pp*G5px*G5xx + 144*G2xx*G4px*G5px*G5xx + 648*G3px*G4px*G5px*G5xx - 720*G4ppx*G4px*G5px*G5xx + 
   72*G2pxx*G4x*G5px*G5xx + 72*G3ppx*G4x*G5px*G5xx - 288*G4pppx*G4x*G5px*G5xx - 48*G2px*G4xx*G5px*G5xx - 384*G3pp*G4xx*G5px*G5xx + 864*G4ppp*G4xx*G5px*G5xx - 36*G2pxx*G5p*G5px*G5xx - 
   36*G3ppx*G5p*G5px*G5xx + 144*G4pppx*G5p*G5px*G5xx - 48*G2xx*G5pp*G5px*G5xx - 96*G3px*G5pp*G5px*G5xx + 288*G4ppx*G5pp*G5px*G5xx - 432*G4px*G5ppp*G5px*G5xx + 72*G2x*G5ppx*G5px*G5xx - 
   144*G3p*G5ppx*G5px*G5xx + 144*G4pp*G5ppx*G5px*G5xx + 24*G2px*pow(G5px,2)*G5xx + 84*G3pp*pow(G5px,2)*G5xx - 216*G4ppp*pow(G5px,2)*G5xx + 24*G2px*G4x*G5pxx*G5xx + 12*G3pp*G4x*G5pxx*G5xx - 
   72*G4ppp*G4x*G5pxx*G5xx - 12*G2px*G5p*G5pxx*G5xx - 6*G3pp*G5p*G5pxx*G5xx + 36*G4ppp*G5p*G5pxx*G5xx - 18*G2x*G5pp*G5pxx*G5xx + 36*G3p*G5pp*G5pxx*G5xx - 36*G4pp*G5pp*G5pxx*G5xx + 6*G2xx*G3px*G5x*G5xx + 
   30*pow(G3px,2)*G5x*G5xx - 18*G2x*G3pxx*G5x*G5xx + 36*G3p*G3pxx*G5x*G5xx - 18*G3pp*G3xx*G5x*G5xx - 36*G3pxx*G4pp*G5x*G5xx + 36*G3xx*G4ppp*G5x*G5xx - 36*G2xx*G4ppx*G5x*G5xx - 108*G3px*G4ppx*G5x*G5xx + 
   144*pow(G4ppx,2)*G5x*G5xx + 36*G2x*G4ppxx*G5x*G5xx - 72*G3p*G4ppxx*G5x*G5xx + 72*G4pp*G4ppxx*G5x*G5xx + 72*G2pxx*G4px*G5x*G5xx - 144*G4pppx*G4px*G5x*G5xx - 48*G2ppx*G4xx*G5x*G5xx + 
   48*G3ppp*G4xx*G5x*G5xx - 12*G2pxx*G5pp*G5x*G5xx + 12*G3ppx*G5pp*G5x*G5xx + 12*G2xx*G5ppp*G5x*G5xx - 12*G3px*G5ppp*G5x*G5xx - 24*G2px*G5ppx*G5x*G5xx + 24*G3pp*G5ppx*G5x*G5xx + 24*G2ppx*G5px*G5x*G5xx - 
   24*G3ppp*G5px*G5x*G5xx + 9*G2x*G3px*pow(G5xx,2) - 18*G3p*G3px*pow(G5xx,2) + 18*G3px*G4pp*pow(G5xx,2) - 18*G2x*G4ppx*pow(G5xx,2) + 36*G3p*G4ppx*pow(G5xx,2) - 36*G4pp*G4ppx*pow(G5xx,2) - 
   12*G2px*G4px*pow(G5xx,2) - 42*G3pp*G4px*pow(G5xx,2) + 108*G4ppp*G4px*pow(G5xx,2) - 12*G2ppx*G4x*pow(G5xx,2) + 12*G3ppp*G4x*pow(G5xx,2) + 6*G2ppx*G5p*pow(G5xx,2) - 
   6*G3ppp*G5p*pow(G5xx,2) + 6*G2px*G5pp*pow(G5xx,2) - 6*G3pp*G5pp*pow(G5xx,2) - 240*pow(G4px,3)*G5xxx - 120*G3px*G4px*G4x*G5xxx + 336*G4ppx*G4px*G4x*G5xxx + 48*G3ppx*pow(G4x,2)*G5xxx - 
   96*G4pppx*pow(G4x,2)*G5xxx - 48*G2x*G4px*G4xx*G5xxx + 96*G3p*G4px*G4xx*G5xxx - 96*G4pp*G4px*G4xx*G5xxx - 48*G3pp*G4x*G4xx*G5xxx + 96*G4ppp*G4x*G4xx*G5xxx + 60*G3px*G4px*G5p*G5xxx - 
   168*G4ppx*G4px*G5p*G5xxx - 48*G3ppx*G4x*G5p*G5xxx + 96*G4pppx*G4x*G5p*G5xxx + 24*G3pp*G4xx*G5p*G5xxx - 48*G4ppp*G4xx*G5p*G5xxx + 12*G3ppx*pow(G5p,2)*G5xxx - 24*G4pppx*pow(G5p,2)*G5xxx + 
   168*pow(G4px,2)*G5pp*G5xxx + 24*G3px*G4x*G5pp*G5xxx - 48*G4ppx*G4x*G5pp*G5xxx + 24*G2x*G4xx*G5pp*G5xxx - 48*G3p*G4xx*G5pp*G5xxx + 48*G4pp*G4xx*G5pp*G5xxx - 12*G3px*G5p*G5pp*G5xxx + 
   24*G4ppx*G5p*G5pp*G5xxx - 24*G4px*pow(G5pp,2)*G5xxx - 48*G4px*G4x*G5ppp*G5xxx + 24*G4px*G5p*G5ppp*G5xxx - 24*G2x*G4x*G5ppx*G5xxx + 48*G3p*G4x*G5ppx*G5xxx - 48*G4pp*G4x*G5ppx*G5xxx + 
   12*G2x*G5p*G5ppx*G5xxx - 24*G3p*G5p*G5ppx*G5xxx + 24*G4pp*G5p*G5ppx*G5xxx + 36*G2x*G4px*G5px*G5xxx - 72*G3p*G4px*G5px*G5xxx + 72*G4pp*G4px*G5px*G5xxx + 24*G3pp*G4x*G5px*G5xxx - 
   48*G4ppp*G4x*G5px*G5xxx - 12*G3pp*G5p*G5px*G5xxx + 24*G4ppp*G5p*G5px*G5xxx - 12*G2x*G5pp*G5px*G5xxx + 24*G3p*G5pp*G5px*G5xxx - 24*G4pp*G5pp*G5px*G5xxx + 6*G2x*G3px*G5x*G5xxx - 12*G3p*G3px*G5x*G5xxx + 
   12*G3px*G4pp*G5x*G5xxx - 12*G2x*G4ppx*G5x*G5xxx + 24*G3p*G4ppx*G5x*G5xxx - 24*G4pp*G4ppx*G5x*G5xxx + 12*G3pp*G4px*G5x*G5xxx - 24*G4ppp*G4px*G5x*G5xxx + 16*G2ppx*G4x*G5x*G5xxx - 
   16*G3ppp*G4x*G5x*G5xxx - 8*G2ppx*G5p*G5x*G5xxx + 8*G3ppp*G5p*G5x*G5xxx - 3*pow(G3x,2)*
    (36*G4xxx*G5px - 90*G4xx*G5pxx + 32*G5px*G5pxx + 12*G4x*G5pxxx - 6*G5p*G5pxxx + 36*G4pxxx*G5x - 16*G5ppxx*G5x + 33*G5ppx*G5xx + 12*G4px*G5xxx - 6*G5pp*G5xxx) - 1440*pow(G4xx,3)*G5px*pow(H,2) - 
   1008*G4x*G4xx*G4xxx*G5px*pow(H,2) + 1368*G4xx*G4xxx*G5p*G5px*pow(H,2) + 2544*pow(G4xx,2)*pow(G5px,2)*pow(H,2) + 432*G4x*G4xxx*pow(G5px,2)*pow(H,2) - 
   720*G4xxx*G5p*pow(G5px,2)*pow(H,2) - 1536*G4xx*pow(G5px,3)*pow(H,2) + 312*pow(G5px,4)*pow(H,2) + 3456*G4x*pow(G4xx,2)*G5pxx*pow(H,2) + 720*pow(G4x,2)*G4xxx*G5pxx*pow(H,2) - 
   432*G4*G4xx*G4xxx*G5pxx*pow(H,2) - 3024*pow(G4xx,2)*G5p*G5pxx*pow(H,2) - 936*G4x*G4xxx*G5p*G5pxx*pow(H,2) + 288*G4xxx*pow(G5p,2)*G5pxx*pow(H,2) - 3120*G4x*G4xx*G5px*G5pxx*pow(H,2) + 
   288*G4*G4xxx*G5px*G5pxx*pow(H,2) + 2664*G4xx*G5p*G5px*G5pxx*pow(H,2) + 768*G4x*pow(G5px,2)*G5pxx*pow(H,2) - 576*G5p*pow(G5px,2)*G5pxx*pow(H,2) - 
   376*pow(G4x,2)*pow(G5pxx,2)*pow(H,2) + 120*G4*G4xx*pow(G5pxx,2)*pow(H,2) + 436*G4x*G5p*pow(G5pxx,2)*pow(H,2) - 124*pow(G5p,2)*pow(G5pxx,2)*pow(H,2) - 
   96*G4*G5px*pow(G5pxx,2)*pow(H,2) - 864*pow(G4x,2)*G4xx*G5pxxx*pow(H,2) + 288*G4*pow(G4xx,2)*G5pxxx*pow(H,2) + 1152*G4x*G4xx*G5p*G5pxxx*pow(H,2) - 
   360*G4xx*pow(G5p,2)*G5pxxx*pow(H,2) + 448*pow(G4x,2)*G5px*G5pxxx*pow(H,2) - 336*G4*G4xx*G5px*G5pxxx*pow(H,2) - 616*G4x*G5p*G5px*G5pxxx*pow(H,2) + 196*pow(G5p,2)*G5px*G5pxxx*pow(H,2) + 
   96*G4*pow(G5px,2)*G5pxxx*pow(H,2) + 32*G4*G4x*G5pxx*G5pxxx*pow(H,2) - 16*G4*G5p*G5pxx*G5pxxx*pow(H,2) - 5616*G4pxxx*G4x*G4xx*G5x*pow(H,2) - 2160*G4px*G4xx*G4xxx*G5x*pow(H,2) + 
   3672*G4pxxx*G4xx*G5p*G5x*pow(H,2) + 648*G4xx*G4xxx*G5pp*G5x*pow(H,2) - 3360*pow(G4xx,2)*G5ppx*G5x*pow(H,2) - 1896*G4x*G4xxx*G5ppx*G5x*pow(H,2) + 1236*G4xxx*G5p*G5ppx*G5x*pow(H,2) + 
   2712*G4x*G4xx*G5ppxx*G5x*pow(H,2) - 96*G4*G4xxx*G5ppxx*G5x*pow(H,2) - 1692*G4xx*G5p*G5ppxx*G5x*pow(H,2) + 2976*G4pxxx*G4x*G5px*G5x*pow(H,2) - 360*G3xx*G4xx*G5px*G5x*pow(H,2) + 
   1512*G4px*G4xxx*G5px*G5x*pow(H,2) - 1992*G4pxxx*G5p*G5px*G5x*pow(H,2) - 288*G4xxx*G5pp*G5px*G5x*pow(H,2) + 3816*G4xx*G5ppx*G5px*G5x*pow(H,2) - 1488*G4x*G5ppxx*G5px*G5x*pow(H,2) + 
   948*G5p*G5ppxx*G5px*G5x*pow(H,2) + 180*G3xx*pow(G5px,2)*G5x*pow(H,2) - 1104*G5ppx*pow(G5px,2)*G5x*pow(H,2) + 96*G4*G4pxxx*G5pxx*G5x*pow(H,2) + 360*G3xx*G4x*G5pxx*G5x*pow(H,2) - 
   3168*G4px*G4xx*G5pxx*G5x*pow(H,2) - 72*G4p*G4xxx*G5pxx*G5x*pow(H,2) - 234*G3xx*G5p*G5pxx*G5x*pow(H,2) - 132*G4xx*G5pp*G5pxx*G5x*pow(H,2) + 924*G4x*G5ppx*G5pxx*G5x*pow(H,2) - 
   486*G5p*G5ppx*G5pxx*G5x*pow(H,2) + 1632*G4px*G5px*G5pxx*G5x*pow(H,2) + 12*G5pp*G5px*G5pxx*G5x*pow(H,2) + 12*G4p*pow(G5pxx,2)*G5x*pow(H,2) + 1032*G4px*G4x*G5pxxx*G5x*pow(H,2) + 
   120*G4p*G4xx*G5pxxx*G5x*pow(H,2) - 672*G4px*G5p*G5pxxx*G5x*pow(H,2) - 72*G4x*G5pp*G5pxxx*G5x*pow(H,2) + 60*G5p*G5pp*G5pxxx*G5x*pow(H,2) - 48*G4*G5ppx*G5pxxx*G5x*pow(H,2) - 
   72*G4p*G5px*G5pxxx*G5x*pow(H,2) + 1728*G4px*G4pxxx*pow(G5x,2)*pow(H,2) - 756*G3pxx*G4xx*pow(G5x,2)*pow(H,2) + 1224*G4ppxx*G4xx*pow(G5x,2)*pow(H,2) + 
   324*G3px*G4xxx*pow(G5x,2)*pow(H,2) - 1080*G4ppx*G4xxx*pow(G5x,2)*pow(H,2) - 144*G4pxxx*G5pp*pow(G5x,2)*pow(H,2) + 144*G4xxx*G5ppp*pow(G5x,2)*pow(H,2) + 
   144*G4xx*G5pppx*pow(G5x,2)*pow(H,2) - 246*G3xx*G5ppx*pow(G5x,2)*pow(H,2) + 36*pow(G5ppx,2)*pow(G5x,2)*pow(H,2) - 864*G4px*G5ppxx*pow(G5x,2)*pow(H,2) + 
   72*G5pp*G5ppxx*pow(G5x,2)*pow(H,2) + 408*G3pxx*G5px*pow(G5x,2)*pow(H,2) - 744*G4ppxx*G5px*pow(G5x,2)*pow(H,2) - 36*G5pppx*G5px*pow(G5x,2)*pow(H,2) + 
   30*G2xx*G5pxx*pow(G5x,2)*pow(H,2) - 156*G3px*G5pxx*pow(G5x,2)*pow(H,2) + 492*G4ppx*G5pxx*pow(G5x,2)*pow(H,2) - 72*G5ppp*G5pxx*pow(G5x,2)*pow(H,2) - 
   30*G2x*G5pxxx*pow(G5x,2)*pow(H,2) + 60*G3p*G5pxxx*pow(G5x,2)*pow(H,2) - 36*G4pp*G5pxxx*pow(G5x,2)*pow(H,2) - 1008*G4pxxx*pow(G4x,2)*G5xx*pow(H,2) + 720*G4*G4pxxx*G4xx*G5xx*pow(H,2) - 
   2304*G4px*pow(G4xx,2)*G5xx*pow(H,2) - 864*G4px*G4x*G4xxx*G5xx*pow(H,2) + 720*G4p*G4xx*G4xxx*G5xx*pow(H,2) + 1368*G4pxxx*G4x*G5p*G5xx*pow(H,2) + 792*G4px*G4xxx*G5p*G5xx*pow(H,2) - 
   432*G4pxxx*pow(G5p,2)*G5xx*pow(H,2) + 288*G4x*G4xxx*G5pp*G5xx*pow(H,2) - 324*G4xxx*G5p*G5pp*G5xx*pow(H,2) - 2616*G4x*G4xx*G5ppx*G5xx*pow(H,2) + 216*G4*G4xxx*G5ppx*G5xx*pow(H,2) + 
   2028*G4xx*G5p*G5ppx*G5xx*pow(H,2) + 520*pow(G4x,2)*G5ppxx*G5xx*pow(H,2) - 264*G4*G4xx*G5ppxx*G5xx*pow(H,2) - 652*G4x*G5p*G5ppxx*G5xx*pow(H,2) + 196*pow(G5p,2)*G5ppxx*G5xx*pow(H,2) - 
   432*G4*G4pxxx*G5px*G5xx*pow(H,2) - 144*G3xx*G4x*G5px*G5xx*pow(H,2) + 5280*G4px*G4xx*G5px*G5xx*pow(H,2) - 288*G4p*G4xxx*G5px*G5xx*pow(H,2) + 162*G3xx*G5p*G5px*G5xx*pow(H,2) - 
   156*G4xx*G5pp*G5px*G5xx*pow(H,2) + 1380*G4x*G5ppx*G5px*G5xx*pow(H,2) - 1194*G5p*G5ppx*G5px*G5xx*pow(H,2) + 168*G4*G5ppxx*G5px*G5xx*pow(H,2) - 2352*G4px*pow(G5px,2)*G5xx*pow(H,2) + 
   132*G5pp*pow(G5px,2)*G5xx*pow(H,2) - 36*G3xx*G4*G5pxx*G5xx*pow(H,2) - 1128*G4px*G4x*G5pxx*G5xx*pow(H,2) - 768*G4p*G4xx*G5pxx*G5xx*pow(H,2) + 972*G4px*G5p*G5pxx*G5xx*pow(H,2) - 
   120*G4x*G5pp*G5pxx*G5xx*pow(H,2) + 54*G5p*G5pp*G5pxx*G5xx*pow(H,2) + 12*G4*G5ppx*G5pxx*G5xx*pow(H,2) + 408*G4p*G5px*G5pxx*G5xx*pow(H,2) - 144*G4*G4px*G5pxxx*G5xx*pow(H,2) + 
   48*G4p*G4x*G5pxxx*G5xx*pow(H,2) - 24*G4p*G5p*G5pxxx*G5xx*pow(H,2) + 24*G4*G5pp*G5pxxx*G5xx*pow(H,2) - 252*G3xx*G4px*G5x*G5xx*pow(H,2) + 144*G4p*G4pxxx*G5x*G5xx*pow(H,2) - 
   540*G3pxx*G4x*G5x*G5xx*pow(H,2) + 1032*G4ppxx*G4x*G5x*G5xx*pow(H,2) + 720*G3px*G4xx*G5x*G5xx*pow(H,2) + 264*G4ppx*G4xx*G5x*G5xx*pow(H,2) + 72*G4pp*G4xxx*G5x*G5xx*pow(H,2) + 
   360*G3pxx*G5p*G5x*G5xx*pow(H,2) - 552*G4ppxx*G5p*G5x*G5xx*pow(H,2) + 90*G3xx*G5pp*G5x*G5xx*pow(H,2) - 528*G4xx*G5ppp*G5x*G5xx*pow(H,2) + 24*G4x*G5pppx*G5x*G5xx*pow(H,2) - 
   84*G5p*G5pppx*G5x*G5xx*pow(H,2) + 1836*G4px*G5ppx*G5x*G5xx*pow(H,2) - 96*G5pp*G5ppx*G5x*G5xx*pow(H,2) - 48*G4p*G5ppxx*G5x*G5xx*pow(H,2) - 30*G2xx*G5px*G5x*G5xx*pow(H,2) - 
   384*G3px*G5px*G5x*G5xx*pow(H,2) + 12*G4ppx*G5px*G5x*G5xx*pow(H,2) + 228*G5ppp*G5px*G5x*G5xx*pow(H,2) + 96*G2x*G5pxx*G5x*G5xx*pow(H,2) - 192*G3p*G5pxx*G5x*G5xx*pow(H,2) + 
   48*G4pp*G5pxx*G5x*G5xx*pow(H,2) - 48*G2pxx*pow(G5x,2)*G5xx*pow(H,2) + 24*G3ppx*pow(G5x,2)*G5xx*pow(H,2) + 48*G4pppx*pow(G5x,2)*G5xx*pow(H,2) + 36*G3pxx*G4*pow(G5xx,2)*pow(H,2) + 
   36*G3xx*G4p*pow(G5xx,2)*pow(H,2) + 756*pow(G4px,2)*pow(G5xx,2)*pow(H,2) + 102*G3px*G4x*pow(G5xx,2)*pow(H,2) - 36*G4ppx*G4x*pow(G5xx,2)*pow(H,2) - 
   684*G4pp*G4xx*pow(G5xx,2)*pow(H,2) - 114*G3px*G5p*pow(G5xx,2)*pow(H,2) - 72*G4ppx*G5p*pow(G5xx,2)*pow(H,2) - 186*G4px*G5pp*pow(G5xx,2)*pow(H,2) - 
   6*pow(G5pp,2)*pow(G5xx,2)*pow(H,2) - 12*G4x*G5ppp*pow(G5xx,2)*pow(H,2) + 114*G5p*G5ppp*pow(G5xx,2)*pow(H,2) - 36*G4*G5pppx*pow(G5xx,2)*pow(H,2) + 
   18*G4p*G5ppx*pow(G5xx,2)*pow(H,2) - 39*G2x*G5px*pow(G5xx,2)*pow(H,2) + 78*G3p*G5px*pow(G5xx,2)*pow(H,2) + 330*G4pp*G5px*pow(G5xx,2)*pow(H,2) - 24*G2px*G5x*pow(G5xx,2)*pow(H,2) + 
   87*G3pp*G5x*pow(G5xx,2)*pow(H,2) - 78*G4ppp*G5x*pow(G5xx,2)*pow(H,2) - 9*G2p*pow(G5xx,3)*pow(H,2) - 576*G4px*G4x*G4xx*G5xxx*pow(H,2) + 288*G4p*pow(G4xx,2)*G5xxx*pow(H,2) + 
   576*G4px*G4xx*G5p*G5xxx*pow(H,2) + 144*G4x*G4xx*G5pp*G5xxx*pow(H,2) - 216*G4xx*G5p*G5pp*G5xxx*pow(H,2) - 304*pow(G4x,2)*G5ppx*G5xxx*pow(H,2) + 192*G4*G4xx*G5ppx*G5xxx*pow(H,2) + 
   400*G4x*G5p*G5ppx*G5xxx*pow(H,2) - 124*pow(G5p,2)*G5ppx*G5xxx*pow(H,2) - 32*G4*G4x*G5ppxx*G5xxx*pow(H,2) + 16*G4*G5p*G5ppxx*G5xxx*pow(H,2) + 288*G4px*G4x*G5px*G5xxx*pow(H,2) - 
   216*G4p*G4xx*G5px*G5xxx*pow(H,2) - 468*G4px*G5p*G5px*G5xxx*pow(H,2) - 24*G4x*G5pp*G5px*G5xxx*pow(H,2) + 120*G5p*G5pp*G5px*G5xxx*pow(H,2) - 120*G4*G5ppx*G5px*G5xxx*pow(H,2) + 
   24*G4p*pow(G5px,2)*G5xxx*pow(H,2) + 120*G4*G4px*G5pxx*G5xxx*pow(H,2) - 24*G4p*G4x*G5pxx*G5xxx*pow(H,2) + 12*G4p*G5p*G5pxx*G5xxx*pow(H,2) - 24*G4*G5pp*G5pxx*G5xxx*pow(H,2) - 
   96*G4*G4ppxx*G5x*G5xxx*pow(H,2) + 264*pow(G4px,2)*G5x*G5xxx*pow(H,2) + 216*G3px*G4x*G5x*G5xxx*pow(H,2) - 672*G4ppx*G4x*G5x*G5xxx*pow(H,2) + 24*G4pp*G4xx*G5x*G5xxx*pow(H,2) - 
   126*G3px*G5p*G5x*G5xxx*pow(H,2) + 420*G4ppx*G5p*G5x*G5xxx*pow(H,2) - 36*G4px*G5pp*G5x*G5xxx*pow(H,2) - 12*pow(G5pp,2)*G5x*G5xxx*pow(H,2) + 72*G4x*G5ppp*G5x*G5xxx*pow(H,2) - 
   60*G5p*G5ppp*G5x*G5xxx*pow(H,2) + 48*G4*G5pppx*G5x*G5xxx*pow(H,2) + 36*G4p*G5ppx*G5x*G5xxx*pow(H,2) + 12*G2px*pow(G5x,2)*G5xxx*pow(H,2) - 42*G3pp*pow(G5x,2)*G5xxx*pow(H,2) + 
   36*G4ppp*pow(G5x,2)*G5xxx*pow(H,2) - 12*G3px*G4*G5xx*G5xxx*pow(H,2) + 72*G4*G4ppx*G5xx*G5xxx*pow(H,2) - 96*G4p*G4px*G5xx*G5xxx*pow(H,2) + 24*G4pp*G4x*G5xx*G5xxx*pow(H,2) - 
   12*G4pp*G5p*G5xx*G5xxx*pow(H,2) - 24*G4*G5ppp*G5xx*G5xxx*pow(H,2) + 6*G2p*G5x*G5xx*G5xxx*pow(H,2) + 954*G4xx*G5pxx*pow(G5x,2)*pow(H,4) - 512*G5px*G5pxx*pow(G5x,2)*pow(H,4) - 
   504*G4x*G5pxxx*pow(G5x,2)*pow(H,4) + 378*G5p*G5pxxx*pow(G5x,2)*pow(H,4) - 540*G4pxxx*pow(G5x,3)*pow(H,4) + 240*G5ppxx*pow(G5x,3)*pow(H,4) + 420*G4xx*G5px*G5x*G5xx*pow(H,4) - 
   210*pow(G5px,2)*G5x*G5xx*pow(H,4) + 876*G4x*G5pxx*G5x*G5xx*pow(H,4) - 858*G5p*G5pxx*G5x*G5xx*pow(H,4) + 120*G4*G5pxxx*G5x*G5xx*pow(H,4) - 385*G5ppx*pow(G5x,2)*G5xx*pow(H,4) + 
   6*G4x*G5px*pow(G5xx,2)*pow(H,4) + 54*G5p*G5px*pow(G5xx,2)*pow(H,4) - 36*G4*G5pxx*pow(G5xx,2)*pow(H,4) - 174*G4px*G5x*pow(G5xx,2)*pow(H,4) - 159*G5pp*G5x*pow(G5xx,2)*pow(H,4) + 
   66*G4p*pow(G5xx,3)*pow(H,4) + 60*G4x*G5px*G5x*G5xxx*pow(H,4) + 30*G5p*G5px*G5x*G5xxx*pow(H,4) - 84*G4*G5pxx*G5x*G5xxx*pow(H,4) - 72*G4px*pow(G5x,2)*G5xxx*pow(H,4) - 
   36*G5pp*pow(G5x,2)*G5xxx*pow(H,4) - 12*G4*G5px*G5xx*G5xxx*pow(H,4) + 108*G4p*G5x*G5xx*G5xxx*pow(H,4) + 
   144*pow(G4pxx,2)*(4*G4x*G4xx - 2*G4xx*G5p - G5x*(G3x - 6*G4px + 2*G5pp + 4*G5x*pow(H,2))) - 
   12*G4pxx*(-288*G3x*pow(G4xx,2) + 720*G4px*pow(G4xx,2) - 72*G3x*G4x*G4xxx + 240*G4px*G4x*G4xxx + 36*G3x*G4xxx*G5p - 120*G4px*G4xxx*G5p + 12*G4pxxx*pow(-2*G4x + G5p,2) - 72*pow(G4xx,2)*G5pp - 
      48*G4x*G4xxx*G5pp + 24*G4xxx*G5p*G5pp + 96*G4x*G4xx*G5ppx - 48*G4xx*G5p*G5ppx - 16*pow(G4x,2)*G5ppxx + 16*G4x*G5p*G5ppxx - 4*pow(G5p,2)*G5ppxx + 282*G3x*G4xx*G5px - 636*G4px*G4xx*G5px + 
      36*G4xx*G5pp*G5px - 36*G4x*G5ppx*G5px + 18*G5p*G5ppx*G5px - 66*G3x*pow(G5px,2) + 132*G4px*pow(G5px,2) + 34*G3x*G4x*G5pxx - 140*G4px*G4x*G5pxx - 17*G3x*G5p*G5pxx + 70*G4px*G5p*G5pxx + 
      36*G4x*G5pp*G5pxx - 18*G5p*G5pp*G5pxx + 24*G3pxx*G4x*G5x - 12*G2xx*G4xx*G5x - 42*G3px*G4xx*G5x + 60*G4ppx*G4xx*G5x - 12*G2x*G4xxx*G5x + 24*G3p*G4xxx*G5x - 24*G4pp*G4xxx*G5x - 12*G3pxx*G5p*G5x + 
      24*G4xx*G5ppp*G5x - 24*G4x*G5pppx*G5x + 12*G5p*G5pppx*G5x + 3*G3x*G5ppx*G5x + 18*G4px*G5ppx*G5x - 12*G5pp*G5ppx*G5x + 8*G2xx*G5px*G5x + 16*G3px*G5px*G5x - 24*G4ppx*G5px*G5x - 12*G5ppp*G5px*G5x + 
      9*G2x*G5pxx*G5x - 18*G3p*G5pxx*G5x + 18*G4pp*G5pxx*G5x + 2*G2pxx*pow(G5x,2) + 4*G3ppx*pow(G5x,2) - 12*G4pppx*pow(G5x,2) - 21*pow(G3x,2)*G5xx + 96*G3x*G4px*G5xx - 84*pow(G4px,2)*G5xx - 
      4*G2xx*G4x*G5xx - 26*G3px*G4x*G5xx + 36*G4ppx*G4x*G5xx - 24*G2x*G4xx*G5xx + 48*G3p*G4xx*G5xx - 48*G4pp*G4xx*G5xx + 2*G2xx*G5p*G5xx + 13*G3px*G5p*G5xx - 18*G4ppx*G5p*G5xx - 6*G3x*G5pp*G5xx - 
      12*G4px*G5pp*G5xx + 6*pow(G5pp,2)*G5xx + 12*G4x*G5ppp*G5xx - 6*G5p*G5ppp*G5xx + 9*G2x*G5px*G5xx - 18*G3p*G5px*G5xx + 18*G4pp*G5px*G5xx - 4*G2px*G5x*G5xx + G3pp*G5x*G5xx + 6*G4ppp*G5x*G5xx - 
      4*G2x*G4x*G5xxx + 8*G3p*G4x*G5xxx - 8*G4pp*G4x*G5xxx + 2*G2x*G5p*G5xxx - 4*G3p*G5p*G5xxx + 4*G4pp*G5p*G5xxx - 720*pow(G4xx,2)*G5x*pow(H,2) - 288*G4x*G4xxx*G5x*pow(H,2) + 
      180*G4xxx*G5p*G5x*pow(H,2) + 736*G4xx*G5px*G5x*pow(H,2) - 194*pow(G5px,2)*G5x*pow(H,2) + 192*G4x*G5pxx*G5x*pow(H,2) - 103*G5p*G5pxx*G5x*pow(H,2) - 8*G4*G5pxxx*G5x*pow(H,2) - 
      23*G5ppx*pow(G5x,2)*pow(H,2) - 504*G4x*G4xx*G5xx*pow(H,2) + 24*G4*G4xxx*G5xx*pow(H,2) + 444*G4xx*G5p*G5xx*pow(H,2) + 240*G4x*G5px*G5xx*pow(H,2) - 225*G5p*G5px*G5xx*pow(H,2) + 
      2*G4*G5pxx*G5xx*pow(H,2) - 126*G3x*G5x*G5xx*pow(H,2) + 314*G4px*G5x*G5xx*pow(H,2) - 7*G5pp*G5x*G5xx*pow(H,2) + 18*G4p*pow(G5xx,2)*pow(H,2) - 48*pow(G4x,2)*G5xxx*pow(H,2) + 
      24*G4*G4xx*G5xxx*pow(H,2) + 60*G4x*G5p*G5xxx*pow(H,2) - 18*pow(G5p,2)*G5xxx*pow(H,2) - 16*G4*G5px*G5xxx*pow(H,2) + 4*G4p*G5x*G5xxx*pow(H,2) - 99*pow(G5x,2)*G5xx*pow(H,4) - 
      6*G3xx*(12*G4x*G4xx - 6*G4xx*G5p - 8*G4x*G5px + 4*G5p*G5px + 3*G3x*G5x - 10*G4px*G5x + 2*G5pp*G5x + 6*pow(G5x,2)*pow(H,2))) - 
   3*G3x*(-216*G4xx*G4xxx*G5pp + 432*pow(G4xx,2)*G5ppx + 168*G4x*G4xxx*G5ppx - 84*G4xxx*G5p*G5ppx - 184*G4x*G4xx*G5ppxx + 92*G4xx*G5p*G5ppxx + 108*G3xx*G4xx*G5px + 144*G4xxx*G5pp*G5px - 
      528*G4xx*G5ppx*G5px + 104*G4x*G5ppxx*G5px - 52*G5p*G5ppxx*G5px - 60*G3xx*pow(G5px,2) + 156*G5ppx*pow(G5px,2) - 36*G3xx*G4x*G5pxx + 18*G3xx*G5p*G5pxx + 44*G4xx*G5pp*G5pxx - 44*G4x*G5ppx*G5pxx + 
      22*G5p*G5ppx*G5pxx - 40*G5pp*G5px*G5pxx + 8*G4x*G5pp*G5pxxx - 4*G5p*G5pp*G5pxxx + 108*G3pxx*G4xx*G5x - 120*G4ppxx*G4xx*G5x - 36*G3px*G4xxx*G5x + 120*G4ppx*G4xxx*G5x - 24*G4xxx*G5ppp*G5x - 
      48*G4xx*G5pppx*G5x + 42*G3xx*G5ppx*G5x - 24*pow(G5ppx,2)*G5x - 12*G5pp*G5ppxx*G5x - 60*G3pxx*G5px*G5x + 72*G4ppxx*G5px*G5x + 24*G5pppx*G5px*G5x - 6*G2xx*G5pxx*G5x + 8*G3px*G5pxx*G5x - 
      28*G4ppx*G5pxx*G5x + 12*G5ppp*G5pxx*G5x + 2*G2x*G5pxxx*G5x - 4*G3p*G5pxxx*G5x + 4*G4pp*G5pxxx*G5x + 48*G3pxx*G4x*G5xx - 48*G4ppxx*G4x*G5xx - 132*G3px*G4xx*G5xx - 24*G4ppx*G4xx*G5xx - 
      24*G3pxx*G5p*G5xx + 24*G4ppxx*G5p*G5xx - 24*G3xx*G5pp*G5xx + 144*G4xx*G5ppp*G5xx - 24*G4x*G5pppx*G5xx + 12*G5p*G5pppx*G5xx + 36*G5pp*G5ppx*G5xx + 8*G2xx*G5px*G5xx + 76*G3px*G5px*G5xx - 
      24*G4ppx*G5px*G5xx - 72*G5ppp*G5px*G5xx - 6*G2x*G5pxx*G5xx + 12*G3p*G5pxx*G5xx - 12*G4pp*G5pxx*G5xx + 8*G2pxx*G5x*G5xx + 4*G3ppx*G5x*G5xx - 24*G4pppx*G5x*G5xx - 9*G3pp*pow(G5xx,2) + 
      18*G4ppp*pow(G5xx,2) - 64*pow(G4px,2)*G5xxx - 12*G3px*G4x*G5xxx + 40*G4ppx*G4x*G5xxx + 6*G3px*G5p*G5xxx - 20*G4ppx*G5p*G5xxx - 4*pow(G5pp,2)*G5xxx - 8*G4x*G5ppp*G5xxx + 4*G5p*G5ppp*G5xxx + 
      2*G2x*G5px*G5xxx - 4*G3p*G5px*G5xxx + 4*G4pp*G5px*G5xxx + 2*G3pp*G5x*G5xxx - 4*G4ppp*G5x*G5xxx + 60*G4xxx*G5px*G5x*pow(H,2) - 528*G4xx*G5pxx*G5x*pow(H,2) + 260*G5px*G5pxx*G5x*pow(H,2) + 
      132*G4x*G5pxxx*G5x*pow(H,2) - 84*G5p*G5pxxx*G5x*pow(H,2) - 104*G5ppxx*pow(G5x,2)*pow(H,2) + 288*G4xx*G5px*G5xx*pow(H,2) - 168*pow(G5px,2)*G5xx*pow(H,2) - 
      200*G4x*G5pxx*G5xx*pow(H,2) + 166*G5p*G5pxx*G5xx*pow(H,2) - 16*G4*G5pxxx*G5xx*pow(H,2) + 208*G5ppx*G5x*G5xx*pow(H,2) - 11*G5pp*pow(G5xx,2)*pow(H,2) + 8*G4x*G5px*G5xxx*pow(H,2) - 
      22*G5p*G5px*G5xxx*pow(H,2) + 12*G4*G5pxx*G5xxx*pow(H,2) - 6*G5pp*G5x*G5xxx*pow(H,2) - 16*G4p*G5xx*G5xxx*pow(H,2) + 
      24*G4pxxx*(-9*G4xx*G5p + 2*G4x*(9*G4xx - 5*G5px) + 5*G5p*G5px - 8*G4px*G5x + G5pp*G5x + 9*pow(G5x,2)*pow(H,2)) + 
      G4px*(432*G4xx*G4xxx - 432*G4xxx*G5px + 272*G4xx*G5pxx - 48*G5px*G5pxx - 64*G4x*G5pxxx + 32*G5p*G5pxxx + 88*G5ppxx*G5x + 48*G3xx*G5xx - 204*G5ppx*G5xx + 40*G5pp*G5xxx + 
         82*pow(G5xx,2)*pow(H,2) + 36*G5x*G5xxx*pow(H,2))) - 6*G5pxxx*G5x*G5xx*pow(H,2)*rptot + 9*G5pxx*pow(G5xx,2)*pow(H,2)*rptot + 6*G5pxx*G5x*G5xxx*pow(H,2)*rptot - 
   6*G5px*G5xx*G5xxx*pow(H,2)*rptot) - pow(a,9)*pow(dphi,8)*(-9*G3pp*pow(G3x,3) - 18*G3px*pow(G3x,2)*G4pp + 18*pow(G3x,3)*G4ppp + 36*pow(G3x,2)*G4pp*G4ppx + 12*G2px*pow(G3x,2)*G4px + 
   42*G3pp*pow(G3x,2)*G4px + 24*G2xx*G3x*G4pp*G4px + 48*G3px*G3x*G4pp*G4px - 108*pow(G3x,2)*G4ppp*G4px - 144*G3x*G4pp*G4ppx*G4px - 48*G2px*G3x*pow(G4px,2) - 60*G3pp*G3x*pow(G4px,2) - 
   48*G2xx*G4pp*pow(G4px,2) - 24*G3px*G4pp*pow(G4px,2) + 216*G3x*G4ppp*pow(G4px,2) + 144*G4pp*G4ppx*pow(G4px,2) + 48*G2px*pow(G4px,3) + 24*G3pp*pow(G4px,3) - 144*G4ppp*pow(G4px,3) + 
   12*G2xx*G3pp*G3x*G4x - 24*G2px*G3px*G3x*G4x + 12*G3pp*G3px*G3x*G4x + 12*G2ppx*pow(G3x,2)*G4x - 12*G3ppp*pow(G3x,2)*G4x - 24*G2xx*G3px*G4pp*G4x + 24*pow(G3px,2)*G4pp*G4x + 24*G2pxx*G3x*G4pp*G4x - 
   24*G3ppx*G3x*G4pp*G4x - 24*G2xx*G3x*G4ppp*G4x + 24*G3px*G3x*G4ppp*G4x + 48*G2px*G3x*G4ppx*G4x - 48*G3pp*G3x*G4ppx*G4x + 48*G2xx*G4pp*G4ppx*G4x - 48*G3px*G4pp*G4ppx*G4x - 24*G2xx*G3pp*G4px*G4x + 
   48*G2px*G3px*G4px*G4x - 24*G3pp*G3px*G4px*G4x - 48*G2ppx*G3x*G4px*G4x + 48*G3ppp*G3x*G4px*G4x - 48*G2pxx*G4pp*G4px*G4x + 48*G3ppx*G4pp*G4px*G4x + 48*G2xx*G4ppp*G4px*G4x - 48*G3px*G4ppp*G4px*G4x - 
   96*G2px*G4ppx*G4px*G4x + 96*G3pp*G4ppx*G4px*G4x + 48*G2ppx*pow(G4px,2)*G4x - 48*G3ppp*pow(G4px,2)*G4x + 16*G2px*G2pxx*pow(G4x,2) - 16*G2ppx*G2xx*pow(G4x,2) - 16*G2pxx*G3pp*pow(G4x,2) + 
   16*G2xx*G3ppp*pow(G4x,2) - 16*G2px*G3ppx*pow(G4x,2) + 16*G3pp*G3ppx*pow(G4x,2) + 16*G2ppx*G3px*pow(G4x,2) - 16*G3ppp*G3px*pow(G4x,2) - 6*G2xx*G3pp*G3x*G5p + 12*G2px*G3px*G3x*G5p - 
   6*G3pp*G3px*G3x*G5p - 6*G2ppx*pow(G3x,2)*G5p + 6*G3ppp*pow(G3x,2)*G5p + 12*G2xx*G3px*G4pp*G5p - 12*pow(G3px,2)*G4pp*G5p - 12*G2pxx*G3x*G4pp*G5p + 12*G3ppx*G3x*G4pp*G5p + 12*G2xx*G3x*G4ppp*G5p - 
   12*G3px*G3x*G4ppp*G5p - 24*G2px*G3x*G4ppx*G5p + 24*G3pp*G3x*G4ppx*G5p - 24*G2xx*G4pp*G4ppx*G5p + 24*G3px*G4pp*G4ppx*G5p + 12*G2xx*G3pp*G4px*G5p - 24*G2px*G3px*G4px*G5p + 12*G3pp*G3px*G4px*G5p + 
   24*G2ppx*G3x*G4px*G5p - 24*G3ppp*G3x*G4px*G5p + 24*G2pxx*G4pp*G4px*G5p - 24*G3ppx*G4pp*G4px*G5p - 24*G2xx*G4ppp*G4px*G5p + 24*G3px*G4ppp*G4px*G5p + 48*G2px*G4ppx*G4px*G5p - 48*G3pp*G4ppx*G4px*G5p - 
   24*G2ppx*pow(G4px,2)*G5p + 24*G3ppp*pow(G4px,2)*G5p - 16*G2px*G2pxx*G4x*G5p + 16*G2ppx*G2xx*G4x*G5p + 16*G2pxx*G3pp*G4x*G5p - 16*G2xx*G3ppp*G4x*G5p + 16*G2px*G3ppx*G4x*G5p - 
   16*G3pp*G3ppx*G4x*G5p - 16*G2ppx*G3px*G4x*G5p + 16*G3ppp*G3px*G4x*G5p + 4*G2px*G2pxx*pow(G5p,2) - 4*G2ppx*G2xx*pow(G5p,2) - 4*G2pxx*G3pp*pow(G5p,2) + 4*G2xx*G3ppp*pow(G5p,2) - 
   4*G2px*G3ppx*pow(G5p,2) + 4*G3pp*G3ppx*pow(G5p,2) + 4*G2ppx*G3px*pow(G5p,2) - 4*G3ppp*G3px*pow(G5p,2) - 6*G2px*pow(G3x,2)*G5pp + 6*G3pp*pow(G3x,2)*G5pp - 12*G2xx*G3x*G4pp*G5pp + 
   12*G3px*G3x*G4pp*G5pp + 24*G2px*G3x*G4px*G5pp - 24*G3pp*G3x*G4px*G5pp + 24*G2xx*G4pp*G4px*G5pp - 24*G3px*G4pp*G4px*G5pp - 24*G2px*pow(G4px,2)*G5pp + 24*G3pp*pow(G4px,2)*G5pp - 
   108*G3pxx*pow(G3x,2)*G4*pow(H,2) + 108*G3px*G3x*G3xx*G4*pow(H,2) - 108*pow(G3x,2)*G3xx*G4p*pow(H,2) - 360*G3x*G3xx*G4*G4ppx*pow(H,2) + 144*pow(G3x,2)*G4*G4ppxx*pow(H,2) + 
   162*pow(G3x,3)*G4px*pow(H,2) + 576*G3pxx*G3x*G4*G4px*pow(H,2) - 360*G3px*G3xx*G4*G4px*pow(H,2) + 288*G3x*G3xx*G4p*G4px*pow(H,2) + 1008*G3xx*G4*G4ppx*G4px*pow(H,2) - 
   864*G3x*G4*G4ppxx*G4px*pow(H,2) - 1332*pow(G3x,2)*pow(G4px,2)*pow(H,2) - 720*G3pxx*G4*pow(G4px,2)*pow(H,2) + 144*G3xx*G4p*pow(G4px,2)*pow(H,2) + 
   1152*G4*G4ppxx*pow(G4px,2)*pow(H,2) + 3672*G3x*pow(G4px,3)*pow(H,2) - 3312*pow(G4px,4)*pow(H,2) + 144*G2xx*G3x*G4*G4pxx*pow(H,2) + 72*G3px*G3x*G4*G4pxx*pow(H,2) + 
   648*pow(G3x,2)*G4p*G4pxx*pow(H,2) + 288*G3xx*G4*G4pp*G4pxx*pow(H,2) - 144*G3x*G4*G4ppx*G4pxx*pow(H,2) - 480*G2xx*G4*G4px*G4pxx*pow(H,2) + 336*G3px*G4*G4px*G4pxx*pow(H,2) - 
   3168*G3x*G4p*G4px*G4pxx*pow(H,2) - 288*G4*G4ppx*G4px*G4pxx*pow(H,2) + 3168*G4p*pow(G4px,2)*G4pxx*pow(H,2) - 576*G4*G4pp*pow(G4pxx,2)*pow(H,2) - 144*G3x*G4*G4pp*G4pxxx*pow(H,2) + 
   288*G4*G4pp*G4px*G4pxxx*pow(H,2) - 378*G3px*pow(G3x,2)*G4x*pow(H,2) - 288*G3px*G3pxx*G4*G4x*pow(H,2) + 288*G3ppx*G3xx*G4*G4x*pow(H,2) - 288*G3pxx*G3x*G4p*G4x*pow(H,2) + 
   72*G3px*G3xx*G4p*G4x*pow(H,2) + 72*G3x*G3xx*G4pp*G4x*pow(H,2) - 576*G3xx*G4*G4pppx*G4x*pow(H,2) + 396*pow(G3x,2)*G4ppx*G4x*pow(H,2) + 576*G3pxx*G4*G4ppx*G4x*pow(H,2) - 
   432*G3xx*G4p*G4ppx*G4x*pow(H,2) + 192*G2xx*G4*G4ppxx*G4x*pow(H,2) + 384*G3px*G4*G4ppxx*G4x*pow(H,2) + 288*G3x*G4p*G4ppxx*G4x*pow(H,2) - 1152*G4*G4ppx*G4ppxx*G4x*pow(H,2) + 
   144*G2xx*G3x*G4px*G4x*pow(H,2) + 1872*G3px*G3x*G4px*G4x*pow(H,2) + 864*G3pxx*G4p*G4px*G4x*pow(H,2) - 432*G3xx*G4pp*G4px*G4x*pow(H,2) - 2592*G3x*G4ppx*G4px*G4x*pow(H,2) - 
   1152*G4p*G4ppxx*G4px*G4x*pow(H,2) - 288*G2xx*pow(G4px,2)*G4x*pow(H,2) - 2520*G3px*pow(G4px,2)*G4x*pow(H,2) + 4176*G4ppx*pow(G4px,2)*G4x*pow(H,2) - 192*G2pxx*G4*G4pxx*G4x*pow(H,2) - 
   384*G3ppx*G4*G4pxx*G4x*pow(H,2) + 96*G2xx*G4p*G4pxx*G4x*pow(H,2) + 624*G3px*G4p*G4pxx*G4x*pow(H,2) - 432*G3x*G4pp*G4pxx*G4x*pow(H,2) + 1152*G4*G4pppx*G4pxx*G4x*pow(H,2) - 
   864*G4p*G4ppx*G4pxx*G4x*pow(H,2) + 2016*G4pp*G4px*G4pxx*G4x*pow(H,2) - 192*G2px*G4*G4pxxx*G4x*pow(H,2) + 192*G3pp*G4*G4pxxx*G4x*pow(H,2) - 288*G4p*G4pp*G4pxxx*G4x*pow(H,2) - 
   216*G2xx*G3px*pow(G4x,2)*pow(H,2) + 216*pow(G3px,2)*pow(G4x,2)*pow(H,2) + 360*G2pxx*G3x*pow(G4x,2)*pow(H,2) - 360*G3ppx*G3x*pow(G4x,2)*pow(H,2) - 
   144*G2px*G3xx*pow(G4x,2)*pow(H,2) + 432*G3pp*G3xx*pow(G4x,2)*pow(H,2) + 288*G3pxx*G4pp*pow(G4x,2)*pow(H,2) - 288*G3xx*G4ppp*pow(G4x,2)*pow(H,2) + 
   624*G2xx*G4ppx*pow(G4x,2)*pow(H,2) - 912*G3px*G4ppx*pow(G4x,2)*pow(H,2) + 576*pow(G4ppx,2)*pow(G4x,2)*pow(H,2) - 576*G4pp*G4ppxx*pow(G4x,2)*pow(H,2) - 
   912*G2pxx*G4px*pow(G4x,2)*pow(H,2) + 1200*G3ppx*G4px*pow(G4x,2)*pow(H,2) - 576*G4pppx*G4px*pow(G4x,2)*pow(H,2) + 384*G2px*G4pxx*pow(G4x,2)*pow(H,2) - 
   1056*G3pp*G4pxx*pow(G4x,2)*pow(H,2) + 576*G4ppp*G4pxx*pow(G4x,2)*pow(H,2) - 96*G2p*G4pxxx*pow(G4x,2)*pow(H,2) + 72*G2xx*G3px*G4*G4xx*pow(H,2) + 216*pow(G3px,2)*G4*G4xx*pow(H,2) - 
   216*G2pxx*G3x*G4*G4xx*pow(H,2) - 72*G3ppx*G3x*G4*G4xx*pow(H,2) - 144*G3pp*G3xx*G4*G4xx*pow(H,2) - 216*G2xx*G3x*G4p*G4xx*pow(H,2) + 720*G3px*G3x*G4p*G4xx*pow(H,2) + 
   972*pow(G3x,2)*G4pp*G4xx*pow(H,2) - 288*G3pxx*G4*G4pp*G4xx*pow(H,2) - 288*G3xx*G4p*G4pp*G4xx*pow(H,2) + 288*G3xx*G4*G4ppp*G4xx*pow(H,2) + 576*G3x*G4*G4pppx*G4xx*pow(H,2) - 
   336*G2xx*G4*G4ppx*G4xx*pow(H,2) - 816*G3px*G4*G4ppx*G4xx*pow(H,2) + 720*G3x*G4p*G4ppx*G4xx*pow(H,2) + 1152*G4*pow(G4ppx,2)*G4xx*pow(H,2) + 576*G4*G4pp*G4ppxx*G4xx*pow(H,2) + 
   624*G2pxx*G4*G4px*G4xx*pow(H,2) - 48*G3ppx*G4*G4px*G4xx*pow(H,2) + 384*G2xx*G4p*G4px*G4xx*pow(H,2) - 2544*G3px*G4p*G4px*G4xx*pow(H,2) - 3888*G3x*G4pp*G4px*G4xx*pow(H,2) - 
   1152*G4*G4pppx*G4px*G4xx*pow(H,2) + 864*G4p*G4ppx*G4px*G4xx*pow(H,2) + 3312*G4pp*pow(G4px,2)*G4xx*pow(H,2) + 384*G2px*G4*G4pxx*G4xx*pow(H,2) - 96*G3pp*G4*G4pxx*G4xx*pow(H,2) + 
   1728*G4p*G4pp*G4pxx*G4xx*pow(H,2) - 576*G4*G4ppp*G4pxx*G4xx*pow(H,2) + 504*G2px*G3x*G4x*G4xx*pow(H,2) - 1296*G3pp*G3x*G4x*G4xx*pow(H,2) - 144*G2p*G3xx*G4x*G4xx*pow(H,2) - 
   240*G2pxx*G4p*G4x*G4xx*pow(H,2) - 336*G3ppx*G4p*G4x*G4xx*pow(H,2) - 48*G2xx*G4pp*G4x*G4xx*pow(H,2) + 48*G3px*G4pp*G4x*G4xx*pow(H,2) + 576*G3x*G4ppp*G4x*G4xx*pow(H,2) + 
   1152*G4p*G4pppx*G4x*G4xx*pow(H,2) - 576*G4pp*G4ppx*G4x*G4xx*pow(H,2) - 1200*G2px*G4px*G4x*G4xx*pow(H,2) + 2688*G3pp*G4px*G4x*G4xx*pow(H,2) - 576*G4ppp*G4px*G4x*G4xx*pow(H,2) + 
   672*G2p*G4pxx*G4x*G4xx*pow(H,2) + 432*G2p*G3x*pow(G4xx,2)*pow(H,2) - 192*G2ppx*G4*pow(G4xx,2)*pow(H,2) + 192*G3ppp*G4*pow(G4xx,2)*pow(H,2) - 192*G2px*G4p*pow(G4xx,2)*pow(H,2) + 
   1056*G3pp*G4p*pow(G4xx,2)*pow(H,2) + 576*pow(G4pp,2)*pow(G4xx,2)*pow(H,2) - 1728*G4p*G4ppp*pow(G4xx,2)*pow(H,2) - 1056*G2p*G4px*pow(G4xx,2)*pow(H,2) - 
   192*G2pp*G4x*pow(G4xx,2)*pow(H,2) - 72*G3pp*G3x*G4*G4xxx*pow(H,2) + 144*G3px*G4*G4pp*G4xxx*pow(H,2) - 144*G3x*G4p*G4pp*G4xxx*pow(H,2) + 144*G3x*G4*G4ppp*G4xxx*pow(H,2) - 
   288*G4*G4pp*G4ppx*G4xxx*pow(H,2) + 144*G3pp*G4*G4px*G4xxx*pow(H,2) - 288*G4*G4ppp*G4px*G4xxx*pow(H,2) - 72*G2p*G3x*G4x*G4xxx*pow(H,2) + 192*G2ppx*G4*G4x*G4xxx*pow(H,2) - 
   192*G3ppp*G4*G4x*G4xxx*pow(H,2) - 144*G3pp*G4p*G4x*G4xxx*pow(H,2) + 288*pow(G4pp,2)*G4x*G4xxx*pow(H,2) + 288*G4p*G4ppp*G4x*G4xxx*pow(H,2) + 144*G2p*G4px*G4x*G4xxx*pow(H,2) + 
   96*G2pp*pow(G4x,2)*G4xxx*pow(H,2) + 324*G3px*pow(G3x,2)*G5p*pow(H,2) + 144*G3px*G3pxx*G4*G5p*pow(H,2) - 144*G3ppx*G3xx*G4*G5p*pow(H,2) + 144*G3pxx*G3x*G4p*G5p*pow(H,2) - 
   36*G3px*G3xx*G4p*G5p*pow(H,2) - 36*G3x*G3xx*G4pp*G5p*pow(H,2) + 288*G3xx*G4*G4pppx*G5p*pow(H,2) - 252*pow(G3x,2)*G4ppx*G5p*pow(H,2) - 288*G3pxx*G4*G4ppx*G5p*pow(H,2) + 
   216*G3xx*G4p*G4ppx*G5p*pow(H,2) - 96*G2xx*G4*G4ppxx*G5p*pow(H,2) - 192*G3px*G4*G4ppxx*G5p*pow(H,2) - 144*G3x*G4p*G4ppxx*G5p*pow(H,2) + 576*G4*G4ppx*G4ppxx*G5p*pow(H,2) - 
   180*G2xx*G3x*G4px*G5p*pow(H,2) - 1656*G3px*G3x*G4px*G5p*pow(H,2) - 432*G3pxx*G4p*G4px*G5p*pow(H,2) + 72*G3xx*G4pp*G4px*G5p*pow(H,2) + 2088*G3x*G4ppx*G4px*G5p*pow(H,2) + 
   576*G4p*G4ppxx*G4px*G5p*pow(H,2) + 456*G2xx*pow(G4px,2)*G5p*pow(H,2) + 2064*G3px*pow(G4px,2)*G5p*pow(H,2) - 3456*G4ppx*pow(G4px,2)*G5p*pow(H,2) + 96*G2pxx*G4*G4pxx*G5p*pow(H,2) + 
   192*G3ppx*G4*G4pxx*G5p*pow(H,2) - 48*G2xx*G4p*G4pxx*G5p*pow(H,2) - 312*G3px*G4p*G4pxx*G5p*pow(H,2) + 504*G3x*G4pp*G4pxx*G5p*pow(H,2) - 576*G4*G4pppx*G4pxx*G5p*pow(H,2) + 
   432*G4p*G4ppx*G4pxx*G5p*pow(H,2) - 1296*G4pp*G4px*G4pxx*G5p*pow(H,2) + 96*G2px*G4*G4pxxx*G5p*pow(H,2) - 96*G3pp*G4*G4pxxx*G5p*pow(H,2) + 144*G4p*G4pp*G4pxxx*G5p*pow(H,2) + 
   252*G2xx*G3px*G4x*G5p*pow(H,2) - 108*pow(G3px,2)*G4x*G5p*pow(H,2) - 468*G2pxx*G3x*G4x*G5p*pow(H,2) + 324*G3ppx*G3x*G4x*G5p*pow(H,2) + 144*G2px*G3xx*G4x*G5p*pow(H,2) - 
   504*G3pp*G3xx*G4x*G5p*pow(H,2) - 432*G3pxx*G4pp*G4x*G5p*pow(H,2) + 432*G3xx*G4ppp*G4x*G5p*pow(H,2) + 288*G3x*G4pppx*G4x*G5p*pow(H,2) - 792*G2xx*G4ppx*G4x*G5p*pow(H,2) + 
   504*G3px*G4ppx*G4x*G5p*pow(H,2) + 864*G4pp*G4ppxx*G4x*G5p*pow(H,2) + 1224*G2pxx*G4px*G4x*G5p*pow(H,2) - 1224*G3ppx*G4px*G4x*G5p*pow(H,2) - 192*G2px*G4pxx*G4x*G5p*pow(H,2) + 
   1008*G3pp*G4pxx*G4x*G5p*pow(H,2) - 864*G4ppp*G4pxx*G4x*G5p*pow(H,2) + 96*G2p*G4pxxx*G4x*G5p*pow(H,2) - 252*G2px*G3x*G4xx*G5p*pow(H,2) + 1080*G3pp*G3x*G4xx*G5p*pow(H,2) + 
   72*G2p*G3xx*G4xx*G5p*pow(H,2) + 120*G2pxx*G4p*G4xx*G5p*pow(H,2) + 168*G3ppx*G4p*G4xx*G5p*pow(H,2) + 24*G2xx*G4pp*G4xx*G5p*pow(H,2) + 264*G3px*G4pp*G4xx*G5p*pow(H,2) - 
   1152*G3x*G4ppp*G4xx*G5p*pow(H,2) - 576*G4p*G4pppx*G4xx*G5p*pow(H,2) - 288*G4pp*G4ppx*G4xx*G5p*pow(H,2) + 408*G2px*G4px*G4xx*G5p*pow(H,2) - 2016*G3pp*G4px*G4xx*G5p*pow(H,2) + 
   2016*G4ppp*G4px*G4xx*G5p*pow(H,2) - 336*G2p*G4pxx*G4xx*G5p*pow(H,2) - 192*G2ppx*G4x*G4xx*G5p*pow(H,2) + 192*G3ppp*G4x*G4xx*G5p*pow(H,2) + 96*G2pp*pow(G4xx,2)*G5p*pow(H,2) + 
   36*G2p*G3x*G4xxx*G5p*pow(H,2) - 96*G2ppx*G4*G4xxx*G5p*pow(H,2) + 96*G3ppp*G4*G4xxx*G5p*pow(H,2) + 72*G3pp*G4p*G4xxx*G5p*pow(H,2) - 144*pow(G4pp,2)*G4xxx*G5p*pow(H,2) - 
   144*G4p*G4ppp*G4xxx*G5p*pow(H,2) - 72*G2p*G4px*G4xxx*G5p*pow(H,2) - 96*G2pp*G4x*G4xxx*G5p*pow(H,2) - 72*G2xx*G3px*pow(G5p,2)*pow(H,2) + 144*G2pxx*G3x*pow(G5p,2)*pow(H,2) - 
   72*G3ppx*G3x*pow(G5p,2)*pow(H,2) - 36*G2px*G3xx*pow(G5p,2)*pow(H,2) + 144*G3pp*G3xx*pow(G5p,2)*pow(H,2) + 144*G3pxx*G4pp*pow(G5p,2)*pow(H,2) - 144*G3xx*G4ppp*pow(G5p,2)*pow(H,2) - 
   144*G3x*G4pppx*pow(G5p,2)*pow(H,2) + 240*G2xx*G4ppx*pow(G5p,2)*pow(H,2) - 24*G3px*G4ppx*pow(G5p,2)*pow(H,2) - 144*pow(G4ppx,2)*pow(G5p,2)*pow(H,2) - 
   288*G4pp*G4ppxx*pow(G5p,2)*pow(H,2) - 384*G2pxx*G4px*pow(G5p,2)*pow(H,2) + 312*G3ppx*G4px*pow(G5p,2)*pow(H,2) + 144*G4pppx*G4px*pow(G5p,2)*pow(H,2) - 
   240*G3pp*G4pxx*pow(G5p,2)*pow(H,2) + 288*G4ppp*G4pxx*pow(G5p,2)*pow(H,2) - 24*G2p*G4pxxx*pow(G5p,2)*pow(H,2) + 96*G2ppx*G4xx*pow(G5p,2)*pow(H,2) - 
   96*G3ppp*G4xx*pow(G5p,2)*pow(H,2) + 24*G2pp*G4xxx*pow(G5p,2)*pow(H,2) - 27*pow(G3x,3)*G5pp*pow(H,2) - 72*G3pxx*G3x*G4*G5pp*pow(H,2) + 72*G3px*G3xx*G4*G5pp*pow(H,2) + 
   72*G3x*G3xx*G4p*G5pp*pow(H,2) - 144*G3xx*G4*G4ppx*G5pp*pow(H,2) + 144*G3x*G4*G4ppxx*G5pp*pow(H,2) + 306*pow(G3x,2)*G4px*G5pp*pow(H,2) + 144*G3pxx*G4*G4px*G5pp*pow(H,2) - 
   432*G3xx*G4p*G4px*G5pp*pow(H,2) - 288*G4*G4ppxx*G4px*G5pp*pow(H,2) - 1188*G3x*pow(G4px,2)*G5pp*pow(H,2) + 1368*pow(G4px,3)*G5pp*pow(H,2) + 96*G2xx*G4*G4pxx*G5pp*pow(H,2) - 
   240*G3px*G4*G4pxx*G5pp*pow(H,2) + 288*G3x*G4p*G4pxx*G5pp*pow(H,2) + 288*G4*G4ppx*G4pxx*G5pp*pow(H,2) + 72*G3px*G3x*G4x*G5pp*pow(H,2) - 144*G3pxx*G4p*G4x*G5pp*pow(H,2) + 
   288*G3xx*G4pp*G4x*G5pp*pow(H,2) - 288*G3x*G4ppx*G4x*G5pp*pow(H,2) + 288*G4p*G4ppxx*G4x*G5pp*pow(H,2) - 96*G2xx*G4px*G4x*G5pp*pow(H,2) + 96*G3px*G4px*G4x*G5pp*pow(H,2) + 
   288*G4ppx*G4px*G4x*G5pp*pow(H,2) - 576*G4pp*G4pxx*G4x*G5pp*pow(H,2) + 48*G2pxx*pow(G4x,2)*G5pp*pow(H,2) - 48*G3ppx*pow(G4x,2)*G5pp*pow(H,2) - 96*G2pxx*G4*G4xx*G5pp*pow(H,2) + 
   96*G3ppx*G4*G4xx*G5pp*pow(H,2) + 24*G2xx*G4p*G4xx*G5pp*pow(H,2) + 552*G3px*G4p*G4xx*G5pp*pow(H,2) + 288*G3x*G4pp*G4xx*G5pp*pow(H,2) - 1152*G4p*G4ppx*G4xx*G5pp*pow(H,2) - 
   288*G4pp*G4px*G4xx*G5pp*pow(H,2) + 288*G2px*G4x*G4xx*G5pp*pow(H,2) - 384*G3pp*G4x*G4xx*G5pp*pow(H,2) + 96*G2p*pow(G4xx,2)*G5pp*pow(H,2) + 144*G4p*G4pp*G4xxx*G5pp*pow(H,2) + 
   54*G2xx*G3x*G5p*G5pp*pow(H,2) + 54*G3px*G3x*G5p*G5pp*pow(H,2) + 72*G3pxx*G4p*G5p*G5pp*pow(H,2) - 72*G3xx*G4pp*G5p*G5pp*pow(H,2) - 144*G3x*G4ppx*G5p*G5pp*pow(H,2) - 
   144*G4p*G4ppxx*G5p*G5pp*pow(H,2) - 156*G2xx*G4px*G5p*G5pp*pow(H,2) - 132*G3px*G4px*G5p*G5pp*pow(H,2) + 432*G4ppx*G4px*G5p*G5pp*pow(H,2) + 144*G4pp*G4pxx*G5p*G5pp*pow(H,2) - 
   96*G2pxx*G4x*G5p*G5pp*pow(H,2) + 96*G3ppx*G4x*G5p*G5pp*pow(H,2) - 48*G2px*G4xx*G5p*G5pp*pow(H,2) + 96*G3pp*G4xx*G5p*G5pp*pow(H,2) + 36*G2pxx*pow(G5p,2)*G5pp*pow(H,2) - 
   36*G3ppx*pow(G5p,2)*G5pp*pow(H,2) + 18*pow(G3x,2)*pow(G5pp,2)*pow(H,2) + 72*G3xx*G4p*pow(G5pp,2)*pow(H,2) - 72*pow(G4px,2)*pow(G5pp,2)*pow(H,2) - 
   144*G4p*G4pxx*pow(G5pp,2)*pow(H,2) + 48*G2xx*G4x*pow(G5pp,2)*pow(H,2) - 48*G3px*G4x*pow(G5pp,2)*pow(H,2) + 72*G3x*G3xx*G4*G5ppp*pow(H,2) - 144*G3xx*G4*G4px*G5ppp*pow(H,2) - 
   144*G3x*G4*G4pxx*G5ppp*pow(H,2) + 288*G4*G4px*G4pxx*G5ppp*pow(H,2) + 36*pow(G3x,2)*G4x*G5ppp*pow(H,2) + 144*G3xx*G4p*G4x*G5ppp*pow(H,2) - 144*pow(G4px,2)*G4x*G5ppp*pow(H,2) - 
   288*G4p*G4pxx*G4x*G5ppp*pow(H,2) - 48*G2xx*pow(G4x,2)*G5ppp*pow(H,2) + 48*G3px*pow(G4x,2)*G5ppp*pow(H,2) + 96*G2xx*G4*G4xx*G5ppp*pow(H,2) - 96*G3px*G4*G4xx*G5ppp*pow(H,2) - 
   864*G3x*G4p*G4xx*G5ppp*pow(H,2) + 1728*G4p*G4px*G4xx*G5ppp*pow(H,2) - 126*pow(G3x,2)*G5p*G5ppp*pow(H,2) - 72*G3xx*G4p*G5p*G5ppp*pow(H,2) + 432*G3x*G4px*G5p*G5ppp*pow(H,2) - 
   360*pow(G4px,2)*G5p*G5ppp*pow(H,2) + 144*G4p*G4pxx*G5p*G5ppp*pow(H,2) + 96*G2xx*G4x*G5p*G5ppp*pow(H,2) - 96*G3px*G4x*G5p*G5ppp*pow(H,2) - 36*G2xx*pow(G5p,2)*G5ppp*pow(H,2) + 
   36*G3px*pow(G5p,2)*G5ppp*pow(H,2) + 36*pow(G3x,2)*G4*G5pppx*pow(H,2) - 144*G3x*G4*G4px*G5pppx*pow(H,2) + 144*G4*pow(G4px,2)*G5pppx*pow(H,2) - 96*G2xx*G4*G4x*G5pppx*pow(H,2) + 
   96*G3px*G4*G4x*G5pppx*pow(H,2) + 144*G3x*G4p*G4x*G5pppx*pow(H,2) - 288*G4p*G4px*G4x*G5pppx*pow(H,2) + 48*G2xx*G4*G5p*G5pppx*pow(H,2) - 48*G3px*G4*G5p*G5pppx*pow(H,2) - 
   72*G3x*G4p*G5p*G5pppx*pow(H,2) + 144*G4p*G4px*G5p*G5pppx*pow(H,2) - 84*G2xx*G3x*G4*G5ppx*pow(H,2) - 132*G3px*G3x*G4*G5ppx*pow(H,2) - 162*pow(G3x,2)*G4p*G5ppx*pow(H,2) - 
   144*G3xx*G4*G4pp*G5ppx*pow(H,2) + 432*G3x*G4*G4ppx*G5ppx*pow(H,2) + 264*G2xx*G4*G4px*G5ppx*pow(H,2) + 168*G3px*G4*G4px*G5ppx*pow(H,2) + 1080*G3x*G4p*G4px*G5ppx*pow(H,2) - 
   864*G4*G4ppx*G4px*G5ppx*pow(H,2) - 1512*G4p*pow(G4px,2)*G5ppx*pow(H,2) + 288*G4*G4pp*G4pxx*G5ppx*pow(H,2) + 96*G2pxx*G4*G4x*G5ppx*pow(H,2) - 96*G3ppx*G4*G4x*G5ppx*pow(H,2) - 
   72*G2xx*G4p*G4x*G5ppx*pow(H,2) - 360*G3px*G4p*G4x*G5ppx*pow(H,2) + 72*G3x*G4pp*G4x*G5ppx*pow(H,2) + 864*G4p*G4ppx*G4x*G5ppx*pow(H,2) - 432*G4pp*G4px*G4x*G5ppx*pow(H,2) - 
   96*G2px*pow(G4x,2)*G5ppx*pow(H,2) + 144*G3pp*pow(G4x,2)*G5ppx*pow(H,2) - 192*G2px*G4*G4xx*G5ppx*pow(H,2) + 192*G3pp*G4*G4xx*G5ppx*pow(H,2) - 576*G4p*G4pp*G4xx*G5ppx*pow(H,2) - 
   192*G2p*G4x*G4xx*G5ppx*pow(H,2) - 48*G2pxx*G4*G5p*G5ppx*pow(H,2) + 48*G3ppx*G4*G5p*G5ppx*pow(H,2) + 36*G2xx*G4p*G5p*G5ppx*pow(H,2) + 180*G3px*G4p*G5p*G5ppx*pow(H,2) - 
   180*G3x*G4pp*G5p*G5ppx*pow(H,2) - 432*G4p*G4ppx*G5p*G5ppx*pow(H,2) + 504*G4pp*G4px*G5p*G5ppx*pow(H,2) - 48*G3pp*G4x*G5p*G5ppx*pow(H,2) + 96*G2p*G4xx*G5p*G5ppx*pow(H,2) + 
   24*G2px*pow(G5p,2)*G5ppx*pow(H,2) - 12*G3pp*pow(G5p,2)*G5ppx*pow(H,2) - 48*G2xx*G4*G5pp*G5ppx*pow(H,2) + 48*G3px*G4*G5pp*G5ppx*pow(H,2) - 216*G3x*G4p*G5pp*G5ppx*pow(H,2) + 
   432*G4p*G4px*G5pp*G5ppx*pow(H,2) + 72*G3x*G4*G4pp*G5ppxx*pow(H,2) - 144*G4*G4pp*G4px*G5ppxx*pow(H,2) + 96*G2px*G4*G4x*G5ppxx*pow(H,2) - 96*G3pp*G4*G4x*G5ppxx*pow(H,2) + 
   144*G4p*G4pp*G4x*G5ppxx*pow(H,2) + 48*G2p*pow(G4x,2)*G5ppxx*pow(H,2) - 48*G2px*G4*G5p*G5ppxx*pow(H,2) + 48*G3pp*G4*G5p*G5ppxx*pow(H,2) - 72*G4p*G4pp*G5p*G5ppxx*pow(H,2) - 
   48*G2p*G4x*G5p*G5ppxx*pow(H,2) + 12*G2p*pow(G5p,2)*G5ppxx*pow(H,2) - 48*G2xx*G3px*G4*G5px*pow(H,2) - 96*pow(G3px,2)*G4*G5px*pow(H,2) + 120*G2pxx*G3x*G4*G5px*pow(H,2) + 
   24*G3ppx*G3x*G4*G5px*pow(H,2) + 72*G3pp*G3xx*G4*G5px*pow(H,2) + 72*G2xx*G3x*G4p*G5px*pow(H,2) - 432*G3px*G3x*G4p*G5px*pow(H,2) - 486*pow(G3x,2)*G4pp*G5px*pow(H,2) + 
   144*G3pxx*G4*G4pp*G5px*pow(H,2) + 72*G3xx*G4p*G4pp*G5px*pow(H,2) - 144*G3xx*G4*G4ppp*G5px*pow(H,2) - 288*G3x*G4*G4pppx*G5px*pow(H,2) + 192*G2xx*G4*G4ppx*G5px*pow(H,2) + 
   384*G3px*G4*G4ppx*G5px*pow(H,2) - 144*G3x*G4p*G4ppx*G5px*pow(H,2) - 576*G4*pow(G4ppx,2)*G5px*pow(H,2) - 288*G4*G4pp*G4ppxx*G5px*pow(H,2) - 336*G2pxx*G4*G4px*G5px*pow(H,2) + 
   48*G3ppx*G4*G4px*G5px*pow(H,2) - 48*G2xx*G4p*G4px*G5px*pow(H,2) + 1344*G3px*G4p*G4px*G5px*pow(H,2) + 1800*G3x*G4pp*G4px*G5px*pow(H,2) + 576*G4*G4pppx*G4px*G5px*pow(H,2) - 
   864*G4p*G4ppx*G4px*G5px*pow(H,2) - 1368*G4pp*pow(G4px,2)*G5px*pow(H,2) - 192*G2px*G4*G4pxx*G5px*pow(H,2) + 48*G3pp*G4*G4pxx*G5px*pow(H,2) - 720*G4p*G4pp*G4pxx*G5px*pow(H,2) + 
   288*G4*G4ppp*G4pxx*G5px*pow(H,2) - 288*G2px*G3x*G4x*G5px*pow(H,2) + 660*G3pp*G3x*G4x*G5px*pow(H,2) + 72*G2p*G3xx*G4x*G5px*pow(H,2) + 144*G2pxx*G4p*G4x*G5px*pow(H,2) + 
   144*G3ppx*G4p*G4x*G5px*pow(H,2) - 24*G2xx*G4pp*G4x*G5px*pow(H,2) + 96*G3px*G4pp*G4x*G5px*pow(H,2) - 216*G3x*G4ppp*G4x*G5px*pow(H,2) - 576*G4p*G4pppx*G4x*G5px*pow(H,2) + 
   144*G4pp*G4ppx*G4x*G5px*pow(H,2) + 672*G2px*G4px*G4x*G5px*pow(H,2) - 1368*G3pp*G4px*G4x*G5px*pow(H,2) + 144*G4ppp*G4px*G4x*G5px*pow(H,2) - 336*G2p*G4pxx*G4x*G5px*pow(H,2) + 
   48*G2ppx*pow(G4x,2)*G5px*pow(H,2) - 48*G3ppp*pow(G4x,2)*G5px*pow(H,2) - 480*G2p*G3x*G4xx*G5px*pow(H,2) + 192*G2ppx*G4*G4xx*G5px*pow(H,2) - 192*G3ppp*G4*G4xx*G5px*pow(H,2) + 
   96*G2px*G4p*G4xx*G5px*pow(H,2) - 960*G3pp*G4p*G4xx*G5px*pow(H,2) - 576*pow(G4pp,2)*G4xx*G5px*pow(H,2) + 1728*G4p*G4ppp*G4xx*G5px*pow(H,2) + 1152*G2p*G4px*G4xx*G5px*pow(H,2) + 
   192*G2pp*G4x*G4xx*G5px*pow(H,2) + 120*G2px*G3x*G5p*G5px*pow(H,2) - 522*G3pp*G3x*G5p*G5px*pow(H,2) - 36*G2p*G3xx*G5p*G5px*pow(H,2) - 72*G2pxx*G4p*G5p*G5px*pow(H,2) - 
   72*G3ppx*G4p*G5p*G5px*pow(H,2) - 12*G2xx*G4pp*G5p*G5px*pow(H,2) - 168*G3px*G4pp*G5p*G5px*pow(H,2) + 540*G3x*G4ppp*G5p*G5px*pow(H,2) + 288*G4p*G4pppx*G5p*G5px*pow(H,2) + 
   216*G4pp*G4ppx*G5p*G5px*pow(H,2) - 192*G2px*G4px*G5p*G5px*pow(H,2) + 972*G3pp*G4px*G5p*G5px*pow(H,2) - 936*G4ppp*G4px*G5p*G5px*pow(H,2) + 168*G2p*G4pxx*G5p*G5px*pow(H,2) + 
   48*G2ppx*G4x*G5p*G5px*pow(H,2) - 48*G3ppp*G4x*G5p*G5px*pow(H,2) - 96*G2pp*G4xx*G5p*G5px*pow(H,2) - 36*G2ppx*pow(G5p,2)*G5px*pow(H,2) + 36*G3ppp*pow(G5p,2)*G5px*pow(H,2) + 
   48*G2pxx*G4*G5pp*G5px*pow(H,2) - 48*G3ppx*G4*G5pp*G5px*pow(H,2) - 48*G2xx*G4p*G5pp*G5px*pow(H,2) - 240*G3px*G4p*G5pp*G5px*pow(H,2) - 108*G3x*G4pp*G5pp*G5px*pow(H,2) + 
   576*G4p*G4ppx*G5pp*G5px*pow(H,2) + 72*G4pp*G4px*G5pp*G5px*pow(H,2) - 144*G2px*G4x*G5pp*G5px*pow(H,2) + 192*G3pp*G4x*G5pp*G5px*pow(H,2) - 96*G2p*G4xx*G5pp*G5px*pow(H,2) + 
   24*G2px*G5p*G5pp*G5px*pow(H,2) - 48*G3pp*G5p*G5pp*G5px*pow(H,2) - 48*G2xx*G4*G5ppp*G5px*pow(H,2) + 48*G3px*G4*G5ppp*G5px*pow(H,2) + 432*G3x*G4p*G5ppp*G5px*pow(H,2) - 
   864*G4p*G4px*G5ppp*G5px*pow(H,2) + 96*G2px*G4*G5ppx*G5px*pow(H,2) - 96*G3pp*G4*G5ppx*G5px*pow(H,2) + 288*G4p*G4pp*G5ppx*G5px*pow(H,2) + 96*G2p*G4x*G5ppx*G5px*pow(H,2) - 
   48*G2p*G5p*G5ppx*G5px*pow(H,2) + 132*G2p*G3x*pow(G5px,2)*pow(H,2) - 48*G2ppx*G4*pow(G5px,2)*pow(H,2) + 48*G3ppp*G4*pow(G5px,2)*pow(H,2) + 216*G3pp*G4p*pow(G5px,2)*pow(H,2) + 
   144*pow(G4pp,2)*pow(G5px,2)*pow(H,2) - 432*G4p*G4ppp*pow(G5px,2)*pow(H,2) - 312*G2p*G4px*pow(G5px,2)*pow(H,2) - 48*G2pp*G4x*pow(G5px,2)*pow(H,2) + 
   24*G2pp*G5p*pow(G5px,2)*pow(H,2) + 24*G2p*G5pp*pow(G5px,2)*pow(H,2) + 24*G2px*G3x*G4*G5pxx*pow(H,2) + 12*G3pp*G3x*G4*G5pxx*pow(H,2) + 24*G2xx*G4*G4pp*G5pxx*pow(H,2) - 
   96*G3px*G4*G4pp*G5pxx*pow(H,2) + 144*G3x*G4p*G4pp*G5pxx*pow(H,2) - 72*G3x*G4*G4ppp*G5pxx*pow(H,2) + 144*G4*G4pp*G4ppx*G5pxx*pow(H,2) - 48*G2px*G4*G4px*G5pxx*pow(H,2) - 
   24*G3pp*G4*G4px*G5pxx*pow(H,2) - 144*G4p*G4pp*G4px*G5pxx*pow(H,2) + 144*G4*G4ppp*G4px*G5pxx*pow(H,2) + 60*G2p*G3x*G4x*G5pxx*pow(H,2) - 96*G2ppx*G4*G4x*G5pxx*pow(H,2) + 
   96*G3ppp*G4*G4x*G5pxx*pow(H,2) + 48*G2px*G4p*G4x*G5pxx*pow(H,2) + 24*G3pp*G4p*G4x*G5pxx*pow(H,2) - 144*pow(G4pp,2)*G4x*G5pxx*pow(H,2) - 144*G4p*G4ppp*G4x*G5pxx*pow(H,2) - 
   120*G2p*G4px*G4x*G5pxx*pow(H,2) - 48*G2pp*pow(G4x,2)*G5pxx*pow(H,2) - 30*G2p*G3x*G5p*G5pxx*pow(H,2) + 48*G2ppx*G4*G5p*G5pxx*pow(H,2) - 48*G3ppp*G4*G5p*G5pxx*pow(H,2) - 
   24*G2px*G4p*G5p*G5pxx*pow(H,2) - 12*G3pp*G4p*G5p*G5pxx*pow(H,2) + 72*pow(G4pp,2)*G5p*G5pxx*pow(H,2) + 72*G4p*G4ppp*G5p*G5pxx*pow(H,2) + 60*G2p*G4px*G5p*G5pxx*pow(H,2) + 
   48*G2pp*G4x*G5p*G5pxx*pow(H,2) - 12*G2pp*pow(G5p,2)*G5pxx*pow(H,2) - 72*G4p*G4pp*G5pp*G5pxx*pow(H,2) + 36*G2px*pow(G3x,2)*G5x*pow(H,2) - 117*G3pp*pow(G3x,2)*G5x*pow(H,2) - 
   18*G2p*G3x*G3xx*G5x*pow(H,2) + 48*G2xx*G3ppx*G4*G5x*pow(H,2) - 48*G2pxx*G3px*G4*G5x*pow(H,2) - 48*G2px*G3pxx*G4*G5x*pow(H,2) + 48*G3pp*G3pxx*G4*G5x*pow(H,2) + 
   48*G2ppx*G3xx*G4*G5x*pow(H,2) - 48*G3ppp*G3xx*G4*G5x*pow(H,2) + 12*G2xx*G3px*G4p*G5x*pow(H,2) + 60*pow(G3px,2)*G4p*G5x*pow(H,2) - 48*G2pxx*G3x*G4p*G5x*pow(H,2) - 
   24*G3ppx*G3x*G4p*G5x*pow(H,2) - 36*G3pp*G3xx*G4p*G5x*pow(H,2) + 12*G2xx*G3x*G4pp*G5x*pow(H,2) - 48*G3px*G3x*G4pp*G5x*pow(H,2) - 72*G3pxx*G4p*G4pp*G5x*pow(H,2) + 
   72*G3xx*pow(G4pp,2)*G5x*pow(H,2) + 90*pow(G3x,2)*G4ppp*G5x*pow(H,2) + 72*G3xx*G4p*G4ppp*G5x*pow(H,2) - 96*G2xx*G4*G4pppx*G5x*pow(H,2) + 96*G3px*G4*G4pppx*G5x*pow(H,2) + 
   144*G3x*G4p*G4pppx*G5x*pow(H,2) + 96*G2pxx*G4*G4ppx*G5x*pow(H,2) - 96*G3ppx*G4*G4ppx*G5x*pow(H,2) - 72*G2xx*G4p*G4ppx*G5x*pow(H,2) - 216*G3px*G4p*G4ppx*G5x*pow(H,2) + 
   288*G4p*pow(G4ppx,2)*G5x*pow(H,2) + 96*G2px*G4*G4ppxx*G5x*pow(H,2) - 96*G3pp*G4*G4ppxx*G5x*pow(H,2) + 144*G4p*G4pp*G4ppxx*G5x*pow(H,2) - 144*G2px*G3x*G4px*G5x*pow(H,2) + 
   456*G3pp*G3x*G4px*G5x*pow(H,2) + 36*G2p*G3xx*G4px*G5x*pow(H,2) + 144*G2pxx*G4p*G4px*G5x*pow(H,2) - 48*G2xx*G4pp*G4px*G5x*pow(H,2) + 192*G3px*G4pp*G4px*G5x*pow(H,2) - 
   288*G3x*G4ppp*G4px*G5x*pow(H,2) - 288*G4p*G4pppx*G4px*G5x*pow(H,2) - 144*G4pp*G4ppx*G4px*G5x*pow(H,2) + 144*G2px*pow(G4px,2)*G5x*pow(H,2) - 444*G3pp*pow(G4px,2)*G5x*pow(H,2) + 
   216*G4ppp*pow(G4px,2)*G5x*pow(H,2) + 84*G2p*G3x*G4pxx*G5x*pow(H,2) - 96*G2ppx*G4*G4pxx*G5x*pow(H,2) + 96*G3ppp*G4*G4pxx*G5x*pow(H,2) + 96*G2px*G4p*G4pxx*G5x*pow(H,2) - 
   24*G3pp*G4p*G4pxx*G5x*pow(H,2) - 144*pow(G4pp,2)*G4pxx*G5x*pow(H,2) - 144*G4p*G4ppp*G4pxx*G5x*pow(H,2) - 168*G2p*G4px*G4pxx*G5x*pow(H,2) - 48*G2px*G2xx*G4x*G5x*pow(H,2) + 
   156*G2xx*G3pp*G4x*G5x*pow(H,2) + 72*G2px*G3px*G4x*G5x*pow(H,2) - 228*G3pp*G3px*G4x*G5x*pow(H,2) - 48*G2p*G3pxx*G4x*G5x*pow(H,2) + 24*G2ppx*G3x*G4x*G5x*pow(H,2) - 
   24*G3ppp*G3x*G4x*G5x*pow(H,2) + 48*G2pp*G3xx*G4x*G5x*pow(H,2) + 120*G2pxx*G4pp*G4x*G5x*pow(H,2) - 120*G3ppx*G4pp*G4x*G5x*pow(H,2) - 120*G2xx*G4ppp*G4x*G5x*pow(H,2) + 
   120*G3px*G4ppp*G4x*G5x*pow(H,2) - 144*G2px*G4ppx*G4x*G5x*pow(H,2) + 240*G3pp*G4ppx*G4x*G5x*pow(H,2) + 96*G2p*G4ppxx*G4x*G5x*pow(H,2) + 48*G2ppx*G4px*G4x*G5x*pow(H,2) - 
   48*G3ppp*G4px*G4x*G5x*pow(H,2) - 96*G2pp*G4pxx*G4x*G5x*pow(H,2) - 24*G2p*G2xx*G4xx*G5x*pow(H,2) + 72*G2p*G3px*G4xx*G5x*pow(H,2) - 48*G2pp*G3x*G4xx*G5x*pow(H,2) - 
   96*G2ppx*G4p*G4xx*G5x*pow(H,2) + 96*G3ppp*G4p*G4xx*G5x*pow(H,2) + 144*G2px*G4pp*G4xx*G5x*pow(H,2) - 192*G3pp*G4pp*G4xx*G5x*pow(H,2) - 96*G2p*G4ppx*G4xx*G5x*pow(H,2) + 
   96*G2pp*G4px*G4xx*G5x*pow(H,2) + 24*G2px*G2xx*G5p*G5x*pow(H,2) - 90*G2xx*G3pp*G5p*G5x*pow(H,2) - 12*G2px*G3px*G5p*G5x*pow(H,2) + 102*G3pp*G3px*G5p*G5x*pow(H,2) + 
   24*G2p*G3pxx*G5p*G5x*pow(H,2) - 36*G2ppx*G3x*G5p*G5x*pow(H,2) + 36*G3ppp*G3x*G5p*G5x*pow(H,2) - 24*G2pp*G3xx*G5p*G5x*pow(H,2) - 84*G2pxx*G4pp*G5p*G5x*pow(H,2) + 
   84*G3ppx*G4pp*G5p*G5x*pow(H,2) + 84*G2xx*G4ppp*G5p*G5x*pow(H,2) - 84*G3px*G4ppp*G5p*G5x*pow(H,2) + 24*G2px*G4ppx*G5p*G5x*pow(H,2) - 72*G3pp*G4ppx*G5p*G5x*pow(H,2) - 
   48*G2p*G4ppxx*G5p*G5x*pow(H,2) + 24*G2ppx*G4px*G5p*G5x*pow(H,2) - 24*G3ppp*G4px*G5p*G5x*pow(H,2) + 48*G2pp*G4pxx*G5p*G5x*pow(H,2) + 24*G2px*G3x*G5pp*G5x*pow(H,2) - 
   36*G3pp*G3x*G5pp*G5x*pow(H,2) - 24*G2pxx*G4p*G5pp*G5x*pow(H,2) + 24*G3ppx*G4p*G5pp*G5x*pow(H,2) + 36*G2xx*G4pp*G5pp*G5x*pow(H,2) - 36*G3px*G4pp*G5pp*G5x*pow(H,2) - 
   48*G2px*G4px*G5pp*G5x*pow(H,2) + 72*G3pp*G4px*G5pp*G5x*pow(H,2) + 24*G2xx*G4p*G5ppp*G5x*pow(H,2) - 24*G3px*G4p*G5ppp*G5x*pow(H,2) - 24*G2p*G3x*G5ppx*G5x*pow(H,2) - 
   48*G2px*G4p*G5ppx*G5x*pow(H,2) + 48*G3pp*G4p*G5ppx*G5x*pow(H,2) + 48*G2p*G4px*G5ppx*G5x*pow(H,2) + 12*G2p*G2xx*G5px*G5x*pow(H,2) - 36*G2p*G3px*G5px*G5x*pow(H,2) + 
   24*G2pp*G3x*G5px*G5x*pow(H,2) + 48*G2ppx*G4p*G5px*G5x*pow(H,2) - 48*G3ppp*G4p*G5px*G5x*pow(H,2) - 72*G2px*G4pp*G5px*G5x*pow(H,2) + 96*G3pp*G4pp*G5px*G5x*pow(H,2) + 
   48*G2p*G4ppx*G5px*G5x*pow(H,2) - 48*G2pp*G4px*G5px*G5x*pow(H,2) + 4*pow(G2px,2)*pow(G5x,2)*pow(H,2) - 4*G2p*G2pxx*pow(G5x,2)*pow(H,2) + 4*G2pp*G2xx*pow(G5x,2)*pow(H,2) - 
   12*G2px*G3pp*pow(G5x,2)*pow(H,2) + 8*pow(G3pp,2)*pow(G5x,2)*pow(H,2) + 4*G2p*G3ppx*pow(G5x,2)*pow(H,2) - 4*G2pp*G3px*pow(G5x,2)*pow(H,2) + 27*G2p*pow(G3x,2)*G5xx*pow(H,2) - 
   12*G2xx*G3pp*G4*G5xx*pow(H,2) + 24*G2px*G3px*G4*G5xx*pow(H,2) - 12*G3pp*G3px*G4*G5xx*pow(H,2) - 24*G2ppx*G3x*G4*G5xx*pow(H,2) + 24*G3ppp*G3x*G4*G5xx*pow(H,2) - 
   24*G2px*G3x*G4p*G5xx*pow(H,2) + 132*G3pp*G3x*G4p*G5xx*pow(H,2) - 24*G2pxx*G4*G4pp*G5xx*pow(H,2) + 24*G3ppx*G4*G4pp*G5xx*pow(H,2) - 24*G2xx*G4p*G4pp*G5xx*pow(H,2) + 
   96*G3px*G4p*G4pp*G5xx*pow(H,2) + 72*G3x*pow(G4pp,2)*G5xx*pow(H,2) + 24*G2xx*G4*G4ppp*G5xx*pow(H,2) - 24*G3px*G4*G4ppp*G5xx*pow(H,2) - 216*G3x*G4p*G4ppp*G5xx*pow(H,2) - 
   48*G2px*G4*G4ppx*G5xx*pow(H,2) + 48*G3pp*G4*G4ppx*G5xx*pow(H,2) - 144*G4p*G4pp*G4ppx*G5xx*pow(H,2) - 132*G2p*G3x*G4px*G5xx*pow(H,2) + 48*G2ppx*G4*G4px*G5xx*pow(H,2) - 
   48*G3ppp*G4*G4px*G5xx*pow(H,2) - 216*G3pp*G4p*G4px*G5xx*pow(H,2) - 144*pow(G4pp,2)*G4px*G5xx*pow(H,2) + 432*G4p*G4ppp*G4px*G5xx*pow(H,2) + 156*G2p*pow(G4px,2)*G5xx*pow(H,2) - 
   12*G2p*G2xx*G4x*G5xx*pow(H,2) + 36*G2p*G3px*G4x*G5xx*pow(H,2) - 24*G2pp*G3x*G4x*G5xx*pow(H,2) - 48*G2ppx*G4p*G4x*G5xx*pow(H,2) + 48*G3ppp*G4p*G4x*G5xx*pow(H,2) + 
   72*G2px*G4pp*G4x*G5xx*pow(H,2) - 96*G3pp*G4pp*G4x*G5xx*pow(H,2) - 48*G2p*G4ppx*G4x*G5xx*pow(H,2) + 48*G2pp*G4px*G4x*G5xx*pow(H,2) + 6*G2p*G2xx*G5p*G5xx*pow(H,2) - 
   18*G2p*G3px*G5p*G5xx*pow(H,2) + 12*G2pp*G3x*G5p*G5xx*pow(H,2) + 24*G2ppx*G4p*G5p*G5xx*pow(H,2) - 24*G3ppp*G4p*G5p*G5xx*pow(H,2) - 36*G2px*G4pp*G5p*G5xx*pow(H,2) + 
   48*G3pp*G4pp*G5p*G5xx*pow(H,2) + 24*G2p*G4ppx*G5p*G5xx*pow(H,2) - 24*G2pp*G4px*G5p*G5xx*pow(H,2) + 12*G2p*G3x*G5pp*G5xx*pow(H,2) + 24*G2px*G4p*G5pp*G5xx*pow(H,2) - 
   24*G3pp*G4p*G5pp*G5xx*pow(H,2) - 24*G2p*G4px*G5pp*G5xx*pow(H,2) + 6*pow(G2x,2)*(4*G4xx*G5px - 2*pow(G5px,2) - 2*G4x*G5pxx + G5p*G5pxx - 4*G4pxx*G5x + 2*G5ppx*G5x + 2*G4px*G5xx - G5pp*G5xx)*
    pow(H,2) + 24*pow(G3p,2)*(4*G4xx*G5px - 2*pow(G5px,2) - 2*G4x*G5pxx + G5p*G5pxx - 4*G4pxx*G5x + 2*G5ppx*G5x + 2*G4px*G5xx - G5pp*G5xx)*pow(H,2) + 576*pow(G4,2)*G4pxx*G4pxxx*pow(H,4) + 
   2592*G3xx*G4*G4pxx*G4x*pow(H,4) - 7488*G4*pow(G4pxx,2)*G4x*pow(H,4) - 3024*G3x*G4*G4pxxx*G4x*pow(H,4) + 7200*G4*G4px*G4pxxx*G4x*pow(H,4) - 2592*G3x*G4pxx*pow(G4x,2)*pow(H,4) + 
   6912*G4px*G4pxx*pow(G4x,2)*pow(H,4) - 2592*G4p*G4pxxx*pow(G4x,2)*pow(H,4) + 1728*G3pxx*pow(G4x,3)*pow(H,4) - 2880*G4ppxx*pow(G4x,3)*pow(H,4) + 3456*G3x*G4*G4pxx*G4xx*pow(H,4) - 
   6336*G4*G4px*G4pxx*G4xx*pow(H,4) + 1440*G4*G4p*G4pxxx*G4xx*pow(H,4) - 3456*G3pxx*G4*G4x*G4xx*pow(H,4) - 2592*G3xx*G4p*G4x*G4xx*pow(H,4) + 9216*G4*G4ppxx*G4x*G4xx*pow(H,4) + 
   864*G3x*G4px*G4x*G4xx*pow(H,4) - 2880*pow(G4px,2)*G4x*G4xx*pow(H,4) + 13248*G4p*G4pxx*G4x*G4xx*pow(H,4) - 3024*G3px*pow(G4x,2)*G4xx*pow(H,4) - 288*G4ppx*pow(G4x,2)*G4xx*pow(H,4) - 
   7776*G3x*G4p*pow(G4xx,2)*pow(H,4) - 2880*G4*G4ppx*pow(G4xx,2)*pow(H,4) + 20160*G4p*G4px*pow(G4xx,2)*pow(H,4) + 12096*G4pp*G4x*pow(G4xx,2)*pow(H,4) - 
   576*pow(G4,2)*G4ppxx*G4xxx*pow(H,4) + 432*G3x*G4*G4px*G4xxx*pow(H,4) - 2016*G4*pow(G4px,2)*G4xxx*pow(H,4) - 576*G4*G4p*G4pxx*G4xxx*pow(H,4) + 2160*G3px*G4*G4x*G4xxx*pow(H,4) - 
   1728*G3x*G4p*G4x*G4xxx*pow(H,4) - 5472*G4*G4ppx*G4x*G4xxx*pow(H,4) + 2880*G4p*G4px*G4x*G4xxx*pow(H,4) + 864*G4pp*pow(G4x,2)*G4xxx*pow(H,4) + 1440*pow(G4p,2)*G4xx*G4xxx*pow(H,4) + 
   288*G4*G4pp*G4xx*G4xxx*pow(H,4) - 1728*G3xx*G4*G4pxx*G5p*pow(H,4) + 3456*G4*pow(G4pxx,2)*G5p*pow(H,4) + 2160*G3x*G4*G4pxxx*G5p*pow(H,4) - 5472*G4*G4px*G4pxxx*G5p*pow(H,4) - 
   864*G3xx*G4px*G4x*G5p*pow(H,4) + 6480*G3x*G4pxx*G4x*G5p*pow(H,4) - 14400*G4px*G4pxx*G4x*G5p*pow(H,4) + 3312*G4p*G4pxxx*G4x*G5p*pow(H,4) - 3888*G3pxx*pow(G4x,2)*G5p*pow(H,4) + 
   6048*G4ppxx*pow(G4x,2)*G5p*pow(H,4) + 2592*G3pxx*G4*G4xx*G5p*pow(H,4) + 2160*G3xx*G4p*G4xx*G5p*pow(H,4) - 5184*G4*G4ppxx*G4xx*G5p*pow(H,4) - 2808*G3x*G4px*G4xx*G5p*pow(H,4) + 
   8496*pow(G4px,2)*G4xx*G5p*pow(H,4) - 11808*G4p*G4pxx*G4xx*G5p*pow(H,4) + 4968*G3px*G4x*G4xx*G5p*pow(H,4) + 2160*G4ppx*G4x*G4xx*G5p*pow(H,4) - 12096*G4pp*pow(G4xx,2)*G5p*pow(H,4) - 
   1296*G3px*G4*G4xxx*G5p*pow(H,4) + 1512*G3x*G4p*G4xxx*G5p*pow(H,4) + 3744*G4*G4ppx*G4xxx*G5p*pow(H,4) - 2592*G4p*G4px*G4xxx*G5p*pow(H,4) - 720*G4pp*G4x*G4xxx*G5p*pow(H,4) + 
   648*G3xx*G4px*pow(G5p,2)*pow(H,4) - 3456*G3x*G4pxx*pow(G5p,2)*pow(H,4) + 7632*G4px*G4pxx*pow(G5p,2)*pow(H,4) - 1008*G4p*G4pxxx*pow(G5p,2)*pow(H,4) + 
   2808*G3pxx*G4x*pow(G5p,2)*pow(H,4) - 4032*G4ppxx*G4x*pow(G5p,2)*pow(H,4) - 2376*G3px*G4xx*pow(G5p,2)*pow(H,4) - 1440*G4ppx*G4xx*pow(G5p,2)*pow(H,4) + 
   144*G4pp*G4xxx*pow(G5p,2)*pow(H,4) - 648*G3pxx*pow(G5p,3)*pow(H,4) + 864*G4ppxx*pow(G5p,3)*pow(H,4) + 864*G3xx*pow(G4x,2)*G5pp*pow(H,4) - 432*G3xx*G4*G4xx*G5pp*pow(H,4) - 
   2592*G4*G4pxx*G4xx*G5pp*pow(H,4) + 6696*G3x*G4x*G4xx*G5pp*pow(H,4) - 18864*G4px*G4x*G4xx*G5pp*pow(H,4) - 3168*G4p*pow(G4xx,2)*G5pp*pow(H,4) - 648*G3x*G4*G4xxx*G5pp*pow(H,4) + 
   3024*G4*G4px*G4xxx*G5pp*pow(H,4) + 288*G4*G4pxxx*G5p*G5pp*pow(H,4) - 648*G3xx*G4x*G5p*G5pp*pow(H,4) - 1296*G4pxx*G4x*G5p*G5pp*pow(H,4) - 3888*G3x*G4xx*G5p*G5pp*pow(H,4) + 
   10800*G4px*G4xx*G5p*G5pp*pow(H,4) - 72*G4p*G4xxx*G5p*G5pp*pow(H,4) + 432*G4pxx*pow(G5p,2)*G5pp*pow(H,4) + 2016*G4x*G4xx*pow(G5pp,2)*pow(H,4) - 576*G4*G4xxx*pow(G5pp,2)*pow(H,4) - 
   1728*G4xx*G5p*pow(G5pp,2)*pow(H,4) + 576*pow(G4x,2)*G4xx*G5ppp*pow(H,4) + 2880*G4*pow(G4xx,2)*G5ppp*pow(H,4) - 2016*G4x*G4xx*G5p*G5ppp*pow(H,4) - 288*G4*G4xxx*G5p*G5ppp*pow(H,4) + 
   1728*G4xx*pow(G5p,2)*G5ppp*pow(H,4) - 288*pow(G4x,3)*G5pppx*pow(H,4) - 1152*G4*G4x*G4xx*G5pppx*pow(H,4) + 288*pow(G4,2)*G4xxx*G5pppx*pow(H,4) + 864*pow(G4x,2)*G5p*G5pppx*pow(H,4) - 
   792*G4x*pow(G5p,2)*G5pppx*pow(H,4) + 216*pow(G5p,3)*G5pppx*pow(H,4) - 288*pow(G4,2)*G4pxxx*G5ppx*pow(H,4) - 1248*G3xx*G4*G4x*G5ppx*pow(H,4) + 5376*G4*G4pxx*G4x*G5ppx*pow(H,4) + 
   144*G3x*pow(G4x,2)*G5ppx*pow(H,4) - 672*G4px*pow(G4x,2)*G5ppx*pow(H,4) - 1848*G3x*G4*G4xx*G5ppx*pow(H,4) + 3120*G4*G4px*G4xx*G5ppx*pow(H,4) - 2544*G4p*G4x*G4xx*G5ppx*pow(H,4) + 
   432*G4*G4p*G4xxx*G5ppx*pow(H,4) + 912*G3xx*G4*G5p*G5ppx*pow(H,4) - 2112*G4*G4pxx*G5p*G5ppx*pow(H,4) - 1608*G3x*G4x*G5p*G5ppx*pow(H,4) + 4176*G4px*G4x*G5p*G5ppx*pow(H,4) + 
   2136*G4p*G4xx*G5p*G5ppx*pow(H,4) + 1092*G3x*pow(G5p,2)*G5ppx*pow(H,4) - 3000*G4px*pow(G5p,2)*G5ppx*pow(H,4) - 1008*pow(G4x,2)*G5pp*G5ppx*pow(H,4) + 2016*G4*G4xx*G5pp*G5ppx*pow(H,4) + 
   1584*G4x*G5p*G5pp*G5ppx*pow(H,4) - 324*pow(G5p,2)*G5pp*G5ppx*pow(H,4) - 864*G4*G4x*pow(G5ppx,2)*pow(H,4) + 144*G4*G5p*pow(G5ppx,2)*pow(H,4) - 48*G3xx*pow(G4,2)*G5ppxx*pow(H,4) - 
   192*pow(G4,2)*G4pxx*G5ppxx*pow(H,4) + 1608*G3x*G4*G4x*G5ppxx*pow(H,4) - 3984*G4*G4px*G4x*G5ppxx*pow(H,4) + 1104*G4p*pow(G4x,2)*G5ppxx*pow(H,4) - 528*G4*G4p*G4xx*G5ppxx*pow(H,4) - 
   1080*G3x*G4*G5p*G5ppxx*pow(H,4) + 2832*G4*G4px*G5p*G5ppxx*pow(H,4) - 1368*G4p*G4x*G5p*G5ppxx*pow(H,4) + 408*G4p*pow(G5p,2)*G5ppxx*pow(H,4) - 144*G4*G5p*G5pp*G5ppxx*pow(H,4) + 
   144*pow(G4,2)*G5ppx*G5ppxx*pow(H,4) + 252*G3x*G3xx*G4*G5px*pow(H,4) - 1080*G3xx*G4*G4px*G5px*pow(H,4) - 2040*G3x*G4*G4pxx*G5px*pow(H,4) + 5808*G4*G4px*G4pxx*G5px*pow(H,4) - 
   864*G4*G4p*G4pxxx*G5px*pow(H,4) - 558*pow(G3x,2)*G4x*G5px*pow(H,4) + 1680*G3pxx*G4*G4x*G5px*pow(H,4) + 1152*G3xx*G4p*G4x*G5px*pow(H,4) - 5088*G4*G4ppxx*G4x*G5px*pow(H,4) + 
   2088*G3x*G4px*G4x*G5px*pow(H,4) - 1656*pow(G4px,2)*G4x*G5px*pow(H,4) - 5952*G4p*G4pxx*G4x*G5px*pow(H,4) - 120*G2xx*pow(G4x,2)*G5px*pow(H,4) + 1920*G3px*pow(G4x,2)*G5px*pow(H,4) - 
   912*G4ppx*pow(G4x,2)*G5px*pow(H,4) + 120*G2xx*G4*G4xx*G5px*pow(H,4) + 720*G3px*G4*G4xx*G5px*pow(H,4) + 8304*G3x*G4p*G4xx*G5px*pow(H,4) + 1488*G4*G4ppx*G4xx*G5px*pow(H,4) - 
   21600*G4p*G4px*G4xx*G5px*pow(H,4) - 14160*G4pp*G4x*G4xx*G5px*pow(H,4) - 720*pow(G4p,2)*G4xxx*G5px*pow(H,4) + 432*G4*G4pp*G4xxx*G5px*pow(H,4) + 144*pow(G3x,2)*G5p*G5px*pow(H,4) - 
   1344*G3pxx*G4*G5p*G5px*pow(H,4) - 900*G3xx*G4p*G5p*G5px*pow(H,4) + 2976*G4*G4ppxx*G5p*G5px*pow(H,4) + 1200*G3x*G4px*G5p*G5px*pow(H,4) - 4992*pow(G4px,2)*G5p*G5px*pow(H,4) + 
   5928*G4p*G4pxx*G5p*G5px*pow(H,4) + 36*G2xx*G4x*G5p*G5px*pow(H,4) - 2712*G3px*G4x*G5p*G5px*pow(H,4) - 72*G4ppx*G4x*G5p*G5px*pow(H,4) + 13560*G4pp*G4xx*G5p*G5px*pow(H,4) + 
   48*G2xx*pow(G5p,2)*G5px*pow(H,4) + 1272*G3px*pow(G5p,2)*G5px*pow(H,4) + 264*G4ppx*pow(G5p,2)*G5px*pow(H,4) + 576*G3xx*G4*G5pp*G5px*pow(H,4) + 288*G4*G4pxx*G5pp*G5px*pow(H,4) - 
   3864*G3x*G4x*G5pp*G5px*pow(H,4) + 11040*G4px*G4x*G5pp*G5px*pow(H,4) + 4152*G4p*G4xx*G5pp*G5px*pow(H,4) + 2250*G3x*G5p*G5pp*G5px*pow(H,4) - 5796*G4px*G5p*G5pp*G5px*pow(H,4) - 
   1152*G4x*pow(G5pp,2)*G5px*pow(H,4) + 864*G5p*pow(G5pp,2)*G5px*pow(H,4) - 144*pow(G4x,2)*G5ppp*G5px*pow(H,4) - 3168*G4*G4xx*G5ppp*G5px*pow(H,4) + 720*G4x*G5p*G5ppp*G5px*pow(H,4) - 
   756*pow(G5p,2)*G5ppp*G5px*pow(H,4) + 864*G4*G4x*G5pppx*G5px*pow(H,4) - 144*G4*G5p*G5pppx*G5px*pow(H,4) + 828*G3x*G4*G5ppx*G5px*pow(H,4) - 1800*G4*G4px*G5ppx*G5px*pow(H,4) + 
   1080*G4p*G4x*G5ppx*G5px*pow(H,4) - 1260*G4p*G5p*G5ppx*G5px*pow(H,4) - 864*G4*G5pp*G5ppx*G5px*pow(H,4) + 336*G4*G4p*G5ppxx*G5px*pow(H,4) - 96*G2xx*G4*pow(G5px,2)*pow(H,4) - 
   288*G3px*G4*pow(G5px,2)*pow(H,4) - 2184*G3x*G4p*pow(G5px,2)*pow(H,4) - 96*G4*G4ppx*pow(G5px,2)*pow(H,4) + 5424*G4p*G4px*pow(G5px,2)*pow(H,4) + 4200*G4pp*G4x*pow(G5px,2)*pow(H,4) - 
   3756*G4pp*G5p*pow(G5px,2)*pow(H,4) - 1176*G4p*G5pp*pow(G5px,2)*pow(H,4) + 864*G4*G5ppp*pow(G5px,2)*pow(H,4) + 144*pow(G3x,2)*G4*G5pxx*pow(H,4) + 48*G3pxx*pow(G4,2)*G5pxx*pow(H,4) - 
   72*G3xx*G4*G4p*G5pxx*pow(H,4) + 192*pow(G4,2)*G4ppxx*G5pxx*pow(H,4) - 864*G3x*G4*G4px*G5pxx*pow(H,4) + 2016*G4*pow(G4px,2)*G5pxx*pow(H,4) - 48*G4*G4p*G4pxx*G5pxx*pow(H,4) + 
   168*G2xx*G4*G4x*G5pxx*pow(H,4) - 1728*G3px*G4*G4x*G5pxx*pow(H,4) + 1608*G3x*G4p*G4x*G5pxx*pow(H,4) + 3888*G4*G4ppx*G4x*G5pxx*pow(H,4) - 2976*G4p*G4px*G4x*G5pxx*pow(H,4) - 
   1032*pow(G4p,2)*G4xx*G5pxx*pow(H,4) - 432*G4*G4pp*G4xx*G5pxx*pow(H,4) - 120*G2xx*G4*G5p*G5pxx*pow(H,4) + 816*G3px*G4*G5p*G5pxx*pow(H,4) - 1476*G3x*G4p*G5p*G5pxx*pow(H,4) - 
   2064*G4*G4ppx*G5p*G5pxx*pow(H,4) + 3144*G4p*G4px*G5p*G5pxx*pow(H,4) - 144*G4pp*G4x*G5p*G5pxx*pow(H,4) + 132*G3x*G4*G5pp*G5pxx*pow(H,4) - 1272*G4*G4px*G5pp*G5pxx*pow(H,4) - 
   288*G4p*G4x*G5pp*G5pxx*pow(H,4) - 12*G4p*G5p*G5pp*G5pxx*pow(H,4) + 288*G4*pow(G5pp,2)*G5pxx*pow(H,4) + 144*G4*G5p*G5ppp*G5pxx*pow(H,4) - 144*pow(G4,2)*G5pppx*G5pxx*pow(H,4) + 
   24*G4*G4p*G5ppx*G5pxx*pow(H,4) + 576*pow(G4p,2)*G5px*G5pxx*pow(H,4) - 144*G4*G4pp*G5px*G5pxx*pow(H,4) + 48*G3px*pow(G4,2)*G5pxxx*pow(H,4) + 96*G3x*G4*G4p*G5pxxx*pow(H,4) - 
   96*pow(G4,2)*G4ppx*G5pxxx*pow(H,4) - 288*G4*G4p*G4px*G5pxxx*pow(H,4) + 48*pow(G4p,2)*G4x*G5pxxx*pow(H,4) - 96*G4*G4pp*G4x*G5pxxx*pow(H,4) - 24*pow(G4p,2)*G5p*G5pxxx*pow(H,4) + 
   96*G4*G4pp*G5p*G5pxxx*pow(H,4) + 48*G4*G4p*G5pp*G5pxxx*pow(H,4) - 864*G3pxx*G3x*G4*G5x*pow(H,4) + 540*G3px*G3xx*G4*G5x*pow(H,4) - 540*G3x*G3xx*G4p*G5x*pow(H,4) - 
   1512*G3xx*G4*G4ppx*G5x*pow(H,4) + 1872*G3x*G4*G4ppxx*G5x*pow(H,4) + 162*pow(G3x,2)*G4px*G5x*pow(H,4) + 2160*G3pxx*G4*G4px*G5x*pow(H,4) + 1008*G3xx*G4p*G4px*G5x*pow(H,4) - 
   5184*G4*G4ppxx*G4px*G5x*pow(H,4) - 1152*G3x*pow(G4px,2)*G5x*pow(H,4) + 1944*pow(G4px,3)*G5x*pow(H,4) + 432*G2xx*G4*G4pxx*G5x*pow(H,4) - 1944*G3px*G4*G4pxx*G5x*pow(H,4) + 
   2952*G3x*G4p*G4pxx*G5x*pow(H,4) + 4464*G4*G4ppx*G4pxx*G5x*pow(H,4) - 6912*G4p*G4px*G4pxx*G5x*pow(H,4) + 144*pow(G4p,2)*G4pxxx*G5x*pow(H,4) - 432*G4*G4pp*G4pxxx*G5x*pow(H,4) - 
   864*G3px*G3x*G4x*G5x*pow(H,4) - 1368*G3pxx*G4p*G4x*G5x*pow(H,4) + 288*G3xx*G4pp*G4x*G5x*pow(H,4) - 1368*G3x*G4ppx*G4x*G5x*pow(H,4) + 1296*G4p*G4ppxx*G4x*G5x*pow(H,4) + 
   144*G2xx*G4px*G4x*G5x*pow(H,4) + 2448*G3px*G4px*G4x*G5x*pow(H,4) + 1008*G4ppx*G4px*G4x*G5x*pow(H,4) + 1152*G4pp*G4pxx*G4x*G5x*pow(H,4) + 936*G2pxx*pow(G4x,2)*G5x*pow(H,4) - 
   360*G3ppx*pow(G4x,2)*G5x*pow(H,4) - 1152*G4pppx*pow(G4x,2)*G5x*pow(H,4) - 648*G2pxx*G4*G4xx*G5x*pow(H,4) + 936*G3ppx*G4*G4xx*G5x*pow(H,4) - 504*G2xx*G4p*G4xx*G5x*pow(H,4) + 
   864*G3px*G4p*G4xx*G5x*pow(H,4) + 5832*G3x*G4pp*G4xx*G5x*pow(H,4) - 576*G4*G4pppx*G4xx*G5x*pow(H,4) + 4032*G4p*G4ppx*G4xx*G5x*pow(H,4) - 14688*G4pp*G4px*G4xx*G5x*pow(H,4) + 
   1224*G2px*G4x*G4xx*G5x*pow(H,4) - 3168*G3pp*G4x*G4xx*G5x*pow(H,4) + 1440*G4ppp*G4x*G4xx*G5x*pow(H,4) + 720*G2p*pow(G4xx,2)*G5x*pow(H,4) + 288*G2px*G4*G4xxx*G5x*pow(H,4) - 
   792*G3pp*G4*G4xxx*G5x*pow(H,4) - 288*G4p*G4pp*G4xxx*G5x*pow(H,4) + 432*G4*G4ppp*G4xxx*G5x*pow(H,4) - 360*G2p*G4x*G4xxx*G5x*pow(H,4) + 972*G3px*G3x*G5p*G5x*pow(H,4) + 
   864*G3pxx*G4p*G5p*G5x*pow(H,4) - 108*G3xx*G4pp*G5p*G5x*pow(H,4) + 864*G3x*G4ppx*G5p*G5x*pow(H,4) - 720*G4p*G4ppxx*G5p*G5x*pow(H,4) - 252*G2xx*G4px*G5p*G5x*pow(H,4) - 
   2592*G3px*G4px*G5p*G5x*pow(H,4) + 72*G4ppx*G4px*G5p*G5x*pow(H,4) - 504*G4pp*G4pxx*G5p*G5x*pow(H,4) - 1404*G2pxx*G4x*G5p*G5x*pow(H,4) + 252*G3ppx*G4x*G5p*G5x*pow(H,4) + 
   2304*G4pppx*G4x*G5p*G5x*pow(H,4) - 1044*G2px*G4xx*G5p*G5x*pow(H,4) + 3384*G3pp*G4xx*G5p*G5x*pow(H,4) - 2592*G4ppp*G4xx*G5p*G5x*pow(H,4) + 252*G2p*G4xxx*G5p*G5x*pow(H,4) + 
   504*G2pxx*pow(G5p,2)*G5x*pow(H,4) - 1008*G4pppx*pow(G5p,2)*G5x*pow(H,4) + 621*pow(G3x,2)*G5pp*G5x*pow(H,4) - 72*G3pxx*G4*G5pp*G5x*pow(H,4) - 36*G3xx*G4p*G5pp*G5x*pow(H,4) + 
   144*G4*G4ppxx*G5pp*G5x*pow(H,4) - 3240*G3x*G4px*G5pp*G5x*pow(H,4) + 3996*pow(G4px,2)*G5pp*G5x*pow(H,4) - 72*G4p*G4pxx*G5pp*G5x*pow(H,4) + 216*G2xx*G4x*G5pp*G5x*pow(H,4) + 
   648*G3px*G4x*G5pp*G5x*pow(H,4) - 2160*G4ppx*G4x*G5pp*G5x*pow(H,4) + 2592*G4pp*G4xx*G5pp*G5x*pow(H,4) - 54*G2xx*G5p*G5pp*G5x*pow(H,4) - 414*G3px*G5p*G5pp*G5x*pow(H,4) + 
   1296*G4ppx*G5p*G5pp*G5x*pow(H,4) + 324*G3x*pow(G5pp,2)*G5x*pow(H,4) - 720*G4px*pow(G5pp,2)*G5x*pow(H,4) + 72*G3xx*G4*G5ppp*G5x*pow(H,4) - 144*G4*G4pxx*G5ppp*G5x*pow(H,4) + 
   288*G3x*G4x*G5ppp*G5x*pow(H,4) - 144*G4px*G4x*G5ppp*G5x*pow(H,4) - 1728*G4p*G4xx*G5ppp*G5x*pow(H,4) - 540*G3x*G5p*G5ppp*G5x*pow(H,4) + 720*G4px*G5p*G5ppp*G5x*pow(H,4) - 
   72*G3x*G4*G5pppx*G5x*pow(H,4) + 432*G4*G4px*G5pppx*G5x*pow(H,4) + 720*G4p*G4x*G5pppx*G5x*pow(H,4) - 504*G4p*G5p*G5pppx*G5x*pow(H,4) - 220*G2xx*G4*G5ppx*G5x*pow(H,4) + 
   580*G3px*G4*G5ppx*G5x*pow(H,4) - 588*G3x*G4p*G5ppx*G5x*pow(H,4) - 1008*G4*G4ppx*G5ppx*G5x*pow(H,4) + 1680*G4p*G4px*G5ppx*G5x*pow(H,4) - 1224*G4pp*G4x*G5ppx*G5x*pow(H,4) + 
   612*G4pp*G5p*G5ppx*G5x*pow(H,4) + 144*G4p*G5pp*G5ppx*G5x*pow(H,4) - 48*pow(G4p,2)*G5ppxx*G5x*pow(H,4) + 216*G4*G4pp*G5ppxx*G5x*pow(H,4) + 328*G2pxx*G4*G5px*G5x*pow(H,4) - 
   616*G3ppx*G4*G5px*G5x*pow(H,4) + 228*G2xx*G4p*G5px*G5x*pow(H,4) - 432*G3px*G4p*G5px*G5x*pow(H,4) - 3264*G3x*G4pp*G5px*G5x*pow(H,4) + 576*G4*G4pppx*G5px*G5x*pow(H,4) - 
   1896*G4p*G4ppx*G5px*G5x*pow(H,4) + 8184*G4pp*G4px*G5px*G5x*pow(H,4) - 616*G2px*G4x*G5px*G5x*pow(H,4) + 1292*G3pp*G4x*G5px*G5x*pow(H,4) - 360*G4ppp*G4x*G5px*G5x*pow(H,4) - 
   648*G2p*G4xx*G5px*G5x*pow(H,4) + 536*G2px*G5p*G5px*G5x*pow(H,4) - 1534*G3pp*G5p*G5px*G5x*pow(H,4) + 1044*G4ppp*G5p*G5px*G5x*pow(H,4) - 1404*G4pp*G5pp*G5px*G5x*pow(H,4) + 
   792*G4p*G5ppp*G5px*G5x*pow(H,4) + 144*G2p*pow(G5px,2)*G5x*pow(H,4) - 200*G2px*G4*G5pxx*G5x*pow(H,4) + 468*G3pp*G4*G5pxx*G5x*pow(H,4) + 144*G4p*G4pp*G5pxx*G5x*pow(H,4) - 
   216*G4*G4ppp*G5pxx*G5x*pow(H,4) + 188*G2p*G4x*G5pxx*G5x*pow(H,4) - 154*G2p*G5p*G5pxx*G5x*pow(H,4) + 16*G2p*G4*G5pxxx*G5x*pow(H,4) + 144*G2px*G3x*pow(G5x,2)*pow(H,4) - 
   387*G3pp*G3x*pow(G5x,2)*pow(H,4) - 54*G2p*G3xx*pow(G5x,2)*pow(H,4) - 120*G2pxx*G4p*pow(G5x,2)*pow(H,4) - 96*G3ppx*G4p*pow(G5x,2)*pow(H,4) + 12*G2xx*G4pp*pow(G5x,2)*pow(H,4) + 
   186*G3px*G4pp*pow(G5x,2)*pow(H,4) + 198*G3x*G4ppp*pow(G5x,2)*pow(H,4) + 432*G4p*G4pppx*pow(G5x,2)*pow(H,4) - 612*G4pp*G4ppx*pow(G5x,2)*pow(H,4) - 
   300*G2px*G4px*pow(G5x,2)*pow(H,4) + 678*G3pp*G4px*pow(G5x,2)*pow(H,4) - 180*G4ppp*G4px*pow(G5x,2)*pow(H,4) + 156*G2p*G4pxx*pow(G5x,2)*pow(H,4) + 132*G2ppx*G4x*pow(G5x,2)*pow(H,4) - 
   132*G3ppp*G4x*pow(G5x,2)*pow(H,4) - 48*G2pp*G4xx*pow(G5x,2)*pow(H,4) - 126*G2ppx*G5p*pow(G5x,2)*pow(H,4) + 126*G3ppp*G5p*pow(G5x,2)*pow(H,4) + 90*G2px*G5pp*pow(G5x,2)*pow(H,4) - 
   114*G3pp*G5pp*pow(G5x,2)*pow(H,4) - 12*G2p*G5ppx*pow(G5x,2)*pow(H,4) + 12*G2pp*G5px*pow(G5x,2)*pow(H,4) + 108*G3px*G3x*G4*G5xx*pow(H,4) - 630*pow(G3x,2)*G4p*G5xx*pow(H,4) + 
   144*G3pxx*G4*G4p*G5xx*pow(H,4) + 144*G3xx*pow(G4p,2)*G5xx*pow(H,4) + 72*G3xx*G4*G4pp*G5xx*pow(H,4) - 552*G3x*G4*G4ppx*G5xx*pow(H,4) - 24*G2xx*G4*G4px*G5xx*pow(H,4) - 
   24*G3px*G4*G4px*G5xx*pow(H,4) + 2976*G3x*G4p*G4px*G5xx*pow(H,4) + 768*G4*G4ppx*G4px*G5xx*pow(H,4) - 3336*G4p*pow(G4px,2)*G5xx*pow(H,4) - 432*pow(G4p,2)*G4pxx*G5xx*pow(H,4) - 
   336*G4*G4pp*G4pxx*G5xx*pow(H,4) - 216*G2pxx*G4*G4x*G5xx*pow(H,4) + 600*G3ppx*G4*G4x*G5xx*pow(H,4) - 192*G2xx*G4p*G4x*G5xx*pow(H,4) + 288*G3px*G4p*G4x*G5xx*pow(H,4) + 
   1968*G3x*G4pp*G4x*G5xx*pow(H,4) - 768*G4*G4pppx*G4x*G5xx*pow(H,4) + 1152*G4p*G4ppx*G4x*G5xx*pow(H,4) - 5328*G4pp*G4px*G4x*G5xx*pow(H,4) + 120*G2px*pow(G4x,2)*G5xx*pow(H,4) - 
   264*G3pp*pow(G4x,2)*G5xx*pow(H,4) + 48*G4ppp*pow(G4x,2)*G5xx*pow(H,4) - 408*G2px*G4*G4xx*G5xx*pow(H,4) + 192*G3pp*G4*G4xx*G5xx*pow(H,4) - 3984*G4p*G4pp*G4xx*G5xx*pow(H,4) + 
   1248*G4*G4ppp*G4xx*G5xx*pow(H,4) + 240*G2p*G4x*G4xx*G5xx*pow(H,4) + 72*G2p*G4*G4xxx*G5xx*pow(H,4) + 168*G2pxx*G4*G5p*G5xx*pow(H,4) - 216*G3ppx*G4*G5p*G5xx*pow(H,4) + 
   156*G2xx*G4p*G5p*G5xx*pow(H,4) - 312*G3px*G4p*G5p*G5xx*pow(H,4) - 1812*G3x*G4pp*G5p*G5xx*pow(H,4) + 96*G4*G4pppx*G5p*G5xx*pow(H,4) - 1224*G4p*G4ppx*G5p*G5xx*pow(H,4) + 
   4488*G4pp*G4px*G5p*G5xx*pow(H,4) - 324*G2px*G4x*G5p*G5xx*pow(H,4) + 900*G3pp*G4x*G5p*G5xx*pow(H,4) - 504*G4ppp*G4x*G5p*G5xx*pow(H,4) - 552*G2p*G4xx*G5p*G5xx*pow(H,4) + 
   132*G2px*pow(G5p,2)*G5xx*pow(H,4) - 492*G3pp*pow(G5p,2)*G5xx*pow(H,4) + 456*G4ppp*pow(G5p,2)*G5xx*pow(H,4) - 12*G2xx*G4*G5pp*G5xx*pow(H,4) - 324*G3px*G4*G5pp*G5xx*pow(H,4) - 
   372*G3x*G4p*G5pp*G5xx*pow(H,4) + 816*G4*G4ppx*G5pp*G5xx*pow(H,4) + 1080*G4p*G4px*G5pp*G5xx*pow(H,4) + 1176*G4pp*G4x*G5pp*G5xx*pow(H,4) - 816*G4pp*G5p*G5pp*G5xx*pow(H,4) - 
   192*G4p*pow(G5pp,2)*G5xx*pow(H,4) + 384*G3x*G4*G5ppp*G5xx*pow(H,4) - 912*G4*G4px*G5ppp*G5xx*pow(H,4) - 384*G4p*G4x*G5ppp*G5xx*pow(H,4) + 624*G4p*G5p*G5ppp*G5xx*pow(H,4) - 
   144*G4*G4p*G5pppx*G5xx*pow(H,4) - 36*pow(G4p,2)*G5ppx*G5xx*pow(H,4) + 168*G4*G4pp*G5ppx*G5xx*pow(H,4) + 288*G2px*G4*G5px*G5xx*pow(H,4) - 156*G3pp*G4*G5px*G5xx*pow(H,4) + 
   2184*G4p*G4pp*G5px*G5xx*pow(H,4) - 696*G4*G4ppp*G5px*G5xx*pow(H,4) - 60*G2p*G4x*G5px*G5xx*pow(H,4) + 270*G2p*G5p*G5px*G5xx*pow(H,4) - 60*G2p*G4*G5pxx*G5xx*pow(H,4) + 
   84*G2p*G3x*G5x*G5xx*pow(H,4) + 40*G2ppx*G4*G5x*G5xx*pow(H,4) - 40*G3ppp*G4*G5x*G5xx*pow(H,4) - 168*G2px*G4p*G5x*G5xx*pow(H,4) + 504*G3pp*G4p*G5x*G5xx*pow(H,4) + 
   384*pow(G4pp,2)*G5x*G5xx*pow(H,4) - 480*G4p*G4ppp*G5x*G5xx*pow(H,4) - 156*G2p*G4px*G5x*G5xx*pow(H,4) - 8*G2pp*G4x*G5x*G5xx*pow(H,4) + 28*G2pp*G5p*G5x*G5xx*pow(H,4) - 
   12*G2p*G5pp*G5x*G5xx*pow(H,4) + 12*G2pp*G4*pow(G5xx,2)*pow(H,4) - 42*G2p*G4p*pow(G5xx,2)*pow(H,4) - 48*G3ppx*pow(G4,2)*G5xxx*pow(H,4) - 24*G3px*G4*G4p*G5xxx*pow(H,4) + 
   96*G3x*pow(G4p,2)*G5xxx*pow(H,4) - 24*G3x*G4*G4pp*G5xxx*pow(H,4) + 96*pow(G4,2)*G4pppx*G5xxx*pow(H,4) + 144*G4*G4p*G4ppx*G5xxx*pow(H,4) - 240*pow(G4p,2)*G4px*G5xxx*pow(H,4) + 
   240*G4*G4pp*G4px*G5xxx*pow(H,4) + 96*G2px*G4*G4x*G5xxx*pow(H,4) - 240*G3pp*G4*G4x*G5xxx*pow(H,4) - 48*G4p*G4pp*G4x*G5xxx*pow(H,4) + 96*G4*G4ppp*G4x*G5xxx*pow(H,4) - 
   48*G2p*pow(G4x,2)*G5xxx*pow(H,4) + 48*G2p*G4*G4xx*G5xxx*pow(H,4) - 48*G2px*G4*G5p*G5xxx*pow(H,4) + 144*G3pp*G4*G5p*G5xxx*pow(H,4) + 72*G4p*G4pp*G5p*G5xxx*pow(H,4) - 
   96*G4*G4ppp*G5p*G5xxx*pow(H,4) + 72*G2p*G4x*G5p*G5xxx*pow(H,4) - 24*G2p*pow(G5p,2)*G5xxx*pow(H,4) + 24*pow(G4p,2)*G5pp*G5xxx*pow(H,4) - 144*G4*G4pp*G5pp*G5xxx*pow(H,4) - 
   48*G4*G4p*G5ppp*G5xxx*pow(H,4) - 24*G2p*G4*G5px*G5xxx*pow(H,4) - 16*G2pp*G4*G5x*G5xxx*pow(H,4) + 12*G2p*G4p*G5x*G5xxx*pow(H,4) - 9648*pow(G4x,2)*G4xx*G5px*pow(H,6) + 
   4608*G4*pow(G4xx,2)*G5px*pow(H,6) + 2448*G4*G4x*G4xxx*G5px*pow(H,6) + 14040*G4x*G4xx*G5p*G5px*pow(H,6) - 1584*G4*G4xxx*G5p*G5px*pow(H,6) - 4824*G4xx*pow(G5p,2)*G5px*pow(H,6) + 
   5208*pow(G4x,2)*pow(G5px,2)*pow(H,6) - 4968*G4*G4xx*pow(G5px,2)*pow(H,6) - 7308*G4x*G5p*pow(G5px,2)*pow(H,6) + 2472*pow(G5p,2)*pow(G5px,2)*pow(H,6) + 
   1440*G4*pow(G5px,3)*pow(H,6) + 1296*pow(G4x,3)*G5pxx*pow(H,6) - 1296*G4*G4x*G4xx*G5pxx*pow(H,6) - 288*pow(G4,2)*G4xxx*G5pxx*pow(H,6) - 1944*pow(G4x,2)*G5p*G5pxx*pow(H,6) + 
   432*G4*G4xx*G5p*G5pxx*pow(H,6) + 432*G4x*pow(G5p,2)*G5pxx*pow(H,6) + 216*pow(G5p,3)*G5pxx*pow(H,6) - 976*G4*G4x*G5px*G5pxx*pow(H,6) + 512*G4*G5p*G5px*G5pxx*pow(H,6) + 
   256*pow(G4,2)*pow(G5pxx,2)*pow(H,6) - 864*G4*pow(G4x,2)*G5pxxx*pow(H,6) + 288*pow(G4,2)*G4xx*G5pxxx*pow(H,6) + 1440*G4*G4x*G5p*G5pxxx*pow(H,6) - 576*G4*pow(G5p,2)*G5pxxx*pow(H,6) - 
   112*pow(G4,2)*G5px*G5pxxx*pow(H,6) - 5616*G4*G4pxxx*G4x*G5x*pow(H,6) + 3456*G4pxx*pow(G4x,2)*G5x*pow(H,6) - 7776*G4px*G4x*G4xx*G5x*pow(H,6) - 12960*G4p*pow(G4xx,2)*G5x*pow(H,6) + 
   1296*G4*G4px*G4xxx*G5x*pow(H,6) - 4320*G4p*G4x*G4xxx*G5x*pow(H,6) + 4752*G4*G4pxxx*G5p*G5x*pow(H,6) - 1296*G4pxx*G4x*G5p*G5x*pow(H,6) + 5400*G4px*G4xx*G5p*G5x*pow(H,6) + 
   3672*G4p*G4xxx*G5p*G5x*pow(H,6) - 1728*G4pxx*pow(G5p,2)*G5x*pow(H,6) + 19656*G4x*G4xx*G5pp*G5x*pow(H,6) - 1944*G4*G4xxx*G5pp*G5x*pow(H,6) - 15552*G4xx*G5p*G5pp*G5x*pow(H,6) - 
   4800*pow(G4x,2)*G5ppx*G5x*pow(H,6) + 984*G4*G4xx*G5ppx*G5x*pow(H,6) + 5640*G4x*G5p*G5ppx*G5x*pow(H,6) - 1212*pow(G5p,2)*G5ppx*G5x*pow(H,6) + 2808*G4*G4x*G5ppxx*G5x*pow(H,6) - 
   2280*G4*G5p*G5ppxx*G5x*pow(H,6) + 540*G3xx*G4*G5px*G5x*pow(H,6) - 1272*G4*G4pxx*G5px*G5x*pow(H,6) - 3600*G3x*G4x*G5px*G5x*pow(H,6) + 12696*G4px*G4x*G5px*G5x*pow(H,6) + 
   17568*G4p*G4xx*G5px*G5x*pow(H,6) + 2700*G3x*G5p*G5px*G5x*pow(H,6) - 9000*G4px*G5p*G5px*G5x*pow(H,6) - 11040*G4x*G5pp*G5px*G5x*pow(H,6) + 8838*G5p*G5pp*G5px*G5x*pow(H,6) - 
   396*G4*G5ppx*G5px*G5x*pow(H,6) - 5652*G4p*pow(G5px,2)*G5x*pow(H,6) - 648*G3x*G4*G5pxx*G5x*pow(H,6) + 1056*G4*G4px*G5pxx*G5x*pow(H,6) + 1992*G4p*G4x*G5pxx*G5x*pow(H,6) - 
   2340*G4p*G5p*G5pxx*G5x*pow(H,6) + 540*G4*G5pp*G5pxx*G5x*pow(H,6) + 336*G4*G4p*G5pxxx*G5x*pow(H,6) - 756*G3pxx*G4*pow(G5x,2)*pow(H,6) - 648*G3xx*G4p*pow(G5x,2)*pow(H,6) + 
   1440*G4*G4ppxx*pow(G5x,2)*pow(H,6) - 810*G3x*G4px*pow(G5x,2)*pow(H,6) + 1764*pow(G4px,2)*pow(G5x,2)*pow(H,6) + 1296*G4p*G4pxx*pow(G5x,2)*pow(H,6) + 
   810*G3px*G4x*pow(G5x,2)*pow(H,6) - 6660*G4ppx*G4x*pow(G5x,2)*pow(H,6) + 6732*G4pp*G4xx*pow(G5x,2)*pow(H,6) - 432*G3px*G5p*pow(G5x,2)*pow(H,6) + 5580*G4ppx*G5p*pow(G5x,2)*pow(H,6) + 
   2079*G3x*G5pp*pow(G5x,2)*pow(H,6) - 5634*G4px*G5pp*pow(G5x,2)*pow(H,6) + 306*pow(G5pp,2)*pow(G5x,2)*pow(H,6) + 612*G4x*G5ppp*pow(G5x,2)*pow(H,6) - 
   846*G5p*G5ppp*pow(G5x,2)*pow(H,6) + 36*G4*G5pppx*pow(G5x,2)*pow(H,6) + 630*G4p*G5ppx*pow(G5x,2)*pow(H,6) - 3762*G4pp*G5px*pow(G5x,2)*pow(H,6) + 180*G2px*pow(G5x,3)*pow(H,6) - 
   495*G3pp*pow(G5x,3)*pow(H,6) + 270*G4ppp*pow(G5x,3)*pow(H,6) + 288*pow(G4,2)*G4pxxx*G5xx*pow(H,6) + 2304*G4*G4pxx*G4x*G5xx*pow(H,6) - 1872*G4px*pow(G4x,2)*G5xx*pow(H,6) - 
   144*G4*G4px*G4xx*G5xx*pow(H,6) - 10800*G4p*G4x*G4xx*G5xx*pow(H,6) + 432*G4*G4p*G4xxx*G5xx*pow(H,6) - 1440*G4*G4pxx*G5p*G5xx*pow(H,6) + 1872*G4px*G4x*G5p*G5xx*pow(H,6) + 
   8928*G4p*G4xx*G5p*G5xx*pow(H,6) - 432*G4px*pow(G5p,2)*G5xx*pow(H,6) + 4320*pow(G4x,2)*G5pp*G5xx*pow(H,6) - 1728*G4*G4xx*G5pp*G5xx*pow(H,6) - 6156*G4x*G5p*G5pp*G5xx*pow(H,6) + 
   2268*pow(G5p,2)*G5pp*G5xx*pow(H,6) - 656*G4*G4x*G5ppx*G5xx*pow(H,6) + 592*G4*G5p*G5ppx*G5xx*pow(H,6) - 256*pow(G4,2)*G5ppxx*G5xx*pow(H,6) + 1020*G3x*G4*G5px*G5xx*pow(H,6) - 
   2760*G4*G4px*G5px*G5xx*pow(H,6) + 7032*G4p*G4x*G5px*G5xx*pow(H,6) - 5616*G4p*G5p*G5px*G5xx*pow(H,6) + 1116*G4*G5pp*G5px*G5xx*pow(H,6) - 336*G4*G4p*G5pxx*G5xx*pow(H,6) + 
   396*G3px*G4*G5x*G5xx*pow(H,6) - 2304*G3x*G4p*G5x*G5xx*pow(H,6) - 1368*G4*G4ppx*G5x*G5xx*pow(H,6) + 6888*G4p*G4px*G5x*G5xx*pow(H,6) + 4704*G4pp*G4x*G5x*G5xx*pow(H,6) - 
   4356*G4pp*G5p*G5x*G5xx*pow(H,6) - 1560*G4p*G5pp*G5x*G5xx*pow(H,6) + 528*G4*G5ppp*G5x*G5xx*pow(H,6) + 81*G2p*pow(G5x,2)*G5xx*pow(H,6) + 300*pow(G4p,2)*pow(G5xx,2)*pow(H,6) + 
   168*G4*G4pp*pow(G5xx,2)*pow(H,6) - 288*pow(G4,2)*G4pxx*G5xxx*pow(H,6) + 576*G4*G4px*G4x*G5xxx*pow(H,6) - 576*G4p*pow(G4x,2)*G5xxx*pow(H,6) + 288*G4*G4p*G4xx*G5xxx*pow(H,6) - 
   288*G4*G4px*G5p*G5xxx*pow(H,6) + 1008*G4p*G4x*G5p*G5xxx*pow(H,6) - 432*G4p*pow(G5p,2)*G5xxx*pow(H,6) - 720*G4*G4x*G5pp*G5xxx*pow(H,6) + 432*G4*G5p*G5pp*G5xxx*pow(H,6) + 
   112*pow(G4,2)*G5ppx*G5xxx*pow(H,6) - 168*G4*G4p*G5px*G5xxx*pow(H,6) + 216*pow(G4p,2)*G5x*G5xxx*pow(H,6) - 120*G4*G4pp*G5x*G5xxx*pow(H,6) - 3330*G4x*G5px*pow(G5x,2)*pow(H,8) + 
   3060*G5p*G5px*pow(G5x,2)*pow(H,8) - 936*G4*G5pxx*pow(G5x,2)*pow(H,8) - 810*G4px*pow(G5x,3)*pow(H,8) + 1215*G5pp*pow(G5x,3)*pow(H,8) + 1260*G4*G5px*G5x*G5xx*pow(H,8) - 
   1674*G4p*pow(G5x,2)*G5xx*pow(H,8) - 2*G3p*(-24*G2xx*pow(G4px,2) + 72*G4ppx*pow(G4px,2) + 24*G2xx*G4ppx*G4x - 24*G2pxx*G4px*G4x + 24*G3ppx*G4px*G4x + 6*pow(G3px,2)*(2*G4x - G5p) - 
      12*G2xx*G4ppx*G5p + 12*G2pxx*G4px*G5p - 12*G3ppx*G4px*G5p + 12*G2xx*G4px*G5pp + 144*G3xx*G4*G4pxx*pow(H,2) - 288*G4*pow(G4pxx,2)*pow(H,2) + 144*G4*G4px*G4pxxx*pow(H,2) + 
      1920*G4px*G4pxx*G4x*pow(H,2) - 144*G4p*G4pxxx*G4x*pow(H,2) + 288*G3pxx*pow(G4x,2)*pow(H,2) - 672*G4ppxx*pow(G4x,2)*pow(H,2) - 144*G3pxx*G4*G4xx*pow(H,2) - 
      144*G3xx*G4p*G4xx*pow(H,2) + 288*G4*G4ppxx*G4xx*pow(H,2) - 2064*pow(G4px,2)*G4xx*pow(H,2) + 864*G4p*G4pxx*G4xx*pow(H,2) + 1200*G4ppx*G4x*G4xx*pow(H,2) + 
      288*G4pp*pow(G4xx,2)*pow(H,2) - 144*G4*G4ppx*G4xxx*pow(H,2) + 144*G4pp*G4x*G4xxx*pow(H,2) - 72*G3xx*G4px*G5p*pow(H,2) - 1104*G4px*G4pxx*G5p*pow(H,2) + 72*G4p*G4pxxx*G5p*pow(H,2) - 
      360*G3pxx*G4x*G5p*pow(H,2) + 816*G4ppxx*G4x*G5p*pow(H,2) - 888*G4ppx*G4xx*G5p*pow(H,2) - 72*G4pp*G4xxx*G5p*pow(H,2) + 108*G3pxx*pow(G5p,2)*pow(H,2) - 
      240*G4ppxx*pow(G5p,2)*pow(H,2) + 72*G3xx*G4x*G5pp*pow(H,2) - 240*G4pxx*G4x*G5pp*pow(H,2) + 840*G4px*G4xx*G5pp*pow(H,2) + 72*G4p*G4xxx*G5pp*pow(H,2) + 48*G4pxx*G5p*G5pp*pow(H,2) - 
      48*G4xx*pow(G5pp,2)*pow(H,2) - 96*G4x*G4xx*G5ppp*pow(H,2) + 48*G4xx*G5p*G5ppp*pow(H,2) + 48*pow(G4x,2)*G5pppx*pow(H,2) - 48*G4x*G5p*G5pppx*pow(H,2) + 
      12*pow(G5p,2)*G5pppx*pow(H,2) - 72*G3xx*G4*G5ppx*pow(H,2) + 144*G4*G4pxx*G5ppx*pow(H,2) - 912*G4px*G4x*G5ppx*pow(H,2) - 288*G4p*G4xx*G5ppx*pow(H,2) + 600*G4px*G5p*G5ppx*pow(H,2) + 
      48*G4x*G5pp*G5ppx*pow(H,2) - 24*G5p*G5pp*G5ppx*pow(H,2) - 72*G4*G4px*G5ppxx*pow(H,2) + 72*G4p*G4x*G5ppxx*pow(H,2) - 36*G4p*G5p*G5ppxx*pow(H,2) + 72*G3pxx*G4*G5px*pow(H,2) + 
      36*G3xx*G4p*G5px*pow(H,2) - 144*G4*G4ppxx*G5px*pow(H,2) + 1428*pow(G4px,2)*G5px*pow(H,2) - 360*G4p*G4pxx*G5px*pow(H,2) - 12*G2xx*G4x*G5px*pow(H,2) - 696*G4ppx*G4x*G5px*pow(H,2) - 
      240*G4pp*G4xx*G5px*pow(H,2) - 6*G2xx*G5p*G5px*pow(H,2) + 492*G4ppx*G5p*G5px*pow(H,2) - 492*G4px*G5pp*G5px*pow(H,2) + 24*pow(G5pp,2)*G5px*pow(H,2) + 48*G4x*G5ppp*G5px*pow(H,2) - 
      24*G5p*G5ppp*G5px*pow(H,2) + 144*G4p*G5ppx*G5px*pow(H,2) + 48*G4pp*pow(G5px,2)*pow(H,2) + 12*G2xx*G4*G5pxx*pow(H,2) + 72*G4*G4ppx*G5pxx*pow(H,2) - 72*G4p*G4px*G5pxx*pow(H,2) - 
      96*G4pp*G4x*G5pxx*pow(H,2) + 48*G4pp*G5p*G5pxx*pow(H,2) - 36*G4p*G5pp*G5pxx*pow(H,2) - 36*G3pxx*G4p*G5x*pow(H,2) + 36*G3xx*G4pp*G5x*pow(H,2) + 72*G4p*G4ppxx*G5x*pow(H,2) + 
      12*G2xx*G4px*G5x*pow(H,2) - 528*G4ppx*G4px*G5x*pow(H,2) - 120*G4pp*G4pxx*G5x*pow(H,2) + 108*G2pxx*G4x*G5x*pow(H,2) - 156*G3ppx*G4x*G5x*pow(H,2) + 96*G4pppx*G4x*G5x*pow(H,2) + 
      24*G2px*G4xx*G5x*pow(H,2) - 24*G3pp*G4xx*G5x*pow(H,2) - 48*G4ppp*G4xx*G5x*pow(H,2) - 66*G2pxx*G5p*G5x*pow(H,2) + 90*G3ppx*G5p*G5x*pow(H,2) - 48*G4pppx*G5p*G5x*pow(H,2) + 
      6*G2xx*G5pp*G5x*pow(H,2) + 24*G4ppx*G5pp*G5x*pow(H,2) + 24*G4px*G5ppp*G5x*pow(H,2) + 24*G4pp*G5ppx*G5x*pow(H,2) - 12*G2px*G5px*G5x*pow(H,2) + 12*G3pp*G5px*G5x*pow(H,2) + 
      24*G4ppp*G5px*G5x*pow(H,2) - 4*G2ppx*pow(G5x,2)*pow(H,2) + 4*G3ppp*pow(G5x,2)*pow(H,2) - 12*G2pxx*G4*G5xx*pow(H,2) + 12*G3ppx*G4*G5xx*pow(H,2) - 12*G2xx*G4p*G5xx*pow(H,2) - 
      72*G4p*G4ppx*G5xx*pow(H,2) - 48*G4pp*G4px*G5xx*pow(H,2) + 12*G2px*G4x*G5xx*pow(H,2) - 12*G3pp*G4x*G5xx*pow(H,2) - 24*G4ppp*G4x*G5xx*pow(H,2) - 6*G2px*G5p*G5xx*pow(H,2) + 
      6*G3pp*G5p*G5xx*pow(H,2) + 12*G4ppp*G5p*G5xx*pow(H,2) - 12*G4pp*G5pp*G5xx*pow(H,2) - 840*G4x*G4xx*G5px*pow(H,4) + 216*G4*G4xxx*G5px*pow(H,4) + 60*G4xx*G5p*G5px*pow(H,4) + 
      468*G4x*pow(G5px,2)*pow(H,4) - 6*G5p*pow(G5px,2)*pow(H,4) + 48*pow(G4x,2)*G5pxx*pow(H,4) + 312*G4*G4xx*G5pxx*pow(H,4) + 144*G4x*G5p*G5pxx*pow(H,4) - 
      120*pow(G5p,2)*G5pxx*pow(H,4) - 312*G4*G5px*G5pxx*pow(H,4) - 144*G4*G4x*G5pxxx*pow(H,4) + 96*G4*G5p*G5pxxx*pow(H,4) - 504*G4*G4pxxx*G5x*pow(H,4) - 720*G4pxx*G4x*G5x*pow(H,4) + 
      504*G4px*G4xx*G5x*pow(H,4) - 216*G4p*G4xxx*G5x*pow(H,4) + 936*G4pxx*G5p*G5x*pow(H,4) + 756*G4xx*G5pp*G5x*pow(H,4) + 64*G4x*G5ppx*G5x*pow(H,4) - 296*G5p*G5ppx*G5x*pow(H,4) + 
      268*G4*G5ppxx*G5x*pow(H,4) - 180*G4px*G5px*G5x*pow(H,4) - 402*G5pp*G5px*G5x*pow(H,4) + 216*G4p*G5pxx*G5x*pow(H,4) + 90*G4ppx*pow(G5x,2)*pow(H,4) - 24*G5ppp*pow(G5x,2)*pow(H,4) + 
      384*G4*G4pxx*G5xx*pow(H,4) - 1032*G4p*G4xx*G5xx*pow(H,4) - 264*G4px*G5p*G5xx*pow(H,4) + 444*G4x*G5pp*G5xx*pow(H,4) - 162*G5p*G5pp*G5xx*pow(H,4) - 168*G4*G5ppx*G5xx*pow(H,4) + 
      492*G4p*G5px*G5xx*pow(H,4) + 252*G4pp*G5x*G5xx*pow(H,4) + 48*G4*G4px*G5xxx*pow(H,4) - 48*G4p*G4x*G5xxx*pow(H,4) + 48*G4p*G5p*G5xxx*pow(H,4) - 48*G4*G5pp*G5xxx*pow(H,4) - 
      225*G5px*pow(G5x,2)*pow(H,6) + 9*pow(G3x,2)*(2*G4ppx + 5*G5px*pow(H,2)) + 
      6*G3x*(2*G2xx*G4px - 12*G4ppx*G4px + 2*G2pxx*G4x - 2*G3ppx*G4x - G2pxx*G5p + G3ppx*G5p - G2xx*G5pp - 12*G4*G4pxxx*pow(H,2) - 120*G4pxx*G4x*pow(H,2) + 132*G4px*G4xx*pow(H,2) - 
         12*G4p*G4xxx*pow(H,2) + 84*G4pxx*G5p*pow(H,2) - 42*G4xx*G5pp*pow(H,2) + 56*G4x*G5ppx*pow(H,2) - 40*G5p*G5ppx*pow(H,2) + 6*G4*G5ppxx*pow(H,2) - 114*G4px*G5px*pow(H,2) + 
         27*G5pp*G5px*pow(H,2) + 12*G4p*G5pxx*pow(H,2) + 34*G4ppx*G5x*pow(H,2) - 2*G5ppp*G5x*pow(H,2) + 6*G4pp*G5xx*pow(H,2) - 10*G5px*G5x*pow(H,4)) - 
      3*G3px*(3*pow(G3x,2) + 4*pow(G4px,2) + 4*G2xx*G4x + 8*G4ppx*G4x - 2*G2xx*G5p - 4*G4ppx*G5p + 4*G4px*G5pp + 216*G4x*G4xx*pow(H,2) - 24*G4*G4xxx*pow(H,2) - 156*G4xx*G5p*pow(H,2) - 
         128*G4x*G5px*pow(H,2) + 84*G5p*G5px*pow(H,2) + 16*G4*G5pxx*pow(H,2) - 88*G4px*G5x*pow(H,2) + 6*G5pp*G5x*pow(H,2) - 16*G4p*G5xx*pow(H,2) + 33*pow(G5x,2)*pow(H,4) - 
         2*G3x*(4*G4px + G5pp - 18*G5x*pow(H,2)))) + G2x*(-24*G2xx*pow(G4px,2) + 72*G4ppx*pow(G4px,2) + 24*G2xx*G4ppx*G4x - 24*G2pxx*G4px*G4x + 24*G3ppx*G4px*G4x + 6*pow(G3px,2)*(2*G4x - G5p) - 
      12*G2xx*G4ppx*G5p + 12*G2pxx*G4px*G5p - 12*G3ppx*G4px*G5p + 12*G2xx*G4px*G5pp + 144*G3xx*G4*G4pxx*pow(H,2) - 288*G4*pow(G4pxx,2)*pow(H,2) + 144*G4*G4px*G4pxxx*pow(H,2) + 
      1920*G4px*G4pxx*G4x*pow(H,2) - 144*G4p*G4pxxx*G4x*pow(H,2) + 288*G3pxx*pow(G4x,2)*pow(H,2) - 672*G4ppxx*pow(G4x,2)*pow(H,2) - 144*G3pxx*G4*G4xx*pow(H,2) - 
      144*G3xx*G4p*G4xx*pow(H,2) + 288*G4*G4ppxx*G4xx*pow(H,2) - 2064*pow(G4px,2)*G4xx*pow(H,2) + 864*G4p*G4pxx*G4xx*pow(H,2) + 1200*G4ppx*G4x*G4xx*pow(H,2) + 
      288*G4pp*pow(G4xx,2)*pow(H,2) - 144*G4*G4ppx*G4xxx*pow(H,2) + 144*G4pp*G4x*G4xxx*pow(H,2) - 72*G3xx*G4px*G5p*pow(H,2) - 1104*G4px*G4pxx*G5p*pow(H,2) + 72*G4p*G4pxxx*G5p*pow(H,2) - 
      360*G3pxx*G4x*G5p*pow(H,2) + 816*G4ppxx*G4x*G5p*pow(H,2) - 888*G4ppx*G4xx*G5p*pow(H,2) - 72*G4pp*G4xxx*G5p*pow(H,2) + 108*G3pxx*pow(G5p,2)*pow(H,2) - 
      240*G4ppxx*pow(G5p,2)*pow(H,2) + 72*G3xx*G4x*G5pp*pow(H,2) - 240*G4pxx*G4x*G5pp*pow(H,2) + 840*G4px*G4xx*G5pp*pow(H,2) + 72*G4p*G4xxx*G5pp*pow(H,2) + 48*G4pxx*G5p*G5pp*pow(H,2) - 
      48*G4xx*pow(G5pp,2)*pow(H,2) - 96*G4x*G4xx*G5ppp*pow(H,2) + 48*G4xx*G5p*G5ppp*pow(H,2) + 48*pow(G4x,2)*G5pppx*pow(H,2) - 48*G4x*G5p*G5pppx*pow(H,2) + 
      12*pow(G5p,2)*G5pppx*pow(H,2) - 72*G3xx*G4*G5ppx*pow(H,2) + 144*G4*G4pxx*G5ppx*pow(H,2) - 912*G4px*G4x*G5ppx*pow(H,2) - 288*G4p*G4xx*G5ppx*pow(H,2) + 600*G4px*G5p*G5ppx*pow(H,2) + 
      48*G4x*G5pp*G5ppx*pow(H,2) - 24*G5p*G5pp*G5ppx*pow(H,2) - 72*G4*G4px*G5ppxx*pow(H,2) + 72*G4p*G4x*G5ppxx*pow(H,2) - 36*G4p*G5p*G5ppxx*pow(H,2) + 72*G3pxx*G4*G5px*pow(H,2) + 
      36*G3xx*G4p*G5px*pow(H,2) - 144*G4*G4ppxx*G5px*pow(H,2) + 1428*pow(G4px,2)*G5px*pow(H,2) - 360*G4p*G4pxx*G5px*pow(H,2) - 12*G2xx*G4x*G5px*pow(H,2) - 696*G4ppx*G4x*G5px*pow(H,2) - 
      96*G3p*G4xx*G5px*pow(H,2) - 240*G4pp*G4xx*G5px*pow(H,2) - 6*G2xx*G5p*G5px*pow(H,2) + 492*G4ppx*G5p*G5px*pow(H,2) - 492*G4px*G5pp*G5px*pow(H,2) + 24*pow(G5pp,2)*G5px*pow(H,2) + 
      48*G4x*G5ppp*G5px*pow(H,2) - 24*G5p*G5ppp*G5px*pow(H,2) + 144*G4p*G5ppx*G5px*pow(H,2) + 48*G3p*pow(G5px,2)*pow(H,2) + 48*G4pp*pow(G5px,2)*pow(H,2) + 12*G2xx*G4*G5pxx*pow(H,2) + 
      72*G4*G4ppx*G5pxx*pow(H,2) - 72*G4p*G4px*G5pxx*pow(H,2) + 48*G3p*G4x*G5pxx*pow(H,2) - 96*G4pp*G4x*G5pxx*pow(H,2) - 24*G3p*G5p*G5pxx*pow(H,2) + 48*G4pp*G5p*G5pxx*pow(H,2) - 
      36*G4p*G5pp*G5pxx*pow(H,2) - 36*G3pxx*G4p*G5x*pow(H,2) + 36*G3xx*G4pp*G5x*pow(H,2) + 72*G4p*G4ppxx*G5x*pow(H,2) + 12*G2xx*G4px*G5x*pow(H,2) - 528*G4ppx*G4px*G5x*pow(H,2) + 
      96*G3p*G4pxx*G5x*pow(H,2) - 120*G4pp*G4pxx*G5x*pow(H,2) + 108*G2pxx*G4x*G5x*pow(H,2) - 156*G3ppx*G4x*G5x*pow(H,2) + 96*G4pppx*G4x*G5x*pow(H,2) + 24*G2px*G4xx*G5x*pow(H,2) - 
      24*G3pp*G4xx*G5x*pow(H,2) - 48*G4ppp*G4xx*G5x*pow(H,2) - 66*G2pxx*G5p*G5x*pow(H,2) + 90*G3ppx*G5p*G5x*pow(H,2) - 48*G4pppx*G5p*G5x*pow(H,2) + 6*G2xx*G5pp*G5x*pow(H,2) + 
      24*G4ppx*G5pp*G5x*pow(H,2) + 24*G4px*G5ppp*G5x*pow(H,2) - 48*G3p*G5ppx*G5x*pow(H,2) + 24*G4pp*G5ppx*G5x*pow(H,2) - 12*G2px*G5px*G5x*pow(H,2) + 12*G3pp*G5px*G5x*pow(H,2) + 
      24*G4ppp*G5px*G5x*pow(H,2) - 4*G2ppx*pow(G5x,2)*pow(H,2) + 4*G3ppp*pow(G5x,2)*pow(H,2) - 12*G2pxx*G4*G5xx*pow(H,2) + 12*G3ppx*G4*G5xx*pow(H,2) - 12*G2xx*G4p*G5xx*pow(H,2) - 
      72*G4p*G4ppx*G5xx*pow(H,2) - 48*G3p*G4px*G5xx*pow(H,2) - 48*G4pp*G4px*G5xx*pow(H,2) + 12*G2px*G4x*G5xx*pow(H,2) - 12*G3pp*G4x*G5xx*pow(H,2) - 24*G4ppp*G4x*G5xx*pow(H,2) - 
      6*G2px*G5p*G5xx*pow(H,2) + 6*G3pp*G5p*G5xx*pow(H,2) + 12*G4ppp*G5p*G5xx*pow(H,2) + 24*G3p*G5pp*G5xx*pow(H,2) - 12*G4pp*G5pp*G5xx*pow(H,2) - 840*G4x*G4xx*G5px*pow(H,4) + 
      216*G4*G4xxx*G5px*pow(H,4) + 60*G4xx*G5p*G5px*pow(H,4) + 468*G4x*pow(G5px,2)*pow(H,4) - 6*G5p*pow(G5px,2)*pow(H,4) + 48*pow(G4x,2)*G5pxx*pow(H,4) + 312*G4*G4xx*G5pxx*pow(H,4) + 
      144*G4x*G5p*G5pxx*pow(H,4) - 120*pow(G5p,2)*G5pxx*pow(H,4) - 312*G4*G5px*G5pxx*pow(H,4) - 144*G4*G4x*G5pxxx*pow(H,4) + 96*G4*G5p*G5pxxx*pow(H,4) - 504*G4*G4pxxx*G5x*pow(H,4) - 
      720*G4pxx*G4x*G5x*pow(H,4) + 504*G4px*G4xx*G5x*pow(H,4) - 216*G4p*G4xxx*G5x*pow(H,4) + 936*G4pxx*G5p*G5x*pow(H,4) + 756*G4xx*G5pp*G5x*pow(H,4) + 64*G4x*G5ppx*G5x*pow(H,4) - 
      296*G5p*G5ppx*G5x*pow(H,4) + 268*G4*G5ppxx*G5x*pow(H,4) - 180*G4px*G5px*G5x*pow(H,4) - 402*G5pp*G5px*G5x*pow(H,4) + 216*G4p*G5pxx*G5x*pow(H,4) + 90*G4ppx*pow(G5x,2)*pow(H,4) - 
      24*G5ppp*pow(G5x,2)*pow(H,4) + 384*G4*G4pxx*G5xx*pow(H,4) - 1032*G4p*G4xx*G5xx*pow(H,4) - 264*G4px*G5p*G5xx*pow(H,4) + 444*G4x*G5pp*G5xx*pow(H,4) - 162*G5p*G5pp*G5xx*pow(H,4) - 
      168*G4*G5ppx*G5xx*pow(H,4) + 492*G4p*G5px*G5xx*pow(H,4) + 252*G4pp*G5x*G5xx*pow(H,4) + 48*G4*G4px*G5xxx*pow(H,4) - 48*G4p*G4x*G5xxx*pow(H,4) + 48*G4p*G5p*G5xxx*pow(H,4) - 
      48*G4*G5pp*G5xxx*pow(H,4) - 225*G5px*pow(G5x,2)*pow(H,6) + 9*pow(G3x,2)*(2*G4ppx + 5*G5px*pow(H,2)) + 
      6*G3x*(2*G2xx*G4px - 12*G4ppx*G4px + 2*G2pxx*G4x - 2*G3ppx*G4x - G2pxx*G5p + G3ppx*G5p - G2xx*G5pp - 12*G4*G4pxxx*pow(H,2) - 120*G4pxx*G4x*pow(H,2) + 132*G4px*G4xx*pow(H,2) - 
         12*G4p*G4xxx*pow(H,2) + 84*G4pxx*G5p*pow(H,2) - 42*G4xx*G5pp*pow(H,2) + 56*G4x*G5ppx*pow(H,2) - 40*G5p*G5ppx*pow(H,2) + 6*G4*G5ppxx*pow(H,2) - 114*G4px*G5px*pow(H,2) + 
         27*G5pp*G5px*pow(H,2) + 12*G4p*G5pxx*pow(H,2) + 34*G4ppx*G5x*pow(H,2) - 2*G5ppp*G5x*pow(H,2) + 6*G4pp*G5xx*pow(H,2) - 10*G5px*G5x*pow(H,4)) - 
      3*G3px*(3*pow(G3x,2) + 4*pow(G4px,2) + 4*G2xx*G4x + 8*G4ppx*G4x - 2*G2xx*G5p - 4*G4ppx*G5p + 4*G4px*G5pp + 216*G4x*G4xx*pow(H,2) - 24*G4*G4xxx*pow(H,2) - 156*G4xx*G5p*pow(H,2) - 
         128*G4x*G5px*pow(H,2) + 84*G5p*G5px*pow(H,2) + 16*G4*G5pxx*pow(H,2) - 88*G4px*G5x*pow(H,2) + 6*G5pp*G5x*pow(H,2) - 16*G4p*G5xx*pow(H,2) + 33*pow(G5x,2)*pow(H,4) - 
         2*G3x*(4*G4px + G5pp - 18*G5x*pow(H,2)))) - 144*G3xx*G4pxx*G4x*pow(H,2)*rptot + 288*pow(G4pxx,2)*G4x*pow(H,2)*rptot + 72*G3x*G4pxxx*G4x*pow(H,2)*rptot - 
   144*G4px*G4pxxx*G4x*pow(H,2)*rptot + 144*G3xx*G4px*G4xx*pow(H,2)*rptot - 288*G3x*G4pxx*G4xx*pow(H,2)*rptot + 288*G4px*G4pxx*G4xx*pow(H,2)*rptot + 144*G3pxx*G4x*G4xx*pow(H,2)*rptot - 
   288*G4ppxx*G4x*G4xx*pow(H,2)*rptot - 144*G3px*pow(G4xx,2)*pow(H,2)*rptot + 288*G4ppx*pow(G4xx,2)*pow(H,2)*rptot + 72*G3x*G4px*G4xxx*pow(H,2)*rptot - 
   144*pow(G4px,2)*G4xxx*pow(H,2)*rptot - 72*G3px*G4x*G4xxx*pow(H,2)*rptot + 144*G4ppx*G4x*G4xxx*pow(H,2)*rptot + 72*G3xx*G4pxx*G5p*pow(H,2)*rptot - 144*pow(G4pxx,2)*G5p*pow(H,2)*rptot - 
   36*G3x*G4pxxx*G5p*pow(H,2)*rptot + 72*G4px*G4pxxx*G5p*pow(H,2)*rptot - 72*G3pxx*G4xx*G5p*pow(H,2)*rptot + 144*G4ppxx*G4xx*G5p*pow(H,2)*rptot + 36*G3px*G4xxx*G5p*pow(H,2)*rptot - 
   72*G4ppx*G4xxx*G5p*pow(H,2)*rptot - 72*G3xx*G4xx*G5pp*pow(H,2)*rptot + 144*G4pxx*G4xx*G5pp*pow(H,2)*rptot - 36*G3x*G4xxx*G5pp*pow(H,2)*rptot + 72*G4px*G4xxx*G5pp*pow(H,2)*rptot + 
   72*G3xx*G4x*G5ppx*pow(H,2)*rptot - 144*G4pxx*G4x*G5ppx*pow(H,2)*rptot + 144*G3x*G4xx*G5ppx*pow(H,2)*rptot - 288*G4px*G4xx*G5ppx*pow(H,2)*rptot - 36*G3xx*G5p*G5ppx*pow(H,2)*rptot + 
   72*G4pxx*G5p*G5ppx*pow(H,2)*rptot - 36*G3x*G4x*G5ppxx*pow(H,2)*rptot + 72*G4px*G4x*G5ppxx*pow(H,2)*rptot + 18*G3x*G5p*G5ppxx*pow(H,2)*rptot - 36*G4px*G5p*G5ppxx*pow(H,2)*rptot + 
   18*G3x*G3xx*G5px*pow(H,2)*rptot - 108*G3xx*G4px*G5px*pow(H,2)*rptot + 108*G3x*G4pxx*G5px*pow(H,2)*rptot - 72*G4px*G4pxx*G5px*pow(H,2)*rptot - 72*G3pxx*G4x*G5px*pow(H,2)*rptot + 
   144*G4ppxx*G4x*G5px*pow(H,2)*rptot + 24*G2xx*G4xx*G5px*pow(H,2)*rptot + 120*G3px*G4xx*G5px*pow(H,2)*rptot - 288*G4ppx*G4xx*G5px*pow(H,2)*rptot + 36*G3pxx*G5p*G5px*pow(H,2)*rptot - 
   72*G4ppxx*G5p*G5px*pow(H,2)*rptot + 36*G3xx*G5pp*G5px*pow(H,2)*rptot - 72*G4pxx*G5pp*G5px*pow(H,2)*rptot - 72*G3x*G5ppx*G5px*pow(H,2)*rptot + 144*G4px*G5ppx*G5px*pow(H,2)*rptot - 
   12*G2xx*pow(G5px,2)*pow(H,2)*rptot - 24*G3px*pow(G5px,2)*pow(H,2)*rptot + 72*G4ppx*pow(G5px,2)*pow(H,2)*rptot - 9*pow(G3x,2)*G5pxx*pow(H,2)*rptot + 
   36*pow(G4px,2)*G5pxx*pow(H,2)*rptot - 12*G2xx*G4x*G5pxx*pow(H,2)*rptot + 48*G3px*G4x*G5pxx*pow(H,2)*rptot - 72*G4ppx*G4x*G5pxx*pow(H,2)*rptot + 6*G2xx*G5p*G5pxx*pow(H,2)*rptot - 
   24*G3px*G5p*G5pxx*pow(H,2)*rptot + 36*G4ppx*G5p*G5pxx*pow(H,2)*rptot + 18*G3x*G5pp*G5pxx*pow(H,2)*rptot - 36*G4px*G5pp*G5pxx*pow(H,2)*rptot + 18*G3pxx*G3x*G5x*pow(H,2)*rptot - 
   18*G3px*G3xx*G5x*pow(H,2)*rptot + 36*G3xx*G4ppx*G5x*pow(H,2)*rptot - 36*G3x*G4ppxx*G5x*pow(H,2)*rptot - 36*G3pxx*G4px*G5x*pow(H,2)*rptot + 72*G4ppxx*G4px*G5x*pow(H,2)*rptot - 
   24*G2xx*G4pxx*G5x*pow(H,2)*rptot + 60*G3px*G4pxx*G5x*pow(H,2)*rptot - 72*G4ppx*G4pxx*G5x*pow(H,2)*rptot + 24*G2pxx*G4xx*G5x*pow(H,2)*rptot - 24*G3ppx*G4xx*G5x*pow(H,2)*rptot + 
   12*G2xx*G5ppx*G5x*pow(H,2)*rptot - 12*G3px*G5ppx*G5x*pow(H,2)*rptot - 12*G2pxx*G5px*G5x*pow(H,2)*rptot + 12*G3ppx*G5px*G5x*pow(H,2)*rptot - 18*G3px*G3x*G5xx*pow(H,2)*rptot + 
   36*G3x*G4ppx*G5xx*pow(H,2)*rptot + 12*G2xx*G4px*G5xx*pow(H,2)*rptot + 24*G3px*G4px*G5xx*pow(H,2)*rptot - 72*G4ppx*G4px*G5xx*pow(H,2)*rptot + 12*G2pxx*G4x*G5xx*pow(H,2)*rptot - 
   12*G3ppx*G4x*G5xx*pow(H,2)*rptot - 6*G2pxx*G5p*G5xx*pow(H,2)*rptot + 6*G3ppx*G5p*G5xx*pow(H,2)*rptot - 6*G2xx*G5pp*G5xx*pow(H,2)*rptot + 6*G3px*G5pp*G5xx*pow(H,2)*rptot + 
   144*pow(G4xx,2)*G5px*pow(H,4)*rptot - 72*G4x*G4xxx*G5px*pow(H,4)*rptot - 36*G4xxx*G5p*G5px*pow(H,4)*rptot - 216*G4xx*pow(G5px,2)*pow(H,4)*rptot + 72*pow(G5px,3)*pow(H,4)*rptot - 
   240*G4x*G4xx*G5pxx*pow(H,4)*rptot + 72*G4*G4xxx*G5pxx*pow(H,4)*rptot + 264*G4xx*G5p*G5pxx*pow(H,4)*rptot + 192*G4x*G5px*G5pxx*pow(H,4)*rptot - 132*G5p*G5px*G5pxx*pow(H,4)*rptot - 
   36*G4*pow(G5pxx,2)*pow(H,4)*rptot + 48*pow(G4x,2)*G5pxxx*pow(H,4)*rptot - 48*G4*G4xx*G5pxxx*pow(H,4)*rptot - 72*G4x*G5p*G5pxxx*pow(H,4)*rptot + 24*pow(G5p,2)*G5pxxx*pow(H,4)*rptot + 
   24*G4*G5px*G5pxxx*pow(H,4)*rptot + 360*G4pxxx*G4x*G5x*pow(H,4)*rptot - 864*G4pxx*G4xx*G5x*pow(H,4)*rptot + 72*G4px*G4xxx*G5x*pow(H,4)*rptot - 252*G4pxxx*G5p*G5x*pow(H,4)*rptot + 
   36*G4xxx*G5pp*G5x*pow(H,4)*rptot + 360*G4xx*G5ppx*G5x*pow(H,4)*rptot - 180*G4x*G5ppxx*G5x*pow(H,4)*rptot + 126*G5p*G5ppxx*G5x*pow(H,4)*rptot + 504*G4pxx*G5px*G5x*pow(H,4)*rptot - 
   216*G5ppx*G5px*G5x*pow(H,4)*rptot - 48*G3x*G5pxx*G5x*pow(H,4)*rptot + 96*G4px*G5pxx*G5x*pow(H,4)*rptot - 18*G5pp*G5pxx*G5x*pow(H,4)*rptot - 12*G4p*G5pxxx*G5x*pow(H,4)*rptot + 
   54*G3pxx*pow(G5x,2)*pow(H,4)*rptot - 108*G4ppxx*pow(G5x,2)*pow(H,4)*rptot - 72*G4*G4pxxx*G5xx*pow(H,4)*rptot - 336*G4pxx*G4x*G5xx*pow(H,4)*rptot + 336*G4px*G4xx*G5xx*pow(H,4)*rptot - 
   72*G4p*G4xxx*G5xx*pow(H,4)*rptot + 312*G4pxx*G5p*G5xx*pow(H,4)*rptot - 24*G4xx*G5pp*G5xx*pow(H,4)*rptot + 132*G4x*G5ppx*G5xx*pow(H,4)*rptot - 138*G5p*G5ppx*G5xx*pow(H,4)*rptot + 
   36*G4*G5ppxx*G5xx*pow(H,4)*rptot + 24*G3x*G5px*G5xx*pow(H,4)*rptot - 288*G4px*G5px*G5xx*pow(H,4)*rptot + 30*G5pp*G5px*G5xx*pow(H,4)*rptot + 72*G4p*G5pxx*G5xx*pow(H,4)*rptot - 
   60*G3px*G5x*G5xx*pow(H,4)*rptot + 84*G4ppx*G5x*G5xx*pow(H,4)*rptot + 18*G4pp*pow(G5xx,2)*pow(H,4)*rptot + 48*G4*G4pxx*G5xxx*pow(H,4)*rptot - 48*G4p*G4xx*G5xxx*pow(H,4)*rptot - 
   24*G4px*G5p*G5xxx*pow(H,4)*rptot + 24*G4x*G5pp*G5xxx*pow(H,4)*rptot - 24*G4*G5ppx*G5xxx*pow(H,4)*rptot + 12*G4p*G5px*G5xxx*pow(H,4)*rptot + 12*G4pp*G5x*G5xxx*pow(H,4)*rptot + 
   9*G5pxx*pow(G5x,2)*pow(H,6)*rptot - 90*G5px*G5x*G5xx*pow(H,6)*rptot) - 
4*pow(a,14)*pow(dphi,3)*H*(-12*G2px*G2xx*pow(G4,2) + 24*G2xx*G3pp*pow(G4,2) - 24*G2x*G3ppx*pow(G4,2) + 48*G3p*G3ppx*pow(G4,2) + 48*G2px*G3px*pow(G4,2) - 72*G3pp*G3px*pow(G4,2) - 
   24*G2ppx*G3x*pow(G4,2) + 24*G3ppp*G3x*pow(G4,2) + 12*G2pp*G3xx*pow(G4,2) - 6*G2x*G3px*G4*G4p + 12*G3p*G3px*G4*G4p + 60*G2px*G3x*G4*G4p - 84*G3pp*G3x*G4*G4p + 96*G2x*G3x*pow(G4p,2) - 
   192*G3p*G3x*pow(G4p,2) - 24*G3ppx*G4*pow(G4p,2) - 12*G2xx*pow(G4p,3) - 6*G3px*pow(G4p,3) - 78*G2x*G3x*G4*G4pp + 156*G3p*G3x*G4*G4pp - 12*G2xx*G4*G4p*G4pp + 156*G3px*G4*G4p*G4pp + 
   270*G3x*pow(G4p,2)*G4pp - 72*G3x*G4*pow(G4pp,2) - 72*G3x*G4*G4p*G4ppp + 24*G2x*pow(G4,2)*G4pppx - 48*G3p*pow(G4,2)*G4pppx + 72*G4*pow(G4p,2)*G4pppx - 96*G2px*pow(G4,2)*G4ppx + 
   120*G3pp*pow(G4,2)*G4ppx + 108*pow(G4p,3)*G4ppx - 360*G4*G4p*G4pp*G4ppx + 12*pow(G2x,2)*G4*G4px - 48*G2x*G3p*G4*G4px + 48*pow(G3p,2)*G4*G4px + 72*G2ppx*pow(G4,2)*G4px - 
   72*G3ppp*pow(G4,2)*G4px - 192*G2px*G4*G4p*G4px + 228*G3pp*G4*G4p*G4px - 240*G2x*pow(G4p,2)*G4px + 480*G3p*pow(G4p,2)*G4px + 240*G2x*G4*G4pp*G4px - 480*G3p*G4*G4pp*G4px - 
   756*pow(G4p,2)*G4pp*G4px + 216*G4*pow(G4pp,2)*G4px + 216*G4*G4p*G4ppp*G4px - 24*G2pp*pow(G4,2)*G4pxx - 12*G2px*G2x*G4*G4x + 24*G2px*G3p*G4*G4x + 12*G2x*G3pp*G4*G4x - 24*G3p*G3pp*G4*G4x - 
   24*G2pp*G3x*G4*G4x - 12*pow(G2x,2)*G4p*G4x + 48*G2x*G3p*G4p*G4x - 48*pow(G3p,2)*G4p*G4x + 48*G2ppx*G4*G4p*G4x - 48*G3ppp*G4*G4p*G4x + 96*G2px*pow(G4p,2)*G4x - 252*G3pp*pow(G4p,2)*G4x - 
   72*G2px*G4*G4pp*G4x + 96*G3pp*G4*G4pp*G4x - 156*G2x*G4p*G4pp*G4x + 312*G3p*G4p*G4pp*G4x - 144*G4p*pow(G4pp,2)*G4x + 24*G2x*G4*G4ppp*G4x - 48*G3p*G4*G4ppp*G4x + 216*pow(G4p,2)*G4ppp*G4x + 
   96*G2pp*G4*G4px*G4x + 48*G2pp*G4p*pow(G4x,2) - 48*G2pp*G4*G4p*G4xx + 12*G2px*G2x*G4*G5p - 24*G2px*G3p*G4*G5p - 12*G2x*G3pp*G4*G5p + 24*G3p*G3pp*G4*G5p + 12*pow(G2x,2)*G4p*G5p - 
   48*G2x*G3p*G4p*G5p + 48*pow(G3p,2)*G4p*G5p - 48*G2ppx*G4*G4p*G5p + 48*G3ppp*G4*G4p*G5p - 72*G2px*pow(G4p,2)*G5p + 204*G3pp*pow(G4p,2)*G5p + 72*G2px*G4*G4pp*G5p - 96*G3pp*G4*G4pp*G5p + 
   126*G2x*G4p*G4pp*G5p - 252*G3p*G4p*G4pp*G5p + 144*G4p*pow(G4pp,2)*G5p - 24*G2x*G4*G4ppp*G5p + 48*G3p*G4*G4ppp*G5p - 216*pow(G4p,2)*G4ppp*G5p - 24*G2pp*G4*G4px*G5p - 72*G2pp*G4p*G4x*G5p + 
   24*G2pp*G4p*pow(G5p,2) - 12*pow(G2x,2)*G4*G5pp + 48*G2x*G3p*G4*G5pp - 48*pow(G3p,2)*G4*G5pp + 60*G2px*G4*G4p*G5pp - 72*G3pp*G4*G4p*G5pp + 24*G2x*pow(G4p,2)*G5pp - 48*G3p*pow(G4p,2)*G5pp - 
   36*G2x*G4*G4pp*G5pp + 72*G3p*G4*G4pp*G5pp + 108*pow(G4p,2)*G4pp*G5pp - 12*G2x*G4*G4p*G5ppp + 24*G3p*G4*G4p*G5ppp - 36*pow(G4p,3)*G5ppp + 24*G2pp*G4*G4p*G5px - 4*G2pp*G2x*G4*G5x + 
   8*G2pp*G3p*G4*G5x - 6*G2pp*pow(G4p,2)*G5x + 108*pow(G3x,2)*G4*G4p*pow(H,2) - 144*G3pxx*pow(G4,2)*G4p*pow(H,2) - 108*G3xx*G4*pow(G4p,2)*pow(H,2) + 144*G3xx*pow(G4,2)*G4pp*pow(H,2) - 
   216*G3x*pow(G4,2)*G4ppx*pow(H,2) + 288*pow(G4,2)*G4p*G4ppxx*pow(H,2) - 72*G2xx*pow(G4,2)*G4px*pow(H,2) + 360*G3px*pow(G4,2)*G4px*pow(H,2) - 828*G3x*G4*G4p*G4px*pow(H,2) - 
   144*pow(G4,2)*G4ppx*G4px*pow(H,2) + 1008*G4*G4p*pow(G4px,2)*pow(H,2) + 216*G2x*pow(G4,2)*G4pxx*pow(H,2) - 432*G3p*pow(G4,2)*G4pxx*pow(H,2) + 288*G4*pow(G4p,2)*G4pxx*pow(H,2) - 
   288*pow(G4,2)*G4pp*G4pxx*pow(H,2) - 144*G3ppx*pow(G4,2)*G4x*pow(H,2) + 72*G2xx*G4*G4p*G4x*pow(H,2) + 756*G3px*G4*G4p*G4x*pow(H,2) + 1008*G3x*pow(G4p,2)*G4x*pow(H,2) - 
   900*G3x*G4*G4pp*G4x*pow(H,2) + 144*pow(G4,2)*G4pppx*G4x*pow(H,2) - 2448*G4*G4p*G4ppx*G4x*pow(H,2) + 288*G2x*G4*G4px*G4x*pow(H,2) - 576*G3p*G4*G4px*G4x*pow(H,2) - 
   3312*pow(G4p,2)*G4px*G4x*pow(H,2) + 2448*G4*G4pp*G4px*G4x*pow(H,2) - 216*G2px*G4*pow(G4x,2)*pow(H,2) + 360*G3pp*G4*pow(G4x,2)*pow(H,2) - 360*G2x*G4p*pow(G4x,2)*pow(H,2) + 
   720*G3p*G4p*pow(G4x,2)*pow(H,2) - 1944*G4p*G4pp*pow(G4x,2)*pow(H,2) + 144*G4*G4ppp*pow(G4x,2)*pow(H,2) - 216*G2px*pow(G4,2)*G4xx*pow(H,2) + 432*G3pp*pow(G4,2)*G4xx*pow(H,2) + 
   72*G2x*G4*G4p*G4xx*pow(H,2) - 144*G3p*G4*G4p*G4xx*pow(H,2) - 432*pow(G4p,3)*G4xx*pow(H,2) + 144*G4*G4p*G4pp*G4xx*pow(H,2) + 144*G3ppx*pow(G4,2)*G5p*pow(H,2) - 
   72*G2xx*G4*G4p*G5p*pow(H,2) - 468*G3px*G4*G4p*G5p*pow(H,2) - 900*G3x*pow(G4p,2)*G5p*pow(H,2) + 612*G3x*G4*G4pp*G5p*pow(H,2) - 144*pow(G4,2)*G4pppx*G5p*pow(H,2) + 
   1584*G4*G4p*G4ppx*G5p*pow(H,2) - 288*G2x*G4*G4px*G5p*pow(H,2) + 576*G3p*G4*G4px*G5p*pow(H,2) + 2844*pow(G4p,2)*G4px*G5p*pow(H,2) - 1584*G4*G4pp*G4px*G5p*pow(H,2) + 
   432*G2px*G4*G4x*G5p*pow(H,2) - 720*G3pp*G4*G4x*G5p*pow(H,2) + 684*G2x*G4p*G4x*G5p*pow(H,2) - 1368*G3p*G4p*G4x*G5p*pow(H,2) + 3708*G4p*G4pp*G4x*G5p*pow(H,2) - 
   288*G4*G4ppp*G4x*G5p*pow(H,2) - 216*G2px*G4*pow(G5p,2)*pow(H,2) + 360*G3pp*G4*pow(G5p,2)*pow(H,2) - 324*G2x*G4p*pow(G5p,2)*pow(H,2) + 648*G3p*G4p*pow(G5p,2)*pow(H,2) - 
   1764*G4p*G4pp*pow(G5p,2)*pow(H,2) + 144*G4*G4ppp*pow(G5p,2)*pow(H,2) + 72*G2xx*pow(G4,2)*G5pp*pow(H,2) - 216*G3px*pow(G4,2)*G5pp*pow(H,2) + 396*G3x*G4*G4p*G5pp*pow(H,2) + 
   360*pow(G4,2)*G4ppx*G5pp*pow(H,2) - 468*G4*G4p*G4px*G5pp*pow(H,2) - 324*G2x*G4*G4x*G5pp*pow(H,2) + 648*G3p*G4*G4x*G5pp*pow(H,2) + 612*pow(G4p,2)*G4x*G5pp*pow(H,2) + 
   72*G4*G4pp*G4x*G5pp*pow(H,2) + 324*G2x*G4*G5p*G5pp*pow(H,2) - 648*G3p*G4*G5p*G5pp*pow(H,2) - 468*pow(G4p,2)*G5p*G5pp*pow(H,2) - 72*G4*G4pp*G5p*G5pp*pow(H,2) - 
   216*G4*G4p*pow(G5pp,2)*pow(H,2) + 72*G3x*pow(G4,2)*G5ppp*pow(H,2) - 216*pow(G4,2)*G4px*G5ppp*pow(H,2) - 216*G4*G4p*G4x*G5ppp*pow(H,2) + 216*G4*G4p*G5p*G5ppp*pow(H,2) - 
   152*G2x*pow(G4,2)*G5ppx*pow(H,2) + 304*G3p*pow(G4,2)*G5ppx*pow(H,2) - 24*G4*pow(G4p,2)*G5ppx*pow(H,2) + 176*G2px*pow(G4,2)*G5px*pow(H,2) - 328*G3pp*pow(G4,2)*G5px*pow(H,2) - 
   114*G2x*G4*G4p*G5px*pow(H,2) + 228*G3p*G4*G4p*G5px*pow(H,2) + 246*pow(G4p,3)*G5px*pow(H,2) + 84*G4*G4p*G4pp*G5px*pow(H,2) - 24*G2ppx*pow(G4,2)*G5x*pow(H,2) + 
   24*G3ppp*pow(G4,2)*G5x*pow(H,2) + 192*G2px*G4*G4p*G5x*pow(H,2) - 348*G3pp*G4*G4p*G5x*pow(H,2) + 198*G2x*pow(G4p,2)*G5x*pow(H,2) - 396*G3p*pow(G4p,2)*G5x*pow(H,2) - 
   174*G2x*G4*G4pp*G5x*pow(H,2) + 348*G3p*G4*G4pp*G5x*pow(H,2) + 630*pow(G4p,2)*G4pp*G5x*pow(H,2) - 72*G4*pow(G4pp,2)*G5x*pow(H,2) - 72*G4*G4p*G4ppp*G5x*pow(H,2) - 
   24*G2pp*G4*G5p*G5x*pow(H,2) + 28*G2pp*pow(G4,2)*G5xx*pow(H,2) + 1296*pow(G4,2)*G4pxx*G4x*pow(H,4) + 432*G4*G4px*pow(G4x,2)*pow(H,4) - 864*G4p*pow(G4x,3)*pow(H,4) - 
   1296*pow(G4,2)*G4px*G4xx*pow(H,4) + 2160*G4*G4p*G4x*G4xx*pow(H,4) - 1296*pow(G4,2)*G4pxx*G5p*pow(H,4) - 864*G4*G4px*G4x*G5p*pow(H,4) + 2376*G4p*pow(G4x,2)*G5p*pow(H,4) - 
   2160*G4*G4p*G4xx*G5p*pow(H,4) + 432*G4*G4px*pow(G5p,2)*pow(H,4) - 2160*G4p*G4x*pow(G5p,2)*pow(H,4) + 648*G4p*pow(G5p,3)*pow(H,4) - 648*G4*pow(G4x,2)*G5pp*pow(H,4) + 
   1296*pow(G4,2)*G4xx*G5pp*pow(H,4) + 1296*G4*G4x*G5p*G5pp*pow(H,4) - 648*G4*pow(G5p,2)*G5pp*pow(H,4) - 912*pow(G4,2)*G4x*G5ppx*pow(H,4) + 912*pow(G4,2)*G5p*G5ppx*pow(H,4) + 
   1128*pow(G4,2)*G4px*G5px*pow(H,4) - 1476*G4*G4p*G4x*G5px*pow(H,4) + 1764*G4*G4p*G5p*G5px*pow(H,4) - 984*pow(G4,2)*G5pp*G5px*pow(H,4) - 336*pow(G4,2)*G4p*G5pxx*pow(H,4) + 
   432*G3x*G4*G4p*G5x*pow(H,4) - 216*pow(G4,2)*G4ppx*G5x*pow(H,4) - 1260*G4*G4p*G4px*G5x*pow(H,4) + 900*pow(G4p,2)*G4x*G5x*pow(H,4) - 900*G4*G4pp*G4x*G5x*pow(H,4) - 
   792*pow(G4p,2)*G5p*G5x*pow(H,4) + 612*G4*G4pp*G5p*G5x*pow(H,4) + 180*G4*G4p*G5pp*G5x*pow(H,4) + 72*pow(G4,2)*G5ppp*G5x*pow(H,4) - 252*G4*pow(G4p,2)*G5xx*pow(H,4) + 
   336*pow(G4,2)*G4pp*G5xx*pow(H,4) + 324*G4*G4p*pow(G5x,2)*pow(H,6) + 
   G2p*(18*pow(G3x,2)*G4 - 12*G3pxx*pow(G4,2) - 18*G3xx*G4*G4p + 24*pow(G4,2)*G4ppxx + 156*G4*pow(G4px,2) + 84*G4*G4p*G4pxx + 12*G2xx*G4*G4x + 12*G3px*G4*G4x - 96*G4*G4ppx*G4x + 60*G4p*G4px*G4x - 
      12*G2x*pow(G4x,2) + 24*G3p*pow(G4x,2) - 48*G4pp*pow(G4x,2) + 12*G2x*G4*G4xx - 24*G3p*G4*G4xx + 60*pow(G4p,2)*G4xx + 48*G4*G4pp*G4xx - 12*G2xx*G4*G5p + 12*G3px*G4*G5p + 24*G4*G4ppx*G5p - 
      138*G4p*G4px*G5p + 18*G2x*G4x*G5p - 36*G3p*G4x*G5p + 72*G4pp*G4x*G5p - 6*G2x*pow(G5p,2) + 12*G3p*pow(G5p,2) - 24*G4pp*pow(G5p,2) - 48*G4*G4px*G5pp + 24*G4p*G4x*G5pp - 24*G4*G4p*G5ppx - 
      6*G2x*G4*G5px + 12*G3p*G4*G5px - 36*pow(G4p,2)*G5px - 24*G4*G4pp*G5px + 4*G2px*G4*G5x - 8*G3pp*G4*G5x + 3*G2x*G4p*G5x - 6*G3p*G4p*G5x + 12*G4p*G4pp*G5x + 72*pow(G4x,3)*pow(H,2) + 
      360*G4*G4x*G4xx*pow(H,2) - 252*pow(G4x,2)*G5p*pow(H,2) - 360*G4*G4xx*G5p*pow(H,2) + 288*G4x*pow(G5p,2)*pow(H,2) - 108*pow(G5p,3)*pow(H,2) - 240*G4*G4x*G5px*pow(H,2) + 
      264*G4*G5p*G5px*pow(H,2) - 28*pow(G4,2)*G5pxx*pow(H,2) - 246*G4*G4px*G5x*pow(H,2) - 198*G4p*G4x*G5x*pow(H,2) + 216*G4p*G5p*G5x*pow(H,2) + 48*G4*G5pp*G5x*pow(H,2) - 
      42*G4*G4p*G5xx*pow(H,2) + 54*G4*pow(G5x,2)*pow(H,4) - 6*G3x*(19*G4*G4px + G4p*(8*G4x - 11*G5p) - 4*G4*(G5pp + 3*G5x*pow(H,2)))) - 18*pow(G3x,2)*G4p*rptot + 18*G3pxx*G4*G4p*rptot + 
   18*G3xx*pow(G4p,2)*rptot - 18*G3xx*G4*G4pp*rptot + 18*G3x*G4*G4ppx*rptot - 36*G4*G4p*G4ppxx*rptot + 12*G2xx*G4*G4px*rptot - 30*G3px*G4*G4px*rptot + 90*G3x*G4p*G4px*rptot - 
   72*G4p*pow(G4px,2)*rptot + 12*G2x*G4*G4pxx*rptot - 24*G3p*G4*G4pxx*rptot - 72*pow(G4p,2)*G4pxx*rptot + 36*G4*G4pp*G4pxx*rptot + 12*G3ppx*G4*G4x*rptot - 12*G2xx*G4p*G4x*rptot + 
   12*G3px*G4p*G4x*rptot + 72*G3x*G4pp*G4x*rptot + 36*G4p*G4ppx*G4x*rptot - 180*G4pp*G4px*G4x*rptot + 12*G2px*pow(G4x,2)*rptot - 24*G3pp*pow(G4x,2)*rptot - 12*G2px*G4*G4xx*rptot + 
   24*G3pp*G4*G4xx*rptot - 12*G2x*G4p*G4xx*rptot + 24*G3p*G4p*G4xx*rptot - 72*G4p*G4pp*G4xx*rptot - 12*G3ppx*G4*G5p*rptot + 12*G2xx*G4p*G5p*rptot - 30*G3px*G4p*G5p*rptot - 54*G3x*G4pp*G5p*rptot + 
   18*G4p*G4ppx*G5p*rptot - 6*G2x*G4px*G5p*rptot + 12*G3p*G4px*G5p*rptot + 126*G4pp*G4px*G5p*rptot - 18*G2px*G4x*G5p*rptot + 36*G3pp*G4x*G5p*rptot + 6*G2px*pow(G5p,2)*rptot - 
   12*G3pp*pow(G5p,2)*rptot - 12*G2xx*G4*G5pp*rptot + 12*G3px*G4*G5pp*rptot - 18*G3x*G4p*G5pp*rptot + 18*G4p*G4px*G5pp*rptot + 6*G2x*G4x*G5pp*rptot - 12*G3p*G4x*G5pp*rptot - 6*G2x*G4*G5ppx*rptot + 
   12*G3p*G4*G5ppx*rptot + 18*pow(G4p,2)*G5ppx*rptot + 6*G2px*G4*G5px*rptot - 12*G3pp*G4*G5px*rptot + 3*G2x*G4p*G5px*rptot - 6*G3p*G4p*G5px*rptot + 36*G4p*G4pp*G5px*rptot - 3*G2px*G4p*G5x*rptot + 
   6*G3pp*G4p*G5x*rptot + 3*G2x*G4pp*G5x*rptot - 6*G3p*G4pp*G5x*rptot - 216*G4*G4pxx*G4x*pow(H,2)*rptot - 72*G4px*pow(G4x,2)*pow(H,2)*rptot + 216*G4*G4px*G4xx*pow(H,2)*rptot - 
   360*G4p*G4x*G4xx*pow(H,2)*rptot + 216*G4*G4pxx*G5p*pow(H,2)*rptot + 144*G4px*G4x*G5p*pow(H,2)*rptot + 360*G4p*G4xx*G5p*pow(H,2)*rptot - 72*G4px*pow(G5p,2)*pow(H,2)*rptot + 
   108*pow(G4x,2)*G5pp*pow(H,2)*rptot - 216*G4*G4xx*G5pp*pow(H,2)*rptot - 216*G4x*G5p*G5pp*pow(H,2)*rptot + 108*pow(G5p,2)*G5pp*pow(H,2)*rptot + 144*G4*G4x*G5ppx*pow(H,2)*rptot - 
   144*G4*G5p*G5ppx*pow(H,2)*rptot + 36*G3x*G4*G5px*pow(H,2)*rptot - 270*G4*G4px*G5px*pow(H,2)*rptot + 270*G4p*G4x*G5px*pow(H,2)*rptot - 288*G4p*G5p*G5px*pow(H,2)*rptot + 
   144*G4*G5pp*G5px*pow(H,2)*rptot + 42*G4*G4p*G5pxx*pow(H,2)*rptot - 36*G3px*G4*G5x*pow(H,2)*rptot - 72*G3x*G4p*G5x*pow(H,2)*rptot + 126*G4*G4ppx*G5x*pow(H,2)*rptot + 
   252*G4p*G4px*G5x*pow(H,2)*rptot + 162*G4pp*G4x*G5x*pow(H,2)*rptot - 144*G4pp*G5p*G5x*pow(H,2)*rptot - 72*G4p*G5pp*G5x*pow(H,2)*rptot + 42*pow(G4p,2)*G5xx*pow(H,2)*rptot - 
   42*G4*G4pp*G5xx*pow(H,2)*rptot - 54*G4p*pow(G5x,2)*pow(H,4)*rptot + 12*G2pxx*G4*(G2x*G4 - 2*G3p*G4 - pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2) - G4x*rptot + G5p*rptot)) - 
4*pow(a,15)*pow(dphi,2)*(4*pow(G2px,2)*pow(G4,2) - 4*G2ppx*G2x*pow(G4,2) + 4*G2pp*G2xx*pow(G4,2) + 8*G2ppx*G3p*pow(G4,2) + 8*pow(G3pp,2)*pow(G4,2) + 4*G2x*G3ppp*pow(G4,2) - 
   8*G3p*G3ppp*pow(G4,2) - 4*G2pp*G3px*pow(G4,2) - 6*G2x*G3pp*G4*G4p + 12*G3p*G3pp*G4*G4p - 12*G2pp*G3x*G4*G4p + 6*pow(G2x,2)*pow(G4p,2) - 24*G2x*G3p*pow(G4p,2) + 24*pow(G3p,2)*pow(G4p,2) - 
   12*G2ppx*G4*pow(G4p,2) + 12*G3ppp*G4*pow(G4p,2) + 30*G3pp*pow(G4p,3) - 6*pow(G2x,2)*G4*G4pp + 24*G2x*G3p*G4*G4pp - 24*pow(G3p,2)*G4*G4pp - 48*G3pp*G4*G4p*G4pp + 30*G2x*pow(G4p,2)*G4pp - 
   60*G3p*pow(G4p,2)*G4pp - 12*G2x*G4*pow(G4pp,2) + 24*G3p*G4*pow(G4pp,2) + 36*pow(G4p,2)*pow(G4pp,2) - 12*G2x*G4*G4p*G4ppp + 24*G3p*G4*G4p*G4ppp - 36*pow(G4p,3)*G4ppp + 24*G2pp*G4*G4p*G4px - 
   8*G2pp*G2x*G4*G4x + 16*G2pp*G3p*G4*G4x - 12*G2pp*pow(G4p,2)*G4x + 4*G2pp*G2x*G4*G5p - 8*G2pp*G3p*G4*G5p + 6*G2pp*pow(G4p,2)*G5p + 36*G2x*G3px*pow(G4,2)*pow(H,2) - 
   72*G3p*G3px*pow(G4,2)*pow(H,2) + 72*G3pp*G3x*pow(G4,2)*pow(H,2) + 18*G2x*G3x*G4*G4p*pow(H,2) - 36*G3p*G3x*G4*G4p*pow(H,2) - 48*G2pxx*pow(G4,2)*G4p*pow(H,2) + 
   48*G3ppx*pow(G4,2)*G4p*pow(H,2) - 36*G2xx*G4*pow(G4p,2)*pow(H,2) - 162*G3x*pow(G4p,3)*pow(H,2) + 48*G2xx*pow(G4,2)*G4pp*pow(H,2) - 48*G3px*pow(G4,2)*G4pp*pow(H,2) + 
   144*G3x*G4*G4p*G4pp*pow(H,2) - 144*G2x*pow(G4,2)*G4ppx*pow(H,2) + 288*G3p*pow(G4,2)*G4ppx*pow(H,2) + 144*G4*pow(G4p,2)*G4ppx*pow(H,2) - 312*G3pp*pow(G4,2)*G4px*pow(H,2) - 
   132*G2x*G4*G4p*G4px*pow(H,2) + 264*G3p*G4*G4p*G4px*pow(H,2) + 432*pow(G4p,3)*G4px*pow(H,2) - 144*G4*G4p*G4pp*G4px*pow(H,2) - 24*G2ppx*pow(G4,2)*G4x*pow(H,2) + 
   24*G3ppp*pow(G4,2)*G4x*pow(H,2) - 564*G3pp*G4*G4p*G4x*pow(H,2) + 276*G2x*pow(G4p,2)*G4x*pow(H,2) - 552*G3p*pow(G4p,2)*G4x*pow(H,2) - 264*G2x*G4*G4pp*G4x*pow(H,2) + 
   528*G3p*G4*G4pp*G4x*pow(H,2) + 900*pow(G4p,2)*G4pp*G4x*pow(H,2) - 72*G4*pow(G4pp,2)*G4x*pow(H,2) - 72*G4*G4p*G4ppp*G4x*pow(H,2) + 96*G2pp*pow(G4,2)*G4xx*pow(H,2) + 
   24*G2ppx*pow(G4,2)*G5p*pow(H,2) - 24*G3ppp*pow(G4,2)*G5p*pow(H,2) + 468*G3pp*G4*G4p*G5p*pow(H,2) - 258*G2x*pow(G4p,2)*G5p*pow(H,2) + 516*G3p*pow(G4p,2)*G5p*pow(H,2) + 
   216*G2x*G4*G4pp*G5p*pow(H,2) - 432*G3p*G4*G4pp*G5p*pow(H,2) - 828*pow(G4p,2)*G4pp*G5p*pow(H,2) + 72*G4*pow(G4pp,2)*G5p*pow(H,2) + 72*G4*G4p*G4ppp*G5p*pow(H,2) - 
   24*G2pp*G4*G4x*G5p*pow(H,2) + 24*G2pp*G4*pow(G5p,2)*pow(H,2) + 48*G3pp*pow(G4,2)*G5pp*pow(H,2) + 114*G2x*G4*G4p*G5pp*pow(H,2) - 228*G3p*G4*G4p*G5pp*pow(H,2) - 
   18*pow(G4p,3)*G5pp*pow(H,2) - 144*G4*G4p*G4pp*G5pp*pow(H,2) + 12*G2x*pow(G4,2)*G5ppp*pow(H,2) - 24*G3p*pow(G4,2)*G5ppp*pow(H,2) + 36*G4*pow(G4p,2)*G5ppp*pow(H,2) - 
   60*G2pp*pow(G4,2)*G5px*pow(H,2) - 36*G2pp*G4*G4p*G5x*pow(H,2) - 216*G3x*pow(G4,2)*G4px*pow(H,4) + 864*pow(G4,2)*pow(G4px,2)*pow(H,4) - 1152*pow(G4,2)*G4p*G4pxx*pow(H,4) + 
   216*G3px*pow(G4,2)*G4x*pow(H,4) + 540*G3x*G4*G4p*G4x*pow(H,4) - 864*pow(G4,2)*G4ppx*G4x*pow(H,4) - 1656*G4*G4p*G4px*G4x*pow(H,4) + 720*pow(G4p,2)*pow(G4x,2)*pow(H,4) - 
   792*G4*G4pp*pow(G4x,2)*pow(H,4) - 864*G4*pow(G4p,2)*G4xx*pow(H,4) + 1152*pow(G4,2)*G4pp*G4xx*pow(H,4) - 216*G3px*pow(G4,2)*G5p*pow(H,4) - 540*G3x*G4*G4p*G5p*pow(H,4) + 
   864*pow(G4,2)*G4ppx*G5p*pow(H,4) + 1944*G4*G4p*G4px*G5p*pow(H,4) - 1332*pow(G4p,2)*G4x*G5p*pow(H,4) + 1296*G4*G4pp*G4x*G5p*pow(H,4) + 612*pow(G4p,2)*pow(G5p,2)*pow(H,4) - 
   504*G4*G4pp*pow(G5p,2)*pow(H,4) + 216*G3x*pow(G4,2)*G5pp*pow(H,4) - 936*pow(G4,2)*G4px*G5pp*pow(H,4) + 252*G4*G4p*G4x*G5pp*pow(H,4) - 540*G4*G4p*G5p*G5pp*pow(H,4) + 
   72*pow(G4,2)*pow(G5pp,2)*pow(H,4) + 72*pow(G4,2)*G4x*G5ppp*pow(H,4) - 72*pow(G4,2)*G5p*G5ppp*pow(H,4) + 720*pow(G4,2)*G4p*G5ppx*pow(H,4) + 36*G2x*pow(G4,2)*G5px*pow(H,4) - 
   72*G3p*pow(G4,2)*G5px*pow(H,4) + 792*G4*pow(G4p,2)*G5px*pow(H,4) - 720*pow(G4,2)*G4pp*G5px*pow(H,4) + 72*G3pp*pow(G4,2)*G5x*pow(H,4) + 54*G2x*G4*G4p*G5x*pow(H,4) - 
   108*G3p*G4*G4p*G5x*pow(H,4) - 126*pow(G4p,3)*G5x*pow(H,4) + 144*G4*G4p*G4pp*G5x*pow(H,4) + 216*pow(G4,2)*G4x*G5px*pow(H,6) - 216*pow(G4,2)*G5p*G5px*pow(H,6) - 
   216*pow(G4,2)*G4px*G5x*pow(H,6) + 756*G4*G4p*G4x*G5x*pow(H,6) - 756*G4*G4p*G5p*G5x*pow(H,6) + 216*pow(G4,2)*G5pp*G5x*pow(H,6) + 
   G2p*(-4*G2pxx*pow(G4,2) + 4*G3ppx*pow(G4,2) - 6*G2xx*G4*G4p + 18*G3px*G4*G4p + 15*G3x*pow(G4p,2) + 12*G3x*G4*G4pp - 24*G4*G4p*G4ppx - 42*pow(G4p,2)*G4px - 24*G4*G4pp*G4px + 8*G2px*G4*G4x - 
      16*G3pp*G4*G4x + 24*G4p*G4pp*G4x - 4*G2px*G4*G5p + 8*G3pp*G4*G5p - 12*G4p*G4pp*G5p + 6*pow(G4p,2)*G5pp - 96*pow(G4,2)*G4pxx*pow(H,2) + 90*G3x*G4*G4x*pow(H,2) - 300*G4*G4px*G4x*pow(H,2) - 
      132*G4p*pow(G4x,2)*pow(H,2) - 144*G4*G4p*G4xx*pow(H,2) - 90*G3x*G4*G5p*pow(H,2) + 324*G4*G4px*G5p*pow(H,2) + 282*G4p*G4x*G5p*pow(H,2) - 150*G4p*pow(G5p,2)*pow(H,2) + 
      48*G4*G4x*G5pp*pow(H,2) - 72*G4*G5p*G5pp*pow(H,2) + 60*pow(G4,2)*G5ppx*pow(H,2) + 126*G4*G4p*G5px*pow(H,2) + 45*pow(G4p,2)*G5x*pow(H,2) + 36*G4*G4pp*G5x*pow(H,2) + 
      126*G4*G4x*G5x*pow(H,4) - 126*G4*G5p*G5x*pow(H,4) + 3*G2x*(G3x*G4 - 2*G4*G4px + 2*G4p*G4x - G4p*G5p + 3*G4*G5x*pow(H,2)) - 
      6*G3p*(G3x*G4 - 2*G4*G4px + 2*G4p*G4x - G4p*G5p + 3*G4*G5x*pow(H,2))) + 3*G2x*G3px*G4*rptot - 6*G3p*G3px*G4*rptot + 6*G3pp*G3x*G4*rptot - 3*G2x*G3x*G4p*rptot + 6*G3p*G3x*G4p*rptot + 
   6*G2pxx*G4*G4p*rptot - 6*G3ppx*G4*G4p*rptot + 6*G2xx*pow(G4p,2)*rptot - 15*G3px*pow(G4p,2)*rptot - 6*G2xx*G4*G4pp*rptot + 6*G3px*G4*G4pp*rptot - 18*G3x*G4p*G4pp*rptot - 6*G2x*G4*G4ppx*rptot + 
   12*G3p*G4*G4ppx*rptot + 18*pow(G4p,2)*G4ppx*rptot - 12*G3pp*G4*G4px*rptot + 36*G4p*G4pp*G4px*rptot + 12*G3pp*G4p*G4x*rptot + 6*G2x*G4pp*G4x*rptot - 12*G3p*G4pp*G4x*rptot - 6*G3pp*G4p*G5p*rptot - 
   3*G2x*G4pp*G5p*rptot + 6*G3p*G4pp*G5p*rptot + 3*G2x*G4p*G5pp*rptot - 6*G3p*G4p*G5pp*rptot + 54*G3x*G4*G4px*pow(H,2)*rptot - 180*G4*pow(G4px,2)*pow(H,2)*rptot + 
   144*G4*G4p*G4pxx*pow(H,2)*rptot - 54*G3px*G4*G4x*pow(H,2)*rptot - 90*G3x*G4p*G4x*pow(H,2)*rptot + 180*G4*G4ppx*G4x*pow(H,2)*rptot + 324*G4p*G4px*G4x*pow(H,2)*rptot + 
   108*G4pp*pow(G4x,2)*pow(H,2)*rptot + 144*pow(G4p,2)*G4xx*pow(H,2)*rptot - 144*G4*G4pp*G4xx*pow(H,2)*rptot + 54*G3px*G4*G5p*pow(H,2)*rptot + 90*G3x*G4p*G5p*pow(H,2)*rptot - 
   180*G4*G4ppx*G5p*pow(H,2)*rptot - 342*G4p*G4px*G5p*pow(H,2)*rptot - 198*G4pp*G4x*G5p*pow(H,2)*rptot + 90*G4pp*pow(G5p,2)*pow(H,2)*rptot - 54*G3x*G4*G5pp*pow(H,2)*rptot + 
   180*G4*G4px*G5pp*pow(H,2)*rptot - 90*G4p*G4x*G5pp*pow(H,2)*rptot + 108*G4p*G5p*G5pp*pow(H,2)*rptot - 90*G4*G4p*G5ppx*pow(H,2)*rptot + 9*G2x*G4*G5px*pow(H,2)*rptot - 
   18*G3p*G4*G5px*pow(H,2)*rptot - 117*pow(G4p,2)*G5px*pow(H,2)*rptot + 90*G4*G4pp*G5px*pow(H,2)*rptot + 18*G3pp*G4*G5x*pow(H,2)*rptot - 9*G2x*G4p*G5x*pow(H,2)*rptot + 
   18*G3p*G4p*G5x*pow(H,2)*rptot - 54*G4p*G4pp*G5x*pow(H,2)*rptot - 18*G4*G4x*G5px*pow(H,4)*rptot + 18*G4*G5p*G5px*pow(H,4)*rptot + 18*G4*G4px*G5x*pow(H,4)*rptot - 
   126*G4p*G4x*G5x*pow(H,4)*rptot + 126*G4p*G5p*G5x*pow(H,4)*rptot - 18*G4*G5pp*G5x*pow(H,4)*rptot - 
   3*G2px*(4*G3pp*pow(G4,2) - 2*G2x*G4*G4p + 4*G3p*G4*G4p + 4*pow(G4p,3) - 12*G4*G4p*G4pp + 12*G3x*pow(G4,2)*pow(H,2) - 56*pow(G4,2)*G4px*pow(H,2) - 100*G4*G4p*G4x*pow(H,2) + 
      84*G4*G4p*G5p*pow(H,2) + 12*pow(G4,2)*G5pp*pow(H,2) + 12*pow(G4,2)*G5x*pow(H,4) + G3x*G4*rptot - 2*G4*G4px*rptot + 2*G4p*G4x*rptot - G4p*G5p*rptot + 3*G4*G5x*pow(H,2)*rptot)) - 
2*pow(a,13)*pow(dphi,4)*(6*G2p*G2xx*G3x*G4 - 18*G2p*G3px*G3x*G4 + 6*G2pp*pow(G3x,2)*G4 + 8*G2px*G2pxx*pow(G4,2) - 8*G2ppx*G2xx*pow(G4,2) - 8*G2pxx*G3pp*pow(G4,2) + 8*G2xx*G3ppp*pow(G4,2) - 
   8*G2px*G3ppx*pow(G4,2) + 8*G3pp*G3ppx*pow(G4,2) + 8*G2ppx*G3px*pow(G4,2) - 8*G3ppp*G3px*pow(G4,2) - 21*G2p*pow(G3x,2)*G4p + 12*G2xx*G3pp*G4*G4p - 24*G2px*G3px*G4*G4p + 12*G3pp*G3px*G4*G4p + 
   24*G2ppx*G3x*G4*G4p - 24*G3ppp*G3x*G4*G4p + 24*G2px*G3x*pow(G4p,2) - 78*G3pp*G3x*pow(G4p,2) - 36*G2px*G3x*G4*G4pp + 48*G3pp*G3x*G4*G4pp + 24*G2pxx*G4*G4p*G4pp - 24*G3ppx*G4*G4p*G4pp + 
   24*G2xx*pow(G4p,2)*G4pp - 60*G3px*pow(G4p,2)*G4pp - 24*G2xx*G4*pow(G4pp,2) + 24*G3px*G4*pow(G4pp,2) - 72*G3x*G4p*pow(G4pp,2) - 24*G2xx*G4*G4p*G4ppp + 24*G3px*G4*G4p*G4ppp + 
   108*G3x*pow(G4p,2)*G4ppp + 24*G2p*G3x*G4*G4ppx + 48*G2px*G4*G4p*G4ppx - 48*G3pp*G4*G4p*G4ppx + 72*pow(G4p,2)*G4pp*G4ppx - 12*G2p*G2xx*G4*G4px + 36*G2p*G3px*G4*G4px - 24*G2pp*G3x*G4*G4px + 
   108*G2p*G3x*G4p*G4px - 48*G2ppx*G4*G4p*G4px + 48*G3ppp*G4*G4p*G4px - 24*G2px*pow(G4p,2)*G4px + 132*G3pp*pow(G4p,2)*G4px + 72*G2px*G4*G4pp*G4px - 96*G3pp*G4*G4pp*G4px + 144*G4p*pow(G4pp,2)*G4px - 
   216*pow(G4p,2)*G4ppp*G4px - 48*G2p*G4*G4ppx*G4px + 24*G2pp*G4*pow(G4px,2) - 132*G2p*G4p*pow(G4px,2) - 16*pow(G2px,2)*G4*G4x + 16*G2p*G2pxx*G4*G4x - 16*G2pp*G2xx*G4*G4x + 48*G2px*G3pp*G4*G4x - 
   32*pow(G3pp,2)*G4*G4x - 16*G2p*G3ppx*G4*G4x + 16*G2pp*G3px*G4*G4x + 12*G2p*G2xx*G4p*G4x - 36*G2p*G3px*G4p*G4x + 24*G2pp*G3x*G4p*G4x + 24*G2ppx*pow(G4p,2)*G4x - 24*G3ppp*pow(G4p,2)*G4x - 
   24*G2p*G3x*G4pp*G4x - 72*G2px*G4p*G4pp*G4x + 96*G3pp*G4p*G4pp*G4x + 48*G2p*G4p*G4ppx*G4x - 48*G2pp*G4p*G4px*G4x + 48*G2p*G4pp*G4px*G4x - 8*G2p*G2px*pow(G4x,2) + 16*G2p*G3pp*pow(G4x,2) + 
   8*pow(G2px,2)*G4*G5p - 8*G2p*G2pxx*G4*G5p + 8*G2pp*G2xx*G4*G5p - 24*G2px*G3pp*G4*G5p + 16*pow(G3pp,2)*G4*G5p + 8*G2p*G3ppx*G4*G5p - 8*G2pp*G3px*G4*G5p - 6*G2p*G2xx*G4p*G5p + 18*G2p*G3px*G4p*G5p - 
   12*G2pp*G3x*G4p*G5p - 12*G2ppx*pow(G4p,2)*G5p + 12*G3ppp*pow(G4p,2)*G5p + 12*G2p*G3x*G4pp*G5p + 36*G2px*G4p*G4pp*G5p - 48*G3pp*G4p*G4pp*G5p - 24*G2p*G4p*G4ppx*G5p + 24*G2pp*G4p*G4px*G5p - 
   24*G2p*G4pp*G4px*G5p + 8*G2p*G2px*G4x*G5p - 16*G2p*G3pp*G4x*G5p - 2*G2p*G2px*pow(G5p,2) + 4*G2p*G3pp*pow(G5p,2) - 12*G2p*G3x*G4p*G5pp - 12*G2px*pow(G4p,2)*G5pp + 12*G3pp*pow(G4p,2)*G5pp + 
   24*G2p*G4p*G4px*G5pp - 72*G2xx*G3px*pow(G4,2)*pow(H,2) + 216*pow(G3px,2)*pow(G4,2)*pow(H,2) + 72*G2pxx*G3x*pow(G4,2)*pow(H,2) - 216*G3ppx*G3x*pow(G4,2)*pow(H,2) - 
   72*G2px*G3xx*pow(G4,2)*pow(H,2) + 144*G3pp*G3xx*pow(G4,2)*pow(H,2) + 36*G2xx*G3x*G4*G4p*pow(H,2) + 36*G3px*G3x*G4*G4p*pow(H,2) + 414*pow(G3x,2)*pow(G4p,2)*pow(H,2) - 
   72*G3pxx*G4*pow(G4p,2)*pow(H,2) - 72*G3xx*pow(G4p,3)*pow(H,2) - 216*pow(G3x,2)*G4*G4pp*pow(H,2) - 72*G3xx*G4*G4p*G4pp*pow(H,2) + 288*G3x*pow(G4,2)*G4pppx*pow(H,2) + 
   144*G2xx*pow(G4,2)*G4ppx*pow(H,2) - 864*G3px*pow(G4,2)*G4ppx*pow(H,2) - 144*G3x*G4*G4p*G4ppx*pow(H,2) + 864*pow(G4,2)*pow(G4ppx,2)*pow(H,2) - 144*G2pxx*pow(G4,2)*G4px*pow(H,2) + 
   576*G3ppx*pow(G4,2)*G4px*pow(H,2) - 48*G2xx*G4*G4p*G4px*pow(H,2) - 600*G3px*G4*G4p*G4px*pow(H,2) - 2268*G3x*pow(G4p,2)*G4px*pow(H,2) + 1440*G3x*G4*G4pp*G4px*pow(H,2) - 
   864*pow(G4,2)*G4pppx*G4px*pow(H,2) + 1296*G4*G4p*G4ppx*G4px*pow(H,2) + 3312*pow(G4p,2)*pow(G4px,2)*pow(H,2) - 2592*G4*G4pp*pow(G4px,2)*pow(H,2) + 
   384*G2px*pow(G4,2)*G4pxx*pow(H,2) - 576*G3pp*pow(G4,2)*G4pxx*pow(H,2) + 144*pow(G4p,3)*G4pxx*pow(H,2) + 1008*G4*G4p*G4pp*G4pxx*pow(H,2) - 48*G2p*pow(G4,2)*G4pxxx*pow(H,2) - 
   108*G2px*G3x*G4*G4x*pow(H,2) + 36*G3pp*G3x*G4*G4x*pow(H,2) + 72*G2p*G3xx*G4*G4x*pow(H,2) + 312*G2pxx*G4*G4p*G4x*pow(H,2) - 24*G3ppx*G4*G4p*G4x*pow(H,2) + 
   192*G2xx*pow(G4p,2)*G4x*pow(H,2) - 84*G3px*pow(G4p,2)*G4x*pow(H,2) - 168*G2xx*G4*G4pp*G4x*pow(H,2) - 336*G3px*G4*G4pp*G4x*pow(H,2) - 2304*G3x*G4p*G4pp*G4x*pow(H,2) + 
   360*G3x*G4*G4ppp*G4x*pow(H,2) - 576*G4*G4p*G4pppx*G4x*pow(H,2) - 1224*pow(G4p,2)*G4ppx*G4x*pow(H,2) + 1296*G4*G4pp*G4ppx*G4x*pow(H,2) + 120*G2px*G4*G4px*G4x*pow(H,2) + 
   360*G3pp*G4*G4px*G4x*pow(H,2) + 5616*G4p*G4pp*G4px*G4x*pow(H,2) - 1008*G4*G4ppp*G4px*G4x*pow(H,2) + 48*G2p*G4*G4pxx*G4x*pow(H,2) + 36*G2p*G3x*pow(G4x,2)*pow(H,2) - 
   600*G2px*G4p*pow(G4x,2)*pow(H,2) + 1464*G3pp*G4p*pow(G4x,2)*pow(H,2) + 432*pow(G4pp,2)*pow(G4x,2)*pow(H,2) - 720*G4p*G4ppp*pow(G4x,2)*pow(H,2) + 
   120*G2p*G4px*pow(G4x,2)*pow(H,2) - 48*G2pp*pow(G4x,3)*pow(H,2) + 288*G2p*G3x*G4*G4xx*pow(H,2) - 192*G2ppx*pow(G4,2)*G4xx*pow(H,2) + 192*G3ppp*pow(G4,2)*G4xx*pow(H,2) + 
   408*G2px*G4*G4p*G4xx*pow(H,2) - 528*G3pp*G4*G4p*G4xx*pow(H,2) + 1944*pow(G4p,2)*G4pp*G4xx*pow(H,2) - 576*G4*pow(G4pp,2)*G4xx*pow(H,2) - 576*G4*G4p*G4ppp*G4xx*pow(H,2) - 
   912*G2p*G4*G4px*G4xx*pow(H,2) - 192*G2pp*G4*G4x*G4xx*pow(H,2) - 384*G2p*G4p*G4x*G4xx*pow(H,2) + 48*G2pp*pow(G4,2)*G4xxx*pow(H,2) - 72*G2p*G4*G4p*G4xxx*pow(H,2) + 
   180*G2px*G3x*G4*G5p*pow(H,2) - 180*G3pp*G3x*G4*G5p*pow(H,2) - 72*G2p*G3xx*G4*G5p*pow(H,2) - 216*G2pxx*G4*G4p*G5p*pow(H,2) - 72*G3ppx*G4*G4p*G5p*pow(H,2) - 
   156*G2xx*pow(G4p,2)*G5p*pow(H,2) + 84*G3px*pow(G4p,2)*G5p*pow(H,2) + 72*G2xx*G4*G4pp*G5p*pow(H,2) + 432*G3px*G4*G4pp*G5p*pow(H,2) + 2160*G3x*G4p*G4pp*G5p*pow(H,2) - 
   360*G3x*G4*G4ppp*G5p*pow(H,2) + 576*G4*G4p*G4pppx*G5p*pow(H,2) + 1080*pow(G4p,2)*G4ppx*G5p*pow(H,2) - 1296*G4*G4pp*G4ppx*G5p*pow(H,2) - 456*G2px*G4*G4px*G5p*pow(H,2) + 
   264*G3pp*G4*G4px*G5p*pow(H,2) - 5472*G4p*G4pp*G4px*G5p*pow(H,2) + 1008*G4*G4ppp*G4px*G5p*pow(H,2) + 144*G2p*G4*G4pxx*G5p*pow(H,2) - 162*G2p*G3x*G4x*G5p*pow(H,2) + 
   48*G2ppx*G4*G4x*G5p*pow(H,2) - 48*G3ppp*G4*G4x*G5p*pow(H,2) + 900*G2px*G4p*G4x*G5p*pow(H,2) - 2364*G3pp*G4p*G4x*G5p*pow(H,2) - 792*pow(G4pp,2)*G4x*G5p*pow(H,2) + 
   1512*G4p*G4ppp*G4x*G5p*pow(H,2) + 60*G2p*G4px*G4x*G5p*pow(H,2) + 144*G2pp*pow(G4x,2)*G5p*pow(H,2) + 528*G2p*G4p*G4xx*G5p*pow(H,2) + 126*G2p*G3x*pow(G5p,2)*pow(H,2) - 
   48*G2ppx*G4*pow(G5p,2)*pow(H,2) + 48*G3ppp*G4*pow(G5p,2)*pow(H,2) - 324*G2px*G4p*pow(G5p,2)*pow(H,2) + 948*G3pp*G4p*pow(G5p,2)*pow(H,2) + 360*pow(G4pp,2)*pow(G5p,2)*pow(H,2) - 
   792*G4p*G4ppp*pow(G5p,2)*pow(H,2) - 192*G2p*G4px*pow(G5p,2)*pow(H,2) - 132*G2pp*G4x*pow(G5p,2)*pow(H,2) + 36*G2pp*pow(G5p,3)*pow(H,2) - 24*G2pxx*pow(G4,2)*G5pp*pow(H,2) + 
   24*G3ppx*pow(G4,2)*G5pp*pow(H,2) + 12*G2xx*G4*G4p*G5pp*pow(H,2) + 492*G3px*G4*G4p*G5pp*pow(H,2) + 270*G3x*pow(G4p,2)*G5pp*pow(H,2) - 288*G3x*G4*G4pp*G5pp*pow(H,2) - 
   1152*G4*G4p*G4ppx*G5pp*pow(H,2) - 1116*pow(G4p,2)*G4px*G5pp*pow(H,2) + 1008*G4*G4pp*G4px*G5pp*pow(H,2) - 96*G2px*G4*G4x*G5pp*pow(H,2) + 96*G3pp*G4*G4x*G5pp*pow(H,2) - 
   504*G4p*G4pp*G4x*G5pp*pow(H,2) - 96*G2p*pow(G4x,2)*G5pp*pow(H,2) + 192*G2p*G4*G4xx*G5pp*pow(H,2) + 168*G2px*G4*G5p*G5pp*pow(H,2) - 192*G3pp*G4*G5p*G5pp*pow(H,2) + 
   648*G4p*G4pp*G5p*G5pp*pow(H,2) + 144*G2p*G4x*G5p*G5pp*pow(H,2) - 36*G2p*pow(G5p,2)*G5pp*pow(H,2) + 180*pow(G4p,2)*pow(G5pp,2)*pow(H,2) + 24*G2xx*pow(G4,2)*G5ppp*pow(H,2) - 
   24*G3px*pow(G4,2)*G5ppp*pow(H,2) - 216*G3x*G4*G4p*G5ppp*pow(H,2) + 576*G4*G4p*G4px*G5ppp*pow(H,2) + 360*pow(G4p,2)*G4x*G5ppp*pow(H,2) - 396*pow(G4p,2)*G5p*G5ppp*pow(H,2) + 
   72*G4*pow(G4p,2)*G5pppx*pow(H,2) - 144*G2px*pow(G4,2)*G5ppx*pow(H,2) + 168*G3pp*pow(G4,2)*G5ppx*pow(H,2) + 36*pow(G4p,3)*G5ppx*pow(H,2) - 504*G4*G4p*G4pp*G5ppx*pow(H,2) - 
   144*G2p*G4*G4x*G5ppx*pow(H,2) + 24*G2p*G4*G5p*G5ppx*pow(H,2) + 24*G2p*pow(G4,2)*G5ppxx*pow(H,2) - 198*G2p*G3x*G4*G5px*pow(H,2) + 120*G2ppx*pow(G4,2)*G5px*pow(H,2) - 
   120*G3ppp*pow(G4,2)*G5px*pow(H,2) - 288*G2px*G4*G4p*G5px*pow(H,2) + 324*G3pp*G4*G4p*G5px*pow(H,2) - 1188*pow(G4p,2)*G4pp*G5px*pow(H,2) + 360*G4*pow(G4pp,2)*G5px*pow(H,2) + 
   360*G4*G4p*G4ppp*G5px*pow(H,2) + 564*G2p*G4*G4px*G5px*pow(H,2) + 144*G2pp*G4*G4x*G5px*pow(H,2) + 132*G2p*G4p*G4x*G5px*pow(H,2) - 24*G2pp*G4*G5p*G5px*pow(H,2) - 
   258*G2p*G4p*G5p*G5px*pow(H,2) - 96*G2p*G4*G5pp*G5px*pow(H,2) - 24*G2pp*pow(G4,2)*G5pxx*pow(H,2) + 60*G2p*G4*G4p*G5pxx*pow(H,2) + 18*G2p*G2xx*G4*G5x*pow(H,2) - 
   6*G2p*G3px*G4*G5x*pow(H,2) - 12*G2pp*G3x*G4*G5x*pow(H,2) - 90*G2p*G3x*G4p*G5x*pow(H,2) + 72*G2ppx*G4*G4p*G5x*pow(H,2) - 72*G3ppp*G4*G4p*G5x*pow(H,2) + 120*G2px*pow(G4p,2)*G5x*pow(H,2) - 
   330*G3pp*pow(G4p,2)*G5x*pow(H,2) - 108*G2px*G4*G4pp*G5x*pow(H,2) + 144*G3pp*G4*G4pp*G5x*pow(H,2) - 216*G4p*pow(G4pp,2)*G5x*pow(H,2) + 324*pow(G4p,2)*G4ppp*G5x*pow(H,2) - 
   72*G2p*G4*G4ppx*G5x*pow(H,2) + 72*G2pp*G4*G4px*G5x*pow(H,2) + 168*G2p*G4p*G4px*G5x*pow(H,2) + 120*G2pp*G4p*G4x*G5x*pow(H,2) - 120*G2p*G4pp*G4x*G5x*pow(H,2) - 84*G2pp*G4p*G5p*G5x*pow(H,2) + 
   84*G2p*G4pp*G5p*G5x*pow(H,2) + 12*G2p*G4p*G5pp*G5x*pow(H,2) - 24*G2pp*G4*G4p*G5xx*pow(H,2) + 30*G2p*pow(G4p,2)*G5xx*pow(H,2) + 24*G2p*G4*G4pp*G5xx*pow(H,2) - 
   432*G3xx*pow(G4,2)*G4px*pow(H,4) + 864*G3x*pow(G4,2)*G4pxx*pow(H,4) + 288*pow(G4,2)*G4px*G4pxx*pow(H,4) - 576*pow(G4,2)*G4p*G4pxxx*pow(H,4) + 432*G3pxx*pow(G4,2)*G4x*pow(H,4) + 
   432*G3xx*G4*G4p*G4x*pow(H,4) - 1152*pow(G4,2)*G4ppxx*G4x*pow(H,4) + 2160*G3x*G4*G4px*G4x*pow(H,4) - 6048*G4*pow(G4px,2)*G4x*pow(H,4) + 4032*G4*G4p*G4pxx*G4x*pow(H,4) - 
   1512*G3px*G4*pow(G4x,2)*pow(H,4) - 3024*G3x*G4p*pow(G4x,2)*pow(H,4) + 5040*G4*G4ppx*pow(G4x,2)*pow(H,4) + 10512*G4p*G4px*pow(G4x,2)*pow(H,4) + 3024*G4pp*pow(G4x,3)*pow(H,4) - 
   864*G3px*pow(G4,2)*G4xx*pow(H,4) + 1728*G3x*G4*G4p*G4xx*pow(H,4) + 864*pow(G4,2)*G4ppx*G4xx*pow(H,4) - 5328*G4*G4p*G4px*G4xx*pow(H,4) + 5904*pow(G4p,2)*G4x*G4xx*pow(H,4) - 
   5040*G4*G4pp*G4x*G4xx*pow(H,4) - 432*G4*pow(G4p,2)*G4xxx*pow(H,4) + 576*pow(G4,2)*G4pp*G4xxx*pow(H,4) - 432*G3pxx*pow(G4,2)*G5p*pow(H,4) - 432*G3xx*G4*G4p*G5p*pow(H,4) + 
   1152*pow(G4,2)*G4ppxx*G5p*pow(H,4) - 1728*G3x*G4*G4px*G5p*pow(H,4) + 4320*G4*pow(G4px,2)*G5p*pow(H,4) - 1728*G4*G4p*G4pxx*G5p*pow(H,4) + 2592*G3px*G4*G4x*G5p*pow(H,4) + 
   5508*G3x*G4p*G4x*G5p*pow(H,4) - 8352*G4*G4ppx*G4x*G5p*pow(H,4) - 19368*G4p*G4px*G4x*G5p*pow(H,4) - 8280*G4pp*pow(G4x,2)*G5p*pow(H,4) - 5040*pow(G4p,2)*G4xx*G5p*pow(H,4) + 
   2736*G4*G4pp*G4xx*G5p*pow(H,4) - 1080*G3px*G4*pow(G5p,2)*pow(H,4) - 2484*G3x*G4p*pow(G5p,2)*pow(H,4) + 3312*G4*G4ppx*pow(G5p,2)*pow(H,4) + 8712*G4p*G4px*pow(G5p,2)*pow(H,4) + 
   7632*G4pp*G4x*pow(G5p,2)*pow(H,4) - 2376*G4pp*pow(G5p,3)*pow(H,4) + 432*G3xx*pow(G4,2)*G5pp*pow(H,4) - 1728*pow(G4,2)*G4pxx*G5pp*pow(H,4) - 2700*G3x*G4*G4x*G5pp*pow(H,4) + 
   7272*G4*G4px*G4x*G5pp*pow(H,4) - 2880*G4p*pow(G4x,2)*G5pp*pow(H,4) + 1872*G4*G4p*G4xx*G5pp*pow(H,4) + 2268*G3x*G4*G5p*G5pp*pow(H,4) - 5400*G4*G4px*G5p*G5pp*pow(H,4) + 
   5508*G4p*G4x*G5p*G5pp*pow(H,4) - 2484*G4p*pow(G5p,2)*G5pp*pow(H,4) + 288*G4*G4x*pow(G5pp,2)*pow(H,4) - 432*G4*G5p*pow(G5pp,2)*pow(H,4) + 288*G4*pow(G4x,2)*G5ppp*pow(H,4) + 
   576*pow(G4,2)*G4xx*G5ppp*pow(H,4) - 720*G4*G4x*G5p*G5ppp*pow(H,4) + 432*G4*pow(G5p,2)*G5ppp*pow(H,4) + 144*pow(G4,2)*G4x*G5pppx*pow(H,4) - 144*pow(G4,2)*G5p*G5pppx*pow(H,4) - 
   744*G3x*pow(G4,2)*G5ppx*pow(H,4) + 1008*pow(G4,2)*G4px*G5ppx*pow(H,4) - 3120*G4*G4p*G4x*G5ppx*pow(H,4) + 1680*G4*G4p*G5p*G5ppx*pow(H,4) + 504*pow(G4,2)*G5pp*G5ppx*pow(H,4) + 
   288*pow(G4,2)*G4p*G5ppxx*pow(H,4) - 72*G2xx*pow(G4,2)*G5px*pow(H,4) + 960*G3px*pow(G4,2)*G5px*pow(H,4) - 1620*G3x*G4*G4p*G5px*pow(H,4) - 1440*pow(G4,2)*G4ppx*G5px*pow(H,4) + 
   3624*G4*G4p*G4px*G5px*pow(H,4) - 4236*pow(G4p,2)*G4x*G5px*pow(H,4) + 2832*G4*G4pp*G4x*G5px*pow(H,4) + 3444*pow(G4p,2)*G5p*G5px*pow(H,4) - 1392*G4*G4pp*G5p*G5px*pow(H,4) - 
   348*G4*G4p*G5pp*G5px*pow(H,4) - 360*pow(G4,2)*G5ppp*G5px*pow(H,4) + 264*G4*pow(G4p,2)*G5pxx*pow(H,4) - 288*pow(G4,2)*G4pp*G5pxx*pow(H,4) + 72*G2pxx*pow(G4,2)*G5x*pow(H,4) - 
   216*G3ppx*pow(G4,2)*G5x*pow(H,4) + 108*G2xx*G4*G4p*G5x*pow(H,4) + 828*G3px*G4*G4p*G5x*pow(H,4) + 1260*G3x*pow(G4p,2)*G5x*pow(H,4) - 864*G3x*G4*G4pp*G5x*pow(H,4) + 
   288*pow(G4,2)*G4pppx*G5x*pow(H,4) - 2592*G4*G4p*G4ppx*G5x*pow(H,4) - 4140*pow(G4p,2)*G4px*G5x*pow(H,4) + 2448*G4*G4pp*G4px*G5x*pow(H,4) - 468*G2px*G4*G4x*G5x*pow(H,4) + 
   684*G3pp*G4*G4x*G5x*pow(H,4) - 5256*G4p*G4pp*G4x*G5x*pow(H,4) + 504*G4*G4ppp*G4x*G5x*pow(H,4) + 324*G2p*pow(G4x,2)*G5x*pow(H,4) + 576*G2p*G4*G4xx*G5x*pow(H,4) + 
   540*G2px*G4*G5p*G5x*pow(H,4) - 828*G3pp*G4*G5p*G5x*pow(H,4) + 5112*G4p*G4pp*G5p*G5x*pow(H,4) - 504*G4*G4ppp*G5p*G5x*pow(H,4) - 774*G2p*G4x*G5p*G5x*pow(H,4) + 
   450*G2p*pow(G5p,2)*G5x*pow(H,4) + 738*pow(G4p,2)*G5pp*G5x*pow(H,4) - 360*G4*G4p*G5ppp*G5x*pow(H,4) - 402*G2p*G4*G5px*G5x*pow(H,4) + 6*G2pp*G4*pow(G5x,2)*pow(H,4) - 
   153*G2p*G4p*pow(G5x,2)*pow(H,4) - 144*G2px*pow(G4,2)*G5xx*pow(H,4) + 288*G3pp*pow(G4,2)*G5xx*pow(H,4) - 204*pow(G4p,3)*G5xx*pow(H,4) - 24*G4*G4p*G4pp*G5xx*pow(H,4) + 
   204*G2p*G4*G4x*G5xx*pow(H,4) - 204*G2p*G4*G5p*G5xx*pow(H,4) + 360*G4*pow(G4x,2)*G5px*pow(H,6) - 864*pow(G4,2)*G4xx*G5px*pow(H,6) - 1152*G4*G4x*G5p*G5px*pow(H,6) + 
   792*G4*pow(G5p,2)*G5px*pow(H,6) + 744*pow(G4,2)*pow(G5px,2)*pow(H,6) + 864*pow(G4,2)*G4x*G5pxx*pow(H,6) - 864*pow(G4,2)*G5p*G5pxx*pow(H,6) + 864*pow(G4,2)*G4pxx*G5x*pow(H,6) + 
   1296*G4*G4px*G4x*G5x*pow(H,6) - 3456*G4p*pow(G4x,2)*G5x*pow(H,6) + 3456*G4*G4p*G4xx*G5x*pow(H,6) - 864*G4*G4px*G5p*G5x*pow(H,6) + 6156*G4p*G4x*G5p*G5x*pow(H,6) - 
   2700*G4p*pow(G5p,2)*G5x*pow(H,6) - 2052*G4*G4x*G5pp*G5x*pow(H,6) + 1620*G4*G5p*G5pp*G5x*pow(H,6) - 744*pow(G4,2)*G5ppx*G5x*pow(H,6) - 2412*G4*G4p*G5px*G5x*pow(H,6) + 
   486*pow(G4p,2)*pow(G5x,2)*pow(H,6) - 360*G4*G4pp*pow(G5x,2)*pow(H,6) - 864*pow(G4,2)*G4px*G5xx*pow(H,6) + 1224*G4*G4p*G4x*G5xx*pow(H,6) - 1224*G4*G4p*G5p*G5xx*pow(H,6) + 
   864*pow(G4,2)*G5pp*G5xx*pow(H,6) + 6*pow(G2x,2)*(G3px*G4 - G3x*G4p - 2*G4*G4ppx + 2*G4pp*G4x - G4pp*G5p + G4p*G5pp + 3*G4*G5px*pow(H,2) - 3*G4p*G5x*pow(H,2)) + 
   24*pow(G3p,2)*(G3px*G4 - G3x*G4p - 2*G4*G4ppx + 2*G4pp*G4x - G4pp*G5p + G4p*G5pp + 3*G4*G5px*pow(H,2) - 3*G4p*G5x*pow(H,2)) + 6*G2xx*G3px*G4*rptot - 6*pow(G3px,2)*G4*rptot - 
   6*G2pxx*G3x*G4*rptot + 6*G3ppx*G3x*G4*rptot - 6*G2xx*G3x*G4p*rptot + 24*G3px*G3x*G4p*rptot + 9*pow(G3x,2)*G4pp*rptot - 12*G2xx*G4*G4ppx*rptot + 12*G3px*G4*G4ppx*rptot - 36*G3x*G4p*G4ppx*rptot + 
   12*G2pxx*G4*G4px*rptot - 12*G3ppx*G4*G4px*rptot - 36*G3px*G4p*G4px*rptot - 36*G3x*G4pp*G4px*rptot + 72*G4p*G4ppx*G4px*rptot + 36*G4pp*pow(G4px,2)*rptot + 6*G2px*G3x*G4x*rptot - 
   12*G3pp*G3x*G4x*rptot - 12*G2pxx*G4p*G4x*rptot + 12*G3ppx*G4p*G4x*rptot + 12*G2xx*G4pp*G4x*rptot - 12*G3px*G4pp*G4x*rptot - 12*G2px*G4px*G4x*rptot + 24*G3pp*G4px*G4x*rptot - 3*G2px*G3x*G5p*rptot + 
   6*G3pp*G3x*G5p*rptot + 6*G2pxx*G4p*G5p*rptot - 6*G3ppx*G4p*G5p*rptot - 6*G2xx*G4pp*G5p*rptot + 6*G3px*G4pp*G5p*rptot + 6*G2px*G4px*G5p*rptot - 12*G3pp*G4px*G5p*rptot + 6*G2xx*G4p*G5pp*rptot - 
   6*G3px*G4p*G5pp*rptot + 72*G3xx*G4*G4px*pow(H,2)*rptot - 288*G4*G4px*G4pxx*pow(H,2)*rptot + 72*G4*G4p*G4pxxx*pow(H,2)*rptot - 72*G3pxx*G4*G4x*pow(H,2)*rptot - 
   72*G3xx*G4p*G4x*pow(H,2)*rptot + 144*G4*G4ppxx*G4x*pow(H,2)*rptot - 72*G3x*G4px*G4x*pow(H,2)*rptot + 144*pow(G4px,2)*G4x*pow(H,2)*rptot + 144*G4p*G4pxx*G4x*pow(H,2)*rptot + 
   36*G3px*pow(G4x,2)*pow(H,2)*rptot - 216*G4ppx*pow(G4x,2)*pow(H,2)*rptot - 288*G3x*G4p*G4xx*pow(H,2)*rptot + 144*G4*G4ppx*G4xx*pow(H,2)*rptot + 720*G4p*G4px*G4xx*pow(H,2)*rptot + 
   576*G4pp*G4x*G4xx*pow(H,2)*rptot + 72*pow(G4p,2)*G4xxx*pow(H,2)*rptot - 72*G4*G4pp*G4xxx*pow(H,2)*rptot + 72*G3pxx*G4*G5p*pow(H,2)*rptot + 72*G3xx*G4p*G5p*pow(H,2)*rptot - 
   144*G4*G4ppxx*G5p*pow(H,2)*rptot + 18*G3x*G4px*G5p*pow(H,2)*rptot + 36*pow(G4px,2)*G5p*pow(H,2)*rptot - 288*G4p*G4pxx*G5p*pow(H,2)*rptot - 18*G3px*G4x*G5p*pow(H,2)*rptot + 
   252*G4ppx*G4x*G5p*pow(H,2)*rptot - 432*G4pp*G4xx*G5p*pow(H,2)*rptot - 18*G3px*pow(G5p,2)*pow(H,2)*rptot - 36*G4ppx*pow(G5p,2)*pow(H,2)*rptot - 72*G3xx*G4*G5pp*pow(H,2)*rptot + 
   144*G4*G4pxx*G5pp*pow(H,2)*rptot + 162*G3x*G4x*G5pp*pow(H,2)*rptot - 396*G4px*G4x*G5pp*pow(H,2)*rptot - 144*G4p*G4xx*G5pp*pow(H,2)*rptot - 108*G3x*G5p*G5pp*pow(H,2)*rptot + 
   216*G4px*G5p*G5pp*pow(H,2)*rptot + 18*G3x*G4*G5ppx*pow(H,2)*rptot + 36*G4*G4px*G5ppx*pow(H,2)*rptot + 36*G4p*G4x*G5ppx*pow(H,2)*rptot + 54*G4p*G5p*G5ppx*pow(H,2)*rptot - 
   36*G4*G4p*G5ppxx*pow(H,2)*rptot + 18*G2xx*G4*G5px*pow(H,2)*rptot - 36*G3px*G4*G5px*pow(H,2)*rptot + 180*G3x*G4p*G5px*pow(H,2)*rptot - 36*G4*G4ppx*G5px*pow(H,2)*rptot - 
   360*G4p*G4px*G5px*pow(H,2)*rptot - 324*G4pp*G4x*G5px*pow(H,2)*rptot + 234*G4pp*G5p*G5px*pow(H,2)*rptot + 54*G4p*G5pp*G5px*pow(H,2)*rptot - 54*pow(G4p,2)*G5pxx*pow(H,2)*rptot + 
   36*G4*G4pp*G5pxx*pow(H,2)*rptot - 18*G2pxx*G4*G5x*pow(H,2)*rptot + 18*G3ppx*G4*G5x*pow(H,2)*rptot - 18*G2xx*G4p*G5x*pow(H,2)*rptot + 36*G3px*G4p*G5x*pow(H,2)*rptot + 
   90*G3x*G4pp*G5x*pow(H,2)*rptot - 216*G4pp*G4px*G5x*pow(H,2)*rptot + 30*G2px*G4x*G5x*pow(H,2)*rptot - 60*G3pp*G4x*G5x*pow(H,2)*rptot - 21*G2px*G5p*G5x*pow(H,2)*rptot + 
   42*G3pp*G5p*G5x*pow(H,2)*rptot - 6*G2px*G4*G5xx*pow(H,2)*rptot + 12*G3pp*G4*G5xx*pow(H,2)*rptot - 36*G4p*G4pp*G5xx*pow(H,2)*rptot - 108*pow(G4x,2)*G5px*pow(H,4)*rptot + 
   288*G4*G4xx*G5px*pow(H,4)*rptot + 234*G4x*G5p*G5px*pow(H,4)*rptot - 126*pow(G5p,2)*G5px*pow(H,4)*rptot - 198*G4*pow(G5px,2)*pow(H,4)*rptot - 132*G4*G4x*G5pxx*pow(H,4)*rptot + 
   132*G4*G5p*G5pxx*pow(H,4)*rptot - 288*G4*G4pxx*G5x*pow(H,4)*rptot - 216*G4px*G4x*G5x*pow(H,4)*rptot - 576*G4p*G4xx*G5x*pow(H,4)*rptot + 198*G4px*G5p*G5x*pow(H,4)*rptot + 
   342*G4x*G5pp*G5x*pow(H,4)*rptot - 324*G5p*G5pp*G5x*pow(H,4)*rptot + 198*G4*G5ppx*G5x*pow(H,4)*rptot + 432*G4p*G5px*G5x*pow(H,4)*rptot + 117*G4pp*pow(G5x,2)*pow(H,4)*rptot + 
   132*G4*G4px*G5xx*pow(H,4)*rptot - 204*G4p*G4x*G5xx*pow(H,4)*rptot + 204*G4p*G5p*G5xx*pow(H,4)*rptot - 132*G4*G5pp*G5xx*pow(H,4)*rptot + 
   G2x*(6*G3pp*G3x*G4 + 12*G2pxx*G4*G4p - 12*G3ppx*G4*G4p + 12*G2xx*pow(G4p,2) - 30*G3px*pow(G4p,2) - 12*G2xx*G4*G4pp + 24*G3px*G4*G4pp - 48*G3x*G4p*G4pp + 12*G3x*G4*G4ppp + 36*pow(G4p,2)*G4ppx - 
      24*G4*G4pp*G4ppx - 12*G3pp*G4*G4px + 72*G4p*G4pp*G4px - 24*G4*G4ppp*G4px - 6*G2p*G3x*G4x + 16*G2ppx*G4*G4x - 16*G3ppp*G4*G4x + 12*G3pp*G4p*G4x + 24*pow(G4pp,2)*G4x + 24*G4p*G4ppp*G4x + 
      12*G2p*G4px*G4x + 8*G2pp*pow(G4x,2) + 3*G2p*G3x*G5p - 8*G2ppx*G4*G5p + 8*G3ppp*G4*G5p - 6*G3pp*G4p*G5p - 12*pow(G4pp,2)*G5p - 12*G4p*G4ppp*G5p - 6*G2p*G4px*G5p - 8*G2pp*G4x*G5p + 
      2*G2pp*pow(G5p,2) + 12*G4p*G4pp*G5pp + 72*G3pxx*pow(G4,2)*pow(H,2) - 192*pow(G4,2)*G4ppxx*pow(H,2) + 288*G3x*G4*G4px*pow(H,2) - 816*G4*pow(G4px,2)*pow(H,2) - 
      48*G4*G4p*G4pxx*pow(H,2) - 72*G3px*G4*G4x*pow(H,2) - 468*G3x*G4p*G4x*pow(H,2) + 480*G4*G4ppx*G4x*pow(H,2) + 1176*G4p*G4px*G4x*pow(H,2) + 576*G4pp*pow(G4x,2)*pow(H,2) + 
      696*pow(G4p,2)*G4xx*pow(H,2) - 552*G4*G4pp*G4xx*pow(H,2) + 450*G3x*G4p*G5p*pow(H,2) - 192*G4*G4ppx*G5p*pow(H,2) - 1044*G4p*G4px*G5p*pow(H,2) - 888*G4pp*G4x*G5p*pow(H,2) + 
      336*G4pp*pow(G5p,2)*pow(H,2) - 306*G3x*G4*G5pp*pow(H,2) + 900*G4*G4px*G5pp*pow(H,2) - 324*G4p*G4x*G5pp*pow(H,2) + 210*G4p*G5p*G5pp*pow(H,2) - 48*G4*pow(G5pp,2)*pow(H,2) - 
      24*G4*G5p*G5ppp*pow(H,2) + 24*pow(G4,2)*G5pppx*pow(H,2) - 402*pow(G4p,2)*G5px*pow(H,2) + 360*G4*G4pp*G5px*pow(H,2) + 18*G3pp*G4*G5x*pow(H,2) - 204*G4p*G4pp*G5x*pow(H,2) + 
      36*G4*G4ppp*G5x*pow(H,2) - 30*G2p*G4x*G5x*pow(H,2) + 21*G2p*G5p*G5x*pow(H,2) + 6*G2p*G4*G5xx*pow(H,2) + 312*G4*G4x*G5px*pow(H,4) - 384*G4*G5p*G5px*pow(H,4) + 
      144*pow(G4,2)*G5pxx*pow(H,4) + 432*G4*G4px*G5x*pow(H,4) - 1044*G4p*G4x*G5x*pow(H,4) + 990*G4p*G5p*G5x*pow(H,4) - 486*G4*G5pp*G5x*pow(H,4) + 36*G4*G4p*G5xx*pow(H,4) - 
      6*G2px*(G3x*G4 - 2*G4*G4px + 2*G4p*G4x - G4p*G5p + 3*G4*G5x*pow(H,2)) - 24*G3p*(G3px*G4 - G3x*G4p - 2*G4*G4ppx + 2*G4pp*G4x - G4pp*G5p + G4p*G5pp + 3*G4*G5px*pow(H,2) - 3*G4p*G5x*pow(H,2)) + 
      6*G3x*G4px*rptot - 12*pow(G4px,2)*rptot - 6*G3px*G4x*rptot + 12*G4ppx*G4x*rptot + 3*G3px*G5p*rptot - 6*G4ppx*G5p*rptot - 3*G3x*G5pp*rptot + 6*G4px*G5pp*rptot - 6*G4x*G5px*pow(H,2)*rptot - 
      3*G5p*G5px*pow(H,2)*rptot + 6*G4*G5pxx*pow(H,2)*rptot + 6*G4px*G5x*pow(H,2)*rptot + 3*G5pp*G5x*pow(H,2)*rptot - 6*G4p*G5xx*pow(H,2)*rptot) + 
   2*G3p*(-12*G2pxx*G4*G4p + 12*G3ppx*G4*G4p - 12*G2xx*pow(G4p,2) + 30*G3px*pow(G4p,2) + 12*G2xx*G4*G4pp - 24*G3px*G4*G4pp + 48*G3x*G4p*G4pp - 12*G3x*G4*G4ppp - 36*pow(G4p,2)*G4ppx + 
      24*G4*G4pp*G4ppx - 72*G4p*G4pp*G4px + 24*G4*G4ppp*G4px + 6*G2p*G3x*G4x - 16*G2ppx*G4*G4x + 16*G3ppp*G4*G4x - 24*pow(G4pp,2)*G4x - 24*G4p*G4ppp*G4x - 12*G2p*G4px*G4x - 8*G2pp*pow(G4x,2) - 
      3*G2p*G3x*G5p + 8*G2ppx*G4*G5p - 8*G3ppp*G4*G5p + 12*pow(G4pp,2)*G5p + 12*G4p*G4ppp*G5p + 6*G2p*G4px*G5p + 8*G2pp*G4x*G5p - 2*G2pp*pow(G5p,2) - 12*G4p*G4pp*G5pp - 
      72*G3pxx*pow(G4,2)*pow(H,2) + 192*pow(G4,2)*G4ppxx*pow(H,2) - 288*G3x*G4*G4px*pow(H,2) + 816*G4*pow(G4px,2)*pow(H,2) + 48*G4*G4p*G4pxx*pow(H,2) + 72*G3px*G4*G4x*pow(H,2) + 
      468*G3x*G4p*G4x*pow(H,2) - 480*G4*G4ppx*G4x*pow(H,2) - 1176*G4p*G4px*G4x*pow(H,2) - 576*G4pp*pow(G4x,2)*pow(H,2) - 696*pow(G4p,2)*G4xx*pow(H,2) + 552*G4*G4pp*G4xx*pow(H,2) - 
      450*G3x*G4p*G5p*pow(H,2) + 192*G4*G4ppx*G5p*pow(H,2) + 1044*G4p*G4px*G5p*pow(H,2) + 888*G4pp*G4x*G5p*pow(H,2) - 336*G4pp*pow(G5p,2)*pow(H,2) + 306*G3x*G4*G5pp*pow(H,2) - 
      900*G4*G4px*G5pp*pow(H,2) + 324*G4p*G4x*G5pp*pow(H,2) - 210*G4p*G5p*G5pp*pow(H,2) + 48*G4*pow(G5pp,2)*pow(H,2) + 24*G4*G5p*G5ppp*pow(H,2) - 24*pow(G4,2)*G5pppx*pow(H,2) + 
      402*pow(G4p,2)*G5px*pow(H,2) - 360*G4*G4pp*G5px*pow(H,2) + 204*G4p*G4pp*G5x*pow(H,2) - 36*G4*G4ppp*G5x*pow(H,2) + 30*G2p*G4x*G5x*pow(H,2) - 21*G2p*G5p*G5x*pow(H,2) - 
      6*G2p*G4*G5xx*pow(H,2) - 312*G4*G4x*G5px*pow(H,4) + 384*G4*G5p*G5px*pow(H,4) - 144*pow(G4,2)*G5pxx*pow(H,4) - 432*G4*G4px*G5x*pow(H,4) + 1044*G4p*G4x*G5x*pow(H,4) - 
      990*G4p*G5p*G5x*pow(H,4) + 486*G4*G5pp*G5x*pow(H,4) - 36*G4*G4p*G5xx*pow(H,4) + 6*G2px*(G3x*G4 - 2*G4*G4px + 2*G4p*G4x - G4p*G5p + 3*G4*G5x*pow(H,2)) - 
      6*G3pp*(G3x*G4 - 2*G4*G4px + 2*G4p*G4x - G4p*G5p + 3*G4*G5x*pow(H,2)) - 6*G3x*G4px*rptot + 12*pow(G4px,2)*rptot + 6*G3px*G4x*rptot - 12*G4ppx*G4x*rptot - 3*G3px*G5p*rptot + 6*G4ppx*G5p*rptot + 
      3*G3x*G5pp*rptot - 6*G4px*G5pp*rptot + 6*G4x*G5px*pow(H,2)*rptot + 3*G5p*G5px*pow(H,2)*rptot - 6*G4*G5pxx*pow(H,2)*rptot - 6*G4px*G5x*pow(H,2)*rptot - 3*G5pp*G5x*pow(H,2)*rptot + 
      6*G4p*G5xx*pow(H,2)*rptot)) - pow(a,11)*pow(dphi,6)*(-24*G2xx*G3p*G3px*G4 + 24*G3p*pow(G3px,2)*G4 + 24*G2pxx*G3p*G3x*G4 - 12*G2xx*G3pp*G3x*G4 - 24*G3p*G3ppx*G3x*G4 + 24*G2px*G3px*G3x*G4 - 
   12*G3pp*G3px*G3x*G4 - 12*G2ppx*pow(G3x,2)*G4 + 12*G3ppp*pow(G3x,2)*G4 + 24*G2xx*G3p*G3x*G4p - 96*G3p*G3px*G3x*G4p - 12*G2px*pow(G3x,2)*G4p + 66*G3pp*pow(G3x,2)*G4p - 36*G3p*pow(G3x,2)*G4pp + 
   24*G2xx*G3px*G4*G4pp - 24*pow(G3px,2)*G4*G4pp - 24*G2pxx*G3x*G4*G4pp + 24*G3ppx*G3x*G4*G4pp - 24*G2xx*G3x*G4p*G4pp + 96*G3px*G3x*G4p*G4pp + 36*pow(G3x,2)*pow(G4pp,2) + 24*G2xx*G3x*G4*G4ppp - 
   24*G3px*G3x*G4*G4ppp - 108*pow(G3x,2)*G4p*G4ppp + 48*G2xx*G3p*G4*G4ppx - 48*G3p*G3px*G4*G4ppx - 48*G2px*G3x*G4*G4ppx + 48*G3pp*G3x*G4*G4ppx + 144*G3p*G3x*G4p*G4ppx - 48*G2xx*G4*G4pp*G4ppx + 
   48*G3px*G4*G4pp*G4ppx - 144*G3x*G4p*G4pp*G4ppx + 48*pow(G3p,2)*G3x*G4px - 48*G2pxx*G3p*G4*G4px + 24*G2xx*G3pp*G4*G4px + 48*G3p*G3ppx*G4*G4px - 48*G2px*G3px*G4*G4px + 24*G3pp*G3px*G4*G4px + 
   48*G2ppx*G3x*G4*G4px - 48*G3ppp*G3x*G4*G4px + 144*G3p*G3px*G4p*G4px - 216*G3pp*G3x*G4p*G4px + 96*G3p*G3x*G4pp*G4px + 48*G2pxx*G4*G4pp*G4px - 48*G3ppx*G4*G4pp*G4px - 144*G3px*G4p*G4pp*G4px - 
   144*G3x*pow(G4pp,2)*G4px - 48*G2xx*G4*G4ppp*G4px + 48*G3px*G4*G4ppp*G4px + 432*G3x*G4p*G4ppp*G4px + 96*G2px*G4*G4ppx*G4px - 96*G3pp*G4*G4ppx*G4px - 288*G3p*G4p*G4ppx*G4px + 288*G4p*G4pp*G4ppx*G4px - 
   96*pow(G3p,2)*pow(G4px,2) - 48*G2ppx*G4*pow(G4px,2) + 48*G3ppp*G4*pow(G4px,2) + 48*G2px*G4p*pow(G4px,2) + 168*G3pp*G4p*pow(G4px,2) - 48*G3p*G4pp*pow(G4px,2) + 
   144*pow(G4pp,2)*pow(G4px,2) - 432*G4p*G4ppp*pow(G4px,2) - 48*pow(G3p,2)*G3px*G4x - 24*G2px*G3p*G3x*G4x + 24*G3p*G3pp*G3x*G4x - 12*G2pp*pow(G3x,2)*G4x - 32*G2px*G2pxx*G4*G4x + 
   32*G2ppx*G2xx*G4*G4x + 32*G2pxx*G3pp*G4*G4x - 32*G2xx*G3ppp*G4*G4x + 32*G2px*G3ppx*G4*G4x - 32*G3pp*G3ppx*G4*G4x - 32*G2ppx*G3px*G4*G4x + 32*G3ppp*G3px*G4*G4x + 48*G2pxx*G3p*G4p*G4x - 
   24*G2xx*G3pp*G4p*G4x - 48*G3p*G3ppx*G4p*G4x + 48*G2px*G3px*G4p*G4x - 24*G3pp*G3px*G4p*G4x - 48*G2ppx*G3x*G4p*G4x + 48*G3ppp*G3x*G4p*G4x - 48*G2xx*G3p*G4pp*G4x + 96*G3p*G3px*G4pp*G4x + 
   72*G2px*G3x*G4pp*G4x - 96*G3pp*G3x*G4pp*G4x - 48*G2pxx*G4p*G4pp*G4x + 48*G3ppx*G4p*G4pp*G4x + 48*G2xx*pow(G4pp,2)*G4x - 48*G3px*pow(G4pp,2)*G4x + 48*G3p*G3x*G4ppp*G4x + 48*G2xx*G4p*G4ppp*G4x - 
   48*G3px*G4p*G4ppp*G4x + 96*pow(G3p,2)*G4ppx*G4x - 96*G2px*G4p*G4ppx*G4x + 96*G3pp*G4p*G4ppx*G4x - 96*G3p*G4pp*G4ppx*G4x + 48*G2px*G3p*G4px*G4x - 48*G3p*G3pp*G4px*G4x + 48*G2pp*G3x*G4px*G4x + 
   96*G2ppx*G4p*G4px*G4x - 96*G3ppp*G4p*G4px*G4x - 144*G2px*G4pp*G4px*G4x + 192*G3pp*G4pp*G4px*G4x - 96*G3p*G4ppp*G4px*G4x - 48*G2pp*pow(G4px,2)*G4x + 16*pow(G2px,2)*pow(G4x,2) + 
   16*G2pp*G2xx*pow(G4x,2) + 32*G2ppx*G3p*pow(G4x,2) - 48*G2px*G3pp*pow(G4x,2) + 32*pow(G3pp,2)*pow(G4x,2) - 32*G3p*G3ppp*pow(G4x,2) - 16*G2pp*G3px*pow(G4x,2) + 24*pow(G3p,2)*G3px*G5p + 
   12*G2px*G3p*G3x*G5p - 12*G3p*G3pp*G3x*G5p + 6*G2pp*pow(G3x,2)*G5p + 16*G2px*G2pxx*G4*G5p - 16*G2ppx*G2xx*G4*G5p - 16*G2pxx*G3pp*G4*G5p + 16*G2xx*G3ppp*G4*G5p - 16*G2px*G3ppx*G4*G5p + 
   16*G3pp*G3ppx*G4*G5p + 16*G2ppx*G3px*G4*G5p - 16*G3ppp*G3px*G4*G5p - 24*G2pxx*G3p*G4p*G5p + 12*G2xx*G3pp*G4p*G5p + 24*G3p*G3ppx*G4p*G5p - 24*G2px*G3px*G4p*G5p + 12*G3pp*G3px*G4p*G5p + 
   24*G2ppx*G3x*G4p*G5p - 24*G3ppp*G3x*G4p*G5p + 24*G2xx*G3p*G4pp*G5p - 48*G3p*G3px*G4pp*G5p - 36*G2px*G3x*G4pp*G5p + 48*G3pp*G3x*G4pp*G5p + 24*G2pxx*G4p*G4pp*G5p - 24*G3ppx*G4p*G4pp*G5p - 
   24*G2xx*pow(G4pp,2)*G5p + 24*G3px*pow(G4pp,2)*G5p - 24*G3p*G3x*G4ppp*G5p - 24*G2xx*G4p*G4ppp*G5p + 24*G3px*G4p*G4ppp*G5p - 48*pow(G3p,2)*G4ppx*G5p + 48*G2px*G4p*G4ppx*G5p - 
   48*G3pp*G4p*G4ppx*G5p + 48*G3p*G4pp*G4ppx*G5p - 24*G2px*G3p*G4px*G5p + 24*G3p*G3pp*G4px*G5p - 24*G2pp*G3x*G4px*G5p - 48*G2ppx*G4p*G4px*G5p + 48*G3ppp*G4p*G4px*G5p + 72*G2px*G4pp*G4px*G5p - 
   96*G3pp*G4pp*G4px*G5p + 48*G3p*G4ppp*G4px*G5p + 24*G2pp*pow(G4px,2)*G5p - 16*pow(G2px,2)*G4x*G5p - 16*G2pp*G2xx*G4x*G5p - 32*G2ppx*G3p*G4x*G5p + 48*G2px*G3pp*G4x*G5p - 32*pow(G3pp,2)*G4x*G5p + 
   32*G3p*G3ppp*G4x*G5p + 16*G2pp*G3px*G4x*G5p + 4*pow(G2px,2)*pow(G5p,2) + 4*G2pp*G2xx*pow(G5p,2) + 8*G2ppx*G3p*pow(G5p,2) - 12*G2px*G3pp*pow(G5p,2) + 8*pow(G3pp,2)*pow(G5p,2) - 
   8*G3p*G3ppp*pow(G5p,2) - 4*G2pp*G3px*pow(G5p,2) - 24*pow(G3p,2)*G3x*G5pp - 24*G2xx*G3p*G4p*G5pp + 24*G3p*G3px*G4p*G5pp + 24*G2px*G3x*G4p*G5pp - 24*G3pp*G3x*G4p*G5pp + 24*G3p*G3x*G4pp*G5pp + 
   24*G2xx*G4p*G4pp*G5pp - 24*G3px*G4p*G4pp*G5pp + 48*pow(G3p,2)*G4px*G5pp - 48*G2px*G4p*G4px*G5pp + 48*G3pp*G4p*G4px*G5pp - 48*G3p*G4pp*G4px*G5pp + 108*G3px*pow(G3x,2)*G4*pow(H,2) + 
   144*G3px*G3pxx*pow(G4,2)*pow(H,2) - 144*G3ppx*G3xx*pow(G4,2)*pow(H,2) - 270*pow(G3x,3)*G4p*pow(H,2) + 288*G3pxx*G3x*G4*G4p*pow(H,2) - 72*G3px*G3xx*G4*G4p*pow(H,2) + 
   288*G3x*G3xx*pow(G4p,2)*pow(H,2) - 72*G3x*G3xx*G4*G4pp*pow(H,2) + 288*G3xx*pow(G4,2)*G4pppx*pow(H,2) - 288*pow(G3x,2)*G4*G4ppx*pow(H,2) - 288*G3pxx*pow(G4,2)*G4ppx*pow(H,2) + 
   432*G3xx*G4*G4p*G4ppx*pow(H,2) - 96*G2xx*pow(G4,2)*G4ppxx*pow(H,2) - 192*G3px*pow(G4,2)*G4ppxx*pow(H,2) - 288*G3x*G4*G4p*G4ppxx*pow(H,2) + 576*pow(G4,2)*G4ppx*G4ppxx*pow(H,2) + 
   72*G2xx*G3x*G4*G4px*pow(H,2) - 432*G3px*G3x*G4*G4px*pow(H,2) - 288*G3p*G3xx*G4*G4px*pow(H,2) + 1800*pow(G3x,2)*G4p*G4px*pow(H,2) - 864*G3pxx*G4*G4p*G4px*pow(H,2) - 
   720*G3xx*pow(G4p,2)*G4px*pow(H,2) + 720*G3xx*G4*G4pp*G4px*pow(H,2) + 1008*G3x*G4*G4ppx*G4px*pow(H,2) + 1152*G4*G4p*G4ppxx*G4px*pow(H,2) - 336*G2xx*G4*pow(G4px,2)*pow(H,2) + 
   912*G3px*G4*pow(G4px,2)*pow(H,2) - 3672*G3x*G4p*pow(G4px,2)*pow(H,2) - 1440*G4*G4ppx*pow(G4px,2)*pow(H,2) + 1728*G4p*pow(G4px,3)*pow(H,2) - 864*G3p*G3x*G4*G4pxx*pow(H,2) + 
   96*G2pxx*pow(G4,2)*G4pxx*pow(H,2) + 192*G3ppx*pow(G4,2)*G4pxx*pow(H,2) - 96*G2xx*G4*G4p*G4pxx*pow(H,2) - 624*G3px*G4*G4p*G4pxx*pow(H,2) - 1008*G3x*pow(G4p,2)*G4pxx*pow(H,2) - 
   144*G3x*G4*G4pp*G4pxx*pow(H,2) - 576*pow(G4,2)*G4pppx*G4pxx*pow(H,2) + 864*G4*G4p*G4ppx*G4pxx*pow(H,2) + 3264*G3p*G4*G4px*G4pxx*pow(H,2) + 3168*pow(G4p,2)*G4px*G4pxx*pow(H,2) - 
   1440*G4*G4pp*G4px*G4pxx*pow(H,2) + 96*G2px*pow(G4,2)*G4pxxx*pow(H,2) - 96*G3pp*pow(G4,2)*G4pxxx*pow(H,2) - 288*G3p*G4*G4p*G4pxxx*pow(H,2) + 288*G4*G4p*G4pp*G4pxxx*pow(H,2) + 
   360*G2xx*G3px*G4*G4x*pow(H,2) - 648*pow(G3px,2)*G4*G4x*pow(H,2) + 864*G3p*G3pxx*G4*G4x*pow(H,2) - 504*G2pxx*G3x*G4*G4x*pow(H,2) + 792*G3ppx*G3x*G4*G4x*pow(H,2) + 
   288*G2px*G3xx*G4*G4x*pow(H,2) - 720*G3pp*G3xx*G4*G4x*pow(H,2) - 288*G2xx*G3x*G4p*G4x*pow(H,2) + 648*G3px*G3x*G4p*G4x*pow(H,2) + 288*G3p*G3xx*G4p*G4x*pow(H,2) + 
   144*G3pxx*pow(G4p,2)*G4x*pow(H,2) + 1404*pow(G3x,2)*G4pp*G4x*pow(H,2) - 288*G3pxx*G4*G4pp*G4x*pow(H,2) - 144*G3xx*G4p*G4pp*G4x*pow(H,2) + 288*G3xx*G4*G4ppp*G4x*pow(H,2) - 
   576*G3x*G4*G4pppx*G4x*pow(H,2) - 912*G2xx*G4*G4ppx*G4x*pow(H,2) + 2640*G3px*G4*G4ppx*G4x*pow(H,2) + 1008*G3x*G4p*G4ppx*G4x*pow(H,2) - 2304*G4*pow(G4ppx,2)*G4x*pow(H,2) - 
   2112*G3p*G4*G4ppxx*G4x*pow(H,2) + 576*G4*G4pp*G4ppxx*G4x*pow(H,2) + 1200*G2pxx*G4*G4px*G4x*pow(H,2) - 2352*G3ppx*G4*G4px*G4x*pow(H,2) + 480*G2xx*G4p*G4px*G4x*pow(H,2) - 
   1344*G3px*G4p*G4px*G4x*pow(H,2) - 7200*G3x*G4pp*G4px*G4x*pow(H,2) + 2304*G4*G4pppx*G4px*G4x*pow(H,2) - 1728*G4p*G4ppx*G4px*G4x*pow(H,2) + 9360*G4pp*pow(G4px,2)*G4x*pow(H,2) - 
   1152*G2px*G4*G4pxx*G4x*pow(H,2) + 2208*G3pp*G4*G4pxx*G4x*pow(H,2) - 1920*G3p*G4p*G4pxx*G4x*pow(H,2) - 288*G4p*G4pp*G4pxx*G4x*pow(H,2) - 576*G4*G4ppp*G4pxx*G4x*pow(H,2) + 
   288*G3p*G3px*pow(G4x,2)*pow(H,2) + 360*G2px*G3x*pow(G4x,2)*pow(H,2) - 792*G3pp*G3x*pow(G4x,2)*pow(H,2) - 432*G2pxx*G4p*pow(G4x,2)*pow(H,2) - 144*G3ppx*G4p*pow(G4x,2)*pow(H,2) + 
   144*G2xx*G4pp*pow(G4x,2)*pow(H,2) + 576*G3px*G4pp*pow(G4x,2)*pow(H,2) + 144*G3x*G4ppp*pow(G4x,2)*pow(H,2) + 1152*G4p*G4pppx*pow(G4x,2)*pow(H,2) + 
   192*G3p*G4ppx*pow(G4x,2)*pow(H,2) - 2016*G4pp*G4ppx*pow(G4x,2)*pow(H,2) - 720*G2px*G4px*pow(G4x,2)*pow(H,2) + 1200*G3pp*G4px*pow(G4x,2)*pow(H,2) + 
   288*G4ppp*G4px*pow(G4x,2)*pow(H,2) + 96*G2ppx*pow(G4x,3)*pow(H,2) - 96*G3ppp*pow(G4x,3)*pow(H,2) - 720*G3p*G3px*G4*G4xx*pow(H,2) - 504*G2px*G3x*G4*G4xx*pow(H,2) + 
   432*G3pp*G3x*G4*G4xx*pow(H,2) + 2160*G3p*G3x*G4p*G4xx*pow(H,2) + 240*G2pxx*G4*G4p*G4xx*pow(H,2) + 336*G3ppx*G4*G4p*G4xx*pow(H,2) + 240*G2xx*pow(G4p,2)*G4xx*pow(H,2) - 
   168*G3px*pow(G4p,2)*G4xx*pow(H,2) + 48*G2xx*G4*G4pp*G4xx*pow(H,2) - 624*G3px*G4*G4pp*G4xx*pow(H,2) - 4464*G3x*G4p*G4pp*G4xx*pow(H,2) + 1152*G3x*G4*G4ppp*G4xx*pow(H,2) - 
   1152*G4*G4p*G4pppx*G4xx*pow(H,2) + 1248*G3p*G4*G4ppx*G4xx*pow(H,2) - 1872*pow(G4p,2)*G4ppx*G4xx*pow(H,2) + 1728*G4*G4pp*G4ppx*G4xx*pow(H,2) + 1584*G2px*G4*G4px*G4xx*pow(H,2) - 
   1344*G3pp*G4*G4px*G4xx*pow(H,2) - 3648*G3p*G4p*G4px*G4xx*pow(H,2) + 10656*G4p*G4pp*G4px*G4xx*pow(H,2) - 2880*G4*G4ppp*G4px*G4xx*pow(H,2) + 384*G2ppx*G4*G4x*G4xx*pow(H,2) - 
   384*G3ppp*G4*G4x*G4xx*pow(H,2) - 1200*G2px*G4p*G4x*G4xx*pow(H,2) + 3168*G3pp*G4p*G4x*G4xx*pow(H,2) - 3360*G3p*G4pp*G4x*G4xx*pow(H,2) + 2304*pow(G4pp,2)*G4x*G4xx*pow(H,2) - 
   2304*G4p*G4ppp*G4x*G4xx*pow(H,2) + 192*G2pp*G4*pow(G4xx,2)*pow(H,2) - 96*G2ppx*pow(G4,2)*G4xxx*pow(H,2) + 96*G3ppp*pow(G4,2)*G4xxx*pow(H,2) + 144*G3pp*G4*G4p*G4xxx*pow(H,2) - 
   288*G3p*pow(G4p,2)*G4xxx*pow(H,2) + 288*G3p*G4*G4pp*G4xxx*pow(H,2) + 288*pow(G4p,2)*G4pp*G4xxx*pow(H,2) - 288*G4*pow(G4pp,2)*G4xxx*pow(H,2) - 288*G4*G4p*G4ppp*G4xxx*pow(H,2) - 
   192*G2pp*G4*G4x*G4xxx*pow(H,2) - 216*G2xx*G3px*G4*G5p*pow(H,2) + 216*pow(G3px,2)*G4*G5p*pow(H,2) - 576*G3p*G3pxx*G4*G5p*pow(H,2) + 360*G2pxx*G3x*G4*G5p*pow(H,2) - 
   360*G3ppx*G3x*G4*G5p*pow(H,2) - 144*G2px*G3xx*G4*G5p*pow(H,2) + 432*G3pp*G3xx*G4*G5p*pow(H,2) + 252*G2xx*G3x*G4p*G5p*pow(H,2) - 684*G3px*G3x*G4p*G5p*pow(H,2) - 
   288*G3p*G3xx*G4p*G5p*pow(H,2) - 72*G3pxx*pow(G4p,2)*G5p*pow(H,2) - 1188*pow(G3x,2)*G4pp*G5p*pow(H,2) + 288*G3pxx*G4*G4pp*G5p*pow(H,2) + 216*G3xx*G4p*G4pp*G5p*pow(H,2) - 
   288*G3xx*G4*G4ppp*G5p*pow(H,2) + 624*G2xx*G4*G4ppx*G5p*pow(H,2) - 912*G3px*G4*G4ppx*G5p*pow(H,2) - 864*G3x*G4p*G4ppx*G5p*pow(H,2) + 576*G4*pow(G4ppx,2)*G5p*pow(H,2) + 
   1344*G3p*G4*G4ppxx*G5p*pow(H,2) - 576*G4*G4pp*G4ppxx*G5p*pow(H,2) + 576*G3p*G3x*G4px*G5p*pow(H,2) - 912*G2pxx*G4*G4px*G5p*pow(H,2) + 1200*G3ppx*G4*G4px*G5p*pow(H,2) - 
   432*G2xx*G4p*G4px*G5p*pow(H,2) + 1944*G3px*G4p*G4px*G5p*pow(H,2) + 5760*G3x*G4pp*G4px*G5p*pow(H,2) - 576*G4*G4pppx*G4px*G5p*pow(H,2) + 432*G4p*G4ppx*G4px*G5p*pow(H,2) - 
   1632*G3p*pow(G4px,2)*G5p*pow(H,2) - 6768*G4pp*pow(G4px,2)*G5p*pow(H,2) + 384*G2px*G4*G4pxx*G5p*pow(H,2) - 1056*G3pp*G4*G4pxx*G5p*pow(H,2) + 1824*G3p*G4p*G4pxx*G5p*pow(H,2) - 
   720*G4p*G4pp*G4pxx*G5p*pow(H,2) + 576*G4*G4ppp*G4pxx*G5p*pow(H,2) - 720*G3p*G3px*G4x*G5p*pow(H,2) - 612*G2px*G3x*G4x*G5p*pow(H,2) + 1548*G3pp*G3x*G4x*G5p*pow(H,2) + 
   552*G2pxx*G4p*G4x*G5p*pow(H,2) + 312*G3ppx*G4p*G4x*G5p*pow(H,2) - 120*G2xx*G4pp*G4x*G5p*pow(H,2) - 816*G3px*G4pp*G4x*G5p*pow(H,2) - 648*G3x*G4ppp*G4x*G5p*pow(H,2) - 
   1728*G4p*G4pppx*G4x*G5p*pow(H,2) + 576*G3p*G4ppx*G4x*G5p*pow(H,2) + 2736*G4pp*G4ppx*G4x*G5p*pow(H,2) + 1320*G2px*G4px*G4x*G5p*pow(H,2) - 2760*G3pp*G4px*G4x*G5p*pow(H,2) + 
   432*G4ppp*G4px*G4x*G5p*pow(H,2) - 288*G2ppx*pow(G4x,2)*G5p*pow(H,2) + 288*G3ppp*pow(G4x,2)*G5p*pow(H,2) + 792*G2px*G4p*G4xx*G5p*pow(H,2) - 2640*G3pp*G4p*G4xx*G5p*pow(H,2) + 
   2256*G3p*G4pp*G4xx*G5p*pow(H,2) - 1728*pow(G4pp,2)*G4xx*G5p*pow(H,2) + 2880*G4p*G4ppp*G4xx*G5p*pow(H,2) + 192*G2pp*G4x*G4xx*G5p*pow(H,2) + 96*G2pp*G4*G4xxx*G5p*pow(H,2) + 
   360*G3p*G3px*pow(G5p,2)*pow(H,2) + 216*G2px*G3x*pow(G5p,2)*pow(H,2) - 684*G3pp*G3x*pow(G5p,2)*pow(H,2) - 168*G2pxx*G4p*pow(G5p,2)*pow(H,2) - 120*G3ppx*G4p*pow(G5p,2)*pow(H,2) + 
   24*G2xx*G4pp*pow(G5p,2)*pow(H,2) + 192*G3px*G4pp*pow(G5p,2)*pow(H,2) + 504*G3x*G4ppp*pow(G5p,2)*pow(H,2) + 576*G4p*G4pppx*pow(G5p,2)*pow(H,2) - 480*G3p*G4ppx*pow(G5p,2)*pow(H,2) - 
   720*G4pp*G4ppx*pow(G5p,2)*pow(H,2) - 432*G2px*G4px*pow(G5p,2)*pow(H,2) + 1248*G3pp*G4px*pow(G5p,2)*pow(H,2) - 720*G4ppp*G4px*pow(G5p,2)*pow(H,2) + 
   264*G2ppx*G4x*pow(G5p,2)*pow(H,2) - 264*G3ppp*G4x*pow(G5p,2)*pow(H,2) - 96*G2pp*G4xx*pow(G5p,2)*pow(H,2) - 72*G2ppx*pow(G5p,3)*pow(H,2) + 72*G3ppp*pow(G5p,3)*pow(H,2) - 
   108*G2xx*G3x*G4*G5pp*pow(H,2) - 252*G3px*G3x*G4*G5pp*pow(H,2) + 288*G3p*G3xx*G4*G5pp*pow(H,2) - 198*pow(G3x,2)*G4p*G5pp*pow(H,2) + 144*G3pxx*G4*G4p*G5pp*pow(H,2) + 
   72*G3xx*pow(G4p,2)*G5pp*pow(H,2) - 432*G3xx*G4*G4pp*G5pp*pow(H,2) + 864*G3x*G4*G4ppx*G5pp*pow(H,2) - 288*G4*G4p*G4ppxx*G5pp*pow(H,2) + 504*G2xx*G4*G4px*G5pp*pow(H,2) + 
   72*G3px*G4*G4px*G5pp*pow(H,2) + 1080*G3x*G4p*G4px*G5pp*pow(H,2) - 1440*G4*G4ppx*G4px*G5pp*pow(H,2) - 792*G4p*pow(G4px,2)*G5pp*pow(H,2) - 768*G3p*G4*G4pxx*G5pp*pow(H,2) - 
   576*pow(G4p,2)*G4pxx*G5pp*pow(H,2) + 864*G4*G4pp*G4pxx*G5pp*pow(H,2) - 936*G3p*G3x*G4x*G5pp*pow(H,2) - 432*G3px*G4p*G4x*G5pp*pow(H,2) + 1080*G3x*G4pp*G4x*G5pp*pow(H,2) + 
   1152*G4p*G4ppx*G4x*G5pp*pow(H,2) + 2352*G3p*G4px*G4x*G5pp*pow(H,2) - 2736*G4pp*G4px*G4x*G5pp*pow(H,2) + 240*G2px*pow(G4x,2)*G5pp*pow(H,2) - 288*G3pp*pow(G4x,2)*G5pp*pow(H,2) - 
   480*G2px*G4*G4xx*G5pp*pow(H,2) + 576*G3pp*G4*G4xx*G5pp*pow(H,2) - 48*G3p*G4p*G4xx*G5pp*pow(H,2) - 1152*G4p*G4pp*G4xx*G5pp*pow(H,2) + 324*G3p*G3x*G5p*G5pp*pow(H,2) + 
   48*G2pxx*G4*G5p*G5pp*pow(H,2) - 48*G3ppx*G4*G5p*G5pp*pow(H,2) - 12*G2xx*G4p*G5p*G5pp*pow(H,2) - 60*G3px*G4p*G5p*G5pp*pow(H,2) - 792*G3x*G4pp*G5p*G5pp*pow(H,2) - 
   552*G3p*G4px*G5p*G5pp*pow(H,2) + 1728*G4pp*G4px*G5p*G5pp*pow(H,2) - 384*G2px*G4x*G5p*G5pp*pow(H,2) + 480*G3pp*G4x*G5p*G5pp*pow(H,2) + 108*G2px*pow(G5p,2)*G5pp*pow(H,2) - 
   144*G3pp*pow(G5p,2)*G5pp*pow(H,2) - 96*G2xx*G4*pow(G5pp,2)*pow(H,2) + 96*G3px*G4*pow(G5pp,2)*pow(H,2) - 216*G3x*G4p*pow(G5pp,2)*pow(H,2) + 288*G4p*G4px*pow(G5pp,2)*pow(H,2) - 
   96*G3p*G4x*pow(G5pp,2)*pow(H,2) + 180*pow(G3x,2)*G4*G5ppp*pow(H,2) - 144*G3xx*G4*G4p*G5ppp*pow(H,2) - 864*G3x*G4*G4px*G5ppp*pow(H,2) + 1008*G4*pow(G4px,2)*G5ppp*pow(H,2) + 
   288*G4*G4p*G4pxx*G5ppp*pow(H,2) - 432*G3x*G4p*G4x*G5ppp*pow(H,2) + 576*G4p*G4px*G4x*G5ppp*pow(H,2) + 96*G3p*pow(G4x,2)*G5ppp*pow(H,2) - 192*G3p*G4*G4xx*G5ppp*pow(H,2) + 
   864*pow(G4p,2)*G4xx*G5ppp*pow(H,2) - 48*G2xx*G4*G5p*G5ppp*pow(H,2) + 48*G3px*G4*G5p*G5ppp*pow(H,2) + 648*G3x*G4p*G5p*G5ppp*pow(H,2) - 1152*G4p*G4px*G5p*G5ppp*pow(H,2) - 
   192*G3p*G4x*G5p*G5ppp*pow(H,2) + 72*G3p*pow(G5p,2)*G5ppp*pow(H,2) + 48*G2xx*pow(G4,2)*G5pppx*pow(H,2) - 48*G3px*pow(G4,2)*G5pppx*pow(H,2) - 144*G3x*G4*G4p*G5pppx*pow(H,2) + 
   288*G4*G4p*G4px*G5pppx*pow(H,2) + 192*G3p*G4*G4x*G5pppx*pow(H,2) - 144*pow(G4p,2)*G4x*G5pppx*pow(H,2) - 96*G3p*G4*G5p*G5pppx*pow(H,2) + 72*pow(G4p,2)*G5p*G5pppx*pow(H,2) + 
   384*G3p*G3x*G4*G5ppx*pow(H,2) - 48*G2pxx*pow(G4,2)*G5ppx*pow(H,2) + 48*G3ppx*pow(G4,2)*G5ppx*pow(H,2) + 72*G2xx*G4*G4p*G5ppx*pow(H,2) + 360*G3px*G4*G4p*G5ppx*pow(H,2) + 
   108*G3x*pow(G4p,2)*G5ppx*pow(H,2) + 216*G3x*G4*G4pp*G5ppx*pow(H,2) - 864*G4*G4p*G4ppx*G5ppx*pow(H,2) - 1248*G3p*G4*G4px*G5ppx*pow(H,2) - 648*pow(G4p,2)*G4px*G5ppx*pow(H,2) - 
   144*G4*G4pp*G4px*G5ppx*pow(H,2) + 384*G2px*G4*G4x*G5ppx*pow(H,2) - 480*G3pp*G4*G4x*G5ppx*pow(H,2) + 576*G3p*G4p*G4x*G5ppx*pow(H,2) + 432*G4p*G4pp*G4x*G5ppx*pow(H,2) - 
   96*G2px*G4*G5p*G5ppx*pow(H,2) + 144*G3pp*G4*G5p*G5ppx*pow(H,2) - 576*G3p*G4p*G5p*G5ppx*pow(H,2) + 72*G4p*G4pp*G5p*G5ppx*pow(H,2) + 96*G3p*G4*G5pp*G5ppx*pow(H,2) + 
   216*pow(G4p,2)*G5pp*G5ppx*pow(H,2) - 48*G2px*pow(G4,2)*G5ppxx*pow(H,2) + 48*G3pp*pow(G4,2)*G5ppxx*pow(H,2) + 144*G3p*G4*G4p*G5ppxx*pow(H,2) - 144*G4*G4p*G4pp*G5ppxx*pow(H,2) - 
   72*G2xx*G3p*G4*G5px*pow(H,2) + 528*G3p*G3px*G4*G5px*pow(H,2) + 336*G2px*G3x*G4*G5px*pow(H,2) - 276*G3pp*G3x*G4*G5px*pow(H,2) - 1008*G3p*G3x*G4p*G5px*pow(H,2) - 
   144*G2pxx*G4*G4p*G5px*pow(H,2) - 144*G3ppx*G4*G4p*G5px*pow(H,2) - 120*G2xx*pow(G4p,2)*G5px*pow(H,2) + 192*G3px*pow(G4p,2)*G5px*pow(H,2) + 72*G2xx*G4*G4pp*G5px*pow(H,2) + 
   144*G3px*G4*G4pp*G5px*pow(H,2) + 2448*G3x*G4p*G4pp*G5px*pow(H,2) - 648*G3x*G4*G4ppp*G5px*pow(H,2) + 576*G4*G4p*G4pppx*G5px*pow(H,2) - 816*G3p*G4*G4ppx*G5px*pow(H,2) + 
   720*pow(G4p,2)*G4ppx*G5px*pow(H,2) - 720*G4*G4pp*G4ppx*G5px*pow(H,2) - 960*G2px*G4*G4px*G5px*pow(H,2) + 792*G3pp*G4*G4px*G5px*pow(H,2) + 1248*G3p*G4p*G4px*G5px*pow(H,2) - 
   5472*G4p*G4pp*G4px*G5px*pow(H,2) + 1584*G4*G4ppp*G4px*G5px*pow(H,2) - 48*pow(G3p,2)*G4x*G5px*pow(H,2) - 288*G2ppx*G4*G4x*G5px*pow(H,2) + 288*G3ppp*G4*G4x*G5px*pow(H,2) + 
   672*G2px*G4p*G4x*G5px*pow(H,2) - 1608*G3pp*G4p*G4x*G5px*pow(H,2) + 1920*G3p*G4pp*G4x*G5px*pow(H,2) - 1296*pow(G4pp,2)*G4x*G5px*pow(H,2) + 1008*G4p*G4ppp*G4x*G5px*pow(H,2) - 
   48*G2pp*pow(G4x,2)*G5px*pow(H,2) - 192*G2pp*G4*G4xx*G5px*pow(H,2) - 24*pow(G3p,2)*G5p*G5px*pow(H,2) + 48*G2ppx*G4*G5p*G5px*pow(H,2) - 48*G3ppp*G4*G5p*G5px*pow(H,2) - 
   384*G2px*G4p*G5p*G5px*pow(H,2) + 1284*G3pp*G4p*G5p*G5px*pow(H,2) - 1200*G3p*G4pp*G5p*G5px*pow(H,2) + 936*pow(G4pp,2)*G5p*G5px*pow(H,2) - 1368*G4p*G4ppp*G5p*G5px*pow(H,2) - 
   48*G2pp*G4x*G5p*G5px*pow(H,2) + 36*G2pp*pow(G5p,2)*G5px*pow(H,2) + 240*G2px*G4*G5pp*G5px*pow(H,2) - 288*G3pp*G4*G5pp*G5px*pow(H,2) + 168*G3p*G4p*G5pp*G5px*pow(H,2) + 
   504*G4p*G4pp*G5pp*G5px*pow(H,2) + 96*G3p*G4*G5ppp*G5px*pow(H,2) - 432*pow(G4p,2)*G5ppp*G5px*pow(H,2) + 48*G2pp*G4*pow(G5px,2)*pow(H,2) + 48*pow(G3p,2)*G4*G5pxx*pow(H,2) + 
   48*G2ppx*pow(G4,2)*G5pxx*pow(H,2) - 48*G3ppp*pow(G4,2)*G5pxx*pow(H,2) - 48*G2px*G4*G4p*G5pxx*pow(H,2) - 24*G3pp*G4*G4p*G5pxx*pow(H,2) + 216*G3p*pow(G4p,2)*G5pxx*pow(H,2) - 
   192*G3p*G4*G4pp*G5pxx*pow(H,2) - 216*pow(G4p,2)*G4pp*G5pxx*pow(H,2) + 144*G4*pow(G4pp,2)*G5pxx*pow(H,2) + 144*G4*G4p*G4ppp*G5pxx*pow(H,2) + 96*G2pp*G4*G4x*G5pxx*pow(H,2) - 
   48*G2pp*G4*G5p*G5pxx*pow(H,2) + 48*G2px*G2xx*G4*G5x*pow(H,2) + 168*G2pxx*G3p*G4*G5x*pow(H,2) - 132*G2xx*G3pp*G4*G5x*pow(H,2) - 264*G3p*G3ppx*G4*G5x*pow(H,2) - 
   120*G2px*G3px*G4*G5x*pow(H,2) + 252*G3pp*G3px*G4*G5x*pow(H,2) + 24*G2ppx*G3x*G4*G5x*pow(H,2) - 24*G3ppp*G3x*G4*G5x*pow(H,2) - 48*G2pp*G3xx*G4*G5x*pow(H,2) + 72*G2xx*G3p*G4p*G5x*pow(H,2) - 
   312*G3p*G3px*G4p*G5x*pow(H,2) - 192*G2px*G3x*G4p*G5x*pow(H,2) + 564*G3pp*G3x*G4p*G5x*pow(H,2) + 24*G2pxx*pow(G4p,2)*G5x*pow(H,2) + 48*G3ppx*pow(G4p,2)*G5x*pow(H,2) - 
   528*G3p*G3x*G4pp*G5x*pow(H,2) - 72*G2pxx*G4*G4pp*G5x*pow(H,2) + 72*G3ppx*G4*G4pp*G5x*pow(H,2) - 48*G2xx*G4p*G4pp*G5x*pow(H,2) - 24*G3px*G4p*G4pp*G5x*pow(H,2) + 
   360*G3x*pow(G4pp,2)*G5x*pow(H,2) + 72*G2xx*G4*G4ppp*G5x*pow(H,2) - 72*G3px*G4*G4ppp*G5x*pow(H,2) - 504*G3x*G4p*G4ppp*G5x*pow(H,2) + 192*G3p*G4*G4pppx*G5x*pow(H,2) - 
   144*pow(G4p,2)*G4pppx*G5x*pow(H,2) + 240*G2px*G4*G4ppx*G5x*pow(H,2) - 336*G3pp*G4*G4ppx*G5x*pow(H,2) + 432*G3p*G4p*G4ppx*G5x*pow(H,2) + 288*G4p*G4pp*G4ppx*G5x*pow(H,2) + 
   48*pow(G3p,2)*G4px*G5x*pow(H,2) - 144*G2ppx*G4*G4px*G5x*pow(H,2) + 144*G3ppp*G4*G4px*G5x*pow(H,2) + 384*G2px*G4p*G4px*G5x*pow(H,2) - 1104*G3pp*G4p*G4px*G5x*pow(H,2) + 
   1248*G3p*G4pp*G4px*G5x*pow(H,2) - 864*pow(G4pp,2)*G4px*G5x*pow(H,2) + 864*G4p*G4ppp*G4px*G5x*pow(H,2) + 96*G2pp*G4*G4pxx*G5x*pow(H,2) - 120*G2px*G3p*G4x*G5x*pow(H,2) + 
   120*G3p*G3pp*G4x*G5x*pow(H,2) - 24*G2pp*G3x*G4x*G5x*pow(H,2) - 240*G2ppx*G4p*G4x*G5x*pow(H,2) + 240*G3ppp*G4p*G4x*G5x*pow(H,2) + 360*G2px*G4pp*G4x*G5x*pow(H,2) - 
   480*G3pp*G4pp*G4x*G5x*pow(H,2) + 240*G3p*G4ppp*G4x*G5x*pow(H,2) - 48*G2pp*G4px*G4x*G5x*pow(H,2) + 96*G2pp*G4p*G4xx*G5x*pow(H,2) + 84*G2px*G3p*G5p*G5x*pow(H,2) - 
   84*G3p*G3pp*G5p*G5x*pow(H,2) + 36*G2pp*G3x*G5p*G5x*pow(H,2) + 168*G2ppx*G4p*G5p*G5x*pow(H,2) - 168*G3ppp*G4p*G5p*G5x*pow(H,2) - 252*G2px*G4pp*G5p*G5x*pow(H,2) + 
   336*G3pp*G4pp*G5p*G5x*pow(H,2) - 168*G3p*G4ppp*G5p*G5x*pow(H,2) - 24*G2pp*G4px*G5p*G5x*pow(H,2) + 24*pow(G3p,2)*G5pp*G5x*pow(H,2) - 48*G2px*G4p*G5pp*G5x*pow(H,2) + 
   72*G3pp*G4p*G5pp*G5x*pow(H,2) - 72*G3p*G4pp*G5pp*G5x*pow(H,2) - 48*G3p*G4p*G5ppp*G5x*pow(H,2) - 48*G2pp*G4p*G5px*G5x*pow(H,2) - 8*G2pp*G3p*pow(G5x,2)*pow(H,2) + 
   24*G2px*G3p*G4*G5xx*pow(H,2) - 24*G3p*G3pp*G4*G5xx*pow(H,2) + 24*G2pp*G3x*G4*G5xx*pow(H,2) - 48*pow(G3p,2)*G4p*G5xx*pow(H,2) + 48*G2ppx*G4*G4p*G5xx*pow(H,2) - 
   48*G3ppp*G4*G4p*G5xx*pow(H,2) + 48*G2px*pow(G4p,2)*G5xx*pow(H,2) - 156*G3pp*pow(G4p,2)*G5xx*pow(H,2) - 72*G2px*G4*G4pp*G5xx*pow(H,2) + 96*G3pp*G4*G4pp*G5xx*pow(H,2) + 
   192*G3p*G4p*G4pp*G5xx*pow(H,2) - 144*G4p*pow(G4pp,2)*G5xx*pow(H,2) - 48*G3p*G4*G4ppp*G5xx*pow(H,2) + 216*pow(G4p,2)*G4ppp*G5xx*pow(H,2) - 48*G2pp*G4*G4px*G5xx*pow(H,2) + 
   48*G2pp*G4p*G4x*G5xx*pow(H,2) - 24*G2pp*G4p*G5p*G5xx*pow(H,2) - 864*G3xx*pow(G4,2)*G4pxx*pow(H,4) + 4032*pow(G4,2)*pow(G4pxx,2)*pow(H,4) + 864*G3x*pow(G4,2)*G4pxxx*pow(H,4) - 
   1728*pow(G4,2)*G4px*G4pxxx*pow(H,4) + 1728*G3xx*G4*G4px*G4x*pow(H,4) - 2592*G3x*G4*G4pxx*G4x*pow(H,4) + 1152*G4*G4px*G4pxx*G4x*pow(H,4) + 3744*G4*G4p*G4pxxx*G4x*pow(H,4) - 
   2592*G3pxx*G4*pow(G4x,2)*pow(H,4) - 1728*G3xx*G4p*pow(G4x,2)*pow(H,4) + 5184*G4*G4ppxx*pow(G4x,2)*pow(H,4) - 3024*G3x*G4px*pow(G4x,2)*pow(H,4) + 
   7200*pow(G4px,2)*pow(G4x,2)*pow(H,4) + 1728*G4p*G4pxx*pow(G4x,2)*pow(H,4) + 1296*G3px*pow(G4x,3)*pow(H,4) - 7200*G4ppx*pow(G4x,3)*pow(H,4) + 864*G3pxx*pow(G4,2)*G4xx*pow(H,4) + 
   864*G3xx*G4*G4p*G4xx*pow(H,4) - 4032*pow(G4,2)*G4ppxx*G4xx*pow(H,4) + 3888*G3x*G4*G4px*G4xx*pow(H,4) - 11232*G4*pow(G4px,2)*G4xx*pow(H,4) - 2880*G4*G4p*G4pxx*G4xx*pow(H,4) + 
   2160*G3px*G4*G4x*G4xx*pow(H,4) - 15120*G3x*G4p*G4x*G4xx*pow(H,4) - 3168*G4*G4ppx*G4x*G4xx*pow(H,4) + 44064*G4p*G4px*G4x*G4xx*pow(H,4) + 17568*G4pp*pow(G4x,2)*G4xx*pow(H,4) + 
   6336*pow(G4p,2)*pow(G4xx,2)*pow(H,4) - 864*G3px*pow(G4,2)*G4xxx*pow(H,4) + 432*G3x*G4*G4p*G4xxx*pow(H,4) + 1728*pow(G4,2)*G4ppx*G4xxx*pow(H,4) - 576*G4*G4p*G4px*G4xxx*pow(H,4) + 
   2304*pow(G4p,2)*G4x*G4xxx*pow(H,4) - 2016*G4*G4pp*G4x*G4xxx*pow(H,4) - 864*G3xx*G4*G4px*G5p*pow(H,4) + 864*G3x*G4*G4pxx*G5p*pow(H,4) - 1728*G4*G4px*G4pxx*G5p*pow(H,4) - 
   2592*G4*G4p*G4pxxx*G5p*pow(H,4) + 4320*G3pxx*G4*G4x*G5p*pow(H,4) + 3024*G3xx*G4p*G4x*G5p*pow(H,4) - 8064*G4*G4ppxx*G4x*G5p*pow(H,4) + 3888*G3x*G4px*G4x*G5p*pow(H,4) - 
   8352*pow(G4px,2)*G4x*G5p*pow(H,4) - 7488*G4p*G4pxx*G4x*G5p*pow(H,4) - 2376*G3px*pow(G4x,2)*G5p*pow(H,4) + 16560*G4ppx*pow(G4x,2)*G5p*pow(H,4) - 432*G3px*G4*G4xx*G5p*pow(H,4) + 
   13392*G3x*G4p*G4xx*G5p*pow(H,4) + 1440*G4*G4ppx*G4xx*G5p*pow(H,4) - 38736*G4p*G4px*G4xx*G5p*pow(H,4) - 30096*G4pp*G4x*G4xx*G5p*pow(H,4) - 1872*pow(G4p,2)*G4xxx*G5p*pow(H,4) + 
   864*G4*G4pp*G4xxx*G5p*pow(H,4) - 1728*G3pxx*G4*pow(G5p,2)*pow(H,4) - 1296*G3xx*G4p*pow(G5p,2)*pow(H,4) + 2880*G4*G4ppxx*pow(G5p,2)*pow(H,4) - 1080*G3x*G4px*pow(G5p,2)*pow(H,4) + 
   2016*pow(G4px,2)*pow(G5p,2)*pow(H,4) + 4608*G4p*G4pxx*pow(G5p,2)*pow(H,4) + 1080*G3px*G4x*pow(G5p,2)*pow(H,4) - 12384*G4ppx*G4x*pow(G5p,2)*pow(H,4) + 
   13680*G4pp*G4xx*pow(G5p,2)*pow(H,4) + 3024*G4ppx*pow(G5p,3)*pow(H,4) - 288*pow(G4,2)*G4pxxx*G5pp*pow(H,4) - 2160*G3xx*G4*G4x*G5pp*pow(H,4) + 2592*G4*G4pxx*G4x*G5pp*pow(H,4) + 
   6048*G3x*pow(G4x,2)*G5pp*pow(H,4) - 15840*G4px*pow(G4x,2)*G5pp*pow(H,4) - 5616*G3x*G4*G4xx*G5pp*pow(H,4) + 16128*G4*G4px*G4xx*G5pp*pow(H,4) - 10512*G4p*G4x*G4xx*G5pp*pow(H,4) + 
   144*G4*G4p*G4xxx*G5pp*pow(H,4) + 1296*G3xx*G4*G5p*G5pp*pow(H,4) + 864*G4*G4pxx*G5p*G5pp*pow(H,4) - 9396*G3x*G4x*G5p*G5pp*pow(H,4) + 24408*G4px*G4x*G5p*G5pp*pow(H,4) + 
   8640*G4p*G4xx*G5p*G5pp*pow(H,4) + 3564*G3x*pow(G5p,2)*G5pp*pow(H,4) - 9504*G4px*pow(G5p,2)*G5pp*pow(H,4) + 288*pow(G4x,2)*pow(G5pp,2)*pow(H,4) - 576*G4*G4xx*pow(G5pp,2)*pow(H,4) - 
   864*G4x*G5p*pow(G5pp,2)*pow(H,4) + 648*pow(G5p,2)*pow(G5pp,2)*pow(H,4) + 288*pow(G4x,3)*G5ppp*pow(H,4) + 1728*G4*G4x*G4xx*G5ppp*pow(H,4) + 288*pow(G4,2)*G4xxx*G5ppp*pow(H,4) - 
   1152*pow(G4x,2)*G5p*G5ppp*pow(H,4) - 2880*G4*G4xx*G5p*G5ppp*pow(H,4) + 1512*G4x*pow(G5p,2)*G5ppp*pow(H,4) - 648*pow(G5p,3)*G5ppp*pow(H,4) + 1152*pow(G4,2)*G4xx*G5pppx*pow(H,4) - 
   288*G4*G4x*G5p*G5pppx*pow(H,4) + 288*G4*pow(G5p,2)*G5pppx*pow(H,4) + 336*G3xx*pow(G4,2)*G5ppx*pow(H,4) - 3264*pow(G4,2)*G4pxx*G5ppx*pow(H,4) + 2640*G3x*G4*G4x*G5ppx*pow(H,4) - 
   5664*G4*G4px*G4x*G5ppx*pow(H,4) + 2496*G4p*pow(G4x,2)*G5ppx*pow(H,4) + 816*G4*G4p*G4xx*G5ppx*pow(H,4) - 1152*G3x*G4*G5p*G5ppx*pow(H,4) + 3648*G4*G4px*G5p*G5ppx*pow(H,4) - 
   1872*G4p*G4x*G5p*G5ppx*pow(H,4) + 96*G4p*pow(G5p,2)*G5ppx*pow(H,4) + 864*G4*G4x*G5pp*G5ppx*pow(H,4) - 1872*G4*G5p*G5pp*G5ppx*pow(H,4) + 720*pow(G4,2)*pow(G5ppx,2)*pow(H,4) - 
   528*G3x*pow(G4,2)*G5ppxx*pow(H,4) + 1152*pow(G4,2)*G4px*G5ppxx*pow(H,4) - 1680*G4*G4p*G4x*G5ppxx*pow(H,4) + 1104*G4*G4p*G5p*G5ppxx*pow(H,4) + 144*pow(G4,2)*G5pp*G5ppxx*pow(H,4) + 
   828*pow(G3x,2)*G4*G5px*pow(H,4) - 336*G3pxx*pow(G4,2)*G5px*pow(H,4) - 504*G3xx*G4*G4p*G5px*pow(H,4) + 2112*pow(G4,2)*G4ppxx*G5px*pow(H,4) - 6576*G3x*G4*G4px*G5px*pow(H,4) + 
   13296*G4*pow(G4px,2)*G5px*pow(H,4) + 48*G4*G4p*G4pxx*G5px*pow(H,4) + 408*G2xx*G4*G4x*G5px*pow(H,4) - 2256*G3px*G4*G4x*G5px*pow(H,4) + 10248*G3x*G4p*G4x*G5px*pow(H,4) + 
   3792*G4*G4ppx*G4x*G5px*pow(H,4) - 27984*G4p*G4px*G4x*G5px*pow(H,4) + 1248*G3p*pow(G4x,2)*G5px*pow(H,4) - 9696*G4pp*pow(G4x,2)*G5px*pow(H,4) - 3120*G3p*G4*G4xx*G5px*pow(H,4) - 
   7464*pow(G4p,2)*G4xx*G5px*pow(H,4) + 1200*G4*G4pp*G4xx*G5px*pow(H,4) - 264*G2xx*G4*G5p*G5px*pow(H,4) + 336*G3px*G4*G5p*G5px*pow(H,4) - 8628*G3x*G4p*G5p*G5px*pow(H,4) - 
   912*G4*G4ppx*G5p*G5px*pow(H,4) + 24360*G4p*G4px*G5p*G5px*pow(H,4) - 1872*G3p*G4x*G5p*G5px*pow(H,4) + 16560*G4pp*G4x*G5p*G5px*pow(H,4) + 552*G3p*pow(G5p,2)*G5px*pow(H,4) - 
   7584*G4pp*pow(G5p,2)*G5px*pow(H,4) + 3228*G3x*G4*G5pp*G5px*pow(H,4) - 10488*G4*G4px*G5pp*G5px*pow(H,4) + 5712*G4p*G4x*G5pp*G5px*pow(H,4) - 5364*G4p*G5p*G5pp*G5px*pow(H,4) + 
   576*G4*pow(G5pp,2)*G5px*pow(H,4) - 864*G4*G4x*G5ppp*G5px*pow(H,4) + 1584*G4*G5p*G5ppp*G5px*pow(H,4) - 720*pow(G4,2)*G5pppx*G5px*pow(H,4) + 360*G4*G4p*G5ppx*G5px*pow(H,4) + 
   1848*G3p*G4*pow(G5px,2)*pow(H,4) + 2304*pow(G4p,2)*pow(G5px,2)*pow(H,4) - 888*G4*G4pp*pow(G5px,2)*pow(H,4) - 48*G2xx*pow(G4,2)*G5pxx*pow(H,4) + 
   912*G3px*pow(G4,2)*G5pxx*pow(H,4) - 264*G3x*G4*G4p*G5pxx*pow(H,4) - 1824*pow(G4,2)*G4ppx*G5pxx*pow(H,4) - 336*G4*G4p*G4px*G5pxx*pow(H,4) + 960*G3p*G4*G4x*G5pxx*pow(H,4) - 
   1560*pow(G4p,2)*G4x*G5pxx*pow(H,4) + 288*G4*G4pp*G4x*G5pxx*pow(H,4) - 384*G3p*G4*G5p*G5pxx*pow(H,4) + 1296*pow(G4p,2)*G5p*G5pxx*pow(H,4) + 288*G4*G4pp*G5p*G5pxx*pow(H,4) + 
   600*G4*G4p*G5pp*G5pxx*pow(H,4) - 144*pow(G4,2)*G5ppp*G5pxx*pow(H,4) - 96*G3p*pow(G4,2)*G5pxxx*pow(H,4) - 48*G4*pow(G4p,2)*G5pxxx*pow(H,4) - 216*G3px*G3x*G4*G5x*pow(H,4) - 
   1566*pow(G3x,2)*G4p*G5x*pow(H,4) + 1008*G3pxx*G4*G4p*G5x*pow(H,4) + 648*G3xx*pow(G4p,2)*G5x*pow(H,4) - 360*G3xx*G4*G4pp*G5x*pow(H,4) + 1008*G3x*G4*G4ppx*G5x*pow(H,4) - 
   1152*G4*G4p*G4ppxx*G5x*pow(H,4) + 216*G2xx*G4*G4px*G5x*pow(H,4) + 288*G3px*G4*G4px*G5x*pow(H,4) + 9144*G3x*G4p*G4px*G5x*pow(H,4) - 2160*G4*G4ppx*G4px*G5x*pow(H,4) - 
   12744*G4p*pow(G4px,2)*G5x*pow(H,4) + 864*G3p*G4*G4pxx*G5x*pow(H,4) - 1872*pow(G4p,2)*G4pxx*G5x*pow(H,4) - 1296*G4*G4pp*G4pxx*G5x*pow(H,4) - 936*G2pxx*G4*G4x*G5x*pow(H,4) + 
   936*G3ppx*G4*G4x*G5x*pow(H,4) - 720*G2xx*G4p*G4x*G5x*pow(H,4) - 1008*G3px*G4p*G4x*G5x*pow(H,4) + 7344*G3x*G4pp*G4x*G5x*pow(H,4) + 10512*G4p*G4ppx*G4x*G5x*pow(H,4) + 
   1152*G3p*G4px*G4x*G5x*pow(H,4) - 18720*G4pp*G4px*G4x*G5x*pow(H,4) + 1512*G2px*pow(G4x,2)*G5x*pow(H,4) - 3672*G3pp*pow(G4x,2)*G5x*pow(H,4) + 1296*G4ppp*pow(G4x,2)*G5x*pow(H,4) - 
   360*G2px*G4*G4xx*G5x*pow(H,4) - 432*G3pp*G4*G4xx*G5x*pow(H,4) + 5040*G3p*G4p*G4xx*G5x*pow(H,4) - 10800*G4p*G4pp*G4xx*G5x*pow(H,4) + 2304*G4*G4ppp*G4xx*G5x*pow(H,4) + 
   792*G2pxx*G4*G5p*G5x*pow(H,4) - 504*G3ppx*G4*G5p*G5x*pow(H,4) + 612*G2xx*G4p*G5p*G5x*pow(H,4) + 180*G3px*G4p*G5p*G5x*pow(H,4) - 6480*G3x*G4pp*G5p*G5x*pow(H,4) - 
   576*G4*G4pppx*G5p*G5x*pow(H,4) - 7920*G4p*G4ppx*G5p*G5x*pow(H,4) - 288*G3p*G4px*G5p*G5x*pow(H,4) + 16272*G4pp*G4px*G5p*G5x*pow(H,4) - 2556*G2px*G4x*G5p*G5x*pow(H,4) + 
   6660*G3pp*G4x*G5p*G5x*pow(H,4) - 3096*G4ppp*G4x*G5p*G5x*pow(H,4) + 1008*G2px*pow(G5p,2)*G5x*pow(H,4) - 2916*G3pp*pow(G5p,2)*G5x*pow(H,4) + 1800*G4ppp*pow(G5p,2)*G5x*pow(H,4) - 
   324*G2xx*G4*G5pp*G5x*pow(H,4) - 468*G3px*G4*G5pp*G5x*pow(H,4) - 2268*G3x*G4p*G5pp*G5x*pow(H,4) + 1728*G4*G4ppx*G5pp*G5x*pow(H,4) + 6336*G4p*G4px*G5pp*G5x*pow(H,4) - 
   3240*G3p*G4x*G5pp*G5x*pow(H,4) + 1944*G4pp*G4x*G5pp*G5x*pow(H,4) + 2268*G3p*G5p*G5pp*G5x*pow(H,4) - 1944*G4pp*G5p*G5pp*G5x*pow(H,4) - 504*G4p*pow(G5pp,2)*G5x*pow(H,4) + 
   504*G3x*G4*G5ppp*G5x*pow(H,4) - 1152*G4*G4px*G5ppp*G5x*pow(H,4) - 1440*G4p*G4x*G5ppp*G5x*pow(H,4) + 1800*G4p*G5p*G5ppp*G5x*pow(H,4) - 432*G4*G4p*G5pppx*G5x*pow(H,4) - 
   928*G3p*G4*G5ppx*G5x*pow(H,4) - 60*pow(G4p,2)*G5ppx*G5x*pow(H,4) + 1224*G4*G4pp*G5ppx*G5x*pow(H,4) + 160*G2px*G4*G5px*G5x*pow(H,4) + 484*G3pp*G4*G5px*G5x*pow(H,4) - 
   3048*G3p*G4p*G5px*G5x*pow(H,4) + 6024*G4p*G4pp*G5px*G5x*pow(H,4) - 1368*G4*G4ppp*G5px*G5x*pow(H,4) - 12*G2ppx*G4*pow(G5x,2)*pow(H,4) + 12*G3ppp*G4*pow(G5x,2)*pow(H,4) - 
   444*G2px*G4p*pow(G5x,2)*pow(H,4) + 1194*G3pp*G4p*pow(G5x,2)*pow(H,4) - 924*G3p*G4pp*pow(G5x,2)*pow(H,4) + 468*pow(G4pp,2)*pow(G5x,2)*pow(H,4) - 
   828*G4p*G4ppp*pow(G5x,2)*pow(H,4) - 132*G2pp*G4x*pow(G5x,2)*pow(H,4) + 126*G2pp*G5p*pow(G5x,2)*pow(H,4) + 48*G2pxx*pow(G4,2)*G5xx*pow(H,4) - 384*G3ppx*pow(G4,2)*G5xx*pow(H,4) + 
   72*G2xx*G4*G4p*G5xx*pow(H,4) + 48*G3px*G4*G4p*G5xx*pow(H,4) + 1176*G3x*pow(G4p,2)*G5xx*pow(H,4) - 312*G3x*G4*G4pp*G5xx*pow(H,4) + 672*pow(G4,2)*G4pppx*G5xx*pow(H,4) + 
   144*G4*G4p*G4ppx*G5xx*pow(H,4) - 1056*G3p*G4*G4px*G5xx*pow(H,4) - 3336*pow(G4p,2)*G4px*G5xx*pow(H,4) + 1680*G4*G4pp*G4px*G5xx*pow(H,4) + 168*G2px*G4*G4x*G5xx*pow(H,4) - 
   744*G3pp*G4*G4x*G5xx*pow(H,4) + 1776*G3p*G4p*G4x*G5xx*pow(H,4) - 3504*G4p*G4pp*G4x*G5xx*pow(H,4) + 816*G4*G4ppp*G4x*G5xx*pow(H,4) + 120*G2px*G4*G5p*G5xx*pow(H,4) + 
   168*G3pp*G4*G5p*G5xx*pow(H,4) - 1704*G3p*G4p*G5p*G5xx*pow(H,4) + 3528*G4p*G4pp*G5p*G5xx*pow(H,4) - 816*G4*G4ppp*G5p*G5xx*pow(H,4) + 1128*G3p*G4*G5pp*G5xx*pow(H,4) + 
   420*pow(G4p,2)*G5pp*G5xx*pow(H,4) - 720*G4*G4pp*G5pp*G5xx*pow(H,4) - 480*G4*G4p*G5ppp*G5xx*pow(H,4) - 40*G2pp*G4*G5x*G5xx*pow(H,4) - 48*G2px*pow(G4,2)*G5xxx*pow(H,4) + 
   96*G3pp*pow(G4,2)*G5xxx*pow(H,4) - 48*pow(G4p,3)*G5xxx*pow(H,4) - 48*G4*G4p*G4pp*G5xxx*pow(H,4) - 2736*pow(G4x,3)*G5px*pow(H,6) + 10512*G4*G4x*G4xx*G5px*pow(H,6) - 
   864*pow(G4,2)*G4xxx*G5px*pow(H,6) + 7848*pow(G4x,2)*G5p*G5px*pow(H,6) - 8784*G4*G4xx*G5p*G5px*pow(H,6) - 7272*G4x*pow(G5p,2)*G5px*pow(H,6) + 2160*pow(G5p,3)*G5px*pow(H,6) - 
   6216*G4*G4x*pow(G5px,2)*pow(H,6) + 4728*G4*G5p*pow(G5px,2)*pow(H,6) - 3888*G4*pow(G4x,2)*G5pxx*pow(H,6) + 864*pow(G4,2)*G4xx*G5pxx*pow(H,6) + 6048*G4*G4x*G5p*G5pxx*pow(H,6) - 
   2160*G4*pow(G5p,2)*G5pxx*pow(H,6) + 464*pow(G4,2)*G5px*G5pxx*pow(H,6) + 288*pow(G4,2)*G4x*G5pxxx*pow(H,6) - 288*pow(G4,2)*G5p*G5pxxx*pow(H,6) + 864*pow(G4,2)*G4pxxx*G5x*pow(H,6) - 
   11232*G4*G4pxx*G4x*G5x*pow(H,6) - 6480*G4px*pow(G4x,2)*G5x*pow(H,6) + 4752*G4*G4px*G4xx*G5x*pow(H,6) - 23760*G4p*G4x*G4xx*G5x*pow(H,6) + 1296*G4*G4p*G4xxx*G5x*pow(H,6) + 
   9504*G4*G4pxx*G5p*G5x*pow(H,6) + 11664*G4px*G4x*G5p*G5x*pow(H,6) + 20304*G4p*G4xx*G5p*G5x*pow(H,6) - 5400*G4px*pow(G5p,2)*G5x*pow(H,6) + 9936*pow(G4x,2)*G5pp*G5x*pow(H,6) - 
   8208*G4*G4xx*G5pp*G5x*pow(H,6) - 17820*G4x*G5p*G5pp*G5x*pow(H,6) + 8100*pow(G5p,2)*G5pp*G5x*pow(H,6) + 7920*G4*G4x*G5ppx*G5x*pow(H,6) - 6432*G4*G5p*G5ppx*G5x*pow(H,6) - 
   528*pow(G4,2)*G5ppxx*G5x*pow(H,6) + 1800*G3x*G4*G5px*G5x*pow(H,6) - 7392*G4*G4px*G5px*G5x*pow(H,6) + 18288*G4p*G4x*G5px*G5x*pow(H,6) - 15876*G4p*G5p*G5px*G5x*pow(H,6) + 
   4404*G4*G5pp*G5px*G5x*pow(H,6) + 696*G4*G4p*G5pxx*G5x*pow(H,6) - 756*G3px*G4*pow(G5x,2)*pow(H,6) - 2538*G3x*G4p*pow(G5x,2)*pow(H,6) + 2160*G4*G4ppx*pow(G5x,2)*pow(H,6) + 
   9072*G4p*G4px*pow(G5x,2)*pow(H,6) + 7236*G4pp*G4x*pow(G5x,2)*pow(H,6) - 6876*G4pp*G5p*pow(G5x,2)*pow(H,6) - 2430*G4p*G5pp*pow(G5x,2)*pow(H,6) + 468*G4*G5ppp*pow(G5x,2)*pow(H,6) - 
   864*pow(G4,2)*G4pxx*G5xx*pow(H,6) + 3744*G4*G4px*G4x*G5xx*pow(H,6) - 4896*G4p*pow(G4x,2)*G5xx*pow(H,6) + 3744*G4*G4p*G4xx*G5xx*pow(H,6) - 2016*G4*G4px*G5p*G5xx*pow(H,6) + 
   8568*G4p*G4x*G5p*G5xx*pow(H,6) - 3672*G4p*pow(G5p,2)*G5xx*pow(H,6) - 4968*G4*G4x*G5pp*G5xx*pow(H,6) + 3240*G4*G5p*G5pp*G5xx*pow(H,6) + 64*pow(G4,2)*G5ppx*G5xx*pow(H,6) - 
   2832*G4*G4p*G5px*G5xx*pow(H,6) + 1584*pow(G4p,2)*G5x*G5xx*pow(H,6) - 696*G4*G4pp*G5x*G5xx*pow(H,6) - 288*pow(G4,2)*G4px*G5xxx*pow(H,6) + 288*G4*G4p*G4x*G5xxx*pow(H,6) - 
   288*G4*G4p*G5p*G5xxx*pow(H,6) + 288*pow(G4,2)*G5pp*G5xxx*pow(H,6) + 540*G4*G5px*pow(G5x,2)*pow(H,8) - 810*G4p*pow(G5x,3)*pow(H,8) + 
   6*pow(G2x,2)*(-4*pow(G4px,2) - 2*G3px*G4x + 4*G4ppx*G4x + G3px*G5p - 2*G4ppx*G5p + G3x*(2*G4px - G5pp) - 2*G4x*G5px*pow(H,2) - G5p*G5px*pow(H,2) + 2*G4*G5pxx*pow(H,2) + G5pp*G5x*pow(H,2) - 
      2*G4p*G5xx*pow(H,2) + 2*G4px*(G5pp + G5x*pow(H,2))) + G2p*(9*pow(G3x,3) - 120*pow(G4px,3) - 16*G2pxx*pow(G4x,2) + 16*G3ppx*pow(G4x,2) + 16*G2pxx*G4x*G5p - 16*G3ppx*G4x*G5p - 
      4*G2pxx*pow(G5p,2) + 4*G3ppx*pow(G5p,2) + 192*G4*G4pxxx*G4x*pow(H,2) - 144*G3xx*pow(G4x,2)*pow(H,2) + 288*G4pxx*pow(G4x,2)*pow(H,2) + 144*G3xx*G4*G4xx*pow(H,2) - 
      672*G4*G4pxx*G4xx*pow(H,2) - 672*G4p*pow(G4xx,2)*pow(H,2) + 144*G4p*G4x*G4xxx*pow(H,2) - 96*G4*G4pxxx*G5p*pow(H,2) + 216*G3xx*G4x*G5p*pow(H,2) - 624*G4pxx*G4x*G5p*pow(H,2) - 
      72*G4p*G4xxx*G5p*pow(H,2) - 72*G3xx*pow(G5p,2)*pow(H,2) + 240*G4pxx*pow(G5p,2)*pow(H,2) - 192*G4x*G4xx*G5pp*pow(H,2) + 48*pow(G4x,2)*G5ppx*pow(H,2) + 192*G4*G4xx*G5ppx*pow(H,2) + 
      48*G4x*G5p*G5ppx*pow(H,2) - 36*pow(G5p,2)*G5ppx*pow(H,2) - 96*G4*G4x*G5ppxx*pow(H,2) + 48*G4*G5p*G5ppxx*pow(H,2) - 72*G3xx*G4*G5px*pow(H,2) + 336*G4*G4pxx*G5px*pow(H,2) + 
      768*G4p*G4xx*G5px*pow(H,2) + 96*G4x*G5pp*G5px*pow(H,2) - 96*G4*G5ppx*G5px*pow(H,2) - 216*G4p*pow(G5px,2)*pow(H,2) - 120*G4p*G4x*G5pxx*pow(H,2) + 60*G4p*G5p*G5pxx*pow(H,2) + 
      48*G3pxx*G4*G5x*pow(H,2) + 36*G3xx*G4p*G5x*pow(H,2) - 96*G4*G4ppxx*G5x*pow(H,2) - 168*G4p*G4pxx*G5x*pow(H,2) - 60*G2xx*G4x*G5x*pow(H,2) + 84*G3px*G4x*G5x*pow(H,2) + 
      48*G4ppx*G4x*G5x*pow(H,2) - 24*G2x*G4xx*G5x*pow(H,2) + 48*G3p*G4xx*G5x*pow(H,2) - 96*G4pp*G4xx*G5x*pow(H,2) + 42*G2xx*G5p*G5x*pow(H,2) - 78*G3px*G5p*G5x*pow(H,2) + 
      24*G4ppx*G5p*G5x*pow(H,2) + 48*G4p*G5ppx*G5x*pow(H,2) + 12*G2x*G5px*G5x*pow(H,2) - 24*G3p*G5px*G5x*pow(H,2) + 48*G4pp*G5px*G5x*pow(H,2) - 4*G2px*pow(G5x,2)*pow(H,2) + 
      8*G3pp*pow(G5x,2)*pow(H,2) + 12*G2xx*G4*G5xx*pow(H,2) - 36*G3px*G4*G5xx*pow(H,2) + 48*G4*G4ppx*G5xx*pow(H,2) - 12*G2x*G4x*G5xx*pow(H,2) + 24*G3p*G4x*G5xx*pow(H,2) - 
      48*G4pp*G4x*G5xx*pow(H,2) + 6*G2x*G5p*G5xx*pow(H,2) - 12*G3p*G5p*G5xx*pow(H,2) + 24*G4pp*G5p*G5xx*pow(H,2) - 24*G4p*G5pp*G5xx*pow(H,2) + 720*G4x*G4xx*G5x*pow(H,4) + 
      216*G4*G4xxx*G5x*pow(H,4) - 1296*G4xx*G5p*G5x*pow(H,4) - 60*G4x*G5px*G5x*pow(H,4) + 462*G5p*G5px*G5x*pow(H,4) - 68*G4*G5pxx*G5x*pow(H,4) - 66*G5pp*pow(G5x,2)*pow(H,4) + 
      24*pow(G4x,2)*G5xx*pow(H,4) + 624*G4*G4xx*G5xx*pow(H,4) - 252*G4x*G5p*G5xx*pow(H,4) + 228*pow(G5p,2)*G5xx*pow(H,4) - 420*G4*G5px*G5xx*pow(H,4) - 168*G4p*G5x*G5xx*pow(H,4) + 
      48*G4*G4x*G5xxx*pow(H,4) - 48*G4*G5p*G5xxx*pow(H,4) + 135*pow(G5x,3)*pow(H,6) + 12*pow(G4px,2)*(2*G5pp + 13*G5x*pow(H,2)) + pow(G3x,2)*(-66*G4px + 6*G5pp + 45*G5x*pow(H,2)) + 
      6*G4px*(4*G2xx*G4x - 12*G3px*G4x + 16*G4ppx*G4x - 2*G2xx*G5p + 6*G3px*G5p - 8*G4ppx*G5p - 48*G4x*G4xx*pow(H,2) - 24*G4*G4xxx*pow(H,2) + 200*G4xx*G5p*pow(H,2) + 4*G4x*G5px*pow(H,2) - 
         98*G5p*G5px*pow(H,2) + 20*G4*G5pxx*pow(H,2) + 4*G5pp*G5x*pow(H,2) + 36*G4p*G5xx*pow(H,2) - 13*pow(G5x,2)*pow(H,4)) + 
      3*G3x*(52*pow(G4px,2) - 4*G2xx*G4x + 12*G3px*G4x - 16*G4ppx*G4x + 2*G2xx*G5p - 6*G3px*G5p + 8*G4ppx*G5p + 96*G4x*G4xx*pow(H,2) + 24*G4*G4xxx*pow(H,2) - 192*G4xx*G5p*pow(H,2) - 
         28*G4x*G5px*pow(H,2) + 94*G5p*G5px*pow(H,2) - 20*G4*G5pxx*pow(H,2) - 4*G5pp*G5x*pow(H,2) - 28*G4p*G5xx*pow(H,2) + 33*pow(G5x,2)*pow(H,4) - 8*G4px*(G5pp + 7*G5x*pow(H,2)))) - 
   9*G3px*pow(G3x,2)*rptot + 18*pow(G3x,2)*G4ppx*rptot + 12*G2xx*G3x*G4px*rptot + 24*G3px*G3x*G4px*rptot - 72*G3x*G4ppx*G4px*rptot - 24*G2xx*pow(G4px,2)*rptot - 12*G3px*pow(G4px,2)*rptot + 
   72*G4ppx*pow(G4px,2)*rptot - 12*G2xx*G3px*G4x*rptot + 12*pow(G3px,2)*G4x*rptot + 12*G2pxx*G3x*G4x*rptot - 12*G3ppx*G3x*G4x*rptot + 24*G2xx*G4ppx*G4x*rptot - 24*G3px*G4ppx*G4x*rptot - 
   24*G2pxx*G4px*G4x*rptot + 24*G3ppx*G4px*G4x*rptot + 6*G2xx*G3px*G5p*rptot - 6*pow(G3px,2)*G5p*rptot - 6*G2pxx*G3x*G5p*rptot + 6*G3ppx*G3x*G5p*rptot - 12*G2xx*G4ppx*G5p*rptot + 
   12*G3px*G4ppx*G5p*rptot + 12*G2pxx*G4px*G5p*rptot - 12*G3ppx*G4px*G5p*rptot - 6*G2xx*G3x*G5pp*rptot + 6*G3px*G3x*G5pp*rptot + 12*G2xx*G4px*G5pp*rptot - 12*G3px*G4px*G5pp*rptot + 
   144*G3xx*G4*G4pxx*pow(H,2)*rptot - 288*G4*pow(G4pxx,2)*pow(H,2)*rptot - 72*G3x*G4*G4pxxx*pow(H,2)*rptot + 144*G4*G4px*G4pxxx*pow(H,2)*rptot - 288*G3x*G4pxx*G4x*pow(H,2)*rptot + 
   864*G4px*G4pxx*G4x*pow(H,2)*rptot - 144*G4p*G4pxxx*G4x*pow(H,2)*rptot + 144*G3pxx*pow(G4x,2)*pow(H,2)*rptot - 288*G4ppxx*pow(G4x,2)*pow(H,2)*rptot - 144*G3pxx*G4*G4xx*pow(H,2)*rptot - 
   144*G3xx*G4p*G4xx*pow(H,2)*rptot + 288*G4*G4ppxx*G4xx*pow(H,2)*rptot + 288*G3x*G4px*G4xx*pow(H,2)*rptot - 864*pow(G4px,2)*G4xx*pow(H,2)*rptot + 864*G4p*G4pxx*G4xx*pow(H,2)*rptot - 
   288*G3px*G4x*G4xx*pow(H,2)*rptot + 288*G4ppx*G4x*G4xx*pow(H,2)*rptot + 288*G4pp*pow(G4xx,2)*pow(H,2)*rptot + 72*G3px*G4*G4xxx*pow(H,2)*rptot - 72*G3x*G4p*G4xxx*pow(H,2)*rptot - 
   144*G4*G4ppx*G4xxx*pow(H,2)*rptot + 144*G4pp*G4x*G4xxx*pow(H,2)*rptot - 72*G3xx*G4px*G5p*pow(H,2)*rptot + 288*G3x*G4pxx*G5p*pow(H,2)*rptot - 576*G4px*G4pxx*G5p*pow(H,2)*rptot + 
   72*G4p*G4pxxx*G5p*pow(H,2)*rptot - 216*G3pxx*G4x*G5p*pow(H,2)*rptot + 432*G4ppxx*G4x*G5p*pow(H,2)*rptot + 288*G3px*G4xx*G5p*pow(H,2)*rptot - 432*G4ppx*G4xx*G5p*pow(H,2)*rptot - 
   72*G4pp*G4xxx*G5p*pow(H,2)*rptot + 72*G3pxx*pow(G5p,2)*pow(H,2)*rptot - 144*G4ppxx*pow(G5p,2)*pow(H,2)*rptot + 72*G3xx*G4x*G5pp*pow(H,2)*rptot - 144*G4pxx*G4x*G5pp*pow(H,2)*rptot + 
   144*G4px*G4xx*G5pp*pow(H,2)*rptot + 72*G4p*G4xxx*G5pp*pow(H,2)*rptot - 72*G3xx*G4*G5ppx*pow(H,2)*rptot + 144*G4*G4pxx*G5ppx*pow(H,2)*rptot + 108*G3x*G4x*G5ppx*pow(H,2)*rptot - 
   360*G4px*G4x*G5ppx*pow(H,2)*rptot - 288*G4p*G4xx*G5ppx*pow(H,2)*rptot - 126*G3x*G5p*G5ppx*pow(H,2)*rptot + 324*G4px*G5p*G5ppx*pow(H,2)*rptot + 36*G3x*G4*G5ppxx*pow(H,2)*rptot - 
   72*G4*G4px*G5ppxx*pow(H,2)*rptot + 72*G4p*G4x*G5ppxx*pow(H,2)*rptot - 36*G4p*G5p*G5ppxx*pow(H,2)*rptot + 9*pow(G3x,2)*G5px*pow(H,2)*rptot + 72*G3pxx*G4*G5px*pow(H,2)*rptot + 
   36*G3xx*G4p*G5px*pow(H,2)*rptot - 144*G4*G4ppxx*G5px*pow(H,2)*rptot - 252*G3x*G4px*G5px*pow(H,2)*rptot + 612*pow(G4px,2)*G5px*pow(H,2)*rptot - 360*G4p*G4pxx*G5px*pow(H,2)*rptot - 
   12*G2xx*G4x*G5px*pow(H,2)*rptot + 192*G3px*G4x*G5px*pow(H,2)*rptot - 216*G4ppx*G4x*G5px*pow(H,2)*rptot - 48*G3p*G4xx*G5px*pow(H,2)*rptot - 288*G4pp*G4xx*G5px*pow(H,2)*rptot - 
   6*G2xx*G5p*G5px*pow(H,2)*rptot - 156*G3px*G5p*G5px*pow(H,2)*rptot + 252*G4ppx*G5p*G5px*pow(H,2)*rptot + 18*G3x*G5pp*G5px*pow(H,2)*rptot - 108*G4px*G5pp*G5px*pow(H,2)*rptot + 
   144*G4p*G5ppx*G5px*pow(H,2)*rptot + 24*G3p*pow(G5px,2)*pow(H,2)*rptot + 72*G4pp*pow(G5px,2)*pow(H,2)*rptot + 12*G2xx*G4*G5pxx*pow(H,2)*rptot - 48*G3px*G4*G5pxx*pow(H,2)*rptot + 
   72*G3x*G4p*G5pxx*pow(H,2)*rptot + 72*G4*G4ppx*G5pxx*pow(H,2)*rptot - 72*G4p*G4px*G5pxx*pow(H,2)*rptot + 24*G3p*G4x*G5pxx*pow(H,2)*rptot - 72*G4pp*G4x*G5pxx*pow(H,2)*rptot - 
   12*G3p*G5p*G5pxx*pow(H,2)*rptot + 36*G4pp*G5p*G5pxx*pow(H,2)*rptot - 36*G4p*G5pp*G5pxx*pow(H,2)*rptot - 54*G3px*G3x*G5x*pow(H,2)*rptot - 36*G3pxx*G4p*G5x*pow(H,2)*rptot + 
   36*G3xx*G4pp*G5x*pow(H,2)*rptot + 72*G3x*G4ppx*G5x*pow(H,2)*rptot + 72*G4p*G4ppxx*G5x*pow(H,2)*rptot + 12*G2xx*G4px*G5x*pow(H,2)*rptot + 132*G3px*G4px*G5x*pow(H,2)*rptot - 
   216*G4ppx*G4px*G5x*pow(H,2)*rptot + 48*G3p*G4pxx*G5x*pow(H,2)*rptot - 72*G4pp*G4pxx*G5x*pow(H,2)*rptot + 60*G2pxx*G4x*G5x*pow(H,2)*rptot - 60*G3ppx*G4x*G5x*pow(H,2)*rptot + 
   24*G2px*G4xx*G5x*pow(H,2)*rptot - 48*G3pp*G4xx*G5x*pow(H,2)*rptot - 42*G2pxx*G5p*G5x*pow(H,2)*rptot + 42*G3ppx*G5p*G5x*pow(H,2)*rptot + 6*G2xx*G5pp*G5x*pow(H,2)*rptot - 
   6*G3px*G5pp*G5x*pow(H,2)*rptot - 24*G3p*G5ppx*G5x*pow(H,2)*rptot - 12*G2px*G5px*G5x*pow(H,2)*rptot + 24*G3pp*G5px*G5x*pow(H,2)*rptot - 12*G2pxx*G4*G5xx*pow(H,2)*rptot + 
   12*G3ppx*G4*G5xx*pow(H,2)*rptot - 12*G2xx*G4p*G5xx*pow(H,2)*rptot + 48*G3px*G4p*G5xx*pow(H,2)*rptot + 36*G3x*G4pp*G5xx*pow(H,2)*rptot - 72*G4p*G4ppx*G5xx*pow(H,2)*rptot - 
   24*G3p*G4px*G5xx*pow(H,2)*rptot - 72*G4pp*G4px*G5xx*pow(H,2)*rptot + 12*G2px*G4x*G5xx*pow(H,2)*rptot - 24*G3pp*G4x*G5xx*pow(H,2)*rptot - 6*G2px*G5p*G5xx*pow(H,2)*rptot + 
   12*G3pp*G5p*G5xx*pow(H,2)*rptot + 12*G3p*G5pp*G5xx*pow(H,2)*rptot - 720*G4x*G4xx*G5px*pow(H,4)*rptot + 216*G4*G4xxx*G5px*pow(H,4)*rptot + 432*G4xx*G5p*G5px*pow(H,4)*rptot + 
   396*G4x*pow(G5px,2)*pow(H,4)*rptot - 198*G5p*pow(G5px,2)*pow(H,4)*rptot + 120*pow(G4x,2)*G5pxx*pow(H,4)*rptot - 48*G4*G4xx*G5pxx*pow(H,4)*rptot - 108*G4x*G5p*G5pxx*pow(H,4)*rptot - 
   12*pow(G5p,2)*G5pxx*pow(H,4)*rptot - 120*G4*G5px*G5pxx*pow(H,4)*rptot - 48*G4*G4x*G5pxxx*pow(H,4)*rptot + 48*G4*G5p*G5pxxx*pow(H,4)*rptot - 216*G4*G4pxxx*G5x*pow(H,4)*rptot + 
   144*G4pxx*G4x*G5x*pow(H,4)*rptot - 144*G4px*G4xx*G5x*pow(H,4)*rptot - 216*G4p*G4xxx*G5x*pow(H,4)*rptot + 144*G4pxx*G5p*G5x*pow(H,4)*rptot + 720*G4xx*G5pp*G5x*pow(H,4)*rptot - 
   252*G4x*G5ppx*G5x*pow(H,4)*rptot + 54*G5p*G5ppx*G5x*pow(H,4)*rptot + 108*G4*G5ppxx*G5x*pow(H,4)*rptot - 90*G3x*G5px*G5x*pow(H,4)*rptot + 216*G4px*G5px*G5x*pow(H,4)*rptot - 
   378*G5pp*G5px*G5x*pow(H,4)*rptot + 132*G4p*G5pxx*G5x*pow(H,4)*rptot - 9*G3px*pow(G5x,2)*pow(H,4)*rptot - 90*G4ppx*pow(G5x,2)*pow(H,4)*rptot + 48*G4*G4pxx*G5xx*pow(H,4)*rptot - 
   144*G4px*G4x*G5xx*pow(H,4)*rptot - 624*G4p*G4xx*G5xx*pow(H,4)*rptot + 12*G4px*G5p*G5xx*pow(H,4)*rptot + 348*G4x*G5pp*G5xx*pow(H,4)*rptot - 216*G5p*G5pp*G5xx*pow(H,4)*rptot + 
   12*G4*G5ppx*G5xx*pow(H,4)*rptot + 372*G4p*G5px*G5xx*pow(H,4)*rptot + 192*G4pp*G5x*G5xx*pow(H,4)*rptot + 48*G4*G4px*G5xxx*pow(H,4)*rptot - 48*G4p*G4x*G5xxx*pow(H,4)*rptot + 
   48*G4p*G5p*G5xxx*pow(H,4)*rptot - 48*G4*G5pp*G5xxx*pow(H,4)*rptot - 135*G5px*pow(G5x,2)*pow(H,6)*rptot + 
   2*G2x*(-6*pow(G3px,2)*G4 - 6*G2pxx*G3x*G4 + 6*G3ppx*G3x*G4 + 9*pow(G3x,2)*G4pp - 36*G3x*G4p*G4ppx - 24*G3p*G3x*G4px + 12*G2pxx*G4*G4px - 12*G3ppx*G4*G4px - 24*G3x*G4pp*G4px + 72*G4p*G4ppx*G4px + 
      48*G3p*pow(G4px,2) + 12*G4pp*pow(G4px,2) + 6*G2px*G3x*G4x - 6*G3pp*G3x*G4x - 12*G2pxx*G4p*G4x + 12*G3ppx*G4p*G4x - 12*G3x*G4ppp*G4x - 48*G3p*G4ppx*G4x + 24*G4pp*G4ppx*G4x - 12*G2px*G4px*G4x + 
      12*G3pp*G4px*G4x + 24*G4ppp*G4px*G4x - 8*G2ppx*pow(G4x,2) + 8*G3ppp*pow(G4x,2) - 3*G2px*G3x*G5p + 3*G3pp*G3x*G5p + 6*G2pxx*G4p*G5p - 6*G3ppx*G4p*G5p + 6*G3x*G4ppp*G5p + 24*G3p*G4ppx*G5p - 
      12*G4pp*G4ppx*G5p + 6*G2px*G4px*G5p - 6*G3pp*G4px*G5p - 12*G4ppp*G4px*G5p + 8*G2ppx*G4x*G5p - 8*G3ppp*G4x*G5p - 2*G2ppx*pow(G5p,2) + 2*G3ppp*pow(G5p,2) + 12*G3p*G3x*G5pp - 6*G3x*G4pp*G5pp - 
      24*G3p*G4px*G5pp + 12*G4pp*G4px*G5pp + 72*G3xx*G4*G4px*pow(H,2) + 216*G3x*G4*G4pxx*pow(H,2) - 816*G4*G4px*G4pxx*pow(H,2) + 72*G4*G4p*G4pxxx*pow(H,2) - 216*G3pxx*G4*G4x*pow(H,2) - 
      72*G3xx*G4p*G4x*pow(H,2) + 528*G4*G4ppxx*G4x*pow(H,2) + 480*G4p*G4pxx*G4x*pow(H,2) - 48*G4ppx*pow(G4x,2)*pow(H,2) - 540*G3x*G4p*G4xx*pow(H,2) - 312*G4*G4ppx*G4xx*pow(H,2) + 
      912*G4p*G4px*G4xx*pow(H,2) + 840*G4pp*G4x*G4xx*pow(H,2) + 72*pow(G4p,2)*G4xxx*pow(H,2) - 72*G4*G4pp*G4xxx*pow(H,2) + 144*G3pxx*G4*G5p*pow(H,2) + 72*G3xx*G4p*G5p*pow(H,2) - 
      336*G4*G4ppxx*G5p*pow(H,2) - 144*G3x*G4px*G5p*pow(H,2) + 408*pow(G4px,2)*G5p*pow(H,2) - 456*G4p*G4pxx*G5p*pow(H,2) - 144*G4ppx*G4x*G5p*pow(H,2) - 564*G4pp*G4xx*G5p*pow(H,2) + 
      120*G4ppx*pow(G5p,2)*pow(H,2) - 72*G3xx*G4*G5pp*pow(H,2) + 192*G4*G4pxx*G5pp*pow(H,2) + 234*G3x*G4x*G5pp*pow(H,2) - 588*G4px*G4x*G5pp*pow(H,2) + 12*G4p*G4xx*G5pp*pow(H,2) - 
      81*G3x*G5p*G5pp*pow(H,2) + 138*G4px*G5p*G5pp*pow(H,2) + 24*G4x*pow(G5pp,2)*pow(H,2) - 24*pow(G4x,2)*G5ppp*pow(H,2) + 48*G4*G4xx*G5ppp*pow(H,2) + 48*G4x*G5p*G5ppp*pow(H,2) - 
      18*pow(G5p,2)*G5ppp*pow(H,2) - 48*G4*G4x*G5pppx*pow(H,2) + 24*G4*G5p*G5pppx*pow(H,2) - 96*G3x*G4*G5ppx*pow(H,2) + 312*G4*G4px*G5ppx*pow(H,2) - 144*G4p*G4x*G5ppx*pow(H,2) + 
      144*G4p*G5p*G5ppx*pow(H,2) - 24*G4*G5pp*G5ppx*pow(H,2) - 36*G4*G4p*G5ppxx*pow(H,2) + 252*G3x*G4p*G5px*pow(H,2) + 204*G4*G4ppx*G5px*pow(H,2) - 312*G4p*G4px*G5px*pow(H,2) + 
      24*G3p*G4x*G5px*pow(H,2) - 480*G4pp*G4x*G5px*pow(H,2) + 12*G3p*G5p*G5px*pow(H,2) + 300*G4pp*G5p*G5px*pow(H,2) - 42*G4p*G5pp*G5px*pow(H,2) - 24*G4*G5ppp*G5px*pow(H,2) - 
      24*G3p*G4*G5pxx*pow(H,2) - 54*pow(G4p,2)*G5pxx*pow(H,2) + 48*G4*G4pp*G5pxx*pow(H,2) - 42*G2pxx*G4*G5x*pow(H,2) + 66*G3ppx*G4*G5x*pow(H,2) + 132*G3x*G4pp*G5x*pow(H,2) - 
      48*G4*G4pppx*G5x*pow(H,2) - 108*G4p*G4ppx*G5x*pow(H,2) - 24*G3p*G4px*G5x*pow(H,2) - 312*G4pp*G4px*G5x*pow(H,2) + 30*G2px*G4x*G5x*pow(H,2) - 30*G3pp*G4x*G5x*pow(H,2) - 
      60*G4ppp*G4x*G5x*pow(H,2) - 21*G2px*G5p*G5x*pow(H,2) + 21*G3pp*G5p*G5x*pow(H,2) + 42*G4ppp*G5p*G5x*pow(H,2) - 12*G3p*G5pp*G5x*pow(H,2) + 18*G4pp*G5pp*G5x*pow(H,2) + 
      12*G4p*G5ppp*G5x*pow(H,2) + 2*G2pp*pow(G5x,2)*pow(H,2) - 6*G2px*G4*G5xx*pow(H,2) + 6*G3pp*G4*G5xx*pow(H,2) + 24*G3p*G4p*G5xx*pow(H,2) - 48*G4p*G4pp*G5xx*pow(H,2) + 
      12*G4*G4ppp*G5xx*pow(H,2) - 312*pow(G4x,2)*G5px*pow(H,4) + 780*G4*G4xx*G5px*pow(H,4) + 468*G4x*G5p*G5px*pow(H,4) - 138*pow(G5p,2)*G5px*pow(H,4) - 462*G4*pow(G5px,2)*pow(H,4) - 
      240*G4*G4x*G5pxx*pow(H,4) + 96*G4*G5p*G5pxx*pow(H,4) + 24*pow(G4,2)*G5pxxx*pow(H,4) - 216*G4*G4pxx*G5x*pow(H,4) - 288*G4px*G4x*G5x*pow(H,4) - 1260*G4p*G4xx*G5x*pow(H,4) + 
      72*G4px*G5p*G5x*pow(H,4) + 810*G4x*G5pp*G5x*pow(H,4) - 567*G5p*G5pp*G5x*pow(H,4) + 232*G4*G5ppx*G5x*pow(H,4) + 762*G4p*G5px*G5x*pow(H,4) + 231*G4pp*pow(G5x,2)*pow(H,4) + 
      264*G4*G4px*G5xx*pow(H,4) - 444*G4p*G4x*G5xx*pow(H,4) + 426*G4p*G5p*G5xx*pow(H,4) - 282*G4*G5pp*G5xx*pow(H,4) + 
      6*G2xx*(G3px*G4 - G3x*G4p - 2*G4*G4ppx + 2*G4pp*G4x - G4pp*G5p + G4p*G5pp + 3*G4*G5px*pow(H,2) - 3*G4p*G5x*pow(H,2)) + 
      6*G3px*(4*G3x*G4p - 6*G4p*G4px + 4*G3p*G4x - 4*G4pp*G4x - 2*G3p*G5p + 2*G4pp*G5p - G4p*G5pp - 12*pow(G4x,2)*pow(H,2) + 30*G4x*G5p*pow(H,2) - 15*pow(G5p,2)*pow(H,2) + 
         13*G4p*G5x*pow(H,2) + 2*G4*(G4ppx + (15*G4xx - 11*G5px)*pow(H,2))) + 12*G4xx*G5px*pow(H,2)*rptot - 6*pow(G5px,2)*pow(H,2)*rptot - 6*G4x*G5pxx*pow(H,2)*rptot + 
      3*G5p*G5pxx*pow(H,2)*rptot - 12*G4pxx*G5x*pow(H,2)*rptot + 6*G5ppx*G5x*pow(H,2)*rptot + 6*G4px*G5xx*pow(H,2)*rptot - 3*G5pp*G5xx*pow(H,2)*rptot)) + 
pow(a,8)*pow(dphi,9)*H*(72*G2x*G3xx*pow(G4px,2) - 144*G3p*G3xx*pow(G4px,2) + 144*G3xx*G4pp*pow(G4px,2) - 240*G2xx*pow(G4px,3) + 720*G4ppx*pow(G4px,3) - 72*G2x*G3xx*G4ppx*G4x + 
   144*G3p*G3xx*G4ppx*G4x - 144*G3xx*G4pp*G4ppx*G4x + 72*G2x*G3pxx*G4px*G4x - 144*G3p*G3pxx*G4px*G4x + 72*G3pp*G3xx*G4px*G4x + 144*G3pxx*G4pp*G4px*G4x - 144*G3xx*G4ppp*G4px*G4x + 
   336*G2xx*G4ppx*G4px*G4x - 576*pow(G4ppx,2)*G4px*G4x - 144*G2x*G4ppxx*G4px*G4x + 288*G3p*G4ppxx*G4px*G4x - 288*G4pp*G4ppxx*G4px*G4x - 240*G2pxx*pow(G4px,2)*G4x + 96*G3ppx*pow(G4px,2)*G4x + 
   288*G4pppx*pow(G4px,2)*G4x + 48*G2x*G2xx*G4pxx*G4x - 96*G2xx*G3p*G4pxx*G4x + 96*G2xx*G4pp*G4pxx*G4x + 144*G2x*G4ppx*G4pxx*G4x - 288*G3p*G4ppx*G4pxx*G4x + 288*G4pp*G4ppx*G4pxx*G4x - 
   192*G2px*G4px*G4pxx*G4x + 48*G3pp*G4px*G4pxx*G4x + 288*G4ppp*G4px*G4pxx*G4x + 48*G2xx*G3ppx*pow(G4x,2) - 48*G2px*G3pxx*pow(G4x,2) + 48*G3pp*G3pxx*pow(G4x,2) + 48*G2ppx*G3xx*pow(G4x,2) - 
   48*G3ppp*G3xx*pow(G4x,2) - 96*G2xx*G4pppx*pow(G4x,2) + 96*G2pxx*G4ppx*pow(G4x,2) - 96*G3ppx*G4ppx*pow(G4x,2) + 96*G2px*G4ppxx*pow(G4x,2) - 96*G3pp*G4ppxx*pow(G4x,2) - 
   96*G2ppx*G4pxx*pow(G4x,2) + 96*G3ppp*G4pxx*pow(G4x,2) - 48*G2x*G2xx*G4px*G4xx + 96*G2xx*G3p*G4px*G4xx - 96*G2xx*G4pp*G4px*G4xx + 288*G2x*G4ppx*G4px*G4xx - 576*G3p*G4ppx*G4px*G4xx + 
   576*G4pp*G4ppx*G4px*G4xx + 192*G2px*pow(G4px,2)*G4xx + 240*G3pp*pow(G4px,2)*G4xx - 864*G4ppp*pow(G4px,2)*G4xx - 48*G2pxx*G2x*G4x*G4xx + 96*G2pxx*G3p*G4x*G4xx - 48*G2xx*G3pp*G4x*G4xx + 
   48*G2x*G3ppx*G4x*G4xx - 96*G3p*G3ppx*G4x*G4xx - 96*G2pxx*G4pp*G4x*G4xx + 96*G3ppx*G4pp*G4x*G4xx + 96*G2xx*G4ppp*G4x*G4xx - 192*G2px*G4ppx*G4x*G4xx + 192*G3pp*G4ppx*G4x*G4xx + 192*G2ppx*G4px*G4x*G4xx - 
   192*G3ppp*G4px*G4x*G4xx + 36*G2x*G3xx*G4ppx*G5p - 72*G3p*G3xx*G4ppx*G5p + 72*G3xx*G4pp*G4ppx*G5p - 36*G2x*G3pxx*G4px*G5p + 72*G3p*G3pxx*G4px*G5p - 36*G3pp*G3xx*G4px*G5p - 72*G3pxx*G4pp*G4px*G5p + 
   72*G3xx*G4ppp*G4px*G5p - 168*G2xx*G4ppx*G4px*G5p + 288*pow(G4ppx,2)*G4px*G5p + 72*G2x*G4ppxx*G4px*G5p - 144*G3p*G4ppxx*G4px*G5p + 144*G4pp*G4ppxx*G4px*G5p + 120*G2pxx*pow(G4px,2)*G5p - 
   48*G3ppx*pow(G4px,2)*G5p - 144*G4pppx*pow(G4px,2)*G5p - 24*G2x*G2xx*G4pxx*G5p + 48*G2xx*G3p*G4pxx*G5p - 48*G2xx*G4pp*G4pxx*G5p - 72*G2x*G4ppx*G4pxx*G5p + 144*G3p*G4ppx*G4pxx*G5p - 
   144*G4pp*G4ppx*G4pxx*G5p + 96*G2px*G4px*G4pxx*G5p - 24*G3pp*G4px*G4pxx*G5p - 144*G4ppp*G4px*G4pxx*G5p - 48*G2xx*G3ppx*G4x*G5p + 48*G2px*G3pxx*G4x*G5p - 48*G3pp*G3pxx*G4x*G5p - 48*G2ppx*G3xx*G4x*G5p + 
   48*G3ppp*G3xx*G4x*G5p + 96*G2xx*G4pppx*G4x*G5p - 96*G2pxx*G4ppx*G4x*G5p + 96*G3ppx*G4ppx*G4x*G5p - 96*G2px*G4ppxx*G4x*G5p + 96*G3pp*G4ppxx*G4x*G5p + 96*G2ppx*G4pxx*G4x*G5p - 96*G3ppp*G4pxx*G4x*G5p + 
   24*G2pxx*G2x*G4xx*G5p - 48*G2pxx*G3p*G4xx*G5p + 24*G2xx*G3pp*G4xx*G5p - 24*G2x*G3ppx*G4xx*G5p + 48*G3p*G3ppx*G4xx*G5p + 48*G2pxx*G4pp*G4xx*G5p - 48*G3ppx*G4pp*G4xx*G5p - 48*G2xx*G4ppp*G4xx*G5p + 
   96*G2px*G4ppx*G4xx*G5p - 96*G3pp*G4ppx*G4xx*G5p - 96*G2ppx*G4px*G4xx*G5p + 96*G3ppp*G4px*G4xx*G5p + 12*G2xx*G3ppx*pow(G5p,2) - 12*G2px*G3pxx*pow(G5p,2) + 12*G3pp*G3pxx*pow(G5p,2) + 
   12*G2ppx*G3xx*pow(G5p,2) - 12*G3ppp*G3xx*pow(G5p,2) - 24*G2xx*G4pppx*pow(G5p,2) + 24*G2pxx*G4ppx*pow(G5p,2) - 24*G3ppx*G4ppx*pow(G5p,2) + 24*G2px*G4ppxx*pow(G5p,2) - 
   24*G3pp*G4ppxx*pow(G5p,2) - 24*G2ppx*G4pxx*pow(G5p,2) + 24*G3ppp*G4pxx*pow(G5p,2) - 36*G2x*G3xx*G4px*G5pp + 72*G3p*G3xx*G4px*G5pp - 72*G3xx*G4pp*G4px*G5pp + 168*G2xx*pow(G4px,2)*G5pp - 
   288*G4ppx*pow(G4px,2)*G5pp + 72*G2x*G4px*G4pxx*G5pp - 144*G3p*G4px*G4pxx*G5pp + 144*G4pp*G4px*G4pxx*G5pp - 48*G2xx*G4ppx*G4x*G5pp + 48*G2pxx*G4px*G4x*G5pp - 48*G3ppx*G4px*G4x*G5pp + 
   24*G2x*G2xx*G4xx*G5pp - 48*G2xx*G3p*G4xx*G5pp + 48*G2xx*G4pp*G4xx*G5pp - 96*G2px*G4px*G4xx*G5pp + 96*G3pp*G4px*G4xx*G5pp + 24*G2xx*G4ppx*G5p*G5pp - 24*G2pxx*G4px*G5p*G5pp + 24*G3ppx*G4px*G5p*G5pp - 
   24*G2xx*G4px*pow(G5pp,2) + 144*pow(G4px,3)*G5ppp - 48*G2xx*G4px*G4x*G5ppp + 24*G2xx*G4px*G5p*G5ppp - 72*G2x*pow(G4px,2)*G5ppx + 144*G3p*pow(G4px,2)*G5ppx - 144*G4pp*pow(G4px,2)*G5ppx - 
   24*G2x*G2xx*G4x*G5ppx + 48*G2xx*G3p*G4x*G5ppx - 48*G2xx*G4pp*G4x*G5ppx + 96*G2px*G4px*G4x*G5ppx - 96*G3pp*G4px*G4x*G5ppx + 12*G2x*G2xx*G5p*G5ppx - 24*G2xx*G3p*G5p*G5ppx + 24*G2xx*G4pp*G5p*G5ppx - 
   48*G2px*G4px*G5p*G5ppx + 48*G3pp*G4px*G5p*G5ppx + 36*G2x*G2xx*G4px*G5px - 72*G2xx*G3p*G4px*G5px + 72*G2xx*G4pp*G4px*G5px - 144*G2x*G4ppx*G4px*G5px + 288*G3p*G4ppx*G4px*G5px - 
   288*G4pp*G4ppx*G4px*G5px - 120*G2px*pow(G4px,2)*G5px - 96*G3pp*pow(G4px,2)*G5px + 432*G4ppp*pow(G4px,2)*G5px + 24*G2pxx*G2x*G4x*G5px - 48*G2pxx*G3p*G4x*G5px + 24*G2xx*G3pp*G4x*G5px - 
   24*G2x*G3ppx*G4x*G5px + 48*G3p*G3ppx*G4x*G5px + 48*G2pxx*G4pp*G4x*G5px - 48*G3ppx*G4pp*G4x*G5px - 48*G2xx*G4ppp*G4x*G5px + 96*G2px*G4ppx*G4x*G5px - 96*G3pp*G4ppx*G4x*G5px - 96*G2ppx*G4px*G4x*G5px + 
   96*G3ppp*G4px*G4x*G5px - 12*G2pxx*G2x*G5p*G5px + 24*G2pxx*G3p*G5p*G5px - 12*G2xx*G3pp*G5p*G5px + 12*G2x*G3ppx*G5p*G5px - 24*G3p*G3ppx*G5p*G5px - 24*G2pxx*G4pp*G5p*G5px + 24*G3ppx*G4pp*G5p*G5px + 
   24*G2xx*G4ppp*G5p*G5px - 48*G2px*G4ppx*G5p*G5px + 48*G3pp*G4ppx*G5p*G5px + 48*G2ppx*G4px*G5p*G5px - 48*G3ppp*G4px*G5p*G5px - 12*G2x*G2xx*G5pp*G5px + 24*G2xx*G3p*G5pp*G5px - 24*G2xx*G4pp*G5pp*G5px + 
   48*G2px*G4px*G5pp*G5px - 48*G3pp*G4px*G5pp*G5px - 12*G2x*G2xx*G4ppx*G5x + 24*G2xx*G3p*G4ppx*G5x - 24*G2xx*G4pp*G4ppx*G5x + 12*G2pxx*G2x*G4px*G5x - 24*G2pxx*G3p*G4px*G5x + 12*G2xx*G3pp*G4px*G5x - 
   12*G2x*G3ppx*G4px*G5x + 24*G3p*G3ppx*G4px*G5x + 24*G2pxx*G4pp*G4px*G5x - 24*G3ppx*G4pp*G4px*G5x - 24*G2xx*G4ppp*G4px*G5x + 48*G2px*G4ppx*G4px*G5x - 48*G3pp*G4ppx*G4px*G5x - 
   24*G2ppx*pow(G4px,2)*G5x + 24*G3ppp*pow(G4px,2)*G5x - 16*G2px*G2pxx*G4x*G5x + 16*G2ppx*G2xx*G4x*G5x + 16*G2pxx*G3pp*G4x*G5x - 16*G2xx*G3ppp*G4x*G5x + 16*G2px*G3ppx*G4x*G5x - 
   16*G3pp*G3ppx*G4x*G5x + 8*G2px*G2pxx*G5p*G5x - 8*G2ppx*G2xx*G5p*G5x - 8*G2pxx*G3pp*G5p*G5x + 8*G2xx*G3ppp*G5p*G5x - 8*G2px*G3ppx*G5p*G5x + 8*G3pp*G3ppx*G5p*G5x + 1440*G3xx*G4*G4px*G4pxx*pow(H,2) - 
   1728*G4*G4px*pow(G4pxx,2)*pow(H,2) + 1440*G4*pow(G4px,2)*G4pxxx*pow(H,2) - 576*G3xx*G4*G4ppxx*G4x*pow(H,2) + 864*G3xx*pow(G4px,2)*G4x*pow(H,2) + 576*G3pxx*G4*G4pxx*G4x*pow(H,2) - 
   288*G3xx*G4p*G4pxx*G4x*pow(H,2) + 10368*pow(G4px,2)*G4pxx*G4x*pow(H,2) - 576*G4p*pow(G4pxx,2)*G4x*pow(H,2) - 1152*G4*G4ppx*G4pxxx*G4x*pow(H,2) - 1728*G4p*G4px*G4pxxx*G4x*pow(H,2) - 
   1872*G3xx*G4ppx*pow(G4x,2)*pow(H,2) + 2736*G3pxx*G4px*pow(G4x,2)*pow(H,2) - 6048*G4ppxx*G4px*pow(G4x,2)*pow(H,2) + 576*G2xx*G4pxx*pow(G4x,2)*pow(H,2) + 
   4320*G4ppx*G4pxx*pow(G4x,2)*pow(H,2) - 576*G2x*G4pxxx*pow(G4x,2)*pow(H,2) + 1152*G3p*G4pxxx*pow(G4x,2)*pow(H,2) - 576*G4pp*G4pxxx*pow(G4x,2)*pow(H,2) + 
   1008*G3xx*G4*G4ppx*G4xx*pow(H,2) - 1872*G3pxx*G4*G4px*G4xx*pow(H,2) - 1152*G3xx*G4p*G4px*G4xx*pow(H,2) + 2592*G4*G4ppxx*G4px*G4xx*pow(H,2) - 12096*pow(G4px,3)*G4xx*pow(H,2) - 
   288*G2xx*G4*G4pxx*G4xx*pow(H,2) + 1440*G4*G4ppx*G4pxx*G4xx*pow(H,2) + 10944*G4p*G4px*G4pxx*G4xx*pow(H,2) + 288*G2x*G4*G4pxxx*G4xx*pow(H,2) - 576*G3p*G4*G4pxxx*G4xx*pow(H,2) + 
   576*G4*G4pp*G4pxxx*G4xx*pow(H,2) + 720*G3pxx*G4p*G4x*G4xx*pow(H,2) + 144*G3xx*G4pp*G4x*G4xx*pow(H,2) - 288*G4p*G4ppxx*G4x*G4xx*pow(H,2) - 576*G2xx*G4px*G4x*G4xx*pow(H,2) + 
   6912*G4ppx*G4px*G4x*G4xx*pow(H,2) + 2592*G2x*G4pxx*G4x*G4xx*pow(H,2) - 5184*G3p*G4pxx*G4x*G4xx*pow(H,2) + 864*G4pp*G4pxx*G4x*G4xx*pow(H,2) - 864*G2pxx*pow(G4x,2)*G4xx*pow(H,2) + 
   864*G3ppx*pow(G4x,2)*G4xx*pow(H,2) + 288*G2pxx*G4*pow(G4xx,2)*pow(H,2) + 288*G3ppx*G4*pow(G4xx,2)*pow(H,2) + 288*G2xx*G4p*pow(G4xx,2)*pow(H,2) - 
   1152*G4*G4pppx*pow(G4xx,2)*pow(H,2) - 2304*G4p*G4ppx*pow(G4xx,2)*pow(H,2) - 1440*G2x*G4px*pow(G4xx,2)*pow(H,2) + 2880*G3p*G4px*pow(G4xx,2)*pow(H,2) + 
   6912*G4pp*G4px*pow(G4xx,2)*pow(H,2) - 864*G2px*G4x*pow(G4xx,2)*pow(H,2) + 2304*G3pp*G4x*pow(G4xx,2)*pow(H,2) - 1152*G4ppp*G4x*pow(G4xx,2)*pow(H,2) - 576*G2p*pow(G4xx,3)*pow(H,2) - 
   2016*G4*G4ppx*G4px*G4xxx*pow(H,2) - 288*G4p*pow(G4px,2)*G4xxx*pow(H,2) - 288*G2x*G4*G4pxx*G4xxx*pow(H,2) + 576*G3p*G4*G4pxx*G4xxx*pow(H,2) - 576*G4*G4pp*G4pxx*G4xxx*pow(H,2) - 
   576*G3ppx*G4*G4x*G4xxx*pow(H,2) + 1152*G4*G4pppx*G4x*G4xxx*pow(H,2) + 864*G4p*G4ppx*G4x*G4xxx*pow(H,2) + 864*G4pp*G4px*G4x*G4xxx*pow(H,2) + 288*G2px*pow(G4x,2)*G4xxx*pow(H,2) - 
   864*G3pp*pow(G4x,2)*G4xxx*pow(H,2) + 576*G4ppp*pow(G4x,2)*G4xxx*pow(H,2) + 288*G3pp*G4*G4xx*G4xxx*pow(H,2) + 288*G2x*G4p*G4xx*G4xxx*pow(H,2) - 576*G3p*G4p*G4xx*G4xxx*pow(H,2) + 
   576*G4p*G4pp*G4xx*G4xxx*pow(H,2) - 576*G4*G4ppp*G4xx*G4xxx*pow(H,2) + 288*G2p*G4x*G4xx*G4xxx*pow(H,2) + 288*G3xx*G4*G4ppxx*G5p*pow(H,2) - 1368*G3xx*pow(G4px,2)*G5p*pow(H,2) - 
   288*G3pxx*G4*G4pxx*G5p*pow(H,2) + 144*G3xx*G4p*G4pxx*G5p*pow(H,2) - 7200*pow(G4px,2)*G4pxx*G5p*pow(H,2) + 288*G4p*pow(G4pxx,2)*G5p*pow(H,2) + 576*G4*G4ppx*G4pxxx*G5p*pow(H,2) + 
   864*G4p*G4px*G4pxxx*G5p*pow(H,2) + 2376*G3xx*G4ppx*G4x*G5p*pow(H,2) - 3672*G3pxx*G4px*G4x*G5p*pow(H,2) + 7344*G4ppxx*G4px*G4x*G5p*pow(H,2) - 720*G2xx*G4pxx*G4x*G5p*pow(H,2) - 
   3600*G4ppx*G4pxx*G4x*G5p*pow(H,2) + 720*G2x*G4pxxx*G4x*G5p*pow(H,2) - 1440*G3p*G4pxxx*G4x*G5p*pow(H,2) + 864*G4pp*G4pxxx*G4x*G5p*pow(H,2) - 360*G3pxx*G4p*G4xx*G5p*pow(H,2) - 
   72*G3xx*G4pp*G4xx*G5p*pow(H,2) + 144*G4p*G4ppxx*G4xx*G5p*pow(H,2) + 576*G2xx*G4px*G4xx*G5p*pow(H,2) - 5760*G4ppx*G4px*G4xx*G5p*pow(H,2) - 1872*G2x*G4pxx*G4xx*G5p*pow(H,2) + 
   3744*G3p*G4pxx*G4xx*G5p*pow(H,2) - 1584*G4pp*G4pxx*G4xx*G5p*pow(H,2) + 1152*G2pxx*G4x*G4xx*G5p*pow(H,2) - 576*G3ppx*G4x*G4xx*G5p*pow(H,2) - 1152*G4pppx*G4x*G4xx*G5p*pow(H,2) + 
   432*G2px*pow(G4xx,2)*G5p*pow(H,2) - 2016*G3pp*pow(G4xx,2)*G5p*pow(H,2) + 2304*G4ppp*pow(G4xx,2)*G5p*pow(H,2) + 288*G3ppx*G4*G4xxx*G5p*pow(H,2) - 576*G4*G4pppx*G4xxx*G5p*pow(H,2) - 
   432*G4p*G4ppx*G4xxx*G5p*pow(H,2) + 144*G2x*G4px*G4xxx*G5p*pow(H,2) - 288*G3p*G4px*G4xxx*G5p*pow(H,2) - 144*G4pp*G4px*G4xxx*G5p*pow(H,2) - 288*G2px*G4x*G4xxx*G5p*pow(H,2) + 
   1008*G3pp*G4x*G4xxx*G5p*pow(H,2) - 864*G4ppp*G4x*G4xxx*G5p*pow(H,2) - 144*G2p*G4xx*G4xxx*G5p*pow(H,2) - 720*G3xx*G4ppx*pow(G5p,2)*pow(H,2) + 1152*G3pxx*G4px*pow(G5p,2)*pow(H,2) - 
   2160*G4ppxx*G4px*pow(G5p,2)*pow(H,2) + 216*G2xx*G4pxx*pow(G5p,2)*pow(H,2) + 720*G4ppx*G4pxx*pow(G5p,2)*pow(H,2) - 216*G2x*G4pxxx*pow(G5p,2)*pow(H,2) + 
   432*G3p*G4pxxx*pow(G5p,2)*pow(H,2) - 288*G4pp*G4pxxx*pow(G5p,2)*pow(H,2) - 360*G2pxx*G4xx*pow(G5p,2)*pow(H,2) + 72*G3ppx*G4xx*pow(G5p,2)*pow(H,2) + 
   576*G4pppx*G4xx*pow(G5p,2)*pow(H,2) + 72*G2px*G4xxx*pow(G5p,2)*pow(H,2) - 288*G3pp*G4xxx*pow(G5p,2)*pow(H,2) + 288*G4ppp*G4xxx*pow(G5p,2)*pow(H,2) - 288*G3xx*G4*G4pxx*G5pp*pow(H,2) + 
   576*G4*pow(G4pxx,2)*G5pp*pow(H,2) - 288*G4*G4px*G4pxxx*G5pp*pow(H,2) + 288*G3xx*G4px*G4x*G5pp*pow(H,2) - 1152*G4px*G4pxx*G4x*G5pp*pow(H,2) + 288*G4p*G4pxxx*G4x*G5pp*pow(H,2) - 
   144*G3pxx*pow(G4x,2)*G5pp*pow(H,2) + 288*G4ppxx*pow(G4x,2)*G5pp*pow(H,2) + 288*G3pxx*G4*G4xx*G5pp*pow(H,2) - 72*G3xx*G4p*G4xx*G5pp*pow(H,2) - 576*G4*G4ppxx*G4xx*G5pp*pow(H,2) + 
   3888*pow(G4px,2)*G4xx*G5pp*pow(H,2) - 1584*G4p*G4pxx*G4xx*G5pp*pow(H,2) + 144*G2xx*G4x*G4xx*G5pp*pow(H,2) + 1152*G4ppx*G4x*G4xx*G5pp*pow(H,2) + 432*G2x*pow(G4xx,2)*G5pp*pow(H,2) - 
   864*G3p*pow(G4xx,2)*G5pp*pow(H,2) - 576*G4pp*pow(G4xx,2)*G5pp*pow(H,2) + 288*G4*G4ppx*G4xxx*G5pp*pow(H,2) + 864*G4p*G4px*G4xxx*G5pp*pow(H,2) - 144*G2x*G4x*G4xxx*G5pp*pow(H,2) + 
   288*G3p*G4x*G4xxx*G5pp*pow(H,2) - 576*G4pp*G4x*G4xxx*G5pp*pow(H,2) + 468*G3xx*G4px*G5p*G5pp*pow(H,2) + 216*G4px*G4pxx*G5p*G5pp*pow(H,2) - 144*G4p*G4pxxx*G5p*G5pp*pow(H,2) + 
   288*G3pxx*G4x*G5p*G5pp*pow(H,2) - 576*G4ppxx*G4x*G5p*G5pp*pow(H,2) - 216*G2xx*G4xx*G5p*G5pp*pow(H,2) + 576*G4ppx*G4xx*G5p*G5pp*pow(H,2) + 144*G4pp*G4xxx*G5p*G5pp*pow(H,2) - 
   108*G3pxx*pow(G5p,2)*G5pp*pow(H,2) + 216*G4ppxx*pow(G5p,2)*G5pp*pow(H,2) - 144*G3xx*G4x*pow(G5pp,2)*pow(H,2) + 288*G4pxx*G4x*pow(G5pp,2)*pow(H,2) - 
   144*G4p*G4xxx*pow(G5pp,2)*pow(H,2) + 144*G3xx*pow(G4x,2)*G5ppp*pow(H,2) - 288*G4pxx*pow(G4x,2)*G5ppp*pow(H,2) - 288*G3xx*G4*G4xx*G5ppp*pow(H,2) + 576*G4*G4pxx*G4xx*G5ppp*pow(H,2) + 
   1728*G4p*pow(G4xx,2)*G5ppp*pow(H,2) + 288*G4*G4px*G4xxx*G5ppp*pow(H,2) - 288*G4p*G4x*G4xxx*G5ppp*pow(H,2) - 288*G3xx*G4x*G5p*G5ppp*pow(H,2) + 576*G4pxx*G4x*G5p*G5ppp*pow(H,2) - 
   1728*G4px*G4xx*G5p*G5ppp*pow(H,2) + 144*G4p*G4xxx*G5p*G5ppp*pow(H,2) + 108*G3xx*pow(G5p,2)*G5ppp*pow(H,2) - 216*G4pxx*pow(G5p,2)*G5ppp*pow(H,2) + 288*G3xx*G4*G4x*G5pppx*pow(H,2) - 
   576*G4*G4pxx*G4x*G5pppx*pow(H,2) + 288*G4px*pow(G4x,2)*G5pppx*pow(H,2) + 576*G4*G4px*G4xx*G5pppx*pow(H,2) - 576*G4p*G4x*G4xx*G5pppx*pow(H,2) - 144*G3xx*G4*G5p*G5pppx*pow(H,2) + 
   288*G4*G4pxx*G5p*G5pppx*pow(H,2) + 288*G4p*G4xx*G5p*G5pppx*pow(H,2) - 72*G4px*pow(G5p,2)*G5pppx*pow(H,2) - 792*G3xx*G4*G4px*G5ppx*pow(H,2) + 432*G4*G4px*G4pxx*G5ppx*pow(H,2) - 
   288*G3pxx*G4*G4x*G5ppx*pow(H,2) + 216*G3xx*G4p*G4x*G5ppx*pow(H,2) + 576*G4*G4ppxx*G4x*G5ppx*pow(H,2) - 5424*pow(G4px,2)*G4x*G5ppx*pow(H,2) + 720*G4p*G4pxx*G4x*G5ppx*pow(H,2) - 
   304*G2xx*pow(G4x,2)*G5ppx*pow(H,2) - 576*G4ppx*pow(G4x,2)*G5ppx*pow(H,2) + 192*G2xx*G4*G4xx*G5ppx*pow(H,2) - 1728*G4*G4ppx*G4xx*G5ppx*pow(H,2) - 3456*G4p*G4px*G4xx*G5ppx*pow(H,2) - 
   1200*G2x*G4x*G4xx*G5ppx*pow(H,2) + 2400*G3p*G4x*G4xx*G5ppx*pow(H,2) - 288*G4pp*G4x*G4xx*G5ppx*pow(H,2) + 144*G2x*G4*G4xxx*G5ppx*pow(H,2) - 288*G3p*G4*G4xxx*G5ppx*pow(H,2) + 
   288*G4*G4pp*G4xxx*G5ppx*pow(H,2) + 144*G3pxx*G4*G5p*G5ppx*pow(H,2) - 108*G3xx*G4p*G5p*G5ppx*pow(H,2) - 288*G4*G4ppxx*G5p*G5ppx*pow(H,2) + 4440*pow(G4px,2)*G5p*G5ppx*pow(H,2) - 
   360*G4p*G4pxx*G5p*G5ppx*pow(H,2) + 400*G2xx*G4x*G5p*G5ppx*pow(H,2) - 288*G4ppx*G4x*G5p*G5ppx*pow(H,2) + 888*G2x*G4xx*G5p*G5ppx*pow(H,2) - 1776*G3p*G4xx*G5p*G5ppx*pow(H,2) + 
   720*G4pp*G4xx*G5p*G5ppx*pow(H,2) - 124*G2xx*pow(G5p,2)*G5ppx*pow(H,2) + 288*G4ppx*pow(G5p,2)*G5ppx*pow(H,2) + 144*G3xx*G4*G5pp*G5ppx*pow(H,2) - 288*G4*G4pxx*G5pp*G5ppx*pow(H,2) + 
   144*G4px*G4x*G5pp*G5ppx*pow(H,2) + 864*G4p*G4xx*G5pp*G5ppx*pow(H,2) - 504*G4px*G5p*G5pp*G5ppx*pow(H,2) + 288*G4*G4px*pow(G5ppx,2)*pow(H,2) - 288*G4p*G4x*pow(G5ppx,2)*pow(H,2) + 
   144*G4p*G5p*pow(G5ppx,2)*pow(H,2) - 672*G4*pow(G4px,2)*G5ppxx*pow(H,2) - 32*G2xx*G4*G4x*G5ppxx*pow(H,2) + 576*G4*G4ppx*G4x*G5ppxx*pow(H,2) + 768*G4p*G4px*G4x*G5ppxx*pow(H,2) + 
   304*G2x*pow(G4x,2)*G5ppxx*pow(H,2) - 608*G3p*pow(G4x,2)*G5ppxx*pow(H,2) + 288*G4pp*pow(G4x,2)*G5ppxx*pow(H,2) - 144*G2x*G4*G4xx*G5ppxx*pow(H,2) + 288*G3p*G4*G4xx*G5ppxx*pow(H,2) - 
   288*G4*G4pp*G4xx*G5ppxx*pow(H,2) + 16*G2xx*G4*G5p*G5ppxx*pow(H,2) - 288*G4*G4ppx*G5p*G5ppxx*pow(H,2) - 384*G4p*G4px*G5p*G5ppxx*pow(H,2) - 376*G2x*G4x*G5p*G5ppxx*pow(H,2) + 
   752*G3p*G4x*G5p*G5ppxx*pow(H,2) - 432*G4pp*G4x*G5p*G5ppxx*pow(H,2) + 112*G2x*pow(G5p,2)*G5ppxx*pow(H,2) - 224*G3p*pow(G5p,2)*G5ppxx*pow(H,2) + 144*G4pp*pow(G5p,2)*G5ppxx*pow(H,2) + 
   144*G4*G4px*G5pp*G5ppxx*pow(H,2) - 144*G4p*G4x*G5pp*G5ppxx*pow(H,2) + 72*G4p*G5p*G5pp*G5ppxx*pow(H,2) - 576*G3xx*G4*G4ppx*G5px*pow(H,2) + 1008*G3pxx*G4*G4px*G5px*pow(H,2) + 
   144*G3xx*G4p*G4px*G5px*pow(H,2) - 1440*G4*G4ppxx*G4px*G5px*pow(H,2) + 8328*pow(G4px,3)*G5px*pow(H,2) + 192*G2xx*G4*G4pxx*G5px*pow(H,2) - 576*G4*G4ppx*G4pxx*G5px*pow(H,2) - 
   5184*G4p*G4px*G4pxx*G5px*pow(H,2) - 144*G2x*G4*G4pxxx*G5px*pow(H,2) + 288*G3p*G4*G4pxxx*G5px*pow(H,2) - 288*G4*G4pp*G4pxxx*G5px*pow(H,2) + 36*G2x*G3xx*G4x*G5px*pow(H,2) - 
   72*G3p*G3xx*G4x*G5px*pow(H,2) - 432*G3pxx*G4p*G4x*G5px*pow(H,2) + 72*G3xx*G4pp*G4x*G5px*pow(H,2) + 288*G4p*G4ppxx*G4x*G5px*pow(H,2) + 288*G2xx*G4px*G4x*G5px*pow(H,2) - 
   4272*G4ppx*G4px*G4x*G5px*pow(H,2) - 1560*G2x*G4pxx*G4x*G5px*pow(H,2) + 3120*G3p*G4pxx*G4x*G5px*pow(H,2) - 1008*G4pp*G4pxx*G4x*G5px*pow(H,2) + 448*G2pxx*pow(G4x,2)*G5px*pow(H,2) - 
   592*G3ppx*pow(G4x,2)*G5px*pow(H,2) + 288*G4pppx*pow(G4x,2)*G5px*pow(H,2) - 336*G2pxx*G4*G4xx*G5px*pow(H,2) - 240*G3ppx*G4*G4xx*G5px*pow(H,2) - 216*G2xx*G4p*G4xx*G5px*pow(H,2) + 
   1152*G4*G4pppx*G4xx*G5px*pow(H,2) + 1440*G4p*G4ppx*G4xx*G5px*pow(H,2) + 2520*G2x*G4px*G4xx*G5px*pow(H,2) - 5040*G3p*G4px*G4xx*G5px*pow(H,2) - 6336*G4pp*G4px*G4xx*G5px*pow(H,2) + 
   1008*G2px*G4x*G4xx*G5px*pow(H,2) - 2352*G3pp*G4x*G4xx*G5px*pow(H,2) + 864*G4ppp*G4x*G4xx*G5px*pow(H,2) + 960*G2p*pow(G4xx,2)*G5px*pow(H,2) - 144*G3pp*G4*G4xxx*G5px*pow(H,2) - 
   72*G2x*G4p*G4xxx*G5px*pow(H,2) + 144*G3p*G4p*G4xxx*G5px*pow(H,2) - 144*G4p*G4pp*G4xxx*G5px*pow(H,2) + 288*G4*G4ppp*G4xxx*G5px*pow(H,2) - 144*G2p*G4x*G4xxx*G5px*pow(H,2) + 
   18*G2x*G3xx*G5p*G5px*pow(H,2) - 36*G3p*G3xx*G5p*G5px*pow(H,2) + 216*G3pxx*G4p*G5p*G5px*pow(H,2) + 36*G3xx*G4pp*G5p*G5px*pow(H,2) - 144*G4p*G4ppxx*G5p*G5px*pow(H,2) - 
   468*G2xx*G4px*G5p*G5px*pow(H,2) + 3720*G4ppx*G4px*G5p*G5px*pow(H,2) + 996*G2x*G4pxx*G5p*G5px*pow(H,2) - 1992*G3p*G4pxx*G5p*G5px*pow(H,2) + 936*G4pp*G4pxx*G5p*G5px*pow(H,2) - 
   616*G2pxx*G4x*G5p*G5px*pow(H,2) + 472*G3ppx*G4x*G5p*G5px*pow(H,2) + 288*G4pppx*G4x*G5p*G5px*pow(H,2) - 408*G2px*G4xx*G5p*G5px*pow(H,2) + 1944*G3pp*G4xx*G5p*G5px*pow(H,2) - 
   2160*G4ppp*G4xx*G5p*G5px*pow(H,2) + 72*G2p*G4xxx*G5p*G5px*pow(H,2) + 196*G2pxx*pow(G5p,2)*G5px*pow(H,2) - 88*G3ppx*pow(G5p,2)*G5px*pow(H,2) - 216*G4pppx*pow(G5p,2)*G5px*pow(H,2) - 
   144*G3pxx*G4*G5pp*G5px*pow(H,2) + 144*G3xx*G4p*G5pp*G5px*pow(H,2) + 288*G4*G4ppxx*G5pp*G5px*pow(H,2) - 2616*pow(G4px,2)*G5pp*G5px*pow(H,2) + 576*G4p*G4pxx*G5pp*G5px*pow(H,2) - 
   24*G2xx*G4x*G5pp*G5px*pow(H,2) - 432*G4ppx*G4x*G5pp*G5px*pow(H,2) - 576*G2x*G4xx*G5pp*G5px*pow(H,2) + 1152*G3p*G4xx*G5pp*G5px*pow(H,2) + 432*G4pp*G4xx*G5pp*G5px*pow(H,2) + 
   120*G2xx*G5p*G5pp*G5px*pow(H,2) - 360*G4ppx*G5p*G5pp*G5px*pow(H,2) + 72*G4px*pow(G5pp,2)*G5px*pow(H,2) + 144*G3xx*G4*G5ppp*G5px*pow(H,2) - 288*G4*G4pxx*G5ppp*G5px*pow(H,2) + 
   144*G4px*G4x*G5ppp*G5px*pow(H,2) - 1728*G4p*G4xx*G5ppp*G5px*pow(H,2) + 792*G4px*G5p*G5ppp*G5px*pow(H,2) - 288*G4*G4px*G5pppx*G5px*pow(H,2) + 288*G4p*G4x*G5pppx*G5px*pow(H,2) - 
   144*G4p*G5p*G5pppx*G5px*pow(H,2) - 120*G2xx*G4*G5ppx*G5px*pow(H,2) + 864*G4*G4ppx*G5ppx*G5px*pow(H,2) + 2016*G4p*G4px*G5ppx*G5px*pow(H,2) + 696*G2x*G4x*G5ppx*G5px*pow(H,2) - 
   1392*G3p*G4x*G5ppx*G5px*pow(H,2) + 288*G4pp*G4x*G5ppx*G5px*pow(H,2) - 492*G2x*G5p*G5ppx*G5px*pow(H,2) + 984*G3p*G5p*G5ppx*G5px*pow(H,2) - 432*G4pp*G5p*G5ppx*G5px*pow(H,2) - 
   432*G4p*G5pp*G5ppx*G5px*pow(H,2) + 72*G2x*G4*G5ppxx*G5px*pow(H,2) - 144*G3p*G4*G5ppxx*G5px*pow(H,2) + 144*G4*G4pp*G5ppxx*G5px*pow(H,2) + 96*G2pxx*G4*pow(G5px,2)*pow(H,2) + 
   48*G3ppx*G4*pow(G5px,2)*pow(H,2) + 24*G2xx*G4p*pow(G5px,2)*pow(H,2) - 288*G4*G4pppx*pow(G5px,2)*pow(H,2) - 144*G4p*G4ppx*pow(G5px,2)*pow(H,2) - 948*G2x*G4px*pow(G5px,2)*pow(H,2) + 
   1896*G3p*G4px*pow(G5px,2)*pow(H,2) + 1368*G4pp*G4px*pow(G5px,2)*pow(H,2) - 288*G2px*G4x*pow(G5px,2)*pow(H,2) + 600*G3pp*G4x*pow(G5px,2)*pow(H,2) - 
   144*G4ppp*G4x*pow(G5px,2)*pow(H,2) - 528*G2p*G4xx*pow(G5px,2)*pow(H,2) + 96*G2px*G5p*pow(G5px,2)*pow(H,2) - 468*G3pp*G5p*pow(G5px,2)*pow(H,2) + 504*G4ppp*G5p*pow(G5px,2)*pow(H,2) + 
   180*G2x*G5pp*pow(G5px,2)*pow(H,2) - 360*G3p*G5pp*pow(G5px,2)*pow(H,2) - 72*G4pp*G5pp*pow(G5px,2)*pow(H,2) + 432*G4p*G5ppp*pow(G5px,2)*pow(H,2) + 96*G2p*pow(G5px,3)*pow(H,2) - 
   36*G2x*G3xx*G4*G5pxx*pow(H,2) + 72*G3p*G3xx*G4*G5pxx*pow(H,2) - 72*G3xx*G4*G4pp*G5pxx*pow(H,2) + 120*G2xx*G4*G4px*G5pxx*pow(H,2) + 624*G4*G4ppx*G4px*G5pxx*pow(H,2) - 
   696*G4p*pow(G4px,2)*G5pxx*pow(H,2) + 216*G2x*G4*G4pxx*G5pxx*pow(H,2) - 432*G3p*G4*G4pxx*G5pxx*pow(H,2) + 432*G4*G4pp*G4pxx*G5pxx*pow(H,2) + 32*G2pxx*G4*G4x*G5pxx*pow(H,2) + 
   256*G3ppx*G4*G4x*G5pxx*pow(H,2) - 24*G2xx*G4p*G4x*G5pxx*pow(H,2) - 576*G4*G4pppx*G4x*G5pxx*pow(H,2) - 48*G4p*G4ppx*G4x*G5pxx*pow(H,2) - 456*G2x*G4px*G4x*G5pxx*pow(H,2) + 
   912*G3p*G4px*G4x*G5pxx*pow(H,2) - 720*G4pp*G4px*G4x*G5pxx*pow(H,2) - 176*G2px*pow(G4x,2)*G5pxx*pow(H,2) + 480*G3pp*pow(G4x,2)*G5pxx*pow(H,2) - 288*G4ppp*pow(G4x,2)*G5pxx*pow(H,2) - 
   96*G2px*G4*G4xx*G5pxx*pow(H,2) - 48*G3pp*G4*G4xx*G5pxx*pow(H,2) - 288*G2x*G4p*G4xx*G5pxx*pow(H,2) + 576*G3p*G4p*G4xx*G5pxx*pow(H,2) - 576*G4p*G4pp*G4xx*G5pxx*pow(H,2) + 
   288*G4*G4ppp*G4xx*G5pxx*pow(H,2) - 240*G2p*G4x*G4xx*G5pxx*pow(H,2) - 16*G2pxx*G4*G5p*G5pxx*pow(H,2) - 128*G3ppx*G4*G5p*G5pxx*pow(H,2) + 12*G2xx*G4p*G5p*G5pxx*pow(H,2) + 
   288*G4*G4pppx*G5p*G5pxx*pow(H,2) + 24*G4p*G4ppx*G5p*G5pxx*pow(H,2) + 228*G2x*G4px*G5p*G5pxx*pow(H,2) - 456*G3p*G4px*G5p*G5pxx*pow(H,2) + 360*G4pp*G4px*G5p*G5pxx*pow(H,2) + 
   128*G2px*G4x*G5p*G5pxx*pow(H,2) - 504*G3pp*G4x*G5p*G5pxx*pow(H,2) + 432*G4ppp*G4x*G5p*G5pxx*pow(H,2) + 120*G2p*G4xx*G5p*G5pxx*pow(H,2) - 20*G2px*pow(G5p,2)*G5pxx*pow(H,2) + 
   132*G3pp*pow(G5p,2)*G5pxx*pow(H,2) - 144*G4ppp*pow(G5p,2)*G5pxx*pow(H,2) - 24*G2xx*G4*G5pp*G5pxx*pow(H,2) - 144*G4*G4ppx*G5pp*G5pxx*pow(H,2) - 240*G4p*G4px*G5pp*G5pxx*pow(H,2) + 
   96*G2x*G4x*G5pp*G5pxx*pow(H,2) - 192*G3p*G4x*G5pp*G5pxx*pow(H,2) + 288*G4pp*G4x*G5pp*G5pxx*pow(H,2) - 12*G2x*G5p*G5pp*G5pxx*pow(H,2) + 24*G3p*G5p*G5pp*G5pxx*pow(H,2) - 
   72*G4pp*G5p*G5pp*G5pxx*pow(H,2) + 72*G4p*pow(G5pp,2)*G5pxx*pow(H,2) - 144*G4*G4px*G5ppp*G5pxx*pow(H,2) + 144*G4p*G4x*G5ppp*G5pxx*pow(H,2) - 72*G4p*G5p*G5ppp*G5pxx*pow(H,2) - 
   72*G2x*G4*G5ppx*G5pxx*pow(H,2) + 144*G3p*G4*G5ppx*G5pxx*pow(H,2) - 144*G4*G4pp*G5ppx*G5pxx*pow(H,2) + 48*G2px*G4*G5px*G5pxx*pow(H,2) + 24*G3pp*G4*G5px*G5pxx*pow(H,2) + 
   108*G2x*G4p*G5px*G5pxx*pow(H,2) - 216*G3p*G4p*G5px*G5pxx*pow(H,2) + 216*G4p*G4pp*G5px*G5pxx*pow(H,2) - 144*G4*G4ppp*G5px*G5pxx*pow(H,2) + 120*G2p*G4x*G5px*G5pxx*pow(H,2) - 
   60*G2p*G5p*G5px*G5pxx*pow(H,2) - 24*G2x*G4*G4px*G5pxxx*pow(H,2) + 48*G3p*G4*G4px*G5pxxx*pow(H,2) - 48*G4*G4pp*G4px*G5pxxx*pow(H,2) + 32*G2px*G4*G4x*G5pxxx*pow(H,2) - 
   32*G3pp*G4*G4x*G5pxxx*pow(H,2) + 24*G2x*G4p*G4x*G5pxxx*pow(H,2) - 48*G3p*G4p*G4x*G5pxxx*pow(H,2) + 48*G4p*G4pp*G4x*G5pxxx*pow(H,2) + 16*G2p*pow(G4x,2)*G5pxxx*pow(H,2) - 
   16*G2px*G4*G5p*G5pxxx*pow(H,2) + 16*G3pp*G4*G5p*G5pxxx*pow(H,2) - 12*G2x*G4p*G5p*G5pxxx*pow(H,2) + 24*G3p*G4p*G5p*G5pxxx*pow(H,2) - 24*G4p*G4pp*G5p*G5pxxx*pow(H,2) - 
   16*G2p*G4x*G5p*G5pxxx*pow(H,2) + 4*G2p*pow(G5p,2)*G5pxxx*pow(H,2) - 144*G3ppx*G3xx*G4*G5x*pow(H,2) + 288*G3xx*G4*G4pppx*G5x*pow(H,2) - 288*G3pxx*G4*G4ppx*G5x*pow(H,2) + 
   216*G3xx*G4p*G4ppx*G5x*pow(H,2) - 96*G2xx*G4*G4ppxx*G5x*pow(H,2) + 576*G4*G4ppx*G4ppxx*G5x*pow(H,2) - 36*G2x*G3xx*G4px*G5x*pow(H,2) + 72*G3p*G3xx*G4px*G5x*pow(H,2) - 
   432*G3pxx*G4p*G4px*G5x*pow(H,2) + 144*G3xx*G4pp*G4px*G5x*pow(H,2) + 576*G4p*G4ppxx*G4px*G5x*pow(H,2) + 264*G2xx*pow(G4px,2)*G5x*pow(H,2) - 2664*G4ppx*pow(G4px,2)*G5x*pow(H,2) + 
   96*G2pxx*G4*G4pxx*G5x*pow(H,2) + 192*G3ppx*G4*G4pxx*G5x*pow(H,2) - 48*G2xx*G4p*G4pxx*G5x*pow(H,2) - 576*G4*G4pppx*G4pxx*G5x*pow(H,2) + 432*G4p*G4ppx*G4pxx*G5x*pow(H,2) - 
   1032*G2x*G4px*G4pxx*G5x*pow(H,2) + 2064*G3p*G4px*G4pxx*G5x*pow(H,2) - 1152*G4pp*G4px*G4pxx*G5x*pow(H,2) + 96*G2px*G4*G4pxxx*G5x*pow(H,2) - 96*G3pp*G4*G4pxxx*G5x*pow(H,2) + 
   72*G2x*G4p*G4pxxx*G5x*pow(H,2) - 144*G3p*G4p*G4pxxx*G5x*pow(H,2) + 144*G4p*G4pp*G4pxxx*G5x*pow(H,2) - 324*G2x*G3pxx*G4x*G5x*pow(H,2) + 648*G3p*G3pxx*G4x*G5x*pow(H,2) + 
   144*G2px*G3xx*G4x*G5x*pow(H,2) - 468*G3pp*G3xx*G4x*G5x*pow(H,2) - 360*G3pxx*G4pp*G4x*G5x*pow(H,2) + 360*G3xx*G4ppp*G4x*G5x*pow(H,2) - 672*G2xx*G4ppx*G4x*G5x*pow(H,2) - 
   288*pow(G4ppx,2)*G4x*G5x*pow(H,2) + 744*G2x*G4ppxx*G4x*G5x*pow(H,2) - 1488*G3p*G4ppxx*G4x*G5x*pow(H,2) + 720*G4pp*G4ppxx*G4x*G5x*pow(H,2) + 1032*G2pxx*G4px*G4x*G5x*pow(H,2) - 
   1176*G3ppx*G4px*G4x*G5x*pow(H,2) + 288*G4pppx*G4px*G4x*G5x*pow(H,2) - 288*G2px*G4pxx*G4x*G5x*pow(H,2) + 1032*G3pp*G4pxx*G4x*G5x*pow(H,2) - 720*G4ppp*G4pxx*G4x*G5x*pow(H,2) + 
   96*G2p*G4pxxx*G4x*G5x*pow(H,2) + 72*G2p*G3xx*G4xx*G5x*pow(H,2) + 120*G2pxx*G4p*G4xx*G5x*pow(H,2) + 168*G3ppx*G4p*G4xx*G5x*pow(H,2) + 24*G2xx*G4pp*G4xx*G5x*pow(H,2) - 
   576*G4p*G4pppx*G4xx*G5x*pow(H,2) - 744*G2x*G4ppx*G4xx*G5x*pow(H,2) + 1488*G3p*G4ppx*G4xx*G5x*pow(H,2) + 504*G2px*G4px*G4xx*G5x*pow(H,2) - 1680*G3pp*G4px*G4xx*G5x*pow(H,2) + 
   1152*G4ppp*G4px*G4xx*G5x*pow(H,2) - 336*G2p*G4pxx*G4xx*G5x*pow(H,2) - 96*G2ppx*G4x*G4xx*G5x*pow(H,2) + 96*G3ppp*G4x*G4xx*G5x*pow(H,2) + 96*G2pp*pow(G4xx,2)*G5x*pow(H,2) - 
   96*G2ppx*G4*G4xxx*G5x*pow(H,2) + 96*G3ppp*G4*G4xxx*G5x*pow(H,2) + 72*G3pp*G4p*G4xxx*G5x*pow(H,2) - 72*G2x*G4pp*G4xxx*G5x*pow(H,2) + 144*G3p*G4pp*G4xxx*G5x*pow(H,2) - 
   144*pow(G4pp,2)*G4xxx*G5x*pow(H,2) - 144*G4p*G4ppp*G4xxx*G5x*pow(H,2) - 72*G2p*G4px*G4xxx*G5x*pow(H,2) - 96*G2pp*G4x*G4xxx*G5x*pow(H,2) + 198*G2x*G3pxx*G5p*G5x*pow(H,2) - 
   396*G3p*G3pxx*G5p*G5x*pow(H,2) - 72*G2px*G3xx*G5p*G5x*pow(H,2) + 270*G3pp*G3xx*G5p*G5x*pow(H,2) + 252*G3pxx*G4pp*G5p*G5x*pow(H,2) - 252*G3xx*G4ppp*G5p*G5x*pow(H,2) + 
   420*G2xx*G4ppx*G5p*G5x*pow(H,2) - 144*pow(G4ppx,2)*G5p*G5x*pow(H,2) - 444*G2x*G4ppxx*G5p*G5x*pow(H,2) + 888*G3p*G4ppxx*G5p*G5x*pow(H,2) - 504*G4pp*G4ppxx*G5p*G5x*pow(H,2) - 
   672*G2pxx*G4px*G5p*G5x*pow(H,2) + 600*G3ppx*G4px*G5p*G5x*pow(H,2) + 144*G4pppx*G4px*G5p*G5x*pow(H,2) + 48*G2px*G4pxx*G5p*G5x*pow(H,2) - 492*G3pp*G4pxx*G5p*G5x*pow(H,2) + 
   504*G4ppp*G4pxx*G5p*G5x*pow(H,2) - 48*G2p*G4pxxx*G5p*G5x*pow(H,2) + 144*G2ppx*G4xx*G5p*G5x*pow(H,2) - 144*G3ppp*G4xx*G5p*G5x*pow(H,2) + 48*G2pp*G4xxx*G5p*G5x*pow(H,2) - 
   18*G2x*G3xx*G5pp*G5x*pow(H,2) + 36*G3p*G3xx*G5pp*G5x*pow(H,2) + 72*G3pxx*G4p*G5pp*G5x*pow(H,2) - 108*G3xx*G4pp*G5pp*G5x*pow(H,2) - 144*G4p*G4ppxx*G5pp*G5x*pow(H,2) - 
   36*G2xx*G4px*G5pp*G5x*pow(H,2) + 144*G4ppx*G4px*G5pp*G5x*pow(H,2) + 84*G2x*G4pxx*G5pp*G5x*pow(H,2) - 168*G3p*G4pxx*G5pp*G5x*pow(H,2) + 216*G4pp*G4pxx*G5pp*G5x*pow(H,2) - 
   72*G2pxx*G4x*G5pp*G5x*pow(H,2) + 72*G3ppx*G4x*G5pp*G5x*pow(H,2) - 96*G2px*G4xx*G5pp*G5x*pow(H,2) + 144*G3pp*G4xx*G5pp*G5x*pow(H,2) + 60*G2pxx*G5p*G5pp*G5x*pow(H,2) - 
   60*G3ppx*G5p*G5pp*G5x*pow(H,2) - 12*G2xx*pow(G5pp,2)*G5x*pow(H,2) - 72*G3xx*G4p*G5ppp*G5x*pow(H,2) - 144*pow(G4px,2)*G5ppp*G5x*pow(H,2) + 144*G4p*G4pxx*G5ppp*G5x*pow(H,2) + 
   72*G2xx*G4x*G5ppp*G5x*pow(H,2) + 48*G2x*G4xx*G5ppp*G5x*pow(H,2) - 96*G3p*G4xx*G5ppp*G5x*pow(H,2) - 60*G2xx*G5p*G5ppp*G5x*pow(H,2) + 48*G2xx*G4*G5pppx*G5x*pow(H,2) + 
   144*G4p*G4px*G5pppx*G5x*pow(H,2) - 48*G2x*G4x*G5pppx*G5x*pow(H,2) + 96*G3p*G4x*G5pppx*G5x*pow(H,2) + 24*G2x*G5p*G5pppx*G5x*pow(H,2) - 48*G3p*G5p*G5pppx*G5x*pow(H,2) - 
   48*G2pxx*G4*G5ppx*G5x*pow(H,2) + 48*G3ppx*G4*G5ppx*G5x*pow(H,2) + 36*G2xx*G4p*G5ppx*G5x*pow(H,2) - 432*G4p*G4ppx*G5ppx*G5x*pow(H,2) + 528*G2x*G4px*G5ppx*G5x*pow(H,2) - 
   1056*G3p*G4px*G5ppx*G5x*pow(H,2) + 360*G4pp*G4px*G5ppx*G5x*pow(H,2) + 48*G2px*G4x*G5ppx*G5x*pow(H,2) - 96*G3pp*G4x*G5ppx*G5x*pow(H,2) + 96*G2p*G4xx*G5ppx*G5x*pow(H,2) + 
   24*G2px*G5p*G5ppx*G5x*pow(H,2) - 24*G2x*G5pp*G5ppx*G5x*pow(H,2) + 48*G3p*G5pp*G5ppx*G5x*pow(H,2) - 48*G2px*G4*G5ppxx*G5x*pow(H,2) + 48*G3pp*G4*G5ppxx*G5x*pow(H,2) - 
   36*G2x*G4p*G5ppxx*G5x*pow(H,2) + 72*G3p*G4p*G5ppxx*G5x*pow(H,2) - 72*G4p*G4pp*G5ppxx*G5x*pow(H,2) - 48*G2p*G4x*G5ppxx*G5x*pow(H,2) + 24*G2p*G5p*G5ppxx*G5x*pow(H,2) - 
   36*G2p*G3xx*G5px*G5x*pow(H,2) - 72*G2pxx*G4p*G5px*G5x*pow(H,2) - 72*G3ppx*G4p*G5px*G5x*pow(H,2) + 288*G4p*G4pppx*G5px*G5x*pow(H,2) + 420*G2x*G4ppx*G5px*G5x*pow(H,2) - 
   840*G3p*G4ppx*G5px*G5x*pow(H,2) + 72*G4pp*G4ppx*G5px*G5x*pow(H,2) - 264*G2px*G4px*G5px*G5x*pow(H,2) + 828*G3pp*G4px*G5px*G5x*pow(H,2) - 504*G4ppp*G4px*G5px*G5x*pow(H,2) + 
   168*G2p*G4pxx*G5px*G5x*pow(H,2) - 96*G2pp*G4xx*G5px*G5x*pow(H,2) - 48*G2ppx*G5p*G5px*G5x*pow(H,2) + 48*G3ppp*G5p*G5px*G5x*pow(H,2) + 48*G2px*G5pp*G5px*G5x*pow(H,2) - 
   72*G3pp*G5pp*G5px*G5x*pow(H,2) - 24*G2x*G5ppp*G5px*G5x*pow(H,2) + 48*G3p*G5ppp*G5px*G5x*pow(H,2) - 48*G2p*G5ppx*G5px*G5x*pow(H,2) + 24*G2pp*pow(G5px,2)*G5x*pow(H,2) + 
   6*pow(G2x,2)*G5pxx*G5x*pow(H,2) - 24*G2x*G3p*G5pxx*G5x*pow(H,2) + 24*pow(G3p,2)*G5pxx*G5x*pow(H,2) + 48*G2ppx*G4*G5pxx*G5x*pow(H,2) - 48*G3ppp*G4*G5pxx*G5x*pow(H,2) - 
   24*G2px*G4p*G5pxx*G5x*pow(H,2) - 12*G3pp*G4p*G5pxx*G5x*pow(H,2) + 48*G2x*G4pp*G5pxx*G5x*pow(H,2) - 96*G3p*G4pp*G5pxx*G5x*pow(H,2) + 72*pow(G4pp,2)*G5pxx*G5x*pow(H,2) + 
   72*G4p*G4ppp*G5pxx*G5x*pow(H,2) + 60*G2p*G4px*G5pxx*G5x*pow(H,2) + 48*G2pp*G4x*G5pxx*G5x*pow(H,2) - 24*G2pp*G5p*G5pxx*G5x*pow(H,2) - 30*G2pxx*G2x*pow(G5x,2)*pow(H,2) + 
   12*G2px*G2xx*pow(G5x,2)*pow(H,2) + 60*G2pxx*G3p*pow(G5x,2)*pow(H,2) - 42*G2xx*G3pp*pow(G5x,2)*pow(H,2) + 42*G2x*G3ppx*pow(G5x,2)*pow(H,2) - 84*G3p*G3ppx*pow(G5x,2)*pow(H,2) + 
   12*G2p*G3pxx*pow(G5x,2)*pow(H,2) - 12*G2pp*G3xx*pow(G5x,2)*pow(H,2) - 36*G2pxx*G4pp*pow(G5x,2)*pow(H,2) + 36*G3ppx*G4pp*pow(G5x,2)*pow(H,2) + 36*G2xx*G4ppp*pow(G5x,2)*pow(H,2) - 
   24*G2x*G4pppx*pow(G5x,2)*pow(H,2) + 48*G3p*G4pppx*pow(G5x,2)*pow(H,2) + 24*G2px*G4ppx*pow(G5x,2)*pow(H,2) - 48*G3pp*G4ppx*pow(G5x,2)*pow(H,2) - 24*G2p*G4ppxx*pow(G5x,2)*pow(H,2) + 
   24*G2pp*G4pxx*pow(G5x,2)*pow(H,2) + 36*G2x*G3pxx*G4*G5xx*pow(H,2) - 72*G3p*G3pxx*G4*G5xx*pow(H,2) + 36*G3pp*G3xx*G4*G5xx*pow(H,2) + 36*G2x*G3xx*G4p*G5xx*pow(H,2) - 
   72*G3p*G3xx*G4p*G5xx*pow(H,2) + 72*G3pxx*G4*G4pp*G5xx*pow(H,2) + 72*G3xx*G4p*G4pp*G5xx*pow(H,2) - 72*G3xx*G4*G4ppp*G5xx*pow(H,2) + 72*G2xx*G4*G4ppx*G5xx*pow(H,2) - 
   288*G4*pow(G4ppx,2)*G5xx*pow(H,2) - 72*G2x*G4*G4ppxx*G5xx*pow(H,2) + 144*G3p*G4*G4ppxx*G5xx*pow(H,2) - 144*G4*G4pp*G4ppxx*G5xx*pow(H,2) - 144*G2pxx*G4*G4px*G5xx*pow(H,2) - 
   96*G2xx*G4p*G4px*G5xx*pow(H,2) + 288*G4*G4pppx*G4px*G5xx*pow(H,2) - 144*G4p*G4ppx*G4px*G5xx*pow(H,2) + 528*G2x*pow(G4px,2)*G5xx*pow(H,2) - 1056*G3p*pow(G4px,2)*G5xx*pow(H,2) - 
   744*G4pp*pow(G4px,2)*G5xx*pow(H,2) - 96*G2px*G4*G4pxx*G5xx*pow(H,2) + 24*G3pp*G4*G4pxx*G5xx*pow(H,2) - 216*G2x*G4p*G4pxx*G5xx*pow(H,2) + 432*G3p*G4p*G4pxx*G5xx*pow(H,2) - 
   432*G4p*G4pp*G4pxx*G5xx*pow(H,2) + 144*G4*G4ppp*G4pxx*G5xx*pow(H,2) + 36*G2p*G3xx*G4x*G5xx*pow(H,2) + 48*G2pxx*G4p*G4x*G5xx*pow(H,2) + 96*G3ppx*G4p*G4x*G5xx*pow(H,2) + 
   24*G2xx*G4pp*G4x*G5xx*pow(H,2) - 288*G4p*G4pppx*G4x*G5xx*pow(H,2) - 312*G2x*G4ppx*G4x*G5xx*pow(H,2) + 624*G3p*G4ppx*G4x*G5xx*pow(H,2) + 96*G4pp*G4ppx*G4x*G5xx*pow(H,2) + 
   288*G2px*G4px*G4x*G5xx*pow(H,2) - 624*G3pp*G4px*G4x*G5xx*pow(H,2) + 96*G4ppp*G4px*G4x*G5xx*pow(H,2) - 168*G2p*G4pxx*G4x*G5xx*pow(H,2) + 16*G2ppx*pow(G4x,2)*G5xx*pow(H,2) - 
   16*G3ppp*pow(G4x,2)*G5xx*pow(H,2) + 96*G2ppx*G4*G4xx*G5xx*pow(H,2) - 96*G3ppp*G4*G4xx*G5xx*pow(H,2) + 96*G2px*G4p*G4xx*G5xx*pow(H,2) - 528*G3pp*G4p*G4xx*G5xx*pow(H,2) - 
   144*G2x*G4pp*G4xx*G5xx*pow(H,2) + 288*G3p*G4pp*G4xx*G5xx*pow(H,2) - 288*pow(G4pp,2)*G4xx*G5xx*pow(H,2) + 864*G4p*G4ppp*G4xx*G5xx*pow(H,2) + 528*G2p*G4px*G4xx*G5xx*pow(H,2) + 
   96*G2pp*G4x*G4xx*G5xx*pow(H,2) - 18*G2p*G3xx*G5p*G5xx*pow(H,2) - 24*G2pxx*G4p*G5p*G5xx*pow(H,2) - 48*G3ppx*G4p*G5p*G5xx*pow(H,2) - 12*G2xx*G4pp*G5p*G5xx*pow(H,2) + 
   144*G4p*G4pppx*G5p*G5xx*pow(H,2) + 228*G2x*G4ppx*G5p*G5xx*pow(H,2) - 456*G3p*G4ppx*G5p*G5xx*pow(H,2) + 96*G4pp*G4ppx*G5p*G5xx*pow(H,2) - 96*G2px*G4px*G5p*G5xx*pow(H,2) + 
   480*G3pp*G4px*G5p*G5xx*pow(H,2) - 480*G4ppp*G4px*G5p*G5xx*pow(H,2) + 84*G2p*G4pxx*G5p*G5xx*pow(H,2) + 32*G2ppx*G4x*G5p*G5xx*pow(H,2) - 32*G3ppp*G4x*G5p*G5xx*pow(H,2) - 
   48*G2pp*G4xx*G5p*G5xx*pow(H,2) - 20*G2ppx*pow(G5p,2)*G5xx*pow(H,2) + 20*G3ppp*pow(G5p,2)*G5xx*pow(H,2) + 24*G2pxx*G4*G5pp*G5xx*pow(H,2) - 24*G3ppx*G4*G5pp*G5xx*pow(H,2) + 
   288*G4p*G4ppx*G5pp*G5xx*pow(H,2) - 216*G2x*G4px*G5pp*G5xx*pow(H,2) + 432*G3p*G4px*G5pp*G5xx*pow(H,2) + 48*G4pp*G4px*G5pp*G5xx*pow(H,2) - 72*G2px*G4x*G5pp*G5xx*pow(H,2) + 
   96*G3pp*G4x*G5pp*G5xx*pow(H,2) - 48*G2p*G4xx*G5pp*G5xx*pow(H,2) + 12*G2px*G5p*G5pp*G5xx*pow(H,2) - 24*G3pp*G5p*G5pp*G5xx*pow(H,2) + 12*G2x*pow(G5pp,2)*G5xx*pow(H,2) - 
   24*G3p*pow(G5pp,2)*G5xx*pow(H,2) - 24*G2xx*G4*G5ppp*G5xx*pow(H,2) - 432*G4p*G4px*G5ppp*G5xx*pow(H,2) + 24*G2x*G4x*G5ppp*G5xx*pow(H,2) - 48*G3p*G4x*G5ppp*G5xx*pow(H,2) - 
   12*G2x*G5p*G5ppp*G5xx*pow(H,2) + 24*G3p*G5p*G5ppp*G5xx*pow(H,2) + 48*G2px*G4*G5ppx*G5xx*pow(H,2) - 48*G3pp*G4*G5ppx*G5xx*pow(H,2) + 72*G2x*G4p*G5ppx*G5xx*pow(H,2) - 
   144*G3p*G4p*G5ppx*G5xx*pow(H,2) + 144*G4p*G4pp*G5ppx*G5xx*pow(H,2) + 48*G2p*G4x*G5ppx*G5xx*pow(H,2) - 24*G2p*G5p*G5ppx*G5xx*pow(H,2) - 6*pow(G2x,2)*G5px*G5xx*pow(H,2) + 
   24*G2x*G3p*G5px*G5xx*pow(H,2) - 24*pow(G3p,2)*G5px*G5xx*pow(H,2) - 48*G2ppx*G4*G5px*G5xx*pow(H,2) + 48*G3ppp*G4*G5px*G5xx*pow(H,2) - 24*G2px*G4p*G5px*G5xx*pow(H,2) + 
   240*G3pp*G4p*G5px*G5xx*pow(H,2) + 60*G2x*G4pp*G5px*G5xx*pow(H,2) - 120*G3p*G4pp*G5px*G5xx*pow(H,2) + 144*pow(G4pp,2)*G5px*G5xx*pow(H,2) - 432*G4p*G4ppp*G5px*G5xx*pow(H,2) - 
   288*G2p*G4px*G5px*G5xx*pow(H,2) - 48*G2pp*G4x*G5px*G5xx*pow(H,2) + 24*G2pp*G5p*G5px*G5xx*pow(H,2) + 24*G2p*G5pp*G5px*G5xx*pow(H,2) - 6*G2px*G2x*G5x*G5xx*pow(H,2) + 
   6*G2p*G2xx*G5x*G5xx*pow(H,2) + 12*G2px*G3p*G5x*G5xx*pow(H,2) + 6*G2x*G3pp*G5x*G5xx*pow(H,2) - 12*G3p*G3pp*G5x*G5xx*pow(H,2) + 24*G2ppx*G4p*G5x*G5xx*pow(H,2) - 
   24*G3ppp*G4p*G5x*G5xx*pow(H,2) - 36*G2px*G4pp*G5x*G5xx*pow(H,2) + 48*G3pp*G4pp*G5x*G5xx*pow(H,2) + 12*G2x*G4ppp*G5x*G5xx*pow(H,2) - 24*G3p*G4ppp*G5x*G5xx*pow(H,2) + 
   24*G2p*G4ppx*G5x*G5xx*pow(H,2) - 24*G2pp*G4px*G5x*G5xx*pow(H,2) + 24*G2x*G4*G4ppx*G5xxx*pow(H,2) - 48*G3p*G4*G4ppx*G5xxx*pow(H,2) + 48*G4*G4pp*G4ppx*G5xxx*pow(H,2) - 
   24*G3pp*G4*G4px*G5xxx*pow(H,2) + 48*G4*G4ppp*G4px*G5xxx*pow(H,2) - 32*G2ppx*G4*G4x*G5xxx*pow(H,2) + 32*G3ppp*G4*G4x*G5xxx*pow(H,2) + 24*G3pp*G4p*G4x*G5xxx*pow(H,2) - 
   24*G2x*G4pp*G4x*G5xxx*pow(H,2) + 48*G3p*G4pp*G4x*G5xxx*pow(H,2) - 48*pow(G4pp,2)*G4x*G5xxx*pow(H,2) - 48*G4p*G4ppp*G4x*G5xxx*pow(H,2) - 24*G2p*G4px*G4x*G5xxx*pow(H,2) - 
   16*G2pp*pow(G4x,2)*G5xxx*pow(H,2) + 16*G2ppx*G4*G5p*G5xxx*pow(H,2) - 16*G3ppp*G4*G5p*G5xxx*pow(H,2) - 12*G3pp*G4p*G5p*G5xxx*pow(H,2) + 12*G2x*G4pp*G5p*G5xxx*pow(H,2) - 
   24*G3p*G4pp*G5p*G5xxx*pow(H,2) + 24*pow(G4pp,2)*G5p*G5xxx*pow(H,2) + 24*G4p*G4ppp*G5p*G5xxx*pow(H,2) + 12*G2p*G4px*G5p*G5xxx*pow(H,2) + 16*G2pp*G4x*G5p*G5xxx*pow(H,2) - 
   4*G2pp*pow(G5p,2)*G5xxx*pow(H,2) - 12*G2x*G4p*G5pp*G5xxx*pow(H,2) + 24*G3p*G4p*G5pp*G5xxx*pow(H,2) - 24*G4p*G4pp*G5pp*G5xxx*pow(H,2) - 3456*G4pxxx*pow(G4x,3)*pow(H,4) + 
   6912*G4*G4pxxx*G4x*G4xx*pow(H,4) + 12096*G4pxx*pow(G4x,2)*G4xx*pow(H,4) - 3456*G4*G4pxx*pow(G4xx,2)*pow(H,4) - 3456*G4px*G4x*pow(G4xx,2)*pow(H,4) + 6912*G4p*pow(G4xx,3)*pow(H,4) - 
   5184*G4*G4pxx*G4x*G4xxx*pow(H,4) + 5184*G4p*G4x*G4xx*G4xxx*pow(H,4) + 7776*G4pxxx*pow(G4x,2)*G5p*pow(H,4) - 5184*G4*G4pxxx*G4xx*G5p*pow(H,4) - 23328*G4pxx*G4x*G4xx*G5p*pow(H,4) + 
   4320*G4px*pow(G4xx,2)*G5p*pow(H,4) + 3456*G4*G4pxx*G4xxx*G5p*pow(H,4) + 1728*G4px*G4x*G4xxx*G5p*pow(H,4) - 4320*G4p*G4xx*G4xxx*G5p*pow(H,4) - 5616*G4pxxx*G4x*pow(G5p,2)*pow(H,4) + 
   11232*G4pxx*G4xx*pow(G5p,2)*pow(H,4) - 1296*G4px*G4xxx*pow(G5p,2)*pow(H,4) + 1296*G4pxxx*pow(G5p,3)*pow(H,4) - 7776*G4x*pow(G4xx,2)*G5pp*pow(H,4) - 
   1728*pow(G4x,2)*G4xxx*G5pp*pow(H,4) + 864*G4*G4xx*G4xxx*G5pp*pow(H,4) + 5184*pow(G4xx,2)*G5p*G5pp*pow(H,4) + 1296*G4x*G4xxx*G5p*G5pp*pow(H,4) - 3168*pow(G4x,2)*G4xx*G5ppx*pow(H,4) + 
   2400*G4*pow(G4xx,2)*G5ppx*pow(H,4) + 2496*G4*G4x*G4xxx*G5ppx*pow(H,4) + 7296*G4x*G4xx*G5p*G5ppx*pow(H,4) - 1824*G4*G4xxx*G5p*G5ppx*pow(H,4) - 3720*G4xx*pow(G5p,2)*G5ppx*pow(H,4) + 
   1632*pow(G4x,3)*G5ppxx*pow(H,4) - 3840*G4*G4x*G4xx*G5ppxx*pow(H,4) + 96*pow(G4,2)*G4xxx*G5ppxx*pow(H,4) - 3600*pow(G4x,2)*G5p*G5ppxx*pow(H,4) + 2592*G4*G4xx*G5p*G5ppxx*pow(H,4) + 
   2544*G4x*pow(G5p,2)*G5ppxx*pow(H,4) - 576*pow(G5p,3)*G5ppxx*pow(H,4) - 3360*G4*G4pxxx*G4x*G5px*pow(H,4) + 360*G3xx*pow(G4x,2)*G5px*pow(H,4) - 7632*G4pxx*pow(G4x,2)*G5px*pow(H,4) - 
   360*G3xx*G4*G4xx*G5px*pow(H,4) + 2544*G4*G4pxx*G4xx*G5px*pow(H,4) - 3168*G4px*G4x*G4xx*G5px*pow(H,4) - 11424*G4p*pow(G4xx,2)*G5px*pow(H,4) + 2160*G4*G4px*G4xxx*G5px*pow(H,4) - 
   2304*G4p*G4x*G4xxx*G5px*pow(H,4) + 2688*G4*G4pxxx*G5p*G5px*pow(H,4) - 108*G3xx*G4x*G5p*G5px*pow(H,4) + 12648*G4pxx*G4x*G5p*G5px*pow(H,4) - 3720*G4px*G4xx*G5p*G5px*pow(H,4) + 
   1800*G4p*G4xxx*G5p*G5px*pow(H,4) - 144*G3xx*pow(G5p,2)*G5px*pow(H,4) - 5784*G4pxx*pow(G5p,2)*G5px*pow(H,4) + 9408*G4x*G4xx*G5pp*G5px*pow(H,4) - 1152*G4*G4xxx*G5pp*G5px*pow(H,4) - 
   6408*G4xx*G5p*G5pp*G5px*pow(H,4) + 2016*pow(G4x,2)*G5ppx*G5px*pow(H,4) - 1584*G4*G4xx*G5ppx*G5px*pow(H,4) - 4032*G4x*G5p*G5ppx*G5px*pow(H,4) + 2088*pow(G5p,2)*G5ppx*G5px*pow(H,4) + 
   1968*G4*G4x*G5ppxx*G5px*pow(H,4) - 1392*G4*G5p*G5ppxx*G5px*pow(H,4) + 288*G3xx*G4*pow(G5px,2)*pow(H,4) - 768*G4*G4pxx*pow(G5px,2)*pow(H,4) + 2784*G4px*G4x*pow(G5px,2)*pow(H,4) + 
   6360*G4p*G4xx*pow(G5px,2)*pow(H,4) + 1092*G4px*G5p*pow(G5px,2)*pow(H,4) - 2904*G4x*G5pp*pow(G5px,2)*pow(H,4) + 1872*G5p*G5pp*pow(G5px,2)*pow(H,4) + 
   264*G4*G5ppx*pow(G5px,2)*pow(H,4) - 1152*G4p*pow(G5px,3)*pow(H,4) - 96*pow(G4,2)*G4pxxx*G5pxx*pow(H,4) - 504*G3xx*G4*G4x*G5pxx*pow(H,4) + 4752*G4*G4pxx*G4x*G5pxx*pow(H,4) - 
   720*G4px*pow(G4x,2)*G5pxx*pow(H,4) + 1440*G4*G4px*G4xx*G5pxx*pow(H,4) - 4848*G4p*G4x*G4xx*G5pxx*pow(H,4) + 144*G4*G4p*G4xxx*G5pxx*pow(H,4) + 360*G3xx*G4*G5p*G5pxx*pow(H,4) - 
   2544*G4*G4pxx*G5p*G5pxx*pow(H,4) + 1992*G4px*G4x*G5p*G5pxx*pow(H,4) + 4248*G4p*G4xx*G5p*G5pxx*pow(H,4) - 1224*G4px*pow(G5p,2)*G5pxx*pow(H,4) + 480*pow(G4x,2)*G5pp*G5pxx*pow(H,4) + 
   336*G4*G4xx*G5pp*G5pxx*pow(H,4) - 48*G4x*G5p*G5pp*G5pxx*pow(H,4) - 108*pow(G5p,2)*G5pp*G5pxx*pow(H,4) - 1920*G4*G4x*G5ppx*G5pxx*pow(H,4) + 1008*G4*G5p*G5ppx*G5pxx*pow(H,4) - 
   1968*G4*G4px*G5px*G5pxx*pow(H,4) + 2184*G4p*G4x*G5px*G5pxx*pow(H,4) - 2016*G4p*G5p*G5px*G5pxx*pow(H,4) + 264*G4*G5pp*G5px*G5pxx*pow(H,4) - 24*G4*G4p*pow(G5pxx,2)*pow(H,4) - 
   96*pow(G4,2)*G4pxx*G5pxxx*pow(H,4) - 1200*G4*G4px*G4x*G5pxxx*pow(H,4) + 432*G4p*pow(G4x,2)*G5pxxx*pow(H,4) - 240*G4*G4p*G4xx*G5pxxx*pow(H,4) + 912*G4*G4px*G5p*G5pxxx*pow(H,4) - 
   552*G4p*G4x*G5p*G5pxxx*pow(H,4) + 168*G4p*pow(G5p,2)*G5pxxx*pow(H,4) - 48*G4*G5p*G5pp*G5pxxx*pow(H,4) + 48*pow(G4,2)*G5ppx*G5pxxx*pow(H,4) + 144*G4*G4p*G5px*G5pxxx*pow(H,4) - 
   1296*G3xx*G4*G4pxx*G5x*pow(H,4) + 3168*G4*pow(G4pxx,2)*G5x*pow(H,4) - 4320*G4*G4px*G4pxxx*G5x*pow(H,4) - 432*G3xx*G4px*G4x*G5x*pow(H,4) - 9216*G4px*G4pxx*G4x*G5x*pow(H,4) + 
   2736*G4p*G4pxxx*G4x*G5x*pow(H,4) - 2808*G3pxx*pow(G4x,2)*G5x*pow(H,4) + 4464*G4ppxx*pow(G4x,2)*G5x*pow(H,4) + 1944*G3pxx*G4*G4xx*G5x*pow(H,4) + 1512*G3xx*G4p*G4xx*G5x*pow(H,4) - 
   4464*G4*G4ppxx*G4xx*G5x*pow(H,4) + 3744*pow(G4px,2)*G4xx*G5x*pow(H,4) - 7920*G4p*G4pxx*G4xx*G5x*pow(H,4) + 2016*G4ppx*G4x*G4xx*G5x*pow(H,4) - 8640*G4pp*pow(G4xx,2)*G5x*pow(H,4) + 
   3024*G4*G4ppx*G4xxx*G5x*pow(H,4) - 2016*G4p*G4px*G4xxx*G5x*pow(H,4) - 576*G4pp*G4x*G4xxx*G5x*pow(H,4) + 756*G3xx*G4px*G5p*G5x*pow(H,4) + 10296*G4px*G4pxx*G5p*G5x*pow(H,4) - 
   1728*G4p*G4pxxx*G5p*G5x*pow(H,4) + 4212*G3pxx*G4x*G5p*G5x*pow(H,4) - 6120*G4ppxx*G4x*G5p*G5x*pow(H,4) - 2808*G4ppx*G4xx*G5p*G5x*pow(H,4) + 216*G4pp*G4xxx*G5p*G5x*pow(H,4) - 
   1512*G3pxx*pow(G5p,2)*G5x*pow(H,4) + 2016*G4ppxx*pow(G5p,2)*G5x*pow(H,4) + 144*G4*G4pxxx*G5pp*G5x*pow(H,4) - 648*G3xx*G4x*G5pp*G5x*pow(H,4) - 720*G4pxx*G4x*G5pp*G5x*pow(H,4) + 
   10152*G4px*G4xx*G5pp*G5x*pow(H,4) + 72*G4p*G4xxx*G5pp*G5x*pow(H,4) + 162*G3xx*G5p*G5pp*G5x*pow(H,4) + 684*G4pxx*G5p*G5pp*G5x*pow(H,4) - 1296*G4xx*pow(G5pp,2)*G5x*pow(H,4) - 
   1152*G4x*G4xx*G5ppp*G5x*pow(H,4) - 144*G4*G4xxx*G5ppp*G5x*pow(H,4) + 2160*G4xx*G5p*G5ppp*G5x*pow(H,4) + 576*pow(G4x,2)*G5pppx*G5x*pow(H,4) + 288*G4*G4xx*G5pppx*G5x*pow(H,4) - 
   1152*G4x*G5p*G5pppx*G5x*pow(H,4) + 504*pow(G5p,2)*G5pppx*G5x*pow(H,4) + 660*G3xx*G4*G5ppx*G5x*pow(H,4) - 2184*G4*G4pxx*G5ppx*G5x*pow(H,4) + 1848*G4px*G4x*G5ppx*G5x*pow(H,4) + 
   1272*G4p*G4xx*G5ppx*G5x*pow(H,4) - 3540*G4px*G5p*G5ppx*G5x*pow(H,4) + 1224*G4x*G5pp*G5ppx*G5x*pow(H,4) - 684*G5p*G5pp*G5ppx*G5x*pow(H,4) + 288*G4*pow(G5ppx,2)*G5x*pow(H,4) + 
   2304*G4*G4px*G5ppxx*G5x*pow(H,4) - 1128*G4p*G4x*G5ppxx*G5x*pow(H,4) + 696*G4p*G5p*G5ppxx*G5x*pow(H,4) - 72*G4*G5pp*G5ppxx*G5x*pow(H,4) - 984*G3pxx*G4*G5px*G5x*pow(H,4) - 
   684*G3xx*G4p*G5px*G5x*pow(H,4) + 2544*G4*G4ppxx*G5px*G5x*pow(H,4) - 1308*pow(G4px,2)*G5px*G5x*pow(H,4) + 3912*G4p*G4pxx*G5px*G5x*pow(H,4) + 60*G2xx*G4x*G5px*G5x*pow(H,4) - 
   48*G4ppx*G4x*G5px*G5x*pow(H,4) + 240*G2x*G4xx*G5px*G5x*pow(H,4) - 480*G3p*G4xx*G5px*G5x*pow(H,4) + 9816*G4pp*G4xx*G5px*G5x*pow(H,4) + 30*G2xx*G5p*G5px*G5x*pow(H,4) + 
   732*G4ppx*G5p*G5px*G5x*pow(H,4) - 5676*G4px*G5pp*G5px*G5x*pow(H,4) + 684*pow(G5pp,2)*G5px*G5x*pow(H,4) + 360*G4x*G5ppp*G5px*G5x*pow(H,4) - 900*G5p*G5ppp*G5px*G5x*pow(H,4) - 
   288*G4*G5pppx*G5px*G5x*pow(H,4) - 684*G4p*G5ppx*G5px*G5x*pow(H,4) - 120*G2x*pow(G5px,2)*G5x*pow(H,4) + 240*G3p*pow(G5px,2)*G5x*pow(H,4) - 2784*G4pp*pow(G5px,2)*G5x*pow(H,4) - 
   84*G2xx*G4*G5pxx*G5x*pow(H,4) - 1896*G4*G4ppx*G5pxx*G5x*pow(H,4) + 2208*G4p*G4px*G5pxx*G5x*pow(H,4) + 48*G2x*G4x*G5pxx*G5x*pow(H,4) - 96*G3p*G4x*G5pxx*G5x*pow(H,4) - 
   144*G4pp*G4x*G5pxx*G5x*pow(H,4) - 156*G2x*G5p*G5pxx*G5x*pow(H,4) + 312*G3p*G5p*G5pxx*G5x*pow(H,4) + 72*G4pp*G5p*G5pxx*G5x*pow(H,4) + 12*G4p*G5pp*G5pxx*G5x*pow(H,4) + 
   72*G4*G5ppp*G5pxx*G5x*pow(H,4) + 84*G2x*G4*G5pxxx*G5x*pow(H,4) - 168*G3p*G4*G5pxxx*G5x*pow(H,4) - 24*pow(G4p,2)*G5pxxx*G5x*pow(H,4) + 72*G4*G4pp*G5pxxx*G5x*pow(H,4) + 
   360*G3pxx*G4p*pow(G5x,2)*pow(H,4) - 36*G3xx*G4pp*pow(G5x,2)*pow(H,4) - 288*G4p*G4ppxx*pow(G5x,2)*pow(H,4) - 72*G2xx*G4px*pow(G5x,2)*pow(H,4) - 252*G4ppx*G4px*pow(G5x,2)*pow(H,4) + 
   324*G2x*G4pxx*pow(G5x,2)*pow(H,4) - 648*G3p*G4pxx*pow(G5x,2)*pow(H,4) - 288*G4pp*G4pxx*pow(G5x,2)*pow(H,4) - 504*G2pxx*G4x*pow(G5x,2)*pow(H,4) + 108*G3ppx*G4x*pow(G5x,2)*pow(H,4) + 
   792*G4pppx*G4x*pow(G5x,2)*pow(H,4) - 396*G2px*G4xx*pow(G5x,2)*pow(H,4) + 1188*G3pp*G4xx*pow(G5x,2)*pow(H,4) - 792*G4ppp*G4xx*pow(G5x,2)*pow(H,4) + 
   108*G2p*G4xxx*pow(G5x,2)*pow(H,4) + 378*G2pxx*G5p*pow(G5x,2)*pow(H,4) - 756*G4pppx*G5p*pow(G5x,2)*pow(H,4) - 36*G2xx*G5pp*pow(G5x,2)*pow(H,4) + 576*G4ppx*G5pp*pow(G5x,2)*pow(H,4) + 
   180*G4px*G5ppp*pow(G5x,2)*pow(H,4) - 216*G4p*G5pppx*pow(G5x,2)*pow(H,4) - 82*G2x*G5ppx*pow(G5x,2)*pow(H,4) + 164*G3p*G5ppx*pow(G5x,2)*pow(H,4) + 
   288*G4pp*G5ppx*pow(G5x,2)*pow(H,4) + 202*G2px*G5px*pow(G5x,2)*pow(H,4) - 518*G3pp*G5px*pow(G5x,2)*pow(H,4) + 288*G4ppp*G5px*pow(G5x,2)*pow(H,4) - 62*G2p*G5pxx*pow(G5x,2)*pow(H,4) - 
   30*G2ppx*pow(G5x,3)*pow(H,4) + 30*G3ppp*pow(G5x,3)*pow(H,4) + 72*G3xx*G4*G4px*G5xx*pow(H,4) + 1488*G4*G4px*G4pxx*G5xx*pow(H,4) - 288*G4*G4p*G4pxxx*G5xx*pow(H,4) + 
   648*G3pxx*G4*G4x*G5xx*pow(H,4) + 576*G3xx*G4p*G4x*G5xx*pow(H,4) - 2064*G4*G4ppxx*G4x*G5xx*pow(H,4) + 1152*pow(G4px,2)*G4x*G5xx*pow(H,4) - 2976*G4p*G4pxx*G4x*G5xx*pow(H,4) - 
   912*G4ppx*pow(G4x,2)*G5xx*pow(H,4) + 1632*G4*G4ppx*G4xx*G5xx*pow(H,4) - 9744*G4p*G4px*G4xx*G5xx*pow(H,4) - 5136*G4pp*G4x*G4xx*G5xx*pow(H,4) - 288*pow(G4p,2)*G4xxx*G5xx*pow(H,4) - 
   144*G4*G4pp*G4xxx*G5xx*pow(H,4) - 504*G3pxx*G4*G5p*G5xx*pow(H,4) - 468*G3xx*G4p*G5p*G5xx*pow(H,4) + 1104*G4*G4ppxx*G5p*G5xx*pow(H,4) - 2400*pow(G4px,2)*G5p*G5xx*pow(H,4) + 
   2568*G4p*G4pxx*G5p*G5xx*pow(H,4) + 504*G4ppx*G4x*G5p*G5xx*pow(H,4) + 5448*G4pp*G4xx*G5p*G5xx*pow(H,4) + 120*G4ppx*pow(G5p,2)*G5xx*pow(H,4) + 36*G3xx*G4*G5pp*G5xx*pow(H,4) + 
   696*G4*G4pxx*G5pp*G5xx*pow(H,4) + 4272*G4px*G4x*G5pp*G5xx*pow(H,4) + 1560*G4p*G4xx*G5pp*G5xx*pow(H,4) - 2328*G4px*G5p*G5pp*G5xx*pow(H,4) - 552*G4x*pow(G5pp,2)*G5xx*pow(H,4) + 
   432*G5p*pow(G5pp,2)*G5xx*pow(H,4) - 96*pow(G4x,2)*G5ppp*G5xx*pow(H,4) - 1536*G4*G4xx*G5ppp*G5xx*pow(H,4) + 408*G4x*G5p*G5ppp*G5xx*pow(H,4) - 396*pow(G5p,2)*G5ppp*G5xx*pow(H,4) + 
   384*G4*G4x*G5pppx*G5xx*pow(H,4) - 48*G4*G5p*G5pppx*G5xx*pow(H,4) - 864*G4*G4px*G5ppx*G5xx*pow(H,4) + 600*G4p*G4x*G5ppx*G5xx*pow(H,4) - 444*G4p*G5p*G5ppx*G5xx*pow(H,4) - 
   456*G4*G5pp*G5ppx*G5xx*pow(H,4) + 96*G4*G4p*G5ppxx*G5xx*pow(H,4) - 12*G2xx*G4*G5px*G5xx*pow(H,4) - 456*G4*G4ppx*G5px*G5xx*pow(H,4) + 5160*G4p*G4px*G5px*G5xx*pow(H,4) + 
   216*G2x*G4x*G5px*G5xx*pow(H,4) - 432*G3p*G4x*G5px*G5xx*pow(H,4) + 3096*G4pp*G4x*G5px*G5xx*pow(H,4) - 12*G2x*G5p*G5px*G5xx*pow(H,4) + 24*G3p*G5p*G5px*G5xx*pow(H,4) - 
   3072*G4pp*G5p*G5px*G5xx*pow(H,4) - 1008*G4p*G5pp*G5px*G5xx*pow(H,4) + 840*G4*G5ppp*G5px*G5xx*pow(H,4) - 84*G2x*G4*G5pxx*G5xx*pow(H,4) + 168*G3p*G4*G5pxx*G5xx*pow(H,4) + 
   204*pow(G4p,2)*G5pxx*G5xx*pow(H,4) + 120*G4*G4pp*G5pxx*G5xx*pow(H,4) + 120*G2pxx*G4*G5x*G5xx*pow(H,4) - 240*G3ppx*G4*G5x*G5xx*pow(H,4) + 108*G2xx*G4p*G5x*G5xx*pow(H,4) + 
   240*G4*G4pppx*G5x*G5xx*pow(H,4) - 1008*G4p*G4ppx*G5x*G5xx*pow(H,4) - 132*G2x*G4px*G5x*G5xx*pow(H,4) + 264*G3p*G4px*G5x*G5xx*pow(H,4) + 3432*G4pp*G4px*G5x*G5xx*pow(H,4) - 
   204*G2px*G4x*G5x*G5xx*pow(H,4) + 528*G3pp*G4x*G5x*G5xx*pow(H,4) - 240*G4ppp*G4x*G5x*G5xx*pow(H,4) - 336*G2p*G4xx*G5x*G5xx*pow(H,4) + 204*G2px*G5p*G5x*G5xx*pow(H,4) - 
   690*G3pp*G5p*G5x*G5xx*pow(H,4) + 564*G4ppp*G5p*G5x*G5xx*pow(H,4) - 192*G2x*G5pp*G5x*G5xx*pow(H,4) + 384*G3p*G5pp*G5x*G5xx*pow(H,4) - 684*G4pp*G5pp*G5x*G5xx*pow(H,4) + 
   408*G4p*G5ppp*G5x*G5xx*pow(H,4) + 150*G2p*G5px*G5x*G5xx*pow(H,4) + 8*G2pp*pow(G5x,2)*G5xx*pow(H,4) + 48*G2px*G4*pow(G5xx,2)*pow(H,4) - 12*G3pp*G4*pow(G5xx,2)*pow(H,4) + 
   132*G2x*G4p*pow(G5xx,2)*pow(H,4) - 264*G3p*G4p*pow(G5xx,2)*pow(H,4) + 492*G4p*G4pp*pow(G5xx,2)*pow(H,4) - 168*G4*G4ppp*pow(G5xx,2)*pow(H,4) - 24*G2p*G4x*pow(G5xx,2)*pow(H,4) + 
   66*G2p*G5p*pow(G5xx,2)*pow(H,4) + 96*pow(G4,2)*G4ppxx*G5xxx*pow(H,4) + 336*G4*pow(G4px,2)*G5xxx*pow(H,4) + 96*G4*G4p*G4pxx*G5xxx*pow(H,4) + 912*G4*G4ppx*G4x*G5xxx*pow(H,4) - 
   480*G4p*G4px*G4x*G5xxx*pow(H,4) - 144*G4pp*pow(G4x,2)*G5xxx*pow(H,4) - 240*pow(G4p,2)*G4xx*G5xxx*pow(H,4) - 48*G4*G4pp*G4xx*G5xxx*pow(H,4) - 624*G4*G4ppx*G5p*G5xxx*pow(H,4) + 
   432*G4p*G4px*G5p*G5xxx*pow(H,4) + 120*G4pp*G4x*G5p*G5xxx*pow(H,4) - 24*G4pp*pow(G5p,2)*G5xxx*pow(H,4) - 504*G4*G4px*G5pp*G5xxx*pow(H,4) + 12*G4p*G5p*G5pp*G5xxx*pow(H,4) + 
   96*G4*pow(G5pp,2)*G5xxx*pow(H,4) + 48*G4*G5p*G5ppp*G5xxx*pow(H,4) - 48*pow(G4,2)*G5pppx*G5xxx*pow(H,4) - 72*G4*G4p*G5ppx*G5xxx*pow(H,4) - 36*G2x*G4*G5px*G5xxx*pow(H,4) + 
   72*G3p*G4*G5px*G5xxx*pow(H,4) + 120*pow(G4p,2)*G5px*G5xxx*pow(H,4) - 72*G4*G4pp*G5px*G5xxx*pow(H,4) - 48*G2px*G4*G5x*G5xxx*pow(H,4) + 132*G3pp*G4*G5x*G5xxx*pow(H,4) + 
   36*G2x*G4p*G5x*G5xxx*pow(H,4) - 72*G3p*G4p*G5x*G5xxx*pow(H,4) + 48*G4p*G4pp*G5x*G5xxx*pow(H,4) - 72*G4*G4ppp*G5x*G5xxx*pow(H,4) + 60*G2p*G4x*G5x*G5xxx*pow(H,4) - 
   42*G2p*G5p*G5x*G5xxx*pow(H,4) - 12*G2p*G4*G5xx*G5xxx*pow(H,4) + 10080*G4x*G4xx*G5px*G5x*pow(H,6) - 1080*G4*G4xxx*G5px*G5x*pow(H,6) - 7560*G4xx*G5p*G5px*G5x*pow(H,6) - 
   5340*G4x*pow(G5px,2)*G5x*pow(H,6) + 3930*G5p*pow(G5px,2)*G5x*pow(H,6) - 1512*pow(G4x,2)*G5pxx*G5x*pow(H,6) + 792*G4*G4xx*G5pxx*G5x*pow(H,6) + 1152*G4x*G5p*G5pxx*G5x*pow(H,6) + 
   180*pow(G5p,2)*G5pxx*G5x*pow(H,6) + 176*G4*G5px*G5pxx*G5x*pow(H,6) + 936*G4*G4x*G5pxxx*G5x*pow(H,6) - 792*G4*G5p*G5pxxx*G5x*pow(H,6) + 1512*G4*G4pxxx*pow(G5x,2)*pow(H,6) - 
   864*G4pxx*G4x*pow(G5x,2)*pow(H,6) + 2160*G4px*G4xx*pow(G5x,2)*pow(H,6) + 1296*G4p*G4xxx*pow(G5x,2)*pow(H,6) - 540*G4pxx*G5p*pow(G5x,2)*pow(H,6) - 
   5616*G4xx*G5pp*pow(G5x,2)*pow(H,6) + 2184*G4x*G5ppx*pow(G5x,2)*pow(H,6) - 1224*G5p*G5ppx*pow(G5x,2)*pow(H,6) - 744*G4*G5ppxx*pow(G5x,2)*pow(H,6) - 
   3594*G4px*G5px*pow(G5x,2)*pow(H,6) + 3216*G5pp*G5px*pow(G5x,2)*pow(H,6) - 690*G4p*G5pxx*pow(G5x,2)*pow(H,6) + 1350*G4ppx*pow(G5x,3)*pow(H,6) - 180*G5ppp*pow(G5x,3)*pow(H,6) + 
   2280*pow(G4x,2)*G5px*G5xx*pow(H,6) - 1704*G4*G4xx*G5px*G5xx*pow(H,6) - 3096*G4x*G5p*G5px*G5xx*pow(H,6) + 996*pow(G5p,2)*G5px*G5xx*pow(H,6) + 924*G4*pow(G5px,2)*G5xx*pow(H,6) - 
   168*G4*G4x*G5pxx*G5xx*pow(H,6) + 168*G4*G5p*G5pxx*G5xx*pow(H,6) - 48*pow(G4,2)*G5pxxx*G5xx*pow(H,6) - 864*G4*G4pxx*G5x*G5xx*pow(H,6) + 1440*G4px*G4x*G5x*G5xx*pow(H,6) + 
   5976*G4p*G4xx*G5x*G5xx*pow(H,6) - 900*G4px*G5p*G5x*G5xx*pow(H,6) - 4428*G4x*G5pp*G5x*G5xx*pow(H,6) + 3402*G5p*G5pp*G5x*G5xx*pow(H,6) + 316*G4*G5ppx*G5x*G5xx*pow(H,6) - 
   3972*G4p*G5px*G5x*G5xx*pow(H,6) - 1470*G4pp*pow(G5x,2)*G5xx*pow(H,6) + 168*G4*G4px*pow(G5xx,2)*pow(H,6) + 1272*G4p*G4x*pow(G5xx,2)*pow(H,6) - 1020*G4p*G5p*pow(G5xx,2)*pow(H,6) + 
   84*G4*G5pp*pow(G5xx,2)*pow(H,6) - 408*G4*G4x*G5px*G5xxx*pow(H,6) + 264*G4*G5p*G5px*G5xxx*pow(H,6) + 48*pow(G4,2)*G5pxx*G5xxx*pow(H,6) - 216*G4*G4px*G5x*G5xxx*pow(H,6) + 
   720*G4p*G4x*G5x*G5xxx*pow(H,6) - 612*G4p*G5p*G5x*G5xxx*pow(H,6) + 324*G4*G5pp*G5x*G5xxx*pow(H,6) - 72*G4*G4p*G5xx*G5xxx*pow(H,6) + 675*G5px*pow(G5x,3)*pow(H,8) - 
   9*pow(G3x,3)*(2*G4ppx + 2*G5ppp + 7*G5px*pow(H,2)) + 6*pow(G3px,2)*(6*G3x*G4x - 4*G4px*G4x - 3*G3x*G5p + 2*G4px*G5p - 4*G4x*G5pp + 2*G5p*G5pp - G2x*G5x + 2*G3p*G5x - 2*G4pp*G5x - 
      24*G4x*G5x*pow(H,2) + 3*G5p*G5x*pow(H,2) - 10*G4*G5xx*pow(H,2)) - 3*pow(G3x,2)*
    (12*G2xx*G4px - 84*G4ppx*G4px - 12*G2x*G4pxx + 24*G3p*G4pxx - 24*G4pp*G4pxx + 12*G2pxx*G4x - 24*G4pppx*G4x - 36*G3pp*G4xx + 72*G4ppp*G4xx - 6*G2pxx*G5p + 12*G4pppx*G5p - 6*G2xx*G5pp + 24*G4ppx*G5pp - 
      36*G4px*G5ppp + 6*G2x*G5ppx - 12*G3p*G5ppx + 12*G4pp*G5ppx + 2*G2px*G5px + 16*G3pp*G5px - 36*G4ppp*G5px + 2*G2ppx*G5x - 2*G3ppp*G5x - 72*G4*G4pxxx*pow(H,2) - 576*G4pxx*G4x*pow(H,2) + 
      504*G4px*G4xx*pow(H,2) - 72*G4p*G4xxx*pow(H,2) + 468*G4pxx*G5p*pow(H,2) - 72*G4xx*G5pp*pow(H,2) + 248*G4x*G5ppx*pow(H,2) - 196*G5p*G5ppx*pow(H,2) + 32*G4*G5ppxx*pow(H,2) - 
      466*G4px*G5px*pow(H,2) + 50*G5pp*G5px*pow(H,2) + 74*G4p*G5pxx*pow(H,2) + 66*G4ppx*G5x*pow(H,2) + 24*G5ppp*G5x*pow(H,2) + 78*G4pp*G5xx*pow(H,2) - 75*G5px*G5x*pow(H,4)) - 
   288*G4pxxx*G4x*G4xx*pow(H,2)*rptot + 576*G4pxx*pow(G4xx,2)*pow(H,2)*rptot + 288*G4pxx*G4x*G4xxx*pow(H,2)*rptot - 288*G4px*G4xx*G4xxx*pow(H,2)*rptot + 144*G4pxxx*G4xx*G5p*pow(H,2)*rptot - 
   144*G4pxx*G4xxx*G5p*pow(H,2)*rptot + 144*G4xx*G4xxx*G5pp*pow(H,2)*rptot - 288*pow(G4xx,2)*G5ppx*pow(H,2)*rptot - 144*G4x*G4xxx*G5ppx*pow(H,2)*rptot + 72*G4xxx*G5p*G5ppx*pow(H,2)*rptot + 
   144*G4x*G4xx*G5ppxx*pow(H,2)*rptot - 72*G4xx*G5p*G5ppxx*pow(H,2)*rptot + 144*G4pxxx*G4x*G5px*pow(H,2)*rptot - 72*G3xx*G4xx*G5px*pow(H,2)*rptot - 432*G4pxx*G4xx*G5px*pow(H,2)*rptot + 
   216*G4px*G4xxx*G5px*pow(H,2)*rptot - 72*G4pxxx*G5p*G5px*pow(H,2)*rptot - 72*G4xxx*G5pp*G5px*pow(H,2)*rptot + 288*G4xx*G5ppx*G5px*pow(H,2)*rptot - 72*G4x*G5ppxx*G5px*pow(H,2)*rptot + 
   36*G5p*G5ppxx*G5px*pow(H,2)*rptot + 36*G3xx*pow(G5px,2)*pow(H,2)*rptot + 72*G4pxx*pow(G5px,2)*pow(H,2)*rptot - 72*G5ppx*pow(G5px,2)*pow(H,2)*rptot + 36*G3xx*G4x*G5pxx*pow(H,2)*rptot - 
   216*G4pxx*G4x*G5pxx*pow(H,2)*rptot - 18*G3xx*G5p*G5pxx*pow(H,2)*rptot + 108*G4pxx*G5p*G5pxx*pow(H,2)*rptot - 72*G4xx*G5pp*G5pxx*pow(H,2)*rptot + 72*G4x*G5ppx*G5pxx*pow(H,2)*rptot - 
   36*G5p*G5ppx*G5pxx*pow(H,2)*rptot - 36*G4px*G5px*G5pxx*pow(H,2)*rptot + 36*G5pp*G5px*G5pxx*pow(H,2)*rptot + 24*G4px*G4x*G5pxxx*pow(H,2)*rptot - 12*G4px*G5p*G5pxxx*pow(H,2)*rptot + 
   72*G3xx*G4pxx*G5x*pow(H,2)*rptot - 144*pow(G4pxx,2)*G5x*pow(H,2)*rptot + 72*G4px*G4pxxx*G5x*pow(H,2)*rptot - 72*G3pxx*G4xx*G5x*pow(H,2)*rptot + 144*G4ppxx*G4xx*G5x*pow(H,2)*rptot - 
   72*G4ppx*G4xxx*G5x*pow(H,2)*rptot - 36*G3xx*G5ppx*G5x*pow(H,2)*rptot + 72*G4pxx*G5ppx*G5x*pow(H,2)*rptot - 36*G4px*G5ppxx*G5x*pow(H,2)*rptot + 36*G3pxx*G5px*G5x*pow(H,2)*rptot - 
   72*G4ppxx*G5px*G5x*pow(H,2)*rptot + 6*G2xx*G5pxx*G5x*pow(H,2)*rptot + 36*G4ppx*G5pxx*G5x*pow(H,2)*rptot - 36*G3xx*G4px*G5xx*pow(H,2)*rptot - 72*G4px*G4pxx*G5xx*pow(H,2)*rptot - 
   36*G3pxx*G4x*G5xx*pow(H,2)*rptot + 72*G4ppxx*G4x*G5xx*pow(H,2)*rptot - 144*G4ppx*G4xx*G5xx*pow(H,2)*rptot + 18*G3pxx*G5p*G5xx*pow(H,2)*rptot - 36*G4ppxx*G5p*G5xx*pow(H,2)*rptot + 
   18*G3xx*G5pp*G5xx*pow(H,2)*rptot - 36*G4pxx*G5pp*G5xx*pow(H,2)*rptot + 72*G4px*G5ppx*G5xx*pow(H,2)*rptot - 6*G2xx*G5px*G5xx*pow(H,2)*rptot + 72*G4ppx*G5px*G5xx*pow(H,2)*rptot - 
   6*G2pxx*G5x*G5xx*pow(H,2)*rptot + 6*G3ppx*G5x*G5xx*pow(H,2)*rptot + 24*pow(G4px,2)*G5xxx*pow(H,2)*rptot - 24*G4ppx*G4x*G5xxx*pow(H,2)*rptot + 12*G4ppx*G5p*G5xxx*pow(H,2)*rptot - 
   12*G4px*G5pp*G5xxx*pow(H,2)*rptot + 192*G4xx*G5pxx*G5x*pow(H,4)*rptot - 114*G5px*G5pxx*G5x*pow(H,4)*rptot - 60*G4x*G5pxxx*G5x*pow(H,4)*rptot + 42*G5p*G5pxxx*G5x*pow(H,4)*rptot - 
   108*G4pxxx*pow(G5x,2)*pow(H,4)*rptot + 54*G5ppxx*pow(G5x,2)*pow(H,4)*rptot - 96*G4xx*G5px*G5xx*pow(H,4)*rptot + 66*pow(G5px,2)*G5xx*pow(H,4)*rptot + 72*G4x*G5pxx*G5xx*pow(H,4)*rptot - 
   72*G5p*G5pxx*G5xx*pow(H,4)*rptot + 12*G4*G5pxxx*G5xx*pow(H,4)*rptot + 240*G4pxx*G5x*G5xx*pow(H,4)*rptot - 102*G5ppx*G5x*G5xx*pow(H,4)*rptot - 48*G4px*pow(G5xx,2)*pow(H,4)*rptot + 
   6*G5pp*pow(G5xx,2)*pow(H,4)*rptot + 12*G4x*G5px*G5xxx*pow(H,4)*rptot + 6*G5p*G5px*G5xxx*pow(H,4)*rptot - 12*G4*G5pxx*G5xxx*pow(H,4)*rptot - 12*G4px*G5x*G5xxx*pow(H,4)*rptot - 
   6*G5pp*G5x*G5xxx*pow(H,4)*rptot + 12*G4p*G5xx*G5xxx*pow(H,4)*rptot + G3px*
    (27*pow(G3x,3) - 264*pow(G4px,3) + 36*G2x*G3xx*G4x - 72*G3p*G3xx*G4x + 72*G3xx*G4pp*G4x - 120*G2x*G4pxx*G4x + 240*G3p*G4pxx*G4x - 240*G4pp*G4pxx*G4x - 48*G2pxx*pow(G4x,2) + 
      96*G4pppx*pow(G4x,2) + 96*G2px*G4x*G4xx - 48*G3pp*G4x*G4xx - 96*G4ppp*G4x*G4xx - 18*G2x*G3xx*G5p + 36*G3p*G3xx*G5p - 36*G3xx*G4pp*G5p + 60*G2x*G4pxx*G5p - 120*G3p*G4pxx*G5p + 120*G4pp*G4pxx*G5p + 
      48*G2pxx*G4x*G5p - 96*G4pppx*G4x*G5p - 48*G2px*G4xx*G5p + 24*G3pp*G4xx*G5p + 48*G4ppp*G4xx*G5p - 12*G2pxx*pow(G5p,2) + 24*G4pppx*pow(G5p,2) + 24*G2xx*G4x*G5pp + 48*G4ppx*G4x*G5pp - 
      24*G2x*G4xx*G5pp + 48*G3p*G4xx*G5pp - 48*G4pp*G4xx*G5pp - 12*G2xx*G5p*G5pp - 24*G4ppx*G5p*G5pp + 24*G2x*G4x*G5ppx - 48*G3p*G4x*G5ppx + 48*G4pp*G4x*G5ppx - 12*G2x*G5p*G5ppx + 24*G3p*G5p*G5ppx - 
      24*G4pp*G5p*G5ppx - 48*G2px*G4x*G5px + 24*G3pp*G4x*G5px + 48*G4ppp*G4x*G5px + 24*G2px*G5p*G5px - 12*G3pp*G5p*G5px - 24*G4ppp*G5p*G5px + 12*G2x*G5pp*G5px - 24*G3p*G5pp*G5px + 24*G4pp*G5pp*G5px + 
      6*G2x*G2xx*G5x - 12*G2xx*G3p*G5x + 12*G2xx*G4pp*G5x + 12*G2x*G4ppx*G5x - 24*G3p*G4ppx*G5x + 24*G4pp*G4ppx*G5x - 16*G2ppx*G4x*G5x + 16*G3ppp*G4x*G5x + 8*G2ppx*G5p*G5x - 8*G3ppp*G5p*G5x + 
      576*G4*G4pxxx*G4x*pow(H,2) + 648*G3xx*pow(G4x,2)*pow(H,2) - 1872*G4pxx*pow(G4x,2)*pow(H,2) - 216*G3xx*G4*G4xx*pow(H,2) - 1008*G4*G4pxx*G4xx*pow(H,2) - 
      864*G4p*pow(G4xx,2)*pow(H,2) - 144*G4p*G4x*G4xxx*pow(H,2) - 288*G4*G4pxxx*G5p*pow(H,2) - 756*G3xx*G4x*G5p*pow(H,2) + 1368*G4pxx*G4x*G5p*pow(H,2) + 72*G4p*G4xxx*G5p*pow(H,2) + 
      216*G3xx*pow(G5p,2)*pow(H,2) - 216*G4pxx*pow(G5p,2)*pow(H,2) - 432*G4x*G4xx*G5pp*pow(H,2) - 144*G4*G4xxx*G5pp*pow(H,2) - 216*G4xx*G5p*G5pp*pow(H,2) + 
      448*pow(G4x,2)*G5ppx*pow(H,2) + 672*G4*G4xx*G5ppx*pow(H,2) - 112*G4x*G5p*G5ppx*pow(H,2) - 56*pow(G5p,2)*G5ppx*pow(H,2) - 256*G4*G4x*G5ppxx*pow(H,2) + 128*G4*G5p*G5ppxx*pow(H,2) + 
      144*G3xx*G4*G5px*pow(H,2) + 384*G4*G4pxx*G5px*pow(H,2) + 1224*G4p*G4xx*G5px*pow(H,2) + 96*G4x*G5pp*G5px*pow(H,2) + 132*G5p*G5pp*G5px*pow(H,2) - 312*G4*G5ppx*G5px*pow(H,2) - 
      384*G4p*pow(G5px,2)*pow(H,2) - 96*G4p*G4x*G5pxx*pow(H,2) + 48*G4p*G5p*G5pxx*pow(H,2) + 96*G4*G5pp*G5pxx*pow(H,2) + 144*G3pxx*G4*G5x*pow(H,2) - 36*G3xx*G4p*G5x*pow(H,2) - 
      192*G4*G4ppxx*G5x*pow(H,2) - 312*G4p*G4pxx*G5x*pow(H,2) + 216*G2xx*G4x*G5x*pow(H,2) + 672*G4ppx*G4x*G5x*pow(H,2) + 396*G2x*G4xx*G5x*pow(H,2) - 792*G3p*G4xx*G5x*pow(H,2) + 
      120*G4pp*G4xx*G5x*pow(H,2) - 126*G2xx*G5p*G5x*pow(H,2) - 132*G4ppx*G5p*G5x*pow(H,2) + 12*pow(G5pp,2)*G5x*pow(H,2) - 72*G4x*G5ppp*G5x*pow(H,2) + 60*G5p*G5ppp*G5x*pow(H,2) - 
      48*G4*G5pppx*G5x*pow(H,2) + 180*G4p*G5ppx*G5x*pow(H,2) - 222*G2x*G5px*G5x*pow(H,2) + 444*G3p*G5px*G5x*pow(H,2) - 108*G4pp*G5px*G5x*pow(H,2) - 12*G2px*pow(G5x,2)*pow(H,2) + 
      54*G3pp*pow(G5x,2)*pow(H,2) - 36*G4ppp*pow(G5x,2)*pow(H,2) - 12*G2xx*G4*G5xx*pow(H,2) + 216*G4*G4ppx*G5xx*pow(H,2) + 168*G2x*G4x*G5xx*pow(H,2) - 336*G3p*G4x*G5xx*pow(H,2) - 
      120*G2x*G5p*G5xx*pow(H,2) + 240*G3p*G5p*G5xx*pow(H,2) - 72*G4pp*G5p*G5xx*pow(H,2) - 144*G4p*G5pp*G5xx*pow(H,2) + 24*G4*G5ppp*G5xx*pow(H,2) - 18*G2p*G5x*G5xx*pow(H,2) - 
      12*G2x*G4*G5xxx*pow(H,2) + 24*G3p*G4*G5xxx*pow(H,2) - 24*G4*G4pp*G5xxx*pow(H,2) + 3240*G4x*G4xx*G5x*pow(H,4) - 1080*G4*G4xxx*G5x*pow(H,4) - 3024*G4xx*G5p*G5x*pow(H,4) - 
      1884*G4x*G5px*G5x*pow(H,4) + 1644*G5p*G5px*G5x*pow(H,4) + 768*G4*G5pxx*G5x*pow(H,4) - 180*G5pp*pow(G5x,2)*pow(H,4) + 1008*pow(G4x,2)*G5xx*pow(H,4) + 72*G4*G4xx*G5xx*pow(H,4) - 
      1440*G4x*G5p*G5xx*pow(H,4) + 612*pow(G5p,2)*G5xx*pow(H,4) - 240*G4*G5px*G5xx*pow(H,4) - 156*G4p*G5x*G5xx*pow(H,4) - 360*G4*G4x*G5xxx*pow(H,4) + 216*G4*G5p*G5xxx*pow(H,4) - 
      135*pow(G5x,3)*pow(H,6) - 9*pow(G3x,2)*(22*G4px - 2*G5pp - 27*G5x*pow(H,2)) + pow(G4px,2)*(-24*G5pp + 1644*G5x*pow(H,2)) + 
      3*G3x*(140*pow(G4px,2) + 12*G2xx*G4x - 56*G4ppx*G4x + 24*G2x*G4xx - 48*G3p*G4xx + 48*G4pp*G4xx - 6*G2xx*G5p + 28*G4ppx*G5p - 8*G4px*G5pp - 4*pow(G5pp,2) - 8*G4x*G5ppp + 4*G5p*G5ppp - 
         10*G2x*G5px + 20*G3p*G5px - 20*G4pp*G5px + 4*G2px*G5x - 2*G3pp*G5x - 4*G4ppp*G5x + 720*G4x*G4xx*pow(H,2) - 72*G4*G4xxx*pow(H,2) - 648*G4xx*G5p*pow(H,2) - 372*G4x*G5px*pow(H,2) + 
         348*G5p*G5px*pow(H,2) + 16*G4*G5pxx*pow(H,2) - 420*G4px*G5x*pow(H,2) + 6*G5pp*G5x*pow(H,2) - 52*G4p*G5xx*pow(H,2) + 99*pow(G5x,2)*pow(H,4)) - 
      6*G4px*(20*G2xx*G4x - 40*G4ppx*G4x + 16*G2x*G4xx - 32*G3p*G4xx + 32*G4pp*G4xx - 10*G2xx*G5p + 20*G4ppx*G5p - 4*pow(G5pp,2) - 8*G4x*G5ppp + 4*G5p*G5ppp - 6*G2x*G5px + 12*G3p*G5px - 12*G4pp*G5px + 
         4*G2px*G5x - 2*G3pp*G5x - 4*G4ppp*G5x + 888*G4x*G4xx*pow(H,2) - 120*G4*G4xxx*pow(H,2) - 876*G4xx*G5p*pow(H,2) - 500*G4x*G5px*pow(H,2) + 472*G5p*G5px*pow(H,2) + 48*G4*G5pxx*pow(H,2) + 
         18*G5pp*G5x*pow(H,2) - 100*G4p*G5xx*pow(H,2) + 141*pow(G5x,2)*pow(H,4)) + 36*G4xxx*G5x*pow(H,2)*rptot - 24*G5pxx*G5x*pow(H,2)*rptot + 72*G4xx*G5xx*pow(H,2)*rptot - 
      30*G5px*G5xx*pow(H,2)*rptot + 12*G4x*G5xxx*pow(H,2)*rptot - 6*G5p*G5xxx*pow(H,2)*rptot) - 
   3*G3x*(24*G3xx*G4pp*G4px - 64*G2xx*pow(G4px,2) + 264*G4ppx*pow(G4px,2) + 48*G4pp*G4px*G4pxx + 12*G3pp*G3xx*G4x + 24*G3pxx*G4pp*G4x - 24*G3xx*G4ppp*G4x + 40*G2xx*G4ppx*G4x - 96*pow(G4ppx,2)*G4x - 
      48*G4pp*G4ppxx*G4x - 64*G2pxx*G4px*G4x + 16*G3ppx*G4px*G4x + 96*G4pppx*G4px*G4x - 32*G2px*G4pxx*G4x + 8*G3pp*G4pxx*G4x + 48*G4ppp*G4pxx*G4x + 96*G4pp*G4ppx*G4xx + 32*G2px*G4px*G4xx + 
      112*G3pp*G4px*G4xx - 288*G4ppp*G4px*G4xx + 32*G2ppx*G4x*G4xx - 32*G3ppp*G4x*G4xx - 6*G3pp*G3xx*G5p - 12*G3pxx*G4pp*G5p + 12*G3xx*G4ppp*G5p - 20*G2xx*G4ppx*G5p + 48*pow(G4ppx,2)*G5p + 
      24*G4pp*G4ppxx*G5p + 32*G2pxx*G4px*G5p - 8*G3ppx*G4px*G5p - 48*G4pppx*G4px*G5p + 16*G2px*G4pxx*G5p - 4*G3pp*G4pxx*G5p - 24*G4ppp*G4pxx*G5p - 16*G2ppx*G4xx*G5p + 16*G3ppp*G4xx*G5p - 
      12*G3xx*G4pp*G5pp + 40*G2xx*G4px*G5pp - 96*G4ppx*G4px*G5pp + 24*G4pp*G4pxx*G5pp + 8*G2pxx*G4x*G5pp - 8*G3ppx*G4x*G5pp - 16*G2px*G4xx*G5pp + 16*G3pp*G4xx*G5pp - 4*G2pxx*G5p*G5pp + 4*G3ppx*G5p*G5pp - 
      4*G2xx*pow(G5pp,2) + 72*pow(G4px,2)*G5ppp - 8*G2xx*G4x*G5ppp + 4*G2xx*G5p*G5ppp - 48*G4pp*G4px*G5ppx + 16*G2px*G4x*G5ppx - 16*G3pp*G4x*G5ppx - 8*G2px*G5p*G5ppx + 8*G3pp*G5p*G5ppx + 
      4*G2xx*G4pp*G5px - 48*G4pp*G4ppx*G5px - 24*G2px*G4px*G5px - 48*G3pp*G4px*G5px + 144*G4ppp*G4px*G5px - 16*G2ppx*G4x*G5px + 16*G3ppp*G4x*G5px + 8*G2ppx*G5p*G5px - 8*G3ppp*G5p*G5px + 
      8*G2px*G5pp*G5px - 8*G3pp*G5pp*G5px + 2*G2xx*G3pp*G5x + 4*G2pxx*G4pp*G5x - 4*G3ppx*G4pp*G5x - 4*G2xx*G4ppp*G5x + 8*G2px*G4ppx*G5x - 8*G3pp*G4ppx*G5x - 8*G2ppx*G4px*G5x + 8*G3ppp*G4px*G5x + 
      144*G3xx*G4*G4pxx*pow(H,2) - 96*G4*pow(G4pxx,2)*pow(H,2) + 384*G4*G4px*G4pxxx*pow(H,2) + 144*G3xx*G4px*G4x*pow(H,2) + 2688*G4px*G4pxx*G4x*pow(H,2) - 192*G4p*G4pxxx*G4x*pow(H,2) + 
      360*G3pxx*pow(G4x,2)*pow(H,2) - 720*G4ppxx*pow(G4x,2)*pow(H,2) - 216*G3pxx*G4*G4xx*pow(H,2) - 216*G3xx*G4p*G4xx*pow(H,2) + 240*G4*G4ppxx*G4xx*pow(H,2) - 
      2832*pow(G4px,2)*G4xx*pow(H,2) + 1296*G4p*G4pxx*G4xx*pow(H,2) + 624*G4ppx*G4x*G4xx*pow(H,2) + 1152*G4pp*pow(G4xx,2)*pow(H,2) - 240*G4*G4ppx*G4xxx*pow(H,2) + 
      192*G4p*G4px*G4xxx*pow(H,2) + 48*G4pp*G4x*G4xxx*pow(H,2) - 180*G3xx*G4px*G5p*pow(H,2) - 2136*G4px*G4pxx*G5p*pow(H,2) + 96*G4p*G4pxxx*G5p*pow(H,2) - 468*G3pxx*G4x*G5p*pow(H,2) + 
      840*G4ppxx*G4x*G5p*pow(H,2) - 312*G4ppx*G4xx*G5p*pow(H,2) - 24*G4pp*G4xxx*G5p*pow(H,2) + 144*G3pxx*pow(G5p,2)*pow(H,2) - 240*G4ppxx*pow(G5p,2)*pow(H,2) - 
      48*G4*G4pxxx*G5pp*pow(H,2) + 600*G4px*G4xx*G5pp*pow(H,2) + 48*G4p*G4xxx*G5pp*pow(H,2) + 54*G3xx*G5p*G5pp*pow(H,2) + 36*G4pxx*G5p*G5pp*pow(H,2) + 48*G4xx*pow(G5pp,2)*pow(H,2) + 
      96*G4x*G4xx*G5ppp*pow(H,2) + 48*G4*G4xxx*G5ppp*pow(H,2) - 336*G4xx*G5p*G5ppp*pow(H,2) + 96*G4*G4xx*G5pppx*pow(H,2) + 48*G4x*G5p*G5pppx*pow(H,2) - 24*pow(G5p,2)*G5pppx*pow(H,2) - 
      84*G3xx*G4*G5ppx*pow(H,2) - 24*G4*G4pxx*G5ppx*pow(H,2) - 1304*G4px*G4x*G5ppx*pow(H,2) - 288*G4p*G4xx*G5ppx*pow(H,2) + 1084*G4px*G5p*G5ppx*pow(H,2) - 24*G4x*G5pp*G5ppx*pow(H,2) - 
      60*G5p*G5pp*G5ppx*pow(H,2) + 48*G4*pow(G5ppx,2)*pow(H,2) - 176*G4*G4px*G5ppxx*pow(H,2) + 80*G4p*G4x*G5ppxx*pow(H,2) - 40*G4p*G5p*G5ppxx*pow(H,2) + 24*G4*G5pp*G5ppxx*pow(H,2) + 
      120*G3pxx*G4*G5px*pow(H,2) + 72*G3xx*G4p*G5px*pow(H,2) - 144*G4*G4ppxx*G5px*pow(H,2) + 2140*pow(G4px,2)*G5px*pow(H,2) - 672*G4p*G4pxx*G5px*pow(H,2) + 8*G2xx*G4x*G5px*pow(H,2) - 
      376*G4ppx*G4x*G5px*pow(H,2) - 1152*G4pp*G4xx*G5px*pow(H,2) - 22*G2xx*G5p*G5px*pow(H,2) + 260*G4ppx*G5p*G5px*pow(H,2) - 440*G4px*G5pp*G5px*pow(H,2) - 12*pow(G5pp,2)*G5px*pow(H,2) - 
      24*G4x*G5ppp*G5px*pow(H,2) + 156*G5p*G5ppp*G5px*pow(H,2) - 48*G4*G5pppx*G5px*pow(H,2) + 192*G4p*G5ppx*G5px*pow(H,2) + 276*G4pp*pow(G5px,2)*pow(H,2) + 12*G2xx*G4*G5pxx*pow(H,2) + 
      56*G4*G4ppx*G5pxx*pow(H,2) - 312*G4p*G4px*G5pxx*pow(H,2) - 48*G4pp*G4x*G5pxx*pow(H,2) + 48*G4pp*G5p*G5pxx*pow(H,2) + 8*G4p*G5pp*G5pxx*pow(H,2) - 24*G4*G5ppp*G5pxx*pow(H,2) - 
      8*G4*G4pp*G5pxxx*pow(H,2) - 48*G3pxx*G4p*G5x*pow(H,2) + 12*G3xx*G4pp*G5x*pow(H,2) + 48*G4p*G4ppxx*G5x*pow(H,2) + 36*G2xx*G4px*G5x*pow(H,2) - 528*G4ppx*G4px*G5x*pow(H,2) - 
      120*G4pp*G4pxx*G5x*pow(H,2) + 132*G2pxx*G4x*G5x*pow(H,2) - 108*G3ppx*G4x*G5x*pow(H,2) - 48*G4pppx*G4x*G5x*pow(H,2) + 84*G2px*G4xx*G5x*pow(H,2) - 288*G3pp*G4xx*G5x*pow(H,2) + 
      240*G4ppp*G4xx*G5x*pow(H,2) - 12*G2p*G4xxx*G5x*pow(H,2) - 84*G2pxx*G5p*G5x*pow(H,2) + 48*G3ppx*G5p*G5x*pow(H,2) + 72*G4pppx*G5p*G5x*pow(H,2) - 6*G2xx*G5pp*G5x*pow(H,2) - 
      72*G4px*G5ppp*G5x*pow(H,2) + 24*G4p*G5pppx*G5x*pow(H,2) + 36*G4pp*G5ppx*G5x*pow(H,2) - 44*G2px*G5px*G5x*pow(H,2) + 142*G3pp*G5px*G5x*pow(H,2) - 108*G4ppp*G5px*G5x*pow(H,2) + 
      10*G2p*G5pxx*G5x*pow(H,2) + 4*G2ppx*pow(G5x,2)*pow(H,2) - 4*G3ppp*pow(G5x,2)*pow(H,2) - 16*G2pxx*G4*G5xx*pow(H,2) - 8*G3ppx*G4*G5xx*pow(H,2) - 16*G2xx*G4p*G5xx*pow(H,2) + 
      48*G4*G4pppx*G5xx*pow(H,2) + 72*G4p*G4ppx*G5xx*pow(H,2) - 304*G4pp*G4px*G5xx*pow(H,2) + 40*G2px*G4x*G5xx*pow(H,2) - 100*G3pp*G4x*G5xx*pow(H,2) + 40*G4ppp*G4x*G5xx*pow(H,2) + 
      72*G2p*G4xx*G5xx*pow(H,2) - 20*G2px*G5p*G5xx*pow(H,2) + 86*G3pp*G5p*G5xx*pow(H,2) - 92*G4ppp*G5p*G5xx*pow(H,2) + 20*G4pp*G5pp*G5xx*pow(H,2) - 72*G4p*G5ppp*G5xx*pow(H,2) - 
      40*G2p*G5px*G5xx*pow(H,2) - 4*G2pp*G5x*G5xx*pow(H,2) - 4*G3pp*G4*G5xxx*pow(H,2) - 8*G4p*G4pp*G5xxx*pow(H,2) + 8*G4*G4ppp*G5xxx*pow(H,2) - 4*G2p*G4x*G5xxx*pow(H,2) + 
      2*G2p*G5p*G5xxx*pow(H,2) - 984*G4x*G4xx*G5px*pow(H,4) + 168*G4*G4xxx*G5px*pow(H,4) + 240*G4xx*G5p*G5px*pow(H,4) + 512*G4x*pow(G5px,2)*pow(H,4) - 98*G5p*pow(G5px,2)*pow(H,4) - 
      96*pow(G4x,2)*G5pxx*pow(H,4) + 264*G4*G4xx*G5pxx*pow(H,4) + 408*G4x*G5p*G5pxx*pow(H,4) - 252*pow(G5p,2)*G5pxx*pow(H,4) - 208*G4*G5px*G5pxx*pow(H,4) - 168*G4*G4x*G5pxxx*pow(H,4) + 
      120*G4*G5p*G5pxxx*pow(H,4) - 576*G4*G4pxxx*G5x*pow(H,4) - 1296*G4pxx*G4x*G5x*pow(H,4) + 360*G4px*G4xx*G5x*pow(H,4) - 360*G4p*G4xxx*G5x*pow(H,4) + 1512*G4pxx*G5p*G5x*pow(H,4) + 
      1224*G4xx*G5pp*G5x*pow(H,4) + 208*G4x*G5ppx*G5x*pow(H,4) - 420*G5p*G5ppx*G5x*pow(H,4) + 296*G4*G5ppxx*G5x*pow(H,4) + 100*G4px*G5px*G5x*pow(H,4) - 706*G5pp*G5px*G5x*pow(H,4) + 
      344*G4p*G5pxx*G5x*pow(H,4) - 150*G4ppx*pow(G5x,2)*pow(H,4) + 54*G5ppp*pow(G5x,2)*pow(H,4) + 288*G4*G4pxx*G5xx*pow(H,4) + 144*G4px*G4x*G5xx*pow(H,4) - 1224*G4p*G4xx*G5xx*pow(H,4) - 
      276*G4px*G5p*G5xx*pow(H,4) + 480*G4x*G5pp*G5xx*pow(H,4) - 270*G5p*G5pp*G5xx*pow(H,4) - 164*G4*G5ppx*G5xx*pow(H,4) + 644*G4p*G5px*G5xx*pow(H,4) + 448*G4pp*G5x*G5xx*pow(H,4) + 
      24*G4*G4px*G5xxx*pow(H,4) - 96*G4p*G4x*G5xxx*pow(H,4) + 84*G4p*G5p*G5xxx*pow(H,4) - 36*G4*G5pp*G5xxx*pow(H,4) - 345*G5px*pow(G5x,2)*pow(H,6) + 
      2*G2x*(6*G3xx*G4px + 12*G4px*G4pxx + 6*G3pxx*G4x - 12*G4ppxx*G4x + 24*G4ppx*G4xx - 3*G3pxx*G5p + 6*G4ppxx*G5p - 3*G3xx*G5pp + 6*G4pxx*G5pp - 12*G4px*G5ppx + G2xx*G5px - 12*G4ppx*G5px + G2pxx*G5x - 
         G3ppx*G5x + 54*G4xx*G5px*pow(H,2) - 31*pow(G5px,2)*pow(H,2) - 28*G4x*G5pxx*pow(H,2) + 20*G5p*G5pxx*pow(H,2) - 2*G4*G5pxxx*pow(H,2) - 72*G4pxx*G5x*pow(H,2) + 
         34*G5ppx*G5x*pow(H,2) + 34*G4px*G5xx*pow(H,2) - 11*G5pp*G5xx*pow(H,2) - 2*G4p*G5xxx*pow(H,2)) - 
      4*G3p*(6*G3xx*G4px + 12*G4px*G4pxx + 6*G3pxx*G4x - 12*G4ppxx*G4x + 24*G4ppx*G4xx - 3*G3pxx*G5p + 6*G4ppxx*G5p - 3*G3xx*G5pp + 6*G4pxx*G5pp - 12*G4px*G5ppx + G2xx*G5px - 12*G4ppx*G5px + G2pxx*G5x - 
         G3ppx*G5x + 54*G4xx*G5px*pow(H,2) - 31*pow(G5px,2)*pow(H,2) - 28*G4x*G5pxx*pow(H,2) + 20*G5p*G5pxx*pow(H,2) - 2*G4*G5pxxx*pow(H,2) - 72*G4pxx*G5x*pow(H,2) + 
         34*G5ppx*G5x*pow(H,2) + 34*G4px*G5xx*pow(H,2) - 11*G5pp*G5xx*pow(H,2) - 2*G4p*G5xxx*pow(H,2)) + 12*G4xxx*G5px*pow(H,2)*rptot - 24*G4xx*G5pxx*pow(H,2)*rptot + 
      6*G5px*G5pxx*pow(H,2)*rptot + 4*G4x*G5pxxx*pow(H,2)*rptot - 2*G5p*G5pxxx*pow(H,2)*rptot + 12*G4pxxx*G5x*pow(H,2)*rptot - 6*G5ppxx*G5x*pow(H,2)*rptot - 24*G4pxx*G5xx*pow(H,2)*rptot + 
      12*G5ppx*G5xx*pow(H,2)*rptot + 4*G4px*G5xxx*pow(H,2)*rptot - 2*G5pp*G5xxx*pow(H,2)*rptot)) + 
pow(a,7)*pow(dphi,10)*pow(H,2)*(-288*pow(G4px,3)*G4pxx - 720*G3pxx*pow(G4px,2)*G4x + 1152*G4ppxx*pow(G4px,2)*G4x - 480*G2xx*G4px*G4pxx*G4x + 336*G3px*G4px*G4pxx*G4x - 
   288*G4ppx*G4px*G4pxx*G4x - 288*G2x*pow(G4pxx,2)*G4x + 576*G3p*pow(G4pxx,2)*G4x - 576*G4pp*pow(G4pxx,2)*G4x + 144*G2x*G4px*G4pxxx*G4x - 288*G3p*G4px*G4pxxx*G4x + 288*G4pp*G4px*G4pxxx*G4x - 
   144*G3px*G3pxx*pow(G4x,2) + 288*G3pxx*G4ppx*pow(G4x,2) + 96*G2xx*G4ppxx*pow(G4x,2) + 192*G3px*G4ppxx*pow(G4x,2) - 576*G4ppx*G4ppxx*pow(G4x,2) - 96*G2pxx*G4pxx*pow(G4x,2) - 
   192*G3ppx*G4pxx*pow(G4x,2) + 576*G4pppx*G4pxx*pow(G4x,2) - 96*G2px*G4pxxx*pow(G4x,2) + 96*G3pp*G4pxxx*pow(G4x,2) + 624*G2xx*pow(G4px,2)*G4xx + 1608*G3px*pow(G4px,2)*G4xx - 
   2736*G4ppx*pow(G4px,2)*G4xx - 288*G2x*G4px*G4pxx*G4xx + 576*G3p*G4px*G4pxx*G4xx - 576*G4pp*G4px*G4pxx*G4xx + 72*G2xx*G3px*G4x*G4xx + 216*pow(G3px,2)*G4x*G4xx - 144*G2x*G3pxx*G4x*G4xx + 
   288*G3p*G3pxx*G4x*G4xx - 288*G3pxx*G4pp*G4x*G4xx - 336*G2xx*G4ppx*G4x*G4xx - 816*G3px*G4ppx*G4x*G4xx + 1152*pow(G4ppx,2)*G4x*G4xx + 288*G2x*G4ppxx*G4x*G4xx - 576*G3p*G4ppxx*G4x*G4xx + 
   576*G4pp*G4ppxx*G4x*G4xx + 624*G2pxx*G4px*G4x*G4xx - 48*G3ppx*G4px*G4x*G4xx - 1152*G4pppx*G4px*G4x*G4xx + 384*G2px*G4pxx*G4x*G4xx - 96*G3pp*G4pxx*G4x*G4xx - 576*G4ppp*G4pxx*G4x*G4xx + 
   144*G2x*G3px*pow(G4xx,2) - 288*G3p*G3px*pow(G4xx,2) + 288*G3px*G4pp*pow(G4xx,2) - 288*G2x*G4ppx*pow(G4xx,2) + 576*G3p*G4ppx*pow(G4xx,2) - 576*G4pp*G4ppx*pow(G4xx,2) - 
   192*G2px*G4px*pow(G4xx,2) - 672*G3pp*G4px*pow(G4xx,2) + 1728*G4ppp*G4px*pow(G4xx,2) - 192*G2ppx*G4x*pow(G4xx,2) + 192*G3ppp*G4x*pow(G4xx,2) + 144*G2x*pow(G4px,2)*G4xxx - 
   288*G3p*pow(G4px,2)*G4xxx + 288*G4pp*pow(G4px,2)*G4xxx + 72*G2x*G3px*G4x*G4xxx - 144*G3p*G3px*G4x*G4xxx + 144*G3px*G4pp*G4x*G4xxx - 144*G2x*G4ppx*G4x*G4xxx + 288*G3p*G4ppx*G4x*G4xxx - 
   288*G4pp*G4ppx*G4x*G4xxx + 144*G3pp*G4px*G4x*G4xxx - 288*G4ppp*G4px*G4x*G4xxx + 96*G2ppx*pow(G4x,2)*G4xxx - 96*G3ppp*pow(G4x,2)*G4xxx + 360*G3pxx*pow(G4px,2)*G5p - 576*G4ppxx*pow(G4px,2)*G5p + 
   240*G2xx*G4px*G4pxx*G5p - 168*G3px*G4px*G4pxx*G5p + 144*G4ppx*G4px*G4pxx*G5p + 144*G2x*pow(G4pxx,2)*G5p - 288*G3p*pow(G4pxx,2)*G5p + 288*G4pp*pow(G4pxx,2)*G5p - 72*G2x*G4px*G4pxxx*G5p + 
   144*G3p*G4px*G4pxxx*G5p - 144*G4pp*G4px*G4pxxx*G5p + 144*G3px*G3pxx*G4x*G5p - 288*G3pxx*G4ppx*G4x*G5p - 96*G2xx*G4ppxx*G4x*G5p - 192*G3px*G4ppxx*G4x*G5p + 576*G4ppx*G4ppxx*G4x*G5p + 
   96*G2pxx*G4pxx*G4x*G5p + 192*G3ppx*G4pxx*G4x*G5p - 576*G4pppx*G4pxx*G4x*G5p + 96*G2px*G4pxxx*G4x*G5p - 96*G3pp*G4pxxx*G4x*G5p - 36*G2xx*G3px*G4xx*G5p - 108*pow(G3px,2)*G4xx*G5p + 
   72*G2x*G3pxx*G4xx*G5p - 144*G3p*G3pxx*G4xx*G5p + 144*G3pxx*G4pp*G4xx*G5p + 168*G2xx*G4ppx*G4xx*G5p + 408*G3px*G4ppx*G4xx*G5p - 576*pow(G4ppx,2)*G4xx*G5p - 144*G2x*G4ppxx*G4xx*G5p + 
   288*G3p*G4ppxx*G4xx*G5p - 288*G4pp*G4ppxx*G4xx*G5p - 312*G2pxx*G4px*G4xx*G5p + 24*G3ppx*G4px*G4xx*G5p + 576*G4pppx*G4px*G4xx*G5p - 192*G2px*G4pxx*G4xx*G5p + 48*G3pp*G4pxx*G4xx*G5p + 
   288*G4ppp*G4pxx*G4xx*G5p + 96*G2ppx*pow(G4xx,2)*G5p - 96*G3ppp*pow(G4xx,2)*G5p - 36*G2x*G3px*G4xxx*G5p + 72*G3p*G3px*G4xxx*G5p - 72*G3px*G4pp*G4xxx*G5p + 72*G2x*G4ppx*G4xxx*G5p - 
   144*G3p*G4ppx*G4xxx*G5p + 144*G4pp*G4ppx*G4xxx*G5p - 72*G3pp*G4px*G4xxx*G5p + 144*G4ppp*G4px*G4xxx*G5p - 96*G2ppx*G4x*G4xxx*G5p + 96*G3ppp*G4x*G4xxx*G5p - 36*G3px*G3pxx*pow(G5p,2) + 
   72*G3pxx*G4ppx*pow(G5p,2) + 24*G2xx*G4ppxx*pow(G5p,2) + 48*G3px*G4ppxx*pow(G5p,2) - 144*G4ppx*G4ppxx*pow(G5p,2) - 24*G2pxx*G4pxx*pow(G5p,2) - 48*G3ppx*G4pxx*pow(G5p,2) + 
   144*G4pppx*G4pxx*pow(G5p,2) - 24*G2px*G4pxxx*pow(G5p,2) + 24*G3pp*G4pxxx*pow(G5p,2) + 24*G2ppx*G4xxx*pow(G5p,2) - 24*G3ppp*G4xxx*pow(G5p,2) - 576*pow(G4px,2)*G4pxx*G5pp + 
   144*G3pxx*G4px*G4x*G5pp - 288*G4ppxx*G4px*G4x*G5pp + 96*G2xx*G4pxx*G4x*G5pp - 240*G3px*G4pxx*G4x*G5pp + 288*G4ppx*G4pxx*G4x*G5pp - 408*G2xx*G4px*G4xx*G5pp - 168*G3px*G4px*G4xx*G5pp + 
   1152*G4ppx*G4px*G4xx*G5pp - 144*G2x*G4pxx*G4xx*G5pp + 288*G3p*G4pxx*G4xx*G5pp - 288*G4pp*G4pxx*G4xx*G5pp - 96*G2pxx*G4x*G4xx*G5pp + 96*G3ppx*G4x*G4xx*G5pp + 96*G2px*pow(G4xx,2)*G5pp - 
   96*G3pp*pow(G4xx,2)*G5pp - 72*G2x*G4px*G4xxx*G5pp + 144*G3p*G4px*G4xxx*G5pp - 144*G4pp*G4px*G4xxx*G5pp - 72*G3pxx*G4px*G5p*G5pp + 144*G4ppxx*G4px*G5p*G5pp - 48*G2xx*G4pxx*G5p*G5pp + 
   120*G3px*G4pxx*G5p*G5pp - 144*G4ppx*G4pxx*G5p*G5pp + 48*G2pxx*G4xx*G5p*G5pp - 48*G3ppx*G4xx*G5p*G5pp + 144*G4px*G4pxx*pow(G5pp,2) + 48*G2xx*G4xx*pow(G5pp,2) - 48*G3px*G4xx*pow(G5pp,2) + 
   288*G4px*G4pxx*G4x*G5ppp - 864*pow(G4px,2)*G4xx*G5ppp + 96*G2xx*G4x*G4xx*G5ppp - 96*G3px*G4x*G4xx*G5ppp - 144*G4px*G4pxx*G5p*G5ppp - 48*G2xx*G4xx*G5p*G5ppp + 48*G3px*G4xx*G5p*G5ppp + 
   144*pow(G4px,2)*G4x*G5pppx - 48*G2xx*pow(G4x,2)*G5pppx + 48*G3px*pow(G4x,2)*G5pppx - 72*pow(G4px,2)*G5p*G5pppx + 48*G2xx*G4x*G5p*G5pppx - 48*G3px*G4x*G5p*G5pppx - 12*G2xx*pow(G5p,2)*G5pppx + 
   12*G3px*pow(G5p,2)*G5pppx + 9*pow(G3x,3)*(12*G4pxx - 5*G5ppx) + 792*pow(G4px,3)*G5ppx + 264*G2xx*G4px*G4x*G5ppx + 168*G3px*G4px*G4x*G5ppx - 864*G4ppx*G4px*G4x*G5ppx + 144*G2x*G4pxx*G4x*G5ppx - 
   288*G3p*G4pxx*G4x*G5ppx + 288*G4pp*G4pxx*G4x*G5ppx + 48*G2pxx*pow(G4x,2)*G5ppx - 48*G3ppx*pow(G4x,2)*G5ppx + 288*G2x*G4px*G4xx*G5ppx - 576*G3p*G4px*G4xx*G5ppx + 576*G4pp*G4px*G4xx*G5ppx - 
   192*G2px*G4x*G4xx*G5ppx + 192*G3pp*G4x*G4xx*G5ppx - 132*G2xx*G4px*G5p*G5ppx - 84*G3px*G4px*G5p*G5ppx + 432*G4ppx*G4px*G5p*G5ppx - 72*G2x*G4pxx*G5p*G5ppx + 144*G3p*G4pxx*G5p*G5ppx - 
   144*G4pp*G4pxx*G5p*G5ppx - 48*G2pxx*G4x*G5p*G5ppx + 48*G3ppx*G4x*G5p*G5ppx + 96*G2px*G4xx*G5p*G5ppx - 96*G3pp*G4xx*G5p*G5ppx + 12*G2pxx*pow(G5p,2)*G5ppx - 12*G3ppx*pow(G5p,2)*G5ppx - 
   216*pow(G4px,2)*G5pp*G5ppx - 48*G2xx*G4x*G5pp*G5ppx + 48*G3px*G4x*G5pp*G5ppx + 24*G2xx*G5p*G5pp*G5ppx - 24*G3px*G5p*G5pp*G5ppx - 72*G2x*G4px*G4x*G5ppxx + 144*G3p*G4px*G4x*G5ppxx - 
   144*G4pp*G4px*G4x*G5ppxx + 48*G2px*pow(G4x,2)*G5ppxx - 48*G3pp*pow(G4x,2)*G5ppxx + 36*G2x*G4px*G5p*G5ppxx - 72*G3p*G4px*G5p*G5ppxx + 72*G4pp*G4px*G5p*G5ppxx - 48*G2px*G4x*G5p*G5ppxx + 
   48*G3pp*G4x*G5p*G5ppxx + 12*G2px*pow(G5p,2)*G5ppxx - 12*G3pp*pow(G5p,2)*G5ppxx - 456*G2xx*pow(G4px,2)*G5px - 768*G3px*pow(G4px,2)*G5px + 1584*G4ppx*pow(G4px,2)*G5px + 
   72*G2x*G4px*G4pxx*G5px - 144*G3p*G4px*G4pxx*G5px + 144*G4pp*G4px*G4pxx*G5px - 48*G2xx*G3px*G4x*G5px - 96*pow(G3px,2)*G4x*G5px + 72*G2x*G3pxx*G4x*G5px - 144*G3p*G3pxx*G4x*G5px + 
   144*G3pxx*G4pp*G4x*G5px + 192*G2xx*G4ppx*G4x*G5px + 384*G3px*G4ppx*G4x*G5px - 576*pow(G4ppx,2)*G4x*G5px - 144*G2x*G4ppxx*G4x*G5px + 288*G3p*G4ppxx*G4x*G5px - 288*G4pp*G4ppxx*G4x*G5px - 
   336*G2pxx*G4px*G4x*G5px + 48*G3ppx*G4px*G4x*G5px + 576*G4pppx*G4px*G4x*G5px - 192*G2px*G4pxx*G4x*G5px + 48*G3pp*G4pxx*G4x*G5px + 288*G4ppp*G4pxx*G4x*G5px - 24*G2x*G2xx*G4xx*G5px + 
   48*G2xx*G3p*G4xx*G5px - 120*G2x*G3px*G4xx*G5px + 240*G3p*G3px*G4xx*G5px - 48*G2xx*G4pp*G4xx*G5px - 240*G3px*G4pp*G4xx*G5px + 288*G2x*G4ppx*G4xx*G5px - 576*G3p*G4ppx*G4xx*G5px + 
   576*G4pp*G4ppx*G4xx*G5px + 288*G2px*G4px*G4xx*G5px + 576*G3pp*G4px*G4xx*G5px - 1728*G4ppp*G4px*G4xx*G5px + 192*G2ppx*G4x*G4xx*G5px - 192*G3ppp*G4x*G4xx*G5px + 24*G2xx*G3px*G5p*G5px + 
   48*pow(G3px,2)*G5p*G5px - 36*G2x*G3pxx*G5p*G5px + 72*G3p*G3pxx*G5p*G5px - 72*G3pxx*G4pp*G5p*G5px - 96*G2xx*G4ppx*G5p*G5px - 192*G3px*G4ppx*G5p*G5px + 288*pow(G4ppx,2)*G5p*G5px + 
   72*G2x*G4ppxx*G5p*G5px - 144*G3p*G4ppxx*G5p*G5px + 144*G4pp*G4ppxx*G5p*G5px + 168*G2pxx*G4px*G5p*G5px - 24*G3ppx*G4px*G5p*G5px - 288*G4pppx*G4px*G5p*G5px + 96*G2px*G4pxx*G5p*G5px - 
   24*G3pp*G4pxx*G5p*G5px - 144*G4ppp*G4pxx*G5p*G5px - 96*G2ppx*G4xx*G5p*G5px + 96*G3ppp*G4xx*G5p*G5px + 240*G2xx*G4px*G5pp*G5px + 48*G3px*G4px*G5pp*G5px - 576*G4ppx*G4px*G5pp*G5px + 
   72*G2x*G4pxx*G5pp*G5px - 144*G3p*G4pxx*G5pp*G5px + 144*G4pp*G4pxx*G5pp*G5px + 48*G2pxx*G4x*G5pp*G5px - 48*G3ppx*G4x*G5pp*G5px - 96*G2px*G4xx*G5pp*G5px + 96*G3pp*G4xx*G5pp*G5px - 
   24*G2pxx*G5p*G5pp*G5px + 24*G3ppx*G5p*G5pp*G5px - 24*G2xx*pow(G5pp,2)*G5px + 24*G3px*pow(G5pp,2)*G5px + 432*pow(G4px,2)*G5ppp*G5px - 48*G2xx*G4x*G5ppp*G5px + 48*G3px*G4x*G5ppp*G5px + 
   24*G2xx*G5p*G5ppp*G5px - 24*G3px*G5p*G5ppp*G5px - 144*G2x*G4px*G5ppx*G5px + 288*G3p*G4px*G5ppx*G5px - 288*G4pp*G4px*G5ppx*G5px + 96*G2px*G4x*G5ppx*G5px - 96*G3pp*G4x*G5ppx*G5px - 
   48*G2px*G5p*G5ppx*G5px + 48*G3pp*G5p*G5ppx*G5px + 12*G2x*G2xx*pow(G5px,2) - 24*G2xx*G3p*pow(G5px,2) + 24*G2x*G3px*pow(G5px,2) - 48*G3p*G3px*pow(G5px,2) + 24*G2xx*G4pp*pow(G5px,2) + 
   48*G3px*G4pp*pow(G5px,2) - 72*G2x*G4ppx*pow(G5px,2) + 144*G3p*G4ppx*pow(G5px,2) - 144*G4pp*G4ppx*pow(G5px,2) - 96*G2px*G4px*pow(G5px,2) - 120*G3pp*G4px*pow(G5px,2) + 
   432*G4ppp*G4px*pow(G5px,2) - 48*G2ppx*G4x*pow(G5px,2) + 48*G3ppp*G4x*pow(G5px,2) + 24*G2ppx*G5p*pow(G5px,2) - 24*G3ppp*G5p*pow(G5px,2) + 24*G2px*G5pp*pow(G5px,2) - 
   24*G3pp*G5pp*pow(G5px,2) - 36*G2x*pow(G4px,2)*G5pxx + 72*G3p*pow(G4px,2)*G5pxx - 72*G4pp*pow(G4px,2)*G5pxx + 12*G2x*G2xx*G4x*G5pxx - 24*G2xx*G3p*G4x*G5pxx - 48*G2x*G3px*G4x*G5pxx + 
   96*G3p*G3px*G4x*G5pxx + 24*G2xx*G4pp*G4x*G5pxx - 96*G3px*G4pp*G4x*G5pxx + 72*G2x*G4ppx*G4x*G5pxx - 144*G3p*G4ppx*G4x*G5pxx + 144*G4pp*G4ppx*G4x*G5pxx - 48*G2px*G4px*G4x*G5pxx - 
   24*G3pp*G4px*G4x*G5pxx + 144*G4ppp*G4px*G4x*G5pxx - 48*G2ppx*pow(G4x,2)*G5pxx + 48*G3ppp*pow(G4x,2)*G5pxx - 6*G2x*G2xx*G5p*G5pxx + 12*G2xx*G3p*G5p*G5pxx + 24*G2x*G3px*G5p*G5pxx - 
   48*G3p*G3px*G5p*G5pxx - 12*G2xx*G4pp*G5p*G5pxx + 48*G3px*G4pp*G5p*G5pxx - 36*G2x*G4ppx*G5p*G5pxx + 72*G3p*G4ppx*G5p*G5pxx - 72*G4pp*G4ppx*G5p*G5pxx + 24*G2px*G4px*G5p*G5pxx + 12*G3pp*G4px*G5p*G5pxx - 
   72*G4ppp*G4px*G5p*G5pxx + 48*G2ppx*G4x*G5p*G5pxx - 48*G3ppp*G4x*G5p*G5pxx - 12*G2ppx*pow(G5p,2)*G5pxx + 12*G3ppp*pow(G5p,2)*G5pxx + 36*G2x*G4px*G5pp*G5pxx - 72*G3p*G4px*G5pp*G5pxx + 
   72*G4pp*G4px*G5pp*G5pxx - 60*G2xx*G3px*G4px*G5x - 12*pow(G3px,2)*G4px*G5x + 36*G2x*G3pxx*G4px*G5x - 72*G3p*G3pxx*G4px*G5x + 72*G3pxx*G4pp*G4px*G5x + 168*G2xx*G4ppx*G4px*G5x + 
   120*G3px*G4ppx*G4px*G5x - 288*pow(G4ppx,2)*G4px*G5x - 72*G2x*G4ppxx*G4px*G5x + 144*G3p*G4ppxx*G4px*G5x - 144*G4pp*G4ppxx*G4px*G5x - 120*G2pxx*pow(G4px,2)*G5x + 48*G3ppx*pow(G4px,2)*G5x + 
   144*G4pppx*pow(G4px,2)*G5x + 24*G2x*G2xx*G4pxx*G5x - 48*G2xx*G3p*G4pxx*G5x - 60*G2x*G3px*G4pxx*G5x + 120*G3p*G3px*G4pxx*G5x + 48*G2xx*G4pp*G4pxx*G5x - 120*G3px*G4pp*G4pxx*G5x + 
   72*G2x*G4ppx*G4pxx*G5x - 144*G3p*G4ppx*G4pxx*G5x + 144*G4pp*G4ppx*G4pxx*G5x - 96*G2px*G4px*G4pxx*G5x + 24*G3pp*G4px*G4pxx*G5x + 144*G4ppp*G4px*G4pxx*G5x + 48*G2xx*G3ppx*G4x*G5x - 
   48*G2pxx*G3px*G4x*G5x - 48*G2px*G3pxx*G4x*G5x + 48*G3pp*G3pxx*G4x*G5x - 96*G2xx*G4pppx*G4x*G5x + 96*G3px*G4pppx*G4x*G5x + 96*G2pxx*G4ppx*G4x*G5x - 96*G3ppx*G4ppx*G4x*G5x + 96*G2px*G4ppxx*G4x*G5x - 
   96*G3pp*G4ppxx*G4x*G5x - 96*G2ppx*G4pxx*G4x*G5x + 96*G3ppp*G4pxx*G4x*G5x - 24*G2pxx*G2x*G4xx*G5x + 48*G2pxx*G3p*G4xx*G5x - 24*G2xx*G3pp*G4xx*G5x + 24*G2x*G3ppx*G4xx*G5x - 48*G3p*G3ppx*G4xx*G5x + 
   48*G2px*G3px*G4xx*G5x - 24*G3pp*G3px*G4xx*G5x - 48*G2pxx*G4pp*G4xx*G5x + 48*G3ppx*G4pp*G4xx*G5x + 48*G2xx*G4ppp*G4xx*G5x - 48*G3px*G4ppp*G4xx*G5x - 96*G2px*G4ppx*G4xx*G5x + 96*G3pp*G4ppx*G4xx*G5x + 
   96*G2ppx*G4px*G4xx*G5x - 96*G3ppp*G4px*G4xx*G5x - 24*G2xx*G3ppx*G5p*G5x + 24*G2pxx*G3px*G5p*G5x + 24*G2px*G3pxx*G5p*G5x - 24*G3pp*G3pxx*G5p*G5x + 48*G2xx*G4pppx*G5p*G5x - 48*G3px*G4pppx*G5p*G5x - 
   48*G2pxx*G4ppx*G5p*G5x + 48*G3ppx*G4ppx*G5p*G5x - 48*G2px*G4ppxx*G5p*G5x + 48*G3pp*G4ppxx*G5p*G5x + 48*G2ppx*G4pxx*G5p*G5x - 48*G3ppp*G4pxx*G5p*G5x + 12*G2xx*G3px*G5pp*G5x - 
   12*pow(G3px,2)*G5pp*G5x - 24*G2xx*G4ppx*G5pp*G5x + 24*G3px*G4ppx*G5pp*G5x + 24*G2pxx*G4px*G5pp*G5x - 24*G3ppx*G4px*G5pp*G5x - 24*G2xx*G4px*G5ppp*G5x + 24*G3px*G4px*G5ppp*G5x - 
   12*G2x*G2xx*G5ppx*G5x + 24*G2xx*G3p*G5ppx*G5x + 12*G2x*G3px*G5ppx*G5x - 24*G3p*G3px*G5ppx*G5x - 24*G2xx*G4pp*G5ppx*G5x + 24*G3px*G4pp*G5ppx*G5x + 48*G2px*G4px*G5ppx*G5x - 48*G3pp*G4px*G5ppx*G5x + 
   12*G2pxx*G2x*G5px*G5x - 24*G2pxx*G3p*G5px*G5x + 12*G2xx*G3pp*G5px*G5x - 12*G2x*G3ppx*G5px*G5x + 24*G3p*G3ppx*G5px*G5x - 24*G2px*G3px*G5px*G5x + 12*G3pp*G3px*G5px*G5x + 24*G2pxx*G4pp*G5px*G5x - 
   24*G3ppx*G4pp*G5px*G5x - 24*G2xx*G4ppp*G5px*G5x + 24*G3px*G4ppp*G5px*G5x + 48*G2px*G4ppx*G5px*G5x - 48*G3pp*G4ppx*G5px*G5x - 48*G2ppx*G4px*G5px*G5x + 48*G3ppp*G4px*G5px*G5x - 
   4*G2px*G2pxx*pow(G5x,2) + 4*G2ppx*G2xx*pow(G5x,2) + 4*G2pxx*G3pp*pow(G5x,2) - 4*G2xx*G3ppp*pow(G5x,2) + 4*G2px*G3ppx*pow(G5x,2) - 4*G3pp*G3ppx*pow(G5x,2) - 4*G2ppx*G3px*pow(G5x,2) + 
   4*G3ppp*G3px*pow(G5x,2) - 12*G2x*G2xx*G4px*G5xx + 24*G2xx*G3p*G4px*G5xx - 24*G2x*G3px*G4px*G5xx + 48*G3p*G3px*G4px*G5xx - 24*G2xx*G4pp*G4px*G5xx - 48*G3px*G4pp*G4px*G5xx + 72*G2x*G4ppx*G4px*G5xx - 
   144*G3p*G4ppx*G4px*G5xx + 144*G4pp*G4ppx*G4px*G5xx + 48*G2px*pow(G4px,2)*G5xx + 60*G3pp*pow(G4px,2)*G5xx - 216*G4ppp*pow(G4px,2)*G5xx - 12*G2pxx*G2x*G4x*G5xx + 24*G2pxx*G3p*G4x*G5xx - 
   12*G2xx*G3pp*G4x*G5xx + 12*G2x*G3ppx*G4x*G5xx - 24*G3p*G3ppx*G4x*G5xx + 24*G2px*G3px*G4x*G5xx - 12*G3pp*G3px*G4x*G5xx - 24*G2pxx*G4pp*G4x*G5xx + 24*G3ppx*G4pp*G4x*G5xx + 24*G2xx*G4ppp*G4x*G5xx - 
   24*G3px*G4ppp*G4x*G5xx - 48*G2px*G4ppx*G4x*G5xx + 48*G3pp*G4ppx*G4x*G5xx + 48*G2ppx*G4px*G4x*G5xx - 48*G3ppp*G4px*G4x*G5xx + 6*G2pxx*G2x*G5p*G5xx - 12*G2pxx*G3p*G5p*G5xx + 6*G2xx*G3pp*G5p*G5xx - 
   6*G2x*G3ppx*G5p*G5xx + 12*G3p*G3ppx*G5p*G5xx - 12*G2px*G3px*G5p*G5xx + 6*G3pp*G3px*G5p*G5xx + 12*G2pxx*G4pp*G5p*G5xx - 12*G3ppx*G4pp*G5p*G5xx - 12*G2xx*G4ppp*G5p*G5xx + 12*G3px*G4ppp*G5p*G5xx + 
   24*G2px*G4ppx*G5p*G5xx - 24*G3pp*G4ppx*G5p*G5xx - 24*G2ppx*G4px*G5p*G5xx + 24*G3ppp*G4px*G5p*G5xx + 6*G2x*G2xx*G5pp*G5xx - 12*G2xx*G3p*G5pp*G5xx - 6*G2x*G3px*G5pp*G5xx + 12*G3p*G3px*G5pp*G5xx + 
   12*G2xx*G4pp*G5pp*G5xx - 12*G3px*G4pp*G5pp*G5xx - 24*G2px*G4px*G5pp*G5xx + 24*G3pp*G4px*G5pp*G5xx + 1152*G4*G4pxx*G4pxxx*G4x*pow(H,2) - 3456*pow(G4pxx,2)*pow(G4x,2)*pow(H,2) + 
   5472*G4px*G4pxxx*pow(G4x,2)*pow(H,2) - 576*G4*pow(G4pxx,2)*G4xx*pow(H,2) - 3744*G4*G4px*G4pxxx*G4xx*pow(H,2) - 23616*G4px*G4pxx*G4x*G4xx*pow(H,2) + 1440*G4p*G4pxxx*G4x*G4xx*pow(H,2) - 
   2592*G3pxx*pow(G4x,2)*G4xx*pow(H,2) + 5184*G4ppxx*pow(G4x,2)*G4xx*pow(H,2) + 864*G3pxx*G4*pow(G4xx,2)*pow(H,2) - 576*G4*G4ppxx*pow(G4xx,2)*pow(H,2) + 
   12672*pow(G4px,2)*pow(G4xx,2)*pow(H,2) - 5184*G4p*G4pxx*pow(G4xx,2)*pow(H,2) + 2592*G3px*G4x*pow(G4xx,2)*pow(H,2) - 1152*G4ppx*G4x*pow(G4xx,2)*pow(H,2) - 
   4032*G4pp*pow(G4xx,3)*pow(H,2) + 2880*G4*G4px*G4pxx*G4xxx*pow(H,2) - 1152*G4*G4ppxx*G4x*G4xxx*pow(H,2) + 1728*pow(G4px,2)*G4x*G4xxx*pow(H,2) - 576*G4p*G4pxx*G4x*G4xxx*pow(H,2) + 
   1296*G3px*pow(G4x,2)*G4xxx*pow(H,2) - 3744*G4ppx*pow(G4x,2)*G4xxx*pow(H,2) - 432*G3px*G4*G4xx*G4xxx*pow(H,2) + 2016*G4*G4ppx*G4xx*G4xxx*pow(H,2) - 2304*G4p*G4px*G4xx*G4xxx*pow(H,2) + 
   288*G4pp*G4x*G4xx*G4xxx*pow(H,2) - 576*G4*G4pxx*G4pxxx*G5p*pow(H,2) + 3168*pow(G4pxx,2)*G4x*G5p*pow(H,2) - 7344*G4px*G4pxxx*G4x*G5p*pow(H,2) + 20448*G4px*G4pxx*G4xx*G5p*pow(H,2) - 
   720*G4p*G4pxxx*G4xx*G5p*pow(H,2) + 3456*G3pxx*G4x*G4xx*G5p*pow(H,2) - 5760*G4ppxx*G4x*G4xx*G5p*pow(H,2) - 2592*G3px*pow(G4xx,2)*G5p*pow(H,2) - 288*G4ppx*pow(G4xx,2)*G5p*pow(H,2) + 
   576*G4*G4ppxx*G4xxx*G5p*pow(H,2) - 2736*pow(G4px,2)*G4xxx*G5p*pow(H,2) + 288*G4p*G4pxx*G4xxx*G5p*pow(H,2) - 1512*G3px*G4x*G4xxx*G5p*pow(H,2) + 4752*G4ppx*G4x*G4xxx*G5p*pow(H,2) - 
   144*G4pp*G4xx*G4xxx*G5p*pow(H,2) - 720*pow(G4pxx,2)*pow(G5p,2)*pow(H,2) + 2304*G4px*G4pxxx*pow(G5p,2)*pow(H,2) - 1080*G3pxx*G4xx*pow(G5p,2)*pow(H,2) + 
   1584*G4ppxx*G4xx*pow(G5p,2)*pow(H,2) + 432*G3px*G4xxx*pow(G5p,2)*pow(H,2) - 1440*G4ppx*G4xxx*pow(G5p,2)*pow(H,2) - 288*G4pxxx*pow(G4x,2)*G5pp*pow(H,2) + 
   576*G4*G4pxxx*G4xx*G5pp*pow(H,2) - 864*G4pxx*G4x*G4xx*G5pp*pow(H,2) - 2304*G4px*pow(G4xx,2)*G5pp*pow(H,2) - 576*G4*G4pxx*G4xxx*G5pp*pow(H,2) + 576*G4px*G4x*G4xxx*G5pp*pow(H,2) - 
   144*G4p*G4xx*G4xxx*G5pp*pow(H,2) + 576*G4pxxx*G4x*G5p*G5pp*pow(H,2) - 432*G4pxx*G4xx*G5p*G5pp*pow(H,2) + 936*G4px*G4xxx*G5p*G5pp*pow(H,2) - 216*G4pxxx*pow(G5p,2)*G5pp*pow(H,2) - 
   288*pow(G4xx,2)*pow(G5pp,2)*pow(H,2) - 288*G4x*G4xxx*pow(G5pp,2)*pow(H,2) - 576*G4x*pow(G4xx,2)*G5ppp*pow(H,2) + 288*pow(G4x,2)*G4xxx*G5ppp*pow(H,2) - 
   576*G4*G4xx*G4xxx*G5ppp*pow(H,2) + 2016*pow(G4xx,2)*G5p*G5ppp*pow(H,2) - 576*G4x*G4xxx*G5p*G5ppp*pow(H,2) + 216*G4xxx*pow(G5p,2)*G5ppp*pow(H,2) - 576*G4*pow(G4xx,2)*G5pppx*pow(H,2) + 
   576*G4*G4x*G4xxx*G5pppx*pow(H,2) - 576*G4x*G4xx*G5p*G5pppx*pow(H,2) - 288*G4*G4xxx*G5p*G5pppx*pow(H,2) + 288*G4xx*pow(G5p,2)*G5pppx*pow(H,2) - 576*G4*G4pxxx*G4x*G5ppx*pow(H,2) + 
   2112*G4pxx*pow(G4x,2)*G5ppx*pow(H,2) + 1152*G4*G4pxx*G4xx*G5ppx*pow(H,2) + 11760*G4px*G4x*G4xx*G5ppx*pow(H,2) + 864*G4p*pow(G4xx,2)*G5ppx*pow(H,2) - 1584*G4*G4px*G4xxx*G5ppx*pow(H,2) + 
   432*G4p*G4x*G4xxx*G5ppx*pow(H,2) + 288*G4*G4pxxx*G5p*G5ppx*pow(H,2) - 1536*G4pxx*G4x*G5p*G5ppx*pow(H,2) - 10200*G4px*G4xx*G5p*G5ppx*pow(H,2) - 216*G4p*G4xxx*G5p*G5ppx*pow(H,2) + 
   240*G4pxx*pow(G5p,2)*G5ppx*pow(H,2) + 288*G4x*G4xx*G5pp*G5ppx*pow(H,2) + 288*G4*G4xxx*G5pp*G5ppx*pow(H,2) + 720*G4xx*G5p*G5pp*G5ppx*pow(H,2) - 144*pow(G4x,2)*pow(G5ppx,2)*pow(H,2) - 
   576*G4*G4xx*pow(G5ppx,2)*pow(H,2) - 144*G4x*G5p*pow(G5ppx,2)*pow(H,2) + 108*pow(G5p,2)*pow(G5ppx,2)*pow(H,2) - 384*G4*G4pxx*G4x*G5ppxx*pow(H,2) - 
   2832*G4px*pow(G4x,2)*G5ppxx*pow(H,2) + 1680*G4*G4px*G4xx*G5ppxx*pow(H,2) - 528*G4p*G4x*G4xx*G5ppxx*pow(H,2) + 192*G4*G4pxx*G5p*G5ppxx*pow(H,2) + 3672*G4px*G4x*G5p*G5ppxx*pow(H,2) + 
   264*G4p*G4xx*G5p*G5ppxx*pow(H,2) - 1128*G4px*pow(G5p,2)*G5ppxx*pow(H,2) + 144*pow(G4x,2)*G5pp*G5ppxx*pow(H,2) - 288*G4*G4xx*G5pp*G5ppxx*pow(H,2) - 288*G4x*G5p*G5pp*G5ppxx*pow(H,2) + 
   108*pow(G5p,2)*G5pp*G5ppxx*pow(H,2) + 288*G4*G4x*G5ppx*G5ppxx*pow(H,2) - 144*G4*G5p*G5ppx*G5ppxx*pow(H,2) + 2016*G4*G4px*G4pxxx*G5px*pow(H,2) + 13440*G4px*G4pxx*G4x*G5px*pow(H,2) - 
   864*G4p*G4pxxx*G4x*G5px*pow(H,2) + 1344*G3pxx*pow(G4x,2)*G5px*pow(H,2) - 2976*G4ppxx*pow(G4x,2)*G5px*pow(H,2) - 1008*G3pxx*G4*G4xx*G5px*pow(H,2) + 864*G4*G4ppxx*G4xx*G5px*pow(H,2) - 
   19848*pow(G4px,2)*G4xx*G5px*pow(H,2) + 5904*G4p*G4pxx*G4xx*G5px*pow(H,2) - 168*G2xx*G4x*G4xx*G5px*pow(H,2) - 2448*G3px*G4x*G4xx*G5px*pow(H,2) + 1488*G4ppx*G4x*G4xx*G5px*pow(H,2) - 
   576*G2x*pow(G4xx,2)*G5px*pow(H,2) + 1152*G3p*pow(G4xx,2)*G5px*pow(H,2) + 6048*G4pp*pow(G4xx,2)*G5px*pow(H,2) + 288*G3px*G4*G4xxx*G5px*pow(H,2) - 1152*G4*G4ppx*G4xxx*G5px*pow(H,2) + 
   288*G4p*G4px*G4xxx*G5px*pow(H,2) + 72*G2x*G4x*G4xxx*G5px*pow(H,2) - 144*G3p*G4x*G4xxx*G5px*pow(H,2) + 144*G4pp*G4x*G4xxx*G5px*pow(H,2) - 10536*G4px*G4pxx*G5p*G5px*pow(H,2) + 
   432*G4p*G4pxxx*G5p*G5px*pow(H,2) - 1848*G3pxx*G4x*G5p*G5px*pow(H,2) + 3408*G4ppxx*G4x*G5p*G5px*pow(H,2) + 228*G2xx*G4xx*G5p*G5px*pow(H,2) + 2808*G3px*G4xx*G5p*G5px*pow(H,2) - 
   744*G4ppx*G4xx*G5p*G5px*pow(H,2) + 36*G2x*G4xxx*G5p*G5px*pow(H,2) - 72*G3p*G4xxx*G5p*G5px*pow(H,2) + 72*G4pp*G4xxx*G5p*G5px*pow(H,2) + 588*G3pxx*pow(G5p,2)*G5px*pow(H,2) - 
   960*G4ppxx*pow(G5p,2)*G5px*pow(H,2) - 288*G4*G4pxxx*G5pp*G5px*pow(H,2) - 144*G4pxx*G4x*G5pp*G5px*pow(H,2) + 3768*G4px*G4xx*G5pp*G5px*pow(H,2) + 288*G4p*G4xxx*G5pp*G5px*pow(H,2) + 
   288*G4pxx*G5p*G5pp*G5px*pow(H,2) + 144*G4xx*pow(G5pp,2)*G5px*pow(H,2) + 288*G4x*G4xx*G5ppp*G5px*pow(H,2) + 288*G4*G4xxx*G5ppp*G5px*pow(H,2) - 1872*G4xx*G5p*G5ppp*G5px*pow(H,2) + 
   144*pow(G4x,2)*G5pppx*G5px*pow(H,2) + 576*G4*G4xx*G5pppx*G5px*pow(H,2) + 144*G4x*G5p*G5pppx*G5px*pow(H,2) - 108*pow(G5p,2)*G5pppx*G5px*pow(H,2) - 432*G4*G4pxx*G5ppx*G5px*pow(H,2) - 
   6696*G4px*G4x*G5ppx*G5px*pow(H,2) - 1440*G4p*G4xx*G5ppx*G5px*pow(H,2) + 5796*G4px*G5p*G5ppx*G5px*pow(H,2) - 432*G5p*G5pp*G5ppx*G5px*pow(H,2) + 288*G4*pow(G5ppx,2)*G5px*pow(H,2) - 
   912*G4*G4px*G5ppxx*G5px*pow(H,2) + 336*G4p*G4x*G5ppxx*G5px*pow(H,2) - 168*G4p*G5p*G5ppxx*G5px*pow(H,2) + 144*G4*G5pp*G5ppxx*G5px*pow(H,2) + 288*G3pxx*G4*pow(G5px,2)*pow(H,2) - 
   288*G4*G4ppxx*pow(G5px,2)*pow(H,2) + 7248*pow(G4px,2)*pow(G5px,2)*pow(H,2) - 1584*G4p*G4pxx*pow(G5px,2)*pow(H,2) + 72*G2xx*G4x*pow(G5px,2)*pow(H,2) + 
   624*G3px*G4x*pow(G5px,2)*pow(H,2) - 528*G4ppx*G4x*pow(G5px,2)*pow(H,2) + 672*G2x*G4xx*pow(G5px,2)*pow(H,2) - 1344*G3p*G4xx*pow(G5px,2)*pow(H,2) - 
   2880*G4pp*G4xx*pow(G5px,2)*pow(H,2) - 120*G2xx*G5p*pow(G5px,2)*pow(H,2) - 768*G3px*G5p*pow(G5px,2)*pow(H,2) + 480*G4ppx*G5p*pow(G5px,2)*pow(H,2) - 
   1416*G4px*G5pp*pow(G5px,2)*pow(H,2) + 432*G5p*G5ppp*pow(G5px,2)*pow(H,2) - 144*G4*G5pppx*pow(G5px,2)*pow(H,2) + 504*G4p*G5ppx*pow(G5px,2)*pow(H,2) - 192*G2x*pow(G5px,3)*pow(H,2) + 
   384*G3p*pow(G5px,3)*pow(H,2) + 432*G4pp*pow(G5px,3)*pow(H,2) - 1680*G4*G4px*G4pxx*G5pxx*pow(H,2) + 96*G3pxx*G4*G4x*G5pxx*pow(H,2) + 384*G4*G4ppxx*G4x*G5pxx*pow(H,2) + 
   1992*pow(G4px,2)*G4x*G5pxx*pow(H,2) - 48*G4p*G4pxx*G4x*G5pxx*pow(H,2) + 120*G2xx*pow(G4x,2)*G5pxx*pow(H,2) - 816*G3px*pow(G4x,2)*G5pxx*pow(H,2) + 
   2064*G4ppx*pow(G4x,2)*G5pxx*pow(H,2) - 72*G2xx*G4*G4xx*G5pxx*pow(H,2) - 96*G3px*G4*G4xx*G5pxx*pow(H,2) - 240*G4*G4ppx*G4xx*G5pxx*pow(H,2) + 3312*G4p*G4px*G4xx*G5pxx*pow(H,2) + 
   600*G2x*G4x*G4xx*G5pxx*pow(H,2) - 1200*G3p*G4x*G4xx*G5pxx*pow(H,2) + 144*G4pp*G4x*G4xx*G5pxx*pow(H,2) - 72*G2x*G4*G4xxx*G5pxx*pow(H,2) + 144*G3p*G4*G4xxx*G5pxx*pow(H,2) - 
   144*G4*G4pp*G4xxx*G5pxx*pow(H,2) - 48*G3pxx*G4*G5p*G5pxx*pow(H,2) - 192*G4*G4ppxx*G5p*G5pxx*pow(H,2) - 984*pow(G4px,2)*G5p*G5pxx*pow(H,2) + 24*G4p*G4pxx*G5p*G5pxx*pow(H,2) - 
   156*G2xx*G4x*G5p*G5pxx*pow(H,2) + 768*G3px*G4x*G5p*G5pxx*pow(H,2) - 2184*G4ppx*G4x*G5p*G5pxx*pow(H,2) - 444*G2x*G4xx*G5p*G5pxx*pow(H,2) + 888*G3p*G4xx*G5p*G5pxx*pow(H,2) - 
   360*G4pp*G4xx*G5p*G5pxx*pow(H,2) + 48*G2xx*pow(G5p,2)*G5pxx*pow(H,2) - 180*G3px*pow(G5p,2)*G5pxx*pow(H,2) + 576*G4ppx*pow(G5p,2)*G5pxx*pow(H,2) + 432*G4*G4pxx*G5pp*G5pxx*pow(H,2) - 
   432*G4px*G4x*G5pp*G5pxx*pow(H,2) - 312*G4p*G4xx*G5pp*G5pxx*pow(H,2) - 204*G4px*G5p*G5pp*G5pxx*pow(H,2) + 144*G4x*pow(G5pp,2)*G5pxx*pow(H,2) - 144*pow(G4x,2)*G5ppp*G5pxx*pow(H,2) + 
   288*G4*G4xx*G5ppp*G5pxx*pow(H,2) + 288*G4x*G5p*G5ppp*G5pxx*pow(H,2) - 108*pow(G5p,2)*G5ppp*G5pxx*pow(H,2) - 288*G4*G4x*G5pppx*G5pxx*pow(H,2) + 144*G4*G5p*G5pppx*G5pxx*pow(H,2) + 
   552*G4*G4px*G5ppx*G5pxx*pow(H,2) + 24*G4p*G4x*G5ppx*G5pxx*pow(H,2) - 12*G4p*G5p*G5ppx*G5pxx*pow(H,2) - 144*G4*G5pp*G5ppx*G5pxx*pow(H,2) + 48*G2xx*G4*G5px*G5pxx*pow(H,2) + 
   192*G4*G4ppx*G5px*G5pxx*pow(H,2) - 1344*G4p*G4px*G5px*G5pxx*pow(H,2) - 384*G2x*G4x*G5px*G5pxx*pow(H,2) + 768*G3p*G4x*G5px*G5pxx*pow(H,2) - 288*G4pp*G4x*G5px*G5pxx*pow(H,2) + 
   228*G2x*G5p*G5px*G5pxx*pow(H,2) - 456*G3p*G5p*G5px*G5pxx*pow(H,2) + 216*G4pp*G5p*G5px*G5pxx*pow(H,2) + 48*G4p*G5pp*G5px*G5pxx*pow(H,2) - 144*G4*G5ppp*G5px*G5pxx*pow(H,2) + 
   36*G2x*G4*pow(G5pxx,2)*pow(H,2) - 72*G3p*G4*pow(G5pxx,2)*pow(H,2) + 72*G4*G4pp*pow(G5pxx,2)*pow(H,2) + 240*G4*pow(G4px,2)*G5pxxx*pow(H,2) + 96*G3px*G4*G4x*G5pxxx*pow(H,2) - 
   192*G4*G4ppx*G4x*G5pxxx*pow(H,2) - 288*G4p*G4px*G4x*G5pxxx*pow(H,2) - 96*G2x*pow(G4x,2)*G5pxxx*pow(H,2) + 192*G3p*pow(G4x,2)*G5pxxx*pow(H,2) - 96*G4pp*pow(G4x,2)*G5pxxx*pow(H,2) + 
   48*G2x*G4*G4xx*G5pxxx*pow(H,2) - 96*G3p*G4*G4xx*G5pxxx*pow(H,2) + 96*G4*G4pp*G4xx*G5pxxx*pow(H,2) - 48*G3px*G4*G5p*G5pxxx*pow(H,2) + 96*G4*G4ppx*G5p*G5pxxx*pow(H,2) + 
   144*G4p*G4px*G5p*G5pxxx*pow(H,2) + 120*G2x*G4x*G5p*G5pxxx*pow(H,2) - 240*G3p*G4x*G5p*G5pxxx*pow(H,2) + 144*G4pp*G4x*G5p*G5pxxx*pow(H,2) - 36*G2x*pow(G5p,2)*G5pxxx*pow(H,2) + 
   72*G3p*pow(G5p,2)*G5pxxx*pow(H,2) - 48*G4pp*pow(G5p,2)*G5pxxx*pow(H,2) - 48*G4*G4px*G5pp*G5pxxx*pow(H,2) + 48*G4p*G4x*G5pp*G5pxxx*pow(H,2) - 24*G4p*G5p*G5pp*G5pxxx*pow(H,2) - 
   24*G2x*G4*G5px*G5pxxx*pow(H,2) + 48*G3p*G4*G5px*G5pxxx*pow(H,2) - 48*G4*G4pp*G5px*G5pxxx*pow(H,2) + 288*G3pxx*G4*G4pxx*G5x*pow(H,2) + 6192*pow(G4px,2)*G4pxx*G5x*pow(H,2) - 
   288*G4p*pow(G4pxx,2)*G5x*pow(H,2) + 288*G3px*G4*G4pxxx*G5x*pow(H,2) - 576*G4*G4ppx*G4pxxx*G5x*pow(H,2) - 864*G4p*G4px*G4pxxx*G5x*pow(H,2) + 3096*G3pxx*G4px*G4x*G5x*pow(H,2) - 
   6480*G4ppxx*G4px*G4x*G5x*pow(H,2) + 576*G2xx*G4pxx*G4x*G5x*pow(H,2) - 1440*G3px*G4pxx*G4x*G5x*pow(H,2) + 3744*G4ppx*G4pxx*G4x*G5x*pow(H,2) - 648*G2x*G4pxxx*G4x*G5x*pow(H,2) + 
   1296*G3p*G4pxxx*G4x*G5x*pow(H,2) - 720*G4pp*G4pxxx*G4x*G5x*pow(H,2) + 360*G3pxx*G4p*G4xx*G5x*pow(H,2) - 144*G4p*G4ppxx*G4xx*G5x*pow(H,2) - 360*G2xx*G4px*G4xx*G5x*pow(H,2) - 
   3816*G3px*G4px*G4xx*G5x*pow(H,2) + 4176*G4ppx*G4px*G4xx*G5x*pow(H,2) + 1584*G2x*G4pxx*G4xx*G5x*pow(H,2) - 3168*G3p*G4pxx*G4xx*G5x*pow(H,2) + 1008*G4pp*G4pxx*G4xx*G5x*pow(H,2) - 
   936*G2pxx*G4x*G4xx*G5x*pow(H,2) + 648*G3ppx*G4x*G4xx*G5x*pow(H,2) + 576*G4pppx*G4x*G4xx*G5x*pow(H,2) - 432*G2px*pow(G4xx,2)*G5x*pow(H,2) + 1584*G3pp*pow(G4xx,2)*G5x*pow(H,2) - 
   1440*G4ppp*pow(G4xx,2)*G5x*pow(H,2) - 288*G3ppx*G4*G4xxx*G5x*pow(H,2) - 72*G3px*G4p*G4xxx*G5x*pow(H,2) + 576*G4*G4pppx*G4xxx*G5x*pow(H,2) + 432*G4p*G4ppx*G4xxx*G5x*pow(H,2) - 
   72*G2x*G4px*G4xxx*G5x*pow(H,2) + 144*G3p*G4px*G4xxx*G5x*pow(H,2) + 288*G4pp*G4px*G4xxx*G5x*pow(H,2) + 288*G2px*G4x*G4xxx*G5x*pow(H,2) - 936*G3pp*G4x*G4xxx*G5x*pow(H,2) + 
   720*G4ppp*G4x*G4xxx*G5x*pow(H,2) + 144*G2p*G4xx*G4xxx*G5x*pow(H,2) - 2016*G3pxx*G4px*G5p*G5x*pow(H,2) + 3888*G4ppxx*G4px*G5p*G5x*pow(H,2) - 360*G2xx*G4pxx*G5p*G5x*pow(H,2) + 
   468*G3px*G4pxx*G5p*G5x*pow(H,2) - 1512*G4ppx*G4pxx*G5p*G5x*pow(H,2) + 396*G2x*G4pxxx*G5p*G5x*pow(H,2) - 792*G3p*G4pxxx*G5p*G5x*pow(H,2) + 504*G4pp*G4pxxx*G5p*G5x*pow(H,2) + 
   612*G2pxx*G4xx*G5p*G5x*pow(H,2) - 180*G3ppx*G4xx*G5p*G5x*pow(H,2) - 864*G4pppx*G4xx*G5p*G5x*pow(H,2) - 144*G2px*G4xxx*G5p*G5x*pow(H,2) + 540*G3pp*G4xxx*G5p*G5x*pow(H,2) - 
   504*G4ppp*G4xxx*G5p*G5x*pow(H,2) - 504*G4px*G4pxx*G5pp*G5x*pow(H,2) + 144*G4p*G4pxxx*G5pp*G5x*pow(H,2) - 216*G3pxx*G4x*G5pp*G5x*pow(H,2) + 432*G4ppxx*G4x*G5pp*G5x*pow(H,2) + 
   108*G2xx*G4xx*G5pp*G5x*pow(H,2) + 36*G3px*G4xx*G5pp*G5x*pow(H,2) - 36*G2x*G4xxx*G5pp*G5x*pow(H,2) + 72*G3p*G4xxx*G5pp*G5x*pow(H,2) - 216*G4pp*G4xxx*G5pp*G5x*pow(H,2) + 
   180*G3pxx*G5p*G5pp*G5x*pow(H,2) - 360*G4ppxx*G5p*G5pp*G5x*pow(H,2) + 72*G4pxx*pow(G5pp,2)*G5x*pow(H,2) - 432*G4pxx*G4x*G5ppp*G5x*pow(H,2) + 864*G4px*G4xx*G5ppp*G5x*pow(H,2) - 
   144*G4p*G4xxx*G5ppp*G5x*pow(H,2) + 360*G4pxx*G5p*G5ppp*G5x*pow(H,2) - 288*G4*G4pxx*G5pppx*G5x*pow(H,2) + 144*G4px*G4x*G5pppx*G5x*pow(H,2) - 288*G4p*G4xx*G5pppx*G5x*pow(H,2) + 
   72*G4px*G5p*G5pppx*G5x*pow(H,2) - 144*G3pxx*G4*G5ppx*G5x*pow(H,2) + 288*G4*G4ppxx*G5ppx*G5x*pow(H,2) - 3468*pow(G4px,2)*G5ppx*G5x*pow(H,2) + 360*G4p*G4pxx*G5ppx*G5x*pow(H,2) - 
   316*G2xx*G4x*G5ppx*G5x*pow(H,2) + 244*G3px*G4x*G5ppx*G5x*pow(H,2) - 144*G4ppx*G4x*G5ppx*G5x*pow(H,2) - 744*G2x*G4xx*G5ppx*G5x*pow(H,2) + 1488*G3p*G4xx*G5ppx*G5x*pow(H,2) - 
   432*G4pp*G4xx*G5ppx*G5x*pow(H,2) + 206*G2xx*G5p*G5ppx*G5x*pow(H,2) + 46*G3px*G5p*G5ppx*G5x*pow(H,2) - 360*G4ppx*G5p*G5ppx*G5x*pow(H,2) + 288*G4px*G5pp*G5ppx*G5x*pow(H,2) - 
   144*G4p*pow(G5ppx,2)*G5x*pow(H,2) - 16*G2xx*G4*G5ppxx*G5x*pow(H,2) - 128*G3px*G4*G5ppxx*G5x*pow(H,2) + 288*G4*G4ppx*G5ppxx*G5x*pow(H,2) + 384*G4p*G4px*G5ppxx*G5x*pow(H,2) + 
   340*G2x*G4x*G5ppxx*G5x*pow(H,2) - 680*G3p*G4x*G5ppxx*G5x*pow(H,2) + 360*G4pp*G4x*G5ppxx*G5x*pow(H,2) - 206*G2x*G5p*G5ppxx*G5x*pow(H,2) + 412*G3p*G5p*G5ppxx*G5x*pow(H,2) - 
   252*G4pp*G5p*G5ppxx*G5x*pow(H,2) - 72*G4p*G5pp*G5ppxx*G5x*pow(H,2) - 216*G3pxx*G4p*G5px*G5x*pow(H,2) + 144*G4p*G4ppxx*G5px*G5x*pow(H,2) + 252*G2xx*G4px*G5px*G5x*pow(H,2) + 
   2112*G3px*G4px*G5px*G5x*pow(H,2) - 2712*G4ppx*G4px*G5px*G5x*pow(H,2) - 888*G2x*G4pxx*G5px*G5x*pow(H,2) + 1776*G3p*G4pxx*G5px*G5x*pow(H,2) - 720*G4pp*G4pxx*G5px*G5x*pow(H,2) + 
   496*G2pxx*G4x*G5px*G5x*pow(H,2) - 496*G3ppx*G4x*G5px*G5x*pow(H,2) + 456*G2px*G4xx*G5px*G5x*pow(H,2) - 1560*G3pp*G4xx*G5px*G5x*pow(H,2) + 1296*G4ppp*G4xx*G5px*G5x*pow(H,2) - 
   72*G2p*G4xxx*G5px*G5x*pow(H,2) - 332*G2pxx*G5p*G5px*G5x*pow(H,2) + 188*G3ppx*G5p*G5px*G5x*pow(H,2) + 288*G4pppx*G5p*G5px*G5x*pow(H,2) - 48*G2xx*G5pp*G5px*G5x*pow(H,2) - 
   60*G3px*G5pp*G5px*G5x*pow(H,2) + 72*G4ppx*G5pp*G5px*G5x*pow(H,2) - 360*G4px*G5ppp*G5px*G5x*pow(H,2) + 144*G4p*G5pppx*G5px*G5x*pow(H,2) + 420*G2x*G5ppx*G5px*G5x*pow(H,2) - 
   840*G3p*G5ppx*G5px*G5x*pow(H,2) + 288*G4pp*G5ppx*G5px*G5x*pow(H,2) - 120*G2px*pow(G5px,2)*G5x*pow(H,2) + 384*G3pp*pow(G5px,2)*G5x*pow(H,2) - 288*G4ppp*pow(G5px,2)*G5x*pow(H,2) + 
   16*G2pxx*G4*G5pxx*G5x*pow(H,2) + 128*G3ppx*G4*G5pxx*G5x*pow(H,2) - 12*G2xx*G4p*G5pxx*G5x*pow(H,2) - 48*G3px*G4p*G5pxx*G5x*pow(H,2) - 288*G4*G4pppx*G5pxx*G5x*pow(H,2) - 
   24*G4p*G4ppx*G5pxx*G5x*pow(H,2) - 228*G2x*G4px*G5pxx*G5x*pow(H,2) + 456*G3p*G4px*G5pxx*G5x*pow(H,2) - 360*G4pp*G4px*G5pxx*G5x*pow(H,2) - 152*G2px*G4x*G5pxx*G5x*pow(H,2) + 
   492*G3pp*G4x*G5pxx*G5x*pow(H,2) - 360*G4ppp*G4x*G5pxx*G5x*pow(H,2) - 120*G2p*G4xx*G5pxx*G5x*pow(H,2) + 52*G2px*G5p*G5pxx*G5x*pow(H,2) - 258*G3pp*G5p*G5pxx*G5x*pow(H,2) + 
   252*G4ppp*G5p*G5pxx*G5x*pow(H,2) + 30*G2x*G5pp*G5pxx*G5x*pow(H,2) - 60*G3p*G5pp*G5pxx*G5x*pow(H,2) + 108*G4pp*G5pp*G5pxx*G5x*pow(H,2) + 72*G4p*G5ppp*G5pxx*G5x*pow(H,2) + 
   60*G2p*G5px*G5pxx*G5x*pow(H,2) + 16*G2px*G4*G5pxxx*G5x*pow(H,2) - 16*G3pp*G4*G5pxxx*G5x*pow(H,2) + 12*G2x*G4p*G5pxxx*G5x*pow(H,2) - 24*G3p*G4p*G5pxxx*G5x*pow(H,2) + 
   24*G4p*G4pp*G5pxxx*G5x*pow(H,2) + 16*G2p*G4x*G5pxxx*G5x*pow(H,2) - 8*G2p*G5p*G5pxxx*G5x*pow(H,2) + 54*G2xx*G3px*pow(G5x,2)*pow(H,2) - 18*pow(G3px,2)*pow(G5x,2)*pow(H,2) - 
   90*G2x*G3pxx*pow(G5x,2)*pow(H,2) + 180*G3p*G3pxx*pow(G5x,2)*pow(H,2) - 108*G3pxx*G4pp*pow(G5x,2)*pow(H,2) - 180*G2xx*G4ppx*pow(G5x,2)*pow(H,2) + 
   108*G3px*G4ppx*pow(G5x,2)*pow(H,2) + 204*G2x*G4ppxx*pow(G5x,2)*pow(H,2) - 408*G3p*G4ppxx*pow(G5x,2)*pow(H,2) + 216*G4pp*G4ppxx*pow(G5x,2)*pow(H,2) + 
   288*G2pxx*G4px*pow(G5x,2)*pow(H,2) - 288*G3ppx*G4px*pow(G5x,2)*pow(H,2) - 48*G2px*G4pxx*pow(G5x,2)*pow(H,2) + 252*G3pp*G4pxx*pow(G5x,2)*pow(H,2) - 
   216*G4ppp*G4pxx*pow(G5x,2)*pow(H,2) + 24*G2p*G4pxxx*pow(G5x,2)*pow(H,2) - 48*G2ppx*G4xx*pow(G5x,2)*pow(H,2) + 48*G3ppp*G4xx*pow(G5x,2)*pow(H,2) - 
   24*G2pp*G4xxx*pow(G5x,2)*pow(H,2) - 24*G2pxx*G5pp*pow(G5x,2)*pow(H,2) + 24*G3ppx*G5pp*pow(G5x,2)*pow(H,2) + 24*G2xx*G5ppp*pow(G5x,2)*pow(H,2) - 24*G3px*G5ppp*pow(G5x,2)*pow(H,2) - 
   12*G2x*G5pppx*pow(G5x,2)*pow(H,2) + 24*G3p*G5pppx*pow(G5x,2)*pow(H,2) - 12*G3pp*G5ppx*pow(G5x,2)*pow(H,2) - 12*G2p*G5ppxx*pow(G5x,2)*pow(H,2) + 12*G2ppx*G5px*pow(G5x,2)*pow(H,2) - 
   12*G3ppp*G5px*pow(G5x,2)*pow(H,2) + 12*G2pp*G5pxx*pow(G5x,2)*pow(H,2) - 432*G3pxx*G4*G4px*G5xx*pow(H,2) + 576*G4*G4ppxx*G4px*G5xx*pow(H,2) - 3048*pow(G4px,3)*G5xx*pow(H,2) - 
   48*G2xx*G4*G4pxx*G5xx*pow(H,2) - 312*G3px*G4*G4pxx*G5xx*pow(H,2) + 432*G4*G4ppx*G4pxx*G5xx*pow(H,2) + 2592*G4p*G4px*G4pxx*G5xx*pow(H,2) + 72*G2x*G4*G4pxxx*G5xx*pow(H,2) - 
   144*G3p*G4*G4pxxx*G5xx*pow(H,2) + 144*G4*G4pp*G4pxxx*G5xx*pow(H,2) + 144*G3pxx*G4p*G4x*G5xx*pow(H,2) - 144*G2xx*G4px*G4x*G5xx*pow(H,2) - 1272*G3px*G4px*G4x*G5xx*pow(H,2) + 
   1776*G4ppx*G4px*G4x*G5xx*pow(H,2) + 672*G2x*G4pxx*G4x*G5xx*pow(H,2) - 1344*G3p*G4pxx*G4x*G5xx*pow(H,2) + 240*G4pp*G4pxx*G4x*G5xx*pow(H,2) - 168*G2pxx*pow(G4x,2)*G5xx*pow(H,2) + 
   216*G3ppx*pow(G4x,2)*G5xx*pow(H,2) - 96*G4pppx*pow(G4x,2)*G5xx*pow(H,2) + 120*G2pxx*G4*G4xx*G5xx*pow(H,2) + 168*G3ppx*G4*G4xx*G5xx*pow(H,2) + 120*G2xx*G4p*G4xx*G5xx*pow(H,2) - 
   336*G3px*G4p*G4xx*G5xx*pow(H,2) - 576*G4*G4pppx*G4xx*G5xx*pow(H,2) - 1296*G4p*G4ppx*G4xx*G5xx*pow(H,2) - 744*G2x*G4px*G4xx*G5xx*pow(H,2) + 1488*G3p*G4px*G4xx*G5xx*pow(H,2) + 
   3216*G4pp*G4px*G4xx*G5xx*pow(H,2) - 408*G2px*G4x*G4xx*G5xx*pow(H,2) + 1056*G3pp*G4x*G4xx*G5xx*pow(H,2) - 480*G4ppp*G4x*G4xx*G5xx*pow(H,2) - 432*G2p*pow(G4xx,2)*G5xx*pow(H,2) + 
   72*G3pp*G4*G4xxx*G5xx*pow(H,2) + 72*G2x*G4p*G4xxx*G5xx*pow(H,2) - 144*G3p*G4p*G4xxx*G5xx*pow(H,2) + 144*G4p*G4pp*G4xxx*G5xx*pow(H,2) - 144*G4*G4ppp*G4xxx*G5xx*pow(H,2) + 
   72*G2p*G4x*G4xxx*G5xx*pow(H,2) - 72*G3pxx*G4p*G5p*G5xx*pow(H,2) + 132*G2xx*G4px*G5p*G5xx*pow(H,2) + 1260*G3px*G4px*G5p*G5xx*pow(H,2) - 1392*G4ppx*G4px*G5p*G5xx*pow(H,2) - 
   480*G2x*G4pxx*G5p*G5xx*pow(H,2) + 960*G3p*G4pxx*G5p*G5xx*pow(H,2) - 408*G4pp*G4pxx*G5p*G5xx*pow(H,2) + 228*G2pxx*G4x*G5p*G5xx*pow(H,2) - 132*G3ppx*G4x*G5p*G5xx*pow(H,2) - 
   192*G4pppx*G4x*G5p*G5xx*pow(H,2) + 204*G2px*G4xx*G5p*G5xx*pow(H,2) - 960*G3pp*G4xx*G5p*G5xx*pow(H,2) + 1104*G4ppp*G4xx*G5p*G5xx*pow(H,2) - 36*G2p*G4xxx*G5p*G5xx*pow(H,2) - 
   72*G2pxx*pow(G5p,2)*G5xx*pow(H,2) + 12*G3ppx*pow(G5p,2)*G5xx*pow(H,2) + 120*G4pppx*pow(G5p,2)*G5xx*pow(H,2) + 72*G3pxx*G4*G5pp*G5xx*pow(H,2) - 144*G4*G4ppxx*G5pp*G5xx*pow(H,2) + 
   1068*pow(G4px,2)*G5pp*G5xx*pow(H,2) - 432*G4p*G4pxx*G5pp*G5xx*pow(H,2) + 48*G2xx*G4x*G5pp*G5xx*pow(H,2) - 96*G3px*G4x*G5pp*G5xx*pow(H,2) + 240*G4ppx*G4x*G5pp*G5xx*pow(H,2) + 
   228*G2x*G4xx*G5pp*G5xx*pow(H,2) - 456*G3p*G4xx*G5pp*G5xx*pow(H,2) - 240*G4pp*G4xx*G5pp*G5xx*pow(H,2) - 54*G2xx*G5p*G5pp*G5xx*pow(H,2) - 66*G3px*G5p*G5pp*G5xx*pow(H,2) + 
   168*G4ppx*G5p*G5pp*G5xx*pow(H,2) - 24*G4px*pow(G5pp,2)*G5xx*pow(H,2) + 144*G4*G4pxx*G5ppp*G5xx*pow(H,2) - 48*G4px*G4x*G5ppp*G5xx*pow(H,2) + 864*G4p*G4xx*G5ppp*G5xx*pow(H,2) - 
   408*G4px*G5p*G5ppp*G5xx*pow(H,2) + 144*G4*G4px*G5pppx*G5xx*pow(H,2) - 144*G4p*G4x*G5pppx*G5xx*pow(H,2) + 72*G4p*G5p*G5pppx*G5xx*pow(H,2) + 36*G2xx*G4*G5ppx*G5xx*pow(H,2) + 
   180*G3px*G4*G5ppx*G5xx*pow(H,2) - 432*G4*G4ppx*G5ppx*G5xx*pow(H,2) - 792*G4p*G4px*G5ppx*G5xx*pow(H,2) - 312*G2x*G4x*G5ppx*G5xx*pow(H,2) + 624*G3p*G4x*G5ppx*G5xx*pow(H,2) - 
   120*G4pp*G4x*G5ppx*G5xx*pow(H,2) + 228*G2x*G5p*G5ppx*G5xx*pow(H,2) - 456*G3p*G5p*G5ppx*G5xx*pow(H,2) + 204*G4pp*G5p*G5ppx*G5xx*pow(H,2) + 216*G4p*G5pp*G5ppx*G5xx*pow(H,2) - 
   36*G2x*G4*G5ppxx*G5xx*pow(H,2) + 72*G3p*G4*G5ppxx*G5xx*pow(H,2) - 72*G4*G4pp*G5ppxx*G5xx*pow(H,2) - 72*G2pxx*G4*G5px*G5xx*pow(H,2) - 72*G3ppx*G4*G5px*G5xx*pow(H,2) - 
   48*G2xx*G4p*G5px*G5xx*pow(H,2) + 264*G3px*G4p*G5px*G5xx*pow(H,2) + 288*G4*G4pppx*G5px*G5xx*pow(H,2) + 432*G4p*G4ppx*G5px*G5xx*pow(H,2) + 648*G2x*G4px*G5px*G5xx*pow(H,2) - 
   1296*G3p*G4px*G5px*G5xx*pow(H,2) - 1440*G4pp*G4px*G5px*G5xx*pow(H,2) + 240*G2px*G4x*G5px*G5xx*pow(H,2) - 540*G3pp*G4x*G5px*G5xx*pow(H,2) + 168*G4ppp*G4x*G5px*G5xx*pow(H,2) + 
   480*G2p*G4xx*G5px*G5xx*pow(H,2) - 96*G2px*G5p*G5px*G5xx*pow(H,2) + 462*G3pp*G5p*G5px*G5xx*pow(H,2) - 516*G4ppp*G5p*G5px*G5xx*pow(H,2) - 150*G2x*G5pp*G5px*G5xx*pow(H,2) + 
   300*G3p*G5pp*G5px*G5xx*pow(H,2) + 84*G4pp*G5pp*G5px*G5xx*pow(H,2) - 432*G4p*G5ppp*G5px*G5xx*pow(H,2) - 132*G2p*pow(G5px,2)*G5xx*pow(H,2) - 24*G2px*G4*G5pxx*G5xx*pow(H,2) - 
   12*G3pp*G4*G5pxx*G5xx*pow(H,2) - 72*G2x*G4p*G5pxx*G5xx*pow(H,2) + 144*G3p*G4p*G5pxx*G5xx*pow(H,2) - 144*G4p*G4pp*G5pxx*G5xx*pow(H,2) + 72*G4*G4ppp*G5pxx*G5xx*pow(H,2) - 
   60*G2p*G4x*G5pxx*G5xx*pow(H,2) + 30*G2p*G5p*G5pxx*G5xx*pow(H,2) + 102*G2x*G3px*G5x*G5xx*pow(H,2) - 204*G3p*G3px*G5x*G5xx*pow(H,2) + 24*G2pxx*G4p*G5x*G5xx*pow(H,2) + 
   48*G3ppx*G4p*G5x*G5xx*pow(H,2) + 12*G2xx*G4pp*G5x*G5xx*pow(H,2) + 36*G3px*G4pp*G5x*G5xx*pow(H,2) - 144*G4p*G4pppx*G5x*G5xx*pow(H,2) - 192*G2x*G4ppx*G5x*G5xx*pow(H,2) + 
   384*G3p*G4ppx*G5x*G5xx*pow(H,2) - 24*G4pp*G4ppx*G5x*G5xx*pow(H,2) + 120*G2px*G4px*G5x*G5xx*pow(H,2) - 396*G3pp*G4px*G5x*G5xx*pow(H,2) + 264*G4ppp*G4px*G5x*G5xx*pow(H,2) - 
   84*G2p*G4pxx*G5x*G5xx*pow(H,2) - 8*G2ppx*G4x*G5x*G5xx*pow(H,2) + 8*G3ppp*G4x*G5x*G5xx*pow(H,2) + 48*G2pp*G4xx*G5x*G5xx*pow(H,2) + 28*G2ppx*G5p*G5x*G5xx*pow(H,2) - 
   28*G3ppp*G5p*G5x*G5xx*pow(H,2) - 24*G2px*G5pp*G5x*G5xx*pow(H,2) + 36*G3pp*G5pp*G5x*G5xx*pow(H,2) + 12*G2x*G5ppp*G5x*G5xx*pow(H,2) - 24*G3p*G5ppp*G5x*G5xx*pow(H,2) + 
   24*G2p*G5ppx*G5x*G5xx*pow(H,2) - 24*G2pp*G5px*G5x*G5xx*pow(H,2) + 12*G2ppx*G4*pow(G5xx,2)*pow(H,2) - 12*G3ppp*G4*pow(G5xx,2)*pow(H,2) + 12*G2px*G4p*pow(G5xx,2)*pow(H,2) - 
   66*G3pp*G4p*pow(G5xx,2)*pow(H,2) - 18*G2x*G4pp*pow(G5xx,2)*pow(H,2) + 36*G3p*G4pp*pow(G5xx,2)*pow(H,2) - 36*pow(G4pp,2)*pow(G5xx,2)*pow(H,2) + 
   108*G4p*G4ppp*pow(G5xx,2)*pow(H,2) + 66*G2p*G4px*pow(G5xx,2)*pow(H,2) + 12*G2pp*G4x*pow(G5xx,2)*pow(H,2) - 6*G2pp*G5p*pow(G5xx,2)*pow(H,2) - 6*G2p*G5pp*pow(G5xx,2)*pow(H,2) + 
   120*G3px*G4*G4px*G5xxx*pow(H,2) - 336*G4*G4ppx*G4px*G5xxx*pow(H,2) - 48*G4p*pow(G4px,2)*G5xxx*pow(H,2) - 48*G2x*G4*G4pxx*G5xxx*pow(H,2) + 96*G3p*G4*G4pxx*G5xxx*pow(H,2) - 
   96*G4*G4pp*G4pxx*G5xxx*pow(H,2) - 96*G3ppx*G4*G4x*G5xxx*pow(H,2) - 24*G3px*G4p*G4x*G5xxx*pow(H,2) + 192*G4*G4pppx*G4x*G5xxx*pow(H,2) + 144*G4p*G4ppx*G4x*G5xxx*pow(H,2) + 
   144*G4pp*G4px*G4x*G5xxx*pow(H,2) + 48*G2px*pow(G4x,2)*G5xxx*pow(H,2) - 144*G3pp*pow(G4x,2)*G5xxx*pow(H,2) + 96*G4ppp*pow(G4x,2)*G5xxx*pow(H,2) + 48*G3pp*G4*G4xx*G5xxx*pow(H,2) + 
   48*G2x*G4p*G4xx*G5xxx*pow(H,2) - 96*G3p*G4p*G4xx*G5xxx*pow(H,2) + 96*G4p*G4pp*G4xx*G5xxx*pow(H,2) - 96*G4*G4ppp*G4xx*G5xxx*pow(H,2) + 48*G2p*G4x*G4xx*G5xxx*pow(H,2) + 
   48*G3ppx*G4*G5p*G5xxx*pow(H,2) + 12*G3px*G4p*G5p*G5xxx*pow(H,2) - 96*G4*G4pppx*G5p*G5xxx*pow(H,2) - 72*G4p*G4ppx*G5p*G5xxx*pow(H,2) + 24*G2x*G4px*G5p*G5xxx*pow(H,2) - 
   48*G3p*G4px*G5p*G5xxx*pow(H,2) - 24*G4pp*G4px*G5p*G5xxx*pow(H,2) - 48*G2px*G4x*G5p*G5xxx*pow(H,2) + 168*G3pp*G4x*G5p*G5xxx*pow(H,2) - 144*G4ppp*G4x*G5p*G5xxx*pow(H,2) - 
   24*G2p*G4xx*G5p*G5xxx*pow(H,2) + 12*G2px*pow(G5p,2)*G5xxx*pow(H,2) - 48*G3pp*pow(G5p,2)*G5xxx*pow(H,2) + 48*G4ppp*pow(G5p,2)*G5xxx*pow(H,2) - 24*G3px*G4*G5pp*G5xxx*pow(H,2) + 
   48*G4*G4ppx*G5pp*G5xxx*pow(H,2) + 144*G4p*G4px*G5pp*G5xxx*pow(H,2) - 24*G2x*G4x*G5pp*G5xxx*pow(H,2) + 48*G3p*G4x*G5pp*G5xxx*pow(H,2) - 96*G4pp*G4x*G5pp*G5xxx*pow(H,2) + 
   24*G4pp*G5p*G5pp*G5xxx*pow(H,2) - 24*G4p*pow(G5pp,2)*G5xxx*pow(H,2) + 48*G4*G4px*G5ppp*G5xxx*pow(H,2) - 48*G4p*G4x*G5ppp*G5xxx*pow(H,2) + 24*G4p*G5p*G5ppp*G5xxx*pow(H,2) + 
   24*G2x*G4*G5ppx*G5xxx*pow(H,2) - 48*G3p*G4*G5ppx*G5xxx*pow(H,2) + 48*G4*G4pp*G5ppx*G5xxx*pow(H,2) - 24*G3pp*G4*G5px*G5xxx*pow(H,2) - 12*G2x*G4p*G5px*G5xxx*pow(H,2) + 
   24*G3p*G4p*G5px*G5xxx*pow(H,2) - 24*G4p*G4pp*G5px*G5xxx*pow(H,2) + 48*G4*G4ppp*G5px*G5xxx*pow(H,2) - 24*G2p*G4x*G5px*G5xxx*pow(H,2) + 12*G2p*G5p*G5px*G5xxx*pow(H,2) - 
   16*G2ppx*G4*G5x*G5xxx*pow(H,2) + 16*G3ppp*G4*G5x*G5xxx*pow(H,2) + 12*G3pp*G4p*G5x*G5xxx*pow(H,2) - 12*G2x*G4pp*G5x*G5xxx*pow(H,2) + 24*G3p*G4pp*G5x*G5xxx*pow(H,2) - 
   24*pow(G4pp,2)*G5x*G5xxx*pow(H,2) - 24*G4p*G4ppp*G5x*G5xxx*pow(H,2) - 12*G2p*G4px*G5x*G5xxx*pow(H,2) - 16*G2pp*G4x*G5x*G5xxx*pow(H,2) + 8*G2pp*G5p*G5x*G5xxx*pow(H,2) + 
   2880*G4x*pow(G4xx,2)*G5px*pow(H,4) + 720*pow(G4x,2)*G4xxx*G5px*pow(H,4) - 720*G4*G4xx*G4xxx*G5px*pow(H,4) - 576*pow(G4xx,2)*G5p*G5px*pow(H,4) - 216*G4x*G4xxx*G5p*G5px*pow(H,4) - 
   288*G4xxx*pow(G5p,2)*G5px*pow(H,4) - 2904*G4x*G4xx*pow(G5px,2)*pow(H,4) + 576*G4*G4xxx*pow(G5px,2)*pow(H,4) + 420*G4xx*G5p*pow(G5px,2)*pow(H,4) + 768*G4x*pow(G5px,3)*pow(H,4) - 
   48*G5p*pow(G5px,3)*pow(H,4) + 2160*pow(G4x,2)*G4xx*G5pxx*pow(H,4) - 864*G4*pow(G4xx,2)*G5pxx*pow(H,4) - 1008*G4*G4x*G4xxx*G5pxx*pow(H,4) - 4968*G4x*G4xx*G5p*G5pxx*pow(H,4) + 
   720*G4*G4xxx*G5p*G5pxx*pow(H,4) + 2592*G4xx*pow(G5p,2)*G5pxx*pow(H,4) - 1616*pow(G4x,2)*G5px*G5pxx*pow(H,4) + 912*G4*G4xx*G5px*G5pxx*pow(H,4) + 2744*G4x*G5p*G5px*G5pxx*pow(H,4) - 
   1244*pow(G5p,2)*G5px*G5pxx*pow(H,4) - 384*G4*pow(G5px,2)*G5pxx*pow(H,4) + 632*G4*G4x*pow(G5pxx,2)*pow(H,4) - 376*G4*G5p*pow(G5pxx,2)*pow(H,4) - 576*pow(G4x,3)*G5pxxx*pow(H,4) + 
   1152*G4*G4x*G4xx*G5pxxx*pow(H,4) + 1296*pow(G4x,2)*G5p*G5pxxx*pow(H,4) - 864*G4*G4xx*G5p*G5pxxx*pow(H,4) - 936*G4x*pow(G5p,2)*G5pxxx*pow(H,4) + 216*pow(G5p,3)*G5pxxx*pow(H,4) - 
   560*G4*G4x*G5px*G5pxxx*pow(H,4) + 448*G4*G5p*G5px*G5pxxx*pow(H,4) - 16*pow(G4,2)*G5pxx*G5pxxx*pow(H,4) - 5616*G4pxxx*pow(G4x,2)*G5x*pow(H,4) + 3888*G4*G4pxxx*G4xx*G5x*pow(H,4) + 
   14688*G4pxx*G4x*G4xx*G5x*pow(H,4) - 1728*G4px*pow(G4xx,2)*G5x*pow(H,4) - 2592*G4*G4pxx*G4xxx*G5x*pow(H,4) - 864*G4px*G4x*G4xxx*G5x*pow(H,4) + 3024*G4p*G4xx*G4xxx*G5x*pow(H,4) + 
   8424*G4pxxx*G4x*G5p*G5x*pow(H,4) - 14688*G4pxx*G4xx*G5p*G5x*pow(H,4) + 1512*G4px*G4xxx*G5p*G5x*pow(H,4) - 3024*G4pxxx*pow(G5p,2)*G5x*pow(H,4) - 4752*pow(G4xx,2)*G5pp*G5x*pow(H,4) - 
   1296*G4x*G4xxx*G5pp*G5x*pow(H,4) + 324*G4xxx*G5p*G5pp*G5x*pow(H,4) - 4008*G4x*G4xx*G5ppx*G5x*pow(H,4) + 1320*G4*G4xxx*G5ppx*G5x*pow(H,4) + 4500*G4xx*G5p*G5ppx*G5x*pow(H,4) + 
   2616*pow(G4x,2)*G5ppxx*G5x*pow(H,4) - 2040*G4*G4xx*G5ppxx*G5x*pow(H,4) - 3828*G4x*G5p*G5ppxx*G5x*pow(H,4) + 1344*pow(G5p,2)*G5ppxx*G5x*pow(H,4) - 1968*G4*G4pxxx*G5px*G5x*pow(H,4) - 
   8376*G4pxx*G4x*G5px*G5x*pow(H,4) - 768*G4px*G4xx*G5px*G5x*pow(H,4) - 1368*G4p*G4xxx*G5px*G5x*pow(H,4) + 7740*G4pxx*G5p*G5px*G5x*pow(H,4) + 5772*G4xx*G5pp*G5px*G5x*pow(H,4) + 
   2340*G4x*G5ppx*G5px*G5x*pow(H,4) - 2538*G5p*G5ppx*G5px*G5x*pow(H,4) + 1080*G4*G5ppxx*G5px*G5x*pow(H,4) + 780*G4px*pow(G5px,2)*G5x*pow(H,4) - 1716*G5pp*pow(G5px,2)*G5x*pow(H,4) + 
   2136*G4*G4pxx*G5pxx*G5x*pow(H,4) - 1032*G4px*G4x*G5pxx*G5x*pow(H,4) - 2904*G4p*G4xx*G5pxx*G5x*pow(H,4) + 1560*G4px*G5p*G5pxx*G5x*pow(H,4) + 192*G4x*G5pp*G5pxx*G5x*pow(H,4) + 
   78*G5p*G5pp*G5pxx*G5x*pow(H,4) - 876*G4*G5ppx*G5pxx*G5x*pow(H,4) + 1392*G4p*G5px*G5pxx*G5x*pow(H,4) - 720*G4*G4px*G5pxxx*G5x*pow(H,4) + 456*G4p*G4x*G5pxxx*G5x*pow(H,4) - 
   288*G4p*G5p*G5pxxx*G5x*pow(H,4) + 24*G4*G5pp*G5pxxx*G5x*pow(H,4) - 3312*G4px*G4pxx*pow(G5x,2)*pow(H,4) + 720*G4p*G4pxxx*pow(G5x,2)*pow(H,4) - 1512*G3pxx*G4x*pow(G5x,2)*pow(H,4) + 
   2232*G4ppxx*G4x*pow(G5x,2)*pow(H,4) + 918*G3px*G4xx*pow(G5x,2)*pow(H,4) + 1260*G4ppx*G4xx*pow(G5x,2)*pow(H,4) - 72*G4pp*G4xxx*pow(G5x,2)*pow(H,4) + 
   1134*G3pxx*G5p*pow(G5x,2)*pow(H,4) - 1512*G4ppxx*G5p*pow(G5x,2)*pow(H,4) - 252*G4pxx*G5pp*pow(G5x,2)*pow(H,4) - 648*G4xx*G5ppp*pow(G5x,2)*pow(H,4) + 
   396*G4x*G5pppx*pow(G5x,2)*pow(H,4) - 378*G5p*G5pppx*pow(G5x,2)*pow(H,4) + 918*G4px*G5ppx*pow(G5x,2)*pow(H,4) + 306*G5pp*G5ppx*pow(G5x,2)*pow(H,4) - 
   288*G4p*G5ppxx*pow(G5x,2)*pow(H,4) - 516*G3px*G5px*pow(G5x,2)*pow(H,4) - 360*G4ppx*G5px*pow(G5x,2)*pow(H,4) + 252*G5ppp*G5px*pow(G5x,2)*pow(H,4) + 
   45*G2x*G5pxx*pow(G5x,2)*pow(H,4) - 90*G3p*G5pxx*pow(G5x,2)*pow(H,4) - 54*G4pp*G5pxx*pow(G5x,2)*pow(H,4) - 90*G2pxx*pow(G5x,3)*pow(H,4) + 180*G4pppx*pow(G5x,3)*pow(H,4) + 
   1296*G4*G4pxxx*G4x*G5xx*pow(H,4) + 3744*G4pxx*pow(G4x,2)*G5xx*pow(H,4) - 1440*G4*G4pxx*G4xx*G5xx*pow(H,4) - 2592*G4px*G4x*G4xx*G5xx*pow(H,4) + 4608*G4p*pow(G4xx,2)*G5xx*pow(H,4) + 
   144*G4*G4px*G4xxx*G5xx*pow(H,4) + 1152*G4p*G4x*G4xxx*G5xx*pow(H,4) - 1008*G4*G4pxxx*G5p*G5xx*pow(H,4) - 6336*G4pxx*G4x*G5p*G5xx*pow(H,4) + 2520*G4px*G4xx*G5p*G5xx*pow(H,4) - 
   936*G4p*G4xxx*G5p*G5xx*pow(H,4) + 2808*G4pxx*pow(G5p,2)*G5xx*pow(H,4) - 2808*G4x*G4xx*G5pp*G5xx*pow(H,4) + 72*G4*G4xxx*G5pp*G5xx*pow(H,4) + 1944*G4xx*G5p*G5pp*G5xx*pow(H,4) - 
   1312*pow(G4x,2)*G5ppx*G5xx*pow(H,4) + 1176*G4*G4xx*G5ppx*G5xx*pow(H,4) + 2296*G4x*G5p*G5ppx*G5xx*pow(H,4) - 1000*pow(G5p,2)*G5ppx*G5xx*pow(H,4) - 776*G4*G4x*G5ppxx*G5xx*pow(H,4) + 
   520*G4*G5p*G5ppxx*G5xx*pow(H,4) + 360*G4*G4pxx*G5px*G5xx*pow(H,4) - 72*G4px*G4x*G5px*G5xx*pow(H,4) - 5064*G4p*G4xx*G5px*G5xx*pow(H,4) - 1308*G4px*G5p*G5px*G5xx*pow(H,4) + 
   1824*G4x*G5pp*G5px*G5xx*pow(H,4) - 1266*G5p*G5pp*G5px*G5xx*pow(H,4) - 372*G4*G5ppx*G5px*G5xx*pow(H,4) + 1416*G4p*pow(G5px,2)*G5xx*pow(H,4) + 312*G4*G4px*G5pxx*G5xx*pow(H,4) - 
   1104*G4p*G4x*G5pxx*G5xx*pow(H,4) + 936*G4p*G5p*G5pxx*G5xx*pow(H,4) + 132*G4*G5pp*G5pxx*G5xx*pow(H,4) - 48*G4*G4p*G5pxxx*G5xx*pow(H,4) + 360*G3pxx*G4*G5x*G5xx*pow(H,4) - 
   960*G4*G4ppxx*G5x*G5xx*pow(H,4) + 1128*pow(G4px,2)*G5x*G5xx*pow(H,4) - 1704*G4p*G4pxx*G5x*G5xx*pow(H,4) + 1008*G3px*G4x*G5x*G5xx*pow(H,4) - 456*G4ppx*G4x*G5x*G5xx*pow(H,4) - 
   3792*G4pp*G4xx*G5x*G5xx*pow(H,4) - 810*G3px*G5p*G5x*G5xx*pow(H,4) - 228*G4ppx*G5p*G5x*G5xx*pow(H,4) + 2268*G4px*G5pp*G5x*G5xx*pow(H,4) - 336*pow(G5pp,2)*G5x*G5xx*pow(H,4) - 
   216*G4x*G5ppp*G5x*G5xx*pow(H,4) + 480*G5p*G5ppp*G5x*G5xx*pow(H,4) + 120*G4*G5pppx*G5x*G5xx*pow(H,4) + 264*G4p*G5ppx*G5x*G5xx*pow(H,4) + 60*G2x*G5px*G5x*G5xx*pow(H,4) - 
   120*G3p*G5px*G5x*G5xx*pow(H,4) + 2184*G4pp*G5px*G5x*G5xx*pow(H,4) - 72*G2px*pow(G5x,2)*G5xx*pow(H,4) + 225*G3pp*pow(G5x,2)*G5xx*pow(H,4) - 162*G4ppp*pow(G5x,2)*G5xx*pow(H,4) + 
   24*G3px*G4*pow(G5xx,2)*pow(H,4) + 216*G4*G4ppx*pow(G5xx,2)*pow(H,4) - 1176*G4p*G4px*pow(G5xx,2)*pow(H,4) - 516*G4pp*G4x*pow(G5xx,2)*pow(H,4) + 600*G4pp*G5p*pow(G5xx,2)*pow(H,4) + 
   198*G4p*G5pp*pow(G5xx,2)*pow(H,4) - 204*G4*G5ppp*pow(G5xx,2)*pow(H,4) - 39*G2p*G5x*pow(G5xx,2)*pow(H,4) - 864*G4*G4pxx*G4x*G5xxx*pow(H,4) + 864*G4p*G4x*G4xx*G5xxx*pow(H,4) + 
   576*G4*G4pxx*G5p*G5xxx*pow(H,4) + 288*G4px*G4x*G5p*G5xxx*pow(H,4) - 720*G4p*G4xx*G5p*G5xxx*pow(H,4) - 216*G4px*pow(G5p,2)*G5xxx*pow(H,4) - 288*pow(G4x,2)*G5pp*G5xxx*pow(H,4) + 
   144*G4*G4xx*G5pp*G5xxx*pow(H,4) + 216*G4x*G5p*G5pp*G5xxx*pow(H,4) + 416*G4*G4x*G5ppx*G5xxx*pow(H,4) - 304*G4*G5p*G5ppx*G5xxx*pow(H,4) + 16*pow(G4,2)*G5ppxx*G5xxx*pow(H,4) + 
   360*G4*G4px*G5px*G5xxx*pow(H,4) - 384*G4p*G4x*G5px*G5xxx*pow(H,4) + 300*G4p*G5p*G5px*G5xxx*pow(H,4) - 192*G4*G5pp*G5px*G5xxx*pow(H,4) + 24*G4*G4p*G5pxx*G5xxx*pow(H,4) - 
   180*G3px*G4*G5x*G5xxx*pow(H,4) + 504*G4*G4ppx*G5x*G5xxx*pow(H,4) - 336*G4p*G4px*G5x*G5xxx*pow(H,4) - 96*G4pp*G4x*G5x*G5xxx*pow(H,4) + 36*G4pp*G5p*G5x*G5xxx*pow(H,4) + 
   12*G4p*G5pp*G5x*G5xxx*pow(H,4) - 24*G4*G5ppp*G5x*G5xxx*pow(H,4) + 18*G2p*pow(G5x,2)*G5xxx*pow(H,4) - 48*pow(G4p,2)*G5xx*G5xxx*pow(H,4) - 24*G4*G4pp*G5xx*G5xxx*pow(H,4) + 
   2790*G4xx*G5px*pow(G5x,2)*pow(H,6) - 1470*pow(G5px,2)*pow(G5x,2)*pow(H,6) - 522*G4x*G5pxx*pow(G5x,2)*pow(H,6) + 54*G5p*G5pxx*pow(G5x,2)*pow(H,6) + 
   252*G4*G5pxxx*pow(G5x,2)*pow(H,6) + 345*G5ppx*pow(G5x,3)*pow(H,6) + 2220*G4x*G5px*G5x*G5xx*pow(H,6) - 1590*G5p*G5px*G5x*G5xx*pow(H,6) - 36*G4*G5pxx*G5x*G5xx*pow(H,6) + 
   378*G4px*pow(G5x,2)*G5xx*pow(H,6) - 1215*G5pp*pow(G5x,2)*G5xx*pow(H,6) - 120*G4*G5px*pow(G5xx,2)*pow(H,6) + 666*G4p*G5x*pow(G5xx,2)*pow(H,6) - 180*G4*G5px*G5x*G5xxx*pow(H,6) + 
   216*G4p*pow(G5x,2)*G5xxx*pow(H,6) - 3*pow(G3x,2)*(36*G3pxx*G4x - 48*G4ppxx*G4x - 90*G3px*G4xx + 36*G4ppx*G4xx - 18*G3pxx*G5p + 24*G4ppxx*G5p + 18*G3xx*(2*G4px - G5pp) + 72*G4xx*G5ppp - 
      12*G4x*G5pppx + 6*G5p*G5pppx + 18*G5pp*G5ppx + 6*G2xx*G5px + 48*G3px*G5px - 36*G4ppx*G5px - 36*G5ppp*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + 6*G2pxx*G5x - 12*G4pppx*G5x - 9*G3pp*G5xx + 
      18*G4ppp*G5xx + 198*G4xx*G5px*pow(H,2) - 116*pow(G5px,2)*pow(H,2) - 138*G4x*G5pxx*pow(H,2) + 114*G5p*G5pxx*pow(H,2) - 12*G4*G5pxxx*pow(H,2) - 360*G4pxx*G5x*pow(H,2) + 
      151*G5ppx*G5x*pow(H,2) - 21*G5pp*G5xx*pow(H,2) - 12*G4p*G5xxx*pow(H,2) + 18*G4px*(12*G4pxx - 7*G5ppx + 7*G5xx*pow(H,2))) - 144*G4xx*G4xxx*G5px*pow(H,2)*rptot + 
   72*G4xxx*pow(G5px,2)*pow(H,2)*rptot + 144*pow(G4xx,2)*G5pxx*pow(H,2)*rptot + 72*G4x*G4xxx*G5pxx*pow(H,2)*rptot - 36*G4xxx*G5p*G5pxx*pow(H,2)*rptot - 72*G4xx*G5px*G5pxx*pow(H,2)*rptot - 
   36*G4x*pow(G5pxx,2)*pow(H,2)*rptot + 18*G5p*pow(G5pxx,2)*pow(H,2)*rptot - 48*G4x*G4xx*G5pxxx*pow(H,2)*rptot + 24*G4xx*G5p*G5pxxx*pow(H,2)*rptot + 24*G4x*G5px*G5pxxx*pow(H,2)*rptot - 
   12*G5p*G5px*G5pxxx*pow(H,2)*rptot - 144*G4pxxx*G4xx*G5x*pow(H,2)*rptot + 144*G4pxx*G4xxx*G5x*pow(H,2)*rptot - 72*G4xxx*G5ppx*G5x*pow(H,2)*rptot + 72*G4xx*G5ppxx*G5x*pow(H,2)*rptot + 
   72*G4pxxx*G5px*G5x*pow(H,2)*rptot - 36*G5ppxx*G5px*G5x*pow(H,2)*rptot - 108*G4pxx*G5pxx*G5x*pow(H,2)*rptot + 36*G5ppx*G5pxx*G5x*pow(H,2)*rptot + 12*G4px*G5pxxx*G5x*pow(H,2)*rptot - 
   72*G4pxxx*G4x*G5xx*pow(H,2)*rptot + 288*G4pxx*G4xx*G5xx*pow(H,2)*rptot - 72*G4px*G4xxx*G5xx*pow(H,2)*rptot + 36*G4pxxx*G5p*G5xx*pow(H,2)*rptot + 36*G4xxx*G5pp*G5xx*pow(H,2)*rptot - 
   144*G4xx*G5ppx*G5xx*pow(H,2)*rptot + 36*G4x*G5ppxx*G5xx*pow(H,2)*rptot - 18*G5p*G5ppxx*G5xx*pow(H,2)*rptot - 108*G4pxx*G5px*G5xx*pow(H,2)*rptot + 72*G5ppx*G5px*G5xx*pow(H,2)*rptot - 
   18*G5pp*G5pxx*G5xx*pow(H,2)*rptot - 18*G3pxx*G5x*G5xx*pow(H,2)*rptot + 36*G4ppxx*G5x*G5xx*pow(H,2)*rptot + 9*G3px*pow(G5xx,2)*pow(H,2)*rptot - 18*G4ppx*pow(G5xx,2)*pow(H,2)*rptot + 
   48*G4pxx*G4x*G5xxx*pow(H,2)*rptot - 48*G4px*G4xx*G5xxx*pow(H,2)*rptot - 24*G4pxx*G5p*G5xxx*pow(H,2)*rptot + 24*G4xx*G5pp*G5xxx*pow(H,2)*rptot - 24*G4x*G5ppx*G5xxx*pow(H,2)*rptot + 
   12*G5p*G5ppx*G5xxx*pow(H,2)*rptot + 36*G4px*G5px*G5xxx*pow(H,2)*rptot - 12*G5pp*G5px*G5xxx*pow(H,2)*rptot + 6*G3px*G5x*G5xxx*pow(H,2)*rptot - 12*G4ppx*G5x*G5xxx*pow(H,2)*rptot - 
   18*G5pxxx*pow(G5x,2)*pow(H,4)*rptot + 54*G5pxx*G5x*G5xx*pow(H,4)*rptot - 15*G5px*pow(G5xx,2)*pow(H,4)*rptot - 
   6*G3xx*(120*pow(G4px,3) + 48*G3p*G4pxx*G4x - 48*G4pp*G4pxx*G4x - 24*G3ppx*pow(G4x,2) + 48*G4pppx*pow(G4x,2) + 24*G3pp*G4x*G4xx - 48*G4ppp*G4x*G4xx - 24*G3p*G4pxx*G5p + 24*G4pp*G4pxx*G5p + 
      24*G3ppx*G4x*G5p - 48*G4pppx*G4x*G5p - 12*G3pp*G4xx*G5p + 24*G4ppp*G4xx*G5p - 6*G3ppx*pow(G5p,2) + 12*G4pppx*pow(G5p,2) - 12*G3px*G4x*G5pp + 24*G4ppx*G4x*G5pp + 24*G3p*G4xx*G5pp - 
      24*G4pp*G4xx*G5pp + 6*G3px*G5p*G5pp - 12*G4ppx*G5p*G5pp - 24*G3p*G4x*G5ppx + 24*G4pp*G4x*G5ppx + 12*G3p*G5p*G5ppx - 12*G4pp*G5p*G5ppx - 12*G3pp*G4x*G5px + 24*G4ppp*G4x*G5px + 6*G3pp*G5p*G5px - 
      12*G4ppp*G5p*G5px - 12*G3p*G5pp*G5px + 12*G4pp*G5pp*G5px + 6*G3p*G3px*G5x - 6*G3px*G4pp*G5x - 12*G3p*G4ppx*G5x + 12*G4pp*G4ppx*G5x - 8*G2ppx*G4x*G5x + 8*G3ppp*G4x*G5x + 4*G2ppx*G5p*G5x - 
      4*G3ppp*G5p*G5x - 3*G2x*(8*G4pxx*G4x - 4*G4pxx*G5p + 4*G4xx*G5pp - 4*G4x*G5ppx + 2*G5p*G5ppx - 2*G5pp*G5px + G3px*G5x - 2*G4ppx*G5x) - 288*G4pxx*pow(G4x,2)*pow(H,2) + 
      144*G4*G4pxx*G4xx*pow(H,2) - 144*G4p*pow(G4xx,2)*pow(H,2) + 360*G4pxx*G4x*G5p*pow(H,2) - 108*G4pxx*pow(G5p,2)*pow(H,2) - 72*G4x*G4xx*G5pp*pow(H,2) + 108*G4xx*G5p*G5pp*pow(H,2) + 
      152*pow(G4x,2)*G5ppx*pow(H,2) - 96*G4*G4xx*G5ppx*pow(H,2) - 200*G4x*G5p*G5ppx*pow(H,2) + 62*pow(G5p,2)*G5ppx*pow(H,2) + 16*G4*G4x*G5ppxx*pow(H,2) - 8*G4*G5p*G5ppxx*pow(H,2) - 
      96*G4*G4pxx*G5px*pow(H,2) + 108*G4p*G4xx*G5px*pow(H,2) + 12*G4x*G5pp*G5px*pow(H,2) - 60*G5p*G5pp*G5px*pow(H,2) + 60*G4*G5ppx*G5px*pow(H,2) - 12*G4p*pow(G5px,2)*pow(H,2) + 
      12*G4p*G4x*G5pxx*pow(H,2) - 6*G4p*G5p*G5pxx*pow(H,2) + 12*G4*G5pp*G5pxx*pow(H,2) + 48*G4*G4ppxx*G5x*pow(H,2) + 24*G4p*G4pxx*G5x*pow(H,2) - 108*G3px*G4x*G5x*pow(H,2) + 
      336*G4ppx*G4x*G5x*pow(H,2) - 12*G4pp*G4xx*G5x*pow(H,2) + 63*G3px*G5p*G5x*pow(H,2) - 210*G4ppx*G5p*G5x*pow(H,2) + 6*pow(G5pp,2)*G5x*pow(H,2) - 36*G4x*G5ppp*G5x*pow(H,2) + 
      30*G5p*G5ppp*G5x*pow(H,2) - 24*G4*G5pppx*G5x*pow(H,2) - 18*G4p*G5ppx*G5x*pow(H,2) - 6*G2px*pow(G5x,2)*pow(H,2) + 21*G3pp*pow(G5x,2)*pow(H,2) - 18*G4ppp*pow(G5x,2)*pow(H,2) + 
      6*G3px*G4*G5xx*pow(H,2) - 36*G4*G4ppx*G5xx*pow(H,2) - 12*G4pp*G4x*G5xx*pow(H,2) + 6*G4pp*G5p*G5xx*pow(H,2) + 12*G4*G5ppp*G5xx*pow(H,2) - 3*G2p*G5x*G5xx*pow(H,2) - 
      30*G4x*G5px*G5x*pow(H,4) - 15*G5p*G5px*G5x*pow(H,4) + 42*G4*G5pxx*G5x*pow(H,4) + 18*G5pp*pow(G5x,2)*pow(H,4) + 6*G4*G5px*G5xx*pow(H,4) - 54*G4p*G5x*G5xx*pow(H,4) - 
      12*pow(G4px,2)*(7*G5pp + 11*G5x*pow(H,2)) + 6*G4px*(10*G3px*G4x - 28*G4ppx*G4x + 4*G2x*G4xx - 8*G3p*G4xx + 8*G4pp*G4xx - 5*G3px*G5p + 14*G4ppx*G5p + 2*pow(G5pp,2) + 4*G4x*G5ppp - 
         2*G5p*G5ppp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px - G3pp*G5x + 2*G4ppp*G5x + 48*G4x*G4xx*pow(H,2) - 48*G4xx*G5p*pow(H,2) - 24*G4x*G5px*pow(H,2) + 39*G5p*G5px*pow(H,2) - 
         10*G4*G5pxx*pow(H,2) + 3*G5pp*G5x*pow(H,2) + 8*G4p*G5xx*pow(H,2) + 6*pow(G5x,2)*pow(H,4)) - 3*G5pxx*G5x*pow(H,2)*rptot + 3*G5px*G5xx*pow(H,2)*rptot) + 
   3*G3x*(48*G2xx*G4pxx*G4x + 24*G3px*G4pxx*G4x - 48*G4ppx*G4pxx*G4x - 24*G2x*G4pxxx*G4x + 48*G3p*G4pxxx*G4x - 48*G4pp*G4pxxx*G4x + 96*G2x*G4pxx*G4xx - 192*G3p*G4pxx*G4xx + 192*G4pp*G4pxx*G4xx - 
      72*G2pxx*G4x*G4xx - 24*G3ppx*G4x*G4xx + 192*G4pppx*G4x*G4xx + 144*G3pp*pow(G4xx,2) - 288*G4ppp*pow(G4xx,2) - 24*G3pp*G4x*G4xxx + 48*G4ppp*G4x*G4xxx - 24*G2xx*G4pxx*G5p - 12*G3px*G4pxx*G5p + 
      24*G4ppx*G4pxx*G5p + 12*G2x*G4pxxx*G5p - 24*G3p*G4pxxx*G5p + 24*G4pp*G4pxxx*G5p + 36*G2pxx*G4xx*G5p + 12*G3ppx*G4xx*G5p - 96*G4pppx*G4xx*G5p + 12*G3pp*G4xxx*G5p - 24*G4ppp*G4xxx*G5p - 
      24*G3pxx*G4x*G5pp + 48*G4ppxx*G4x*G5pp + 36*G2xx*G4xx*G5pp + 60*G3px*G4xx*G5pp - 192*G4ppx*G4xx*G5pp + 12*G2x*G4xxx*G5pp - 24*G3p*G4xxx*G5pp + 24*G4pp*G4xxx*G5pp + 12*G3pxx*G5p*G5pp - 
      24*G4ppxx*G5p*G5pp - 24*G4pxx*pow(G5pp,2) - 48*G4pxx*G4x*G5ppp + 24*G4pxx*G5p*G5ppp - 28*G2xx*G4x*G5ppx - 44*G3px*G4x*G5ppx + 144*G4ppx*G4x*G5ppx - 48*G2x*G4xx*G5ppx + 96*G3p*G4xx*G5ppx - 
      96*G4pp*G4xx*G5ppx + 14*G2xx*G5p*G5ppx + 22*G3px*G5p*G5ppx - 72*G4ppx*G5p*G5ppx + 12*G2x*G4x*G5ppxx - 24*G3p*G4x*G5ppxx + 24*G4pp*G4x*G5ppxx - 6*G2x*G5p*G5ppxx + 12*G3p*G5p*G5ppxx - 
      12*G4pp*G5p*G5ppxx - 36*G2x*G4pxx*G5px + 72*G3p*G4pxx*G5px - 72*G4pp*G4pxx*G5px + 40*G2pxx*G4x*G5px + 8*G3ppx*G4x*G5px - 96*G4pppx*G4x*G5px - 16*G2px*G4xx*G5px - 128*G3pp*G4xx*G5px + 
      288*G4ppp*G4xx*G5px - 20*G2pxx*G5p*G5px - 4*G3ppx*G5p*G5px + 48*G4pppx*G5p*G5px - 24*G2xx*G5pp*G5px - 24*G3px*G5pp*G5px + 96*G4ppx*G5pp*G5px + 24*G2x*G5ppx*G5px - 48*G3p*G5ppx*G5px + 
      48*G4pp*G5ppx*G5px + 8*G2px*pow(G5px,2) + 28*G3pp*pow(G5px,2) - 72*G4ppp*pow(G5px,2) + 8*G2px*G4x*G5pxx + 4*G3pp*G4x*G5pxx - 24*G4ppp*G4x*G5pxx - 4*G2px*G5p*G5pxx - 2*G3pp*G5p*G5pxx + 
      12*G4ppp*G5p*G5pxx - 6*G2x*G5pp*G5pxx + 12*G3p*G5pp*G5pxx - 12*G4pp*G5pp*G5pxx + 6*G2xx*G3px*G5x + 6*pow(G3px,2)*G5x - 6*G2x*G3pxx*G5x + 12*G3p*G3pxx*G5x - 12*G3pxx*G4pp*G5x - 20*G2xx*G4ppx*G5x - 
      28*G3px*G4ppx*G5x + 48*pow(G4ppx,2)*G5x + 12*G2x*G4ppxx*G5x - 24*G3p*G4ppxx*G5x + 24*G4pp*G4ppxx*G5x + 16*G2px*G4pxx*G5x - 4*G3pp*G4pxx*G5x - 24*G4ppp*G4pxx*G5x - 16*G2ppx*G4xx*G5x + 
      16*G3ppp*G4xx*G5x - 4*G2pxx*G5pp*G5x + 4*G3ppx*G5pp*G5x + 4*G2xx*G5ppp*G5x - 4*G3px*G5ppp*G5x - 8*G2px*G5ppx*G5x + 8*G3pp*G5ppx*G5x + 8*G2ppx*G5px*G5x - 8*G3ppp*G5px*G5x + 6*G2x*G3px*G5xx - 
      12*G3p*G3px*G5xx + 12*G3px*G4pp*G5xx - 12*G2x*G4ppx*G5xx + 24*G3p*G4ppx*G5xx - 24*G4pp*G4ppx*G5xx - 8*G2ppx*G4x*G5xx + 8*G3ppp*G4x*G5xx + 4*G2ppx*G5p*G5xx - 4*G3ppp*G5p*G5xx + 4*G2px*G5pp*G5xx - 
      4*G3pp*G5pp*G5xx - 720*G4pxxx*pow(G4x,2)*pow(H,2) + 432*G4*G4pxxx*G4xx*pow(H,2) + 3456*G4pxx*G4x*G4xx*pow(H,2) - 288*G4*G4pxx*G4xxx*pow(H,2) + 432*G4p*G4xx*G4xxx*pow(H,2) + 
      936*G4pxxx*G4x*G5p*pow(H,2) - 2880*G4pxx*G4xx*G5p*pow(H,2) - 288*G4pxxx*pow(G5p,2)*pow(H,2) + 144*pow(G4xx,2)*G5pp*pow(H,2) - 108*G4xxx*G5p*G5pp*pow(H,2) - 
      1480*G4x*G4xx*G5ppx*pow(H,2) + 168*G4*G4xxx*G5ppx*pow(H,2) + 1172*G4xx*G5p*G5ppx*pow(H,2) + 360*pow(G4x,2)*G5ppxx*pow(H,2) - 184*G4*G4xx*G5ppxx*pow(H,2) - 
      452*G4x*G5p*G5ppxx*pow(H,2) + 136*pow(G5p,2)*G5ppxx*pow(H,2) - 240*G4*G4pxxx*G5px*pow(H,2) - 1808*G4pxx*G4x*G5px*pow(H,2) - 144*G4p*G4xxx*G5px*pow(H,2) + 
      1468*G4pxx*G5p*G5px*pow(H,2) - 220*G4xx*G5pp*G5px*pow(H,2) + 804*G4x*G5ppx*G5px*pow(H,2) - 666*G5p*G5ppx*G5px*pow(H,2) + 104*G4*G5ppxx*G5px*pow(H,2) + 92*G5pp*pow(G5px,2)*pow(H,2) + 
      136*G4*G4pxx*G5pxx*pow(H,2) - 448*G4p*G4xx*G5pxx*pow(H,2) + 22*G5p*G5pp*G5pxx*pow(H,2) - 44*G4*G5ppx*G5pxx*pow(H,2) + 208*G4p*G5px*G5pxx*pow(H,2) + 32*G4p*G4x*G5pxxx*pow(H,2) - 
      16*G4p*G5p*G5pxxx*pow(H,2) + 8*G4*G5pp*G5pxxx*pow(H,2) + 96*G4p*G4pxxx*G5x*pow(H,2) - 396*G3pxx*G4x*G5x*pow(H,2) + 744*G4ppxx*G4x*G5x*pow(H,2) + 468*G3px*G4xx*G5x*pow(H,2) - 
      240*G4ppx*G4xx*G5x*pow(H,2) - 24*G4pp*G4xxx*G5x*pow(H,2) + 252*G3pxx*G5p*G5x*pow(H,2) - 432*G4ppxx*G5p*G5x*pow(H,2) + 36*G4pxx*G5pp*G5x*pow(H,2) - 192*G4xx*G5ppp*G5x*pow(H,2) + 
      24*G4x*G5pppx*G5x*pow(H,2) - 36*G5p*G5pppx*G5x*pow(H,2) - 24*G5pp*G5ppx*G5x*pow(H,2) - 40*G4p*G5ppxx*G5x*pow(H,2) - 10*G2xx*G5px*G5x*pow(H,2) - 252*G3px*G5px*G5x*pow(H,2) + 
      188*G4ppx*G5px*G5x*pow(H,2) + 84*G5ppp*G5px*G5x*pow(H,2) + 34*G2x*G5pxx*G5x*pow(H,2) - 68*G3p*G5pxx*G5x*pow(H,2) + 36*G4pp*G5pxx*G5x*pow(H,2) - 36*G2pxx*pow(G5x,2)*pow(H,2) + 
      24*G3ppx*pow(G5x,2)*pow(H,2) + 24*G4pppx*pow(G5x,2)*pow(H,2) + 48*G3pxx*G4*G5xx*pow(H,2) - 48*G4*G4ppxx*G5xx*pow(H,2) - 288*G4p*G4pxx*G5xx*pow(H,2) + 168*G3px*G4x*G5xx*pow(H,2) - 
      160*G4ppx*G4x*G5xx*pow(H,2) - 552*G4pp*G4xx*G5xx*pow(H,2) - 150*G3px*G5p*G5xx*pow(H,2) + 68*G4ppx*G5p*G5xx*pow(H,2) - 8*pow(G5pp,2)*G5xx*pow(H,2) - 16*G4x*G5ppp*G5xx*pow(H,2) + 
      80*G5p*G5ppp*G5xx*pow(H,2) - 24*G4*G5pppx*G5xx*pow(H,2) + 60*G4p*G5ppx*G5xx*pow(H,2) - 28*G2x*G5px*G5xx*pow(H,2) + 56*G3p*G5px*G5xx*pow(H,2) + 272*G4pp*G5px*G5xx*pow(H,2) - 
      20*G2px*G5x*G5xx*pow(H,2) + 68*G3pp*G5x*G5xx*pow(H,2) - 56*G4ppp*G5x*G5xx*pow(H,2) - 9*G2p*pow(G5xx,2)*pow(H,2) - 12*G3px*G4*G5xxx*pow(H,2) + 40*G4*G4ppx*G5xxx*pow(H,2) - 
      8*G4pp*G4x*G5xxx*pow(H,2) + 4*G4pp*G5p*G5xxx*pow(H,2) - 8*G4p*G5pp*G5xxx*pow(H,2) - 8*G4*G5ppp*G5xxx*pow(H,2) + 2*G2p*G5x*G5xxx*pow(H,2) + 420*G4xx*G5px*G5x*pow(H,4) - 
      210*pow(G5px,2)*G5x*pow(H,4) + 204*G4x*G5pxx*G5x*pow(H,4) - 312*G5p*G5pxx*G5x*pow(H,4) + 96*G4*G5pxxx*G5x*pow(H,4) + 468*G4pxx*pow(G5x,2)*pow(H,4) - 
      101*G5ppx*pow(G5x,2)*pow(H,4) + 196*G4x*G5px*G5xx*pow(H,4) - 26*G5p*G5px*G5xx*pow(H,4) - 68*G4*G5pxx*G5xx*pow(H,4) - 264*G5pp*G5x*G5xx*pow(H,4) + 142*G4p*pow(G5xx,2)*pow(H,4) - 
      28*G4*G5px*G5xxx*pow(H,4) + 60*G4p*G5x*G5xxx*pow(H,4) + 4*pow(G4px,2)*(84*G4pxx - 81*G5ppx + 178*G5xx*pow(H,2)) - 
      4*G4px*(72*G4ppxx*G4x + 18*G2xx*G4xx + 120*G3px*G4xx - 132*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx - 36*G4ppxx*G5p + 24*G3pxx*(-2*G4x + G5p) - 24*G4pxx*G5pp - 72*G4xx*G5ppp + 
         12*G4x*G5pppx - 6*G5p*G5pppx - 18*G5pp*G5ppx - 18*G2xx*G5px - 60*G3px*G5px + 84*G4ppx*G5px + 36*G5ppp*G5px - 8*G2pxx*G5x + 2*G3ppx*G5x + 12*G4pppx*G5x + 2*G2px*G5xx + 7*G3pp*G5xx - 
         18*G4ppp*G5xx + 360*pow(G4xx,2)*pow(H,2) + 72*G4x*G4xxx*pow(H,2) - 90*G4xxx*G5p*pow(H,2) - 698*G4xx*G5px*pow(H,2) + 288*pow(G5px,2)*pow(H,2) + 140*G4x*G5pxx*pow(H,2) - 
         104*G5p*G5pxx*pow(H,2) + 16*G4*G5pxxx*pow(H,2) + 426*G4pxx*G5x*pow(H,2) - 208*G5ppx*G5x*pow(H,2) + 43*G5pp*G5xx*pow(H,2) + 8*G4p*G5xxx*pow(H,2) + 30*G5x*G5xx*pow(H,4)) + 
      6*G3xx*(32*pow(G4px,2) + 6*G3px*G4x - 20*G4ppx*G4x - 3*G3px*G5p + 10*G4ppx*G5p + 2*pow(G5pp,2) + 4*G4x*G5ppp - 2*G5p*G5ppp - G2x*G5px + 2*G3p*G5px - 2*G4pp*G5px - G3pp*G5x + 2*G4ppp*G5x - 
         4*G4x*G5px*pow(H,2) + 11*G5p*G5px*pow(H,2) - 6*G4*G5pxx*pow(H,2) + 3*G5pp*G5x*pow(H,2) + 8*G4p*G5xx*pow(H,2) - 2*G4px*(10*G5pp + 9*G5x*pow(H,2))) - 2*G5pxxx*G5x*pow(H,2)*rptot + 
      6*G5pxx*G5xx*pow(H,2)*rptot - 2*G5px*G5xxx*pow(H,2)*rptot)) + 8*pow(a,17)*
 (-2*G2pp*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
   G2p*(2*G2px*pow(G4,2) - 4*G3pp*pow(G4,2) + 3*(G2x*G4*G4p - 2*G3p*G4*G4p + pow(G4p,3) + 4*G4*G4p*G4pp + 4*pow(G4,2)*G4px*pow(H,2) + 6*G4*G4p*G4x*pow(H,2) - 6*G4*G4p*G5p*pow(H,2) - 
         4*pow(G4,2)*G5pp*pow(H,2))) + 3*(6*G2x*G4*pow(G4p,2)*pow(H,2) - 12*G3p*G4*pow(G4p,2)*pow(H,2) - 6*pow(G4p,4)*pow(H,2) - 8*G2x*pow(G4,2)*G4pp*pow(H,2) + 
      16*G3p*pow(G4,2)*G4pp*pow(H,2) + 24*G4*pow(G4p,2)*G4pp*pow(H,2) + 48*pow(G4,2)*G4p*G4px*pow(H,4) + 36*G4*pow(G4p,2)*G4x*pow(H,4) - 48*pow(G4,2)*G4pp*G4x*pow(H,4) - 
      36*G4*pow(G4p,2)*G5p*pow(H,4) + 48*pow(G4,2)*G4pp*G5p*pow(H,4) - 48*pow(G4,2)*G4p*G5pp*pow(H,4) + G2px*G4*G4p*(8*G4*pow(H,2) - rptot) - G2x*pow(G4p,2)*rptot + 
      2*G3p*pow(G4p,2)*rptot + G2x*G4*G4pp*rptot - 2*G3p*G4*G4pp*rptot - 3*pow(G4p,2)*G4pp*rptot - 6*G4*G4p*G4px*pow(H,2)*rptot - 6*pow(G4p,2)*G4x*pow(H,2)*rptot + 
      6*G4*G4pp*G4x*pow(H,2)*rptot + 6*pow(G4p,2)*G5p*pow(H,2)*rptot - 6*G4*G4pp*G5p*pow(H,2)*rptot + 6*G4*G4p*G5pp*pow(H,2)*rptot + 2*G3pp*G4*G4p*(-8*G4*pow(H,2) + rptot))) - 
2*pow(a,12)*pow(dphi,5)*H*(-36*G2px*pow(G3x,2)*G4 + 36*G3pp*pow(G3x,2)*G4 + 18*G2p*G3x*G3xx*G4 - 24*G2xx*G3ppx*pow(G4,2) + 24*G2pxx*G3px*pow(G4,2) + 24*G2px*G3pxx*pow(G4,2) - 
   24*G3pp*G3pxx*pow(G4,2) - 24*G2ppx*G3xx*pow(G4,2) + 24*G3ppp*G3xx*pow(G4,2) - 12*G2xx*G3px*G4*G4p - 60*pow(G3px,2)*G4*G4p + 48*G2pxx*G3x*G4*G4p + 24*G3ppx*G3x*G4*G4p + 36*G3pp*G3xx*G4*G4p + 
   48*G2xx*G3x*pow(G4p,2) - 66*G3px*G3x*pow(G4p,2) - 12*G2xx*G3x*G4*G4pp - 60*G3px*G3x*G4*G4pp - 306*pow(G3x,2)*G4p*G4pp + 72*G3pxx*G4*G4p*G4pp + 72*G3xx*pow(G4p,2)*G4pp - 
   72*G3xx*G4*pow(G4pp,2) + 72*pow(G3x,2)*G4*G4ppp - 72*G3xx*G4*G4p*G4ppp + 48*G2xx*pow(G4,2)*G4pppx - 48*G3px*pow(G4,2)*G4pppx - 144*G3x*G4*G4p*G4pppx - 48*G2pxx*pow(G4,2)*G4ppx + 
   48*G3ppx*pow(G4,2)*G4ppx + 72*G2xx*G4*G4p*G4ppx + 216*G3px*G4*G4p*G4ppx - 180*G3x*pow(G4p,2)*G4ppx + 216*G3x*G4*G4pp*G4ppx - 288*G4*G4p*pow(G4ppx,2) - 48*G2px*pow(G4,2)*G4ppxx + 
   48*G3pp*pow(G4,2)*G4ppxx - 144*G4*G4p*G4pp*G4ppxx + 216*G2px*G3x*G4*G4px - 204*G3pp*G3x*G4*G4px - 36*G2p*G3xx*G4*G4px - 144*G2pxx*G4*G4p*G4px - 120*G2xx*pow(G4p,2)*G4px + 
   300*G3px*pow(G4p,2)*G4px + 120*G2xx*G4*G4pp*G4px - 48*G3px*G4*G4pp*G4px + 1440*G3x*G4p*G4pp*G4px - 360*G3x*G4*G4ppp*G4px + 288*G4*G4p*G4pppx*G4px + 72*pow(G4p,2)*G4ppx*G4px - 
   288*G4*G4pp*G4ppx*G4px - 288*G2px*G4*pow(G4px,2) + 264*G3pp*G4*pow(G4px,2) - 1512*G4p*G4pp*pow(G4px,2) + 432*G4*G4ppp*pow(G4px,2) - 84*G2p*G3x*G4*G4pxx + 48*G2ppx*pow(G4,2)*G4pxx - 
   48*G3ppp*pow(G4,2)*G4pxx - 96*G2px*G4*G4p*G4pxx + 24*G3pp*G4*G4p*G4pxx - 288*pow(G4p,2)*G4pp*G4pxx + 144*G4*pow(G4pp,2)*G4pxx + 144*G4*G4p*G4ppp*G4pxx + 168*G2p*G4*G4px*G4pxx + 
   18*G2p*pow(G3x,2)*G4x + 48*G2px*G2xx*G4*G4x - 120*G2xx*G3pp*G4*G4x - 144*G2px*G3px*G4*G4x + 264*G3pp*G3px*G4*G4x + 48*G2p*G3pxx*G4*G4x + 48*G2ppx*G3x*G4*G4x - 48*G3ppp*G3x*G4*G4x - 
   48*G2pp*G3xx*G4*G4x - 168*G2px*G3x*G4p*G4x + 432*G3pp*G3x*G4p*G4x + 36*G2p*G3xx*G4p*G4x + 24*G2pxx*pow(G4p,2)*G4x + 48*G3ppx*pow(G4p,2)*G4x - 48*G2pxx*G4*G4pp*G4x + 48*G3ppx*G4*G4pp*G4x - 
   24*G2xx*G4p*G4pp*G4x - 120*G3px*G4p*G4pp*G4x + 288*G3x*pow(G4pp,2)*G4x + 48*G2xx*G4*G4ppp*G4x - 48*G3px*G4*G4ppp*G4x - 288*G3x*G4p*G4ppp*G4x - 144*pow(G4p,2)*G4pppx*G4x + 288*G2px*G4*G4ppx*G4x - 
   384*G3pp*G4*G4ppx*G4x + 432*G4p*G4pp*G4ppx*G4x - 96*G2p*G4*G4ppxx*G4x - 36*G2p*G3x*G4px*G4x - 192*G2ppx*G4*G4px*G4x + 192*G3ppp*G4*G4px*G4x + 384*G2px*G4p*G4px*G4x - 888*G3pp*G4p*G4px*G4x - 
   720*pow(G4pp,2)*G4px*G4x + 432*G4p*G4ppp*G4px*G4x + 96*G2pp*G4*G4pxx*G4x - 168*G2p*G4p*G4pxx*G4x - 24*G2p*G2xx*pow(G4x,2) + 24*G2p*G3px*pow(G4x,2) - 96*G2ppx*G4p*pow(G4x,2) + 
   96*G3ppp*G4p*pow(G4x,2) + 144*G2px*G4pp*pow(G4x,2) - 192*G3pp*G4pp*pow(G4x,2) + 48*G2p*G4ppx*pow(G4x,2) - 48*G2pp*G4px*pow(G4x,2) + 24*G2p*G2xx*G4*G4xx - 72*G2p*G3px*G4*G4xx + 
   48*G2pp*G3x*G4*G4xx - 168*G2p*G3x*G4p*G4xx + 96*G2ppx*G4*G4p*G4xx - 96*G3ppp*G4*G4p*G4xx + 96*G2px*pow(G4p,2)*G4xx - 312*G3pp*pow(G4p,2)*G4xx - 144*G2px*G4*G4pp*G4xx + 192*G3pp*G4*G4pp*G4xx - 
   288*G4p*pow(G4pp,2)*G4xx + 432*pow(G4p,2)*G4ppp*G4xx + 96*G2p*G4*G4ppx*G4xx - 96*G2pp*G4*G4px*G4xx + 432*G2p*G4p*G4px*G4xx + 96*G2pp*G4p*G4x*G4xx - 96*G2p*G4pp*G4x*G4xx - 36*G2p*pow(G3x,2)*G5p - 
   24*G2px*G2xx*G4*G5p + 72*G2xx*G3pp*G4*G5p + 48*G2px*G3px*G4*G5p - 120*G3pp*G3px*G4*G5p - 24*G2p*G3pxx*G4*G5p + 24*G2pp*G3xx*G4*G5p + 108*G2px*G3x*G4p*G5p - 348*G3pp*G3x*G4p*G5p - 18*G2p*G3xx*G4p*G5p - 
   12*G2pxx*pow(G4p,2)*G5p - 24*G3ppx*pow(G4p,2)*G5p + 48*G2pxx*G4*G4pp*G5p - 48*G3ppx*G4*G4pp*G5p + 36*G2xx*G4p*G4pp*G5p - 36*G3px*G4p*G4pp*G5p - 216*G3x*pow(G4pp,2)*G5p - 48*G2xx*G4*G4ppp*G5p + 
   48*G3px*G4*G4ppp*G5p + 360*G3x*G4p*G4ppp*G5p + 72*pow(G4p,2)*G4pppx*G5p - 96*G2px*G4*G4ppx*G5p + 144*G3pp*G4*G4ppx*G5p - 72*G4p*G4pp*G4ppx*G5p + 48*G2p*G4*G4ppxx*G5p + 150*G2p*G3x*G4px*G5p + 
   48*G2ppx*G4*G4px*G5p - 48*G3ppp*G4*G4px*G5p - 192*G2px*G4p*G4px*G5p + 660*G3pp*G4p*G4px*G5p + 504*pow(G4pp,2)*G4px*G5p - 648*G4p*G4ppp*G4px*G5p - 156*G2p*pow(G4px,2)*G5p - 48*G2pp*G4*G4pxx*G5p + 
   84*G2p*G4p*G4pxx*G5p + 36*G2p*G2xx*G4x*G5p - 60*G2p*G3px*G4x*G5p + 24*G2pp*G3x*G4x*G5p + 144*G2ppx*G4p*G4x*G5p - 144*G3ppp*G4p*G4x*G5p - 216*G2px*G4pp*G4x*G5p + 288*G3pp*G4pp*G4x*G5p - 
   48*G2pp*G4p*G4xx*G5p + 48*G2p*G4pp*G4xx*G5p - 12*G2p*G2xx*pow(G5p,2) + 24*G2p*G3px*pow(G5p,2) - 12*G2pp*G3x*pow(G5p,2) - 48*G2ppx*G4p*pow(G5p,2) + 48*G3ppp*G4p*pow(G5p,2) + 
   72*G2px*G4pp*pow(G5p,2) - 96*G3pp*G4pp*pow(G5p,2) - 12*G2p*G4ppx*pow(G5p,2) + 12*G2pp*G4px*pow(G5p,2) - 60*G2px*G3x*G4*G5pp + 72*G3pp*G3x*G4*G5pp + 24*G2pxx*G4*G4p*G5pp - 
   24*G3ppx*G4*G4p*G5pp + 12*G2xx*pow(G4p,2)*G5pp - 84*G3px*pow(G4p,2)*G5pp - 72*G2xx*G4*G4pp*G5pp + 72*G3px*G4*G4pp*G5pp - 144*G3x*G4p*G4pp*G5pp + 144*pow(G4p,2)*G4ppx*G5pp + 
   120*G2px*G4*G4px*G5pp - 144*G3pp*G4*G4px*G5pp + 216*G4p*G4pp*G4px*G5pp - 24*G2p*G3x*G4x*G5pp - 72*G2px*G4p*G4x*G5pp + 96*G3pp*G4p*G4x*G5pp + 48*G2p*G4px*G4x*G5pp - 48*G2p*G4p*G4xx*G5pp + 
   12*G2px*G4p*G5p*G5pp - 24*G3pp*G4p*G5p*G5pp - 24*G2xx*G4*G4p*G5ppp + 24*G3px*G4*G4p*G5ppp + 108*G3x*pow(G4p,2)*G5ppp - 216*pow(G4p,2)*G4px*G5ppp + 24*G2p*G3x*G4*G5ppx + 48*G2px*G4*G4p*G5ppx - 
   48*G3pp*G4*G4p*G5ppx + 72*pow(G4p,2)*G4pp*G5ppx - 48*G2p*G4*G4px*G5ppx + 48*G2p*G4p*G4x*G5ppx - 24*G2p*G4p*G5p*G5ppx - 12*G2p*G2xx*G4*G5px + 36*G2p*G3px*G4*G5px - 24*G2pp*G3x*G4*G5px + 
   96*G2p*G3x*G4p*G5px - 48*G2ppx*G4*G4p*G5px + 48*G3ppp*G4*G4p*G5px - 36*G2px*pow(G4p,2)*G5px + 144*G3pp*pow(G4p,2)*G5px + 72*G2px*G4*G4pp*G5px - 96*G3pp*G4*G4pp*G5px + 144*G4p*pow(G4pp,2)*G5px - 
   216*pow(G4p,2)*G4ppp*G5px - 48*G2p*G4*G4ppx*G5px + 48*G2pp*G4*G4px*G5px - 240*G2p*G4p*G4px*G5px - 48*G2pp*G4p*G4x*G5px + 48*G2p*G4pp*G4x*G5px + 24*G2pp*G4p*G5p*G5px - 24*G2p*G4pp*G5p*G5px + 
   24*G2p*G4p*G5pp*G5px - 8*pow(G2px,2)*G4*G5x + 8*G2p*G2pxx*G4*G5x - 8*G2pp*G2xx*G4*G5x + 24*G2px*G3pp*G4*G5x - 16*pow(G3pp,2)*G4*G5x - 8*G2p*G3ppx*G4*G5x + 8*G2pp*G3px*G4*G5x + 6*G2p*G2xx*G4p*G5x - 
   18*G2p*G3px*G4p*G5x + 12*G2pp*G3x*G4p*G5x + 12*G2ppx*pow(G4p,2)*G5x - 12*G3ppp*pow(G4p,2)*G5x - 12*G2p*G3x*G4pp*G5x - 36*G2px*G4p*G4pp*G5x + 48*G3pp*G4p*G4pp*G5x + 24*G2p*G4p*G4ppx*G5x - 
   24*G2pp*G4p*G4px*G5x + 24*G2p*G4pp*G4px*G5x - 8*G2p*G2px*G4x*G5x + 16*G2p*G3pp*G4x*G5x + 4*G2p*G2px*G5p*G5x - 8*G2p*G3pp*G5p*G5x + 
   6*pow(G2x,2)*(4*G4*G4pxx - 4*G4p*G4xx - 2*G4px*G5p + 2*G4x*G5pp - 2*G4*G5ppx + G4p*G5px + G4pp*G5x) + 
   24*pow(G3p,2)*(4*G4*G4pxx - 4*G4p*G4xx - 2*G4px*G5p + 2*G4x*G5pp - 2*G4*G5ppx + G4p*G5px + G4pp*G5x) + 216*G3pxx*G3x*pow(G4,2)*pow(H,2) - 216*G3px*G3xx*pow(G4,2)*pow(H,2) + 
   108*G3x*G3xx*G4*G4p*pow(H,2) + 432*G3xx*pow(G4,2)*G4ppx*pow(H,2) - 720*G3x*pow(G4,2)*G4ppxx*pow(H,2) + 432*pow(G3x,2)*G4*G4px*pow(H,2) - 432*G3pxx*pow(G4,2)*G4px*pow(H,2) - 
   144*G3xx*G4*G4p*G4px*pow(H,2) + 1728*pow(G4,2)*G4ppxx*G4px*pow(H,2) - 2376*G3x*G4*pow(G4px,2)*pow(H,2) + 3600*G4*pow(G4px,3)*pow(H,2) - 144*G2xx*pow(G4,2)*G4pxx*pow(H,2) + 
   1440*G3px*pow(G4,2)*G4pxx*pow(H,2) - 360*G3x*G4*G4p*G4pxx*pow(H,2) - 2880*pow(G4,2)*G4ppx*G4pxx*pow(H,2) - 864*G4*G4p*G4px*G4pxx*pow(H,2) - 144*G4*pow(G4p,2)*G4pxxx*pow(H,2) - 
   108*G3px*G3x*G4*G4x*pow(H,2) - 1188*pow(G3x,2)*G4p*G4x*pow(H,2) + 936*G3pxx*G4*G4p*G4x*pow(H,2) + 576*G3xx*pow(G4p,2)*G4x*pow(H,2) - 504*G3xx*G4*G4pp*G4x*pow(H,2) + 
   1008*G3x*G4*G4ppx*G4x*pow(H,2) - 1296*G4*G4p*G4ppxx*G4x*pow(H,2) + 288*G2xx*G4*G4px*G4x*pow(H,2) - 432*G3px*G4*G4px*G4x*pow(H,2) + 6912*G3x*G4p*G4px*G4x*pow(H,2) - 
   1296*G4*G4ppx*G4px*G4x*pow(H,2) - 9360*G4p*pow(G4px,2)*G4x*pow(H,2) - 1728*pow(G4p,2)*G4pxx*G4x*pow(H,2) - 432*G4*G4pp*G4pxx*G4x*pow(H,2) - 432*G2pxx*G4*pow(G4x,2)*pow(H,2) + 
   432*G3ppx*G4*pow(G4x,2)*pow(H,2) - 288*G2xx*G4p*pow(G4x,2)*pow(H,2) - 504*G3px*G4p*pow(G4x,2)*pow(H,2) + 2952*G3x*G4pp*pow(G4x,2)*pow(H,2) + 4320*G4p*G4ppx*pow(G4x,2)*pow(H,2) - 
   7488*G4pp*G4px*pow(G4x,2)*pow(H,2) + 432*G2px*pow(G4x,3)*pow(H,2) - 1008*G3pp*pow(G4x,3)*pow(H,2) + 288*G4ppp*pow(G4x,3)*pow(H,2) + 144*G2pxx*pow(G4,2)*G4xx*pow(H,2) - 
   720*G3ppx*pow(G4,2)*G4xx*pow(H,2) + 144*G2xx*G4*G4p*G4xx*pow(H,2) + 144*G3px*G4*G4p*G4xx*pow(H,2) + 2448*G3x*pow(G4p,2)*G4xx*pow(H,2) - 864*G3x*G4*G4pp*G4xx*pow(H,2) + 
   1152*pow(G4,2)*G4pppx*G4xx*pow(H,2) - 144*G4*G4p*G4ppx*G4xx*pow(H,2) - 6912*pow(G4p,2)*G4px*G4xx*pow(H,2) + 3600*G4*G4pp*G4px*G4xx*pow(H,2) - 720*G3pp*G4*G4x*G4xx*pow(H,2) - 
   7488*G4p*G4pp*G4x*G4xx*pow(H,2) + 1440*G4*G4ppp*G4x*G4xx*pow(H,2) + 144*G2p*pow(G4x,2)*G4xx*pow(H,2) + 576*G2p*G4*pow(G4xx,2)*pow(H,2) - 144*G2px*pow(G4,2)*G4xxx*pow(H,2) + 
   288*G3pp*pow(G4,2)*G4xxx*pow(H,2) - 144*pow(G4p,3)*G4xxx*pow(H,2) - 144*G4*G4p*G4pp*G4xxx*pow(H,2) + 144*G2p*G4*G4x*G4xxx*pow(H,2) + 108*G3px*G3x*G4*G5p*pow(H,2) + 
   1080*pow(G3x,2)*G4p*G5p*pow(H,2) - 648*G3pxx*G4*G4p*G5p*pow(H,2) - 468*G3xx*pow(G4p,2)*G5p*pow(H,2) + 216*G3xx*G4*G4pp*G5p*pow(H,2) - 576*G3x*G4*G4ppx*G5p*pow(H,2) + 
   720*G4*G4p*G4ppxx*G5p*pow(H,2) - 144*G2xx*G4*G4px*G5p*pow(H,2) - 288*G3px*G4*G4px*G5p*pow(H,2) - 6084*G3x*G4p*G4px*G5p*pow(H,2) + 1584*G4*G4ppx*G4px*G5p*pow(H,2) + 
   8352*G4p*pow(G4px,2)*G5p*pow(H,2) + 1440*pow(G4p,2)*G4pxx*G5p*pow(H,2) + 1008*G4*G4pp*G4pxx*G5p*pow(H,2) + 720*G2pxx*G4*G4x*G5p*pow(H,2) - 576*G3ppx*G4*G4x*G5p*pow(H,2) + 
   504*G2xx*G4p*G4x*G5p*pow(H,2) + 252*G3px*G4p*G4x*G5p*pow(H,2) - 5004*G3x*G4pp*G4x*G5p*pow(H,2) - 288*G4*G4pppx*G4x*G5p*pow(H,2) - 6192*G4p*G4ppx*G4x*G5p*pow(H,2) + 
   12528*G4pp*G4px*G4x*G5p*pow(H,2) - 1080*G2px*pow(G4x,2)*G5p*pow(H,2) + 2664*G3pp*pow(G4x,2)*G5p*pow(H,2) - 1008*G4ppp*pow(G4x,2)*G5p*pow(H,2) + 432*G2px*G4*G4xx*G5p*pow(H,2) - 
   144*G3pp*G4*G4xx*G5p*pow(H,2) + 7344*G4p*G4pp*G4xx*G5p*pow(H,2) - 1440*G4*G4ppp*G4xx*G5p*pow(H,2) - 648*G2p*G4x*G4xx*G5p*pow(H,2) - 144*G2p*G4*G4xxx*G5p*pow(H,2) - 
   288*G2pxx*G4*pow(G5p,2)*pow(H,2) + 144*G3ppx*G4*pow(G5p,2)*pow(H,2) - 216*G2xx*G4p*pow(G5p,2)*pow(H,2) + 108*G3px*G4p*pow(G5p,2)*pow(H,2) + 2196*G3x*G4pp*pow(G5p,2)*pow(H,2) + 
   288*G4*G4pppx*pow(G5p,2)*pow(H,2) + 2304*G4p*G4ppx*pow(G5p,2)*pow(H,2) - 5472*G4pp*G4px*pow(G5p,2)*pow(H,2) + 864*G2px*G4x*pow(G5p,2)*pow(H,2) - 2304*G3pp*G4x*pow(G5p,2)*pow(H,2) + 
   1152*G4ppp*G4x*pow(G5p,2)*pow(H,2) + 504*G2p*G4xx*pow(G5p,2)*pow(H,2) - 216*G2px*pow(G5p,3)*pow(H,2) + 648*G3pp*pow(G5p,3)*pow(H,2) - 432*G4ppp*pow(G5p,3)*pow(H,2) - 
   540*pow(G3x,2)*G4*G5pp*pow(H,2) - 72*G3pxx*pow(G4,2)*G5pp*pow(H,2) + 36*G3xx*G4*G4p*G5pp*pow(H,2) + 144*pow(G4,2)*G4ppxx*G5pp*pow(H,2) + 2988*G3x*G4*G4px*G5pp*pow(H,2) - 
   4680*G4*pow(G4px,2)*G5pp*pow(H,2) + 1368*G4*G4p*G4pxx*G5pp*pow(H,2) - 360*G2xx*G4*G4x*G5pp*pow(H,2) - 1800*G3x*G4p*G4x*G5pp*pow(H,2) + 720*G4*G4ppx*G4x*G5pp*pow(H,2) + 
   4608*G4p*G4px*G4x*G5pp*pow(H,2) + 576*G4pp*pow(G4x,2)*G5pp*pow(H,2) + 864*pow(G4p,2)*G4xx*G5pp*pow(H,2) - 1152*G4*G4pp*G4xx*G5pp*pow(H,2) + 216*G2xx*G4*G5p*G5pp*pow(H,2) + 
   432*G3px*G4*G5p*G5pp*pow(H,2) + 1404*G3x*G4p*G5p*G5pp*pow(H,2) - 1440*G4*G4ppx*G5p*G5pp*pow(H,2) - 4140*G4p*G4px*G5p*G5pp*pow(H,2) - 1224*G4pp*G4x*G5p*G5pp*pow(H,2) + 
   648*G4pp*pow(G5p,2)*G5pp*pow(H,2) - 72*G3x*G4*pow(G5pp,2)*pow(H,2) + 432*G4*G4px*pow(G5pp,2)*pow(H,2) - 216*G4p*G4x*pow(G5pp,2)*pow(H,2) + 432*G4p*G5p*pow(G5pp,2)*pow(H,2) + 
   72*G3xx*pow(G4,2)*G5ppp*pow(H,2) - 144*pow(G4,2)*G4pxx*G5ppp*pow(H,2) + 216*G3x*G4*G4x*G5ppp*pow(H,2) - 432*G4*G4px*G4x*G5ppp*pow(H,2) - 432*G4p*pow(G4x,2)*G5ppp*pow(H,2) - 
   864*G4*G4p*G4xx*G5ppp*pow(H,2) - 360*G3x*G4*G5p*G5ppp*pow(H,2) + 864*G4*G4px*G5p*G5ppp*pow(H,2) + 1080*G4p*G4x*G5p*G5ppp*pow(H,2) - 648*G4p*pow(G5p,2)*G5ppp*pow(H,2) + 
   144*G3x*pow(G4,2)*G5pppx*pow(H,2) - 432*pow(G4,2)*G4px*G5pppx*pow(H,2) - 288*G4*G4p*G4x*G5pppx*pow(H,2) + 288*G4*G4p*G5p*G5pppx*pow(H,2) + 56*G2xx*pow(G4,2)*G5ppx*pow(H,2) - 
   560*G3px*pow(G4,2)*G5ppx*pow(H,2) + 48*G3x*G4*G4p*G5ppx*pow(H,2) + 1152*pow(G4,2)*G4ppx*G5ppx*pow(H,2) + 696*G4*G4p*G4px*G5ppx*pow(H,2) + 48*pow(G4p,2)*G4x*G5ppx*pow(H,2) + 
   864*G4*G4pp*G4x*G5ppx*pow(H,2) - 24*pow(G4p,2)*G5p*G5ppx*pow(H,2) - 864*G4*G4pp*G5p*G5ppx*pow(H,2) - 792*G4*G4p*G5pp*G5ppx*pow(H,2) + 48*G4*pow(G4p,2)*G5ppxx*pow(H,2) - 
   56*G2pxx*pow(G4,2)*G5px*pow(H,2) + 416*G3ppx*pow(G4,2)*G5px*pow(H,2) - 84*G2xx*G4*G4p*G5px*pow(H,2) - 360*G3px*G4*G4p*G5px*pow(H,2) - 1446*G3x*pow(G4p,2)*G5px*pow(H,2) + 
   636*G3x*G4*G4pp*G5px*pow(H,2) - 720*pow(G4,2)*G4pppx*G5px*pow(H,2) + 600*G4*G4p*G4ppx*G5px*pow(H,2) + 4308*pow(G4p,2)*G4px*G5px*pow(H,2) - 2784*G4*G4pp*G4px*G5px*pow(H,2) - 
   104*G2px*G4*G4x*G5px*pow(H,2) + 688*G3pp*G4*G4x*G5px*pow(H,2) + 4080*G4p*G4pp*G4x*G5px*pow(H,2) - 864*G4*G4ppp*G4x*G5px*pow(H,2) + 48*G2p*pow(G4x,2)*G5px*pow(H,2) - 
   792*G2p*G4*G4xx*G5px*pow(H,2) - 248*G2px*G4*G5p*G5px*pow(H,2) - 32*G3pp*G4*G5p*G5px*pow(H,2) - 4164*G4p*G4pp*G5p*G5px*pow(H,2) + 864*G4*G4ppp*G5p*G5px*pow(H,2) + 
   144*G2p*G4x*G5p*G5px*pow(H,2) - 204*G2p*pow(G5p,2)*G5px*pow(H,2) - 708*pow(G4p,2)*G5pp*G5px*pow(H,2) + 792*G4*G4pp*G5pp*G5px*pow(H,2) + 504*G4*G4p*G5ppp*G5px*pow(H,2) + 
   252*G2p*G4*pow(G5px,2)*pow(H,2) + 136*G2px*pow(G4,2)*G5pxx*pow(H,2) - 216*G3pp*pow(G4,2)*G5pxx*pow(H,2) + 84*pow(G4p,3)*G5pxx*pow(H,2) + 288*G4*G4p*G4pp*G5pxx*pow(H,2) - 
   8*G2p*G4*G4x*G5pxx*pow(H,2) + 64*G2p*G4*G5p*G5pxx*pow(H,2) - 8*G2p*pow(G4,2)*G5pxxx*pow(H,2) - 108*G2px*G3x*G4*G5x*pow(H,2) + 72*G3pp*G3x*G4*G5x*pow(H,2) + 54*G2p*G3xx*G4*G5x*pow(H,2) + 
   168*G2pxx*G4*G4p*G5x*pow(H,2) + 48*G3ppx*G4*G4p*G5x*pow(H,2) + 108*G2xx*pow(G4p,2)*G5x*pow(H,2) - 18*G3px*pow(G4p,2)*G5x*pow(H,2) - 60*G2xx*G4*G4pp*G5x*pow(H,2) - 
   300*G3px*G4*G4pp*G5x*pow(H,2) - 1620*G3x*G4p*G4pp*G5x*pow(H,2) + 288*G3x*G4*G4ppp*G5x*pow(H,2) - 432*G4*G4p*G4pppx*G5x*pow(H,2) - 900*pow(G4p,2)*G4ppx*G5x*pow(H,2) + 
   936*G4*G4pp*G4ppx*G5x*pow(H,2) + 240*G2px*G4*G4px*G5x*pow(H,2) + 12*G3pp*G4*G4px*G5x*pow(H,2) + 4032*G4p*G4pp*G4px*G5x*pow(H,2) - 792*G4*G4ppp*G4px*G5x*pow(H,2) - 
   60*G2p*G4*G4pxx*G5x*pow(H,2) + 90*G2p*G3x*G4x*G5x*pow(H,2) - 732*G2px*G4p*G4x*G5x*pow(H,2) + 1860*G3pp*G4p*G4x*G5x*pow(H,2) + 648*pow(G4pp,2)*G4x*G5x*pow(H,2) - 
   1080*G4p*G4ppp*G4x*G5x*pow(H,2) + 48*G2p*G4px*G4x*G5x*pow(H,2) - 96*G2pp*pow(G4x,2)*G5x*pow(H,2) - 48*G2pp*G4*G4xx*G5x*pow(H,2) - 360*G2p*G4p*G4xx*G5x*pow(H,2) - 
   162*G2p*G3x*G5p*G5x*pow(H,2) + 48*G2ppx*G4*G5p*G5x*pow(H,2) - 48*G3ppp*G4*G5p*G5x*pow(H,2) + 540*G2px*G4p*G5p*G5x*pow(H,2) - 1512*G3pp*G4p*G5p*G5x*pow(H,2) - 
   576*pow(G4pp,2)*G5p*G5x*pow(H,2) + 1152*G4p*G4ppp*G5p*G5x*pow(H,2) + 198*G2p*G4px*G5p*G5x*pow(H,2) + 192*G2pp*G4x*G5p*G5x*pow(H,2) - 84*G2pp*pow(G5p,2)*G5x*pow(H,2) - 
   108*G2px*G4*G5pp*G5x*pow(H,2) + 120*G3pp*G4*G5pp*G5x*pow(H,2) - 432*G4p*G4pp*G5pp*G5x*pow(H,2) - 120*G2p*G4x*G5pp*G5x*pow(H,2) + 72*G2p*G5p*G5pp*G5x*pow(H,2) + 
   288*pow(G4p,2)*G5ppp*G5x*pow(H,2) - 48*G2p*G4*G5ppx*G5x*pow(H,2) + 48*G2pp*G4*G5px*G5x*pow(H,2) + 162*G2p*G4p*G5px*G5x*pow(H,2) + 36*G2pp*G4p*pow(G5x,2)*pow(H,2) - 
   36*G2p*G4pp*pow(G5x,2)*pow(H,2) + 78*G2p*G3x*G4*G5xx*pow(H,2) - 56*G2ppx*pow(G4,2)*G5xx*pow(H,2) + 56*G3ppp*pow(G4,2)*G5xx*pow(H,2) + 96*G2px*G4*G4p*G5xx*pow(H,2) - 
   108*G3pp*G4*G4p*G5xx*pow(H,2) + 492*pow(G4p,2)*G4pp*G5xx*pow(H,2) - 168*G4*pow(G4pp,2)*G5xx*pow(H,2) - 168*G4*G4p*G4ppp*G5xx*pow(H,2) - 240*G2p*G4*G4px*G5xx*pow(H,2) - 
   64*G2pp*G4*G4x*G5xx*pow(H,2) - 84*G2p*G4p*G4x*G5xx*pow(H,2) + 8*G2pp*G4*G5p*G5xx*pow(H,2) + 126*G2p*G4p*G5p*G5xx*pow(H,2) + 48*G2p*G4*G5pp*G5xx*pow(H,2) + 
   8*G2pp*pow(G4,2)*G5xxx*pow(H,2) - 12*G2p*G4*G4p*G5xxx*pow(H,2) + 864*pow(G4,2)*G4pxxx*G4x*pow(H,4) - 6048*G4*G4pxx*pow(G4x,2)*pow(H,4) - 1728*G4px*pow(G4x,3)*pow(H,4) + 
   6912*G4*G4px*G4x*G4xx*pow(H,4) - 9504*G4p*pow(G4x,2)*G4xx*pow(H,4) + 3456*G4*G4p*pow(G4xx,2)*pow(H,4) - 864*pow(G4,2)*G4px*G4xxx*pow(H,4) + 864*G4*G4p*G4x*G4xxx*pow(H,4) - 
   864*pow(G4,2)*G4pxxx*G5p*pow(H,4) + 9504*G4*G4pxx*G4x*G5p*pow(H,4) + 4752*G4px*pow(G4x,2)*G5p*pow(H,4) - 4320*G4*G4px*G4xx*G5p*pow(H,4) + 16848*G4p*G4x*G4xx*G5p*pow(H,4) - 
   864*G4*G4p*G4xxx*G5p*pow(H,4) - 3456*G4*G4pxx*pow(G5p,2)*pow(H,4) - 4320*G4px*G4x*pow(G5p,2)*pow(H,4) - 7344*G4p*G4xx*pow(G5p,2)*pow(H,4) + 1296*G4px*pow(G5p,3)*pow(H,4) + 
   2592*pow(G4x,3)*G5pp*pow(H,4) - 9072*G4*G4x*G4xx*G5pp*pow(H,4) + 864*pow(G4,2)*G4xxx*G5pp*pow(H,4) - 7128*pow(G4x,2)*G5p*G5pp*pow(H,4) + 6480*G4*G4xx*G5p*G5pp*pow(H,4) + 
   6480*G4x*pow(G5p,2)*G5pp*pow(H,4) - 1944*pow(G5p,3)*G5pp*pow(H,4) + 4176*G4*pow(G4x,2)*G5ppx*pow(H,4) - 816*pow(G4,2)*G4xx*G5ppx*pow(H,4) - 6528*G4*G4x*G5p*G5ppx*pow(H,4) + 
   2352*G4*pow(G5p,2)*G5ppx*pow(H,4) - 480*pow(G4,2)*G4x*G5ppxx*pow(H,4) + 480*pow(G4,2)*G5p*G5ppxx*pow(H,4) - 216*G3xx*pow(G4,2)*G5px*pow(H,4) + 1536*pow(G4,2)*G4pxx*G5px*pow(H,4) + 
   1476*G3x*G4*G4x*G5px*pow(H,4) - 8160*G4*G4px*G4x*G5px*pow(H,4) + 7296*G4p*pow(G4x,2)*G5px*pow(H,4) - 5616*G4*G4p*G4xx*G5px*pow(H,4) - 1476*G3x*G4*G5p*G5px*pow(H,4) + 
   5904*G4*G4px*G5p*G5px*pow(H,4) - 13116*G4p*G4x*G5p*G5px*pow(H,4) + 5676*G4p*pow(G5p,2)*G5px*pow(H,4) + 5208*G4*G4x*G5pp*G5px*pow(H,4) - 3240*G4*G5p*G5pp*G5px*pow(H,4) - 
   144*pow(G4,2)*G5ppx*G5px*pow(H,4) + 1764*G4*G4p*pow(G5px,2)*pow(H,4) + 360*G3x*pow(G4,2)*G5pxx*pow(H,4) - 96*pow(G4,2)*G4px*G5pxx*pow(H,4) + 960*G4*G4p*G4x*G5pxx*pow(H,4) - 
   288*G4*G4p*G5p*G5pxx*pow(H,4) - 648*pow(G4,2)*G5pp*G5pxx*pow(H,4) - 96*pow(G4,2)*G4p*G5pxxx*pow(H,4) + 216*G3pxx*pow(G4,2)*G5x*pow(H,4) + 324*G3xx*G4*G4p*G5x*pow(H,4) - 
   720*pow(G4,2)*G4ppxx*G5x*pow(H,4) + 1080*G3x*G4*G4px*G5x*pow(H,4) - 2952*G4*pow(G4px,2)*G5x*pow(H,4) + 1800*G4*G4p*G4pxx*G5x*pow(H,4) - 1404*G3px*G4*G4x*G5x*pow(H,4) - 
   3996*G3x*G4p*G4x*G5x*pow(H,4) + 4608*G4*G4ppx*G4x*G5x*pow(H,4) + 13896*G4p*G4px*G4x*G5x*pow(H,4) + 5760*G4pp*pow(G4x,2)*G5x*pow(H,4) + 3312*pow(G4p,2)*G4xx*G5x*pow(H,4) - 
   1728*G4*G4pp*G4xx*G5x*pow(H,4) + 1404*G3px*G4*G5p*G5x*pow(H,4) + 3564*G3x*G4p*G5p*G5x*pow(H,4) - 4176*G4*G4ppx*G5p*G5x*pow(H,4) - 12636*G4p*G4px*G5p*G5x*pow(H,4) - 
   10620*G4pp*G4x*G5p*G5x*pow(H,4) + 5004*G4pp*pow(G5p,2)*G5x*pow(H,4) - 1512*G3x*G4*G5pp*G5x*pow(H,4) + 3924*G4*G4px*G5pp*G5x*pow(H,4) - 3708*G4p*G4x*G5pp*G5x*pow(H,4) + 
   3528*G4p*G5p*G5pp*G5x*pow(H,4) + 216*G4*pow(G5pp,2)*G5x*pow(H,4) + 504*G4*G4x*G5ppp*G5x*pow(H,4) - 648*G4*G5p*G5ppp*G5x*pow(H,4) + 144*pow(G4,2)*G5pppx*G5x*pow(H,4) - 
   1512*G4*G4p*G5ppx*G5x*pow(H,4) - 2358*pow(G4p,2)*G5px*G5x*pow(H,4) + 972*G4*G4pp*G5px*G5x*pow(H,4) - 144*G2px*G4*pow(G5x,2)*pow(H,4) + 180*G3pp*G4*pow(G5x,2)*pow(H,4) - 
   1818*G4p*G4pp*pow(G5x,2)*pow(H,4) + 216*G4*G4ppp*pow(G5x,2)*pow(H,4) + 252*G2p*G4x*pow(G5x,2)*pow(H,4) - 306*G2p*G5p*pow(G5x,2)*pow(H,4) - 360*G3px*pow(G4,2)*G5xx*pow(H,4) + 
   468*G3x*G4*G4p*G5xx*pow(H,4) + 576*pow(G4,2)*G4ppx*G5xx*pow(H,4) - 1272*G4*G4p*G4px*G5xx*pow(H,4) + 1488*pow(G4p,2)*G4x*G5xx*pow(H,4) - 1248*G4*G4pp*G4x*G5xx*pow(H,4) - 
   1236*pow(G4p,2)*G5p*G5xx*pow(H,4) + 576*G4*G4pp*G5p*G5xx*pow(H,4) + 372*G4*G4p*G5pp*G5xx*pow(H,4) + 168*pow(G4,2)*G5ppp*G5xx*pow(H,4) + 162*G2p*G4*G5x*G5xx*pow(H,4) - 
   72*G4*pow(G4p,2)*G5xxx*pow(H,4) + 96*pow(G4,2)*G4pp*G5xxx*pow(H,4) + 900*G4*G4x*G5px*G5x*pow(H,6) - 900*G4*G5p*G5px*G5x*pow(H,6) + 360*pow(G4,2)*G5pxx*G5x*pow(H,6) + 
   216*G4*G4px*pow(G5x,2)*pow(H,6) - 2160*G4p*G4x*pow(G5x,2)*pow(H,6) + 1836*G4p*G5p*pow(G5x,2)*pow(H,6) - 540*G4*G5pp*pow(G5x,2)*pow(H,6) - 360*pow(G4,2)*G5px*G5xx*pow(H,6) + 
   972*G4*G4p*G5x*G5xx*pow(H,6) - 18*G3pxx*G3x*G4*rptot + 18*G3px*G3xx*G4*rptot - 18*G3x*G3xx*G4p*rptot - 36*G3xx*G4*G4ppx*rptot + 36*G3x*G4*G4ppxx*rptot + 18*pow(G3x,2)*G4px*rptot + 
   36*G3pxx*G4*G4px*rptot - 72*G4*G4ppxx*G4px*rptot - 108*G3x*pow(G4px,2)*rptot + 144*pow(G4px,3)*rptot + 24*G2xx*G4*G4pxx*rptot - 60*G3px*G4*G4pxx*rptot + 108*G3x*G4p*G4pxx*rptot + 
   72*G4*G4ppx*G4pxx*rptot - 144*G4p*G4px*G4pxx*rptot - 36*G3px*G3x*G4x*rptot - 36*G3pxx*G4p*G4x*rptot + 36*G3xx*G4pp*G4x*rptot + 36*G3x*G4ppx*G4x*rptot + 72*G4p*G4ppxx*G4x*rptot + 
   108*G3px*G4px*G4x*rptot - 144*G4ppx*G4px*G4x*rptot - 72*G4pp*G4pxx*G4x*rptot + 24*G2pxx*pow(G4x,2)*rptot - 24*G3ppx*pow(G4x,2)*rptot - 24*G2pxx*G4*G4xx*rptot + 24*G3ppx*G4*G4xx*rptot - 
   24*G2xx*G4p*G4xx*rptot + 96*G3px*G4p*G4xx*rptot + 72*G3x*G4pp*G4xx*rptot - 144*G4p*G4ppx*G4xx*rptot - 144*G4pp*G4px*G4xx*rptot + 24*G2px*G4x*G4xx*rptot - 48*G3pp*G4x*G4xx*rptot + 
   36*G3px*G3x*G5p*rptot + 18*G3pxx*G4p*G5p*rptot - 18*G3xx*G4pp*G5p*rptot - 54*G3x*G4ppx*G5p*rptot - 36*G4p*G4ppxx*G5p*rptot - 12*G2xx*G4px*G5p*rptot - 78*G3px*G4px*G5p*rptot + 
   144*G4ppx*G4px*G5p*rptot + 36*G4pp*G4pxx*G5p*rptot - 36*G2pxx*G4x*G5p*rptot + 36*G3ppx*G4x*G5p*rptot - 12*G2px*G4xx*G5p*rptot + 24*G3pp*G4xx*G5p*rptot + 12*G2pxx*pow(G5p,2)*rptot - 
   12*G3ppx*pow(G5p,2)*rptot + 18*G3xx*G4p*G5pp*rptot + 18*G3x*G4px*G5pp*rptot - 36*pow(G4px,2)*G5pp*rptot - 36*G4p*G4pxx*G5pp*rptot + 12*G2xx*G4x*G5pp*rptot - 12*G3px*G4x*G5pp*rptot - 
   12*G2xx*G4*G5ppx*rptot + 12*G3px*G4*G5ppx*rptot - 36*G3x*G4p*G5ppx*rptot + 72*G4p*G4px*G5ppx*rptot + 12*G2pxx*G4*G5px*rptot - 12*G3ppx*G4*G5px*rptot + 6*G2xx*G4p*G5px*rptot - 42*G3px*G4p*G5px*rptot - 
   36*G3x*G4pp*G5px*rptot + 72*G4p*G4ppx*G5px*rptot + 72*G4pp*G4px*G5px*rptot - 12*G2px*G4x*G5px*rptot + 24*G3pp*G4x*G5px*rptot + 6*G2px*G5p*G5px*rptot - 12*G3pp*G5p*G5px*rptot + 3*G2px*G3x*G5x*rptot - 
   6*G3pp*G3x*G5x*rptot - 6*G2pxx*G4p*G5x*rptot + 6*G3ppx*G4p*G5x*rptot + 6*G2xx*G4pp*G5x*rptot - 6*G3px*G4pp*G5x*rptot - 6*G2px*G4px*G5x*rptot + 12*G3pp*G4px*G5x*rptot - 
   144*G4*G4pxxx*G4x*pow(H,2)*rptot + 144*G4pxx*pow(G4x,2)*pow(H,2)*rptot - 288*G4px*G4x*G4xx*pow(H,2)*rptot - 576*G4p*pow(G4xx,2)*pow(H,2)*rptot + 144*G4*G4px*G4xxx*pow(H,2)*rptot - 
   144*G4p*G4x*G4xxx*pow(H,2)*rptot + 144*G4*G4pxxx*G5p*pow(H,2)*rptot - 72*G4pxx*G4x*G5p*pow(H,2)*rptot + 72*G4px*G4xx*G5p*pow(H,2)*rptot + 144*G4p*G4xxx*G5p*pow(H,2)*rptot - 
   72*G4pxx*pow(G5p,2)*pow(H,2)*rptot + 648*G4x*G4xx*G5pp*pow(H,2)*rptot - 144*G4*G4xxx*G5pp*pow(H,2)*rptot - 432*G4xx*G5p*G5pp*pow(H,2)*rptot - 144*pow(G4x,2)*G5ppx*pow(H,2)*rptot + 
   72*G4*G4xx*G5ppx*pow(H,2)*rptot + 144*G4x*G5p*G5ppx*pow(H,2)*rptot + 72*G4*G4x*G5ppxx*pow(H,2)*rptot - 72*G4*G5p*G5ppxx*pow(H,2)*rptot + 54*G3xx*G4*G5px*pow(H,2)*rptot - 
   180*G4*G4pxx*G5px*pow(H,2)*rptot - 90*G3x*G4x*G5px*pow(H,2)*rptot + 360*G4px*G4x*G5px*pow(H,2)*rptot + 720*G4p*G4xx*G5px*pow(H,2)*rptot + 54*G3x*G5p*G5px*pow(H,2)*rptot - 
   90*G4px*G5p*G5px*pow(H,2)*rptot - 360*G4x*G5pp*G5px*pow(H,2)*rptot + 216*G5p*G5pp*G5px*pow(H,2)*rptot - 198*G4p*pow(G5px,2)*pow(H,2)*rptot - 6*G3x*G4*G5pxx*pow(H,2)*rptot - 
   96*G4*G4px*G5pxx*pow(H,2)*rptot + 60*G4p*G4x*G5pxx*pow(H,2)*rptot - 102*G4p*G5p*G5pxx*pow(H,2)*rptot + 72*G4*G5pp*G5pxx*pow(H,2)*rptot + 12*G4*G4p*G5pxxx*pow(H,2)*rptot - 
   54*G3pxx*G4*G5x*pow(H,2)*rptot - 54*G3xx*G4p*G5x*pow(H,2)*rptot + 108*G4*G4ppxx*G5x*pow(H,2)*rptot - 18*G3x*G4px*G5x*pow(H,2)*rptot + 180*G4p*G4pxx*G5x*pow(H,2)*rptot + 
   18*G3px*G4x*G5x*pow(H,2)*rptot - 216*G4ppx*G4x*G5x*pow(H,2)*rptot + 360*G4pp*G4xx*G5x*pow(H,2)*rptot + 18*G3px*G5p*G5x*pow(H,2)*rptot + 90*G4ppx*G5p*G5x*pow(H,2)*rptot + 
   90*G3x*G5pp*G5x*pow(H,2)*rptot - 198*G4px*G5pp*G5x*pow(H,2)*rptot - 18*G4p*G5ppx*G5x*pow(H,2)*rptot - 198*G4pp*G5px*G5x*pow(H,2)*rptot + 9*G2px*pow(G5x,2)*pow(H,2)*rptot - 
   18*G3pp*pow(G5x,2)*pow(H,2)*rptot + 6*G3px*G4*G5xx*pow(H,2)*rptot - 78*G3x*G4p*G5xx*pow(H,2)*rptot + 24*G4*G4ppx*G5xx*pow(H,2)*rptot + 180*G4p*G4px*G5xx*pow(H,2)*rptot + 
   156*G4pp*G4x*G5xx*pow(H,2)*rptot - 114*G4pp*G5p*G5xx*pow(H,2)*rptot - 30*G4p*G5pp*G5xx*pow(H,2)*rptot + 12*pow(G4p,2)*G5xxx*pow(H,2)*rptot - 12*G4*G4pp*G5xxx*pow(H,2)*rptot - 
   180*G4x*G5px*G5x*pow(H,4)*rptot + 180*G5p*G5px*G5x*pow(H,4)*rptot - 90*G4*G5pxx*G5x*pow(H,4)*rptot - 72*G4px*pow(G5x,2)*pow(H,4)*rptot + 126*G5pp*pow(G5x,2)*pow(H,4)*rptot + 
   90*G4*G5px*G5xx*pow(H,4)*rptot - 162*G4p*G5x*G5xx*pow(H,4)*rptot - 2*G3p*
    (-72*pow(G3x,2)*G4p + 36*G3pxx*G4*G4p + 36*G3xx*pow(G4p,2) - 36*G3xx*G4*G4pp - 72*G4*G4p*G4ppxx + 24*G2xx*G4*G4px + 312*G4*G4ppx*G4px - 48*G4p*pow(G4px,2) - 144*pow(G4p,2)*G4pxx + 
      120*G4*G4pp*G4pxx - 72*G2pxx*G4*G4x + 120*G3ppx*G4*G4x - 24*G2xx*G4p*G4x - 96*G4*G4pppx*G4x - 144*G4p*G4ppx*G4x - 576*G4pp*G4px*G4x + 24*G2px*pow(G4x,2) - 24*G3pp*pow(G4x,2) - 
      48*G4ppp*pow(G4x,2) - 24*G2px*G4*G4xx + 24*G3pp*G4*G4xx - 192*G4p*G4pp*G4xx + 48*G4*G4ppp*G4xx - 24*G2p*G4x*G4xx + 48*G2pxx*G4*G5p - 72*G3ppx*G4*G5p + 24*G2xx*G4p*G5p + 48*G4*G4pppx*G5p + 
      144*G4p*G4ppx*G5p + 336*G4pp*G4px*G5p - 36*G2px*G4x*G5p + 36*G3pp*G4x*G5p + 72*G4ppp*G4x*G5p + 12*G2p*G4xx*G5p + 12*G2px*pow(G5p,2) - 12*G3pp*pow(G5p,2) - 24*G4ppp*pow(G5p,2) - 
      24*G2xx*G4*G5pp - 24*G4*G4ppx*G5pp - 108*G4p*G4px*G5pp + 48*G4pp*G4x*G5pp - 12*G4pp*G5p*G5pp + 12*G4p*pow(G5pp,2) - 24*G4*G4px*G5ppp + 24*G4p*G4x*G5ppp - 12*G4p*G5p*G5ppp + 
      36*pow(G4p,2)*G5ppx - 24*G4*G4pp*G5ppx + 12*G2px*G4*G5px - 12*G3pp*G4*G5px + 84*G4p*G4pp*G5px - 24*G4*G4ppp*G5px + 12*G2p*G4x*G5px - 6*G2p*G5p*G5px + 8*G2ppx*G4*G5x - 8*G3ppp*G4*G5x - 
      6*G2px*G4p*G5x + 6*G3pp*G4p*G5x + 12*pow(G4pp,2)*G5x + 12*G4p*G4ppp*G5x + 6*G2p*G4px*G5x + 8*G2pp*G4x*G5x - 4*G2pp*G5p*G5x + 144*pow(G4,2)*G4pxxx*pow(H,2) - 576*G4*G4pxx*G4x*pow(H,2) - 
      288*G4px*pow(G4x,2)*pow(H,2) + 1008*G4*G4px*G4xx*pow(H,2) - 1728*G4p*G4x*G4xx*pow(H,2) + 144*G4*G4pxx*G5p*pow(H,2) + 288*G4px*G4x*G5p*pow(H,2) + 1656*G4p*G4xx*G5p*pow(H,2) + 
      648*pow(G4x,2)*G5pp*pow(H,2) - 1080*G4*G4xx*G5pp*pow(H,2) - 972*G4x*G5p*G5pp*pow(H,2) + 324*pow(G5p,2)*G5pp*pow(H,2) + 512*G4*G4x*G5ppx*pow(H,2) - 208*G4*G5p*G5ppx*pow(H,2) - 
      80*pow(G4,2)*G5ppxx*pow(H,2) - 1224*G4*G4px*G5px*pow(H,2) + 1092*G4p*G4x*G5px*pow(H,2) - 978*G4p*G5p*G5px*pow(H,2) + 672*G4*G5pp*G5px*pow(H,2) + 216*G4*G4ppx*G5x*pow(H,2) + 
      816*G4p*G4px*G5x*pow(H,2) + 732*G4pp*G4x*G5x*pow(H,2) - 558*G4pp*G5p*G5x*pow(H,2) - 186*G4p*G5pp*G5x*pow(H,2) + 12*G4*G5ppp*G5x*pow(H,2) - 9*G2p*pow(G5x,2)*pow(H,2) + 
      180*pow(G4p,2)*G5xx*pow(H,2) - 144*G4*G4pp*G5xx*pow(H,2) + 270*G4*G5px*G5x*pow(H,4) - 378*G4p*pow(G5x,2)*pow(H,4) - 24*G4pxx*G4x*rptot + 24*G4px*G4xx*rptot + 12*G4pxx*G5p*rptot - 
      12*G4xx*G5pp*rptot + 12*G4x*G5ppx*rptot - 6*G5p*G5ppx*rptot - 18*G4px*G5px*rptot + 6*G5pp*G5px*rptot + 6*G4ppx*G5x*rptot + 
      3*G3x*(76*G4p*G4px + 76*G4pp*G4x - 50*G4pp*G5p + 4*G4p*G5pp - G2p*G5x - 114*G4p*G5x*pow(H,2) + G4*(-32*G4ppx + 4*G5ppp + 74*G5px*pow(H,2)) + G5px*rptot) + 
      3*G3px*(18*G3x*G4 - 64*G4*G4px + 36*G4p*G4x - 34*G4p*G5p + 12*G4*G5pp - 6*G4*G5x*pow(H,2) - G5x*rptot)) + 
   G2x*(-72*pow(G3x,2)*G4p + 36*G3pxx*G4*G4p + 36*G3xx*pow(G4p,2) - 36*G3xx*G4*G4pp - 72*G4*G4p*G4ppxx + 24*G2xx*G4*G4px + 312*G4*G4ppx*G4px - 48*G4p*pow(G4px,2) - 96*G3p*G4*G4pxx - 
      144*pow(G4p,2)*G4pxx + 120*G4*G4pp*G4pxx - 72*G2pxx*G4*G4x + 120*G3ppx*G4*G4x - 24*G2xx*G4p*G4x - 96*G4*G4pppx*G4x - 144*G4p*G4ppx*G4x - 576*G4pp*G4px*G4x + 24*G2px*pow(G4x,2) - 
      24*G3pp*pow(G4x,2) - 48*G4ppp*pow(G4x,2) - 24*G2px*G4*G4xx + 24*G3pp*G4*G4xx + 96*G3p*G4p*G4xx - 192*G4p*G4pp*G4xx + 48*G4*G4ppp*G4xx - 24*G2p*G4x*G4xx + 48*G2pxx*G4*G5p - 72*G3ppx*G4*G5p + 
      24*G2xx*G4p*G5p + 48*G4*G4pppx*G5p + 144*G4p*G4ppx*G5p + 48*G3p*G4px*G5p + 336*G4pp*G4px*G5p - 36*G2px*G4x*G5p + 36*G3pp*G4x*G5p + 72*G4ppp*G4x*G5p + 12*G2p*G4xx*G5p + 12*G2px*pow(G5p,2) - 
      12*G3pp*pow(G5p,2) - 24*G4ppp*pow(G5p,2) - 24*G2xx*G4*G5pp - 24*G4*G4ppx*G5pp - 108*G4p*G4px*G5pp - 48*G3p*G4x*G5pp + 48*G4pp*G4x*G5pp - 12*G4pp*G5p*G5pp + 12*G4p*pow(G5pp,2) - 
      24*G4*G4px*G5ppp + 24*G4p*G4x*G5ppp - 12*G4p*G5p*G5ppp + 48*G3p*G4*G5ppx + 36*pow(G4p,2)*G5ppx - 24*G4*G4pp*G5ppx + 12*G2px*G4*G5px - 12*G3pp*G4*G5px - 24*G3p*G4p*G5px + 84*G4p*G4pp*G5px - 
      24*G4*G4ppp*G5px + 12*G2p*G4x*G5px - 6*G2p*G5p*G5px + 8*G2ppx*G4*G5x - 8*G3ppp*G4*G5x - 6*G2px*G4p*G5x + 6*G3pp*G4p*G5x - 24*G3p*G4pp*G5x + 12*pow(G4pp,2)*G5x + 12*G4p*G4ppp*G5x + 
      6*G2p*G4px*G5x + 8*G2pp*G4x*G5x - 4*G2pp*G5p*G5x + 144*pow(G4,2)*G4pxxx*pow(H,2) - 576*G4*G4pxx*G4x*pow(H,2) - 288*G4px*pow(G4x,2)*pow(H,2) + 1008*G4*G4px*G4xx*pow(H,2) - 
      1728*G4p*G4x*G4xx*pow(H,2) + 144*G4*G4pxx*G5p*pow(H,2) + 288*G4px*G4x*G5p*pow(H,2) + 1656*G4p*G4xx*G5p*pow(H,2) + 648*pow(G4x,2)*G5pp*pow(H,2) - 1080*G4*G4xx*G5pp*pow(H,2) - 
      972*G4x*G5p*G5pp*pow(H,2) + 324*pow(G5p,2)*G5pp*pow(H,2) + 512*G4*G4x*G5ppx*pow(H,2) - 208*G4*G5p*G5ppx*pow(H,2) - 80*pow(G4,2)*G5ppxx*pow(H,2) - 1224*G4*G4px*G5px*pow(H,2) + 
      1092*G4p*G4x*G5px*pow(H,2) - 978*G4p*G5p*G5px*pow(H,2) + 672*G4*G5pp*G5px*pow(H,2) + 216*G4*G4ppx*G5x*pow(H,2) + 816*G4p*G4px*G5x*pow(H,2) + 732*G4pp*G4x*G5x*pow(H,2) - 
      558*G4pp*G5p*G5x*pow(H,2) - 186*G4p*G5pp*G5x*pow(H,2) + 12*G4*G5ppp*G5x*pow(H,2) - 9*G2p*pow(G5x,2)*pow(H,2) + 180*pow(G4p,2)*G5xx*pow(H,2) - 144*G4*G4pp*G5xx*pow(H,2) + 
      270*G4*G5px*G5x*pow(H,4) - 378*G4p*pow(G5x,2)*pow(H,4) - 24*G4pxx*G4x*rptot + 24*G4px*G4xx*rptot + 12*G4pxx*G5p*rptot - 12*G4xx*G5pp*rptot + 12*G4x*G5ppx*rptot - 6*G5p*G5ppx*rptot - 
      18*G4px*G5px*rptot + 6*G5pp*G5px*rptot + 6*G4ppx*G5x*rptot + 3*G3x*(76*G4p*G4px + 76*G4pp*G4x - 50*G4pp*G5p + 4*G4p*G5pp - G2p*G5x - 114*G4p*G5x*pow(H,2) + 
         G4*(-32*G4ppx + 4*G5ppp + 74*G5px*pow(H,2)) + G5px*rptot) + 3*G3px*(18*G3x*G4 - 64*G4*G4px + 36*G4p*G4x - 34*G4p*G5p + 12*G4*G5pp - 6*G4*G5x*pow(H,2) - G5x*rptot))) + 
2*pow(a,10)*pow(dphi,7)*H*(18*G2x*G3pxx*G3x*G4 - 36*G3p*G3pxx*G3x*G4 + 18*G2pxx*pow(G3x,2)*G4 + 18*G3pp*G3x*G3xx*G4 + 18*G2x*G3x*G3xx*G4p - 36*G3p*G3x*G3xx*G4p - 45*pow(G3x,3)*G4pp + 
   36*G3pxx*G3x*G4*G4pp + 36*G3x*G3xx*G4p*G4pp - 36*G3x*G3xx*G4*G4ppp - 36*pow(G3x,2)*G4*G4pppx + 36*G2x*G3xx*G4*G4ppx - 72*G3p*G3xx*G4*G4ppx - 18*pow(G3x,2)*G4p*G4ppx + 72*G3xx*G4*G4pp*G4ppx - 
   144*G3x*G4*pow(G4ppx,2) - 36*G2x*G3x*G4*G4ppxx + 72*G3p*G3x*G4*G4ppxx - 72*G3x*G4*G4pp*G4ppxx - 54*G2x*pow(G3x,2)*G4px + 108*G3p*pow(G3x,2)*G4px - 36*G2x*G3pxx*G4*G4px + 72*G3p*G3pxx*G4*G4px - 
   96*G2pxx*G3x*G4*G4px + 24*G3ppx*G3x*G4*G4px - 36*G3pp*G3xx*G4*G4px + 270*pow(G3x,2)*G4pp*G4px - 72*G3pxx*G4*G4pp*G4px + 72*G3xx*G4*G4ppp*G4px + 144*G3x*G4*G4pppx*G4px - 216*G3x*G4p*G4ppx*G4px + 
   288*G4*pow(G4ppx,2)*G4px + 72*G2x*G4*G4ppxx*G4px - 144*G3p*G4*G4ppxx*G4px + 144*G4*G4pp*G4ppxx*G4px + 276*G2x*G3x*pow(G4px,2) - 552*G3p*G3x*pow(G4px,2) + 120*G2pxx*G4*pow(G4px,2) - 
   48*G3ppx*G4*pow(G4px,2) - 468*G3x*G4pp*pow(G4px,2) - 144*G4*G4pppx*pow(G4px,2) + 504*G4p*G4ppx*pow(G4px,2) - 336*G2x*pow(G4px,3) + 672*G3p*pow(G4px,3) + 216*G4pp*pow(G4px,3) - 
   48*G2px*G3x*G4*G4pxx + 12*G3pp*G3x*G4*G4pxx - 108*G2x*G3x*G4p*G4pxx + 216*G3p*G3x*G4p*G4pxx - 216*G3x*G4p*G4pp*G4pxx + 72*G3x*G4*G4ppp*G4pxx - 72*G2x*G4*G4ppx*G4pxx + 144*G3p*G4*G4ppx*G4pxx - 
   144*G4*G4pp*G4ppx*G4pxx + 96*G2px*G4*G4px*G4pxx - 24*G3pp*G4*G4px*G4pxx + 144*G2x*G4p*G4px*G4pxx - 288*G3p*G4p*G4px*G4pxx + 288*G4p*G4pp*G4px*G4pxx - 144*G4*G4ppp*G4px*G4pxx - 
   36*G2px*pow(G3x,2)*G4x + 90*G3pp*pow(G3x,2)*G4x + 18*G2p*G3x*G3xx*G4x + 48*G2px*G3pxx*G4*G4x - 48*G3pp*G3pxx*G4*G4x - 48*G2ppx*G3xx*G4*G4x + 48*G3ppp*G3xx*G4*G4x + 36*G2x*G3pxx*G4p*G4x - 
   72*G3p*G3pxx*G4p*G4x + 48*G2pxx*G3x*G4p*G4x + 24*G3ppx*G3x*G4p*G4x + 36*G3pp*G3xx*G4p*G4x - 36*G2x*G3xx*G4pp*G4x + 72*G3p*G3xx*G4pp*G4x + 72*G3pxx*G4p*G4pp*G4x - 72*G3xx*pow(G4pp,2)*G4x - 
   36*pow(G3x,2)*G4ppp*G4x - 72*G3xx*G4p*G4ppp*G4x - 144*G3x*G4p*G4pppx*G4x - 168*G2x*G3x*G4ppx*G4x + 336*G3p*G3x*G4ppx*G4x - 96*G2pxx*G4*G4ppx*G4x + 96*G3ppx*G4*G4ppx*G4x + 72*G3x*G4pp*G4ppx*G4x - 
   288*G4p*pow(G4ppx,2)*G4x - 96*G2px*G4*G4ppxx*G4x + 96*G3pp*G4*G4ppxx*G4x - 72*G2x*G4p*G4ppxx*G4x + 144*G3p*G4p*G4ppxx*G4x - 144*G4p*G4pp*G4ppxx*G4x + 168*G2px*G3x*G4px*G4x - 372*G3pp*G3x*G4px*G4x - 
   36*G2p*G3xx*G4px*G4x - 144*G2pxx*G4p*G4px*G4x + 72*G3x*G4ppp*G4px*G4x + 288*G4p*G4pppx*G4px*G4x + 456*G2x*G4ppx*G4px*G4x - 912*G3p*G4ppx*G4px*G4x - 192*G2px*pow(G4px,2)*G4x + 
   384*G3pp*pow(G4px,2)*G4x + 24*pow(G2x,2)*G4pxx*G4x - 96*G2x*G3p*G4pxx*G4x + 96*pow(G3p,2)*G4pxx*G4x - 84*G2p*G3x*G4pxx*G4x + 96*G2ppx*G4*G4pxx*G4x - 96*G3ppp*G4*G4pxx*G4x - 
   96*G2px*G4p*G4pxx*G4x + 24*G3pp*G4p*G4pxx*G4x + 120*G2x*G4pp*G4pxx*G4x - 240*G3p*G4pp*G4pxx*G4x + 144*pow(G4pp,2)*G4pxx*G4x + 144*G4p*G4ppp*G4pxx*G4x + 168*G2p*G4px*G4pxx*G4x - 
   48*G2pxx*G2x*pow(G4x,2) + 96*G2pxx*G3p*pow(G4x,2) + 72*G2x*G3ppx*pow(G4x,2) - 144*G3p*G3ppx*pow(G4x,2) + 24*G2p*G3pxx*pow(G4x,2) - 24*G2pp*G3xx*pow(G4x,2) - 48*G2pxx*G4pp*pow(G4x,2) + 
   48*G3ppx*G4pp*pow(G4x,2) - 48*G2x*G4pppx*pow(G4x,2) + 96*G3p*G4pppx*pow(G4x,2) + 96*G2px*G4ppx*pow(G4x,2) - 144*G3pp*G4ppx*pow(G4x,2) - 48*G2p*G4ppxx*pow(G4x,2) - 
   48*G2ppx*G4px*pow(G4x,2) + 48*G3ppp*G4px*pow(G4x,2) + 48*G2pp*G4pxx*pow(G4x,2) - 54*G2p*pow(G3x,2)*G4xx + 24*G2pxx*G2x*G4*G4xx - 48*G2pxx*G3p*G4*G4xx - 24*G2x*G3ppx*G4*G4xx + 
   48*G3p*G3ppx*G4*G4xx + 48*G2ppx*G3x*G4*G4xx - 48*G3ppp*G3x*G4*G4xx + 48*G2px*G3x*G4p*G4xx - 264*G3pp*G3x*G4p*G4xx - 72*G2x*G3x*G4pp*G4xx + 144*G3p*G3x*G4pp*G4xx + 48*G2pxx*G4*G4pp*G4xx - 
   48*G3ppx*G4*G4pp*G4xx - 144*G3x*pow(G4pp,2)*G4xx + 432*G3x*G4p*G4ppp*G4xx + 96*G2px*G4*G4ppx*G4xx - 96*G3pp*G4*G4ppx*G4xx + 144*G2x*G4p*G4ppx*G4xx - 288*G3p*G4p*G4ppx*G4xx + 
   288*G4p*G4pp*G4ppx*G4xx - 24*pow(G2x,2)*G4px*G4xx + 96*G2x*G3p*G4px*G4xx - 96*pow(G3p,2)*G4px*G4xx + 264*G2p*G3x*G4px*G4xx - 96*G2ppx*G4*G4px*G4xx + 96*G3ppp*G4*G4px*G4xx + 
   432*G3pp*G4p*G4px*G4xx + 96*G2x*G4pp*G4px*G4xx - 192*G3p*G4pp*G4px*G4xx + 288*pow(G4pp,2)*G4px*G4xx - 864*G4p*G4ppp*G4px*G4xx - 312*G2p*pow(G4px,2)*G4xx - 24*G2px*G2x*G4x*G4xx + 
   48*G2px*G3p*G4x*G4xx + 24*G2x*G3pp*G4x*G4xx - 48*G3p*G3pp*G4x*G4xx + 48*G2pp*G3x*G4x*G4xx + 96*G2ppx*G4p*G4x*G4xx - 96*G3ppp*G4p*G4x*G4xx - 144*G2px*G4pp*G4x*G4xx + 192*G3pp*G4pp*G4x*G4xx + 
   48*G2x*G4ppp*G4x*G4xx - 96*G3p*G4ppp*G4x*G4xx + 96*G2p*G4ppx*G4x*G4xx - 96*G2pp*G4px*G4x*G4xx + 18*G2px*pow(G3x,2)*G5p - 72*G3pp*pow(G3x,2)*G5p - 9*G2p*G3x*G3xx*G5p - 24*G2px*G3pxx*G4*G5p + 
   24*G3pp*G3pxx*G4*G5p + 24*G2ppx*G3xx*G4*G5p - 24*G3ppp*G3xx*G4*G5p - 18*G2x*G3pxx*G4p*G5p + 36*G3p*G3pxx*G4p*G5p - 24*G2pxx*G3x*G4p*G5p - 12*G3ppx*G3x*G4p*G5p - 18*G3pp*G3xx*G4p*G5p + 
   18*G2x*G3xx*G4pp*G5p - 36*G3p*G3xx*G4pp*G5p - 36*G3pxx*G4p*G4pp*G5p + 36*G3xx*pow(G4pp,2)*G5p + 72*pow(G3x,2)*G4ppp*G5p + 36*G3xx*G4p*G4ppp*G5p + 72*G3x*G4p*G4pppx*G5p + 120*G2x*G3x*G4ppx*G5p - 
   240*G3p*G3x*G4ppx*G5p + 48*G2pxx*G4*G4ppx*G5p - 48*G3ppx*G4*G4ppx*G5p + 36*G3x*G4pp*G4ppx*G5p + 144*G4p*pow(G4ppx,2)*G5p + 48*G2px*G4*G4ppxx*G5p - 48*G3pp*G4*G4ppxx*G5p + 36*G2x*G4p*G4ppxx*G5p - 
   72*G3p*G4p*G4ppxx*G5p + 72*G4p*G4pp*G4ppxx*G5p - 60*G2px*G3x*G4px*G5p + 270*G3pp*G3x*G4px*G5p + 18*G2p*G3xx*G4px*G5p + 72*G2pxx*G4p*G4px*G5p - 252*G3x*G4ppp*G4px*G5p - 144*G4p*G4pppx*G4px*G5p - 
   300*G2x*G4ppx*G4px*G5p + 600*G3p*G4ppx*G4px*G5p - 144*G4pp*G4ppx*G4px*G5p + 48*G2px*pow(G4px,2)*G5p - 252*G3pp*pow(G4px,2)*G5p + 216*G4ppp*pow(G4px,2)*G5p - 12*pow(G2x,2)*G4pxx*G5p + 
   48*G2x*G3p*G4pxx*G5p - 48*pow(G3p,2)*G4pxx*G5p + 42*G2p*G3x*G4pxx*G5p - 48*G2ppx*G4*G4pxx*G5p + 48*G3ppp*G4*G4pxx*G5p + 48*G2px*G4p*G4pxx*G5p - 12*G3pp*G4p*G4pxx*G5p - 60*G2x*G4pp*G4pxx*G5p + 
   120*G3p*G4pp*G4pxx*G5p - 72*pow(G4pp,2)*G4pxx*G5p - 72*G4p*G4ppp*G4pxx*G5p - 84*G2p*G4px*G4pxx*G5p + 60*G2pxx*G2x*G4x*G5p - 120*G2pxx*G3p*G4x*G5p - 84*G2x*G3ppx*G4x*G5p + 168*G3p*G3ppx*G4x*G5p - 
   24*G2p*G3pxx*G4x*G5p + 24*G2ppx*G3x*G4x*G5p - 24*G3ppp*G3x*G4x*G5p + 24*G2pp*G3xx*G4x*G5p + 72*G2pxx*G4pp*G4x*G5p - 72*G3ppx*G4pp*G4x*G5p + 48*G2x*G4pppx*G4x*G5p - 96*G3p*G4pppx*G4x*G5p - 
   48*G2px*G4ppx*G4x*G5p + 96*G3pp*G4ppx*G4x*G5p + 48*G2p*G4ppxx*G4x*G5p - 48*G2pp*G4pxx*G4x*G5p + 12*G2px*G2x*G4xx*G5p - 24*G2px*G3p*G4xx*G5p - 12*G2x*G3pp*G4xx*G5p + 24*G3p*G3pp*G4xx*G5p - 
   24*G2pp*G3x*G4xx*G5p - 48*G2ppx*G4p*G4xx*G5p + 48*G3ppp*G4p*G4xx*G5p + 72*G2px*G4pp*G4xx*G5p - 96*G3pp*G4pp*G4xx*G5p - 24*G2x*G4ppp*G4xx*G5p + 48*G3p*G4ppp*G4xx*G5p - 48*G2p*G4ppx*G4xx*G5p + 
   48*G2pp*G4px*G4xx*G5p - 18*G2pxx*G2x*pow(G5p,2) + 36*G2pxx*G3p*pow(G5p,2) + 24*G2x*G3ppx*pow(G5p,2) - 48*G3p*G3ppx*pow(G5p,2) + 6*G2p*G3pxx*pow(G5p,2) - 12*G2ppx*G3x*pow(G5p,2) + 
   12*G3ppp*G3x*pow(G5p,2) - 6*G2pp*G3xx*pow(G5p,2) - 24*G2pxx*G4pp*pow(G5p,2) + 24*G3ppx*G4pp*pow(G5p,2) - 12*G2x*G4pppx*pow(G5p,2) + 24*G3p*G4pppx*pow(G5p,2) - 12*G3pp*G4ppx*pow(G5p,2) - 
   12*G2p*G4ppxx*pow(G5p,2) + 12*G2ppx*G4px*pow(G5p,2) - 12*G3ppp*G4px*pow(G5p,2) + 12*G2pp*G4pxx*pow(G5p,2) + 18*G2x*pow(G3x,2)*G5pp - 36*G3p*pow(G3x,2)*G5pp + 12*G2pxx*G3x*G4*G5pp - 
   12*G3ppx*G3x*G4*G5pp - 18*G2x*G3xx*G4p*G5pp + 36*G3p*G3xx*G4p*G5pp - 18*pow(G3x,2)*G4pp*G5pp - 36*G3xx*G4p*G4pp*G5pp + 144*G3x*G4p*G4ppx*G5pp - 114*G2x*G3x*G4px*G5pp + 228*G3p*G3x*G4px*G5pp - 
   24*G2pxx*G4*G4px*G5pp + 24*G3ppx*G4*G4px*G5pp + 36*G3x*G4pp*G4px*G5pp - 288*G4p*G4ppx*G4px*G5pp + 156*G2x*pow(G4px,2)*G5pp - 312*G3p*pow(G4px,2)*G5pp + 36*G2x*G4p*G4pxx*G5pp - 
   72*G3p*G4p*G4pxx*G5pp + 72*G4p*G4pp*G4pxx*G5pp - 36*G2px*G3x*G4x*G5pp + 48*G3pp*G3x*G4x*G5pp + 24*G2pxx*G4p*G4x*G5pp - 24*G3ppx*G4p*G4x*G5pp - 24*G2x*G4ppx*G4x*G5pp + 48*G3p*G4ppx*G4x*G5pp + 
   72*G2px*G4px*G4x*G5pp - 96*G3pp*G4px*G4x*G5pp + 12*pow(G2x,2)*G4xx*G5pp - 48*G2x*G3p*G4xx*G5pp + 48*pow(G3p,2)*G4xx*G5pp - 24*G2p*G3x*G4xx*G5pp - 48*G2px*G4p*G4xx*G5pp + 48*G3pp*G4p*G4xx*G5pp + 
   24*G2x*G4pp*G4xx*G5pp - 48*G3p*G4pp*G4xx*G5pp + 48*G2p*G4px*G4xx*G5pp + 6*G2px*G3x*G5p*G5pp - 12*G3pp*G3x*G5p*G5pp - 12*G2pxx*G4p*G5p*G5pp + 12*G3ppx*G4p*G5p*G5pp + 12*G2x*G4ppx*G5p*G5pp - 
   24*G3p*G4ppx*G5p*G5pp - 12*G2px*G4px*G5p*G5pp + 24*G3pp*G4px*G5p*G5pp + 6*G2x*G3x*pow(G5pp,2) - 12*G3p*G3x*pow(G5pp,2) - 12*G2x*G4px*pow(G5pp,2) + 24*G3p*G4px*pow(G5pp,2) + 
   54*pow(G3x,2)*G4p*G5ppp - 216*G3x*G4p*G4px*G5ppp + 216*G4p*pow(G4px,2)*G5ppp + 12*G2x*G3x*G4x*G5ppp - 24*G3p*G3x*G4x*G5ppp - 24*G2x*G4px*G4x*G5ppp + 48*G3p*G4px*G4x*G5ppp - 6*G2x*G3x*G5p*G5ppp + 
   12*G3p*G3x*G5p*G5ppp + 12*G2x*G4px*G5p*G5ppp - 24*G3p*G4px*G5p*G5ppp + 24*G2px*G3x*G4*G5ppx - 24*G3pp*G3x*G4*G5ppx + 36*G2x*G3x*G4p*G5ppx - 72*G3p*G3x*G4p*G5ppx + 72*G3x*G4p*G4pp*G5ppx - 
   48*G2px*G4*G4px*G5ppx + 48*G3pp*G4*G4px*G5ppx - 72*G2x*G4p*G4px*G5ppx + 144*G3p*G4p*G4px*G5ppx - 144*G4p*G4pp*G4px*G5ppx - 12*pow(G2x,2)*G4x*G5ppx + 48*G2x*G3p*G4x*G5ppx - 
   48*pow(G3p,2)*G4x*G5ppx + 24*G2p*G3x*G4x*G5ppx + 48*G2px*G4p*G4x*G5ppx - 48*G3pp*G4p*G4x*G5ppx - 24*G2x*G4pp*G4x*G5ppx + 48*G3p*G4pp*G4x*G5ppx - 48*G2p*G4px*G4x*G5ppx + 6*pow(G2x,2)*G5p*G5ppx - 
   24*G2x*G3p*G5p*G5ppx + 24*pow(G3p,2)*G5p*G5ppx - 12*G2p*G3x*G5p*G5ppx - 24*G2px*G4p*G5p*G5ppx + 24*G3pp*G4p*G5p*G5ppx + 12*G2x*G4pp*G5p*G5ppx - 24*G3p*G4pp*G5p*G5ppx + 24*G2p*G4px*G5p*G5ppx - 
   3*pow(G2x,2)*G3x*G5px + 12*G2x*G3p*G3x*G5px - 12*pow(G3p,2)*G3x*G5px + 30*G2p*pow(G3x,2)*G5px - 12*G2pxx*G2x*G4*G5px + 24*G2pxx*G3p*G4*G5px + 12*G2x*G3ppx*G4*G5px - 24*G3p*G3ppx*G4*G5px - 
   24*G2ppx*G3x*G4*G5px + 24*G3ppp*G3x*G4*G5px - 12*G2px*G3x*G4p*G5px + 120*G3pp*G3x*G4p*G5px + 30*G2x*G3x*G4pp*G5px - 60*G3p*G3x*G4pp*G5px - 24*G2pxx*G4*G4pp*G5px + 24*G3ppx*G4*G4pp*G5px + 
   72*G3x*pow(G4pp,2)*G5px - 216*G3x*G4p*G4ppp*G5px - 48*G2px*G4*G4ppx*G5px + 48*G3pp*G4*G4ppx*G5px - 72*G2x*G4p*G4ppx*G5px + 144*G3p*G4p*G4ppx*G5px - 144*G4p*G4pp*G4ppx*G5px + 
   18*pow(G2x,2)*G4px*G5px - 72*G2x*G3p*G4px*G5px + 72*pow(G3p,2)*G4px*G5px - 144*G2p*G3x*G4px*G5px + 48*G2ppx*G4*G4px*G5px - 48*G3ppp*G4*G4px*G5px - 24*G2px*G4p*G4px*G5px - 192*G3pp*G4p*G4px*G5px - 
   36*G2x*G4pp*G4px*G5px + 72*G3p*G4pp*G4px*G5px - 144*pow(G4pp,2)*G4px*G5px + 432*G4p*G4ppp*G4px*G5px + 168*G2p*pow(G4px,2)*G5px + 12*G2px*G2x*G4x*G5px - 24*G2px*G3p*G4x*G5px - 
   12*G2x*G3pp*G4x*G5px + 24*G3p*G3pp*G4x*G5px - 24*G2pp*G3x*G4x*G5px - 48*G2ppx*G4p*G4x*G5px + 48*G3ppp*G4p*G4x*G5px + 72*G2px*G4pp*G4x*G5px - 96*G3pp*G4pp*G4x*G5px - 24*G2x*G4ppp*G4x*G5px + 
   48*G3p*G4ppp*G4x*G5px - 48*G2p*G4ppx*G4x*G5px + 48*G2pp*G4px*G4x*G5px - 6*G2px*G2x*G5p*G5px + 12*G2px*G3p*G5p*G5px + 6*G2x*G3pp*G5p*G5px - 12*G3p*G3pp*G5p*G5px + 12*G2pp*G3x*G5p*G5px + 
   24*G2ppx*G4p*G5p*G5px - 24*G3ppp*G4p*G5p*G5px - 36*G2px*G4pp*G5p*G5px + 48*G3pp*G4pp*G5p*G5px + 12*G2x*G4ppp*G5p*G5px - 24*G3p*G4ppp*G5p*G5px + 24*G2p*G4ppx*G5p*G5px - 24*G2pp*G4px*G5p*G5px - 
   6*pow(G2x,2)*G5pp*G5px + 24*G2x*G3p*G5pp*G5px - 24*pow(G3p,2)*G5pp*G5px + 12*G2p*G3x*G5pp*G5px + 24*G2px*G4p*G5pp*G5px - 24*G3pp*G4p*G5pp*G5px - 12*G2x*G4pp*G5pp*G5px + 24*G3p*G4pp*G5pp*G5px - 
   24*G2p*G4px*G5pp*G5px - 3*G2px*G2x*G3x*G5x + 6*G2px*G3p*G3x*G5x + 3*G2x*G3pp*G3x*G5x - 6*G3p*G3pp*G3x*G5x + 3*G2pp*pow(G3x,2)*G5x + 8*G2px*G2pxx*G4*G5x - 8*G2pxx*G3pp*G4*G5x - 8*G2px*G3ppx*G4*G5x + 
   8*G3pp*G3ppx*G4*G5x + 6*G2pxx*G2x*G4p*G5x - 12*G2pxx*G3p*G4p*G5x - 6*G2x*G3ppx*G4p*G5x + 12*G3p*G3ppx*G4p*G5x + 12*G2ppx*G3x*G4p*G5x - 12*G3ppp*G3x*G4p*G5x - 18*G2px*G3x*G4pp*G5x + 
   24*G3pp*G3x*G4pp*G5x + 12*G2pxx*G4p*G4pp*G5x - 12*G3ppx*G4p*G4pp*G5x + 6*G2x*G3x*G4ppp*G5x - 12*G3p*G3x*G4ppp*G5x - 6*pow(G2x,2)*G4ppx*G5x + 24*G2x*G3p*G4ppx*G5x - 24*pow(G3p,2)*G4ppx*G5x + 
   12*G2p*G3x*G4ppx*G5x + 24*G2px*G4p*G4ppx*G5x - 24*G3pp*G4p*G4ppx*G5x - 12*G2x*G4pp*G4ppx*G5x + 24*G3p*G4pp*G4ppx*G5x + 6*G2px*G2x*G4px*G5x - 12*G2px*G3p*G4px*G5x - 6*G2x*G3pp*G4px*G5x + 
   12*G3p*G3pp*G4px*G5x - 12*G2pp*G3x*G4px*G5x - 24*G2ppx*G4p*G4px*G5x + 24*G3ppp*G4p*G4px*G5x + 36*G2px*G4pp*G4px*G5x - 48*G3pp*G4pp*G4px*G5x - 12*G2x*G4ppp*G4px*G5x + 24*G3p*G4ppp*G4px*G5x - 
   24*G2p*G4ppx*G4px*G5x + 12*G2pp*pow(G4px,2)*G5x - 8*pow(G2px,2)*G4x*G5x + 8*G2p*G2pxx*G4x*G5x + 8*G2ppx*G2x*G4x*G5x - 16*G2ppx*G3p*G4x*G5x + 24*G2px*G3pp*G4x*G5x - 16*pow(G3pp,2)*G4x*G5x - 
   8*G2x*G3ppp*G4x*G5x + 16*G3p*G3ppp*G4x*G5x - 8*G2p*G3ppx*G4x*G5x + 4*pow(G2px,2)*G5p*G5x - 4*G2p*G2pxx*G5p*G5x - 4*G2ppx*G2x*G5p*G5x + 8*G2ppx*G3p*G5p*G5x - 12*G2px*G3pp*G5p*G5x + 
   8*pow(G3pp,2)*G5p*G5x + 4*G2x*G3ppp*G5p*G5x - 8*G3p*G3ppp*G5p*G5x + 4*G2p*G3ppx*G5p*G5x + 144*G3xx*pow(G4,2)*G4ppxx*pow(H,2) - 108*G3x*G3xx*G4*G4px*pow(H,2) + 
   504*G3xx*G4*pow(G4px,2)*pow(H,2) - 324*pow(G3x,2)*G4*G4pxx*pow(H,2) - 144*G3pxx*pow(G4,2)*G4pxx*pow(H,2) + 144*G3xx*G4*G4p*G4pxx*pow(H,2) + 1656*G3x*G4*G4px*G4pxx*pow(H,2) - 
   3168*G4*pow(G4px,2)*G4pxx*pow(H,2) + 288*G4*G4p*pow(G4pxx,2)*pow(H,2) - 288*G3x*G4*G4p*G4pxxx*pow(H,2) + 288*pow(G4,2)*G4ppx*G4pxxx*pow(H,2) + 864*G4*G4p*G4px*G4pxxx*pow(H,2) + 
   756*G3pxx*G3x*G4*G4x*pow(H,2) + 432*G3x*G3xx*G4p*G4x*pow(H,2) + 1368*G3xx*G4*G4ppx*G4x*pow(H,2) - 1800*G3x*G4*G4ppxx*G4x*pow(H,2) - 1800*G3pxx*G4*G4px*G4x*pow(H,2) - 
   720*G3xx*G4p*G4px*G4x*pow(H,2) + 4752*G4*G4ppxx*G4px*G4x*pow(H,2) + 144*G3x*pow(G4px,2)*G4x*pow(H,2) - 288*pow(G4px,3)*G4x*pow(H,2) - 2304*G3x*G4p*G4pxx*G4x*pow(H,2) - 
   5040*G4*G4ppx*G4pxx*G4x*pow(H,2) + 4608*G4p*G4px*G4pxx*G4x*pow(H,2) + 432*G2x*G4*G4pxxx*G4x*pow(H,2) - 864*G3p*G4*G4pxxx*G4x*pow(H,2) - 144*pow(G4p,2)*G4pxxx*G4x*pow(H,2) + 
   288*G4*G4pp*G4pxxx*G4x*pow(H,2) + 648*G3pxx*G4p*pow(G4x,2)*pow(H,2) - 216*G3xx*G4pp*pow(G4x,2)*pow(H,2) + 576*G3x*G4ppx*pow(G4x,2)*pow(H,2) - 720*G4p*G4ppxx*pow(G4x,2)*pow(H,2) - 
   432*G4ppx*G4px*pow(G4x,2)*pow(H,2) + 144*G2x*G4pxx*pow(G4x,2)*pow(H,2) - 288*G3p*G4pxx*pow(G4x,2)*pow(H,2) - 432*G4pp*G4pxx*pow(G4x,2)*pow(H,2) - 288*G2pxx*pow(G4x,3)*pow(H,2) + 
   144*G3ppx*pow(G4x,3)*pow(H,2) + 288*G4pppx*pow(G4x,3)*pow(H,2) + 1296*pow(G3x,2)*G4p*G4xx*pow(H,2) - 360*G3pxx*G4*G4p*G4xx*pow(H,2) - 360*G3xx*pow(G4p,2)*G4xx*pow(H,2) - 
   72*G3xx*G4*G4pp*G4xx*pow(H,2) + 936*G3x*G4*G4ppx*G4xx*pow(H,2) + 144*G4*G4p*G4ppxx*G4xx*pow(H,2) - 6120*G3x*G4p*G4px*G4xx*pow(H,2) - 1152*G4*G4ppx*G4px*G4xx*pow(H,2) + 
   6912*G4p*pow(G4px,2)*G4xx*pow(H,2) - 720*G2x*G4*G4pxx*G4xx*pow(H,2) + 1440*G3p*G4*G4pxx*G4xx*pow(H,2) + 1152*pow(G4p,2)*G4pxx*G4xx*pow(H,2) + 720*G4*G4pp*G4pxx*G4xx*pow(H,2) + 
   576*G2pxx*G4*G4x*G4xx*pow(H,2) - 1152*G3ppx*G4*G4x*G4xx*pow(H,2) - 4320*G3x*G4pp*G4x*G4xx*pow(H,2) + 1152*G4*G4pppx*G4x*G4xx*pow(H,2) - 2448*G4p*G4ppx*G4x*G4xx*pow(H,2) + 
   11376*G4pp*G4px*G4x*G4xx*pow(H,2) - 432*G2px*pow(G4x,2)*G4xx*pow(H,2) + 1008*G3pp*pow(G4x,2)*G4xx*pow(H,2) - 288*G4ppp*pow(G4x,2)*G4xx*pow(H,2) + 432*G2px*G4*pow(G4xx,2)*pow(H,2) - 
   288*G3pp*G4*pow(G4xx,2)*pow(H,2) + 1008*G2x*G4p*pow(G4xx,2)*pow(H,2) - 2016*G3p*G4p*pow(G4xx,2)*pow(H,2) + 4032*G4p*G4pp*pow(G4xx,2)*pow(H,2) - 1152*G4*G4ppp*pow(G4xx,2)*pow(H,2) - 
   288*G2p*G4x*pow(G4xx,2)*pow(H,2) + 144*G3ppx*pow(G4,2)*G4xxx*pow(H,2) - 288*G3x*pow(G4p,2)*G4xxx*pow(H,2) + 72*G3x*G4*G4pp*G4xxx*pow(H,2) - 288*pow(G4,2)*G4pppx*G4xxx*pow(H,2) - 
   432*G4*G4p*G4ppx*G4xxx*pow(H,2) - 144*G2x*G4*G4px*G4xxx*pow(H,2) + 288*G3p*G4*G4px*G4xxx*pow(H,2) + 720*pow(G4p,2)*G4px*G4xxx*pow(H,2) - 720*G4*G4pp*G4px*G4xxx*pow(H,2) - 
   288*G2px*G4*G4x*G4xxx*pow(H,2) + 720*G3pp*G4*G4x*G4xxx*pow(H,2) + 144*G2x*G4p*G4x*G4xxx*pow(H,2) - 288*G3p*G4p*G4x*G4xxx*pow(H,2) + 144*G4p*G4pp*G4x*G4xxx*pow(H,2) - 
   288*G4*G4ppp*G4x*G4xxx*pow(H,2) + 144*G2p*pow(G4x,2)*G4xxx*pow(H,2) - 144*G2p*G4*G4xx*G4xxx*pow(H,2) - 540*G3pxx*G3x*G4*G5p*pow(H,2) - 378*G3x*G3xx*G4p*G5p*pow(H,2) - 
   936*G3xx*G4*G4ppx*G5p*pow(H,2) + 1080*G3x*G4*G4ppxx*G5p*pow(H,2) + 216*pow(G3x,2)*G4px*G5p*pow(H,2) + 1368*G3pxx*G4*G4px*G5p*pow(H,2) + 648*G3xx*G4p*G4px*G5p*pow(H,2) - 
   3024*G4*G4ppxx*G4px*G5p*pow(H,2) - 1332*G3x*pow(G4px,2)*G5p*pow(H,2) + 2088*pow(G4px,3)*G5p*pow(H,2) + 2124*G3x*G4p*G4pxx*G5p*pow(H,2) + 2160*G4*G4ppx*G4pxx*G5p*pow(H,2) - 
   5040*G4p*G4px*G4pxx*G5p*pow(H,2) - 288*G2x*G4*G4pxxx*G5p*pow(H,2) + 576*G3p*G4*G4pxxx*G5p*pow(H,2) + 72*pow(G4p,2)*G4pxxx*G5p*pow(H,2) - 288*G4*G4pp*G4pxxx*G5p*pow(H,2) - 
   828*G3pxx*G4p*G4x*G5p*pow(H,2) + 180*G3xx*G4pp*G4x*G5p*pow(H,2) - 648*G3x*G4ppx*G4x*G5p*pow(H,2) + 792*G4p*G4ppxx*G4x*G5p*pow(H,2) + 216*G4ppx*G4px*G4x*G5p*pow(H,2) - 
   576*G2x*G4pxx*G4x*G5p*pow(H,2) + 1152*G3p*G4pxx*G4x*G5p*pow(H,2) + 648*G4pp*G4pxx*G4x*G5p*pow(H,2) + 648*G2pxx*pow(G4x,2)*G5p*pow(H,2) - 216*G3ppx*pow(G4x,2)*G5p*pow(H,2) - 
   864*G4pppx*pow(G4x,2)*G5p*pow(H,2) - 432*G2pxx*G4*G4xx*G5p*pow(H,2) + 432*G3ppx*G4*G4xx*G5p*pow(H,2) + 3888*G3x*G4pp*G4xx*G5p*pow(H,2) + 2376*G4p*G4ppx*G4xx*G5p*pow(H,2) + 
   504*G2x*G4px*G4xx*G5p*pow(H,2) - 1008*G3p*G4px*G4xx*G5p*pow(H,2) - 9576*G4pp*G4px*G4xx*G5p*pow(H,2) + 864*G2px*G4x*G4xx*G5p*pow(H,2) - 2376*G3pp*G4x*G4xx*G5p*pow(H,2) + 
   1296*G4ppp*G4x*G4xx*G5p*pow(H,2) + 576*G2p*pow(G4xx,2)*G5p*pow(H,2) + 144*G2px*G4*G4xxx*G5p*pow(H,2) - 432*G3pp*G4*G4xxx*G5p*pow(H,2) - 144*G2x*G4p*G4xxx*G5p*pow(H,2) + 
   288*G3p*G4p*G4xxx*G5p*pow(H,2) - 216*G4p*G4pp*G4xxx*G5p*pow(H,2) + 288*G4*G4ppp*G4xxx*G5p*pow(H,2) - 216*G2p*G4x*G4xxx*G5p*pow(H,2) + 252*G3pxx*G4p*pow(G5p,2)*pow(H,2) - 
   36*G3xx*G4pp*pow(G5p,2)*pow(H,2) + 180*G3x*G4ppx*pow(G5p,2)*pow(H,2) - 216*G4p*G4ppxx*pow(G5p,2)*pow(H,2) + 288*G4ppx*G4px*pow(G5p,2)*pow(H,2) + 324*G2x*G4pxx*pow(G5p,2)*pow(H,2) - 
   648*G3p*G4pxx*pow(G5p,2)*pow(H,2) - 72*G4pp*G4pxx*pow(G5p,2)*pow(H,2) - 468*G2pxx*G4x*pow(G5p,2)*pow(H,2) + 72*G3ppx*G4x*pow(G5p,2)*pow(H,2) + 792*G4pppx*G4x*pow(G5p,2)*pow(H,2) - 
   324*G2px*G4xx*pow(G5p,2)*pow(H,2) + 1152*G3pp*G4xx*pow(G5p,2)*pow(H,2) - 1008*G4ppp*G4xx*pow(G5p,2)*pow(H,2) + 72*G2p*G4xxx*pow(G5p,2)*pow(H,2) + 108*G2pxx*pow(G5p,3)*pow(H,2) - 
   216*G4pppx*pow(G5p,3)*pow(H,2) + 162*G3x*G3xx*G4*G5pp*pow(H,2) - 756*G3xx*G4*G4px*G5pp*pow(H,2) + 108*G3x*G4*G4pxx*G5pp*pow(H,2) + 936*G4*G4px*G4pxx*G5pp*pow(H,2) - 
   144*G4*G4p*G4pxxx*G5pp*pow(H,2) - 594*pow(G3x,2)*G4x*G5pp*pow(H,2) + 3168*G3x*G4px*G4x*G5pp*pow(H,2) - 4248*pow(G4px,2)*G4x*G5pp*pow(H,2) + 576*G4p*G4pxx*G4x*G5pp*pow(H,2) + 
   864*G4ppx*pow(G4x,2)*G5pp*pow(H,2) + 792*G3x*G4p*G4xx*G5pp*pow(H,2) - 1728*G4*G4ppx*G4xx*G5pp*pow(H,2) - 2376*G4p*G4px*G4xx*G5pp*pow(H,2) - 864*G2x*G4x*G4xx*G5pp*pow(H,2) + 
   1728*G3p*G4x*G4xx*G5pp*pow(H,2) - 2160*G4pp*G4x*G4xx*G5pp*pow(H,2) + 144*G2x*G4*G4xxx*G5pp*pow(H,2) - 288*G3p*G4*G4xxx*G5pp*pow(H,2) - 72*pow(G4p,2)*G4xxx*G5pp*pow(H,2) + 
   432*G4*G4pp*G4xxx*G5pp*pow(H,2) + 324*pow(G3x,2)*G5p*G5pp*pow(H,2) - 72*G3pxx*G4*G5p*G5pp*pow(H,2) + 18*G3xx*G4p*G5p*G5pp*pow(H,2) + 144*G4*G4ppxx*G5p*G5pp*pow(H,2) - 
   1674*G3x*G4px*G5p*G5pp*pow(H,2) + 1908*pow(G4px,2)*G5p*G5pp*pow(H,2) + 108*G4p*G4pxx*G5p*G5pp*pow(H,2) - 1368*G4ppx*G4x*G5p*G5pp*pow(H,2) + 324*G2x*G4xx*G5p*G5pp*pow(H,2) - 
   648*G3p*G4xx*G5p*G5pp*pow(H,2) + 1584*G4pp*G4xx*G5p*G5pp*pow(H,2) + 324*G4ppx*pow(G5p,2)*G5pp*pow(H,2) + 144*G3xx*G4*pow(G5pp,2)*pow(H,2) - 288*G4*G4pxx*pow(G5pp,2)*pow(H,2) - 
   252*G3x*G4x*pow(G5pp,2)*pow(H,2) + 648*G4px*G4x*pow(G5pp,2)*pow(H,2) + 432*G4p*G4xx*pow(G5pp,2)*pow(H,2) + 216*G3x*G5p*pow(G5pp,2)*pow(H,2) - 432*G4px*G5p*pow(G5pp,2)*pow(H,2) - 
   72*G3x*pow(G4x,2)*G5ppp*pow(H,2) - 720*G3x*G4*G4xx*G5ppp*pow(H,2) + 1728*G4*G4px*G4xx*G5ppp*pow(H,2) + 864*G4p*G4x*G4xx*G5ppp*pow(H,2) + 144*G4*G4p*G4xxx*G5ppp*pow(H,2) + 
   72*G3xx*G4*G5p*G5ppp*pow(H,2) - 144*G4*G4pxx*G5p*G5ppp*pow(H,2) + 252*G3x*G4x*G5p*G5ppp*pow(H,2) - 216*G4px*G4x*G5p*G5ppp*pow(H,2) - 1296*G4p*G4xx*G5p*G5ppp*pow(H,2) - 
   216*G3x*pow(G5p,2)*G5ppp*pow(H,2) + 324*G4px*pow(G5p,2)*G5ppp*pow(H,2) - 72*G3xx*pow(G4,2)*G5pppx*pow(H,2) + 144*pow(G4,2)*G4pxx*G5pppx*pow(H,2) + 144*G3x*G4*G4x*G5pppx*pow(H,2) - 
   576*G4*G4px*G4x*G5pppx*pow(H,2) - 288*G4p*pow(G4x,2)*G5pppx*pow(H,2) + 288*G4*G4p*G4xx*G5pppx*pow(H,2) + 144*G4*G4px*G5p*G5pppx*pow(H,2) + 432*G4p*G4x*G5p*G5pppx*pow(H,2) - 
   144*G4p*pow(G5p,2)*G5pppx*pow(H,2) + 156*pow(G3x,2)*G4*G5ppx*pow(H,2) + 72*G3pxx*pow(G4,2)*G5ppx*pow(H,2) - 108*G3xx*G4*G4p*G5ppx*pow(H,2) - 144*pow(G4,2)*G4ppxx*G5ppx*pow(H,2) - 
   660*G3x*G4*G4px*G5ppx*pow(H,2) + 984*G4*pow(G4px,2)*G5ppx*pow(H,2) - 360*G4*G4p*G4pxx*G5ppx*pow(H,2) + 480*G3x*G4p*G4x*G5ppx*pow(H,2) + 1440*G4*G4ppx*G4x*G5ppx*pow(H,2) - 
   1032*G4p*G4px*G4x*G5ppx*pow(H,2) + 64*G2x*pow(G4x,2)*G5ppx*pow(H,2) - 128*G3p*pow(G4x,2)*G5ppx*pow(H,2) + 576*G4pp*pow(G4x,2)*G5ppx*pow(H,2) + 312*G2x*G4*G4xx*G5ppx*pow(H,2) - 
   624*G3p*G4*G4xx*G5ppx*pow(H,2) - 432*G4*G4pp*G4xx*G5ppx*pow(H,2) - 456*G3x*G4p*G5p*G5ppx*pow(H,2) - 288*G4*G4ppx*G5p*G5ppx*pow(H,2) + 1380*G4p*G4px*G5p*G5ppx*pow(H,2) + 
   128*G2x*G4x*G5p*G5ppx*pow(H,2) - 256*G3p*G4x*G5p*G5ppx*pow(H,2) - 720*G4pp*G4x*G5p*G5ppx*pow(H,2) - 116*G2x*pow(G5p,2)*G5ppx*pow(H,2) + 232*G3p*pow(G5p,2)*G5ppx*pow(H,2) + 
   144*G4pp*pow(G5p,2)*G5ppx*pow(H,2) - 252*G3x*G4*G5pp*G5ppx*pow(H,2) + 360*G4*G4px*G5pp*G5ppx*pow(H,2) - 360*G4p*G4x*G5pp*G5ppx*pow(H,2) - 36*G4p*G5p*G5pp*G5ppx*pow(H,2) + 
   144*G4*G4p*pow(G5ppx,2)*pow(H,2) + 120*G3x*G4*G4p*G5ppxx*pow(H,2) - 144*pow(G4,2)*G4ppx*G5ppxx*pow(H,2) - 384*G4*G4p*G4px*G5ppxx*pow(H,2) - 232*G2x*G4*G4x*G5ppxx*pow(H,2) + 
   464*G3p*G4*G4x*G5ppxx*pow(H,2) + 48*pow(G4p,2)*G4x*G5ppxx*pow(H,2) - 144*G4*G4pp*G4x*G5ppxx*pow(H,2) + 152*G2x*G4*G5p*G5ppxx*pow(H,2) - 304*G3p*G4*G5p*G5ppxx*pow(H,2) - 
   24*pow(G4p,2)*G5p*G5ppxx*pow(H,2) + 144*G4*G4pp*G5p*G5ppxx*pow(H,2) + 72*G4*G4p*G5pp*G5ppxx*pow(H,2) - 54*G2x*G3xx*G4*G5px*pow(H,2) + 108*G3p*G3xx*G4*G5px*pow(H,2) - 
   681*pow(G3x,2)*G4p*G5px*pow(H,2) + 216*G3pxx*G4*G4p*G5px*pow(H,2) + 180*G3xx*pow(G4p,2)*G5px*pow(H,2) - 108*G3xx*G4*G4pp*G5px*pow(H,2) - 348*G3x*G4*G4ppx*G5px*pow(H,2) - 
   144*G4*G4p*G4ppxx*G5px*pow(H,2) + 3132*G3x*G4p*G4px*G5px*pow(H,2) + 552*G4*G4ppx*G4px*G5px*pow(H,2) - 3108*G4p*pow(G4px,2)*G5px*pow(H,2) + 564*G2x*G4*G4pxx*G5px*pow(H,2) - 
   1128*G3p*G4*G4pxx*G5px*pow(H,2) - 720*pow(G4p,2)*G4pxx*G5px*pow(H,2) + 72*G4*G4pp*G4pxx*G5px*pow(H,2) + 114*G2x*G3x*G4x*G5px*pow(H,2) - 228*G3p*G3x*G4x*G5px*pow(H,2) - 
   280*G2pxx*G4*G4x*G5px*pow(H,2) + 712*G3ppx*G4*G4x*G5px*pow(H,2) + 2472*G3x*G4pp*G4x*G5px*pow(H,2) - 864*G4*G4pppx*G4x*G5px*pow(H,2) + 1320*G4p*G4ppx*G4x*G5px*pow(H,2) - 
   288*G2x*G4px*G4x*G5px*pow(H,2) + 576*G3p*G4px*G4x*G5px*pow(H,2) - 6600*G4pp*G4px*G4x*G5px*pow(H,2) + 200*G2px*pow(G4x,2)*G5px*pow(H,2) - 352*G3pp*pow(G4x,2)*G5px*pow(H,2) - 
   600*G2px*G4*G4xx*G5px*pow(H,2) + 408*G3pp*G4*G4xx*G5px*pow(H,2) - 972*G2x*G4p*G4xx*G5px*pow(H,2) + 1944*G3p*G4p*G4xx*G5px*pow(H,2) - 4464*G4p*G4pp*G4xx*G5px*pow(H,2) + 
   1296*G4*G4ppp*G4xx*G5px*pow(H,2) + 168*G2p*G4x*G4xx*G5px*pow(H,2) + 72*G2p*G4*G4xxx*G5px*pow(H,2) - 3*G2x*G3x*G5p*G5px*pow(H,2) + 6*G3p*G3x*G5p*G5px*pow(H,2) + 
   224*G2pxx*G4*G5p*G5px*pow(H,2) - 296*G3ppx*G4*G5p*G5px*pow(H,2) - 2154*G3x*G4pp*G5p*G5px*pow(H,2) + 144*G4*G4pppx*G5p*G5px*pow(H,2) - 1020*G4p*G4ppx*G5p*G5px*pow(H,2) - 
   324*G2x*G4px*G5p*G5px*pow(H,2) + 648*G3p*G4px*G5p*G5px*pow(H,2) + 5208*G4pp*G4px*G5p*G5px*pow(H,2) - 452*G2px*G4x*G5p*G5px*pow(H,2) + 1048*G3pp*G4x*G5p*G5px*pow(H,2) - 
   432*G4ppp*G4x*G5p*G5px*pow(H,2) - 564*G2p*G4xx*G5p*G5px*pow(H,2) + 164*G2px*pow(G5p,2)*G5px*pow(H,2) - 532*G3pp*pow(G5p,2)*G5px*pow(H,2) + 432*G4ppp*pow(G5p,2)*G5px*pow(H,2) - 
   492*G3x*G4p*G5pp*G5px*pow(H,2) + 792*G4*G4ppx*G5pp*G5px*pow(H,2) + 1128*G4p*G4px*G5pp*G5px*pow(H,2) + 492*G2x*G4x*G5pp*G5px*pow(H,2) - 984*G3p*G4x*G5pp*G5px*pow(H,2) + 
   1224*G4pp*G4x*G5pp*G5px*pow(H,2) - 156*G2x*G5p*G5pp*G5px*pow(H,2) + 312*G3p*G5p*G5pp*G5px*pow(H,2) - 828*G4pp*G5p*G5pp*G5px*pow(H,2) - 180*G4p*pow(G5pp,2)*G5px*pow(H,2) + 
   396*G3x*G4*G5ppp*G5px*pow(H,2) - 936*G4*G4px*G5ppp*G5px*pow(H,2) - 360*G4p*G4x*G5ppp*G5px*pow(H,2) + 612*G4p*G5p*G5ppp*G5px*pow(H,2) - 144*G4*G4p*G5pppx*G5px*pow(H,2) - 
   204*G2x*G4*G5ppx*G5px*pow(H,2) + 408*G3p*G4*G5ppx*G5px*pow(H,2) + 72*pow(G4p,2)*G5ppx*G5px*pow(H,2) + 144*G4*G4pp*G5ppx*G5px*pow(H,2) + 192*G2px*G4*pow(G5px,2)*pow(H,2) - 
   132*G3pp*G4*pow(G5px,2)*pow(H,2) + 210*G2x*G4p*pow(G5px,2)*pow(H,2) - 420*G3p*G4p*pow(G5px,2)*pow(H,2) + 1188*G4p*G4pp*pow(G5px,2)*pow(H,2) - 360*G4*G4ppp*pow(G5px,2)*pow(H,2) - 
   12*G2p*G4x*pow(G5px,2)*pow(H,2) + 138*G2p*G5p*pow(G5px,2)*pow(H,2) - 48*G2x*G3x*G4*G5pxx*pow(H,2) + 96*G3p*G3x*G4*G5pxx*pow(H,2) - 8*G2pxx*pow(G4,2)*G5pxx*pow(H,2) - 
   64*G3ppx*pow(G4,2)*G5pxx*pow(H,2) + 210*G3x*pow(G4p,2)*G5pxx*pow(H,2) + 144*pow(G4,2)*G4pppx*G5pxx*pow(H,2) + 24*G4*G4p*G4ppx*G5pxx*pow(H,2) + 228*G2x*G4*G4px*G5pxx*pow(H,2) - 
   456*G3p*G4*G4px*G5pxx*pow(H,2) - 588*pow(G4p,2)*G4px*G5pxx*pow(H,2) + 360*G4*G4pp*G4px*G5pxx*pow(H,2) + 224*G2px*G4*G4x*G5pxx*pow(H,2) - 456*G3pp*G4*G4x*G5pxx*pow(H,2) - 
   144*G2x*G4p*G4x*G5pxx*pow(H,2) + 288*G3p*G4p*G4x*G5pxx*pow(H,2) + 144*G4*G4ppp*G4x*G5pxx*pow(H,2) - 64*G2p*pow(G4x,2)*G5pxx*pow(H,2) + 120*G2p*G4*G4xx*G5pxx*pow(H,2) - 
   88*G2px*G4*G5p*G5pxx*pow(H,2) + 240*G3pp*G4*G5p*G5pxx*pow(H,2) + 144*G2x*G4p*G5p*G5pxx*pow(H,2) - 288*G3p*G4p*G5p*G5pxx*pow(H,2) + 144*G4p*G4pp*G5p*G5pxx*pow(H,2) - 
   144*G4*G4ppp*G5p*G5pxx*pow(H,2) + 124*G2p*G4x*G5p*G5pxx*pow(H,2) - 46*G2p*pow(G5p,2)*G5pxx*pow(H,2) - 84*G2x*G4*G5pp*G5pxx*pow(H,2) + 168*G3p*G4*G5pp*G5pxx*pow(H,2) + 
   84*pow(G4p,2)*G5pp*G5pxx*pow(H,2) - 216*G4*G4pp*G5pp*G5pxx*pow(H,2) - 72*G4*G4p*G5ppp*G5pxx*pow(H,2) - 60*G2p*G4*G5px*G5pxx*pow(H,2) - 8*G2px*pow(G4,2)*G5pxxx*pow(H,2) + 
   8*G3pp*pow(G4,2)*G5pxxx*pow(H,2) - 12*G2x*G4*G4p*G5pxxx*pow(H,2) + 24*G3p*G4*G4p*G5pxxx*pow(H,2) - 24*G4*G4p*G4pp*G5pxxx*pow(H,2) - 16*G2p*G4*G4x*G5pxxx*pow(H,2) + 
   8*G2p*G4*G5p*G5pxxx*pow(H,2) + 126*G2x*G3pxx*G4*G5x*pow(H,2) - 252*G3p*G3pxx*G4*G5x*pow(H,2) + 144*G2pxx*G3x*G4*G5x*pow(H,2) - 180*G3ppx*G3x*G4*G5x*pow(H,2) - 
   72*G2px*G3xx*G4*G5x*pow(H,2) + 198*G3pp*G3xx*G4*G5x*pow(H,2) + 54*G2x*G3xx*G4p*G5x*pow(H,2) - 108*G3p*G3xx*G4p*G5x*pow(H,2) - 36*G3pxx*pow(G4p,2)*G5x*pow(H,2) - 
   459*pow(G3x,2)*G4pp*G5x*pow(H,2) + 108*G3pxx*G4*G4pp*G5x*pow(H,2) + 72*G3xx*G4p*G4pp*G5x*pow(H,2) - 108*G3xx*G4*G4ppp*G5x*pow(H,2) + 72*G3x*G4*G4pppx*G5x*pow(H,2) - 
   396*G3x*G4p*G4ppx*G5x*pow(H,2) + 432*G4*pow(G4ppx,2)*G5x*pow(H,2) - 300*G2x*G4*G4ppxx*G5x*pow(H,2) + 600*G3p*G4*G4ppxx*G5x*pow(H,2) - 216*G4*G4pp*G4ppxx*G5x*pow(H,2) - 
   72*G2x*G3x*G4px*G5x*pow(H,2) + 144*G3p*G3x*G4px*G5x*pow(H,2) - 360*G2pxx*G4*G4px*G5x*pow(H,2) + 576*G3ppx*G4*G4px*G5x*pow(H,2) + 2268*G3x*G4pp*G4px*G5x*pow(H,2) - 
   432*G4*G4pppx*G4px*G5x*pow(H,2) + 432*G4p*G4ppx*G4px*G5x*pow(H,2) + 204*G2x*pow(G4px,2)*G5x*pow(H,2) - 408*G3p*pow(G4px,2)*G5x*pow(H,2) - 2772*G4pp*pow(G4px,2)*G5x*pow(H,2) + 
   240*G2px*G4*G4pxx*G5x*pow(H,2) - 540*G3pp*G4*G4pxx*G5x*pow(H,2) - 348*G2x*G4p*G4pxx*G5x*pow(H,2) + 696*G3p*G4p*G4pxx*G5x*pow(H,2) - 144*G4p*G4pp*G4pxx*G5x*pow(H,2) + 
   216*G4*G4ppp*G4pxx*G5x*pow(H,2) - 48*G2p*G4*G4pxxx*G5x*pow(H,2) - 234*G2px*G3x*G4x*G5x*pow(H,2) + 558*G3pp*G3x*G4x*G5x*pow(H,2) + 90*G2p*G3xx*G4x*G5x*pow(H,2) + 
   228*G2pxx*G4p*G4x*G5x*pow(H,2) + 132*G3ppx*G4p*G4x*G5x*pow(H,2) - 180*G3x*G4ppp*G4x*G5x*pow(H,2) - 720*G4p*G4pppx*G4x*G5x*pow(H,2) - 48*G2x*G4ppx*G4x*G5x*pow(H,2) + 
   96*G3p*G4ppx*G4x*G5x*pow(H,2) + 1152*G4pp*G4ppx*G4x*G5x*pow(H,2) + 492*G2px*G4px*G4x*G5x*pow(H,2) - 936*G3pp*G4px*G4x*G5x*pow(H,2) - 228*G2p*G4pxx*G4x*G5x*pow(H,2) - 
   96*G2ppx*pow(G4x,2)*G5x*pow(H,2) + 96*G3ppp*pow(G4x,2)*G5x*pow(H,2) - 180*G2p*G3x*G4xx*G5x*pow(H,2) - 48*G2ppx*G4*G4xx*G5x*pow(H,2) + 48*G3ppp*G4*G4xx*G5x*pow(H,2) + 
   348*G2px*G4p*G4xx*G5x*pow(H,2) - 1056*G3pp*G4p*G4xx*G5x*pow(H,2) - 492*G2x*G4pp*G4xx*G5x*pow(H,2) + 984*G3p*G4pp*G4xx*G5x*pow(H,2) - 720*pow(G4pp,2)*G4xx*G5x*pow(H,2) + 
   1008*G4p*G4ppp*G4xx*G5x*pow(H,2) + 336*G2p*G4px*G4xx*G5x*pow(H,2) + 48*G2pp*G4x*G4xx*G5x*pow(H,2) + 48*G2pp*G4*G4xxx*G5x*pow(H,2) - 36*G2p*G4p*G4xxx*G5x*pow(H,2) + 
   180*G2px*G3x*G5p*G5x*pow(H,2) - 522*G3pp*G3x*G5p*G5x*pow(H,2) - 63*G2p*G3xx*G5p*G5x*pow(H,2) - 144*G2pxx*G4p*G5p*G5x*pow(H,2) - 108*G3ppx*G4p*G5p*G5x*pow(H,2) + 
   324*G3x*G4ppp*G5p*G5x*pow(H,2) + 504*G4p*G4pppx*G5p*G5x*pow(H,2) + 156*G2x*G4ppx*G5p*G5x*pow(H,2) - 312*G3p*G4ppx*G5p*G5x*pow(H,2) - 684*G4pp*G4ppx*G5p*G5x*pow(H,2) - 
   372*G2px*G4px*G5p*G5x*pow(H,2) + 942*G3pp*G4px*G5p*G5x*pow(H,2) - 396*G4ppp*G4px*G5p*G5x*pow(H,2) + 198*G2p*G4pxx*G5p*G5x*pow(H,2) + 192*G2ppx*G4x*G5p*G5x*pow(H,2) - 
   192*G3ppp*G4x*G5p*G5x*pow(H,2) - 72*G2pp*G4xx*G5p*G5x*pow(H,2) - 84*G2ppx*pow(G5p,2)*G5x*pow(H,2) + 84*G3ppp*pow(G5p,2)*G5x*pow(H,2) - 99*G2x*G3x*G5pp*G5x*pow(H,2) + 
   198*G3p*G3x*G5pp*G5x*pow(H,2) + 12*G2pxx*G4*G5pp*G5x*pow(H,2) - 12*G3ppx*G4*G5pp*G5x*pow(H,2) - 324*G3x*G4pp*G5pp*G5x*pow(H,2) - 144*G4p*G4ppx*G5pp*G5x*pow(H,2) + 
   216*G2x*G4px*G5pp*G5x*pow(H,2) - 432*G3p*G4px*G5pp*G5x*pow(H,2) + 756*G4pp*G4px*G5pp*G5x*pow(H,2) - 156*G2px*G4x*G5pp*G5x*pow(H,2) + 192*G3pp*G4x*G5pp*G5x*pow(H,2) + 
   24*G2p*G4xx*G5pp*G5x*pow(H,2) + 102*G2px*G5p*G5pp*G5x*pow(H,2) - 132*G3pp*G5p*G5pp*G5x*pow(H,2) - 6*G2x*pow(G5pp,2)*G5x*pow(H,2) + 12*G3p*pow(G5pp,2)*G5x*pow(H,2) + 
   216*G3x*G4p*G5ppp*G5x*pow(H,2) - 360*G4p*G4px*G5ppp*G5x*pow(H,2) + 36*G2x*G4x*G5ppp*G5x*pow(H,2) - 72*G3p*G4x*G5ppp*G5x*pow(H,2) - 30*G2x*G5p*G5ppp*G5x*pow(H,2) + 
   60*G3p*G5p*G5ppp*G5x*pow(H,2) + 24*G2x*G4*G5pppx*G5x*pow(H,2) - 48*G3p*G4*G5pppx*G5x*pow(H,2) + 36*pow(G4p,2)*G5pppx*G5x*pow(H,2) - 72*G2px*G4*G5ppx*G5x*pow(H,2) + 
   96*G3pp*G4*G5ppx*G5x*pow(H,2) + 108*G2x*G4p*G5ppx*G5x*pow(H,2) - 216*G3p*G4p*G5ppx*G5x*pow(H,2) - 36*G4p*G4pp*G5ppx*G5x*pow(H,2) - 24*G2p*G5p*G5ppx*G5x*pow(H,2) + 
   24*G2p*G4*G5ppxx*G5x*pow(H,2) + 81*G2p*G3x*G5px*G5x*pow(H,2) + 48*G2ppx*G4*G5px*G5x*pow(H,2) - 48*G3ppp*G4*G5px*G5x*pow(H,2) - 180*G2px*G4p*G5px*G5x*pow(H,2) + 
   522*G3pp*G4p*G5px*G5x*pow(H,2) + 270*G2x*G4pp*G5px*G5x*pow(H,2) - 540*G3p*G4pp*G5px*G5x*pow(H,2) + 396*pow(G4pp,2)*G5px*G5x*pow(H,2) - 468*G4p*G4ppp*G5px*G5x*pow(H,2) - 
   150*G2p*G4px*G5px*G5x*pow(H,2) + 24*G2pp*G5p*G5px*G5x*pow(H,2) - 12*G2p*G5pp*G5px*G5x*pow(H,2) - 24*G2pp*G4*G5pxx*G5x*pow(H,2) + 30*G2p*G4p*G5pxx*G5x*pow(H,2) - 
   9*G2px*G2x*pow(G5x,2)*pow(H,2) + 18*G2px*G3p*pow(G5x,2)*pow(H,2) + 9*G2x*G3pp*pow(G5x,2)*pow(H,2) - 18*G3p*G3pp*pow(G5x,2)*pow(H,2) + 6*G2pp*G3x*pow(G5x,2)*pow(H,2) + 
   36*G2ppx*G4p*pow(G5x,2)*pow(H,2) - 36*G3ppp*G4p*pow(G5x,2)*pow(H,2) - 54*G2px*G4pp*pow(G5x,2)*pow(H,2) + 72*G3pp*G4pp*pow(G5x,2)*pow(H,2) + 18*G2x*G4ppp*pow(G5x,2)*pow(H,2) - 
   36*G3p*G4ppp*pow(G5x,2)*pow(H,2) + 60*G2px*G3x*G4*G5xx*pow(H,2) - 42*G3pp*G3x*G4*G5xx*pow(H,2) - 18*G2p*G3xx*G4*G5xx*pow(H,2) + 138*G2x*G3x*G4p*G5xx*pow(H,2) - 
   276*G3p*G3x*G4p*G5xx*pow(H,2) - 24*G2pxx*G4*G4p*G5xx*pow(H,2) - 48*G3ppx*G4*G4p*G5xx*pow(H,2) + 552*G3x*G4p*G4pp*G5xx*pow(H,2) - 156*G3x*G4*G4ppp*G5xx*pow(H,2) + 
   144*G4*G4p*G4pppx*G5xx*pow(H,2) + 84*G2x*G4*G4ppx*G5xx*pow(H,2) - 168*G3p*G4*G4ppx*G5xx*pow(H,2) + 252*pow(G4p,2)*G4ppx*G5xx*pow(H,2) - 192*G4*G4pp*G4ppx*G5xx*pow(H,2) - 
   192*G2px*G4*G4px*G5xx*pow(H,2) + 144*G3pp*G4*G4px*G5xx*pow(H,2) - 228*G2x*G4p*G4px*G5xx*pow(H,2) + 456*G3p*G4p*G4px*G5xx*pow(H,2) - 1296*G4p*G4pp*G4px*G5xx*pow(H,2) + 
   384*G4*G4ppp*G4px*G5xx*pow(H,2) + 84*G2p*G4*G4pxx*G5xx*pow(H,2) - 30*G2p*G3x*G4x*G5xx*pow(H,2) - 64*G2ppx*G4*G4x*G5xx*pow(H,2) + 64*G3ppp*G4*G4x*G5xx*pow(H,2) + 
   144*G2px*G4p*G4x*G5xx*pow(H,2) - 372*G3pp*G4p*G4x*G5xx*pow(H,2) - 216*G2x*G4pp*G4x*G5xx*pow(H,2) + 432*G3p*G4pp*G4x*G5xx*pow(H,2) - 312*pow(G4pp,2)*G4x*G5xx*pow(H,2) + 
   264*G4p*G4ppp*G4x*G5xx*pow(H,2) + 24*G2p*G4px*G4x*G5xx*pow(H,2) - 8*G2pp*pow(G4x,2)*G5xx*pow(H,2) - 48*G2pp*G4*G4xx*G5xx*pow(H,2) + 168*G2p*G4p*G4xx*G5xx*pow(H,2) + 
   69*G2p*G3x*G5p*G5xx*pow(H,2) + 8*G2ppx*G4*G5p*G5xx*pow(H,2) - 8*G3ppp*G4*G5p*G5xx*pow(H,2) - 96*G2px*G4p*G5p*G5xx*pow(H,2) + 318*G3pp*G4p*G5p*G5xx*pow(H,2) + 
   144*G2x*G4pp*G5p*G5xx*pow(H,2) - 288*G3p*G4pp*G5p*G5xx*pow(H,2) + 228*pow(G4pp,2)*G5p*G5xx*pow(H,2) - 348*G4p*G4ppp*G5p*G5xx*pow(H,2) - 144*G2p*G4px*G5p*G5xx*pow(H,2) - 
   16*G2pp*G4x*G5p*G5xx*pow(H,2) + 10*G2pp*pow(G5p,2)*G5xx*pow(H,2) + 60*G2px*G4*G5pp*G5xx*pow(H,2) - 72*G3pp*G4*G5pp*G5xx*pow(H,2) - 6*G2x*G4p*G5pp*G5xx*pow(H,2) + 
   12*G3p*G4p*G5pp*G5xx*pow(H,2) + 132*G4p*G4pp*G5pp*G5xx*pow(H,2) + 24*G2p*G4x*G5pp*G5xx*pow(H,2) - 12*G2x*G4*G5ppp*G5xx*pow(H,2) + 24*G3p*G4*G5ppp*G5xx*pow(H,2) - 
   108*pow(G4p,2)*G5ppp*G5xx*pow(H,2) - 24*G2p*G4*G5ppx*G5xx*pow(H,2) + 24*G2pp*G4*G5px*G5xx*pow(H,2) - 96*G2p*G4p*G5px*G5xx*pow(H,2) + 3*G2p*G2x*G5x*G5xx*pow(H,2) - 
   6*G2p*G3p*G5x*G5xx*pow(H,2) - 12*G2pp*G4p*G5x*G5xx*pow(H,2) + 12*G2p*G4pp*G5x*G5xx*pow(H,2) - 6*G2p*G3x*G4*G5xxx*pow(H,2) + 8*G2ppx*pow(G4,2)*G5xxx*pow(H,2) - 
   8*G3ppp*pow(G4,2)*G5xxx*pow(H,2) - 12*G3pp*G4*G4p*G5xxx*pow(H,2) - 12*G2x*pow(G4p,2)*G5xxx*pow(H,2) + 24*G3p*pow(G4p,2)*G5xxx*pow(H,2) + 12*G2x*G4*G4pp*G5xxx*pow(H,2) - 
   24*G3p*G4*G4pp*G5xxx*pow(H,2) - 24*pow(G4p,2)*G4pp*G5xxx*pow(H,2) + 24*G4*pow(G4pp,2)*G5xxx*pow(H,2) + 24*G4*G4p*G4ppp*G5xxx*pow(H,2) + 12*G2p*G4*G4px*G5xxx*pow(H,2) + 
   16*G2pp*G4*G4x*G5xxx*pow(H,2) - 12*G2p*G4p*G4x*G5xxx*pow(H,2) - 8*G2pp*G4*G5p*G5xxx*pow(H,2) + 6*G2p*G4p*G5p*G5xxx*pow(H,2) + 2592*G4*G4pxxx*pow(G4x,2)*pow(H,4) - 
   1728*G4pxx*pow(G4x,3)*pow(H,4) - 864*pow(G4,2)*G4pxxx*G4xx*pow(H,4) - 864*G4*G4pxx*G4x*G4xx*pow(H,4) + 4320*G4px*pow(G4x,2)*G4xx*pow(H,4) - 864*G4*G4px*pow(G4xx,2)*pow(H,4) + 
   11232*G4p*G4x*pow(G4xx,2)*pow(H,4) + 864*pow(G4,2)*G4pxx*G4xxx*pow(H,4) - 1728*G4*G4px*G4x*G4xxx*pow(H,4) + 1728*G4p*pow(G4x,2)*G4xxx*pow(H,4) - 864*G4*G4p*G4xx*G4xxx*pow(H,4) - 
   4320*G4*G4pxxx*G4x*G5p*pow(H,4) + 2160*G4pxx*pow(G4x,2)*G5p*pow(H,4) + 864*G4*G4pxx*G4xx*G5p*pow(H,4) - 5184*G4px*G4x*G4xx*G5p*pow(H,4) - 9504*G4p*pow(G4xx,2)*G5p*pow(H,4) + 
   864*G4*G4px*G4xxx*G5p*pow(H,4) - 3024*G4p*G4x*G4xxx*G5p*pow(H,4) + 1728*G4*G4pxxx*pow(G5p,2)*pow(H,4) + 216*G4pxx*G4x*pow(G5p,2)*pow(H,4) + 1512*G4px*G4xx*pow(G5p,2)*pow(H,4) + 
   1296*G4p*G4xxx*pow(G5p,2)*pow(H,4) - 648*G4pxx*pow(G5p,3)*pow(H,4) - 9072*pow(G4x,2)*G4xx*G5pp*pow(H,4) + 2592*G4*pow(G4xx,2)*G5pp*pow(H,4) + 2160*G4*G4x*G4xxx*G5pp*pow(H,4) + 
   13608*G4x*G4xx*G5p*G5pp*pow(H,4) - 1296*G4*G4xxx*G5p*G5pp*pow(H,4) - 5184*G4xx*pow(G5p,2)*G5pp*pow(H,4) + 1776*pow(G4x,3)*G5ppx*pow(H,4) - 960*G4*G4x*G4xx*G5ppx*pow(H,4) - 
   336*pow(G4,2)*G4xxx*G5ppx*pow(H,4) - 3240*pow(G4x,2)*G5p*G5ppx*pow(H,4) + 144*G4*G4xx*G5p*G5ppx*pow(H,4) + 1608*G4x*pow(G5p,2)*G5ppx*pow(H,4) - 144*pow(G5p,3)*G5ppx*pow(H,4) - 
   1296*G4*pow(G4x,2)*G5ppxx*pow(H,4) + 624*pow(G4,2)*G4xx*G5ppxx*pow(H,4) + 2112*G4*G4x*G5p*G5ppxx*pow(H,4) - 816*G4*pow(G5p,2)*G5ppxx*pow(H,4) + 336*pow(G4,2)*G4pxxx*G5px*pow(H,4) - 
   612*G3xx*G4*G4x*G5px*pow(H,4) + 2616*G4*G4pxx*G4x*G5px*pow(H,4) + 1584*G3x*pow(G4x,2)*G5px*pow(H,4) - 6168*G4px*pow(G4x,2)*G5px*pow(H,4) - 2232*G3x*G4*G4xx*G5px*pow(H,4) + 
   6888*G4*G4px*G4xx*G5px*pow(H,4) - 15096*G4p*G4x*G4xx*G5px*pow(H,4) + 504*G4*G4p*G4xxx*G5px*pow(H,4) + 396*G3xx*G4*G5p*G5px*pow(H,4) - 1080*G4*G4pxx*G5p*G5px*pow(H,4) - 
   2430*G3x*G4x*G5p*G5px*pow(H,4) + 8256*G4px*G4x*G5p*G5px*pow(H,4) + 12288*G4p*G4xx*G5p*G5px*pow(H,4) + 846*G3x*pow(G5p,2)*G5px*pow(H,4) - 2652*G4px*pow(G5p,2)*G5px*pow(H,4) + 
   5064*pow(G4x,2)*G5pp*G5px*pow(H,4) - 3000*G4*G4xx*G5pp*G5px*pow(H,4) - 7524*G4x*G5p*G5pp*G5px*pow(H,4) + 2952*pow(G5p,2)*G5pp*G5px*pow(H,4) - 144*G4*G5p*G5ppx*G5px*pow(H,4) - 
   288*pow(G4,2)*G5ppxx*G5px*pow(H,4) + 1242*G3x*G4*pow(G5px,2)*pow(H,4) - 3876*G4*G4px*pow(G5px,2)*pow(H,4) + 4728*G4p*G4x*pow(G5px,2)*pow(H,4) - 3846*G4p*G5p*pow(G5px,2)*pow(H,4) + 
   1032*G4*G5pp*pow(G5px,2)*pow(H,4) + 72*G3xx*pow(G4,2)*G5pxx*pow(H,4) - 1104*pow(G4,2)*G4pxx*G5pxx*pow(H,4) + 648*G3x*G4*G4x*G5pxx*pow(H,4) - 552*G4*G4px*G4x*G5pxx*pow(H,4) - 
   624*G4p*pow(G4x,2)*G5pxx*pow(H,4) + 600*G4*G4p*G4xx*G5pxx*pow(H,4) - 288*G3x*G4*G5p*G5pxx*pow(H,4) + 456*G4*G4px*G5p*G5pxx*pow(H,4) + 1728*G4p*G4x*G5p*G5pxx*pow(H,4) - 
   936*G4p*pow(G5p,2)*G5pxx*pow(H,4) - 912*G4*G4x*G5pp*G5pxx*pow(H,4) + 264*G4*G5p*G5pp*G5pxx*pow(H,4) + 456*pow(G4,2)*G5ppx*G5pxx*pow(H,4) - 168*G4*G4p*G5px*G5pxx*pow(H,4) - 
   72*G3x*pow(G4,2)*G5pxxx*pow(H,4) + 144*pow(G4,2)*G4px*G5pxxx*pow(H,4) - 312*G4*G4p*G4x*G5pxxx*pow(H,4) + 216*G4*G4p*G5p*G5pxxx*pow(H,4) + 24*pow(G4,2)*G5pp*G5pxxx*pow(H,4) - 
   324*G3xx*G4*G4px*G5x*pow(H,4) + 648*G3x*G4*G4pxx*G5x*pow(H,4) - 1080*G4*G4px*G4pxx*G5x*pow(H,4) - 1008*G4*G4p*G4pxxx*G5x*pow(H,4) + 1404*G3pxx*G4*G4x*G5x*pow(H,4) + 
   1080*G3xx*G4p*G4x*G5x*pow(H,4) - 2808*G4*G4ppxx*G4x*G5x*pow(H,4) + 1512*G3x*G4px*G4x*G5x*pow(H,4) - 3456*pow(G4px,2)*G4x*G5x*pow(H,4) - 1728*G4p*G4pxx*G4x*G5x*pow(H,4) + 
   5832*G4ppx*pow(G4x,2)*G5x*pow(H,4) + 4752*G3x*G4p*G4xx*G5x*pow(H,4) + 792*G4*G4ppx*G4xx*G5x*pow(H,4) - 14400*G4p*G4px*G4xx*G5x*pow(H,4) - 10584*G4pp*G4x*G4xx*G5x*pow(H,4) - 
   648*pow(G4p,2)*G4xxx*G5x*pow(H,4) + 360*G4*G4pp*G4xxx*G5x*pow(H,4) - 1188*G3pxx*G4*G5p*G5x*pow(H,4) - 918*G3xx*G4p*G5p*G5x*pow(H,4) + 2088*G4*G4ppxx*G5p*G5x*pow(H,4) - 
   972*G3x*G4px*G5p*G5x*pow(H,4) + 1980*pow(G4px,2)*G5p*G5x*pow(H,4) + 2628*G4p*G4pxx*G5p*G5x*pow(H,4) - 9360*G4ppx*G4x*G5p*G5x*pow(H,4) + 9720*G4pp*G4xx*G5p*G5x*pow(H,4) + 
   3636*G4ppx*pow(G5p,2)*G5x*pow(H,4) + 486*G3xx*G4*G5pp*G5x*pow(H,4) + 36*G4*G4pxx*G5pp*G5x*pow(H,4) - 3510*G3x*G4x*G5pp*G5x*pow(H,4) + 9324*G4px*G4x*G5pp*G5x*pow(H,4) + 
   3456*G4p*G4xx*G5pp*G5x*pow(H,4) + 2754*G3x*G5p*G5pp*G5x*pow(H,4) - 7362*G4px*G5p*G5pp*G5x*pow(H,4) - 324*G4x*pow(G5pp,2)*G5x*pow(H,4) + 432*G5p*pow(G5pp,2)*G5x*pow(H,4) - 
   360*pow(G4x,2)*G5ppp*G5x*pow(H,4) - 1008*G4*G4xx*G5ppp*G5x*pow(H,4) + 972*G4x*G5p*G5ppp*G5x*pow(H,4) - 648*pow(G5p,2)*G5ppp*G5x*pow(H,4) + 144*G4*G5p*G5pppx*G5x*pow(H,4) - 
   636*G3x*G4*G5ppx*G5x*pow(H,4) + 1692*G4*G4px*G5ppx*G5x*pow(H,4) - 1200*G4p*G4x*G5ppx*G5x*pow(H,4) + 444*G4p*G5p*G5ppx*G5x*pow(H,4) - 540*G4*G5pp*G5ppx*G5x*pow(H,4) + 
   432*G4*G4p*G5ppxx*G5x*pow(H,4) - 3168*G3x*G4p*G5px*G5x*pow(H,4) - 684*G4*G4ppx*G5px*G5x*pow(H,4) + 9192*G4p*G4px*G5px*G5x*pow(H,4) + 390*G2x*G4x*G5px*G5x*pow(H,4) - 
   780*G3p*G4x*G5px*G5x*pow(H,4) + 5880*G4pp*G4x*G5px*G5x*pow(H,4) - 255*G2x*G5p*G5px*G5x*pow(H,4) + 510*G3p*G5p*G5px*G5x*pow(H,4) - 5394*G4pp*G5p*G5px*G5x*pow(H,4) - 
   2046*G4p*G5pp*G5px*G5x*pow(H,4) + 540*G4*G5ppp*G5px*G5x*pow(H,4) + 108*G2x*G4*G5pxx*G5x*pow(H,4) - 216*G3p*G4*G5pxx*G5x*pow(H,4) + 438*pow(G4p,2)*G5pxx*G5x*pow(H,4) + 
   72*G4*G4pp*G5pxx*G5x*pow(H,4) + 126*G2pxx*G4*pow(G5x,2)*pow(H,4) - 108*G3ppx*G4*pow(G5x,2)*pow(H,4) - 1179*G3x*G4pp*pow(G5x,2)*pow(H,4) - 36*G4*G4pppx*pow(G5x,2)*pow(H,4) - 
   1674*G4p*G4ppx*pow(G5x,2)*pow(H,4) + 54*G2x*G4px*pow(G5x,2)*pow(H,4) - 108*G3p*G4px*pow(G5x,2)*pow(H,4) + 3006*G4pp*G4px*pow(G5x,2)*pow(H,4) - 450*G2px*G4x*pow(G5x,2)*pow(H,4) + 
   1152*G3pp*G4x*pow(G5x,2)*pow(H,4) - 504*G4ppp*G4x*pow(G5x,2)*pow(H,4) - 198*G2p*G4xx*pow(G5x,2)*pow(H,4) + 378*G2px*G5p*pow(G5x,2)*pow(H,4) - 1062*G3pp*G5p*pow(G5x,2)*pow(H,4) + 
   612*G4ppp*G5p*pow(G5x,2)*pow(H,4) - 243*G2x*G5pp*pow(G5x,2)*pow(H,4) + 486*G3p*G5pp*pow(G5x,2)*pow(H,4) - 378*G4pp*G5pp*pow(G5x,2)*pow(H,4) + 306*G4p*G5ppp*pow(G5x,2)*pow(H,4) + 
   57*G2p*G5px*pow(G5x,2)*pow(H,4) + 15*G2pp*pow(G5x,3)*pow(H,4) - 72*G3pxx*pow(G4,2)*G5xx*pow(H,4) - 108*G3xx*G4*G4p*G5xx*pow(H,4) + 480*pow(G4,2)*G4ppxx*G5xx*pow(H,4) - 
   396*G3x*G4*G4px*G5xx*pow(H,4) + 1248*G4*pow(G4px,2)*G5xx*pow(H,4) + 408*G4*G4p*G4pxx*G5xx*pow(H,4) + 1872*G3x*G4p*G4x*G5xx*pow(H,4) + 1320*G4*G4ppx*G4x*G5xx*pow(H,4) - 
   5280*G4p*G4px*G4x*G5xx*pow(H,4) - 2016*G4pp*pow(G4x,2)*G5xx*pow(H,4) - 1416*pow(G4p,2)*G4xx*G5xx*pow(H,4) - 312*G4*G4pp*G4xx*G5xx*pow(H,4) - 1638*G3x*G4p*G5p*G5xx*pow(H,4) - 
   744*G4*G4ppx*G5p*G5xx*pow(H,4) + 4644*G4p*G4px*G5p*G5xx*pow(H,4) + 3408*G4pp*G4x*G5p*G5xx*pow(H,4) - 1560*G4pp*pow(G5p,2)*G5xx*pow(H,4) + 630*G3x*G4*G5pp*G5xx*pow(H,4) - 
   1944*G4*G4px*G5pp*G5xx*pow(H,4) + 1152*G4p*G4x*G5pp*G5xx*pow(H,4) - 966*G4p*G5p*G5pp*G5xx*pow(H,4) + 120*G4*pow(G5pp,2)*G5xx*pow(H,4) - 216*G4*G4x*G5ppp*G5xx*pow(H,4) + 
   384*G4*G5p*G5ppp*G5xx*pow(H,4) - 168*pow(G4,2)*G5pppx*G5xx*pow(H,4) - 156*G4*G4p*G5ppx*G5xx*pow(H,4) - 204*G2x*G4*G5px*G5xx*pow(H,4) + 408*G3p*G4*G5px*G5xx*pow(H,4) + 
   822*pow(G4p,2)*G5px*G5xx*pow(H,4) - 24*G4*G4pp*G5px*G5xx*pow(H,4) + 162*G3pp*G4*G5x*G5xx*pow(H,4) + 324*G2x*G4p*G5x*G5xx*pow(H,4) - 648*G3p*G4p*G5x*G5xx*pow(H,4) + 
   1284*G4p*G4pp*G5x*G5xx*pow(H,4) - 324*G4*G4ppp*G5x*G5xx*pow(H,4) - 60*G2p*G4x*G5x*G5xx*pow(H,4) + 141*G2p*G5p*G5x*G5xx*pow(H,4) - 42*G2p*G4*pow(G5xx,2)*pow(H,4) - 
   36*G3x*G4*G4p*G5xxx*pow(H,4) - 144*pow(G4,2)*G4ppx*G5xxx*pow(H,4) + 48*G4*G4p*G4px*G5xxx*pow(H,4) - 192*pow(G4p,2)*G4x*G5xxx*pow(H,4) + 168*G4*G4pp*G4x*G5xxx*pow(H,4) + 
   156*pow(G4p,2)*G5p*G5xxx*pow(H,4) - 72*G4*G4pp*G5p*G5xxx*pow(H,4) - 12*G4*G4p*G5pp*G5xxx*pow(H,4) - 24*pow(G4,2)*G5ppp*G5xxx*pow(H,4) - 18*G2p*G4*G5x*G5xxx*pow(H,4) + 
   2700*pow(G4x,2)*G5px*G5x*pow(H,6) - 2520*G4*G4xx*G5px*G5x*pow(H,6) - 4950*G4x*G5p*G5px*G5x*pow(H,6) + 2250*pow(G5p,2)*G5px*G5x*pow(H,6) + 1410*G4*pow(G5px,2)*G5x*pow(H,6) + 
   1872*G4*G4x*G5pxx*G5x*pow(H,6) - 1512*G4*G5p*G5pxx*G5x*pow(H,6) - 72*pow(G4,2)*G5pxxx*G5x*pow(H,6) + 1404*G4*G4pxx*pow(G5x,2)*pow(H,6) + 1944*G4px*G4x*pow(G5x,2)*pow(H,6) + 
   3456*G4p*G4xx*pow(G5x,2)*pow(H,6) - 1836*G4px*G5p*pow(G5x,2)*pow(H,6) - 3024*G4x*G5pp*pow(G5x,2)*pow(H,6) + 2754*G5p*G5pp*pow(G5x,2)*pow(H,6) - 960*G4*G5ppx*pow(G5x,2)*pow(H,6) - 
   2727*G4p*G5px*pow(G5x,2)*pow(H,6) - 765*G4pp*pow(G5x,3)*pow(H,6) - 1464*G4*G4x*G5px*G5xx*pow(H,6) + 1104*G4*G5p*G5px*G5xx*pow(H,6) - 540*G4*G4px*G5x*G5xx*pow(H,6) + 
   2988*G4p*G4x*G5x*G5xx*pow(H,6) - 2502*G4p*G5p*G5x*G5xx*pow(H,6) + 1026*G4*G5pp*G5x*G5xx*pow(H,6) - 252*G4*G4p*pow(G5xx,2)*pow(H,6) + 72*pow(G4,2)*G5px*G5xxx*pow(H,6) - 
   108*G4*G4p*G5x*G5xxx*pow(H,6) - 18*G3x*G3xx*G4px*rptot + 36*G3xx*pow(G4px,2)*rptot + 18*pow(G3x,2)*G4pxx*rptot - 36*G3x*G4px*G4pxx*rptot - 18*G3pxx*G3x*G4x*rptot - 36*G3xx*G4ppx*G4x*rptot + 
   36*G3x*G4ppxx*G4x*rptot + 36*G3pxx*G4px*G4x*rptot - 72*G4ppxx*G4px*G4x*rptot + 72*G4ppx*G4pxx*G4x*rptot - 72*G3x*G4ppx*G4xx*rptot + 144*G4ppx*G4px*G4xx*rptot - 24*G2pxx*G4x*G4xx*rptot + 
   24*G3ppx*G4x*G4xx*rptot + 9*G3pxx*G3x*G5p*rptot + 18*G3xx*G4ppx*G5p*rptot - 18*G3x*G4ppxx*G5p*rptot - 18*G3pxx*G4px*G5p*rptot + 36*G4ppxx*G4px*G5p*rptot - 36*G4ppx*G4pxx*G5p*rptot + 
   12*G2pxx*G4xx*G5p*rptot - 12*G3ppx*G4xx*G5p*rptot + 9*G3x*G3xx*G5pp*rptot - 18*G3xx*G4px*G5pp*rptot - 18*G3x*G4pxx*G5pp*rptot + 36*G4px*G4pxx*G5pp*rptot - 9*pow(G3x,2)*G5ppx*rptot + 
   36*G3x*G4px*G5ppx*rptot - 36*pow(G4px,2)*G5ppx*rptot + 36*G3x*G4ppx*G5px*rptot - 72*G4ppx*G4px*G5px*rptot + 12*G2pxx*G4x*G5px*rptot - 12*G3ppx*G4x*G5px*rptot - 6*G2pxx*G5p*G5px*rptot + 
   6*G3ppx*G5p*G5px*rptot - 3*G2pxx*G3x*G5x*rptot + 3*G3ppx*G3x*G5x*rptot + 6*G2pxx*G4px*G5x*rptot - 6*G3ppx*G4px*G5x*rptot - 144*G4pxxx*pow(G4x,2)*pow(H,2)*rptot + 
   144*G4*G4pxxx*G4xx*pow(H,2)*rptot + 576*G4pxx*G4x*G4xx*pow(H,2)*rptot - 288*G4px*pow(G4xx,2)*pow(H,2)*rptot - 144*G4*G4pxx*G4xxx*pow(H,2)*rptot + 144*G4p*G4xx*G4xxx*pow(H,2)*rptot + 
   216*G4pxxx*G4x*G5p*pow(H,2)*rptot - 576*G4pxx*G4xx*G5p*pow(H,2)*rptot + 72*G4px*G4xxx*G5p*pow(H,2)*rptot - 72*G4pxxx*pow(G5p,2)*pow(H,2)*rptot - 72*G4x*G4xxx*G5pp*pow(H,2)*rptot - 
   216*G4x*G4xx*G5ppx*pow(H,2)*rptot + 72*G4*G4xxx*G5ppx*pow(H,2)*rptot + 252*G4xx*G5p*G5ppx*pow(H,2)*rptot + 72*pow(G4x,2)*G5ppxx*pow(H,2)*rptot - 72*G4*G4xx*G5ppxx*pow(H,2)*rptot - 
   108*G4x*G5p*G5ppxx*pow(H,2)*rptot + 36*pow(G5p,2)*G5ppxx*pow(H,2)*rptot - 72*G4*G4pxxx*G5px*pow(H,2)*rptot + 18*G3xx*G4x*G5px*pow(H,2)*rptot - 396*G4pxx*G4x*G5px*pow(H,2)*rptot - 
   36*G3x*G4xx*G5px*pow(H,2)*rptot + 504*G4px*G4xx*G5px*pow(H,2)*rptot - 36*G4p*G4xxx*G5px*pow(H,2)*rptot + 9*G3xx*G5p*G5px*pow(H,2)*rptot + 306*G4pxx*G5p*G5px*pow(H,2)*rptot - 
   36*G4xx*G5pp*G5px*pow(H,2)*rptot + 144*G4x*G5ppx*G5px*pow(H,2)*rptot - 144*G5p*G5ppx*G5px*pow(H,2)*rptot + 36*G4*G5ppxx*G5px*pow(H,2)*rptot + 27*G3x*pow(G5px,2)*pow(H,2)*rptot - 
   198*G4px*pow(G5px,2)*pow(H,2)*rptot + 18*G5pp*pow(G5px,2)*pow(H,2)*rptot - 18*G3xx*G4*G5pxx*pow(H,2)*rptot + 108*G4*G4pxx*G5pxx*pow(H,2)*rptot + 30*G3x*G4x*G5pxx*pow(H,2)*rptot - 
   96*G4px*G4x*G5pxx*pow(H,2)*rptot - 144*G4p*G4xx*G5pxx*pow(H,2)*rptot - 33*G3x*G5p*G5pxx*pow(H,2)*rptot + 48*G4px*G5p*G5pxx*pow(H,2)*rptot + 36*G4x*G5pp*G5pxx*pow(H,2)*rptot - 
   36*G4*G5ppx*G5pxx*pow(H,2)*rptot + 54*G4p*G5px*G5pxx*pow(H,2)*rptot + 6*G3x*G4*G5pxxx*pow(H,2)*rptot - 12*G4*G4px*G5pxxx*pow(H,2)*rptot + 12*G4p*G4x*G5pxxx*pow(H,2)*rptot - 
   6*G4p*G5p*G5pxxx*pow(H,2)*rptot - 18*G3xx*G4px*G5x*pow(H,2)*rptot + 108*G3x*G4pxx*G5x*pow(H,2)*rptot - 252*G4px*G4pxx*G5x*pow(H,2)*rptot + 36*G4p*G4pxxx*G5x*pow(H,2)*rptot - 
   90*G3pxx*G4x*G5x*pow(H,2)*rptot + 180*G4ppxx*G4x*G5x*pow(H,2)*rptot - 144*G4ppx*G4xx*G5x*pow(H,2)*rptot - 36*G4pp*G4xxx*G5x*pow(H,2)*rptot + 63*G3pxx*G5p*G5x*pow(H,2)*rptot - 
   126*G4ppxx*G5p*G5x*pow(H,2)*rptot - 9*G3xx*G5pp*G5x*pow(H,2)*rptot + 18*G4pxx*G5pp*G5x*pow(H,2)*rptot - 45*G3x*G5ppx*G5x*pow(H,2)*rptot + 126*G4px*G5ppx*G5x*pow(H,2)*rptot - 
   18*G4p*G5ppxx*G5x*pow(H,2)*rptot + 90*G4ppx*G5px*G5x*pow(H,2)*rptot + 3*G2x*G5pxx*G5x*pow(H,2)*rptot - 6*G3p*G5pxx*G5x*pow(H,2)*rptot + 18*G4pp*G5pxx*G5x*pow(H,2)*rptot - 
   9*G2pxx*pow(G5x,2)*pow(H,2)*rptot + 9*G3ppx*pow(G5x,2)*pow(H,2)*rptot + 18*G3pxx*G4*G5xx*pow(H,2)*rptot + 18*G3xx*G4p*G5xx*pow(H,2)*rptot - 36*G4*G4ppxx*G5xx*pow(H,2)*rptot - 
   42*G3x*G4px*G5xx*pow(H,2)*rptot + 120*pow(G4px,2)*G5xx*pow(H,2)*rptot - 108*G4p*G4pxx*G5xx*pow(H,2)*rptot - 48*G4ppx*G4x*G5xx*pow(H,2)*rptot - 72*G4pp*G4xx*G5xx*pow(H,2)*rptot + 
   60*G4ppx*G5p*G5xx*pow(H,2)*rptot + 3*G3x*G5pp*G5xx*pow(H,2)*rptot - 24*G4px*G5pp*G5xx*pow(H,2)*rptot + 36*G4p*G5ppx*G5xx*pow(H,2)*rptot - 3*G2x*G5px*G5xx*pow(H,2)*rptot + 
   6*G3p*G5px*G5xx*pow(H,2)*rptot + 36*G4pp*G5px*G5xx*pow(H,2)*rptot - 3*G2px*G5x*G5xx*pow(H,2)*rptot + 6*G3pp*G5x*G5xx*pow(H,2)*rptot + 6*G3x*G4p*G5xxx*pow(H,2)*rptot + 
   12*G4*G4ppx*G5xxx*pow(H,2)*rptot - 12*G4pp*G4x*G5xxx*pow(H,2)*rptot + 6*G4pp*G5p*G5xxx*pow(H,2)*rptot - 6*G4p*G5pp*G5xxx*pow(H,2)*rptot + 180*G4xx*G5px*G5x*pow(H,4)*rptot - 
   90*pow(G5px,2)*G5x*pow(H,4)*rptot - 48*G4x*G5pxx*G5x*pow(H,4)*rptot + 3*G5p*G5pxx*G5x*pow(H,4)*rptot + 18*G4*G5pxxx*G5x*pow(H,4)*rptot + 18*G4pxx*pow(G5x,2)*pow(H,4)*rptot + 
   18*G5ppx*pow(G5x,2)*pow(H,4)*rptot + 96*G4x*G5px*G5xx*pow(H,4)*rptot - 51*G5p*G5px*G5xx*pow(H,4)*rptot + 12*G4px*G5x*G5xx*pow(H,4)*rptot - 93*G5pp*G5x*G5xx*pow(H,4)*rptot + 
   42*G4p*pow(G5xx,2)*pow(H,4)*rptot - 18*G4*G5px*G5xxx*pow(H,4)*rptot + 18*G4p*G5x*G5xxx*pow(H,4)*rptot - 
   3*pow(G3px,2)*(6*G3x*G4 + 20*G4p*G4x - 10*G4p*G5p - 2*G4*(2*G4px + 2*G5pp + 21*G5x*pow(H,2)) + G5x*rptot) + 
   G3px*(-63*pow(G3x,2)*G4p - 36*G3xx*G4*G4pp + 84*G3x*G4*G4ppx + 372*G3x*G4p*G4px - 120*G4*G4ppx*G4px - 444*G4p*pow(G4px,2) + 120*G4*G4pp*G4pxx + 48*G2pxx*G4*G4x + 12*G3x*G4pp*G4x - 
      96*G4*G4pppx*G4x + 216*G4p*G4ppx*G4x - 144*G4pp*G4px*G4x - 48*G2px*pow(G4x,2) + 120*G3pp*pow(G4x,2) - 48*G4ppp*pow(G4x,2) - 48*G2px*G4*G4xx + 24*G3pp*G4*G4xx - 192*G4p*G4pp*G4xx + 
      48*G4*G4ppp*G4xx - 72*G2p*G4x*G4xx - 24*G2pxx*G4*G5p - 42*G3x*G4pp*G5p + 48*G4*G4pppx*G5p - 108*G4p*G4ppx*G5p + 120*G4pp*G4px*G5p + 24*G2px*G4x*G5p - 108*G3pp*G4x*G5p + 72*G4ppp*G4x*G5p + 
      36*G2p*G4xx*G5p + 24*G3pp*pow(G5p,2) - 24*G4ppp*pow(G5p,2) - 60*G3x*G4p*G5pp - 24*G4*G4ppx*G5pp + 72*G4p*G4px*G5pp + 48*G4pp*G4x*G5pp - 12*G4pp*G5p*G5pp + 12*G4p*pow(G5pp,2) + 
      12*G3x*G4*G5ppp - 24*G4*G4px*G5ppp + 24*G4p*G4x*G5ppp - 12*G4p*G5p*G5ppp - 24*G4*G4pp*G5ppx + 24*G2px*G4*G5px - 12*G3pp*G4*G5px + 84*G4p*G4pp*G5px - 24*G4*G4ppp*G5px + 36*G2p*G4x*G5px - 
      18*G2p*G5p*G5px + 3*pow(G2x,2)*G5x + 12*pow(G3p,2)*G5x - 9*G2p*G3x*G5x + 8*G2ppx*G4*G5x - 8*G3ppp*G4*G5x - 12*G2px*G4p*G5x + 6*G3pp*G4p*G5x + 12*pow(G4pp,2)*G5x + 12*G4p*G4ppp*G5x + 
      18*G2p*G4px*G5x + 8*G2pp*G4x*G5x - 4*G2pp*G5p*G5x - 144*pow(G4,2)*G4pxxx*pow(H,2) - 540*G3xx*G4*G4x*pow(H,2) + 2376*G4*G4pxx*G4x*pow(H,2) + 324*G3x*pow(G4x,2)*pow(H,2) - 
      1008*G4px*pow(G4x,2)*pow(H,2) - 216*G3x*G4*G4xx*pow(H,2) + 72*G4*G4px*G4xx*pow(H,2) - 720*G4p*G4x*G4xx*pow(H,2) + 72*G4*G4p*G4xxx*pow(H,2) + 324*G3xx*G4*G5p*pow(H,2) - 
      936*G4*G4pxx*G5p*pow(H,2) - 702*G3x*G4x*G5p*pow(H,2) + 1800*G4px*G4x*G5p*pow(H,2) + 792*G4p*G4xx*G5p*pow(H,2) + 378*G3x*pow(G5p,2)*pow(H,2) - 972*G4px*pow(G5p,2)*pow(H,2) - 
      216*pow(G4x,2)*G5pp*pow(H,2) + 648*G4*G4xx*G5pp*pow(H,2) + 432*G4x*G5p*G5pp*pow(H,2) - 108*pow(G5p,2)*G5pp*pow(H,2) - 784*G4*G4x*G5ppx*pow(H,2) + 224*G4*G5p*G5ppx*pow(H,2) + 
      64*pow(G4,2)*G5ppxx*pow(H,2) + 72*G3x*G4*G5px*pow(H,2) - 168*G4*G4px*G5px*pow(H,2) + 252*G4p*G4x*G5px*pow(H,2) - 432*G4p*G5p*G5px*pow(H,2) - 228*G4*G5pp*G5px*pow(H,2) + 
      48*G4*G4p*G5pxx*pow(H,2) - 216*G3x*G4p*G5x*pow(H,2) - 540*G4*G4ppx*G5x*pow(H,2) + 600*G4p*G4px*G5x*pow(H,2) - 348*G4pp*G4x*G5x*pow(H,2) + 198*G4pp*G5p*G5x*pow(H,2) + 
      30*G4p*G5pp*G5x*pow(H,2) + 12*G4*G5ppp*G5x*pow(H,2) - 15*G2p*pow(G5x,2)*pow(H,2) + 6*pow(G4p,2)*G5xx*pow(H,2) + 72*G4*G4pp*G5xx*pow(H,2) - 864*pow(G4x,2)*G5x*pow(H,4) - 
      216*G4*G4xx*G5x*pow(H,4) + 1026*G4x*G5p*G5x*pow(H,4) - 162*pow(G5p,2)*G5x*pow(H,4) + 240*G4*G5px*G5x*pow(H,4) + 135*G4p*pow(G5x,2)*pow(H,4) - 576*G4*G4x*G5xx*pow(H,4) + 
      216*G4*G5p*G5xx*pow(H,4) + 72*pow(G4,2)*G5xxx*pow(H,4) + 6*G3p*(6*G3xx*G4 - 30*G3x*G4x + 80*G4px*G4x + 32*G4p*G4xx + 21*G3x*G5p - 48*G4px*G5p - 8*G4x*G5pp + 2*G5p*G5pp - 14*G4p*G5px - 
         4*G4pp*G5x - 42*G4x*G5x*pow(H,2) + 45*G5p*G5x*pow(H,2) + 4*G4*(-5*G4pxx + G5ppx + 4*G5xx*pow(H,2))) - 
      3*G2x*(6*G3xx*G4 - 30*G3x*G4x + 80*G4px*G4x + 32*G4p*G4xx + 21*G3x*G5p - 48*G4px*G5p - 8*G4x*G5pp + 2*G5p*G5pp - 14*G4p*G5px + 4*G3p*G5x - 4*G4pp*G5x - 42*G4x*G5x*pow(H,2) + 
         45*G5p*G5x*pow(H,2) + 4*G4*(-5*G4pxx + G5ppx + 4*G5xx*pow(H,2))) + 18*G3xx*G4x*rptot - 60*G4pxx*G4x*rptot + 36*G3x*G4xx*rptot - 48*G4px*G4xx*rptot - 9*G3xx*G5p*rptot + 30*G4pxx*G5p*rptot - 
      12*G4xx*G5pp*rptot + 12*G4x*G5ppx*rptot - 6*G5p*G5ppx*rptot - 15*G3x*G5px*rptot + 18*G4px*G5px*rptot + 6*G5pp*G5px*rptot + 6*G4ppx*G5x*rptot + 108*G4xx*G5x*pow(H,2)*rptot - 
      63*G5px*G5x*pow(H,2)*rptot + 42*G4x*G5xx*pow(H,2)*rptot - 39*G5p*G5xx*pow(H,2)*rptot - 6*G4*G5xxx*pow(H,2)*rptot) - 
   G2xx*(-18*pow(G3x,2)*G4p + 168*G4*G4ppx*G4px + 24*G4p*pow(G4px,2) + 24*G2x*G4*G4pxx - 48*G3p*G4*G4pxx + 48*G4*G4pp*G4pxx + 48*G3ppx*G4*G4x - 96*G4*G4pppx*G4x - 72*G4p*G4ppx*G4x - 
      72*G4pp*G4px*G4x - 24*G2px*pow(G4x,2) + 72*G3pp*pow(G4x,2) - 48*G4ppp*pow(G4x,2) - 24*G3pp*G4*G4xx - 24*G2x*G4p*G4xx + 48*G3p*G4p*G4xx - 48*G4p*G4pp*G4xx + 48*G4*G4ppp*G4xx - 
      24*G2p*G4x*G4xx - 24*G3ppx*G4*G5p + 48*G4*G4pppx*G5p + 36*G4p*G4ppx*G5p - 12*G2x*G4px*G5p + 24*G3p*G4px*G5p + 12*G4pp*G4px*G5p + 24*G2px*G4x*G5p - 84*G3pp*G4x*G5p + 72*G4ppp*G4x*G5p + 
      12*G2p*G4xx*G5p - 6*G2px*pow(G5p,2) + 24*G3pp*pow(G5p,2) - 24*G4ppp*pow(G5p,2) - 24*G4*G4ppx*G5pp - 72*G4p*G4px*G5pp + 12*G2x*G4x*G5pp - 24*G3p*G4x*G5pp + 48*G4pp*G4x*G5pp - 
      12*G4pp*G5p*G5pp + 12*G4p*pow(G5pp,2) - 24*G4*G4px*G5ppp + 24*G4p*G4x*G5ppp - 12*G4p*G5p*G5ppp - 12*G2x*G4*G5ppx + 24*G3p*G4*G5ppx - 24*G4*G4pp*G5ppx + 12*G3pp*G4*G5px + 6*G2x*G4p*G5px - 
      12*G3p*G4p*G5px + 12*G4p*G4pp*G5px - 24*G4*G4ppp*G5px + 12*G2p*G4x*G5px - 6*G2p*G5p*G5px + 8*G2ppx*G4*G5x - 8*G3ppp*G4*G5x - 6*G3pp*G4p*G5x + 6*G2x*G4pp*G5x - 12*G3p*G4pp*G5x + 
      12*pow(G4pp,2)*G5x + 12*G4p*G4ppp*G5x + 6*G2p*G4px*G5x + 8*G2pp*G4x*G5x - 4*G2pp*G5p*G5x + 432*G4*G4pxx*G4x*pow(H,2) - 432*G4p*G4x*G4xx*pow(H,2) - 288*G4*G4pxx*G5p*pow(H,2) - 
      144*G4px*G4x*G5p*pow(H,2) + 360*G4p*G4xx*G5p*pow(H,2) + 108*G4px*pow(G5p,2)*pow(H,2) + 144*pow(G4x,2)*G5pp*pow(H,2) - 72*G4*G4xx*G5pp*pow(H,2) - 108*G4x*G5p*G5pp*pow(H,2) - 
      208*G4*G4x*G5ppx*pow(H,2) + 152*G4*G5p*G5ppx*pow(H,2) - 8*pow(G4,2)*G5ppxx*pow(H,2) - 180*G4*G4px*G5px*pow(H,2) + 192*G4p*G4x*G5px*pow(H,2) - 150*G4p*G5p*G5px*pow(H,2) + 
      96*G4*G5pp*G5px*pow(H,2) - 12*G4*G4p*G5pxx*pow(H,2) - 252*G4*G4ppx*G5x*pow(H,2) + 168*G4p*G4px*G5x*pow(H,2) + 48*G4pp*G4x*G5x*pow(H,2) - 18*G4pp*G5p*G5x*pow(H,2) - 
      6*G4p*G5pp*G5x*pow(H,2) + 12*G4*G5ppp*G5x*pow(H,2) - 9*G2p*pow(G5x,2)*pow(H,2) + 24*pow(G4p,2)*G5xx*pow(H,2) + 12*G4*G4pp*G5xx*pow(H,2) + 90*G4*G5px*G5x*pow(H,4) - 
      108*G4p*pow(G5x,2)*pow(H,4) - 24*G4pxx*G4x*rptot + 24*G4px*G4xx*rptot + 12*G4pxx*G5p*rptot - 12*G4xx*G5pp*rptot + 12*G4x*G5ppx*rptot - 6*G5p*G5ppx*rptot - 18*G4px*G5px*rptot + 
      6*G5pp*G5px*rptot + 6*G4ppx*G5x*rptot + 3*G3x*(16*G4p*G4px + 4*G4pp*G4x - 2*G4pp*G5p + 4*G4p*G5pp - G2p*G5x - 30*G4p*G5x*pow(H,2) + G4*(-20*G4ppx + 4*G5ppp + 14*G5px*pow(H,2)) + G5px*rptot) + 
      3*G3px*(6*G3x*G4 - 20*G4*G4px + 4*G4p*G4x - 2*G4p*G5p + 4*G4*G5pp + 30*G4*G5x*pow(H,2) - G5x*rptot))))/
(a*pow(H,2)*pow(-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
  4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
  2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
     18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
  2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
     14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
     3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
   (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
     2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*
   (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
     18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4)),2));

   // Milestone spih
   *spih = (4*pow(a,10)*G4*G4p*H + ddphi*pow(dphi,8)*(-pow(G5xx,2) + G5x*G5xxx)*pow(H,3) + 
   a*(ddphi*pow(dphi,7)*(4*G4xxx*G5x - 2*G5pxx*G5x - 6*G4xx*G5xx + 3*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,2) + pow(dphi,9)*(pow(G5xx,2) - G5x*G5xxx)*pow(H,4)) - 
   2*pow(a,9)*dphi*(G2x*G4 + 2*(-(G3p*G4) + pow(G4p,2) + G4*G4x*pow(H,2) - G4*G5p*pow(H,2))) - 
   2*pow(a,8)*(4*dH*dphi*G4*(G4x - G5p) + 4*ddphi*G4*(G4x - G5p)*H + pow(dphi,2)*H*(G3x*G4 + G4p*G5p - 2*G4*G5pp - 3*G4*G5x*pow(H,2))) - 
   pow(a,7)*dphi*(12*dH*dphi*G4*G5x*H + 4*ddphi*(G3x*G4 - 3*G4*G4px - G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) + 
      pow(dphi,2)*(2*G3px*G4 - 2*G3x*G4p - 4*G4*G4ppx - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 2*G4p*G5pp - 12*pow(G4x,2)*pow(H,2) - 20*G4*G4xx*pow(H,2) + 22*G4x*G5p*pow(H,2) - 
         10*pow(G5p,2)*pow(H,2) + 18*G4*G5px*pow(H,2) + 2*G4p*G5x*pow(H,2))) + 
   pow(a,2)*pow(dphi,6)*H*(dphi*H*(dH*G5x*G5xx + dphi*(-4*G4xxx*G5x + 3*G5pxx*G5x + 6*G4xx*G5xx - 4*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*H) + 
      ddphi*(-8*pow(G4xx,2) + 8*G4x*G4xxx - 4*G4xxx*G5p + 8*G4xx*G5px - 2*pow(G5px,2) - 4*G4x*G5pxx + 2*G5p*G5pxx + G3xx*G5x - 2*G4pxx*G5x - G3x*G5xx + 2*G4px*G5xx + G5x*G5xx*pow(H,2))) + 
   pow(a,4)*pow(dphi,4)*(ddphi*H*(-8*G4*G4xxx + 4*G4x*(4*G4xx - 3*G5px) + 2*G5p*G5px + 4*G4*G5pxx - G3x*G5x + 2*G4p*G5xx - 3*pow(G5x,2)*pow(H,2)) + 
      dH*dphi*(8*G4x*G4xx - 4*G4xx*G5p - 4*G4x*G5px + 2*G5p*G5px - G3x*G5x + 2*G4px*G5x + 3*pow(G5x,2)*pow(H,2)) + 
      pow(dphi,2)*H*(12*G4pxx*G4x + 2*G3x*G4xx - 12*G4px*G4xx - 6*G4pxx*G5p + G3xx*(-2*G4x + G5p) + 4*G4xx*G5pp - 4*G4x*G5ppx + 2*G5p*G5ppx - 2*G3x*G5px + 8*G4px*G5px - 2*G5pp*G5px + G3px*G5x - 
         2*G4ppx*G5x + 8*G4xx*G5x*pow(H,2) - 3*G5px*G5x*pow(H,2) - 6*G4x*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2))) - 
   pow(a,3)*pow(dphi,5)*(dphi*H*(2*dH*(-2*G4x + G5p)*G5xx - dphi*(8*pow(G4xx,2) - 8*G4x*G4xxx + 4*G4xxx*G5p - 12*G4xx*G5px + 4*pow(G5px,2) + 6*G4x*G5pxx - 3*G5p*G5pxx - G3xx*G5x + 6*G4pxx*G5x - 
            2*G5ppx*G5x + G3x*G5xx - 4*G4px*G5xx + G5pp*G5xx)*H) + ddphi*(4*G4pxx*G4x + 2*G3x*G4xx - 4*G4px*G4xx - 2*G4pxx*G5p + G3xx*(-2*G4x + G5p) - G3x*G5px + 2*G4px*G5px + 2*G4xx*G5x*pow(H,2) + 
         G5px*G5x*pow(H,2) - 8*G4x*G5xx*pow(H,2) + G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2))) + 
   pow(a,6)*pow(dphi,2)*(2*dH*dphi*(4*pow(G4x,2) - 4*G4*G4xx - 6*G4x*G5p + 2*pow(G5p,2) + 2*G4*G5px + G4p*G5x) + 2*ddphi*(-16*G4*G4xx + 2*G4x*G5p - 2*pow(G5p,2) + 10*G4*G5px + 3*G4p*G5x)*H + 
      pow(dphi,2)*H*(2*G3xx*G4 + 4*G3x*G4x - 12*G4px*G4x + 4*G4p*G4xx - 3*G3x*G5p + 12*G4px*G5p - 2*G5p*G5pp + G2x*G5x - 2*G3p*G5x + 14*G4x*G5x*pow(H,2) - 17*G5p*G5x*pow(H,2) + 
         4*G4*(-3*G4pxx + G5ppx + 3*G5xx*pow(H,2)))) + pow(a,5)*pow(dphi,3)*
    (-2*ddphi*(G3xx*G4 - 2*G4*G4pxx - G3x*G4x + 4*G4px*G4x - 2*G4p*G4xx - G4px*G5p + G4p*G5px + G4x*G5x*pow(H,2) - 4*G5p*G5x*pow(H,2) + 7*G4*G5xx*pow(H,2)) + 
      dphi*(2*dH*(6*G4x*G5x - 3*G5p*G5x - 2*G4*G5xx)*H + dphi*(4*pow(G4px,2) + 2*G3px*G4x - 4*G4ppx*G4x - G3px*G5p + 2*G4ppx*G5p - 4*G4x*G4xx*pow(H,2) + 8*G4*G4xxx*pow(H,2) - 6*G4xx*G5p*pow(H,2) + 
            6*G4x*G5px*pow(H,2) + 3*G5p*G5px*pow(H,2) - 6*G4*G5pxx*pow(H,2) + G5pp*G5x*pow(H,2) + 6*pow(G5x,2)*pow(H,4) + G3x*(-2*G4px + G5pp + 4*G5x*pow(H,2)) - 
            2*G4px*(G5pp + 6*G5x*pow(H,2))))))/
 (a*(-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
     4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
     2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
        18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
     2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
        14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
        3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
      (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
        2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*
      (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
        18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))));

   // Milestone spie
   *spie = (-2*a*H*(4*pow(a,8)*G4*G4p + 8*pow(a,7)*dphi*G4*(G4x - G5p)*H + ddphi*pow(dphi,6)*G5x*G5xx*pow(H,2) - 
   2*pow(a,5)*dphi*(2*dH*dphi*G4*G5x + 4*ddphi*G4*G5x*H + pow(dphi,2)*(-8*G4*G4xx + 2*G4x*G5p - 2*pow(G5p,2) + 6*G4*G5px + 3*G4p*G5x)*H) - 
   a*pow(dphi,5)*G5xx*H*(-4*ddphi*G4x + 2*ddphi*G5p + pow(dphi,2)*G5x*pow(H,2)) + 
   2*pow(a,6)*(4*ddphi*G4*(-G4x + G5p) + pow(dphi,2)*(G3x*G4 - 6*G4*G4px - 4*G4p*G4x + 3*G4p*G5p + 2*G4*G5pp + 5*G4*G5x*pow(H,2))) + 
   pow(a,2)*pow(dphi,4)*(dphi*H*(2*dH*pow(G5x,2) + dphi*(2*G5px*G5x - 4*G4x*G5xx + G5p*G5xx)*H) + 
      ddphi*(8*G4x*G4xx - 4*G4xx*G5p - 4*G4x*G5px + 2*G5p*G5px - G3x*G5x + 2*G4px*G5x + pow(G5x,2)*pow(H,2))) + 
   pow(a,4)*(2*ddphi*pow(dphi,2)*(4*pow(G4x,2) - 4*G4*G4xx - 6*G4x*G5p + 2*pow(G5p,2) + 2*G4*G5px + G4p*G5x) + 
      pow(dphi,4)*(8*G4px*G4x - G3x*G5p - 2*G4px*G5p - 4*G4x*G5pp + 2*G5p*G5pp - 5*G5p*G5x*pow(H,2) + 6*G4*G5xx*pow(H,2))) + 
   pow(a,3)*pow(dphi,3)*(2*dH*dphi*(2*G4x - G5p)*G5x + 4*ddphi*(2*G4x*G5x - G5p*G5x - G4*G5xx)*H + pow(dphi,2)*H*(-8*G4x*(G4xx - G5px) - 2*G5p*G5px + G5x*(G3x + 2*G4px - 2*G5pp + G5x*pow(H,2))))))/
(-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
 4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
 2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
    18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
 2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
    14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
    3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
  (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 2*G4x*G5xx*pow(H,2) + 
    5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*(3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
    24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4)));

   // Milestone spim
   *spim = (3*pow(a,6)*(-2*pow(a,4)*G4p + 4*pow(a,3)*dphi*(G4x - G5p)*H + 2*a*pow(dphi,3)*(2*G4xx - G5px)*H + pow(dphi,4)*G5xx*pow(H,2) + pow(a,2)*pow(dphi,2)*(G3x - 2*G4px + 3*G5x*pow(H,2))))/
   (3*pow(-2*pow(a,4)*G4p + 4*pow(a,3)*dphi*(G4x - G5p)*H + 2*a*pow(dphi,3)*(2*G4xx - G5px)*H + pow(dphi,4)*G5xx*pow(H,2) + pow(a,2)*pow(dphi,2)*(G3x - 2*G4px + 3*G5x*pow(H,2)),2) + 
     2*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H)*(3*a*pow(dphi,4)*(2*G4xxx - G5pxx)*pow(H,2) + pow(dphi,5)*G5xxx*pow(H,3) + 
        pow(a,3)*pow(dphi,2)*(G2xx - G3px + 3*(8*G4xx - 5*G5px)*pow(H,2)) + 6*pow(a,4)*dphi*H*(G3x - 3*G4px + G5x*pow(H,2)) + pow(a,2)*pow(dphi,3)*H*(3*G3xx - 6*G4pxx + 7*G5xx*pow(H,2)) + 
        pow(a,5)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))));

   // Milestone s00
   *s00 = (-2*a*(6*pow(a,4)*dphi*G4pp*H + 3*a*pow(dphi,4)*(-2*G4pxx + G5ppx)*pow(H,2) - pow(dphi,5)*G5pxx*pow(H,3) + pow(a,5)*(G2p + 6*G4p*pow(H,2)) + 
   pow(a,3)*pow(dphi,2)*(-G2px + G3pp + 3*(-4*G4px + 3*G5pp)*pow(H,2)) + pow(a,2)*pow(dphi,3)*H*(-3*G3px + 6*G4ppx - 5*G5px*pow(H,2))))/
(2*pow(a,4)*dphi*G4p + 4*pow(a,5)*G4*H + 2*pow(a,3)*pow(dphi,2)*(-4*G4x + 3*G5p)*H + 2*a*pow(dphi,4)*(-2*G4xx + G5px)*H - pow(dphi,5)*G5xx*pow(H,2) - 
 pow(a,2)*pow(dphi,3)*(G3x - 2*G4px + 5*G5x*pow(H,2)));   

   // Milestone s00k
   *s00k = (-4*pow(a,4)*G4p + 8*pow(a,3)*dphi*(G4x - G5p)*H + 4*a*pow(dphi,3)*(2*G4xx - G5px)*H + 2*pow(dphi,4)*G5xx*pow(H,2) + 2*pow(a,2)*pow(dphi,2)*(G3x - 2*G4px + 3*G5x*pow(H,2)))/
   (2*pow(a,4)*dphi*G4p + 4*pow(a,5)*G4*H + 2*pow(a,3)*pow(dphi,2)*(-4*G4x + 3*G5p)*H + 2*a*pow(dphi,4)*(-2*G4xx + G5px)*H - pow(dphi,5)*G5xx*pow(H,2) - 
     pow(a,2)*pow(dphi,3)*(G3x - 2*G4px + 5*G5x*pow(H,2)));   

   // Milestone s00p
   *s00p = (-2*(6*pow(a,6)*G4p*H + 3*a*pow(dphi,5)*(-2*G4xxx + G5pxx)*pow(H,2) - pow(dphi,6)*G5xxx*pow(H,3) + pow(a,3)*pow(dphi,3)*(-G2xx + G3px + 3*(-12*G4xx + 7*G5px)*pow(H,2)) - 
   3*pow(a,4)*pow(dphi,2)*H*(3*G3x - 8*G4px + 5*G5x*pow(H,2)) + pow(a,2)*pow(dphi,4)*H*(-3*G3xx + 6*G4pxx - 10*G5xx*pow(H,2)) - pow(a,5)*dphi*(G2x - 2*(G3p + 9*(-G4x + G5p)*pow(H,2)))))/
(a*(2*pow(a,4)*dphi*G4p + 4*pow(a,5)*G4*H + 2*pow(a,3)*pow(dphi,2)*(-4*G4x + 3*G5p)*H + 2*a*pow(dphi,4)*(-2*G4xx + G5px)*H - pow(dphi,5)*G5xx*pow(H,2) - 
   pow(a,2)*pow(dphi,3)*(G3x - 2*G4px + 5*G5x*pow(H,2))));

   // Milestone s0i
   *s0i = (-2*pow(a,4)*G4p*H + 2*a*pow(dphi,3)*(3*G4xx - 2*G5px)*pow(H,2) + pow(dphi,4)*G5xx*pow(H,3) + pow(a,3)*dphi*(G2x - 2*G3p + 2*G4pp + 6*G4x*pow(H,2) - 6*G5p*pow(H,2)) + 
   pow(a,2)*pow(dphi,2)*H*(3*G3x - 10*G4px + 2*G5pp + 3*G5x*pow(H,2)))/(4*pow(a,3)*G4 + 2*a*pow(dphi,2)*(-2*G4x + G5p) - 2*pow(dphi,3)*G5x*H);

   // Milestone s0ip
   *s0ip = (2*pow(a,4)*G4p + 4*pow(a,3)*dphi*(-G4x + G5p)*H + 2*a*pow(dphi,3)*(-2*G4xx + G5px)*H - pow(dphi,4)*G5xx*pow(H,2) - pow(a,2)*pow(dphi,2)*(G3x - 2*G4px + 3*G5x*pow(H,2)))/
   (2.*a*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H));

   // Milestone sii
   *sii = (-3*(6*pow(a,8)*G4p*pow(H,2) + pow(ddphi,2)*pow(dphi,4)*G5xxx*pow(H,2) + 2*a*ddphi*pow(dphi,3)*H*(2*ddphi*G4xxx - ddphi*G5pxx - pow(dphi,2)*G5xxx*pow(H,2)) + 
   pow(a,7)*(-2*dH*G4p + 2*dphi*H*(G4pp + 4*(-G4x + G5p)*pow(H,2))) + pow(a,2)*pow(dphi,2)*
    (dddphi*dphi*G5xx*pow(H,2) + pow(dphi,4)*G5xxx*pow(H,4) + ddphi*dphi*H*(4*dH*G5xx + dphi*(-8*G4xxx + 5*G5pxx)*H) + pow(ddphi,2)*(G3xx - 2*G4pxx + 6*G5xx*pow(H,2))) - 
   pow(a,6)*(2*ddphi*G4pp + dphi*H*(-6*dH*G5p + dphi*H*(G3x + 4*G4px - 4*G5pp + 3*G5x*pow(H,2)))) - 
   pow(a,5)*(2*ddH*dphi*(-2*G4x + G5p) + H*(4*dddphi*(-G4x + G5p) + pow(dphi,3)*(G3px - 2*G4ppx + (-4*G4xx + 7*G5px)*pow(H,2)) + 2*ddphi*dphi*(G3x - 6*G4px + 2*G5pp + 3*G5x*pow(H,2))) + 
      dH*(-4*ddphi*G4x + 6*ddphi*G5p + pow(dphi,2)*(G3x - 6*G4px + 2*G5pp + 9*G5x*pow(H,2)))) + 
   pow(a,3)*dphi*(4*pow(ddphi,2)*(3*G4xx - 2*G5px)*H + dphi*H*(dddphi*(4*G4xx - 2*G5px) + pow(dphi,2)*H*(-5*dH*G5xx + 4*dphi*G4xxx*H - 3*dphi*G5pxx*H)) - 
      2*ddphi*dphi*(dH*(-4*G4xx + 2*G5px) + dphi*H*(G3xx - 4*G4pxx + G5ppx + 6*G5xx*pow(H,2)))) + 
   pow(a,4)*(pow(ddphi,2)*(G3x - 4*G4px + 3*G5x*pow(H,2)) + ddphi*dphi*(10*dH*G5x*H + dphi*(G3px - 2*G4ppx + (-24*G4xx + 19*G5px)*pow(H,2))) + 
      dphi*(dddphi*(G3x - 2*G4px + 3*G5x*pow(H,2)) + dphi*(2*pow(dH,2)*G5x + 4*dH*dphi*(-3*G4xx + 2*G5px)*H + H*(2*ddH*G5x + pow(dphi,2)*H*(G3xx - 6*G4pxx + 2*G5ppx + 4*G5xx*pow(H,2))))))) + 
(12*pow(a,14)*(48*pow(dphi,18)*G4*G5xx*G5xxx*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,12)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
     4*a*pow(dphi,17)*pow(H,11)*(-3*G4*G5xx*G5xxx*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
        2*G4*G5xx*G5xxx*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
        12*G4*G5xx*G5xxx*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
        24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) - 
        4*pow(-3*pow(G5xx,2) + G5x*G5xxx,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))) + 
     576*pow(a,18)*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*rptot - 
     2*pow(a,18)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
      (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*rptot*
      (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 32*pow(a,18)*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
      pow(6*pow(G4p,2) - 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) + 3*G4x*rptot - 3*G5p*rptot,2) + 
     96*pow(a,18)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot*
      (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
     32*pow(a,18)*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
      (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
     4*pow(a,18)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot*
      (4*G2px*G4 - 4*G3pp*G4 + 6*G2x*G4p - 12*G3p*G4p + 12*G4p*G4pp + 4*G2p*G4x - 2*G2p*G5p + 12*G3x*G4*pow(H,2) + 156*G4p*G4x*pow(H,2) - 132*G4p*G5p*pow(H,2) - 12*G4*G5pp*pow(H,2) + 
        12*G4*G5x*pow(H,4) - 3*G3x*rptot + 6*G4px*rptot - 9*G5x*pow(H,2)*rptot) + 
     96*pow(a,18)*G4*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
      (4*G2x*G4p - 8*G3p*G4p + 2*G2p*G4x - 2*G2p*G5p + 48*G4p*G4x*pow(H,2) - 48*G4p*G5p*pow(H,2) - 3*G3x*rptot + 9*G4px*rptot - 3*G5x*pow(H,2)*rptot) + 
     96*pow(a,18)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot)*
      (5*G2x*G4*G4p - 10*G3p*G4*G4p + 15*pow(G4p,3) + 30*G4*G4p*G4x*pow(H,2) - 30*G4*G4p*G5p*pow(H,2) - 3*G3x*G4*rptot + 9*G4*G4px*rptot + 6*G4p*G4x*rptot - 6*G4p*G5p*rptot - 
        3*G4*G5x*pow(H,2)*rptot) - 288*pow(a,18)*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
      (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
     96*pow(a,18)*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
      (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
     48*pow(a,18)*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
      (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
     32*pow(a,18)*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
      (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
      (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
     4*pow(a,18)*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
      (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
        3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 4*pow(a,18)*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
      (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot)*
      (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
        3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 4*pow(a,18)*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
      (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
        2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
        3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
     32*pow(a,18)*pow(H,2)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
      (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
           G4*G5x*pow(H,2)*rptot)) + 4*pow(a,2)*pow(dphi,16)*pow(H,10)*
      (2*G4*G5xx*G5xxx*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        6*G4*G5xx*G5xxx*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        6*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
        4*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        2*G4p*pow(-3*pow(G5xx,2) + G5x*G5xxx,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
        24*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*
         (G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 12*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) + 
        12*G4*G5xx*G5xxx*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*pow(G5xx,2) - G5x*G5xxx)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        (3*pow(G5xx,2) - 2*G5x*G5xxx)*(3*pow(G5xx,2) - G5x*G5xxx)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,3)*pow(dphi,15)*pow(H,9)*(48*G4*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 
        16*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        16*G4p*(6*G4xxx*G5x - 4*G5pxx*G5x - 27*G4xx*G5xx + 15*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*(-3*pow(G5xx,2) + G5x*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        24*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) + 
        48*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) - 12*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) + 
        8*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) + 
        48*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) + 
        24*G4*G5xx*G5xxx*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) + 
        8*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 96*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 16*pow(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        32*(3*pow(G5xx,2) - G5x*G5xxx)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) + 3*pow(3*pow(G5xx,2) - 2*G5x*G5xxx,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot - 
        2*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot - 
        4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*(-3*pow(G5xx,2) + G5x*G5xxx)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot))) + 2*pow(a,5)*pow(dphi,13)*pow(H,7)*
      (12*G4*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) + 4*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) - 12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) - 
        6*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2))) + 4*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        48*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        4*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 
        4*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) + 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 
        24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2))) - 6*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        4*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        24*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 
           3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 16*(3*pow(G5xx,2) - G5x*G5xxx)*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        8*pow(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)),2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 2*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot + 
        6*pow(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot + 
        6*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot - 
        2*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot - 
        (3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot + 12*G4*G5xx*G5xxx*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 24*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - (3*pow(G5xx,2) - 2*G5x*G5xxx)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 
           12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 
           2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(3*pow(G5xx,2) - G5x*G5xxx)*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     2*pow(a,4)*pow(dphi,14)*pow(H,8)*(-4*G4p*pow(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx,2)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) - 
        24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 
        6*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) + 4*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) + 24*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 
           6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 
           96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) + 
        4*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) - 
        12*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) + 
        24*G4*G5xx*G5xxx*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        2*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) + 
        24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 
        12*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2))) + 24*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        8*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 8*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 8*(3*pow(G5xx,2) - G5x*G5xxx)*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        16*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - (-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot + 
        6*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot - 
        2*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*(-3*pow(G5xx,2) + G5x*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot - 
        4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G5xx,2) - G5x*G5xxx)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 
           6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        2*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,6)*pow(dphi,12)*pow(H,6)*(-24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) - 24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2))) + 16*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) - 
        16*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        48*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) + 
        24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 
           3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 
           2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 
        24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2))) + 8*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        24*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 
           3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 
        48*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 8*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)),2) - 24*G4*G5xx*G5xxx*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))) - 
        32*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        16*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot + 
        24*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot - 
        (3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot + 
        12*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot + 48*G4*G5xx*G5xxx*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        48*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        4*G4*G5xx*G5xxx*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 48*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 16*(3*pow(G5xx,2) - G5x*G5xxx)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 12*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        8*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        48*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        4*(3*pow(G5xx,2) - G5x*G5xxx)*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + 
           G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        8*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     2*pow(a,17)*dphi*H*(-24*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (-2*G2xx*G4 + 2*G3px*G4 + 6*G3x*G4p - 12*G4p*G4px + 2*G2x*G4x - 4*G3p*G4x - G2x*G5p + 2*G3p*G5p - 12*pow(G4x,2)*pow(H,2) - 48*G4*G4xx*pow(H,2) + 30*G4x*G5p*pow(H,2) - 
           18*pow(G5p,2)*pow(H,2) + 30*G4*G5px*pow(H,2) + 18*G4p*G5x*pow(H,2))*rptot + 
        864*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*pow(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)),2)*rptot - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 32*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*pow(6*pow(G4p,2) - 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) + 3*G4x*rptot - 3*G5p*rptot,2) + 
        4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        4*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px + 10*G2x*G4x - 20*G3p*G4x + 12*G4pp*G4x - 8*G2x*G5p + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           84*pow(G4x,2)*pow(H,2) + 12*G4*G4xx*pow(H,2) - 156*G4x*G5p*pow(H,2) + 72*pow(G5p,2)*pow(H,2) - 16*G4*G5px*pow(H,2) - 48*G4p*G5x*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        288*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        96*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        12*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 4*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 
        8*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        2*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) + 
        96*pow(H,2)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,16)*pow(dphi,2)*(48*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 14*G4*G5xx*pow(H,2))*rptot + 
        288*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot - (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        8*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        48*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 
           6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        32*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        48*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        16*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 32*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 2*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot,2) - 288*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        96*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        48*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        16*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*pow(H,2)*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) + 
        16*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 16*pow(H,2)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 48*pow(H,2)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,8)*pow(dphi,10)*pow(H,4)*(-48*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2)) - 
        288*G4*G5xx*G5xxx*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2))) - 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2))) + 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 
           6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) - 2*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)),2) + 
        48*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 
           12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        16*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))) - 
        16*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        2*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot + 
        12*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot + 
        24*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 
           3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot + 
        12*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        96*G4*G5xx*G5xxx*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        16*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        32*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        4*G4*G5xx*G5xxx*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 48*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 16*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 48*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        8*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        24*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        48*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        12*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        48*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 16*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - (3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        8*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 
           6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     2*pow(a,7)*pow(dphi,11)*pow(H,5)*(-12*G4*G5xx*G5xxx*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2)) - 
        6*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2))) - 8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 
           6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2))) + 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 
           3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 
           2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2))) - 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2))) - 24*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))) - 
        2*pow(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)),2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) - 
        16*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) + 6*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2),2)*
         rptot + 3*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        (3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        (3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot - 
        (-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot - 
        2*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot + 12*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot + 4*G4*G5xx*G5xxx*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        48*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        12*G4*G5xx*G5xxx*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 6*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        24*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        24*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        4*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 4*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        12*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        2*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 2*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 
           6*G3p*G5xx + 6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 
           2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 
           3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - (3*pow(G5xx,2) - 2*G5x*G5xxx)*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 
              2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 2*pow(a,13)*pow(dphi,5)*H*
      (-24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        24*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot - 
        4*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx - 3*G3xx*G5x + 12*G4pxx*G5x - 3*G5ppx*G5x + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           14*G5x*G5xx*pow(H,2))*rptot + 6*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot - 
        12*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot + 
        144*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot + 6*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(G2x*G5x - 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)),2)*rptot + 
        12*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        4*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        4*pow(H,2)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) + 
        2*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) + 8*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        32*pow(H,2)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 
           3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*
            pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        8*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         pow(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot,2) - 
        12*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        12*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 
        6*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        12*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 144*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 
           4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        48*pow(H,2)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 2*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 8*pow(H,2)*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + 
              (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 
           12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        2*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     2*pow(a,9)*pow(dphi,9)*pow(H,3)*(-24*G4*G5xx*G5xxx*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2) - 
        12*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2)) - 
        288*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2))) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2))) - 
        6*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 
           6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        12*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + 
              (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))) - 8*pow(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 
           4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)),2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2)))) + 
        6*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot - 2*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot + 6*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot + 6*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)),2)*rptot - (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot - 
        2*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot + 
        12*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot + 
        24*G4*G5xx*G5xxx*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) - 
        8*G4*G5xx*G5xxx*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        96*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        32*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 
           6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*pow(H,2)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 
           12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        (3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        48*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        4*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 4*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 12*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 8*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 6*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        4*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        24*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        12*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 
           4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) + 12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) - 4*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 
           6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - (-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 
           4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 48*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 
           6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 
           3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - (3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        2*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     2*pow(a,15)*pow(dphi,3)*H*(6*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2),2)*rptot + 
        12*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
           18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))*rptot - 
        4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot + 
        144*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        32*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        6*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 2*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         pow(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot,2) - 12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        144*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 
           4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        12*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 8*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 12*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 2*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 
           4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 8*pow(H,2)*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + 
              (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 48*pow(H,2)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 2*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 
              2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + pow(a,14)*pow(dphi,4)*
      (48*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (-6*G3xx*G4x + 12*G4pxx*G4x + 12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - G2xx*G5x + G3px*G5x + 12*G4xx*G5x*pow(H,2) - 3*G5px*G5x*pow(H,2) - 
           2*G4x*G5xx*pow(H,2) - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2))*rptot - 
        4*pow(H,2)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (48*G4pxx*G4x + 78*G3x*G4xx - 204*G4px*G4xx - 24*G4pxx*G5p + 6*G3xx*(-2*G4x + G5p) + 24*G4xx*G5pp - 12*G4x*G5ppx + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px - 2*G2xx*G5x + 
           8*G3px*G5x - 12*G4ppx*G5x + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx + 150*G4xx*G5x*pow(H,2) - 74*G5px*G5x*pow(H,2) + 38*G4x*G5xx*pow(H,2) - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2))*
         rptot + 144*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot - 
        48*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot + 
        24*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        4*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        16*pow(H,2)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 
           3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*
            pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        2*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 32*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        8*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        8*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot,2) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 
           3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 16*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 16*pow(H,2)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 
           12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        288*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        96*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 
           4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 32*pow(H,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 48*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        16*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        16*pow(H,2)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 8*pow(H,2)*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + 
              (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 48*pow(H,2)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,11)*pow(dphi,7)*H*(-48*G4*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + 10*pow(G5xx,2)*pow(H,2) + 3*G5x*G5xxx*pow(H,2)) - 
        288*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2))) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) + 
        24*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot - 
        8*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot + 
        24*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot - 
        48*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot + 288*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot + 
        3*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)),2)*rptot - 
        2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + 
           G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*rptot - 4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot - 2*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot + 
        24*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + 
              (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot + 8*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 
        8*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        24*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 16*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        32*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        96*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        16*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        64*pow(H,4)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        8*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        32*pow(H,2)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 
           3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*
            pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        12*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        8*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        48*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        16*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        48*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        96*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)
          + 4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 8*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 
        8*pow(H,2)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 2*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) + 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 
           4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot)) - 24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - 
           G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*pow(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot),2) - 24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        16*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 96*pow(H,4)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        8*pow(H,2)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 
           3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*
            pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 16*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*pow(H,2)*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 2*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + 
     pow(a,10)*pow(dphi,8)*pow(H,2)*(-96*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*
         pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) - 288*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2))) - 8*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)),2) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))) - 
        4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot + 24*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         rptot + 144*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         rptot - 48*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         rptot - (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 
           12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 
           2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot + 12*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 
           2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 
        4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot + 24*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-(G2x*G5x) + 
           2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        16*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 8*G4*G5xx*G5xxx*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) + 96*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) - 32*G4*(3*G4xxx*G5xx - G5pxx*G5xx + 3*G4xx*G5xxx - 2*G5px*G5xxx)*pow(H,4)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        32*G4p*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        96*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        64*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 
           6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 
           12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        16*pow(H,2)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        48*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        24*G4*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        16*G4*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        96*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 
           4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 16*pow(H,2)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - (3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 8*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        24*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        48*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        24*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 
           4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        4*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) + 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 16*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + 
           G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 2*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 12*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 
           12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        8*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        24*G4*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 
           4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        16*(3*pow(G5xx,2) - G5x*G5xxx)*pow(H,4)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 96*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*pow(H,2)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 8*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 8*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 
           6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,2)*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 
           2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + 
           pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 
              3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 
              2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot))) + pow(a,12)*pow(dphi,6)*
      (-48*G4*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2)) - 288*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 
           12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2))) - 
        48*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot + 
        8*(6*G4xxx*G5x - 4*G5pxx*G5x - 27*G4xx*G5xx + 15*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,4)*pow(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2),2)*rptot + 
        288*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*rptot - 2*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*rptot + 
        24*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot - 48*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*rptot + 
        12*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(-(G2x*G5x) + 
           2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*rptot - 
        8*G4p*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 4*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (36*G4xx*G4xxx - 24*G4xxx*G5px - 14*G4xx*G5pxx + 10*G5px*G5pxx + 3*G3xx*G5xx - 3*G5ppx*G5xx + 3*G3x*G5xxx - 10*G4px*G5xxx + 2*G5pp*G5xxx + (10*pow(G5xx,2) + 3*G5x*G5xxx)*pow(H,2))*
         (4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) - 16*pow(H,4)*(60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 
           24*G4px*G5xx + 3*G5pp*G5xx + G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*(4*G4*(G2p + 12*G4p*pow(H,2)) - 6*G4p*rptot) + 
        48*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) - 
        2*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot) - 
        16*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (18*G3x*G4xxx - 60*G4px*G4xxx + 12*G4xxx*G5pp - 12*G4xx*G5ppx + 6*G3xx*(3*G4xx - 2*G5px) - 12*G4pxx*(G4xx - G5px) + 6*G5ppx*G5px - 8*G3x*G5pxx + 28*G4px*G5pxx - 6*G5pp*G5pxx + G2xx*G5xx + 
           2*G3px*G5xx - 6*G4ppx*G5xx + G2x*G5xxx - 2*G3p*G5xxx + 2*G4pp*G5xxx + 18*G4xxx*G5x*pow(H,2) - 6*G5pxx*G5x*pow(H,2) + 96*G4xx*G5xx*pow(H,2) - 56*G5px*G5xx*pow(H,2) + 
           6*G4x*G5xxx*pow(H,2) - 6*G5p*G5xxx*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        96*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        32*G4p*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        32*pow(H,4)*(78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) + 
        8*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*rptot*
         (-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot) - 
        16*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) - 
        4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*rptot*(2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 
           12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot) + 
        8*G4*(-6*G4xxx*G5x + 4*G5pxx*G5x + 27*G4xx*G5xx - 15*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        24*G4*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) + 
        48*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(2*G4p*(G2p + 12*G4p*pow(H,2)) + (G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))*rptot) - 
        48*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        48*G4*pow(H,4)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot) + 
        16*G4*pow(H,4)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-4*G2x*G4p + 8*G3p*G4p + 2*G2p*G5p + 48*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 24*G4p*pow(H,2)) + 3*(G3x - 3*G4px + G5x*pow(H,2))*rptot)
          - 4*G4p*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-30*G3xx*G4px + 48*G4px*G4pxx + 6*G2xx*G4xx + 6*G3px*G4xx - 24*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx + 6*G3xx*G5pp - 12*G4pxx*G5pp + 6*G4px*G5ppx - 4*G2xx*G5px - 2*G3px*G5px + 
           12*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx + 
           (216*pow(G4xx,2) + 36*G4x*G4xxx - 36*G4xxx*G5p - 250*G4xx*G5px + 74*pow(G5px,2) - 14*G4x*G5pxx + 14*G5p*G5pxx + 9*G3xx*G5x - 9*G5ppx*G5x - 112*G4px*G5xx + 11*G5pp*G5xx - 2*G4p*G5xxx)*
            pow(H,2) + 45*G5x*G5xx*pow(H,4) + 3*G3x*(3*G3xx - 4*G4pxx - G5ppx + 13*G5xx*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 16*pow(H,2)*(6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 
           4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 
           2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 30*pow(G5x,2)*pow(H,4) + 
           3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*(3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) - 2*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*rptot*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot) + 48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        4*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 
           2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) - 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 
           2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*(pow(G2x,2) + 4*pow(G3p,2) - G2p*G3x - 2*G2px*G4p + 2*G3pp*G4p - 4*G3p*G4pp + 2*G2p*G4px + 
           6*(18*pow(G4x,2) - 36*G4x*G5p + 18*pow(G5p,2) - 11*G4p*G5x)*pow(H,4) + 2*G2x*(-2*G3p + G4pp + 12*(G4x - G5p)*pow(H,2)) + G2xx*rptot - G3px*rptot - 
           3*pow(H,2)*(14*G3x*G4p + 16*G3p*G4x - 4*G4pp*G4x - 16*G3p*G5p + 4*G4pp*G5p - 2*G4p*(16*G4px + G5pp) + G2p*G5x - 8*G4xx*rptot + 5*G5px*rptot)) + 
        24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) - 
        48*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot)) + 
        12*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*(24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 
           4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot)) - 
        (G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 
           15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*rptot*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 16*pow(H,2)*
         (3*G3x*G4*G4p + 9*pow(G4p,2)*(-G4x + G5p) + G4*G4p*(-9*G4px + 3*G5x*pow(H,2)) - G4*(G4x - G5p)*(G2x - 2*(G3p + 3*(-G4x + G5p)*pow(H,2))))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 4*G4*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (-(G2x*G2xx) + 2*G2xx*G3p + G2x*G3px - 2*G3p*G3px - G2px*G3x + G3pp*G3x - 2*G2xx*G4pp + 2*G3px*G4pp + 2*G2px*G4px - 2*G3pp*G4px + 
           (-324*G4x*G4xx + 324*G4xx*G5p + 178*G4x*G5px - 178*G5p*G5px - 72*G3x*G5x + 186*G4px*G5x - 3*G5pp*G5x + 32*G4p*G5xx)*pow(H,4) - 45*pow(G5x,2)*pow(H,6) + 
           pow(H,2)*(-27*pow(G3x,2) + 6*G3xx*G4p - 216*pow(G4px,2) - 6*(G2xx + G3px - 4*G4ppx)*G4x - 42*G2x*G4xx + 84*G3p*G4xx - 48*G4pp*G4xx + 6*G2xx*G5p + 6*G3px*G5p - 24*G4ppx*G5p + 
              3*G3x*(50*G4px - 3*G5pp) + 30*G4px*G5pp - 6*G4p*G5ppx + 25*G2x*G5px - 50*G3p*G5px + 30*G4pp*G5px - 3*G2px*G5x + 3*G3pp*G5x + G2p*G5xx - 6*G4xxx*rptot + 3*G5pxx*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (24*G3p*G3x + 2*G2xx*G4p + 4*G3px*G4p - 12*G3x*G4pp - 12*G4p*G4ppx - 68*G3p*G4px + 36*G4pp*G4px - 4*G2px*G4x + 4*G3pp*G4x + 4*G2p*G4xx + 4*G2px*G5p - 4*G3pp*G5p + 4*G3p*G5pp - 2*G2p*G5px - 
           144*(G4x - G5p)*G5x*pow(H,4) - 2*G2x*(6*G3x - 17*G4px + G5pp + 9*G5x*pow(H,2)) - 3*G3xx*rptot + 6*G4pxx*rptot + 
           pow(H,2)*(-108*G3x*(G4x - G5p) + 4*(33*G4p*G4xx + 69*G4px*(G4x - G5p) - 17*G4p*G5px + 9*G3p*G5x - 3*G4pp*G5x) - 7*G5xx*rptot))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot)) - 2*G4p*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         pow(2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot),2) - 24*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        8*G4*pow(H,2)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 
           2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + 
           pow(H,2)*(-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) + 
        24*G4*pow(H,2)*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(3*G2x*G3xx - 6*G3p*G3xx + 6*G3xx*G4pp - 10*G2xx*G4px + 4*G3px*G4px + 12*G4ppx*G4px - 6*G2x*G4pxx + 12*G3p*G4pxx - 12*G4pp*G4pxx + 4*G2px*G4xx - 
           4*G3pp*G4xx + 2*G2xx*G5pp - 2*G3px*G5pp - 2*G2px*G5px + 2*G3pp*G5px + 6*(33*G4xx*G5x - 18*G5px*G5x + 13*(G4x - G5p)*G5xx)*pow(H,4) + 
           G3x*(3*G2xx - 6*G4ppx + 2*(81*G4xx - 47*G5px)*pow(H,2)) + pow(H,2)*
            (-456*G4px*G4xx - 12*G4p*G4xxx - 18*G3xx*G5p + 12*G4pxx*G5p + 36*G4xx*G5pp + 12*G5p*G5ppx + 6*G4x*(3*G3xx - 2*(G4pxx + G5ppx)) + 272*G4px*G5px - 24*G5pp*G5px + 4*G4p*G5pxx + 
              3*(G2xx + 2*G3px - 6*G4ppx)*G5x + 11*G2x*G5xx - 22*G3p*G5xx + 14*G4pp*G5xx + G5xxx*rptot)) - 
        8*pow(H,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
           6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2))*
         (6*pow(G3x,2) + 36*pow(G4px,2) - 2*G2xx*G4x + 8*G3px*G4x - 12*G4ppx*G4x + 6*G2x*G4xx - 12*G3p*G4xx + 12*G4pp*G4xx + G2xx*G5p - 4*G3px*G5p + 6*G4ppx*G5p - 6*G4px*G5pp - 3*G2x*G5px + 
           6*G3p*G5px - 6*G4pp*G5px + G2px*G5x - G3pp*G5x + 2*(54*G4x*G4xx + 6*G4*G4xxx - 57*G4xx*G5p - 25*G4x*G5px + 29*G5p*G5px - 4*G4*G5pxx - 36*G4px*G5x + 3*G5pp*G5x - 6*G4p*G5xx)*pow(H,2) + 
           30*pow(G5x,2)*pow(H,4) + 3*G3x*(-10*G4px + G5pp + 10*G5x*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 48*pow(H,4)*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2)))*
         (78*G3x*G4xx - 204*G4px*G4xx + 6*G3xx*G5p - 24*G4pxx*G5p + 24*G4xx*G5pp + 6*G5p*G5ppx - 42*G3x*G5px + 108*G4px*G5px - 12*G5pp*G5px + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx - 
           46*G5p*G5xx*pow(H,2) + 4*G4*G5xxx*pow(H,2) - 2*G5x*(G2xx - 4*G3px + 6*G4ppx + (-75*G4xx + 37*G5px)*pow(H,2)) - 2*G4x*(6*G3xx - 24*G4pxx + 6*G5ppx - 19*G5xx*pow(H,2)))*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 16*pow(H,4)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)))*
         (60*pow(G4xx,2) - 12*G4x*G4xxx + 6*G4xxx*G5p - 66*G4xx*G5px + 18*pow(G5px,2) + 8*G4x*G5pxx - 4*G5p*G5pxx + 9*G3x*G5xx - 24*G4px*G5xx + 3*G5pp*G5xx + 
           G5x*(-3*G3xx + 12*G4pxx - 3*G5ppx + 14*G5xx*pow(H,2)))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) + 
        8*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,4)*(2*G2p*G4 + 24*G4*G4p*pow(H,2) - 3*G4p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) + 16*pow(H,4)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
           6*G4px*G5xx + 2*G5x*G5xx*pow(H,2))*(-6*pow(G4p,2) + 2*G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2)) - 3*G4x*rptot + 3*G5p*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
           6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2)))*
         (2*G2xx*G4 - 8*G3px*G4 - 18*G3x*G4p + 12*G4*G4ppx + 48*G4p*G4px - 20*G3p*G4x + 12*G4pp*G4x + 2*G2x*(5*G4x - 4*G5p) + 16*G3p*G5p - 12*G4pp*G5p - 6*G4p*G5pp - G2p*G5x + 
           4*(21*pow(G4x,2) + 3*G4*G4xx - 39*G4x*G5p + 18*pow(G5p,2) - 4*G4*G5px - 12*G4p*G5x)*pow(H,2) + 6*G4xx*rptot - 3*G5px*rptot)*
         (-5*G2x*G4*G4p + 10*G3p*G4*G4p + 3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + 
              G4*G5x*pow(H,2)*rptot)) - 4*pow(H,2)*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
           G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2)))*
         (2*(-3*G2x*G4p + 6*G3p*G4p - 6*G4p*G4pp + G2p*G5p + 66*G4p*G5p*pow(H,2) - 2*G4x*(G2p + 39*G4p*pow(H,2))) - 4*G4*(G2px - G3pp + 3*pow(H,2)*(G3x - G5pp + G5x*pow(H,2))) + 
           3*(G3x - 2*G4px + 3*G5x*pow(H,2))*rptot)*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)) - 
        4*pow(H,2)*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2)))*
         (2*(-3*G3p*G3x + 3*G3x*G4pp + 6*G3p*G4px - 6*G4pp*G4px + 2*G2px*G4x - 2*G3pp*G4x - G2px*G5p + G3pp*G5p) + 2*(99*G4x*G5x - 96*G5p*G5x + 8*G4*G5xx)*pow(H,4) + 
           G2x*(3*G3x - 6*G4px + 13*G5x*pow(H,2)) + pow(H,2)*(-204*G4px*G4x - 108*G4p*G4xx + G3x*(90*G4x - 84*G5p) + 204*G4px*G5p + 12*G4x*G5pp - 18*G5p*G5pp + 12*G4*(G3xx - 4*G4pxx + G5ppx) + 
              60*G4p*G5px - 26*G3p*G5x + 18*G4pp*G5x + 3*G5xx*rptot))*(-5*G2x*G4*G4p + 10*G3p*G4*G4p + 
           3*(-5*pow(G4p,3) - 10*G4*G4p*G4x*pow(H,2) + 10*G4*G4p*G5p*pow(H,2) + G3x*G4*rptot - 3*G4*G4px*rptot - 2*G4p*G4x*rptot + 2*G4p*G5p*rptot + G4*G5x*pow(H,2)*rptot)))))/
 ((-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
     4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
     2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 
        30*G4x*G5p*pow(H,2) + 18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
     2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
        14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 
        6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
      (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
        2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + 
     pow(a,4)*pow(dphi,4)*(3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 
        12*G4x*G5px*pow(H,2) + 18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4)))*
   pow(2*a*pow(dphi,7)*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
     2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px + 4*G3p*G4x - 2*G3p*G5p + G2x*(-2*G4x + G5p) + 
        6*(2*pow(G4x,2) + 8*G4*G4xx - 5*G4x*G5p + 3*pow(G5p,2) - 5*G4*G5px - 3*G4p*G5x)*pow(H,2)) + 
     2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 
        6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) + 4*pow(a,8)*(3*pow(G4p,2) + G4*(G2x - 2*G3p + 6*(G4x - G5p)*pow(H,2))) + 
     24*pow(a,7)*dphi*H*(2*G4p*(-G4x + G5p) + G4*(G3x - 3*G4px + G5x*pow(H,2))) + 
     pow(a,4)*pow(dphi,4)*(3*pow(G3x,2) + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 
        6*(4*G4*G4xxx - 8*G4xx*G5p + 2*G4x*G5px + 3*G5p*G5px - 2*G4*G5pxx - 2*G4p*G5xx)*pow(H,2) + 15*pow(G5x,2)*pow(H,4) + 6*G3x*(-2*G4px + G5x*pow(H,2))) + 
     2*pow(a,3)*pow(dphi,5)*H*(12*G3x*G4xx - 24*G4px*G4xx + 3*G3xx*G5p - 6*G4pxx*G5p - 6*G3x*G5px + 12*G4px*G5px - 5*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2) + 
        G5x*(-G2xx + G3px + 3*(4*G4xx - G5px)*pow(H,2)) - 2*G4x*(3*G3xx - 6*G4pxx + G5xx*pow(H,2))) + 
     2*pow(a,5)*pow(dphi,3)*H*(-(G2x*G5x) + 2*(3*G4*(G3xx - 2*G4pxx) + 6*G4px*G4x - 12*G4p*G4xx - 3*G3x*G5p + 3*G4px*G5p + 6*G4p*G5px + G3p*G5x + (9*G4x*G5x - 12*G5p*G5x + 7*G4*G5xx)*pow(H,2))),2))\
 - (3*pow(a,3)*(-4*pow(ddphi,2)*pow(dphi,9)*G4*(pow(G5xxx,2) - G5xx*G5xxxx)*pow(H,5) + 
     2*a*ddphi*pow(dphi,8)*pow(H,4)*(2*ddphi*G4*(6*G4xxxx*G5xx - 3*G5pxxx*G5xx - 10*G4xxx*G5xxx + 5*G5pxx*G5xxx + 4*G4xx*G5xxxx - 2*G5px*G5xxxx) + ddphi*G4p*(3*G5xx*G5xxx - 2*G5x*G5xxxx) + 
        2*pow(dphi,2)*G4*(3*pow(G5xxx,2) - 2*G5xx*G5xxxx)*pow(H,2)) - 2*pow(a,2)*pow(dphi,7)*pow(H,3)*
      (2*pow(dphi,4)*G4*(2*pow(G5xxx,2) - G5xx*G5xxxx)*pow(H,4) + ddphi*dphi*H*
         (-4*dH*G4*G5xx*G5xxx + dphi*(4*G4*(6*G4xxxx*G5xx - 4*G5pxxx*G5xx - 16*G4xxx*G5xxx + 9*G5pxx*G5xxx + 4*G4xx*G5xxxx - 2*G5px*G5xxxx) + G4p*(9*G5xx*G5xxx - 4*G5x*G5xxxx))*H) + 
        2*pow(ddphi,2)*(G4p*(6*G4xxxx*G5x - 3*G5pxxx*G5x - 6*G4xxx*G5xx + 3*G5pxx*G5xx - 6*G4xx*G5xxx + 3*G5px*G5xxx + 2*G4x*G5xxxx - G5p*G5xxxx) + 
           G4*(24*pow(G4xxx,2) + 12*G4xxxx*G5px - 24*G4xxx*G5pxx + 6*pow(G5pxx,2) - 6*G5px*G5pxxx + 12*G4xx*(-2*G4xxxx + G5pxxx) - 3*G3xxx*G5xx + 6*G4pxxx*G5xx + 4*G3xx*G5xxx - 8*G4pxx*G5xxx - 
              G3x*G5xxxx + 2*G4px*G5xxxx + 2*G5xx*G5xxx*pow(H,2) - 3*G5x*G5xxxx*pow(H,2)))) + 
     2*pow(a,3)*pow(dphi,6)*pow(H,2)*(2*pow(dphi,3)*pow(H,3)*(dH*G4*G5xx*G5xxx + 
           dphi*(G4*(6*G4xxxx*G5xx - 5*G5pxxx*G5xx - 22*G4xxx*G5xxx + 14*G5pxx*G5xxx + 4*G4xx*G5xxxx - 2*G5px*G5xxxx) + G4p*(3*G5xx*G5xxx - G5x*G5xxxx))*H) + 
        pow(ddphi,2)*(6*G3xxx*(4*G4*G4xx - 2*G4*G5px - G4p*G5x) + G4p*(48*G4xx*G4xxx - 24*G4x*G4xxxx + 12*G4xxxx*G5p - 24*G4xxx*G5px - 24*G4xx*G5pxx + 12*G5px*G5pxx + 12*G4x*G5pxxx - 6*G5p*G5pxxx + 
              12*G4pxxx*G5x + 3*G3xx*G5xx - 6*G4pxx*G5xx + 3*G3x*G5xxx - 6*G4px*G5xxx + 24*pow(G5xx,2)*pow(H,2) - 17*G5x*G5xxx*pow(H,2)) - 
           2*G4*(24*G4pxxx*G4xx + 18*G3xx*G4xxx - 36*G4pxx*G4xxx - 6*G3x*G4xxxx + 12*G4px*G4xxxx - 12*G4pxxx*G5px - 9*G3xx*G5pxx + 18*G4pxx*G5pxx + 3*G3x*G5pxxx - 6*G4px*G5pxxx - G2xxx*G5xx + 
              G3pxx*G5xx + G2xx*G5xxx - G3px*G5xxx - 18*G4xxxx*G5x*pow(H,2) + 9*G5pxxx*G5x*pow(H,2) + 22*G4xxx*G5xx*pow(H,2) - 8*G5pxx*G5xx*pow(H,2) - 8*G4xx*G5xxx*pow(H,2) - 
              G5px*G5xxx*pow(H,2) - 4*G4x*G5xxxx*pow(H,2) + 4*G5p*G5xxxx*pow(H,2))) + 
        ddphi*dphi*H*(4*dH*(3*G4p*pow(G5xx,2) + 8*G4*G4xx*G5xxx - 4*G4*G5px*G5xxx - 3*G4p*G5x*G5xxx) + 
           dphi*H*(G4p*(24*G4xxxx*G5x - 16*G5pxxx*G5x - 42*G4xxx*G5xx + 27*G5pxx*G5xx - 36*G4xx*G5xxx + 18*G5px*G5xxx + 8*G4x*G5xxxx - 4*G5p*G5xxxx) + 
              2*G4*(84*pow(G4xxx,2) - 48*G4xx*G4xxxx + 24*G4xxxx*G5px - 96*G4xxx*G5pxx + 27*pow(G5pxx,2) + 32*G4xx*G5pxxx - 16*G5px*G5pxxx - 6*G3xxx*G5xx + 24*G4pxxx*G5xx - 6*G5ppxx*G5xx + 
                 14*G3xx*G5xxx - 36*G4pxx*G5xxx + 4*G5ppx*G5xxx - 2*G3x*G5xxxx + 4*G4px*G5xxxx + 21*G5xx*G5xxx*pow(H,2) - 6*G5x*G5xxxx*pow(H,2))))) - 
     2*pow(a,4)*pow(dphi,5)*H*(2*pow(ddphi,2)*(3*pow(G3xx,2)*G4 + 6*G3xxx*G4*G4px + 12*G4*pow(G4pxx,2) - 12*G4*G4px*G4pxxx + 6*G3xxx*G4p*G4x - 12*G4p*G4pxxx*G4x - 4*G2xxx*G4*G4xx + 
           4*G3pxx*G4*G4xx + 12*G4p*G4pxx*G4xx + 4*G2xx*G4*G4xxx - 4*G3px*G4*G4xxx + 12*G4p*G4px*G4xxx - 3*G3xxx*G4p*G5p + 6*G4p*G4pxxx*G5p + 2*G2xxx*G4*G5px - 2*G3pxx*G4*G5px - 6*G4p*G4pxx*G5px - 
           2*G2xx*G4*G5pxx + 2*G3px*G4*G5pxx - 6*G4p*G4px*G5pxx + G2xxx*G4p*G5x - G3pxx*G4p*G5x - 24*G4*G4x*G4xxxx*pow(H,2) + 24*G4*G4xxxx*G5p*pow(H,2) - 24*G4*G4xxx*G5px*pow(H,2) + 
           12*G4*G4xx*G5pxx*pow(H,2) + 6*G4*G5px*G5pxx*pow(H,2) + 12*G4*G4x*G5pxxx*pow(H,2) - 12*G4*G5p*G5pxxx*pow(H,2) - 9*G3xxx*G4*G5x*pow(H,2) + 18*G4*G4pxxx*G5x*pow(H,2) + 
           36*G4p*G4xxx*G5x*pow(H,2) - 21*G4p*G5pxx*G5x*pow(H,2) - 20*G4*G4pxx*G5xx*pow(H,2) - 78*G4p*G4xx*G5xx*pow(H,2) + 42*G4p*G5px*G5xx*pow(H,2) + 20*G4p*G4x*G5xxx*pow(H,2) - 
           7*G4p*G5p*G5xxx*pow(H,2) + 22*G4*pow(G5xx,2)*pow(H,4) - 24*G4*G5x*G5xxx*pow(H,4) + G3xx*(-12*G4*G4pxx - 6*G4p*G4xx + 3*G4p*G5px + 13*G4*G5xx*pow(H,2)) + 
           G3x*(-3*G3xxx*G4 + 6*G4*G4pxxx - 6*G4p*G4xxx + 3*G4p*G5pxx - 4*G4*G5xxx*pow(H,2))) - 
        ddphi*dphi*H*(2*dH*(-6*G4p*(4*G4xxx*G5x - 2*G5pxx*G5x - 6*G4xx*G5xx + 3*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx) + 
              G4*(48*G4xx*G4xxx - 24*G4xxx*G5px - 24*G4xx*G5pxx + 12*G5px*G5pxx - 6*G3xx*G5xx + 12*G4pxx*G5xx + 6*G3x*G5xxx - 12*G4px*G5xxx + 17*pow(G5xx,2)*pow(H,2) + 4*G5x*G5xxx*pow(H,2))) + 
           dphi*H*(12*G3xxx*(-4*G4*G4xx + 2*G4*G5px + G4p*G5x) + G4p*(-168*G4xx*G4xxx + 48*G4x*G4xxxx - 24*G4xxxx*G5p + 84*G4xxx*G5px + 108*G4xx*G5pxx - 54*G5px*G5pxx - 32*G4x*G5pxxx + 16*G5p*G5pxxx - 
                 48*G4pxxx*G5x + 12*G5ppxx*G5x - 15*G3xx*G5xx + 54*G4pxx*G5xx - 12*G5ppx*G5xx - 9*G3x*G5xxx + 18*G4px*G5xxx - 78*pow(G5xx,2)*pow(H,2) + 25*G5x*G5xxx*pow(H,2)) + 
              2*G4*(96*G4pxxx*G4xx + 72*G3xx*G4xxx - 192*G4pxx*G4xxx - 12*G3x*G4xxxx + 24*G4px*G4xxxx + 24*G4xxx*G5ppx - 24*G4xx*G5ppxx - 48*G4pxxx*G5px + 12*G5ppxx*G5px - 42*G3xx*G5pxx + 
                 108*G4pxx*G5pxx - 12*G5ppx*G5pxx + 8*G3x*G5pxxx - 16*G4px*G5pxxx - 2*G2xxx*G5xx + 8*G3pxx*G5xx - 12*G4ppxx*G5xx + 4*G2xx*G5xxx - 6*G3px*G5xxx + 4*G4ppx*G5xxx - 36*G4xxxx*G5x*pow(H,2) + 
                 24*G5pxxx*G5x*pow(H,2) + 146*G4xxx*G5xx*pow(H,2) - 67*G5pxx*G5xx*pow(H,2) + 46*G4xx*G5xxx*pow(H,2) - 46*G5px*G5xxx*pow(H,2) - 8*G4x*G5xxxx*pow(H,2) + 8*G5p*G5xxxx*pow(H,2)))
           ) + 2*pow(dphi,2)*pow(H,2)*(pow(dH,2)*G4*(-6*pow(G5xx,2) + 2*G5x*G5xxx) + 
           dH*dphi*(G4*(-18*G4xxx*G5xx + 3*G5pxx*G5xx + 4*G4xx*G5xxx + 2*G5px*G5xxx) + G4p*(12*pow(G5xx,2) - 7*G5x*G5xxx))*H + 
           pow(dphi,2)*pow(H,2)*(G4p*(6*G4xxxx*G5x - 5*G5pxxx*G5x - 15*G4xxx*G5xx + 12*G5pxx*G5xx - 12*G4xx*G5xxx + 6*G5px*G5xxx + 2*G4x*G5xxxx - G5p*G5xxxx) + 
              G4*(60*pow(G4xxx,2) - 24*G4xx*G4xxxx + 12*G4xxxx*G5px - 78*G4xxx*G5pxx + 24*pow(G5pxx,2) + 20*G4xx*G5pxxx - 10*G5px*G5pxxx - 3*G3xxx*G5xx + 18*G4pxxx*G5xx - 7*G5ppxx*G5xx + 
                 10*G3xx*G5xxx - 34*G4pxx*G5xxx + 8*G5ppx*G5xxx - G3x*G5xxxx + 2*G4px*G5xxxx + 18*G5xx*G5xxx*pow(H,2) - 3*G5x*G5xxxx*pow(H,2))))) - 
     24*pow(a,13)*pow(H,2)*(3*G3x*G4 - 9*G4*G4px - 7*G4p*G4x + 7*G4p*G5p + 3*G4*G5x*pow(H,2))*rptot - 
     2*pow(a,12)*H*(12*(drptot*G4p*(-G4x + G5p) + 12*dH*G4*H*(G3x*G4 - 3*G4*G4px - G4p*G4x + G4p*G5p + G4*G5x*pow(H,2))) + 
        dphi*(2*pow(G2x,2)*G4 + 8*pow(G3p,2)*G4 + 36*pow(G4p,2)*G4pp + 8*G2pp*G4*G4x - 8*G2pp*G4*G5p + 60*pow(G4p,2)*G4x*pow(H,2) + 168*G4*G4pp*G4x*pow(H,2) - 60*pow(G4p,2)*G5p*pow(H,2) - 
           168*G4*G4pp*G5p*pow(H,2) + 216*G4*pow(G4x,2)*pow(H,4) - 432*G4*G4x*G5p*pow(H,4) + 216*G4*pow(G5p,2)*pow(H,4) + 6*G2xx*G4*rptot - 6*G3px*G4*rptot - 21*G3x*G4p*rptot + 
           42*G4p*G4px*rptot + 48*pow(G4x,2)*pow(H,2)*rptot + 144*G4*G4xx*pow(H,2)*rptot - 120*G4x*G5p*pow(H,2)*rptot + 72*pow(G5p,2)*pow(H,2)*rptot - 90*G4*G5px*pow(H,2)*rptot - 
           63*G4p*G5x*pow(H,2)*rptot + G2x*(-8*G3p*G4 + 6*pow(G4p,2) + 12*G4*G4pp + 48*G4*G4x*pow(H,2) - 48*G4*G5p*pow(H,2) - 8*G4x*rptot + 4*G5p*rptot) - 
           4*G3p*(3*pow(G4p,2) + 6*G4*(G4pp + 4*(G4x - G5p)*pow(H,2)) + 2*(-2*G4x + G5p)*rptot))) + 
     pow(a,5)*pow(dphi,4)*(2*pow(ddphi,2)*(-2*G2xx*G3xx*G4 + 2*G3px*G3xx*G4 + 3*G3x*G3xx*G4p - 6*G3xx*G4p*G4px + 4*G2xx*G4*G4pxx - 4*G3px*G4*G4pxx - 6*G3x*G4p*G4pxx + 12*G4p*G4px*G4pxx + 
           24*G3xxx*G4*G4x*pow(H,2) - 48*G4*G4pxxx*G4x*pow(H,2) - 24*G3xx*G4*G4xx*pow(H,2) + 240*G4p*pow(G4xx,2)*pow(H,2) + 24*G3x*G4*G4xxx*pow(H,2) + 24*G4*G4px*G4xxx*pow(H,2) - 
           168*G4p*G4x*G4xxx*pow(H,2) - 24*G3xxx*G4*G5p*pow(H,2) + 48*G4*G4pxxx*G5p*pow(H,2) + 60*G4p*G4xxx*G5p*pow(H,2) + 30*G3xx*G4*G5px*pow(H,2) - 36*G4*G4pxx*G5px*pow(H,2) - 
           264*G4p*G4xx*G5px*pow(H,2) + 72*G4p*pow(G5px,2)*pow(H,2) - 18*G3x*G4*G5pxx*pow(H,2) + 96*G4p*G4x*G5pxx*pow(H,2) - 36*G4p*G5p*G5pxx*pow(H,2) - 27*G3xx*G4p*G5x*pow(H,2) + 
           66*G4p*G4pxx*G5x*pow(H,2) - 8*G2xx*G4*G5xx*pow(H,2) + 6*G3px*G4*G5xx*pow(H,2) + 33*G3x*G4p*G5xx*pow(H,2) - 72*G4p*G4px*G5xx*pow(H,2) - 2*G2x*G4*G5xxx*pow(H,2) + 
           4*G3p*G4*G5xxx*pow(H,2) - 6*pow(G4p,2)*G5xxx*pow(H,2) + 168*G4*G4xxx*G5x*pow(H,4) - 102*G4*G5pxx*G5x*pow(H,4) - 236*G4*G4xx*G5xx*pow(H,4) + 170*G4*G5px*G5xx*pow(H,4) + 
           31*G4p*G5x*G5xx*pow(H,4) + 76*G4*G4x*G5xxx*pow(H,4) - 76*G4*G5p*G5xxx*pow(H,4) + 2*G2xxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) - 
           2*G3pxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2))) + 
        pow(dphi,2)*pow(H,2)*(-12*pow(dH,2)*(G4p*G5x*G5xx + 2*G4*(2*G4xxx*G5x - G5pxx*G5x - 6*G4xx*G5xx + 3*G5px*G5xx)) + 
           4*dH*dphi*H*(G4p*(30*G4xxx*G5x - 21*G5pxx*G5x - 84*G4xx*G5xx + 48*G5px*G5xx + 14*G4x*G5xxx - 7*G5p*G5xxx) + 
              G4*(-36*G4xxx*G5px + 6*G5px*G5pxx + 12*G4xx*(2*G4xxx + G5pxx) + 15*G3xx*G5xx - 6*G4pxx*G5xx - 12*G5ppx*G5xx - 3*G3x*G5xxx - 2*G4px*G5xxx + 4*G5pp*G5xxx + 7*pow(G5xx,2)*pow(H,2) + 
                 9*G5x*G5xxx*pow(H,2))) + dphi*pow(H,2)*(-3*drptot*pow(G5xx,2) + 2*drptot*G5x*G5xxx + 
              4*dphi*(3*G3xxx*(4*G4*G4xx - 2*G4*G5px - G4p*G5x) + G4p*(60*G4xx*G4xxx - 12*G4x*G4xxxx + 6*G4xxxx*G5p - 30*G4xxx*G5px - 48*G4xx*G5pxx + 24*G5px*G5pxx + 10*G4x*G5pxxx - 5*G5p*G5pxxx + 
                    18*G4pxxx*G5x - 7*G5ppxx*G5x + 6*G3xx*G5xx - 33*G4pxx*G5xx + 12*G5ppx*G5xx + 3*G3x*G5xxx - 6*G4px*G5xxx + 21*pow(G5xx,2)*pow(H,2) - G5x*G5xxx*pow(H,2)) + 
                 G4*(-54*G3xx*G4xxx + 192*G4pxx*G4xxx + 6*G3x*G4xxxx - 12*G4px*G4xxxx - 48*G4xxx*G5ppx + 28*G4xx*G5ppxx - 14*G5ppxx*G5px + 36*G4pxxx*(-2*G4xx + G5px) + 36*G3xx*G5pxx - 114*G4pxx*G5pxx + 
                    24*G5ppx*G5pxx - 5*G3x*G5pxxx + 10*G4px*G5pxxx + G2xxx*G5xx - 7*G3pxx*G5xx + 18*G4ppxx*G5xx - 3*G5pppx*G5xx - 3*G2xx*G5xxx + 8*G3px*G5xxx - 14*G4ppx*G5xxx + 2*G5ppp*G5xxx + 
                    18*G4xxxx*G5x*pow(H,2) - 15*G5pxxx*G5x*pow(H,2) - 118*G4xxx*G5xx*pow(H,2) + 66*G5pxx*G5xx*pow(H,2) - 48*G4xx*G5xxx*pow(H,2) + 44*G5px*G5xxx*pow(H,2) + 
                    4*G4x*G5xxxx*pow(H,2) - 4*G5p*G5xxxx*pow(H,2))))) + ddphi*dphi*H*
         (8*dH*(3*G4p*(8*pow(G4xx,2) - 8*G4x*G4xxx + 4*G4xxx*G5p - 8*G4xx*G5px + 2*pow(G5px,2) + 4*G4x*G5pxx - 2*G5p*G5pxx - G3xx*G5x + 2*G4pxx*G5x + G3x*G5xx - 2*G4px*G5xx - G5x*G5xx*pow(H,2)) + 
              G4*(12*G3x*G4xxx - 24*G4px*G4xxx - 6*G3x*G5pxx + 12*G4px*G5pxx - 2*G2xx*G5xx + 2*G3px*G5xx - 6*G4xxx*G5x*pow(H,2) + 3*G5pxx*G5x*pow(H,2) + 68*G4xx*G5xx*pow(H,2) - 
                 34*G5px*G5xx*pow(H,2) + 6*G4x*G5xxx*pow(H,2) - 7*G5p*G5xxx*pow(H,2))) + 
           dphi*H*(60*pow(G3xx,2)*G4 + 48*G3xxx*G4*G4px + 432*G4*pow(G4pxx,2) - 192*G4*G4px*G4pxxx + 48*G3xxx*G4p*G4x - 192*G4p*G4pxxx*G4x - 32*G2xxx*G4*G4xx + 128*G3pxx*G4*G4xx - 
              192*G4*G4ppxx*G4xx + 432*G4p*G4pxx*G4xx + 80*G2xx*G4*G4xxx - 128*G3px*G4*G4xxx + 96*G4*G4ppx*G4xxx + 168*G4p*G4px*G4xxx - 24*G3xxx*G4p*G5p + 96*G4p*G4pxxx*G5p - 96*G4*G4pxx*G5ppx - 
              96*G4p*G4xx*G5ppx + 48*G4*G4px*G5ppxx + 48*G4p*G4x*G5ppxx - 24*G4p*G5p*G5ppxx + 16*G2xxx*G4*G5px - 64*G3pxx*G4*G5px + 96*G4*G4ppxx*G5px - 216*G4p*G4pxx*G5px + 48*G4p*G5ppx*G5px - 
              48*G2xx*G4*G5pxx + 72*G3px*G4*G5pxx - 48*G4*G4ppx*G5pxx - 108*G4p*G4px*G5pxx + 8*G2xxx*G4p*G5x - 32*G3pxx*G4p*G5x + 48*G4p*G4ppxx*G5x + 8*G2pxx*G4*G5xx - 8*G3ppx*G4*G5xx - 6*G2xx*G4p*G5xx + 
              18*G3px*G4p*G5xx - 24*G4p*G4ppx*G5xx + 1488*G4*G4xx*G4xxx*pow(H,2) - 192*G4*G4x*G4xxxx*pow(H,2) + 192*G4*G4xxxx*G5p*pow(H,2) - 1248*G4*G4xxx*G5px*pow(H,2) - 
              616*G4*G4xx*G5pxx*pow(H,2) + 584*G4*G5px*G5pxx*pow(H,2) + 128*G4*G4x*G5pxxx*pow(H,2) - 128*G4*G5p*G5pxxx*pow(H,2) - 72*G3xxx*G4*G5x*pow(H,2) + 288*G4*G4pxxx*G5x*pow(H,2) + 
              180*G4p*G4xxx*G5x*pow(H,2) - 72*G4*G5ppxx*G5x*pow(H,2) - 134*G4p*G5pxx*G5x*pow(H,2) - 600*G4*G4pxx*G5xx*pow(H,2) - 1092*G4p*G4xx*G5xx*pow(H,2) - 8*G4*G5ppx*G5xx*pow(H,2) + 
              630*G4p*G5px*G5xx*pow(H,2) - 232*G4*G4px*G5xxx*pow(H,2) + 136*G4p*G4x*G5xxx*pow(H,2) - 32*G4p*G5p*G5xxx*pow(H,2) + 40*G4*G5pp*G5xxx*pow(H,2) + 456*G4*pow(G5xx,2)*pow(H,4) - 
              108*G4*G5x*G5xxx*pow(H,4) - 4*G3xx*(84*G4*G4pxx + 30*G4p*G4xx - 12*G4*G5ppx - 15*G4p*G5px - 77*G4*G5xx*pow(H,2)) - 
              6*G3x*(4*G3xxx*G4 - 16*G4*G4pxxx + 14*G4p*G4xxx + 4*G4*G5ppxx - 9*G4p*G5pxx - 6*G4*G5xxx*pow(H,2)) + 3*pow(G5xx,2)*pow(H,2)*rptot - 2*G5x*G5xxx*pow(H,2)*rptot))) + 
     2*pow(a,8)*dphi*(12*pow(dH,2)*pow(dphi,2)*H*(3*G4p*(-2*G4x + G5p)*G5x + G4*(16*G4x*G4xx - 16*G4xx*G5p - 8*G4x*G5px + 8*G5p*G5px + G3x*G5x + 7*pow(G5x,2)*pow(H,2))) + 
        2*dH*dphi*(ddphi*(15*pow(G3x,2)*G4 + 12*pow(G4,2)*(2*G4xxx - G5pxx)*pow(H,2) - 6*G3x*(12*G4*G4px + 2*G4p*G4x - G4*G5x*pow(H,2)) - 
              12*G4p*(G4px*(-4*G4x + G5p) + G4p*(2*G4xx - G5px) - (G4x - 4*G5p)*G5x*pow(H,2)) + 
              G4*(84*pow(G4px,2) - 12*G2xx*G4x + 12*G3px*G4x - 8*G2x*G4xx + 16*G3p*G4xx + 10*G2xx*G5p - 10*G3px*G5p + 4*G2x*G5px - 8*G3p*G5px + 240*G4x*G4xx*pow(H,2) - 288*G4xx*G5p*pow(H,2) - 
                 132*G4x*G5px*pow(H,2) + 162*G5p*G5px*pow(H,2) + 36*G4px*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 51*pow(G5x,2)*pow(H,4))) + 
           pow(dphi,2)*(-6*pow(G3x,2)*G4p - 12*G3x*G4*G4ppx + 36*G3x*G4p*G4px + 24*G4*G4ppx*G4px - 48*G4p*pow(G4px,2) + 24*G4p*G4ppx*G4x - 12*G4p*G4ppx*G5p - 6*G3x*G4p*G5pp + 12*G4p*G4px*G5pp + 
              72*G3xx*G4*G4x*pow(H,2) - 48*G4*G4pxx*G4x*pow(H,2) + 78*G3x*G4*G4xx*pow(H,2) - 372*G4*G4px*G4xx*pow(H,2) - 132*G4p*G4x*G4xx*pow(H,2) - 54*G3xx*G4*G5p*pow(H,2) + 
              12*G4*G4pxx*G5p*pow(H,2) + 210*G4p*G4xx*G5p*pow(H,2) + 24*G4*G4xx*G5pp*pow(H,2) - 48*G4*G4x*G5ppx*pow(H,2) + 48*G4*G5p*G5ppx*pow(H,2) - 45*G3x*G4*G5px*pow(H,2) + 
              246*G4*G4px*G5px*pow(H,2) + 30*G4p*G4x*G5px*pow(H,2) - 99*G4p*G5p*G5px*pow(H,2) - 24*G4*G5pp*G5px*pow(H,2) - 51*G3x*G4p*G5x*pow(H,2) - 36*G4*G4ppx*G5x*pow(H,2) + 
              108*G4p*G4px*G5x*pow(H,2) + 10*G2x*G4*G5xx*pow(H,2) - 20*G3p*G4*G5xx*pow(H,2) + 24*pow(G4p,2)*G5xx*pow(H,2) - 6*G4*G4pp*G5xx*pow(H,2) + 354*G4*G4xx*G5x*pow(H,4) - 
              201*G4*G5px*G5x*pow(H,4) - 93*G4p*pow(G5x,2)*pow(H,4) + 80*G4*G4x*G5xx*pow(H,4) - 38*G4*G5p*G5xx*pow(H,4) - 12*pow(G4,2)*G5xxx*pow(H,4) + 
              G2xx*(3*G3x*G4 - 14*G4*G4px + 2*G4p*G4x - G4p*G5p + 4*G4*G5pp + 27*G4*G5x*pow(H,2)) + G3px*(3*G3x*G4 + 7*G4p*(-2*G4x + G5p) + G4*(2*G4px - 4*G5pp - 9*G5x*pow(H,2))))) - 
        H*(-(pow(dphi,3)*drptot*(-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 
                3*G5px*G5x*pow(H,2) + 2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2))) + 
           2*pow(ddphi,2)*(6*pow(G3x,2)*G4 - 6*G3xx*G4*G4p + 72*G4*pow(G4px,2) + 12*G4*G4p*G4pxx - 8*G2xx*G4*G4x + 12*G3px*G4*G4x + 16*G2x*G4*G4xx - 32*G3p*G4*G4xx + 72*pow(G4p,2)*G4xx + 
              8*G2xx*G4*G5p - 12*G3px*G4*G5p - 36*G4p*G4px*G5p - 10*G2x*G4*G5px + 20*G3p*G4*G5px - 42*pow(G4p,2)*G5px + G2x*G4p*G5x - 2*G3p*G4p*G5x - 24*G4*G4x*G4xx*pow(H,2) + 
              24*G4*G4xx*G5p*pow(H,2) + 24*G4*G4x*G5px*pow(H,2) - 24*G4*G5p*G5px*pow(H,2) - 78*G4*G4px*G5x*pow(H,2) - 60*G4p*G4x*G5x*pow(H,2) + 72*G4p*G5p*G5x*pow(H,2) - 
              14*G4*G4p*G5xx*pow(H,2) + 18*G4*pow(G5x,2)*pow(H,4) - 6*G3x*(7*G4*G4px + G4p*G4x - 3*G4p*G5p - 4*G4*G5x*pow(H,2))) + 
           pow(dphi,4)*(2*pow(G2xx,2)*G4 + 12*pow(G3px,2)*G4 + 4*G2pxx*G3x*G4 - 10*G3ppx*G3x*G4 - 6*G2px*G3xx*G4 + 12*G3pp*G3xx*G4 - 12*G3xx*G4*G4ppp + 12*G3x*G4*G4pppx - 42*G3x*G4p*G4ppx - 
              8*G2pxx*G4*G4px + 20*G3ppx*G4*G4px - 24*G4*G4pppx*G4px + 84*G4p*G4ppx*G4px + 12*G2px*G4*G4pxx - 24*G3pp*G4*G4pxx + 24*G4*G4ppp*G4pxx - 8*G2pxx*G4p*G4x + 20*G3ppx*G4p*G4x - 
              24*G4p*G4pppx*G4x - 8*G2ppx*G4*G4xx + 8*G3ppp*G4*G4xx + 12*G2px*G4p*G4xx - 24*G3pp*G4p*G4xx + 24*G4p*G4ppp*G4xx + 4*G2pxx*G4p*G5p - 10*G3ppx*G4p*G5p + 12*G4p*G4pppx*G5p + 6*G3x*G4p*G5ppp - 
              12*G4p*G4px*G5ppp + 4*G2ppx*G4*G5px - 4*G3ppp*G4*G5px - 6*G2px*G4p*G5px + 12*G3pp*G4p*G5px - 12*G4p*G4ppp*G5px + 2*G2ppx*G4p*G5x - 2*G3ppp*G4p*G5x + 78*G3x*G3xx*G4*pow(H,2) - 
              336*G3xx*G4*G4px*pow(H,2) - 252*G3x*G4*G4pxx*pow(H,2) + 1032*G4*G4px*G4pxx*pow(H,2) - 8*G2xxx*G4*G4x*pow(H,2) + 56*G3pxx*G4*G4x*pow(H,2) - 12*G3xx*G4p*G4x*pow(H,2) - 
              144*G4*G4ppxx*G4x*pow(H,2) + 48*G4p*G4pxx*G4x*pow(H,2) - 216*G3x*G4p*G4xx*pow(H,2) + 528*G4*G4ppx*G4xx*pow(H,2) + 696*G4p*G4px*G4xx*pow(H,2) + 32*G2x*G4*G4xxx*pow(H,2) - 
              64*G3p*G4*G4xxx*pow(H,2) + 60*pow(G4p,2)*G4xxx*pow(H,2) + 72*G4*G4pp*G4xxx*pow(H,2) + 8*G2xxx*G4*G5p*pow(H,2) - 56*G3pxx*G4*G5p*pow(H,2) + 30*G3xx*G4p*G5p*pow(H,2) + 
              144*G4*G4ppxx*G5p*pow(H,2) - 156*G4p*G4pxx*G5p*pow(H,2) + 84*G3xx*G4*G5pp*pow(H,2) - 168*G4*G4pxx*G5pp*pow(H,2) - 168*G4p*G4xx*G5pp*pow(H,2) - 72*G4*G4xx*G5ppp*pow(H,2) + 
              24*G4*G4x*G5pppx*pow(H,2) - 24*G4*G5p*G5pppx*pow(H,2) + 50*G3x*G4*G5ppx*pow(H,2) - 196*G4*G4px*G5ppx*pow(H,2) - 4*G4p*G4x*G5ppx*pow(H,2) + 50*G4p*G5p*G5ppx*pow(H,2) + 
              150*G3x*G4p*G5px*pow(H,2) - 348*G4*G4ppx*G5px*pow(H,2) - 432*G4p*G4px*G5px*pow(H,2) + 84*G4p*G5pp*G5px*pow(H,2) + 48*G4*G5ppp*G5px*pow(H,2) - 22*G2x*G4*G5pxx*pow(H,2) + 
              44*G3p*G4*G5pxx*pow(H,2) - 48*pow(G4p,2)*G5pxx*pow(H,2) - 36*G4*G4pp*G5pxx*pow(H,2) + 12*G2pxx*G4*G5x*pow(H,2) - 30*G3ppx*G4*G5x*pow(H,2) + 36*G4*G4pppx*G5x*pow(H,2) - 
              90*G4p*G4ppx*G5x*pow(H,2) + 12*G4p*G5ppp*G5x*pow(H,2) - 16*G2px*G4*G5xx*pow(H,2) + 32*G3pp*G4*G5xx*pow(H,2) - 3*G2x*G4p*G5xx*pow(H,2) + 6*G3p*G4p*G5xx*pow(H,2) - 
              18*G4p*G4pp*G5xx*pow(H,2) - 28*G4*G4ppp*G5xx*pow(H,2) + 1872*G4*pow(G4xx,2)*pow(H,4) - 48*G4*G4x*G4xxx*pow(H,4) + 48*G4*G4xxx*G5p*pow(H,4) - 2736*G4*G4xx*G5px*pow(H,4) + 
              984*G4*pow(G5px,2)*pow(H,4) + 100*G4*G4x*G5pxx*pow(H,4) - 100*G4*G5p*G5pxx*pow(H,4) + 66*G3xx*G4*G5x*pow(H,4) - 84*G4*G4pxx*G5x*pow(H,4) - 492*G4p*G4xx*G5x*pow(H,4) - 
              42*G4*G5ppx*G5x*pow(H,4) + 360*G4p*G5px*G5x*pow(H,4) + 266*G3x*G4*G5xx*pow(H,4) - 1032*G4*G4px*G5xx*pow(H,4) - 146*G4p*G4x*G5xx*pow(H,4) + 172*G4p*G5p*G5xx*pow(H,4) + 
              208*G4*G5pp*G5xx*pow(H,4) + 270*G4*G5x*G5xx*pow(H,6) - G2xx*(14*G3px*G4 + 3*G3x*G4p - 28*G4*G4ppx - 6*G4p*G4px + 4*G4*G5ppp - 120*G4*G4xx*pow(H,2) + 94*G4*G5px*pow(H,2) + 
                 7*G4p*G5x*pow(H,2)) + 2*G3px*(9*G3x*G4p - 14*G4*G4ppx - 18*G4p*G4px + 2*G4*G5ppp - 164*G4*G4xx*pow(H,2) + 114*G4*G5px*pow(H,2) + 22*G4p*G5x*pow(H,2)) + 
              96*pow(G4xx,2)*pow(H,2)*rptot - 48*G4x*G4xxx*pow(H,2)*rptot + 24*G4xxx*G5p*pow(H,2)*rptot - 96*G4xx*G5px*pow(H,2)*rptot + 24*pow(G5px,2)*pow(H,2)*rptot + 
              24*G4x*G5pxx*pow(H,2)*rptot - 12*G5p*G5pxx*pow(H,2)*rptot - 12*G3xx*G5x*pow(H,2)*rptot + 24*G4pxx*G5x*pow(H,2)*rptot + 12*G3x*G5xx*pow(H,2)*rptot - 24*G4px*G5xx*pow(H,2)*rptot + 
              8*G5x*G5xx*pow(H,4)*rptot) + ddphi*pow(dphi,2)*(-22*G2x*G3xx*G4 + 44*G3p*G3xx*G4 + 51*pow(G3x,2)*G4p - 30*G3xx*pow(G4p,2) - 36*G3xx*G4*G4pp + 24*G3x*G4*G4ppx - 270*G3x*G4p*G4px + 
              336*G4p*pow(G4px,2) + 60*G2x*G4*G4pxx - 120*G3p*G4*G4pxx + 108*pow(G4p,2)*G4pxx + 72*G4*G4pp*G4pxx - 16*G2pxx*G4*G4x + 16*G3ppx*G4*G4x - 96*G4p*G4ppx*G4x - 16*G2px*G4*G4xx + 
              32*G3pp*G4*G4xx + 24*G2x*G4p*G4xx - 48*G3p*G4p*G4xx + 72*G4p*G4pp*G4xx + 16*G2pxx*G4*G5p - 16*G3ppx*G4*G5p + 24*G4p*G4ppx*G5p + 30*G3x*G4p*G5pp - 60*G4p*G4px*G5pp - 8*G2x*G4*G5ppx + 
              16*G3p*G4*G5ppx - 24*pow(G4p,2)*G5ppx + 8*G2px*G4*G5px - 16*G3pp*G4*G5px - 12*G2x*G4p*G5px + 24*G3p*G4p*G5px - 36*G4p*G4pp*G5px + 4*G2px*G4p*G5x - 8*G3pp*G4p*G5x + 
              12*G3xx*G4*G4x*pow(H,2) - 408*G4*G4pxx*G4x*pow(H,2) - 1248*G3x*G4*G4xx*pow(H,2) + 4200*G4*G4px*G4xx*pow(H,2) + 744*G4p*G4x*G4xx*pow(H,2) + 120*G4*G4p*G4xxx*pow(H,2) - 
              12*G3xx*G4*G5p*pow(H,2) + 408*G4*G4pxx*G5p*pow(H,2) - 1056*G4p*G4xx*G5p*pow(H,2) - 384*G4*G4xx*G5pp*pow(H,2) + 192*G4*G4x*G5ppx*pow(H,2) - 192*G4*G5p*G5ppx*pow(H,2) + 
              822*G3x*G4*G5px*pow(H,2) - 2796*G4*G4px*G5px*pow(H,2) - 396*G4p*G4x*G5px*pow(H,2) + 624*G4p*G5p*G5px*pow(H,2) + 252*G4*G5pp*G5px*pow(H,2) - 60*G4*G4p*G5pxx*pow(H,2) + 
              204*G3x*G4p*G5x*pow(H,2) + 168*G4*G4ppx*G5x*pow(H,2) - 534*G4p*G4px*G5x*pow(H,2) + 66*G4p*G5pp*G5x*pow(H,2) - 76*G2x*G4*G5xx*pow(H,2) + 152*G3p*G4*G5xx*pow(H,2) - 
              186*pow(G4p,2)*G5xx*pow(H,2) - 84*G4*G4pp*G5xx*pow(H,2) - 1392*G4*G4xx*G5x*pow(H,4) + 846*G4*G5px*G5x*pow(H,4) + 249*G4p*pow(G5x,2)*pow(H,4) - 248*G4*G4x*G5xx*pow(H,4) + 
              248*G4*G5p*G5xx*pow(H,4) + 6*G3xx*G4x*rptot - 12*G4pxx*G4x*rptot - 12*G3x*G4xx*rptot + 24*G4px*G4xx*rptot - 3*G3xx*G5p*rptot + 6*G4pxx*G5p*rptot + 6*G3x*G5px*rptot - 12*G4px*G5px*rptot - 
              12*G4xx*G5x*pow(H,2)*rptot + 3*G5px*G5x*pow(H,2)*rptot + 2*G4x*G5xx*pow(H,2)*rptot + 5*G5p*G5xx*pow(H,2)*rptot - 2*G4*G5xxx*pow(H,2)*rptot + 
              G3px*(26*G3x*G4 - 132*G4*G4px + 52*G4p*G4x - 8*G4p*G5p + 20*G4*G5pp - 30*G4*G5x*pow(H,2) - G5x*rptot) + 
              G2xx*(-30*G3x*G4 + 116*G4*G4px - 20*G4p*G4x + 4*G4p*G5p - 20*G4*G5pp - 30*G4*G5x*pow(H,2) + G5x*rptot)))) + 
     2*pow(a,7)*pow(dphi,2)*(pow(ddphi,2)*(-2*G2x*G3xx*G4 + 4*G3p*G3xx*G4 + 9*pow(G3x,2)*G4p - 6*G3xx*pow(G4p,2) - 42*G3x*G4p*G4px + 48*G4p*pow(G4px,2) + 4*G2x*G4*G4pxx - 8*G3p*G4*G4pxx + 
           12*pow(G4p,2)*G4pxx + 84*G3xx*G4*G4x*pow(H,2) - 216*G4*G4pxx*G4x*pow(H,2) - 132*G3x*G4*G4xx*pow(H,2) + 504*G4*G4px*G4xx*pow(H,2) + 24*G4p*G4x*G4xx*pow(H,2) + 
           24*G4*G4p*G4xxx*pow(H,2) - 84*G3xx*G4*G5p*pow(H,2) + 216*G4*G4pxx*G5p*pow(H,2) - 180*G4p*G4xx*G5p*pow(H,2) + 84*G3x*G4*G5px*pow(H,2) - 324*G4*G4px*G5px*pow(H,2) + 
           12*G4p*G4x*G5px*pow(H,2) + 90*G4p*G5p*G5px*pow(H,2) - 12*G4*G4p*G5pxx*pow(H,2) + 30*G3x*G4p*G5x*pow(H,2) - 54*G4p*G4px*G5x*pow(H,2) - 14*G2x*G4*G5xx*pow(H,2) + 
           28*G3p*G4*G5xx*pow(H,2) - 54*pow(G4p,2)*G5xx*pow(H,2) - 108*G4*G4xx*G5x*pow(H,4) + 60*G4*G5px*G5x*pow(H,4) + 57*G4p*pow(G5x,2)*pow(H,4) + 76*G4*G4x*G5xx*pow(H,4) - 
           76*G4*G5p*G5xx*pow(H,4) + 2*G2xx*(G3x*G4 - 8*G4p*G4x + 4*G4p*G5p + 3*G4*G5x*pow(H,2)) - 2*G3px*(2*G3x*G4 - 2*G4*G4px - 10*G4p*G4x + 5*G4p*G5p + 6*G4*G5x*pow(H,2))) + 
        ddphi*dphi*(4*dH*H*(3*G3xx*G4*(-2*G4x + G5p) + 2*pow(G4,2)*G5xxx*pow(H,2) + 3*G4p*(-16*G4x*G4xx + 12*G4x*G5px - 2*G5p*G5px + G3x*G5x - 2*G4p*G5xx + 3*pow(G5x,2)*pow(H,2)) + 
              G4*(12*G4pxx*G4x + 60*G3x*G4xx - 120*G4px*G4xx - 6*G4pxx*G5p - 36*G3x*G5px + 72*G4px*G5px - 7*G2xx*G5x + 7*G3px*G5x - 2*G2x*G5xx + 4*G3p*G5xx + 60*G4xx*G5x*pow(H,2) - 
                 27*G5px*G5x*pow(H,2) + 54*G4x*G5xx*pow(H,2) - 61*G5p*G5xx*pow(H,2))) + 
           dphi*(2*pow(G2xx,2)*G4 + 6*pow(G3px,2)*G4 + 4*G2pxx*G3x*G4 - 4*G3ppx*G3x*G4 - 12*G3x*G4p*G4ppx - 8*G2pxx*G4*G4px + 8*G3ppx*G4*G4px + 24*G4p*G4ppx*G4px - 8*G2pxx*G4p*G4x + 8*G3ppx*G4p*G4x + 
              4*G2pxx*G4p*G5p - 4*G3ppx*G4p*G5p + 90*G3x*G3xx*G4*pow(H,2) - 372*G3xx*G4*G4px*pow(H,2) - 156*G3x*G4*G4pxx*pow(H,2) + 792*G4*G4px*G4pxx*pow(H,2) - 16*G2xxx*G4*G4x*pow(H,2) + 
              64*G3pxx*G4*G4x*pow(H,2) + 84*G3xx*G4p*G4x*pow(H,2) - 96*G4*G4ppxx*G4x*pow(H,2) - 312*G4p*G4pxx*G4x*pow(H,2) - 438*G3x*G4p*G4xx*pow(H,2) - 96*G4*G4ppx*G4xx*pow(H,2) + 
              1140*G4p*G4px*G4xx*pow(H,2) + 52*G2x*G4*G4xxx*pow(H,2) - 104*G3p*G4*G4xxx*pow(H,2) + 84*pow(G4p,2)*G4xxx*pow(H,2) + 72*G4*G4pp*G4xxx*pow(H,2) + 16*G2xxx*G4*G5p*pow(H,2) - 
              64*G3pxx*G4*G5p*pow(H,2) - 12*G3xx*G4p*G5p*pow(H,2) + 96*G4*G4ppxx*G5p*pow(H,2) + 48*G4p*G4pxx*G5p*pow(H,2) + 60*G3xx*G4*G5pp*pow(H,2) - 120*G4*G4pxx*G5pp*pow(H,2) - 
              120*G4p*G4xx*G5pp*pow(H,2) - 12*G3x*G4*G5ppx*pow(H,2) - 24*G4*G4px*G5ppx*pow(H,2) + 72*G4p*G4x*G5ppx*pow(H,2) - 12*G4p*G5p*G5ppx*pow(H,2) + 261*G3x*G4p*G5px*pow(H,2) + 
              24*G4*G4ppx*G5px*pow(H,2) - 654*G4p*G4px*G5px*pow(H,2) + 60*G4p*G5pp*G5px*pow(H,2) - 30*G2x*G4*G5pxx*pow(H,2) + 60*G3p*G4*G5pxx*pow(H,2) - 54*pow(G4p,2)*G5pxx*pow(H,2) - 
              36*G4*G4pp*G5pxx*pow(H,2) + 12*G2pxx*G4*G5x*pow(H,2) - 12*G3ppx*G4*G5x*pow(H,2) + 36*G4p*G4ppx*G5x*pow(H,2) + 4*G2px*G4*G5xx*pow(H,2) - 8*G3pp*G4*G5xx*pow(H,2) - 
              6*G2x*G4p*G5xx*pow(H,2) + 12*G3p*G4p*G5xx*pow(H,2) - 18*G4p*G4pp*G5xx*pow(H,2) + 2496*G4*pow(G4xx,2)*pow(H,4) - 264*G4*G4x*G4xxx*pow(H,4) + 264*G4*G4xxx*G5p*pow(H,4) - 
              3348*G4*G4xx*G5px*pow(H,4) + 1134*G4*pow(G5px,2)*pow(H,4) + 268*G4*G4x*G5pxx*pow(H,4) - 268*G4*G5p*G5pxx*pow(H,4) + 42*G3xx*G4*G5x*pow(H,4) + 180*G4*G4pxx*G5x*pow(H,4) - 
              786*G4p*G4xx*G5x*pow(H,4) - 132*G4*G5ppx*G5x*pow(H,4) + 447*G4p*G5px*G5x*pow(H,4) + 366*G3x*G4*G5xx*pow(H,4) - 1280*G4*G4px*G5xx*pow(H,4) - 148*G4p*G4x*G5xx*pow(H,4) + 
              284*G4p*G5p*G5xx*pow(H,4) + 116*G4*G5pp*G5xx*pow(H,4) - 20*G4*G4p*G5xxx*pow(H,4) + 318*G4*G5x*G5xx*pow(H,6) + 
              G3px*(9*G3x*G4p - 8*G4*G4ppx - 18*G4p*G4px - 124*G4*G4xx*pow(H,2) + 108*G4*G5px*pow(H,2) - 17*G4p*G5x*pow(H,2)) + 
              G2xx*(-8*G3px*G4 - 3*G3x*G4p + 8*G4*G4ppx + 6*G4p*G4px + 140*G4*G4xx*pow(H,2) - 104*G4*G5px*pow(H,2) + 7*G4p*G5x*pow(H,2)) + 24*pow(G4xx,2)*pow(H,2)*rptot - 
              12*G4x*G4xxx*pow(H,2)*rptot + 6*G4xxx*G5p*pow(H,2)*rptot - 24*G4xx*G5px*pow(H,2)*rptot + 6*pow(G5px,2)*pow(H,2)*rptot + 6*G4x*G5pxx*pow(H,2)*rptot - 
              3*G5p*G5pxx*pow(H,2)*rptot - 3*G3xx*G5x*pow(H,2)*rptot + 6*G4pxx*G5x*pow(H,2)*rptot + 3*G3x*G5xx*pow(H,2)*rptot - 6*G4px*G5xx*pow(H,2)*rptot + 2*G5x*G5xx*pow(H,4)*rptot)) + 
        pow(dphi,2)*(2*dH*dphi*H*(12*G4*G4px*G4pxx - 84*G4p*G4pxx*G4x + 20*G2xx*G4*G4xx + 4*G3px*G4*G4xx - 48*G4*G4ppx*G4xx + 168*G4p*G4px*G4xx + 42*G4p*G4pxx*G5p - 24*G4*G4pxx*G5pp - 
              24*G4p*G4xx*G5pp + 24*G4*G4px*G5ppx + 24*G4p*G4x*G5ppx - 12*G4p*G5p*G5ppx - 14*G2xx*G4*G5px + 2*G3px*G4*G5px + 24*G4*G4ppx*G5px - 96*G4p*G4px*G5px + 12*G4p*G5pp*G5px + G2xx*G4p*G5x - 
              7*G3px*G4p*G5x + 12*G4p*G4ppx*G5x + 120*G4*pow(G4xx,2)*pow(H,2) + 96*G4*G4x*G4xxx*pow(H,2) - 60*G4*G4xxx*G5p*pow(H,2) - 168*G4*G4xx*G5px*pow(H,2) + 
              66*G4*pow(G5px,2)*pow(H,2) - 24*G4*G4x*G5pxx*pow(H,2) + 6*G4*G5p*G5pxx*pow(H,2) - 54*G4*G4pxx*G5x*pow(H,2) - 174*G4p*G4xx*G5x*pow(H,2) - 36*G4*G5ppx*G5x*pow(H,2) + 
              75*G4p*G5px*G5x*pow(H,2) - 104*G4*G4px*G5xx*pow(H,2) - 4*G4p*G4x*G5xx*pow(H,2) + 44*G4p*G5p*G5xx*pow(H,2) + 10*G4*G5pp*G5xx*pow(H,2) + 69*G4*G5x*G5xx*pow(H,4) + 
              G3xx*(-30*G4*G4px + 18*G4p*G4x - 9*G4p*G5p + 12*G4*G5pp + 63*G4*G5x*pow(H,2)) + 3*G3x*(G3xx*G4 + 6*G4*G4pxx - 20*G4p*G4xx - 4*G4*G5ppx + 12*G4p*G5px + 5*G4*G5xx*pow(H,2))) + 
           2*pow(dH,2)*(3*G3x*(4*G4*G4xx - 2*G4*G5px + G4p*G5x) - 3*G4p*(8*G4x*G4xx - 4*G4xx*G5p - 4*G4x*G5px + 2*G5p*G5px + 2*G4px*G5x + 3*pow(G5x,2)*pow(H,2)) - 
              2*G4*(12*G4px*G4xx - 6*G4px*G5px + G2xx*G5x - G3px*G5x - 30*G4xx*G5x*pow(H,2) + 12*G5px*G5x*pow(H,2) - 18*G4x*G5xx*pow(H,2) + 18*G5p*G5xx*pow(H,2))) + 
           dphi*pow(H,2)*(-(drptot*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 
                   2*G5x*G5xx*pow(H,2))) + dphi*(-14*G2xx*G3xx*G4 + 44*G3px*G3xx*G4 + 12*G3x*G3xx*G4p - 84*G3xx*G4*G4ppx + 36*G3x*G4*G4ppxx - 24*G3xx*G4p*G4px - 72*G4*G4ppxx*G4px + 56*G2xx*G4*G4pxx - 
                 116*G3px*G4*G4pxx - 66*G3x*G4p*G4pxx + 168*G4*G4ppx*G4pxx + 132*G4p*G4px*G4pxx - 72*G4p*G4ppxx*G4x - 16*G2pxx*G4*G4xx + 40*G3ppx*G4*G4xx + 12*G2xx*G4p*G4xx - 72*G3px*G4p*G4xx - 
                 48*G4*G4pppx*G4xx + 168*G4p*G4ppx*G4xx + 12*G2px*G4*G4xxx - 24*G3pp*G4*G4xxx + 24*G4*G4ppp*G4xxx + 36*G4p*G4ppxx*G5p + 12*G3xx*G4*G5ppp - 24*G4*G4pxx*G5ppp - 24*G4p*G4xx*G5ppp - 
                 6*G3x*G4*G5pppx + 12*G4*G4px*G5pppx + 12*G4p*G4x*G5pppx - 6*G4p*G5p*G5pppx - 16*G2xx*G4*G5ppx + 16*G3px*G4*G5ppx + 24*G3x*G4p*G5ppx - 48*G4p*G4px*G5ppx + 8*G2pxx*G4*G5px - 
                 20*G3ppx*G4*G5px - 6*G2xx*G4p*G5px + 36*G3px*G4p*G5px + 24*G4*G4pppx*G5px - 84*G4p*G4ppx*G5px + 12*G4p*G5ppp*G5px - 6*G2px*G4*G5pxx + 12*G3pp*G4*G5pxx - 12*G4*G4ppp*G5pxx + 
                 4*G2pxx*G4p*G5x - 10*G3ppx*G4p*G5x + 12*G4p*G4pppx*G5x + 2*G2ppx*G4*G5xx - 2*G3ppp*G4*G5xx - 3*G2px*G4p*G5xx + 6*G3pp*G4p*G5xx - 6*G4p*G4ppp*G5xx + 24*G3xxx*G4*G4x*pow(H,2) - 
                 144*G4*G4pxxx*G4x*pow(H,2) - 360*G3xx*G4*G4xx*pow(H,2) + 1104*G4*G4pxx*G4xx*pow(H,2) + 480*G4p*pow(G4xx,2)*pow(H,2) - 144*G3x*G4*G4xxx*pow(H,2) + 
                 672*G4*G4px*G4xxx*pow(H,2) - 24*G4p*G4x*G4xxx*pow(H,2) - 24*G3xxx*G4*G5p*pow(H,2) + 144*G4*G4pxxx*G5p*pow(H,2) - 48*G4p*G4xxx*G5p*pow(H,2) - 168*G4*G4xxx*G5pp*pow(H,2) - 
                 200*G4*G4xx*G5ppx*pow(H,2) + 56*G4*G4x*G5ppxx*pow(H,2) - 56*G4*G5p*G5ppxx*pow(H,2) + 288*G3xx*G4*G5px*pow(H,2) - 852*G4*G4pxx*G5px*pow(H,2) - 648*G4p*G4xx*G5px*pow(H,2) + 
                 148*G4*G5ppx*G5px*pow(H,2) + 204*G4p*pow(G5px,2)*pow(H,2) + 80*G3x*G4*G5pxx*pow(H,2) - 388*G4*G4px*G5pxx*pow(H,2) + 32*G4p*G4x*G5pxx*pow(H,2) + 32*G4p*G5p*G5pxx*pow(H,2) + 
                 84*G4*G5pp*G5pxx*pow(H,2) + 18*G3xx*G4p*G5x*pow(H,2) + 108*G4*G4ppxx*G5x*pow(H,2) - 90*G4p*G4pxx*G5x*pow(H,2) - 18*G4*G5pppx*G5x*pow(H,2) + 26*G4p*G5ppx*G5x*pow(H,2) - 
                 40*G2xx*G4*G5xx*pow(H,2) + 102*G3px*G4*G5xx*pow(H,2) + 66*G3x*G4p*G5xx*pow(H,2) - 160*G4*G4ppx*G5xx*pow(H,2) - 198*G4p*G4px*G5xx*pow(H,2) + 42*G4p*G5pp*G5xx*pow(H,2) + 
                 22*G4*G5ppp*G5xx*pow(H,2) - 6*G2x*G4*G5xxx*pow(H,2) + 12*G3p*G4*G5xxx*pow(H,2) - 12*pow(G4p,2)*G5xxx*pow(H,2) - 12*G4*G4pp*G5xxx*pow(H,2) - 48*G4*G4xxx*G5x*pow(H,4) - 
                 24*G4*G5pxx*G5x*pow(H,4) - 1116*G4*G4xx*G5xx*pow(H,4) + 826*G4*G5px*G5xx*pow(H,4) + 136*G4p*G5x*G5xx*pow(H,4) + 36*G4*G4x*G5xxx*pow(H,4) - 36*G4*G5p*G5xxx*pow(H,4) + 
                 2*G2xxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) - 14*G3pxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) + 24*G4xxx*G5x*pow(H,2)*rptot - 
                 12*G5pxx*G5x*pow(H,2)*rptot - 48*G4xx*G5xx*pow(H,2)*rptot + 24*G5px*G5xx*pow(H,2)*rptot + 8*G4x*G5xxx*pow(H,2)*rptot - 4*G5p*G5xxx*pow(H,2)*rptot)))) - 
     2*pow(a,6)*pow(dphi,3)*(-2*pow(ddphi,2)*H*(4*G2xxx*G4*G4x - 4*G3pxx*G4*G4x + 72*G4p*G4pxx*G4x - 4*G2xx*G4*G4xx - 108*G4p*G4px*G4xx - 4*G2x*G4*G4xxx + 8*G3p*G4*G4xxx - 12*pow(G4p,2)*G4xxx - 
           4*G2xxx*G4*G5p + 4*G3pxx*G4*G5p - 30*G4p*G4pxx*G5p + 4*G2xx*G4*G5px - 2*G3px*G4*G5px + 60*G4p*G4px*G5px + 2*G2x*G4*G5pxx - 4*G3p*G4*G5pxx + 6*pow(G4p,2)*G5pxx - 4*G2xx*G4p*G5x + 
           5*G3px*G4p*G5x - 168*G4*pow(G4xx,2)*pow(H,2) + 144*G4*G4x*G4xxx*pow(H,2) - 144*G4*G4xxx*G5p*pow(H,2) + 228*G4*G4xx*G5px*pow(H,2) - 78*G4*pow(G5px,2)*pow(H,2) - 
           84*G4*G4x*G5pxx*pow(H,2) + 84*G4*G5p*G5pxx*pow(H,2) - 60*G4*G4pxx*G5x*pow(H,2) + 66*G4p*G4xx*G5x*pow(H,2) - 30*G4p*G5px*G5x*pow(H,2) + 96*G4*G4px*G5xx*pow(H,2) - 
           8*G4p*G4x*G5xx*pow(H,2) - 26*G4p*G5p*G5xx*pow(H,2) + 2*G4*G4p*G5xxx*pow(H,2) + 3*G4*G5x*G5xx*pow(H,4) + 3*G3xx*(2*G4*G4px - 10*G4p*G4x + 4*G4p*G5p + 7*G4*G5x*pow(H,2)) + 
           G3x*(3*G3xx*G4 - 12*G4*G4pxx + 48*G4p*G4xx - 27*G4p*G5px - 23*G4*G5xx*pow(H,2))) + 
        ddphi*dphi*(-4*dH*(12*G4*G4px*G4pxx + 12*G4p*G4pxx*G4x - 4*G2xx*G4*G4xx + 4*G3px*G4*G4xx - 12*G4p*G4px*G4xx - 6*G4p*G4pxx*G5p + 2*G2xx*G4*G5px - 2*G3px*G4*G5px + 6*G4p*G4px*G5px + 
              120*G4*pow(G4xx,2)*pow(H,2) + 12*G4*G4x*G4xxx*pow(H,2) - 18*G4*G4xxx*G5p*pow(H,2) - 132*G4*G4xx*G5px*pow(H,2) + 36*G4*pow(G5px,2)*pow(H,2) - 6*G4*G4x*G5pxx*pow(H,2) + 
              9*G4*G5p*G5pxx*pow(H,2) + 24*G4*G4pxx*G5x*pow(H,2) + 6*G4p*G4xx*G5x*pow(H,2) + 3*G4p*G5px*G5x*pow(H,2) - 30*G4*G4px*G5xx*pow(H,2) - 24*G4p*G4x*G5xx*pow(H,2) + 
              3*G4p*G5p*G5xx*pow(H,2) + 29*G4*G5x*G5xx*pow(H,4) - 3*G3xx*(2*G4*G4px + 2*G4p*G4x - G4p*G5p + 4*G4*G5x*pow(H,2)) + 
              3*G3x*(G3xx*G4 - 2*G4*G4pxx + 2*G4p*G4xx - G4p*G5px + 6*G4*G5xx*pow(H,2))) + 
           dphi*H*(-16*G2xx*G3xx*G4 + 28*G3px*G3xx*G4 + 15*G3x*G3xx*G4p - 24*G3xx*G4*G4ppx + 24*G3x*G4*G4ppxx - 30*G3xx*G4p*G4px - 48*G4*G4ppxx*G4px + 48*G2xx*G4*G4pxx - 72*G3px*G4*G4pxx - 
              54*G3x*G4p*G4pxx + 48*G4*G4ppx*G4pxx + 108*G4p*G4px*G4pxx - 48*G4p*G4ppxx*G4x - 16*G2pxx*G4*G4xx + 16*G3ppx*G4*G4xx + 12*G2xx*G4p*G4xx - 36*G3px*G4p*G4xx + 48*G4p*G4ppx*G4xx + 
              24*G4p*G4ppxx*G5p - 8*G2xx*G4*G5ppx + 8*G3px*G4*G5ppx + 12*G3x*G4p*G5ppx - 24*G4p*G4px*G5ppx + 8*G2pxx*G4*G5px - 8*G3ppx*G4*G5px - 6*G2xx*G4p*G5px + 18*G3px*G4p*G5px - 24*G4p*G4ppx*G5px + 
              4*G2pxx*G4p*G5x - 4*G3ppx*G4p*G5x + 48*G3xxx*G4*G4x*pow(H,2) - 192*G4*G4pxxx*G4x*pow(H,2) - 420*G3xx*G4*G4xx*pow(H,2) + 744*G4*G4pxx*G4xx*pow(H,2) + 
              936*G4p*pow(G4xx,2)*pow(H,2) - 156*G3x*G4*G4xxx*pow(H,2) + 744*G4*G4px*G4xxx*pow(H,2) - 264*G4p*G4x*G4xxx*pow(H,2) - 48*G3xxx*G4*G5p*pow(H,2) + 192*G4*G4pxxx*G5p*pow(H,2) + 
              48*G4p*G4xxx*G5p*pow(H,2) - 120*G4*G4xxx*G5pp*pow(H,2) + 48*G4*G4xx*G5ppx*pow(H,2) + 48*G4*G4x*G5ppxx*pow(H,2) - 48*G4*G5p*G5ppxx*pow(H,2) + 324*G3xx*G4*G5px*pow(H,2) - 
              648*G4*G4pxx*G5px*pow(H,2) - 1104*G4p*G4xx*G5px*pow(H,2) + 318*G4p*pow(G5px,2)*pow(H,2) + 62*G3x*G4*G5pxx*pow(H,2) - 364*G4*G4px*G5pxx*pow(H,2) + 188*G4p*G4x*G5pxx*pow(H,2) - 
              40*G4p*G5p*G5pxx*pow(H,2) + 60*G4*G5pp*G5pxx*pow(H,2) - 27*G3xx*G4p*G5x*pow(H,2) + 72*G4*G4ppxx*G5x*pow(H,2) + 102*G4p*G4pxx*G5x*pow(H,2) - 24*G4p*G5ppx*G5x*pow(H,2) - 
              50*G2xx*G4*G5xx*pow(H,2) + 50*G3px*G4*G5xx*pow(H,2) + 129*G3x*G4p*G5xx*pow(H,2) + 16*G4*G4ppx*G5xx*pow(H,2) - 324*G4p*G4px*G5xx*pow(H,2) + 30*G4p*G5pp*G5xx*pow(H,2) - 
              10*G2x*G4*G5xxx*pow(H,2) + 20*G3p*G4*G5xxx*pow(H,2) - 18*pow(G4p,2)*G5xxx*pow(H,2) - 12*G4*G4pp*G5xxx*pow(H,2) + 84*G4*G4xxx*G5x*pow(H,4) - 138*G4*G5pxx*G5x*pow(H,4) - 
              1484*G4*G4xx*G5xx*pow(H,4) + 1022*G4*G5px*G5xx*pow(H,4) + 197*G4p*G5x*G5xx*pow(H,4) + 100*G4*G4x*G5xxx*pow(H,4) - 100*G4*G5p*G5xxx*pow(H,4) + 
              4*G2xxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) - 16*G3pxx*(G3x*G4 - 2*G4*G4px - 2*G4p*G4x + G4p*G5p + 3*G4*G5x*pow(H,2)) + 6*G4xxx*G5x*pow(H,2)*rptot - 
              3*G5pxx*G5x*pow(H,2)*rptot - 12*G4xx*G5xx*pow(H,2)*rptot + 6*G5px*G5xx*pow(H,2)*rptot + 2*G4x*G5xxx*pow(H,2)*rptot - G5p*G5xxx*pow(H,2)*rptot)) + 
        pow(dphi,2)*H*(-4*pow(dH,2)*(3*G4p*(-2*G4x + G5p)*G5xx + G4*(24*pow(G4xx,2) - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G3xx*G5x + 6*G4pxx*G5x + 3*G3x*G5xx - 6*G4px*G5xx + 
                 11*G5x*G5xx*pow(H,2))) - 2*dH*dphi*H*(G3xx*(6*G4*(6*G4xx - 5*G5px) + 9*G4p*G5x) - 
              G4p*(144*pow(G4xx,2) - 60*G4x*G4xxx + 30*G4xxx*G5p - 168*G4xx*G5px + 48*pow(G5px,2) + 42*G4x*G5pxx - 21*G5p*G5pxx + 42*G4pxx*G5x - 12*G5ppx*G5x + 18*G3x*G5xx - 48*G4px*G5xx + 
                 6*G5pp*G5xx + 32*G5x*G5xx*pow(H,2)) + G4*(-6*G3x*G4xxx - 36*G4px*G4xxx + 24*G4xxx*G5pp - 48*G4xx*G5ppx + 24*G5ppx*G5px + 12*G4pxx*(2*G4xx + G5px) + 9*G3x*G5pxx + 6*G4px*G5pxx - 
                 12*G5pp*G5pxx + 7*G2xx*G5xx - G3px*G5xx - 12*G4ppx*G5xx + 90*G4xxx*G5x*pow(H,2) - 27*G5pxx*G5x*pow(H,2) + 50*G4xx*G5xx*pow(H,2) - 47*G5px*G5xx*pow(H,2) + 8*G4x*G5xxx*pow(H,2) - 
                 2*G5p*G5xxx*pow(H,2))) + dphi*pow(H,2)*(drptot*(-6*G4xxx*G5x + 3*G5pxx*G5x + 12*G4xx*G5xx - 6*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx) + 
              dphi*(24*pow(G3xx,2)*G4 + 12*G3xxx*G4*G4px + 264*G4*pow(G4pxx,2) - 72*G4*G4px*G4pxxx + 12*G3xxx*G4p*G4x - 72*G4p*G4pxxx*G4x - 8*G2xxx*G4*G4xx + 56*G3pxx*G4*G4xx - 144*G4*G4ppxx*G4xx + 
                 264*G4p*G4pxx*G4xx + 32*G2xx*G4*G4xxx - 92*G3px*G4*G4xxx + 168*G4*G4ppx*G4xxx + 60*G4p*G4px*G4xxx - 6*G3xxx*G4p*G5p + 36*G4p*G4pxxx*G5p - 24*G4*G4xxx*G5ppp + 24*G4*G4xx*G5pppx - 
                 96*G4*G4pxx*G5ppx - 96*G4p*G4xx*G5ppx + 28*G4*G4px*G5ppxx + 28*G4p*G4x*G5ppxx - 14*G4p*G5p*G5ppxx + 4*G2xxx*G4*G5px - 28*G3pxx*G4*G5px + 72*G4*G4ppxx*G5px - 132*G4p*G4pxx*G5px - 
                 12*G4*G5pppx*G5px + 48*G4p*G5ppx*G5px + 24*G3xx*G4p*(-2*G4xx + G5px) - 22*G2xx*G4*G5pxx + 52*G3px*G4*G5pxx - 84*G4*G4ppx*G5pxx - 48*G4p*G4px*G5pxx + 12*G4*G5ppp*G5pxx + 2*G2xxx*G4p*G5x - 
                 14*G3pxx*G4p*G5x + 36*G4p*G4ppxx*G5x - 6*G4p*G5pppx*G5x + 4*G2pxx*G4*G5xx - 10*G3ppx*G4*G5xx - 3*G2xx*G4p*G5xx + 18*G3px*G4p*G5xx + 12*G4*G4pppx*G5xx - 42*G4p*G4ppx*G5xx + 
                 6*G4p*G5ppp*G5xx - 2*G2px*G4*G5xxx + 4*G3pp*G4*G5xxx - 4*G4*G4ppp*G5xxx + 672*G4*G4xx*G4xxx*pow(H,2) - 48*G4*G4x*G4xxxx*pow(H,2) + 48*G4*G4xxxx*G5p*pow(H,2) - 
                 564*G4*G4xxx*G5px*pow(H,2) - 368*G4*G4xx*G5pxx*pow(H,2) + 316*G4*G5px*G5pxx*pow(H,2) + 40*G4*G4x*G5pxxx*pow(H,2) - 40*G4*G5p*G5pxxx*pow(H,2) - 18*G3xxx*G4*G5x*pow(H,2) + 
                 108*G4*G4pxxx*G5x*pow(H,2) - 18*G4p*G4xxx*G5x*pow(H,2) - 42*G4*G5ppxx*G5x*pow(H,2) + 8*G4p*G5pxx*G5x*pow(H,2) - 368*G4*G4pxx*G5xx*pow(H,2) - 288*G4p*G4xx*G5xx*pow(H,2) + 
                 66*G4*G5ppx*G5xx*pow(H,2) + 186*G4p*G5px*G5xx*pow(H,2) - 108*G4*G4px*G5xxx*pow(H,2) + 16*G4p*G4x*G5xxx*pow(H,2) + 4*G4p*G5p*G5xxx*pow(H,2) + 28*G4*G5pp*G5xxx*pow(H,2) + 
                 170*G4*pow(G5xx,2)*pow(H,4) - 12*G4*G5x*G5xxx*pow(H,4) + G3xx*G4*(-180*G4pxx + 48*G5ppx + 122*G5xx*pow(H,2)) + 
                 G3x*(-6*G3xxx*G4 + 6*G4p*(-5*G4xxx + 4*G5pxx) + 2*G4*(18*G4pxxx - 7*G5ppxx + 10*G5xxx*pow(H,2))) + 6*pow(G5xx,2)*pow(H,2)*rptot - 4*G5x*G5xxx*pow(H,2)*rptot)))) + 
     2*pow(a,11)*(dphi*(drptot*(3*G3x*G4p - 6*G4p*G4px + 2*G2x*G4x - 4*G3p*G4x - G2x*G5p + 2*G3p*G5p - 12*pow(G4x,2)*pow(H,2) + 30*G4x*G5p*pow(H,2) - 18*pow(G5p,2)*pow(H,2) + 
              9*G4p*G5x*pow(H,2)) + 4*dH*H*(-6*G2xx*pow(G4,2) + 6*G3px*pow(G4,2) + 9*G3x*G4*G4p - 18*G4*G4p*G4px + 22*G2x*G4*G4x - 44*G3p*G4*G4x + 54*pow(G4p,2)*G4x - 12*G4*G4pp*G4x - 
              19*G2x*G4*G5p + 38*G3p*G4*G5p - 45*pow(G4p,2)*G5p + 12*G4*G4pp*G5p + 60*G4*pow(G4x,2)*pow(H,2) - 144*pow(G4,2)*G4xx*pow(H,2) - 102*G4*G4x*G5p*pow(H,2) + 
              42*G4*pow(G5p,2)*pow(H,2) + 90*pow(G4,2)*G5px*pow(H,2) + 27*G4*G4p*G5x*pow(H,2))) - 
        2*pow(dphi,2)*(-4*G3p*G3pp*G4 + G2pp*G3x*G4 + 6*G3pp*pow(G4p,2) + 4*G3p*G4*G4ppp - 6*pow(G4p,2)*G4ppp - 2*G2pp*G4*G4px - 2*G2pp*G4p*G4x + G2pp*G4p*G5p - 32*G3p*G3x*G4*pow(H,2) + 
           24*G3x*pow(G4p,2)*pow(H,2) + 48*G3x*G4*G4pp*pow(H,2) + 120*G3p*G4*G4px*pow(H,2) - 114*pow(G4p,2)*G4px*pow(H,2) - 132*G4*G4pp*G4px*pow(H,2) + 20*G3pp*G4*G4x*pow(H,2) + 
           20*G3p*G4p*G4x*pow(H,2) - 60*G4p*G4pp*G4x*pow(H,2) - 12*G4*G4ppp*G4x*pow(H,2) - 20*G3pp*G4*G5p*pow(H,2) - 16*G3p*G4p*G5p*pow(H,2) + 48*G4p*G4pp*G5p*pow(H,2) + 
           12*G4*G4ppp*G5p*pow(H,2) - 28*G3p*G4*G5pp*pow(H,2) + 42*pow(G4p,2)*G5pp*pow(H,2) + 3*G2pp*G4*G5x*pow(H,2) + 144*G3x*G4*G4x*pow(H,4) - 528*G4*G4px*G4x*pow(H,4) - 
           84*G4p*pow(G4x,2)*pow(H,4) - 144*G3x*G4*G5p*pow(H,4) + 528*G4*G4px*G5p*pow(H,4) + 156*G4p*G4x*G5p*pow(H,4) - 72*G4p*pow(G5p,2)*pow(H,4) + 108*G4*G4x*G5pp*pow(H,4) - 
           108*G4*G5p*G5pp*pow(H,4) - 48*G3p*G4*G5x*pow(H,4) + 36*pow(G4p,2)*G5x*pow(H,4) + 72*G4*G4pp*G5x*pow(H,4) + 192*G4*G4x*G5x*pow(H,6) - 192*G4*G5p*G5x*pow(H,6) - 
           G2px*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 10*G4*G4x*pow(H,2) - 10*G4*G5p*pow(H,2)) + 9*G3xx*G4*pow(H,2)*rptot - 18*G4*G4pxx*pow(H,2)*rptot + 24*G4px*G4x*pow(H,2)*rptot - 
           42*G4p*G4xx*pow(H,2)*rptot - 12*G3x*G5p*pow(H,2)*rptot + 12*G4px*G5p*pow(H,2)*rptot + 21*G4p*G5px*pow(H,2)*rptot + 4*G3p*G5x*pow(H,2)*rptot + 36*G4x*G5x*pow(H,4)*rptot - 
           48*G5p*G5x*pow(H,4)*rptot + 21*G4*G5xx*pow(H,4)*rptot + 2*G2x*(G3pp*G4 - G4*G4ppp + G4*pow(H,2)*(8*G3x - 30*G4px + 7*G5pp + 12*G5x*pow(H,2)) - 
              pow(H,2)*(5*G4p*G4x - 4*G4p*G5p + G5x*rptot))) + ddphi*(4*pow(G2x,2)*G4 + 
           G2x*(-16*G3p*G4 + 12*pow(G4p,2) + 12*G4*G4pp + 80*G4*G4x*pow(H,2) - 80*G4*G5p*pow(H,2) - 2*G4x*rptot + G5p*rptot) + 
           2*(8*pow(G3p,2)*G4 + 36*G4*G4pp*G4x*pow(H,2) - 36*G4*G4pp*G5p*pow(H,2) + 168*G4*pow(G4x,2)*pow(H,4) - 336*G4*G4x*G5p*pow(H,4) + 168*G4*pow(G5p,2)*pow(H,4) + 
              6*pow(G4p,2)*(3*G4pp + 28*(G4x - G5p)*pow(H,2)) + G2xx*G4*rptot - G3px*G4*rptot + 6*pow(G4x,2)*pow(H,2)*rptot + 24*G4*G4xx*pow(H,2)*rptot - 15*G4x*G5p*pow(H,2)*rptot + 
              9*pow(G5p,2)*pow(H,2)*rptot - 15*G4*G5px*pow(H,2)*rptot - G3p*(12*pow(G4p,2) + 12*G4*G4pp + 80*G4*G4x*pow(H,2) - 80*G4*G5p*pow(H,2) - 2*G4x*rptot + G5p*rptot) - 
              3*G4p*(-60*G4*G4px*pow(H,2) + 20*G4*G5x*pow(H,4) - 2*G4px*rptot + 3*G5x*pow(H,2)*rptot + G3x*(20*G4*pow(H,2) + rptot))))) + 
     pow(a,9)*(8*ddphi*dH*dphi*H*(6*G3xx*pow(G4,2) - 6*G4p*(2*G4x*G5p - 2*pow(G5p,2) + 3*G4p*G5x) + 
           G4*(12*G3x*G4x - 12*G4px*G4x - 24*G4p*G4xx - 18*G3x*G5p + 30*G4px*G5p + 12*G4p*G5px - 7*G2x*G5x + 14*G3p*G5x + 66*G4x*G5x*pow(H,2) - 72*G5p*G5x*pow(H,2)) - 
           2*pow(G4,2)*(6*G4pxx - 7*G5xx*pow(H,2))) - 4*pow(ddphi,2)*(G2x*(2*G3x*G4 - 6*G4*G4px + 2*G4p*G4x - G4p*G5p + 6*G4*G5x*pow(H,2)) - 
           2*(G2xx*G4*G4p - G3px*G4*G4p - 6*G3x*pow(G4p,2) + 15*pow(G4p,2)*G4px - 6*G3x*G4*G4x*pow(H,2) + 18*G4*G4px*G4x*pow(H,2) + 18*G4p*pow(G4x,2)*pow(H,2) + 24*G4*G4p*G4xx*pow(H,2) + 
              6*G3x*G4*G5p*pow(H,2) - 18*G4*G4px*G5p*pow(H,2) - 39*G4p*G4x*G5p*pow(H,2) + 21*G4p*pow(G5p,2)*pow(H,2) - 15*G4*G4p*G5px*pow(H,2) - 18*pow(G4p,2)*G5x*pow(H,2) - 
              18*G4*G4x*G5x*pow(H,4) + 18*G4*G5p*G5x*pow(H,4) + G3p*(2*G3x*G4 - 6*G4*G4px + 2*G4p*G4x - G4p*G5p + 6*G4*G5x*pow(H,2)))) + 
        pow(dphi,3)*(-(drptot*(3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
                18*G5p*G5px*pow(H,2) + 6*G3x*G5x*pow(H,2) - 6*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))) + 
           4*dH*H*(12*pow(G3x,2)*G4 + 36*pow(G4,2)*(-2*G4xxx + G5pxx)*pow(H,2) + 
              G4p*(96*G4px*G4x + 72*G4p*G4xx - 132*G4px*G5p + 12*G4x*G5pp + 6*G5p*G5pp - 48*G4p*G5px - 2*G2x*G5x + 4*G3p*G5x + 6*G4pp*G5x - 234*G4x*G5x*pow(H,2) + 249*G5p*G5x*pow(H,2)) + 
              G3x*(-54*G4p*G4x + 57*G4p*G5p + 6*G4*(-15*G4px + G5pp + 21*G5x*pow(H,2))) + 
              2*G4*(90*pow(G4px,2) + 16*G2xx*G4x - 4*G3px*G4x - 24*G4ppx*G4x + 16*G2x*G4xx - 32*G3p*G4xx - 12*G4pp*G4xx - 13*G2xx*G5p + G3px*G5p + 24*G4ppx*G5p - 10*G2x*G5px + 20*G3p*G5px + 
                 6*G4pp*G5px + 228*G4x*G4xx*pow(H,2) - 156*G4xx*G5p*pow(H,2) - 120*G4x*G5px*pow(H,2) + 75*G5p*G5px*pow(H,2) - 15*G5pp*G5x*pow(H,2) + 9*G4p*G5xx*pow(H,2) + 
                 27*pow(G5x,2)*pow(H,4) - 3*G4px*(6*G5pp + 65*G5x*pow(H,2))))) + 
        pow(dphi,2)*(-8*pow(dH,2)*(12*G4*G4px*(G4x - G5p) + 6*G3x*G4*(-G4x + G5p) + 3*G4p*(4*pow(G4x,2) - 6*G4x*G5p + 2*pow(G5p,2) + G4p*G5x) + 
              G4*G5x*(G2x - 2*(G3p + 24*(G4x - G5p)*pow(H,2)))) + ddphi*(40*G3p*G3px*G4 + 8*G2px*G3x*G4 - 16*G3pp*G3x*G4 + 24*G3p*G3x*G4p - 36*G3px*pow(G4p,2) - 24*G3px*G4*G4pp - 36*G3x*G4p*G4pp - 
              32*G3p*G4*G4ppx + 48*pow(G4p,2)*G4ppx - 16*G2px*G4*G4px + 32*G3pp*G4*G4px - 48*G3p*G4p*G4px + 72*G4p*G4pp*G4px - 16*G2px*G4p*G4x + 32*G3pp*G4p*G4x + 8*G2px*G4p*G5p - 16*G3pp*G4p*G5p + 
              312*pow(G3x,2)*G4*pow(H,2) - 120*G3xx*G4*G4p*pow(H,2) - 2040*G3x*G4*G4px*pow(H,2) + 3360*G4*pow(G4px,2)*pow(H,2) + 240*G4*G4p*G4pxx*pow(H,2) + 136*G3px*G4*G4x*pow(H,2) - 
              432*G3x*G4p*G4x*pow(H,2) - 480*G4*G4ppx*G4x*pow(H,2) + 1104*G4p*G4px*G4x*pow(H,2) - 944*G3p*G4*G4xx*pow(H,2) + 1176*pow(G4p,2)*G4xx*pow(H,2) + 576*G4*G4pp*G4xx*pow(H,2) - 
              136*G3px*G4*G5p*pow(H,2) + 528*G3x*G4p*G5p*pow(H,2) + 480*G4*G4ppx*G5p*pow(H,2) - 1440*G4p*G4px*G5p*pow(H,2) + 192*G3x*G4*G5pp*pow(H,2) - 624*G4*G4px*G5pp*pow(H,2) - 
              144*G4p*G4x*G5pp*pow(H,2) + 192*G4p*G5p*G5pp*pow(H,2) + 632*G3p*G4*G5px*pow(H,2) - 756*pow(G4p,2)*G5px*pow(H,2) - 360*G4*G4pp*G5px*pow(H,2) + 24*G2px*G4*G5x*pow(H,2) - 
              48*G3pp*G4*G5x*pow(H,2) + 56*G3p*G4p*G5x*pow(H,2) - 108*G4p*G4pp*G5x*pow(H,2) + 2640*G4*G4x*G4xx*pow(H,4) - 2640*G4*G4xx*G5p*pow(H,4) - 1512*G4*G4x*G5px*pow(H,4) + 
              1512*G4*G5p*G5px*pow(H,4) + 816*G3x*G4*G5x*pow(H,4) - 2520*G4*G4px*G5x*pow(H,4) - 1200*G4p*G4x*G5x*pow(H,4) + 1296*G4p*G5p*G5x*pow(H,4) + 96*G4*G5pp*G5x*pow(H,4) - 
              280*G4*G4p*G5xx*pow(H,4) + 504*G4*pow(G5x,2)*pow(H,6) + 4*G2x*
               (3*G2xx*G4 - 5*G3px*G4 - 3*G3x*G4p + 4*G4*G4ppx + 6*G4p*G4px + 118*G4*G4xx*pow(H,2) - 79*G4*G5px*pow(H,2) - 7*G4p*G5x*pow(H,2)) + 3*pow(G3x,2)*rptot - 12*G3x*G4px*rptot + 
              12*pow(G4px,2)*rptot + 4*G3px*G4x*rptot - 2*G3px*G5p*rptot + 24*G4*G4xxx*pow(H,2)*rptot - 48*G4xx*G5p*pow(H,2)*rptot + 12*G4x*G5px*pow(H,2)*rptot + 18*G5p*G5px*pow(H,2)*rptot - 
              12*G4*G5pxx*pow(H,2)*rptot + 6*G3x*G5x*pow(H,2)*rptot - 12*G4p*G5xx*pow(H,2)*rptot + 15*pow(G5x,2)*pow(H,4)*rptot + 
              G2xx*(-24*G3p*G4 + 12*pow(G4p,2) + 24*G4*G4pp + 40*G4*G4x*pow(H,2) - 40*G4*G5p*pow(H,2) - 4*G4x*rptot + 2*G5p*rptot))) + 
        2*pow(dphi,4)*(G2px*(2*G2xx*G4 - 2*G3px*G4 - 3*G3x*G4p + 6*G4p*G4px + 56*G4*G4xx*pow(H,2) - 34*G4*G5px*pow(H,2) - 11*G4p*G5x*pow(H,2)) - 
           2*(-2*G3pp*G3px*G4 - G2ppx*G3x*G4 + G3ppp*G3x*G4 - 3*G3pp*G3x*G4p + 2*G3px*G4*G4ppp + 3*G3x*G4p*G4ppp + 2*G2ppx*G4*G4px - 2*G3ppp*G4*G4px + 6*G3pp*G4p*G4px - 6*G4p*G4ppp*G4px + 
              2*G2ppx*G4p*G4x - 2*G3ppp*G4p*G4x - G2ppx*G4p*G5p + G3ppp*G4p*G5p - 39*G3px*G3x*G4*pow(H,2) + 7*G2x*G3xx*G4*pow(H,2) - 14*G3p*G3xx*G4*pow(H,2) - 12*pow(G3x,2)*G4p*pow(H,2) + 
              12*G3xx*pow(G4p,2)*pow(H,2) + 18*G3xx*G4*G4pp*pow(H,2) + 66*G3x*G4*G4ppx*pow(H,2) + 136*G3px*G4*G4px*pow(H,2) + 81*G3x*G4p*G4px*pow(H,2) - 216*G4*G4ppx*G4px*pow(H,2) - 
              114*G4p*pow(G4px,2)*pow(H,2) - 28*G2x*G4*G4pxx*pow(H,2) + 56*G3p*G4*G4pxx*pow(H,2) - 66*pow(G4p,2)*G4pxx*pow(H,2) - 36*G4*G4pp*G4pxx*pow(H,2) + 8*G2pxx*G4*G4x*pow(H,2) - 
              20*G3ppx*G4*G4x*pow(H,2) + 26*G3px*G4p*G4x*pow(H,2) + 24*G4*G4pppx*G4x*pow(H,2) - 48*G4p*G4ppx*G4x*pow(H,2) + 56*G3pp*G4*G4xx*pow(H,2) - 6*G2x*G4p*G4xx*pow(H,2) + 
              12*G3p*G4p*G4xx*pow(H,2) - 36*G4p*G4pp*G4xx*pow(H,2) - 48*G4*G4ppp*G4xx*pow(H,2) - 8*G2pxx*G4*G5p*pow(H,2) + 20*G3ppx*G4*G5p*pow(H,2) - 31*G3px*G4p*G5p*pow(H,2) - 
              24*G4*G4pppx*G5p*pow(H,2) + 66*G4p*G4ppx*G5p*pow(H,2) - 14*G3px*G4*G5pp*pow(H,2) - 21*G3x*G4p*G5pp*pow(H,2) + 42*G4p*G4px*G5pp*pow(H,2) - 9*G3x*G4*G5ppp*pow(H,2) + 
              30*G4*G4px*G5ppp*pow(H,2) + 6*G4p*G4x*G5ppp*pow(H,2) - 9*G4p*G5p*G5ppp*pow(H,2) + 8*G2x*G4*G5ppx*pow(H,2) - 16*G3p*G4*G5ppx*pow(H,2) + 24*pow(G4p,2)*G5ppx*pow(H,2) - 
              34*G3pp*G4*G5px*pow(H,2) + 3*G2x*G4p*G5px*pow(H,2) - 6*G3p*G4p*G5px*pow(H,2) + 18*G4p*G4pp*G5px*pow(H,2) + 30*G4*G4ppp*G5px*pow(H,2) - 3*G2ppx*G4*G5x*pow(H,2) + 
              3*G3ppp*G4*G5x*pow(H,2) - 11*G3pp*G4p*G5x*pow(H,2) + 9*G4p*G4ppp*G5x*pow(H,2) + G2pp*G4*G5xx*pow(H,2) + 18*G3xx*G4*G4x*pow(H,4) + 24*G4*G4pxx*G4x*pow(H,4) + 
              450*G3x*G4*G4xx*pow(H,4) - 1716*G4*G4px*G4xx*pow(H,4) - 288*G4p*G4x*G4xx*pow(H,4) - 18*G3xx*G4*G5p*pow(H,4) - 24*G4*G4pxx*G5p*pow(H,4) + 294*G4p*G4xx*G5p*pow(H,4) + 
              360*G4*G4xx*G5pp*pow(H,4) - 44*G4*G4x*G5ppx*pow(H,4) + 44*G4*G5p*G5ppx*pow(H,4) - 327*G3x*G4*G5px*pow(H,4) + 1212*G4*G4px*G5px*pow(H,4) + 216*G4p*G4x*G5px*pow(H,4) - 
              225*G4p*G5p*G5px*pow(H,4) - 222*G4*G5pp*G5px*pow(H,4) - 33*G3px*G4*G5x*pow(H,4) - 60*G3x*G4p*G5x*pow(H,4) + 30*G4*G4ppx*G5x*pow(H,4) + 231*G4p*G4px*G5x*pow(H,4) - 
              69*G4p*G5pp*G5x*pow(H,4) - 3*G4*G5ppp*G5x*pow(H,4) + 23*G2x*G4*G5xx*pow(H,4) - 46*G3p*G4*G5xx*pow(H,4) + 42*pow(G4p,2)*G5xx*pow(H,4) + 54*G4*G4pp*G5xx*pow(H,4) + 
              534*G4*G4xx*G5x*pow(H,6) - 369*G4*G5px*G5x*pow(H,6) - 60*G4p*pow(G5x,2)*pow(H,6) + 114*G4*G4x*G5xx*pow(H,6) - 114*G4*G5p*G5xx*pow(H,6) - 12*G3xx*G4x*pow(H,2)*rptot + 
              24*G4pxx*G4x*pow(H,2)*rptot + 24*G3x*G4xx*pow(H,2)*rptot - 48*G4px*G4xx*pow(H,2)*rptot + 6*G3xx*G5p*pow(H,2)*rptot - 12*G4pxx*G5p*pow(H,2)*rptot - 12*G3x*G5px*pow(H,2)*rptot + 
              24*G4px*G5px*pow(H,2)*rptot + 2*G3px*G5x*pow(H,2)*rptot + 24*G4xx*G5x*pow(H,4)*rptot - 6*G5px*G5x*pow(H,4)*rptot - 4*G4x*G5xx*pow(H,4)*rptot - 10*G5p*G5xx*pow(H,4)*rptot + 
              3*G4*G5xxx*pow(H,4)*rptot + G2xx*(2*G3pp*G4 - 2*G4*G4ppp + G4*pow(H,2)*(13*G3x - 54*G4px + 14*G5pp + 15*G5x*pow(H,2)) + pow(H,2)*(-4*G4p*G4x + 5*G4p*G5p - 2*G5x*rptot))))) + 
     2*pow(a,10)*(96*pow(dH,2)*dphi*G4*pow(G4x - G5p,2)*H + 4*dH*(ddphi*(2*G2xx*pow(G4,2) - 2*G3px*pow(G4,2) - 6*G3x*G4*G4p + 12*G4*G4p*G4px - 6*G2x*G4*G4x + 12*G3p*G4*G4x - 
              12*pow(G4p,2)*G4x + 5*G2x*G4*G5p - 10*G3p*G4*G5p + 12*pow(G4p,2)*G5p + 36*G4*pow(G4x,2)*pow(H,2) + 48*pow(G4,2)*G4xx*pow(H,2) - 78*G4*G4x*G5p*pow(H,2) + 
              42*G4*pow(G5p,2)*pow(H,2) - 30*pow(G4,2)*G5px*pow(H,2) - 18*G4*G4p*G5x*pow(H,2)) + 
           pow(dphi,2)*(G2x*(3*G3x*G4 + G4p*(-2*G4x + G5p) + 2*G4*(-5*G4px + G5pp + 9*G5x*pow(H,2))) - 2*G3p*(3*G3x*G4 + G4p*(-2*G4x + G5p) + 2*G4*(-5*G4px + G5pp + 9*G5x*pow(H,2))) - 
              3*(-2*G3x*pow(G4p,2) + G4p*G4pp*(-2*G4x + G5p) + 2*G4p*(14*pow(G4x,2) - 6*G4*G4xx - 27*G4x*G5p + 13*pow(G5p,2) + 3*G4*G5px)*pow(H,2) + 
                 G3x*G4*(G4pp + 2*(-14*G4x + 11*G5p)*pow(H,2)) + pow(G4p,2)*(8*G4px - 2*G5pp - 15*G5x*pow(H,2)) + G4*G4pp*(-2*G4px + 3*G5x*pow(H,2)) + 
                 2*G4*pow(H,2)*(3*G3xx*G4 - 6*G4*G4pxx + 40*G4px*G4x - 31*G4px*G5p + 4*G4x*G5pp - 4*G5p*G5pp - 13*G4x*G5x*pow(H,2) + 10*G5p*G5x*pow(H,2) + 7*G4*G5xx*pow(H,2))))) + 
        dphi*H*(dphi*drptot*(12*G4p*G4xx + 6*G3x*G5p - 6*G4px*(2*G4x + G5p) - 6*G4p*G5px + G2x*G5x - 2*G3p*G5x - 18*G4x*G5x*pow(H,2) + 24*G5p*G5x*pow(H,2)) - 
           pow(dphi,2)*(28*G3p*G3px*G4 - 14*G2px*G3x*G4 + 28*G3pp*G3x*G4 + 6*G3p*G3x*G4p - 36*G3px*pow(G4p,2) - 12*G3px*G4*G4pp - 18*G3x*G4p*G4pp - 24*G3x*G4*G4ppp - 56*G3p*G4*G4ppx + 
              84*pow(G4p,2)*G4ppx + 40*G2px*G4*G4px - 80*G3pp*G4*G4px - 12*G3p*G4p*G4px + 36*G4p*G4pp*G4px + 72*G4*G4ppp*G4px - 8*G2ppx*G4*G4x + 8*G3ppp*G4*G4x + 16*G2px*G4p*G4x - 32*G3pp*G4p*G4x + 
              24*G4p*G4ppp*G4x + 8*G2pp*G4*G4xx + 8*G2ppx*G4*G5p - 8*G3ppp*G4*G5p - 14*G2px*G4p*G5p + 28*G3pp*G4p*G5p - 24*G4p*G4ppp*G5p + 8*G3p*G4*G5ppp - 12*pow(G4p,2)*G5ppp - 4*G2pp*G4*G5px - 
              2*G2pp*G4p*G5x + 108*pow(G3x,2)*G4*pow(H,2) - 816*G3x*G4*G4px*pow(H,2) + 1464*G4*pow(G4px,2)*pow(H,2) - 52*G3px*G4*G4x*pow(H,2) - 150*G3x*G4p*G4x*pow(H,2) + 
              24*G4*G4ppx*G4x*pow(H,2) + 612*G4p*G4px*G4x*pow(H,2) - 288*G3p*G4*G4xx*pow(H,2) + 240*pow(G4p,2)*G4xx*pow(H,2) + 384*G4*G4pp*G4xx*pow(H,2) + 52*G3px*G4*G5p*pow(H,2) + 
              138*G3x*G4p*G5p*pow(H,2) - 24*G4*G4ppx*G5p*pow(H,2) - 564*G4p*G4px*G5p*pow(H,2) + 180*G3x*G4*G5pp*pow(H,2) - 528*G4*G4px*G5pp*pow(H,2) - 192*G4p*G4x*G5pp*pow(H,2) + 
              180*G4p*G5p*G5pp*pow(H,2) + 212*G3p*G4*G5px*pow(H,2) - 204*pow(G4p,2)*G5px*pow(H,2) - 228*G4*G4pp*G5px*pow(H,2) - 18*G2px*G4*G5x*pow(H,2) + 36*G3pp*G4*G5x*pow(H,2) + 
              26*G3p*G4p*G5x*pow(H,2) - 78*G4p*G4pp*G5x*pow(H,2) - 24*G4*G4ppp*G5x*pow(H,2) + 1008*G4*G4x*G4xx*pow(H,4) - 1008*G4*G4xx*G5p*pow(H,4) - 684*G4*G4x*G5px*pow(H,4) + 
              684*G4*G5p*G5px*pow(H,4) + 288*G3x*G4*G5x*pow(H,4) - 1056*G4*G4px*G5x*pow(H,4) - 282*G4p*G4x*G5x*pow(H,4) + 270*G4p*G5p*G5x*pow(H,4) + 204*G4*G5pp*G5x*pow(H,4) + 
              180*G4*pow(G5x,2)*pow(H,6) + G2x*(4*G2xx*G4 - 14*G3px*G4 - 3*G3x*G4p + 28*G4*G4ppx + 6*G4p*G4px - 4*G4*G5ppp + 144*G4*G4xx*pow(H,2) - 106*G4*G5px*pow(H,2) - 13*G4p*G5x*pow(H,2)) + 
              6*pow(G3x,2)*rptot - 24*G3x*G4px*rptot + 24*pow(G4px,2)*rptot + 8*G3px*G4x*rptot - 4*G3px*G5p*rptot + 36*G4*G4xxx*pow(H,2)*rptot - 96*G4xx*G5p*pow(H,2)*rptot + 
              24*G4x*G5px*pow(H,2)*rptot + 36*G5p*G5px*pow(H,2)*rptot - 18*G4*G5pxx*pow(H,2)*rptot + 12*G3x*G5x*pow(H,2)*rptot - 21*G4p*G5xx*pow(H,2)*rptot + 30*pow(G5x,2)*pow(H,4)*rptot + 
              G2xx*(-8*G3p*G4 + 6*pow(G4p,2) + 12*G4*(G4pp + 2*(G4x - G5p)*pow(H,2)) + 4*(-2*G4x + G5p)*rptot)) + 
           ddphi*(G2x*(54*G3x*G4 - 176*G4*G4px - 16*G4p*G4x + 20*G4p*G5p + 20*G4*G5pp + 78*G4*G5x*pow(H,2) - G5x*rptot) + 
              2*(-10*G2xx*G4*G4p + 10*G3px*G4*G4p + 66*G3x*pow(G4p,2) + 36*G3x*G4*G4pp - 198*pow(G4p,2)*G4px - 108*G4*G4pp*G4px + 8*G2px*G4*G4x - 16*G3pp*G4*G4x - 36*G4p*G4pp*G4x - 8*G2px*G4*G5p + 
                 16*G3pp*G4*G5p + 36*G4p*G4pp*G5p + 30*pow(G4p,2)*G5pp + 210*G3x*G4*G4x*pow(H,2) - 624*G4*G4px*G4x*pow(H,2) - 192*G4p*pow(G4x,2)*pow(H,2) - 240*G4*G4p*G4xx*pow(H,2) - 
                 210*G3x*G4*G5p*pow(H,2) + 624*G4*G4px*G5p*pow(H,2) + 396*G4p*G4x*G5p*pow(H,2) - 204*G4p*pow(G5p,2)*pow(H,2) + 12*G4*G4x*G5pp*pow(H,2) - 12*G4*G5p*G5pp*pow(H,2) + 
                 150*G4*G4p*G5px*pow(H,2) + 144*pow(G4p,2)*G5x*pow(H,2) + 36*G4*G4pp*G5x*pow(H,2) + 282*G4*G4x*G5x*pow(H,4) - 282*G4*G5p*G5x*pow(H,4) + 3*G3xx*G4*rptot - 6*G4*G4pxx*rptot + 
                 6*G4px*G4x*rptot - 12*G4p*G4xx*rptot - 3*G3x*G5p*rptot + 3*G4px*G5p*rptot + 6*G4p*G5px*rptot + 9*G4x*G5x*pow(H,2)*rptot - 12*G5p*G5x*pow(H,2)*rptot + 7*G4*G5xx*pow(H,2)*rptot + 
                 G3p*(-54*G3x*G4 + 176*G4*G4px + 16*G4p*G4x - 20*G4p*G5p - 20*G4*G5pp - 78*G4*G5x*pow(H,2) + G5x*rptot)))))))/
 (-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
   4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
   2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
      18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
   2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
      14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
      3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
    (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
      2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + 
   pow(a,4)*pow(dphi,4)*(3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 
      12*G4x*G5px*pow(H,2) + 18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))))/
(pow(a,3)*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H));
   
   // Milestone siik
   *siik = (-2*(2*pow(a,5)*G4p - ddphi*pow(dphi,3)*G5xx*H - pow(a,2)*dphi*(dH*dphi*G5x + 2*(pow(dphi,2)*(-G4xx + G5px) + ddphi*G5x)*H) + a*(ddphi*pow(dphi,2)*(-2*G4xx + G5px) + pow(dphi,4)*G5xx*pow(H,2)) + 
   pow(a,3)*(2*ddphi*(-G4x + G5p) + pow(dphi,2)*(-2*G4px + G5pp + G5x*pow(H,2)))))/(pow(a,2)*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H));

   // Milestone siip
   *siip = (3*(a*pow(dphi,13)*G5xx*(G5pxx*G5xx + 14*G4xx*G5xxx - 8*G5px*G5xxx)*pow(H,6) + pow(dphi,14)*pow(G5xx,2)*G5xxx*pow(H,7) + 
   pow(a,2)*pow(dphi,12)*pow(H,5)*(2*pow(G5pxx,2)*G5x - 3*G3xx*pow(G5xx,2) + 6*G4pxx*pow(G5xx,2) - 4*G4xxx*(G5pxx*G5x + 3*(-2*G4xx + G5px)*G5xx) + 48*pow(G4xx,2)*G5xxx - 56*G4xx*G5px*G5xxx + 
      16*pow(G5px,2)*G5xxx + 4*G4pxx*G5x*G5xxx - 2*G5ppx*G5x*G5xxx + 6*G3x*G5xx*G5xxx - 16*G4px*G5xx*G5xxx + 2*G5pp*G5xx*G5xxx + 4*pow(G5xx,3)*pow(H,2) + 10*G5x*G5xx*G5xxx*pow(H,2)) + 
   pow(a,3)*pow(dphi,11)*pow(H,4)*(8*pow(G4xx,2)*(18*G4xxx - 5*G5pxx) - 16*pow(G5px,2)*G5pxx + 4*G4x*pow(G5pxx,2) - 2*G5p*pow(G5pxx,2) - 4*G3xx*G5pxx*G5x + 8*G4pxx*G5pxx*G5x + 
      12*G3xx*G5px*G5xx - 36*G4pxx*G5px*G5xx + 6*G5ppx*G5px*G5xx - 2*G3x*G5pxx*G5xx + 4*G4px*G5pxx*G5xx - 2*G2xx*pow(G5xx,2) - G3px*pow(G5xx,2) + 6*G4ppx*pow(G5xx,2) + 8*G4pxx*G4x*G5xxx - 
      4*G4pxx*G5p*G5xxx - 4*G4x*G5ppx*G5xxx + 2*G5p*G5ppx*G5xxx - 24*G3x*G5px*G5xxx + 64*G4px*G5px*G5xxx - 8*G5pp*G5px*G5xxx + 4*G3px*G5x*G5xxx - 8*G4ppx*G5x*G5xxx + G2x*G5xx*G5xxx - 2*G3p*G5xx*G5xxx + 
      2*G4pp*G5xx*G5xxx - 27*G5px*pow(G5xx,2)*pow(H,2) - 44*G5px*G5x*G5xxx*pow(H,2) + 22*G4x*G5xx*G5xxx*pow(H,2) - 18*G5p*G5xx*G5xxx*pow(H,2) + 
      4*G4xxx*(12*pow(G5px,2) - 2*G4x*G5pxx + G5p*G5pxx + 3*G3x*G5xx - 6*G4px*G5xx + 3*G5x*G5xx*pow(H,2)) - 
      2*G4xx*(84*G4xxx*G5px - 26*G5px*G5pxx + 9*G3xx*G5xx - 30*G4pxx*G5xx + 6*G5ppx*G5xx - 21*G3x*G5xxx + 58*G4px*G5xxx - 8*G5pp*G5xxx - 28*pow(G5xx,2)*pow(H,2) - 37*G5x*G5xxx*pow(H,2))) + 
   pow(a,4)*pow(dphi,10)*pow(H,3)*(-384*G4px*G4xx*G4xxx + 48*G4xx*G4xxx*G5pp - 48*pow(G4xx,2)*G5ppx + 216*G4px*G4xxx*G5px - 24*G4xxx*G5pp*G5px + 48*G4xx*G5ppx*G5px - 12*G5ppx*pow(G5px,2) - 
      8*G3xx*G4x*G5pxx + 152*G4px*G4xx*G5pxx + 4*G3xx*G5p*G5pxx - 24*G4xx*G5pp*G5pxx - 88*G4px*G5px*G5pxx + 12*G5pp*G5px*G5pxx + 24*pow(G4pxx,2)*G5x + 12*G3px*G4xxx*G5x - 24*G4ppx*G4xxx*G5x + 
      6*G3xx*G5ppx*G5x - 2*G2xx*G5pxx*G5x - 4*G3px*G5pxx*G5x + 12*G4ppx*G5pxx*G5x + 24*G3xx*G4px*G5xx - 16*G2xx*G4xx*G5xx + 4*G3px*G4xx*G5xx + 24*G4ppx*G4xx*G5xx - 6*G3xx*G5pp*G5xx + 12*G4px*G5ppx*G5xx + 
      10*G2xx*G5px*G5xx - 4*G3px*G5px*G5xx - 12*G4ppx*G5px*G5xx - 2*G2px*pow(G5xx,2) + 2*G3pp*pow(G5xx,2) + 9*pow(G3x,2)*G5xxx + 60*pow(G4px,2)*G5xxx + 8*G3px*G4x*G5xxx - 16*G4ppx*G4x*G5xxx + 
      8*G2x*G4xx*G5xxx - 16*G3p*G4xx*G5xxx + 16*G4pp*G4xx*G5xxx - 4*G3px*G5p*G5xxx + 8*G4ppx*G5p*G5xxx - 12*G4px*G5pp*G5xxx - 4*G2x*G5px*G5xxx + 8*G3p*G5px*G5xxx - 8*G4pp*G5px*G5xxx + 2*G2px*G5x*G5xxx - 
      2*G3pp*G5x*G5xxx + 192*G4xx*G4xxx*G5x*pow(H,2) - 128*G4xxx*G5px*G5x*pow(H,2) - 76*G4xx*G5pxx*G5x*pow(H,2) + 56*G5px*G5pxx*G5x*pow(H,2) + 264*pow(G4xx,2)*G5xx*pow(H,2) + 
      48*G4x*G4xxx*G5xx*pow(H,2) - 36*G4xxx*G5p*G5xx*pow(H,2) - 260*G4xx*G5px*G5xx*pow(H,2) + 64*pow(G5px,2)*G5xx*pow(H,2) - 12*G4x*G5pxx*G5xx*pow(H,2) + 6*G5p*G5pxx*G5xx*pow(H,2) - 
      18*G3xx*G5x*G5xx*pow(H,2) - 14*G5ppx*G5x*G5xx*pow(H,2) - 28*G4px*pow(G5xx,2)*pow(H,2) + 2*G5pp*pow(G5xx,2)*pow(H,2) + 160*G4x*G4xx*G5xxx*pow(H,2) - 128*G4xx*G5p*G5xxx*pow(H,2) - 
      96*G4x*G5px*G5xxx*pow(H,2) + 76*G5p*G5px*G5xxx*pow(H,2) - 104*G4px*G5x*G5xxx*pow(H,2) + 12*G5pp*G5x*G5xxx*pow(H,2) - 8*G4p*G5xx*G5xxx*pow(H,2) + 49*G5x*pow(G5xx,2)*pow(H,4) + 
      21*pow(G5x,2)*G5xxx*pow(H,4) + 4*G4pxx*(24*pow(G4xx,2) - 24*G4xx*G5px + 6*pow(G5px,2) + 4*G4x*G5pxx - 2*G5p*G5pxx - 3*G3xx*G5x - 3*G5ppx*G5x + 6*G3x*G5xx - 18*G4px*G5xx + 3*G5pp*G5xx + 
         16*G5x*G5xx*pow(H,2)) + G3x*(144*G4xx*G4xxx - 84*G4xxx*G5px - 52*G4xx*G5pxx + 32*G5px*G5pxx - 6*G3xx*G5xx - 6*G5ppx*G5xx - 48*G4px*G5xxx + 6*G5pp*G5xxx + 15*pow(G5xx,2)*pow(H,2) + 
         38*G5x*G5xxx*pow(H,2))) - 8*pow(a,14)*H*(G2x*G4*G4p - 2*G3p*G4*G4p + 3*pow(G4p,3) - 2*G2p*G4*G4x + 2*G2p*G4*G5p - 18*G4*G4p*G4x*pow(H,2) + 18*G4*G4p*G5p*pow(H,2) + 3*G4p*G4x*rptot - 
      3*G4p*G5p*rptot) + pow(a,5)*pow(dphi,9)*pow(H,2)*(-24*G2xx*pow(G4xx,2) + 24*G3px*pow(G4xx,2) + 36*pow(G3x,2)*G4xxx - 192*G3x*G4px*G4xxx + 240*pow(G4px,2)*G4xxx + 24*G3px*G4x*G4xxx - 
      48*G4ppx*G4x*G4xxx + 24*G2x*G4xx*G4xxx - 48*G3p*G4xx*G4xxx + 48*G4pp*G4xx*G4xxx + 24*pow(G4pxx,2)*(2*G4x - G5p) - 12*G3px*G4xxx*G5p + 24*G4ppx*G4xxx*G5p + 24*G3x*G4xxx*G5pp - 48*G4px*G4xxx*G5pp - 
      36*G3x*G4xx*G5ppx + 72*G4px*G4xx*G5ppx + 28*G2xx*G4xx*G5px - 28*G3px*G4xx*G5px - 12*G2x*G4xxx*G5px + 24*G3p*G4xxx*G5px - 24*G4pp*G4xxx*G5px + 18*G3x*G5ppx*G5px - 36*G4px*G5ppx*G5px - 
      8*G2xx*pow(G5px,2) + 8*G3px*pow(G5px,2) - 15*pow(G3x,2)*G5pxx + 84*G3x*G4px*G5pxx - 108*pow(G4px,2)*G5pxx - 4*G2xx*G4x*G5pxx - 8*G3px*G4x*G5pxx + 24*G4ppx*G4x*G5pxx - 12*G2x*G4xx*G5pxx + 
      24*G3p*G4xx*G5pxx - 24*G4pp*G4xx*G5pxx + 2*G2xx*G5p*G5pxx + 4*G3px*G5p*G5pxx - 12*G4ppx*G5p*G5pxx - 12*G3x*G5pp*G5pxx + 24*G4px*G5pp*G5pxx + 6*G2x*G5px*G5pxx - 12*G3p*G5px*G5pxx + 
      12*G4pp*G5px*G5pxx + 8*G2px*G4xxx*G5x - 8*G3pp*G4xxx*G5x + 4*G2xx*G5ppx*G5x - 4*G3px*G5ppx*G5x - 4*G2px*G5pxx*G5x + 4*G3pp*G5pxx*G5x - 6*G2xx*G3x*G5xx + 6*G3px*G3x*G5xx + 20*G2xx*G4px*G5xx - 
      20*G3px*G4px*G5xx - 12*G2px*G4xx*G5xx + 12*G3pp*G4xx*G5xx - 4*G2xx*G5pp*G5xx + 4*G3px*G5pp*G5xx + 6*G2px*G5px*G5xx - 6*G3pp*G5px*G5xx + 3*G2x*G3x*G5xxx - 6*G3p*G3x*G5xxx + 6*G3x*G4pp*G5xxx - 
      6*G2x*G4px*G5xxx + 12*G3p*G4px*G5xxx - 12*G4pp*G4px*G5xxx + 4*G2px*G4x*G5xxx - 4*G3pp*G4x*G5xxx - 2*G2px*G5p*G5xxx + 2*G3pp*G5p*G5xxx + 384*pow(G4xx,3)*pow(H,2) + 528*G4x*G4xx*G4xxx*pow(H,2) - 
      408*G4xx*G4xxx*G5p*pow(H,2) - 584*pow(G4xx,2)*G5px*pow(H,2) - 352*G4x*G4xxx*G5px*pow(H,2) + 260*G4xxx*G5p*G5px*pow(H,2) + 308*G4xx*pow(G5px,2)*pow(H,2) - 56*pow(G5px,3)*pow(H,2) - 
      232*G4x*G4xx*G5pxx*pow(H,2) + 8*G4*G4xxx*G5pxx*pow(H,2) + 156*G4xx*G5p*G5pxx*pow(H,2) + 164*G4x*G5px*G5pxx*pow(H,2) - 108*G5p*G5px*G5pxx*pow(H,2) - 4*G4*pow(G5pxx,2)*pow(H,2) + 
      120*G3x*G4xxx*G5x*pow(H,2) - 336*G4px*G4xxx*G5x*pow(H,2) + 48*G4xxx*G5pp*G5x*pow(H,2) - 36*G4xx*G5ppx*G5x*pow(H,2) + 18*G5ppx*G5px*G5x*pow(H,2) - 56*G3x*G5pxx*G5x*pow(H,2) + 
      168*G4px*G5pxx*G5x*pow(H,2) - 24*G5pp*G5pxx*G5x*pow(H,2) + 162*G3x*G4xx*G5xx*pow(H,2) - 320*G4px*G4xx*G5xx*pow(H,2) - 24*G4p*G4xxx*G5xx*pow(H,2) + 4*G4xx*G5pp*G5xx*pow(H,2) - 
      40*G4x*G5ppx*G5xx*pow(H,2) + 26*G5p*G5ppx*G5xx*pow(H,2) - 82*G3x*G5px*G5xx*pow(H,2) + 156*G4px*G5px*G5xx*pow(H,2) - 2*G5pp*G5px*G5xx*pow(H,2) + 4*G4p*G5pxx*G5xx*pow(H,2) - 
      14*G2xx*G5x*G5xx*pow(H,2) + 24*G3px*G5x*G5xx*pow(H,2) - 20*G4ppx*G5x*G5xx*pow(H,2) - 4*G2x*pow(G5xx,2)*pow(H,2) + 8*G3p*pow(G5xx,2)*pow(H,2) + 2*G4pp*pow(G5xx,2)*pow(H,2) + 
      82*G3x*G4x*G5xxx*pow(H,2) - 228*G4px*G4x*G5xxx*pow(H,2) - 52*G4p*G4xx*G5xxx*pow(H,2) - 62*G3x*G5p*G5xxx*pow(H,2) + 172*G4px*G5p*G5xxx*pow(H,2) + 28*G4x*G5pp*G5xxx*pow(H,2) - 
      22*G5p*G5pp*G5xxx*pow(H,2) + 4*G4*G5ppx*G5xxx*pow(H,2) + 32*G4p*G5px*G5xxx*pow(H,2) + 11*G2x*G5x*G5xxx*pow(H,2) - 22*G3p*G5x*G5xxx*pow(H,2) + 14*G4pp*G5x*G5xxx*pow(H,2) + 
      36*G4xxx*pow(G5x,2)*pow(H,4) - 9*G5pxx*pow(G5x,2)*pow(H,4) + 398*G4xx*G5x*G5xx*pow(H,4) - 212*G5px*G5x*G5xx*pow(H,4) + 112*G4x*pow(G5xx,2)*pow(H,4) - 
      84*G5p*pow(G5xx,2)*pow(H,4) + 98*G4x*G5x*G5xxx*pow(H,4) - 86*G5p*G5x*G5xxx*pow(H,4) - 8*G4*G5xx*G5xxx*pow(H,4) + 
      G3xx*(18*G3x*G4xx - 36*G4px*G4xx + 12*G4pxx*(-2*G4x + G5p) + 12*G4x*G5ppx - 6*G5p*G5ppx - 12*G3x*G5px + 24*G4px*G5px - 3*G2x*G5xx + 6*G3p*G5xx - 6*G4pp*G5xx - 30*G4xx*G5x*pow(H,2) + 
         4*G5px*G5x*pow(H,2) - 18*G4x*G5xx*pow(H,2) + 18*G5p*G5xx*pow(H,2)) + 
      2*G4pxx*(18*G3x*G4xx - 36*G4px*G4xx - 12*G4x*G5ppx + 6*G5p*G5ppx - 6*G3x*G5px + 12*G4px*G5px - 4*G2xx*G5x + 4*G3px*G5x + 3*G2x*G5xx - 6*G3p*G5xx + 6*G4pp*G5xx + 66*G4xx*G5x*pow(H,2) - 
         22*G5px*G5x*pow(H,2) + 58*G4x*G5xx*pow(H,2) - 44*G5p*G5xx*pow(H,2) - 4*G4*G5xxx*pow(H,2)) + G5xx*G5xxx*pow(H,2)*rptot) + 
   pow(a,9)*pow(dphi,5)*(4*G2xx*G3px*G4 - 4*pow(G3px,2)*G4 - 4*G2px*G3xx*G4 + 4*G3pp*G3xx*G4 - 12*G3px*G3x*G4p - 12*G3xx*G4p*G4pp - 8*G2xx*G4*G4ppx + 8*G3px*G4*G4ppx + 24*G3x*G4p*G4ppx + 
      24*G3px*G4p*G4px - 12*G3x*G4pp*G4px - 48*G4p*G4ppx*G4px + 24*G4pp*pow(G4px,2) + 8*G2px*G4*G4pxx - 8*G3pp*G4*G4pxx + 24*G4p*G4pp*G4pxx + 4*G2px*G3x*G4x - 4*G3pp*G3x*G4x - 4*G2p*G3xx*G4x + 
      8*G2xx*G4pp*G4x - 8*G3px*G4pp*G4x - 16*G2px*G4px*G4x + 16*G3pp*G4px*G4x + 8*G2p*G4pxx*G4x + 4*G2p*G3x*G4xx + 8*G2px*G4p*G4xx - 8*G3pp*G4p*G4xx - 8*G2p*G4px*G4xx + 2*G2p*G3xx*G5p + 4*G2px*G4px*G5p - 
      4*G3pp*G4px*G5p - 4*G2p*G4pxx*G5p - 2*G2p*G3x*G5px - 4*G2px*G4p*G5px + 4*G3pp*G4p*G5px + 4*G2p*G4px*G5px + 8*pow(G3p,2)*(-2*G4xx + G5px) + pow(G2x,2)*(-4*G4xx + 2*G5px) - 
      48*G3x*G3xx*G4*pow(H,2) + 192*G3xx*G4*G4px*pow(H,2) + 144*G3x*G4*G4pxx*pow(H,2) - 528*G4*G4px*G4pxx*pow(H,2) + 114*pow(G3x,2)*G4x*pow(H,2) - 60*G3xx*G4p*G4x*pow(H,2) - 
      636*G3x*G4px*G4x*pow(H,2) + 1008*pow(G4px,2)*G4x*pow(H,2) - 24*G4p*G4pxx*G4x*pow(H,2) - 32*G2xx*pow(G4x,2)*pow(H,2) + 56*G3px*pow(G4x,2)*pow(H,2) - 48*G4ppx*pow(G4x,2)*pow(H,2) + 
      8*G2xx*G4*G4xx*pow(H,2) - 104*G3px*G4*G4xx*pow(H,2) - 360*G3x*G4p*G4xx*pow(H,2) + 192*G4*G4ppx*G4xx*pow(H,2) + 792*G4p*G4px*G4xx*pow(H,2) + 240*G4pp*G4x*G4xx*pow(H,2) + 
      48*pow(G4p,2)*G4xxx*pow(H,2) - 96*G4*G4pp*G4xxx*pow(H,2) - 66*pow(G3x,2)*G5p*pow(H,2) + 48*G3xx*G4p*G5p*pow(H,2) + 300*G3x*G4px*G5p*pow(H,2) - 432*pow(G4px,2)*G5p*pow(H,2) + 
      48*G4p*G4pxx*G5p*pow(H,2) + 60*G2xx*G4x*G5p*pow(H,2) - 96*G3px*G4x*G5p*pow(H,2) + 72*G4ppx*G4x*G5p*pow(H,2) - 48*G4pp*G4xx*G5p*pow(H,2) - 28*G2xx*pow(G5p,2)*pow(H,2) + 
      40*G3px*pow(G5p,2)*pow(H,2) - 24*G4ppx*pow(G5p,2)*pow(H,2) - 60*G3xx*G4*G5pp*pow(H,2) + 120*G4*G4pxx*G5pp*pow(H,2) + 12*G3x*G4x*G5pp*pow(H,2) - 144*G4px*G4x*G5pp*pow(H,2) + 
      72*G4p*G4xx*G5pp*pow(H,2) + 48*G3x*G5p*G5pp*pow(H,2) - 36*G4px*G5p*G5pp*pow(H,2) - 24*G3x*G4*G5ppx*pow(H,2) + 72*G4*G4px*G5ppx*pow(H,2) + 72*G4p*G4x*G5ppx*pow(H,2) - 
      72*G4p*G5p*G5ppx*pow(H,2) + 12*G2xx*G4*G5px*pow(H,2) + 48*G3px*G4*G5px*pow(H,2) + 192*G3x*G4p*G5px*pow(H,2) - 120*G4*G4ppx*G5px*pow(H,2) - 480*G4p*G4px*G5px*pow(H,2) - 
      192*G4pp*G4x*G5px*pow(H,2) + 72*G4pp*G5p*G5px*pow(H,2) - 12*G4p*G5pp*G5px*pow(H,2) - 12*pow(G4p,2)*G5pxx*pow(H,2) + 48*G4*G4pp*G5pxx*pow(H,2) + 16*G2xx*G4p*G5x*pow(H,2) - 
      16*G3px*G4p*G5x*pow(H,2) - 24*G3x*G4pp*G5x*pow(H,2) + 36*G4pp*G4px*G5x*pow(H,2) - 4*G2px*G4x*G5x*pow(H,2) + 4*G3pp*G4x*G5x*pow(H,2) + 4*G2p*G4xx*G5x*pow(H,2) + 
      16*G2px*G5p*G5x*pow(H,2) - 16*G3pp*G5p*G5x*pow(H,2) + 2*G2p*G5px*G5x*pow(H,2) - 28*G2px*G4*G5xx*pow(H,2) + 28*G3pp*G4*G5xx*pow(H,2) - 36*G4p*G4pp*G5xx*pow(H,2) - 
      16*G2p*G4x*G5xx*pow(H,2) + 2*G2p*G5p*G5xx*pow(H,2) + 4*G2p*G4*G5xxx*pow(H,2) + 1488*pow(G4x,2)*G4xx*pow(H,4) - 960*G4*pow(G4xx,2)*pow(H,4) - 240*G4*G4x*G4xxx*pow(H,4) - 
      2376*G4x*G4xx*G5p*pow(H,4) + 240*G4*G4xxx*G5p*pow(H,4) + 888*G4xx*pow(G5p,2)*pow(H,4) - 864*pow(G4x,2)*G5px*pow(H,4) + 1336*G4*G4xx*G5px*pow(H,4) + 1336*G4x*G5p*G5px*pow(H,4) - 
      472*pow(G5p,2)*G5px*pow(H,4) - 460*G4*pow(G5px,2)*pow(H,4) + 128*G4*G4x*G5pxx*pow(H,4) - 128*G4*G5p*G5pxx*pow(H,4) - 48*G4*G4pxx*G5x*pow(H,4) + 120*G3x*G4x*G5x*pow(H,4) - 
      84*G4px*G4x*G5x*pow(H,4) - 528*G4p*G4xx*G5x*pow(H,4) - 252*G4px*G5p*G5x*pow(H,4) - 108*G4x*G5pp*G5x*pow(H,4) + 144*G5p*G5pp*G5x*pow(H,4) + 24*G4*G5ppx*G5x*pow(H,4) + 
      252*G4p*G5px*G5x*pow(H,4) - 24*G4pp*pow(G5x,2)*pow(H,4) - 280*G3x*G4*G5xx*pow(H,4) + 784*G4*G4px*G5xx*pow(H,4) - 548*G4p*G4x*G5xx*pow(H,4) + 408*G4p*G5p*G5xx*pow(H,4) - 
      84*G4*G5pp*G5xx*pow(H,4) + 40*G4*G4p*G5xxx*pow(H,4) + 246*G4x*pow(G5x,2)*pow(H,6) - 174*G5p*pow(G5x,2)*pow(H,6) - 216*G4*G5x*G5xx*pow(H,6) - 
      2*G3p*(3*pow(G3x,2) - 6*G3xx*G4p - 18*G3x*G4px + 24*pow(G4px,2) + 12*G4p*G4pxx - 4*G3px*G4x + 8*G4ppx*G4x - 8*G4pp*G4xx + 2*G2xx*G5p - 4*G4ppx*G5p + 4*G4pp*G5px + 160*G4x*G4xx*pow(H,2) - 
         72*G4*G4xxx*pow(H,2) - 20*G4xx*G5p*pow(H,2) - 140*G4x*G5px*pow(H,2) + 44*G5p*G5px*pow(H,2) + 40*G4*G5pxx*pow(H,2) - 28*G3x*G5x*pow(H,2) + 58*G4px*G5x*pow(H,2) - 
         14*G4p*G5xx*pow(H,2) - 39*pow(G5x,2)*pow(H,4)) + G2x*(3*pow(G3x,2) - 6*G3xx*G4p + 24*pow(G4px,2) + 12*G4p*G4pxx - 4*G3px*G4x + 8*G4ppx*G4x + 16*G3p*G4xx - 8*G4pp*G4xx + 2*G2xx*G5p - 
         4*G4ppx*G5p - 8*G3p*G5px + 4*G4pp*G5px + 160*G4x*G4xx*pow(H,2) - 72*G4*G4xxx*pow(H,2) - 20*G4xx*G5p*pow(H,2) - 140*G4x*G5px*pow(H,2) + 44*G5p*G5px*pow(H,2) + 40*G4*G5pxx*pow(H,2) + 
         58*G4px*G5x*pow(H,2) - 14*G4p*G5xx*pow(H,2) - 39*pow(G5x,2)*pow(H,4) - 2*G3x*(9*G4px + 14*G5x*pow(H,2))) + 3*G3x*G3xx*rptot - 6*G3xx*G4px*rptot - 6*G3x*G4pxx*rptot + 
      12*G4px*G4pxx*rptot - 4*G2xx*G4xx*rptot + 4*G3px*G4xx*rptot + 2*G2xx*G5px*rptot - 2*G3px*G5px*rptot + 96*pow(G4xx,2)*pow(H,2)*rptot + 24*G4x*G4xxx*pow(H,2)*rptot - 
      24*G4xxx*G5p*pow(H,2)*rptot - 108*G4xx*G5px*pow(H,2)*rptot + 30*pow(G5px,2)*pow(H,2)*rptot - 12*G4x*G5pxx*pow(H,2)*rptot + 12*G5p*G5pxx*pow(H,2)*rptot - 9*G3xx*G5x*pow(H,2)*rptot + 
      18*G4pxx*G5x*pow(H,2)*rptot + 15*G3x*G5xx*pow(H,2)*rptot - 24*G4px*G5xx*pow(H,2)*rptot - 6*G4p*G5xxx*pow(H,2)*rptot + 27*G5x*G5xx*pow(H,4)*rptot) + 
   pow(a,7)*pow(dphi,7)*(-6*G3p*G3x*G3xx + 6*G3x*G3xx*G4pp - 6*pow(G3x,2)*G4ppx + 12*G3p*G3xx*G4px - 12*G3xx*G4pp*G4px + 24*G3x*G4ppx*G4px - 24*G4ppx*pow(G4px,2) + 12*G3p*G3x*G4pxx - 
      12*G3x*G4pp*G4pxx - 24*G3p*G4px*G4pxx + 24*G4pp*G4px*G4pxx + 4*G2px*G3xx*G4x - 4*G3pp*G3xx*G4x + 8*G2xx*G4ppx*G4x - 8*G2px*G4pxx*G4x + 8*G3pp*G4pxx*G4x + 8*G2xx*G3p*G4xx - 4*G2px*G3x*G4xx + 
      4*G3pp*G3x*G4xx - 8*G2xx*G4pp*G4xx + 8*G2px*G4px*G4xx - 8*G3pp*G4px*G4xx + pow(G3px,2)*(4*G4x - 2*G5p) - 2*G2px*G3xx*G5p + 2*G3pp*G3xx*G5p - 4*G2xx*G4ppx*G5p + 4*G2px*G4pxx*G5p - 
      4*G3pp*G4pxx*G5p - 4*G2xx*G3p*G5px + 2*G2px*G3x*G5px - 2*G3pp*G3x*G5px + 4*G2xx*G4pp*G5px - 4*G2px*G4px*G5px + 4*G3pp*G4px*G5px + 24*G3xx*G4*G4pxx*pow(H,2) - 48*G4*pow(G4pxx,2)*pow(H,2) + 
      66*G3x*G3xx*G4x*pow(H,2) - 228*G3xx*G4px*G4x*pow(H,2) - 108*G3x*G4pxx*G4x*pow(H,2) + 456*G4px*G4pxx*G4x*pow(H,2) + 126*pow(G3x,2)*G4xx*pow(H,2) - 36*G3xx*G4p*G4xx*pow(H,2) - 
      540*G3x*G4px*G4xx*pow(H,2) + 576*pow(G4px,2)*G4xx*pow(H,2) - 72*G4p*G4pxx*G4xx*pow(H,2) - 56*G2xx*G4x*G4xx*pow(H,2) - 192*G4ppx*G4x*G4xx*pow(H,2) + 48*G3p*pow(G4xx,2)*pow(H,2) - 
      96*G3x*G4p*G4xxx*pow(H,2) + 48*G4*G4ppx*G4xxx*pow(H,2) + 288*G4p*G4px*G4xxx*pow(H,2) - 192*G3p*G4x*G4xxx*pow(H,2) + 144*G4pp*G4x*G4xxx*pow(H,2) - 42*G3x*G3xx*G5p*pow(H,2) + 
      132*G3xx*G4px*G5p*pow(H,2) + 36*G3x*G4pxx*G5p*pow(H,2) - 192*G4px*G4pxx*G5p*pow(H,2) + 52*G2xx*G4xx*G5p*pow(H,2) + 96*G4ppx*G4xx*G5p*pow(H,2) + 120*G3p*G4xxx*G5p*pow(H,2) - 
      96*G4pp*G4xxx*G5p*pow(H,2) + 60*G3xx*G4x*G5pp*pow(H,2) - 120*G4pxx*G4x*G5pp*pow(H,2) - 36*G3x*G4xx*G5pp*pow(H,2) + 72*G4px*G4xx*G5pp*pow(H,2) - 48*G4p*G4xxx*G5pp*pow(H,2) - 
      30*G3xx*G5p*G5pp*pow(H,2) + 60*G4pxx*G5p*G5pp*pow(H,2) - 12*G3xx*G4*G5ppx*pow(H,2) + 24*G4*G4pxx*G5ppx*pow(H,2) - 12*G3x*G4x*G5ppx*pow(H,2) + 72*G4p*G4xx*G5ppx*pow(H,2) + 
      24*G3x*G5p*G5ppx*pow(H,2) - 36*G4px*G5p*G5ppx*pow(H,2) - 75*pow(G3x,2)*G5px*pow(H,2) + 24*G3xx*G4p*G5px*pow(H,2) + 348*G3x*G4px*G5px*pow(H,2) - 396*pow(G4px,2)*G5px*pow(H,2) + 
      24*G4p*G4pxx*G5px*pow(H,2) + 16*G2xx*G4x*G5px*pow(H,2) + 120*G4ppx*G4x*G5px*pow(H,2) - 32*G3p*G4xx*G5px*pow(H,2) - 24*G4pp*G4xx*G5px*pow(H,2) - 22*G2xx*G5p*G5px*pow(H,2) - 
      60*G4ppx*G5p*G5px*pow(H,2) + 6*G3x*G5pp*G5px*pow(H,2) - 12*G4px*G5pp*G5px*pow(H,2) - 36*G4p*G5ppx*G5px*pow(H,2) + 4*G3p*pow(G5px,2)*pow(H,2) + 12*G4pp*pow(G5px,2)*pow(H,2) + 
      4*G2xx*G4*G5pxx*pow(H,2) + 36*G3x*G4p*G5pxx*pow(H,2) - 24*G4*G4ppx*G5pxx*pow(H,2) - 120*G4p*G4px*G5pxx*pow(H,2) + 104*G3p*G4x*G5pxx*pow(H,2) - 72*G4pp*G4x*G5pxx*pow(H,2) - 
      64*G3p*G5p*G5pxx*pow(H,2) + 48*G4pp*G5p*G5pxx*pow(H,2) + 24*G4p*G5pp*G5pxx*pow(H,2) - 16*G2xx*G3x*G5x*pow(H,2) - 6*G3p*G3xx*G5x*pow(H,2) + 6*G3xx*G4pp*G5x*pow(H,2) + 
      40*G2xx*G4px*G5x*pow(H,2) + 28*G3p*G4pxx*G5x*pow(H,2) - 12*G4pp*G4pxx*G5x*pow(H,2) - 4*G2px*G4xx*G5x*pow(H,2) + 4*G3pp*G4xx*G5x*pow(H,2) - 8*G2p*G4xxx*G5x*pow(H,2) - 
      8*G3p*G5ppx*G5x*pow(H,2) - 2*G2px*G5px*G5x*pow(H,2) + 2*G3pp*G5px*G5x*pow(H,2) + 4*G2p*G5pxx*G5x*pow(H,2) + 6*G3p*G3x*G5xx*pow(H,2) + 4*G2xx*G4p*G5xx*pow(H,2) + 
      6*G3x*G4pp*G5xx*pow(H,2) - 40*G3p*G4px*G5xx*pow(H,2) + 16*G2px*G4x*G5xx*pow(H,2) - 16*G3pp*G4x*G5xx*pow(H,2) + 12*G2p*G4xx*G5xx*pow(H,2) - 2*G2px*G5p*G5xx*pow(H,2) + 
      2*G3pp*G5p*G5xx*pow(H,2) + 8*G3p*G5pp*G5xx*pow(H,2) - 6*G2p*G5px*G5xx*pow(H,2) - 4*G2px*G4*G5xxx*pow(H,2) + 4*G3pp*G4*G5xxx*pow(H,2) + 12*G3p*G4p*G5xxx*pow(H,2) - 
      12*G4p*G4pp*G5xxx*pow(H,2) - 4*G2p*G4x*G5xxx*pow(H,2) + 2*G2p*G5p*G5xxx*pow(H,2) + 1680*G4x*pow(G4xx,2)*pow(H,4) + 384*pow(G4x,2)*G4xxx*pow(H,4) - 240*G4*G4xx*G4xxx*pow(H,4) - 
      1200*pow(G4xx,2)*G5p*pow(H,4) - 648*G4x*G4xxx*G5p*pow(H,4) + 264*G4xxx*pow(G5p,2)*pow(H,4) - 2000*G4x*G4xx*G5px*pow(H,4) + 184*G4*G4xxx*G5px*pow(H,4) + 1332*G4xx*G5p*G5px*pow(H,4) + 
      624*G4x*pow(G5px,2)*pow(H,4) - 394*G5p*pow(G5px,2)*pow(H,4) - 168*pow(G4x,2)*G5pxx*pow(H,4) + 152*G4*G4xx*G5pxx*pow(H,4) + 272*G4x*G5p*G5pxx*pow(H,4) - 
      104*pow(G5p,2)*G5pxx*pow(H,4) - 112*G4*G5px*G5pxx*pow(H,4) - 30*G3xx*G4x*G5x*pow(H,4) + 180*G4pxx*G4x*G5x*pow(H,4) + 240*G3x*G4xx*G5x*pow(H,4) - 420*G4px*G4xx*G5x*pow(H,4) - 
      144*G4p*G4xxx*G5x*pow(H,4) + 30*G3xx*G5p*G5x*pow(H,4) - 156*G4pxx*G5p*G5x*pow(H,4) - 36*G4xx*G5pp*G5x*pow(H,4) - 60*G4x*G5ppx*G5x*pow(H,4) + 48*G5p*G5ppx*G5x*pow(H,4) - 
      160*G3x*G5px*G5x*pow(H,4) + 360*G4px*G5px*G5x*pow(H,4) - 6*G5pp*G5px*G5x*pow(H,4) + 48*G4p*G5pxx*G5x*pow(H,4) - 24*G2xx*pow(G5x,2)*pow(H,4) + 6*G4ppx*pow(G5x,2)*pow(H,4) - 
      56*G4*G4pxx*G5xx*pow(H,4) + 406*G3x*G4x*G5xx*pow(H,4) - 1032*G4px*G4x*G5xx*pow(H,4) - 340*G4p*G4xx*G5xx*pow(H,4) - 266*G3x*G5p*G5xx*pow(H,4) + 640*G4px*G5p*G5xx*pow(H,4) + 
      88*G4x*G5pp*G5xx*pow(H,4) - 46*G5p*G5pp*G5xx*pow(H,4) + 28*G4*G5ppx*G5xx*pow(H,4) + 184*G4p*G5px*G5xx*pow(H,4) - 30*G3p*G5x*G5xx*pow(H,4) + 38*G4pp*G5x*G5xx*pow(H,4) - 
      40*G3x*G4*G5xxx*pow(H,4) + 112*G4*G4px*G5xxx*pow(H,4) - 92*G4p*G4x*G5xxx*pow(H,4) + 72*G4p*G5p*G5xxx*pow(H,4) - 12*G4*G5pp*G5xxx*pow(H,4) + 210*G4xx*pow(G5x,2)*pow(H,6) - 
      93*G5px*pow(G5x,2)*pow(H,6) + 578*G4x*G5x*G5xx*pow(H,6) - 470*G5p*G5x*G5xx*pow(H,6) - 56*G4*pow(G5xx,2)*pow(H,6) - 24*G4*G5x*G5xxx*pow(H,6) + 
      G3px*(3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x - 8*G4ppx*G4x + 4*G2x*G4xx - 8*G3p*G4xx + 8*G4pp*G4xx + 2*G2xx*G5p + 4*G4ppx*G5p - 2*G2x*G5px + 4*G3p*G5px - 4*G4pp*G5px + 
         152*G4x*G4xx*pow(H,2) - 24*G4*G4xxx*pow(H,2) - 100*G4xx*G5p*pow(H,2) - 76*G4x*G5px*pow(H,2) + 52*G5p*G5px*pow(H,2) + 8*G4*G5pxx*pow(H,2) + 16*G3x*G5x*pow(H,2) - 
         40*G4px*G5x*pow(H,2) - 4*G4p*G5xx*pow(H,2) + 21*pow(G5x,2)*pow(H,4)) + 
      G2x*(-6*G3xx*G4px + 12*G4px*G4pxx - 4*G2xx*G4xx + 2*G2xx*G5px - 24*pow(G4xx,2)*pow(H,2) + 96*G4x*G4xxx*pow(H,2) - 60*G4xxx*G5p*pow(H,2) + 16*G4xx*G5px*pow(H,2) - 
         2*pow(G5px,2)*pow(H,2) - 52*G4x*G5pxx*pow(H,2) + 32*G5p*G5pxx*pow(H,2) + 3*G3xx*G5x*pow(H,2) - 14*G4pxx*G5x*pow(H,2) + 4*G5ppx*G5x*pow(H,2) + 20*G4px*G5xx*pow(H,2) - 
         4*G5pp*G5xx*pow(H,2) - 6*G4p*G5xxx*pow(H,2) + 15*G5x*G5xx*pow(H,4) + 3*G3x*(G3xx - 2*G4pxx - G5xx*pow(H,2))) + 24*G4xx*G4xxx*pow(H,2)*rptot - 12*G4xxx*G5px*pow(H,2)*rptot - 
      12*G4xx*G5pxx*pow(H,2)*rptot + 6*G5px*G5pxx*pow(H,2)*rptot - 3*G3xx*G5xx*pow(H,2)*rptot + 6*G4pxx*G5xx*pow(H,2)*rptot + 3*G3x*G5xxx*pow(H,2)*rptot - 6*G4px*G5xxx*pow(H,2)*rptot + 
      7*pow(G5xx,2)*pow(H,4)*rptot + 3*G5x*G5xxx*pow(H,4)*rptot) - 4*pow(a,13)*dphi*
    (pow(G2x,2)*G4 + 4*pow(G3p,2)*G4 - 2*G2p*G3x*G4 + 12*pow(G4p,2)*G4pp + 6*G2p*G4*G4px + 2*G2p*G4p*G4x - 2*G2p*G4p*G5p - 12*G3x*G4*G4p*pow(H,2) + 36*G4*G4p*G4px*pow(H,2) - 
      18*pow(G4p,2)*G4x*pow(H,2) + 24*G4*G4pp*G4x*pow(H,2) + 18*pow(G4p,2)*G5p*pow(H,2) - 24*G4*G4pp*G5p*pow(H,2) - 6*G2p*G4*G5x*pow(H,2) + 60*G4*pow(G4x,2)*pow(H,4) - 
      120*G4*G4x*G5p*pow(H,4) + 60*G4*pow(G5p,2)*pow(H,4) - 60*G4*G4p*G5x*pow(H,4) + 3*G3x*G4p*rptot - 9*G4p*G4px*rptot - 6*pow(G4x,2)*pow(H,2)*rptot + 12*G4x*G5p*pow(H,2)*rptot - 
      6*pow(G5p,2)*pow(H,2)*rptot + 9*G4p*G5x*pow(H,2)*rptot - 2*G3p*(3*pow(G4p,2) + 4*G4*(G4pp + 4*(G4x - G5p)*pow(H,2)) + (G4x - G5p)*rptot) + 
      G2x*(-4*G3p*G4 + 3*pow(G4p,2) + 4*G4*G4pp + 16*G4*G4x*pow(H,2) - 16*G4*G5p*pow(H,2) + G4x*rptot - G5p*rptot)) + 
   2*pow(a,11)*pow(dphi,3)*(-8*G3p*G3px*G4 - 4*G2px*G3x*G4 + 4*G3pp*G3x*G4 + 2*G2p*G3xx*G4 + 6*G3px*pow(G4p,2) + 8*G3px*G4*G4pp + 12*G3x*G4p*G4pp + 8*G3p*G4*G4ppx - 12*pow(G4p,2)*G4ppx + 
      12*G2px*G4*G4px - 12*G3pp*G4*G4px - 12*G3p*G4p*G4px - 12*G4p*G4pp*G4px - 4*G2p*G4*G4pxx - 2*G2p*G3x*G4x + 4*G2px*G4p*G4x - 4*G3pp*G4p*G4x - 8*G3p*G4pp*G4x + 8*G2p*G4px*G4x - 4*G2p*G4p*G4xx + 
      pow(G2x,2)*G5p + 4*pow(G3p,2)*G5p - 4*G2px*G4p*G5p + 4*G3pp*G4p*G5p - 2*G2p*G4px*G5p + 2*G2p*G4p*G5px - 48*pow(G3x,2)*G4*pow(H,2) + 12*G3xx*G4*G4p*pow(H,2) + 336*G3x*G4*G4px*pow(H,2) - 
      576*G4*pow(G4px,2)*pow(H,2) - 24*G4*G4p*G4pxx*pow(H,2) - 16*G3px*G4*G4x*pow(H,2) - 60*G3x*G4p*G4x*pow(H,2) + 24*G4*G4ppx*G4x*pow(H,2) + 36*G4p*G4px*G4x*pow(H,2) + 
      32*G3p*pow(G4x,2)*pow(H,2) - 24*G4pp*pow(G4x,2)*pow(H,2) + 280*G3p*G4*G4xx*pow(H,2) + 108*pow(G4p,2)*G4xx*pow(H,2) - 192*G4*G4pp*G4xx*pow(H,2) + 16*G3px*G4*G5p*pow(H,2) + 
      48*G3x*G4p*G5p*pow(H,2) - 24*G4*G4ppx*G5p*pow(H,2) - 96*G3p*G4x*G5p*pow(H,2) + 72*G4pp*G4x*G5p*pow(H,2) + 64*G3p*pow(G5p,2)*pow(H,2) - 48*G4pp*pow(G5p,2)*pow(H,2) - 
      60*G3x*G4*G5pp*pow(H,2) + 180*G4*G4px*G5pp*pow(H,2) + 108*G4p*G4x*G5pp*pow(H,2) - 108*G4p*G5p*G5pp*pow(H,2) - 192*G3p*G4*G5px*pow(H,2) - 42*pow(G4p,2)*G5px*pow(H,2) + 
      120*G4*G4pp*G5px*pow(H,2) - 12*G2px*G4*G5x*pow(H,2) + 12*G3pp*G4*G5x*pow(H,2) - 52*G3p*G4p*G5x*pow(H,2) + 36*G4p*G4pp*G5x*pow(H,2) + 2*G2p*G4x*G5x*pow(H,2) - 8*G2p*G5p*G5x*pow(H,2) + 
      14*G2p*G4*G5xx*pow(H,2) + 96*pow(G4x,3)*pow(H,4) - 600*G4*G4x*G4xx*pow(H,4) - 228*pow(G4x,2)*G5p*pow(H,4) + 600*G4*G4xx*G5p*pow(H,4) + 168*G4x*pow(G5p,2)*pow(H,4) - 
      36*pow(G5p,3)*pow(H,4) + 392*G4*G4x*G5px*pow(H,4) - 392*G4*G5p*G5px*pow(H,4) - 120*G3x*G4*G5x*pow(H,4) + 336*G4*G4px*G5x*pow(H,4) - 168*G4p*G4x*G5x*pow(H,4) + 
      108*G4p*G5p*G5x*pow(H,4) - 36*G4*G5pp*G5x*pow(H,4) + 140*G4*G4p*G5xx*pow(H,4) - 72*G4*pow(G5x,2)*pow(H,6) + 3*pow(G3x,2)*rptot - 3*G3xx*G4p*rptot - 15*G3x*G4px*rptot + 
      18*pow(G4px,2)*rptot + 6*G4p*G4pxx*rptot + 2*G3px*G4x*rptot + 4*G3p*G4xx*rptot - 2*G3px*G5p*rptot - 2*G3p*G5px*rptot + 60*G4x*G4xx*pow(H,2)*rptot - 60*G4xx*G5p*pow(H,2)*rptot - 
      36*G4x*G5px*pow(H,2)*rptot + 36*G5p*G5px*pow(H,2)*rptot + 9*G4px*G5x*pow(H,2)*rptot - 21*G4p*G5xx*pow(H,2)*rptot + 9*pow(G5x,2)*pow(H,4)*rptot + 
      2*G2xx*(2*G3p*G4 - 4*G4*G4pp + 2*G4*G4x*pow(H,2) - 2*G4*G5p*pow(H,2) - G4x*rptot + G5p*rptot) + 
      G2x*(-2*G2xx*G4 + 4*G3px*G4 - 4*G4*G4ppx + 6*G4p*G4px + 4*G4pp*G4x - 4*G3p*G5p - 16*pow(G4x,2)*pow(H,2) - 140*G4*G4xx*pow(H,2) + 48*G4x*G5p*pow(H,2) - 32*pow(G5p,2)*pow(H,2) + 
         96*G4*G5px*pow(H,2) + 26*G4p*G5x*pow(H,2) - 2*G4xx*rptot + G5px*rptot)) + 
   pow(a,8)*pow(dphi,6)*H*(15*pow(G3x,3) - 240*pow(G4px,3) + 16*G2xx*G4*G4pxx - 16*G3px*G4*G4pxx - 48*G4p*G4px*G4pxx - 24*G3px*G4px*G4x + 48*G4ppx*G4px*G4x - 64*G2x*G4pxx*G4x + 128*G3p*G4pxx*G4x - 
      96*G4pp*G4pxx*G4x + 8*G2xx*G4p*G4xx - 32*G3px*G4p*G4xx + 48*G4p*G4ppx*G4xx + 16*G2x*G4px*G4xx - 32*G3p*G4px*G4xx + 32*G2px*G4x*G4xx - 32*G3pp*G4x*G4xx + 16*G2p*pow(G4xx,2) - 16*G2px*G4*G4xxx + 
      16*G3pp*G4*G4xxx - 24*G2x*G4p*G4xxx + 48*G3p*G4p*G4xxx - 48*G4p*G4pp*G4xxx - 16*G2p*G4x*G4xxx - 20*G2xx*G4px*G5p + 44*G3px*G4px*G5p - 48*G4ppx*G4px*G5p + 32*G2x*G4pxx*G5p - 64*G3p*G4pxx*G5p + 
      48*G4pp*G4pxx*G5p + 8*G2p*G4xxx*G5p + 24*G4p*G4pxx*G5pp + 16*G2xx*G4x*G5pp - 16*G3px*G4x*G5pp - 8*G2x*G4xx*G5pp + 16*G3p*G4xx*G5pp - 4*G2xx*G5p*G5pp + 4*G3px*G5p*G5pp - 8*G2xx*G4*G5ppx + 
      8*G3px*G4*G5ppx - 48*G4p*G4px*G5ppx + 8*G2x*G4x*G5ppx - 16*G3p*G4x*G5ppx - 4*G2x*G5p*G5ppx + 8*G3p*G5p*G5ppx - 4*G2xx*G4p*G5px + 16*G3px*G4p*G5px - 24*G4p*G4ppx*G5px + 4*G2x*G4px*G5px - 
      8*G3p*G4px*G5px + 24*G4pp*G4px*G5px - 24*G2px*G4x*G5px + 24*G3pp*G4x*G5px - 16*G2p*G4xx*G5px + 4*G2px*G5p*G5px - 4*G3pp*G5p*G5px + 4*G2x*G5pp*G5px - 8*G3p*G5pp*G5px + 4*G2p*pow(G5px,2) + 
      8*G2px*G4*G5pxx - 8*G3pp*G4*G5pxx + 12*G2x*G4p*G5pxx - 24*G3p*G4p*G5pxx + 24*G4p*G4pp*G5pxx + 8*G2p*G4x*G5pxx - 4*G2p*G5p*G5pxx - 4*G2x*G2xx*G5x + 8*G2xx*G3p*G5x + 2*G2x*G3px*G5x - 4*G3p*G3px*G5x - 
      4*G2xx*G4pp*G5x + 4*G3px*G4pp*G5x + 4*G2x*G4ppx*G5x - 8*G3p*G4ppx*G5x + 4*G2p*G4pxx*G5x - 2*pow(G2x,2)*G5xx + 8*G2x*G3p*G5xx - 8*pow(G3p,2)*G5xx + 4*G2px*G4p*G5xx - 4*G3pp*G4p*G5xx - 
      4*G2x*G4pp*G5xx + 8*G3p*G4pp*G5xx - 4*G2p*G4px*G5xx + 96*G4*G4pxx*G4xx*pow(H,2) - 2688*G4px*G4x*G4xx*pow(H,2) - 624*G4p*pow(G4xx,2)*pow(H,2) + 528*G4*G4px*G4xxx*pow(H,2) - 
      336*G4p*G4x*G4xxx*pow(H,2) - 48*G4pxx*G4x*G5p*pow(H,2) + 1488*G4px*G4xx*G5p*pow(H,2) + 264*G4p*G4xxx*G5p*pow(H,2) + 48*G4pxx*pow(G5p,2)*pow(H,2) + 144*G4x*G4xx*G5pp*pow(H,2) - 
      96*G4*G4xxx*G5pp*pow(H,2) + 48*G4xx*G5p*G5pp*pow(H,2) - 48*pow(G4x,2)*G5ppx*pow(H,2) + 96*G4x*G5p*G5ppx*pow(H,2) - 48*pow(G5p,2)*G5ppx*pow(H,2) - 128*G4*G4pxx*G5px*pow(H,2) + 
      1888*G4px*G4x*G5px*pow(H,2) + 640*G4p*G4xx*G5px*pow(H,2) - 1036*G4px*G5p*G5px*pow(H,2) - 144*G4x*G5pp*G5px*pow(H,2) + 24*G5p*G5pp*G5px*pow(H,2) - 176*G4p*pow(G5px,2)*pow(H,2) - 
      312*G4*G4px*G5pxx*pow(H,2) + 128*G4p*G4x*G5pxx*pow(H,2) - 92*G4p*G5p*G5pxx*pow(H,2) + 48*G4*G5pp*G5pxx*pow(H,2) - 36*pow(G4px,2)*G5x*pow(H,2) - 48*G4p*G4pxx*G5x*pow(H,2) - 
      64*G2xx*G4x*G5x*pow(H,2) + 88*G3px*G4x*G5x*pow(H,2) - 48*G4ppx*G4x*G5x*pow(H,2) - 40*G2x*G4xx*G5x*pow(H,2) + 80*G3p*G4xx*G5x*pow(H,2) + 58*G2xx*G5p*G5x*pow(H,2) - 
      58*G3px*G5p*G5x*pow(H,2) + 36*G4px*G5pp*G5x*pow(H,2) + 36*G4p*G5ppx*G5x*pow(H,2) - 24*G4pp*G5px*G5x*pow(H,2) - 6*G2px*pow(G5x,2)*pow(H,2) + 6*G3pp*pow(G5x,2)*pow(H,2) + 
      4*G2xx*G4*G5xx*pow(H,2) - 60*G3px*G4*G5xx*pow(H,2) + 112*G4*G4ppx*G5xx*pow(H,2) + 312*G4p*G4px*G5xx*pow(H,2) + 88*G2x*G4x*G5xx*pow(H,2) - 176*G3p*G4x*G5xx*pow(H,2) + 
      128*G4pp*G4x*G5xx*pow(H,2) - 34*G2x*G5p*G5xx*pow(H,2) + 68*G3p*G5p*G5xx*pow(H,2) - 72*G4pp*G5p*G5xx*pow(H,2) - 24*G4p*G5pp*G5xx*pow(H,2) - 2*G2p*G5x*G5xx*pow(H,2) - 
      16*G2x*G4*G5xxx*pow(H,2) + 32*G3p*G4*G5xxx*pow(H,2) + 12*pow(G4p,2)*G5xxx*pow(H,2) - 16*G4*G4pp*G5xxx*pow(H,2) + 1200*G4x*G4xx*G5x*pow(H,4) - 72*G4*G4xxx*G5x*pow(H,4) - 
      864*G4xx*G5p*G5x*pow(H,4) - 644*G4x*G5px*G5x*pow(H,4) + 434*G5p*G5px*G5x*pow(H,4) + 36*G4*G5pxx*G5x*pow(H,4) + 132*G4px*pow(G5x,2)*pow(H,4) - 36*G5pp*pow(G5x,2)*pow(H,4) + 
      664*pow(G4x,2)*G5xx*pow(H,4) - 544*G4*G4xx*G5xx*pow(H,4) - 1092*G4x*G5p*G5xx*pow(H,4) + 428*pow(G5p,2)*G5xx*pow(H,4) + 340*G4*G5px*G5xx*pow(H,4) - 248*G4p*G5x*G5xx*pow(H,4) - 
      64*G4*G4x*G5xxx*pow(H,4) + 64*G4*G5p*G5xxx*pow(H,4) + 45*pow(G5x,3)*pow(H,6) - 3*pow(G3x,2)*(36*G4px + 2*G5pp - G5x*pow(H,2)) + 
      2*G3xx*(12*G2x*G4x - 24*G3p*G4x + 24*G4pp*G4x - 6*G2x*G5p + 12*G3p*G5p - 12*G4pp*G5p - G2p*G5x + 24*pow(G4x,2)*pow(H,2) - 24*G4*G4xx*pow(H,2) - 36*G4x*G5p*pow(H,2) + 
         12*pow(G5p,2)*pow(H,2) + 32*G4*G5px*pow(H,2) + 6*G4p*(6*G4px - G5pp - G5x*pow(H,2))) - 24*G4px*G4xxx*rptot + 12*G4px*G5pxx*rptot - 2*G2xx*G5xx*rptot + 2*G3px*G5xx*rptot + 
      56*G4xx*G5xx*pow(H,2)*rptot - 28*G5px*G5xx*pow(H,2)*rptot + 8*G4x*G5xxx*pow(H,2)*rptot - 8*G5p*G5xxx*pow(H,2)*rptot + 
      G3x*(-24*G3xx*G4p + 276*pow(G4px,2) - 8*G2xx*G4x + 20*G3px*G4x - 24*G4ppx*G4x + 10*G2xx*G5p - 22*G3px*G5p + 24*G4ppx*G5p + 24*G4p*G5ppx - 6*G2x*G5px + 12*G3p*G5px - 12*G4pp*G5px - 2*G2px*G5x + 
         2*G3pp*G5x + 2*G2p*G5xx + 1056*G4x*G4xx*pow(H,2) - 168*G4*G4xxx*pow(H,2) - 672*G4xx*G5p*pow(H,2) - 688*G4x*G5px*pow(H,2) + 414*G5p*G5px*pow(H,2) + 100*G4*G5pxx*pow(H,2) - 
         30*G5pp*G5x*pow(H,2) - 120*G4p*G5xx*pow(H,2) - 15*pow(G5x,2)*pow(H,4) + 12*G4px*(G5pp + 2*G5x*pow(H,2)) + 12*G4xxx*rptot - 6*G5pxx*rptot)) - 
   2*pow(a,10)*pow(dphi,4)*H*(12*G2x*G3xx*G4 - 24*G3p*G3xx*G4 + 21*pow(G3x,2)*G4p - 6*G3xx*pow(G4p,2) + 24*G3xx*G4*G4pp - 84*G3x*G4p*G4px + 108*G4p*pow(G4px,2) - 32*G2x*G4*G4pxx + 
      64*G3p*G4*G4pxx - 12*pow(G4p,2)*G4pxx - 48*G4*G4pp*G4pxx - 8*G2x*G3x*G4x + 16*G3p*G3x*G4x - 12*G3x*G4pp*G4x - 24*G4p*G4ppx*G4x + 48*G2x*G4px*G4x - 96*G3p*G4px*G4x + 72*G4pp*G4px*G4x + 
      32*G2px*G4*G4xx - 32*G3pp*G4*G4xx + 8*G2x*G4p*G4xx - 16*G3p*G4p*G4xx + 16*G2p*G4x*G4xx - 8*G2p*G4*G4xxx - 5*G2x*G3x*G5p + 10*G3p*G3x*G5p - 12*G3x*G4pp*G5p + 24*G4p*G4ppx*G5p - 2*G2x*G4px*G5p + 
      4*G3p*G4px*G5p - 4*G2px*G4x*G5p + 4*G3pp*G4x*G5p + 4*G2px*pow(G5p,2) - 4*G3pp*pow(G5p,2) - 24*G3x*G4p*G5pp + 36*G4p*G4px*G5pp - 8*G2x*G4x*G5pp + 16*G3p*G4x*G5pp + 2*G2x*G5p*G5pp - 
      4*G3p*G5p*G5pp + 4*G2x*G4*G5ppx - 8*G3p*G4*G5ppx + 12*pow(G4p,2)*G5ppx - 20*G2px*G4*G5px + 20*G3pp*G4*G5px - 10*G2x*G4p*G5px + 20*G3p*G4p*G5px - 12*G4p*G4pp*G5px - 12*G2p*G4x*G5px + 
      2*G2p*G5p*G5px + 4*G2p*G4*G5pxx + 2*pow(G2x,2)*G5x - 8*G2x*G3p*G5x + 8*pow(G3p,2)*G5x - G2p*G3x*G5x - 6*G2px*G4p*G5x + 6*G3pp*G4p*G5x + 2*G2x*G4pp*G5x - 4*G3p*G4pp*G5x + 2*G2p*G4p*G5xx + 
      24*G3xx*G4*G4x*pow(H,2) - 48*G4*G4pxx*G4x*pow(H,2) - 96*G3x*pow(G4x,2)*pow(H,2) + 192*G4px*pow(G4x,2)*pow(H,2) + 384*G3x*G4*G4xx*pow(H,2) - 1200*G4*G4px*G4xx*pow(H,2) + 
      624*G4p*G4x*G4xx*pow(H,2) - 72*G4*G4p*G4xxx*pow(H,2) - 24*G3xx*G4*G5p*pow(H,2) + 48*G4*G4pxx*G5p*pow(H,2) + 126*G3x*G4x*G5p*pow(H,2) - 180*G4px*G4x*G5p*pow(H,2) - 
      480*G4p*G4xx*G5p*pow(H,2) - 30*G3x*pow(G5p,2)*pow(H,2) - 12*G4px*pow(G5p,2)*pow(H,2) + 48*pow(G4x,2)*G5pp*pow(H,2) + 192*G4*G4xx*G5pp*pow(H,2) - 120*G4x*G5p*G5pp*pow(H,2) + 
      72*pow(G5p,2)*G5pp*pow(H,2) - 274*G3x*G4*G5px*pow(H,2) + 852*G4*G4px*G5px*pow(H,2) - 320*G4p*G4x*G5px*pow(H,2) + 230*G4p*G5p*G5px*pow(H,2) - 120*G4*G5pp*G5px*pow(H,2) + 
      36*G4*G4p*G5pxx*pow(H,2) + 12*G3x*G4p*G5x*pow(H,2) - 48*G4*G4ppx*G5x*pow(H,2) + 36*G4p*G4px*G5x*pow(H,2) + 32*G2x*G4x*G5x*pow(H,2) - 64*G3p*G4x*G5x*pow(H,2) + 
      24*G4pp*G4x*G5x*pow(H,2) - 53*G2x*G5p*G5x*pow(H,2) + 106*G3p*G5p*G5x*pow(H,2) - 48*G4pp*G5p*G5x*pow(H,2) - 54*G4p*G5pp*G5x*pow(H,2) - 3*G2p*pow(G5x,2)*pow(H,2) + 
      54*G2x*G4*G5xx*pow(H,2) - 108*G3p*G4*G5xx*pow(H,2) - 36*pow(G4p,2)*G5xx*pow(H,2) + 56*G4*G4pp*G5xx*pow(H,2) - 216*pow(G4x,2)*G5x*pow(H,4) + 336*G4*G4xx*G5x*pow(H,4) + 
      318*G4x*G5p*G5x*pow(H,4) - 102*pow(G5p,2)*G5x*pow(H,4) - 210*G4*G5px*G5x*pow(H,4) + 51*G4p*pow(G5x,2)*pow(H,4) + 236*G4*G4x*G5xx*pow(H,4) - 236*G4*G5p*G5xx*pow(H,4) - 
      24*G3x*G4xx*rptot + 48*G4px*G4xx*rptot + 12*G4p*G4xxx*rptot + 15*G3x*G5px*rptot - 30*G4px*G5px*rptot - 6*G4p*G5pxx*rptot + G2x*G5xx*rptot - 2*G3p*G5xx*rptot - 24*G4xx*G5x*pow(H,2)*rptot + 
      12*G5px*G5x*pow(H,2)*rptot - 28*G4x*G5xx*pow(H,2)*rptot + 28*G5p*G5xx*pow(H,2)*rptot + 
      G3px*(-2*G3x*G4 + 20*G4*G4px + 20*G4p*G4x - 18*G4p*G5p - 12*G4*G5pp + 30*G4*G5x*pow(H,2) - 3*G5x*rptot) + 
      G2xx*(2*G3x*G4 - 20*G4*G4px - 8*G4p*G4x + 6*G4p*G5p + 12*G4*G5pp - 6*G4*G5x*pow(H,2) + 3*G5x*rptot)) + 
   pow(a,6)*pow(dphi,8)*H*(G3xx*(60*pow(G4px,2) + 2*G2px*G5x - 2*G3pp*G5x + 48*G4x*G4xx*pow(H,2) - 24*G4xx*G5p*pow(H,2) - 64*G4x*G5px*pow(H,2) + 32*G5p*G5px*pow(H,2) + 
         8*G4*G5pxx*pow(H,2) + 12*G5pp*G5x*pow(H,2) - 27*pow(G5x,2)*pow(H,4) - 12*G4px*(G5pp + 2*G5x*pow(H,2))) + pow(G3x,2)*(9*G3xx - 6*(G4pxx + G5ppx - 5*G5xx*pow(H,2))) + 
      2*G3x*(-6*G2xx*G4xx + 12*G3px*G4xx - 12*G4ppx*G4xx + 6*G2x*G4xxx - 12*G3p*G4xxx + 12*G4pp*G4xxx - 6*G4pxx*G5pp + 3*G2xx*G5px - 6*G3px*G5px + 6*G4ppx*G5px - 3*G2x*G5pxx + 6*G3p*G5pxx - 6*G4pp*G5pxx - 
         G2px*G5xx + G3pp*G5xx + 180*pow(G4xx,2)*pow(H,2) + 156*G4x*G4xxx*pow(H,2) - 114*G4xxx*G5p*pow(H,2) - 196*G4xx*G5px*pow(H,2) + 56*pow(G5px,2)*pow(H,2) - 76*G4x*G5pxx*pow(H,2) + 
         51*G5p*G5pxx*pow(H,2) - 3*G5ppx*G5x*pow(H,2) - 12*G4p*G5xxx*pow(H,2) + 70*G5x*G5xx*pow(H,4) + 3*G3xx*(-8*G4px + G5pp + G5x*pow(H,2)) + 12*G4px*(2*G4pxx + G5ppx - 5*G5xx*pow(H,2))) - 
      2*(-8*G3px*G4pxx*G4x + 8*G2px*pow(G4xx,2) - 8*G3pp*pow(G4xx,2) - 8*G2px*G4x*G4xxx + 8*G3pp*G4x*G4xxx + 4*G3px*G4pxx*G5p + 4*G2px*G4xxx*G5p - 4*G3pp*G4xxx*G5p - 4*G3px*G4xx*G5pp + 
         4*G3px*G4x*G5ppx - 2*G3px*G5p*G5ppx - 8*G2px*G4xx*G5px + 8*G3pp*G4xx*G5px + 2*G3px*G5pp*G5px + 2*G2px*pow(G5px,2) - 2*G3pp*pow(G5px,2) + 4*G2px*G4x*G5pxx - 4*G3pp*G4x*G5pxx - 
         2*G2px*G5p*G5pxx + 2*G3pp*G5p*G5pxx - pow(G3px,2)*G5x + 2*G3px*G4ppx*G5x + 2*G2px*G4pxx*G5x - 2*G3pp*G4pxx*G5x - G2x*G3px*G5xx + 2*G3p*G3px*G5xx - 2*G3px*G4pp*G5xx - 
         48*G4pxx*G4x*G4xx*pow(H,2) + 96*G4p*G4xx*G4xxx*pow(H,2) + 72*G4pxx*G4xx*G5p*pow(H,2) + 24*pow(G4xx,2)*G5pp*pow(H,2) - 72*G4x*G4xxx*G5pp*pow(H,2) + 48*G4xxx*G5p*G5pp*pow(H,2) + 
         48*G4x*G4xx*G5ppx*pow(H,2) - 48*G4xx*G5p*G5ppx*pow(H,2) - 16*G4pxx*G4x*G5px*pow(H,2) - 60*G4p*G4xxx*G5px*pow(H,2) - 16*G4pxx*G5p*G5px*pow(H,2) - 12*G4xx*G5pp*G5px*pow(H,2) - 
         24*G4x*G5ppx*G5px*pow(H,2) + 24*G5p*G5ppx*G5px*pow(H,2) + 8*G4*G4pxx*G5pxx*pow(H,2) - 28*G4p*G4xx*G5pxx*pow(H,2) + 36*G4x*G5pp*G5pxx*pow(H,2) - 24*G5p*G5pp*G5pxx*pow(H,2) + 
         20*G4p*G5px*G5pxx*pow(H,2) - 32*G3px*G4xx*G5x*pow(H,2) + 12*G4ppx*G4xx*G5x*pow(H,2) - 18*G2x*G4xxx*G5x*pow(H,2) + 36*G3p*G4xxx*G5x*pow(H,2) - 24*G4pp*G4xxx*G5x*pow(H,2) + 
         12*G4pxx*G5pp*G5x*pow(H,2) + 18*G3px*G5px*G5x*pow(H,2) - 12*G4ppx*G5px*G5x*pow(H,2) + 10*G2x*G5pxx*G5x*pow(H,2) - 20*G3p*G5pxx*G5x*pow(H,2) + 12*G4pp*G5pxx*G5x*pow(H,2) + 
         12*G4p*G4pxx*G5xx*pow(H,2) - 32*G3px*G4x*G5xx*pow(H,2) + 44*G4ppx*G4x*G5xx*pow(H,2) + 10*G2x*G4xx*G5xx*pow(H,2) - 20*G3p*G4xx*G5xx*pow(H,2) - 8*G4pp*G4xx*G5xx*pow(H,2) + 
         17*G3px*G5p*G5xx*pow(H,2) - 16*G4ppx*G5p*G5xx*pow(H,2) - 6*G4p*G5ppx*G5xx*pow(H,2) - 6*G2x*G5px*G5xx*pow(H,2) + 12*G3p*G5px*G5xx*pow(H,2) + 4*G4pp*G5px*G5xx*pow(H,2) - 
         G2px*G5x*G5xx*pow(H,2) + G3pp*G5x*G5xx*pow(H,2) - G2p*pow(G5xx,2)*pow(H,2) + 4*G3px*G4*G5xxx*pow(H,2) - 8*G4*G4ppx*G5xxx*pow(H,2) - 12*G2x*G4x*G5xxx*pow(H,2) + 
         24*G3p*G4x*G5xxx*pow(H,2) - 16*G4pp*G4x*G5xxx*pow(H,2) + 8*G2x*G5p*G5xxx*pow(H,2) - 16*G3p*G5p*G5xxx*pow(H,2) + 12*G4pp*G5p*G5xxx*pow(H,2) + 6*G4p*G5pp*G5xxx*pow(H,2) + 
         G2p*G5x*G5xxx*pow(H,2) - 276*pow(G4xx,2)*G5x*pow(H,4) - 132*G4x*G4xxx*G5x*pow(H,4) + 114*G4xxx*G5p*G5x*pow(H,4) + 304*G4xx*G5px*G5x*pow(H,4) - 91*pow(G5px,2)*G5x*pow(H,4) + 
         50*G4x*G5pxx*G5x*pow(H,4) - 41*G5p*G5pxx*G5x*pow(H,4) - 39*G4pxx*pow(G5x,2)*pow(H,4) + 6*G5ppx*pow(G5x,2)*pow(H,4) - 500*G4x*G4xx*G5xx*pow(H,4) + 12*G4*G4xxx*G5xx*pow(H,4) + 
         364*G4xx*G5p*G5xx*pow(H,4) + 282*G4x*G5px*G5xx*pow(H,4) - 197*G5p*G5px*G5xx*pow(H,4) - 6*G4*G5pxx*G5xx*pow(H,4) - 12*G5pp*G5x*G5xx*pow(H,4) + 19*G4p*pow(G5xx,2)*pow(H,4) - 
         56*pow(G4x,2)*G5xxx*pow(H,4) + 32*G4*G4xx*G5xxx*pow(H,4) + 96*G4x*G5p*G5xxx*pow(H,4) - 40*pow(G5p,2)*G5xxx*pow(H,4) - 20*G4*G5px*G5xxx*pow(H,4) + 22*G4p*G5x*G5xxx*pow(H,4) - 
         63*pow(G5x,2)*G5xx*pow(H,6) + 12*pow(G4px,2)*(3*G4pxx + G5ppx - 4*G5xx*pow(H,2)) + 
         G2xx*(8*G4pxx*G4x - 4*G4pxx*G5p + 4*G4xx*G5pp - 4*G4x*G5ppx + 2*G5p*G5ppx - 2*G5pp*G5px + G3px*G5x - 2*G4ppx*G5x + G2x*G5xx - 2*G3p*G5xx + 2*G4pp*G5xx + 26*G4xx*G5x*pow(H,2) - 
            12*G5px*G5x*pow(H,2) + 10*G4x*G5xx*pow(H,2) - 9*G5p*G5xx*pow(H,2)) - 
         2*G4px*(10*G2xx*G4xx - 16*G3px*G4xx + 12*G4ppx*G4xx - 6*G2x*G4xxx + 12*G3p*G4xxx - 12*G4pp*G4xxx + 6*G4pxx*G5pp - 5*G2xx*G5px + 8*G3px*G5px - 6*G4ppx*G5px + 3*G2x*G5pxx - 6*G3p*G5pxx + 
            6*G4pp*G5pxx + G2px*G5xx - G3pp*G5xx - 180*pow(G4xx,2)*pow(H,2) - 228*G4x*G4xxx*pow(H,2) + 162*G4xxx*G5p*pow(H,2) + 208*G4xx*G5px*pow(H,2) - 62*pow(G5px,2)*pow(H,2) + 
            116*G4x*G5pxx*pow(H,2) - 77*G5p*G5pxx*pow(H,2) + 12*G4pxx*G5x*pow(H,2) + 3*G5pp*G5xx*pow(H,2) + 18*G4p*G5xxx*pow(H,2) - 80*G5x*G5xx*pow(H,4)) - 4*G4xx*G5xxx*pow(H,2)*rptot + 
         2*G5px*G5xxx*pow(H,2)*rptot)) - 2*pow(a,12)*pow(dphi,2)*H*(G2x*(26*G3x*G4 - 92*G4*G4px - 32*G4p*G4x + 30*G4p*G5p + 12*G4*G5pp + 42*G4*G5x*pow(H,2) + 3*G5x*rptot) - 
      2*(-2*G2xx*G4*G4p + 2*G3px*G4*G4p + 9*G3x*pow(G4p,2) - 24*G3x*G4*G4pp + 72*G4*G4pp*G4px - 4*G2px*G4*G4x + 4*G3pp*G4*G4x + 36*G4p*G4pp*G4x + 16*G2p*G4*G4xx + 4*G2px*G4*G5p - 4*G3pp*G4*G5p - 
         36*G4p*G4pp*G5p - 2*G2p*G4x*G5p + 2*G2p*pow(G5p,2) - 18*pow(G4p,2)*G5pp - 10*G2p*G4*G5px - 3*G2p*G4p*G5x - 66*G3x*G4*G4x*pow(H,2) + 204*G4*G4px*G4x*pow(H,2) - 
         48*G4p*pow(G4x,2)*pow(H,2) + 144*G4*G4p*G4xx*pow(H,2) + 66*G3x*G4*G5p*pow(H,2) - 204*G4*G4px*G5p*pow(H,2) + 78*G4p*G4x*G5p*pow(H,2) - 30*G4p*pow(G5p,2)*pow(H,2) - 
         24*G4*G4x*G5pp*pow(H,2) + 24*G4*G5p*G5pp*pow(H,2) - 90*G4*G4p*G5px*pow(H,2) + 9*pow(G4p,2)*G5x*pow(H,2) - 24*G4*G4pp*G5x*pow(H,2) - 114*G4*G4x*G5x*pow(H,4) + 
         114*G4*G5p*G5x*pow(H,4) + 3*G3x*G4x*rptot - 6*G4px*G4x*rptot - 24*G4p*G4xx*rptot - 3*G3x*G5p*rptot + 6*G4px*G5p*rptot + 15*G4p*G5px*rptot + 12*G4x*G5x*pow(H,2)*rptot - 
         12*G5p*G5x*pow(H,2)*rptot + G3p*(26*G3x*G4 - 92*G4*G4px - 32*G4p*G4x + 30*G4p*G5p + 12*G4*G5pp + 42*G4*G5x*pow(H,2) + 3*G5x*rptot)))))/
(pow(a,2)*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H)*
 (-2*a*pow(dphi,7)*(6*G4xxx*G5x - 3*G5pxx*G5x - 12*G4xx*G5xx + 6*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,3) + pow(dphi,8)*(3*pow(G5xx,2) - 2*G5x*G5xxx)*pow(H,4) + 
   4*pow(a,8)*(G2x*G4 - 2*G3p*G4 + 3*pow(G4p,2) + 6*G4*G4x*pow(H,2) - 6*G4*G5p*pow(H,2)) + 24*pow(a,7)*dphi*H*(G3x*G4 - 3*G4*G4px - 2*G4p*G4x + 2*G4p*G5p + G4*G5x*pow(H,2)) + 
   2*pow(a,6)*pow(dphi,2)*(2*G2xx*G4 - 2*G3px*G4 - 6*G3x*G4p + 12*G4p*G4px - 2*G2x*G4x + 4*G3p*G4x + G2x*G5p - 2*G3p*G5p + 12*pow(G4x,2)*pow(H,2) + 48*G4*G4xx*pow(H,2) - 30*G4x*G5p*pow(H,2) + 
      18*pow(G5p,2)*pow(H,2) - 30*G4*G5px*pow(H,2) - 18*G4p*G5x*pow(H,2)) + 
   2*pow(a,5)*pow(dphi,3)*H*(6*G3xx*G4 - 12*G4*G4pxx + 12*G4px*G4x - 24*G4p*G4xx - 6*G3x*G5p + 6*G4px*G5p + 12*G4p*G5px - G2x*G5x + 2*G3p*G5x + 18*G4x*G5x*pow(H,2) - 24*G5p*G5x*pow(H,2) + 
      14*G4*G5xx*pow(H,2)) + 2*pow(a,2)*pow(dphi,6)*pow(H,2)*(24*pow(G4xx,2) + 6*G4xxx*G5p - 24*G4xx*G5px + 6*pow(G5px,2) - 3*G5p*G5pxx + 6*G4x*(-2*G4xxx + G5pxx) - 3*G3xx*G5x + 6*G4pxx*G5x + 
      3*G3x*G5xx - 6*G4px*G5xx + 2*G5x*G5xx*pow(H,2)) - 2*pow(a,3)*pow(dphi,5)*H*
    (-12*G3x*G4xx + 24*G4px*G4xx + G3xx*(6*G4x - 3*G5p) + 6*G4pxx*(-2*G4x + G5p) + 6*G3x*G5px - 12*G4px*G5px + G2xx*G5x - G3px*G5x - 12*G4xx*G5x*pow(H,2) + 3*G5px*G5x*pow(H,2) + 
      2*G4x*G5xx*pow(H,2) + 5*G5p*G5xx*pow(H,2) - 2*G4*G5xxx*pow(H,2)) + pow(a,4)*pow(dphi,4)*
    (3*pow(G3x,2) - 12*G3x*G4px + 12*pow(G4px,2) - 4*G2xx*G4x + 4*G3px*G4x + 2*G2xx*G5p - 2*G3px*G5p + 24*G4*G4xxx*pow(H,2) - 48*G4xx*G5p*pow(H,2) + 12*G4x*G5px*pow(H,2) + 
      18*G5p*G5px*pow(H,2) - 12*G4*G5pxx*pow(H,2) + 6*G3x*G5x*pow(H,2) - 12*G4p*G5xx*pow(H,2) + 15*pow(G5x,2)*pow(H,4))));

   // Milestone siipp
   *siipp = (3*(-2*pow(a,4)*G4p + 4*pow(a,3)*dphi*(G4x - G5p)*H + 2*a*pow(dphi,3)*(2*G4xx - G5px)*H + pow(dphi,4)*G5xx*pow(H,2) + pow(a,2)*pow(dphi,2)*(G3x - 2*G4px + 3*G5x*pow(H,2))))/
   (a*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H));

   // Milestone sij
   *sij = (-2*pow(a,5)*G4p + ddphi*pow(dphi,3)*G5xx*H + pow(a,2)*dphi*(dH*dphi*G5x + 2*(pow(dphi,2)*(-G4xx + G5px) + ddphi*G5x)*H) - a*pow(dphi,2)*(ddphi*(-2*G4xx + G5px) + pow(dphi,2)*G5xx*pow(H,2)) + 
   pow(a,3)*(2*ddphi*(G4x - G5p) + pow(dphi,2)*(2*G4px - G5pp - G5x*pow(H,2))))/(pow(a,2)*(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H));

   // Milestone sijdot
   *sijdot = (4*pow(a,10)*dphi*(pow(G4p,2) - G4*G4pp) + pow(ddphi,2)*pow(dphi,7)*(pow(G5xx,2) - G5x*G5xxx)*pow(H,2) + 
   a*ddphi*pow(dphi,6)*H*(ddphi*(-2*G4xxx*G5x + G5pxx*G5x + 4*G4xx*G5xx - 2*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx) - 2*pow(dphi,2)*(pow(G5xx,2) - G5x*G5xxx)*pow(H,2)) + 
   pow(a,2)*pow(dphi,5)*(pow(ddphi,2)*(4*pow(G4xx,2) - 4*G4x*G4xxx + 2*G4xxx*G5p - 4*G4xx*G5px + pow(G5px,2) + 2*G4x*G5pxx - G5p*G5pxx) - dddphi*dphi*G5x*G5xx*pow(H,2) + 
      2*ddphi*pow(dphi,2)*(2*G4xxx*G5x - 2*G5pxx*G5x - 4*G4xx*G5xx + 3*G5px*G5xx + 2*G4x*G5xxx - G5p*G5xxx)*pow(H,2) + pow(dphi,4)*(pow(G5xx,2) - G5x*G5xxx)*pow(H,4)) + 
   2*pow(a,8)*(2*dddphi*G4*(G4x - G5p) - 5*dH*pow(dphi,2)*G4*G5x*H + pow(dphi,3)*
       (-4*G4p*G4px + 2*G4pp*G4x - G4pp*G5p + 2*G4p*G5pp + 4*G4p*G5x*pow(H,2) + G4*(2*G4ppx - G5ppp + (6*G4xx - 7*G5px)*pow(H,2))) + 4*ddphi*dphi*(G4p*(-G4x + G5p) + G4*(G4px - G5pp - 2*G5x*pow(H,2))))
     + 4*pow(a,9)*H*(2*ddphi*G4*(-G4x + G5p) + pow(dphi,2)*(G4p*(G4x - G5p) + G4*(-G4px + G5pp + G5x*pow(H,2)))) + 
   pow(a,6)*dphi*(4*pow(ddphi,2)*(pow(G4x,2) + 3*G4*G4xx - 2*G4x*G5p + pow(G5p,2) - 2*G4*G5px) - 2*dddphi*dphi*(2*pow(G4x,2) - 2*G4*G4xx - 3*G4x*G5p + pow(G5p,2) + G4*G5px) + 
      dH*pow(dphi,3)*(8*G4x*G5x - 3*G5p*G5x - 6*G4*G5xx)*H + pow(dphi,4)*(4*pow(G4px,2) - 4*G4ppx*G4x + 2*G4ppx*G5p - 4*G4px*G5pp + pow(G5pp,2) + 2*G4x*G5ppp - G5p*G5ppp - 8*G4x*G4xx*pow(H,2) + 
         4*G4*G4xxx*pow(H,2) + 2*G4xx*G5p*pow(H,2) + 10*G4x*G5px*pow(H,2) - 3*G5p*G5px*pow(H,2) - 6*G4*G5pxx*pow(H,2) - 6*G4px*G5x*pow(H,2) + 2*G5pp*G5x*pow(H,2) + 4*G4p*G5xx*pow(H,2) + 
         pow(G5x,2)*pow(H,4)) + 4*ddphi*pow(dphi,2)*(-2*G4p*G4xx - G4px*G5p + G4x*G5pp + G4p*G5px + 2*G4x*G5x*pow(H,2) + G4*(2*G4pxx - G5ppx - 6*G5xx*pow(H,2)))) + 
   2*pow(a,7)*(ddH*pow(dphi,2)*G4*G5x + dH*(4*ddphi*dphi*G4*G5x + pow(dphi,3)*(-2*G4*G4xx + 3*G4*G5px - 2*G4p*G5x)) + 
      H*(2*pow(ddphi,2)*G4*G5x + 2*dddphi*dphi*G4*G5x + ddphi*pow(dphi,2)*(2*pow(G4x,2) - 16*G4*G4xx - 2*G4x*G5p + 14*G4*G5px - 5*G4p*G5x) + 
         pow(dphi,4)*(-4*G4*G4pxx + 4*G4p*G4xx + G4px*G5p - G4x*G5pp + 3*G4*G5ppx - 4*G4p*G5px + G4pp*G5x - G4x*G5x*pow(H,2) + 5*G4*G5xx*pow(H,2)))) + 
   pow(a,3)*pow(dphi,4)*(pow(ddphi,2)*(4*G4xx*G5x - G5px*G5x - 6*G4x*G5xx + G5p*G5xx + 2*G4*G5xxx)*H + 
      dphi*H*(dddphi*(-2*G4xx*G5x + G5px*G5x - 2*G4x*G5xx + G5p*G5xx) + pow(dphi,2)*H*(dH*G5x*G5xx + dphi*(-2*G4xxx*G5x + 3*G5pxx*G5x + 4*G4xx*G5xx - 4*G5px*G5xx - 2*G4x*G5xxx + G5p*G5xxx)*H)) + 
      ddphi*dphi*(2*dH*(2*G4xx*G5x - G5px*G5x - 2*G4x*G5xx + G5p*G5xx) + dphi*H*(-8*pow(G4xx,2) + 8*G4x*G4xxx - 4*G4xxx*G5p + 12*G4xx*G5px - 4*pow(G5px,2) - 8*G4x*G5pxx + 4*G5p*G5pxx - 4*G4pxx*G5x + 
            2*G5ppx*G5x + 4*G4px*G5xx - 2*G5pp*G5xx + 3*G5x*G5xx*pow(H,2)))) + pow(a,5)*pow(dphi,2)*
    (ddH*pow(dphi,2)*(-2*G4x + G5p)*G5x + dH*dphi*(-4*ddphi*G4x*G5x + 4*ddphi*G4*G5xx + 
         pow(dphi,2)*(4*G4x*G4xx - 2*G4xx*G5p - 6*G4x*G5px + 3*G5p*G5px + 4*G4px*G5x - 2*G5pp*G5x + pow(G5x,2)*pow(H,2))) + 
      H*(2*dddphi*dphi*(-3*G4x*G5x + 2*G5p*G5x + G4*G5xx) + 2*pow(ddphi,2)*(3*G4x*G5x - 4*G5p*G5x + 5*G4*G5xx) + 
         ddphi*pow(dphi,2)*(20*G4x*G4xx - 8*G4*G4xxx - 4*G4xx*G5p - 18*G4x*G5px + 4*G5p*G5px + 8*G4*G5pxx + 6*G4px*G5x - G5pp*G5x - 4*G4p*G5xx - pow(G5x,2)*pow(H,2)) + 
         pow(dphi,4)*(8*G4pxx*G4x - 8*G4px*G4xx - 4*G4pxx*G5p + 4*G4xx*G5pp - 6*G4x*G5ppx + 3*G5p*G5ppx + 8*G4px*G5px - 4*G5pp*G5px - 2*G4ppx*G5x + G5ppp*G5x + 2*G4xx*G5x*pow(H,2) - 
            G5px*G5x*pow(H,2) - 8*G4x*G5xx*pow(H,2) + 3*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2)))) - 
   pow(a,4)*pow(dphi,3)*(2*pow(ddphi,2)*(-2*G4*G4xxx + G4xx*G5p + 2*G4x*(G4xx - G5px) + G4*G5pxx - 2*pow(G5x,2)*pow(H,2)) + 
      ddphi*dphi*(-(dH*pow(G5x,2)*H) + 2*dphi*(4*G4pxx*G4x - 4*G4px*G4xx - 2*G4pxx*G5p + 2*G4xx*G5pp - 2*G4x*G5ppx + G5p*G5ppx + 2*G4px*G5px - G5pp*G5px + G4xx*G5x*pow(H,2) - 9*G4x*G5xx*pow(H,2) + 
            3*G5p*G5xx*pow(H,2) + 2*G4*G5xxx*pow(H,2))) + dphi*(dddphi*(4*G4x*G4xx - 2*G4xx*G5p - 2*G4x*G5px + G5p*G5px + 2*pow(G5x,2)*pow(H,2)) + 
         dphi*(-(pow(dH,2)*pow(G5x,2)) + dH*dphi*(2*G4xx*G5x - G5px*G5x - 6*G4x*G5xx + 3*G5p*G5xx)*H + 
            H*(ddH*pow(G5x,2) + pow(dphi,2)*H*(-4*pow(G4xx,2) + 4*G4x*G4xxx - 2*G4xxx*G5p + 8*G4xx*G5px - 4*pow(G5px,2) - 6*G4x*G5pxx + 3*G5p*G5pxx - 4*G4pxx*G5x + 3*G5ppx*G5x + 4*G4px*G5xx - 
                  2*G5pp*G5xx + G5x*G5xx*pow(H,2)))))))/(pow(a,4)*pow(2*pow(a,3)*G4 + a*pow(dphi,2)*(-2*G4x + G5p) - pow(dphi,3)*G5x*H,2));

}
