#include <math.h>

/** function to compute coefficients of hubble equation*/
void horndeski_hubble_coefficients(
    double a,
    double dphi,
    double rhotot,
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
    double *E0,
    double *E1,
    double *E2,
    double *E3)
{
    double X = 0.5 * dphi * dphi / a / a;

    *E0 = G2 - rhotot + 2 * (-G2x + G3p) * X;
    *E1 = (6 * dphi * (G4p - G3x * X + 2 * G4px * X)) / a;
    *E2 = 6 * (G4 + X * (-4 * G4x + 3 * G5p - 4 * G4xx * X + 2 * G5px * X));
    *E3 = (-2 * dphi * X * (5 * G5x + 2 * G5xx * X)) / a;
}

/** function to compute the eft functions alpha and a few derivatives*/
void horndeski_eft_alphas(
    double a,
    double H,
    double dH,
    double ddH,
    double dphi,
    double ddphi,
    double dddphi,
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
    double *alphaM,
    double *dalphaM,
    double *alphaB,
    double *dalphaB,
    double *alphaK,
    double *dalphaK,
    double *alphaT,
    double *dalphaT,
    double *M,
    double *dM,
    double *ddM)
{
    double X = 0.5 * dphi * dphi / a / a;

    double M2 = 2 * (G4 - 2 * G4x * X + G5p * X - (dphi * G5x * H * X) / a);

    double dM2 = (dphi * (2 * pow(a, 5) * G4p + 2 * pow(a, 4) * dphi * (G4x - G5p) * H - ddphi * pow(dphi, 3) * G5xx * H -
                          pow(a, 2) * dphi * (dH * dphi * G5x - 2 * pow(dphi, 2) * (G4xx - G5px) * H + 3 * ddphi * G5x * H) +
                          a * (ddphi * pow(dphi, 2) * (-2 * G4xx + G5px) + pow(dphi, 4) * G5xx * pow(H, 2)) +
                          pow(a, 3) * (2 * ddphi * (-G4x + G5p) + pow(dphi, 2) * (-2 * G4px + G5pp + 3 * G5x * pow(H, 2))))) /
                 pow(a, 5);

    double ddM2 = (-(pow(ddphi, 2) * pow(dphi, 5) * G5xxx * H) + a * ddphi * pow(dphi, 4) * (ddphi * (-2 * G4xxx + G5pxx) + 2 * pow(dphi, 2) * G5xxx * pow(H, 2)) -
                   pow(a, 2) * pow(dphi, 3) * (7 * pow(ddphi, 2) * G5xx * H + dddphi * dphi * G5xx * H + pow(dphi, 4) * G5xxx * pow(H, 3) + 2 * ddphi * dphi * (dH * G5xx + 2 * dphi * (-G4xxx + G5pxx) * H)) +
                   2 * pow(a, 7) * (ddphi * G4p + pow(dphi, 2) * (G4pp + (-G4x + G5p) * pow(H, 2))) +
                   2 * pow(a, 6) * dphi * (dH * dphi * (G4x - G5p) + 4 * ddphi * (G4x - G5p) * H + pow(dphi, 2) * H * (2 * G4px - 2 * G5pp - 3 * G5x * pow(H, 2))) +
                   pow(a, 5) * (2 * pow(ddphi, 2) * (-G4x + G5p) + 2 * dddphi * dphi * (-G4x + G5p) + 9 * dH * pow(dphi, 3) * G5x * H +
                                pow(dphi, 4) * (-2 * G4ppx + G5ppp + (-8 * G4xx + 11 * G5px) * pow(H, 2)) +
                                ddphi * pow(dphi, 2) * (-6 * G4px + 5 * G5pp + 18 * G5x * pow(H, 2))) +
                   pow(a, 3) * pow(dphi, 2) * (dddphi * dphi * (-2 * G4xx + G5px) + pow(ddphi, 2) * (-8 * G4xx + 5 * G5px) + pow(dphi, 3) * H * (3 * dH * G5xx - 2 * dphi * G4xxx * H + 3 * dphi * G5pxx * H) + 2 * ddphi * pow(dphi, 2) * (-2 * G4pxx + G5ppx + 8 * G5xx * pow(H, 2))) -
                   pow(a, 4) * dphi * (ddH * pow(dphi, 2) * G5x + dH * (pow(dphi, 3) * (-2 * G4xx + 3 * G5px) + 6 * ddphi * dphi * G5x) + H * (ddphi * pow(dphi, 2) * (-20 * G4xx + 19 * G5px) + 6 * pow(ddphi, 2) * G5x + 3 * dddphi * dphi * G5x + pow(dphi, 4) * (-4 * G4pxx + 3 * G5ppx + 7 * G5xx * pow(H, 2))))) /
                  pow(a, 7);

    double aM = (dphi * (2 * pow(a, 5) * G4p + 2 * pow(a, 4) * dphi * (G4x - G5p) * H - ddphi * pow(dphi, 3) * G5xx * H -
                         pow(a, 2) * dphi * (dH * dphi * G5x - 2 * pow(dphi, 2) * (G4xx - G5px) * H + 3 * ddphi * G5x * H) +
                         a * (ddphi * pow(dphi, 2) * (-2 * G4xx + G5px) + pow(dphi, 4) * G5xx * pow(H, 2)) +
                         pow(a, 3) * (2 * ddphi * (-G4x + G5p) + pow(dphi, 2) * (-2 * G4px + G5pp + 3 * G5x * pow(H, 2))))) /
                (pow(a, 3) * H * (2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H));

    double daM = (-4 * pow(a, 11) * dphi * G4 * G4p * pow(H, 2) + pow(ddphi, 2) * pow(dphi, 8) * (-pow(G5xx, 2) + G5x * G5xxx) * pow(H, 3) +
                  a * ddphi * pow(dphi, 7) * pow(H, 2) * (ddphi * (2 * G4xxx * G5x - G5pxx * G5x - 4 * G4xx * G5xx + 2 * G5px * G5xx + 2 * G4x * G5xxx - G5p * G5xxx) + 2 * pow(dphi, 2) * (pow(G5xx, 2) - G5x * G5xxx) * pow(H, 2)) +
                  pow(a, 2) * pow(dphi, 6) * H * (dddphi * dphi * G5x * G5xx * pow(H, 2) + pow(dphi, 4) * (-pow(G5xx, 2) + G5x * G5xxx) * pow(H, 4) - ddphi * dphi * H * (dH * G5x * G5xx + 2 * dphi * (2 * G4xxx * G5x - 2 * G5pxx * G5x - 4 * G4xx * G5xx + 3 * G5px * G5xx + 2 * G4x * G5xxx - G5p * G5xxx) * H) + pow(ddphi, 2) * (-4 * pow(G4xx, 2) + 4 * G4x * G4xxx - 2 * G4xxx * G5p + 4 * G4xx * G5px - pow(G5px, 2) - 2 * G4x * G5pxx + G5p * G5pxx + G5x * G5xx * pow(H, 2))) - 4 * pow(a, 10) * (dH * dphi * G4 * G4p - ddphi * G4 * G4p * H + pow(dphi, 2) * H * (pow(G4p, 2) - G4 * (G4pp + 2 * (-G4x + G5p) * pow(H, 2)))) +
                  2 * pow(a, 9) * dphi * pow(H, 2) * (10 * ddphi * G4 * (G4x - G5p) + pow(dphi, 2) * (-2 * G4p * G4x + 3 * G4p * G5p + G4 * (6 * G4px - 5 * G5pp - 9 * G5x * pow(H, 2)))) +
                  pow(a, 6) * pow(dphi, 2) * (dH * (-2 * ddphi * dphi * (2 * pow(G4x, 2) - 2 * G4 * G4xx - 3 * G4x * G5p + pow(G5p, 2) + G4 * G5px) + pow(dphi, 3) * (-4 * G4px * G4x + 2 * G4px * G5p + 2 * G4x * G5pp - G5p * G5pp - 10 * G4x * G5x * pow(H, 2) + 3 * G5p * G5x * pow(H, 2) + 4 * G4 * G5xx * pow(H, 2))) + H * (-2 * pow(ddphi, 2) * (8 * G4 * G4xx - G4x * G5p + pow(G5p, 2) - 5 * G4 * G5px) + 2 * dddphi * dphi * (2 * pow(G4x, 2) - 2 * G4 * G4xx - 3 * G4x * G5p + pow(G5p, 2) + G4 * G5px) + ddphi * pow(dphi, 2) * (-8 * G4 * G4pxx + 4 * G4px * G4x + 8 * G4p * G4xx + 2 * G4px * G5p - 6 * G4x * G5pp + G5p * G5pp + 4 * G4 * G5ppx - 4 * G4p * G5px - 28 * G4x * G5x * pow(H, 2) + 7 * G5p * G5x * pow(H, 2) + 34 * G4 * G5xx * pow(H, 2)) - pow(dphi, 4) * (4 * pow(G4px, 2) - 4 * G4ppx * G4x + 2 * G4ppx * G5p - 4 * G4px * G5pp + pow(G5pp, 2) + 2 * G4x * G5ppp - G5p * G5ppp - 12 * G4x * G4xx * pow(H, 2) + 4 * G4 * G4xxx * pow(H, 2) + 2 * G4xx * G5p * pow(H, 2) + 18 * G4x * G5px * pow(H, 2) - 5 * G5p * G5px * pow(H, 2) - 6 * G4 * G5pxx * pow(H, 2) - 6 * G4px * G5x * pow(H, 2) + G5pp * G5x * pow(H, 2) + 4 * G4p * G5xx * pow(H, 2)))) - pow(a, 3) * pow(dphi, 5) * H * (dddphi * dphi * (-2 * G4xx * G5x + G5px * G5x - 2 * G4x * G5xx + G5p * G5xx) * H + pow(ddphi, 2) * (4 * G4xx * G5x - G5px * G5x - 10 * G4x * G5xx + 3 * G5p * G5xx + 2 * G4 * G5xxx) * H + pow(dphi, 4) * (-2 * G4xxx * G5x + 3 * G5pxx * G5x + 4 * G4xx * G5xx - 4 * G5px * G5xx - 2 * G4x * G5xxx + G5p * G5xxx) * pow(H, 3) + ddphi * dphi * (dH * (6 * G4xx * G5x - 3 * G5px * G5x - 2 * G4x * G5xx + G5p * G5xx) + dphi * H * (-8 * pow(G4xx, 2) + 8 * G4x * G4xxx - 4 * G4xxx * G5p + 12 * G4xx * G5px - 4 * pow(G5px, 2) - 8 * G4x * G5pxx + 4 * G5p * G5pxx - 4 * G4pxx * G5x + 2 * G5ppx * G5x + 4 * G4px * G5xx - 2 * G5pp * G5xx + 5 * G5x * G5xx * pow(H, 2)))) +
                  2 * pow(a, 8) * (dH * dphi * (2 * ddphi * G4 * (G4x - G5p) + pow(dphi, 2) * (2 * G4 * G4px + 2 * G4p * G4x - G4p * G5p - G4 * G5pp + 7 * G4 * G5x * pow(H, 2))) + H * (2 * pow(ddphi, 2) * G4 * (-G4x + G5p) + 2 * dddphi * dphi * G4 * (-G4x + G5p) + ddphi * pow(dphi, 2) * (-6 * G4 * G4px + 2 * G4p * G4x - 3 * G4p * G5p + 5 * G4 * G5pp + 21 * G4 * G5x * pow(H, 2)) + pow(dphi, 4) * (4 * G4p * G4px - 2 * G4pp * G4x + G4pp * G5p - 2 * G4p * G5pp + 2 * pow(G4x, 2) * pow(H, 2) - 2 * G4x * G5p * pow(H, 2) - 5 * G4p * G5x * pow(H, 2) + G4 * (-2 * G4ppx + G5ppp + (-10 * G4xx + 13 * G5px) * pow(H, 2))))) -
                  pow(a, 7) * dphi * (-2 * pow(dH, 2) * pow(dphi, 2) * G4 * G5x + 2 * dH * dphi * (3 * ddphi * G4 * G5x + pow(dphi, 2) * (G4 * G5px - 3 * G4p * G5x)) * H + H * (2 * ddH * pow(dphi, 2) * G4 * G5x + H * (12 * pow(ddphi, 2) * G4 * G5x + 6 * dddphi * dphi * G4 * G5x + 2 * ddphi * pow(dphi, 2) * (6 * pow(G4x, 2) - 22 * G4 * G4xx - 7 * G4x * G5p + pow(G5p, 2) + 20 * G4 * G5px - 5 * G4p * G5x) + pow(dphi, 4) * (-8 * G4 * G4pxx + 4 * G4px * G4x + 8 * G4p * G4xx + 2 * G4px * G5p - 6 * G4x * G5pp + G5p * G5pp + 6 * G4 * G5ppx - 8 * G4p * G5px + 2 * G4pp * G5x - 10 * G4x * G5x * pow(H, 2) + G5p * G5x * pow(H, 2) + 16 * G4 * G5xx * pow(H, 2))))) -
                  pow(a, 5) * pow(dphi, 3) * (pow(dH, 2) * pow(dphi, 2) * (2 * G4x - G5p) * G5x + dH * dphi * H * (-3 * ddphi * G5p * G5x + 2 * ddphi * G4 * G5xx + pow(dphi, 2) * (-2 * G4x * G5px + G5p * G5px + G5x * (6 * G4px - 3 * G5pp + G5x * pow(H, 2)))) + H * (ddH * pow(dphi, 2) * (-2 * G4x + G5p) * G5x + H * (-2 * pow(ddphi, 2) * (G4x * G5x + 2 * G5p * G5x - 7 * G4 * G5xx) + dddphi * dphi * (-8 * G4x * G5x + 5 * G5p * G5x + 2 * G4 * G5xx) + ddphi * pow(dphi, 2) * (-8 * G4 * G4xxx - 6 * G4xx * G5p + 28 * G4x * (G4xx - G5px) + 8 * G5p * G5px + 8 * G4 * G5pxx + 6 * G4px * G5x - G5pp * G5x - 4 * G4p * G5xx + 3 * pow(G5x, 2) * pow(H, 2)) + pow(dphi, 4) * (8 * G4pxx * G4x - 8 * G4px * G4xx - 4 * G4pxx * G5p + 4 * G4xx * G5pp - 6 * G4x * G5ppx + 3 * G5p * G5ppx + 8 * G4px * G5px - 4 * G5pp * G5px - 2 * G4ppx * G5x + G5ppp * G5x + 2 * G4xx * G5x * pow(H, 2) + G5px * G5x * pow(H, 2) - 12 * G4x * G5xx * pow(H, 2) + 4 * G5p * G5xx * pow(H, 2) + 2 * G4 * G5xxx * pow(H, 2))))) +
                  pow(a, 4) * pow(dphi, 4) * (pow(ddphi, 2) * H * (8 * G4x * G4xx - 4 * G4 * G4xxx - 6 * G4x * G5px + G5p * G5px + 2 * G4 * G5pxx - 3 * pow(G5x, 2) * pow(H, 2)) + ddphi * dphi * (-(dH * (4 * G4x * G4xx - 2 * G4xx * G5p - 2 * G4x * G5px + G5p * G5px + 3 * pow(G5x, 2) * pow(H, 2))) + dphi * H * (8 * G4pxx * G4x - 8 * G4px * G4xx - 4 * G4pxx * G5p + 4 * G4xx * G5pp - 4 * G4x * G5ppx + 2 * G5p * G5ppx + 4 * G4px * G5px - 2 * G5pp * G5px + 2 * G4xx * G5x * pow(H, 2) + 2 * G5px * G5x * pow(H, 2) - 26 * G4x * G5xx * pow(H, 2) + 9 * G5p * G5xx * pow(H, 2) + 4 * G4 * G5xxx * pow(H, 2))) + dphi * H * (dddphi * (4 * G4x * G4xx - 2 * G4xx * G5p - 2 * G4x * G5px + G5p * G5px + 3 * pow(G5x, 2) * pow(H, 2)) + dphi * (-2 * pow(dH, 2) * pow(G5x, 2) + dH * dphi * (4 * G4xx * G5x - 3 * G5px * G5x - 4 * G4x * G5xx + 2 * G5p * G5xx) * H + H * (ddH * pow(G5x, 2) + pow(dphi, 2) * H * (-4 * pow(G4xx, 2) + 4 * G4x * G4xxx - 2 * G4xxx * G5p + 8 * G4xx * G5px - 4 * pow(G5px, 2) - 6 * G4x * G5pxx + 3 * G5p * G5pxx - 4 * G4pxx * G5x + 3 * G5ppx * G5x + 4 * G4px * G5xx - 2 * G5pp * G5xx + 2 * G5x * G5xx * pow(H, 2))))))) /
                 (pow(a, 5) * pow(H, 2) * pow(2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H, 2));

    double aB = (dphi * (-2 * pow(a, 4) * G4p + 4 * pow(a, 3) * dphi * (G4x - G5p) * H + 2 * a * pow(dphi, 3) * (2 * G4xx - G5px) * H + pow(dphi, 4) * G5xx * pow(H, 2) +
                         pow(a, 2) * pow(dphi, 2) * (G3x - 2 * G4px + 3 * G5x * pow(H, 2)))) /
                (pow(a, 2) * H * (2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H));

    double daB = (4 * pow(a, 10) * dphi * G4 * G4p * pow(H, 2) + ddphi * pow(dphi, 9) * (pow(G5xx, 2) - G5x * G5xxx) * pow(H, 4) +
                  a * pow(dphi, 8) * pow(H, 3) * (ddphi * (-4 * G4xxx * G5x + 2 * G5pxx * G5x + 6 * G4xx * G5xx - 3 * G5px * G5xx - 2 * G4x * G5xxx + G5p * G5xxx) + pow(dphi, 2) * (-pow(G5xx, 2) + G5x * G5xxx) * pow(H, 2)) -
                  2 * pow(a, 8) * dphi * pow(H, 2) * (8 * ddphi * G4 * (-G4x + G5p) + pow(dphi, 2) * (3 * G3x * G4 - 12 * G4 * G4px + 4 * G4p * G4x - 3 * G4p * G5p + 4 * G4 * G5pp + 9 * G4 * G5x * pow(H, 2))) -
                  2 * pow(a, 7) * pow(dphi, 2) * (dH * dphi * (G3x * G4 - 2 * G4 * G4px + 2 * G4p * G4x - G4p * G5p - 3 * G4 * G5x * pow(H, 2)) - ddphi * H * (3 * G3x * G4 - 8 * G4 * G4px + G4p * G5p + 9 * G4 * G5x * pow(H, 2)) + pow(dphi, 2) * H * (-(G3px * G4) + G3x * G4p + 2 * G4 * G4ppx - 2 * G4pp * G4x + G4pp * G5p - G4p * G5pp - 4 * pow(G4x, 2) * pow(H, 2) + 20 * G4 * G4xx * pow(H, 2) + 4 * G4x * G5p * pow(H, 2) - 15 * G4 * G5px * pow(H, 2) + G4p * G5x * pow(H, 2))) +
                  pow(a, 2) * (pow(dphi, 9) * (4 * G4xxx * G5x - 3 * G5pxx * G5x - 6 * G4xx * G5xx + 4 * G5px * G5xx + 2 * G4x * G5xxx - G5p * G5xxx) * pow(H, 4) +
                               ddphi * pow(dphi, 7) * pow(H, 2) * (8 * pow(G4xx, 2) - 8 * G4x * G4xxx + 4 * G4xxx * G5p - 8 * G4xx * G5px + 2 * pow(G5px, 2) + 4 * G4x * G5pxx - 2 * G5p * G5pxx - G3xx * G5x + 2 * G4pxx * G5x + G3x * G5xx - 2 * G4px * G5xx - 2 * G5x * G5xx * pow(H, 2))) +
                  pow(a, 4) * pow(dphi, 5) * H * (2 * dH * dphi * (G3x - 2 * G4px) * G5x + 2 * ddphi * (-12 * G4x * G4xx + 4 * G4 * G4xxx + 2 * G4xx * G5p + 8 * G4x * G5px - 2 * G5p * G5px - 2 * G4 * G5pxx + G4px * G5x - G4p * G5xx) * H + pow(dphi, 2) * H * (-2 * G3x * G4xx + 12 * G4px * G4xx + G3xx * (2 * G4x - G5p) + 6 * G4pxx * (-2 * G4x + G5p) - 4 * G4xx * G5pp + 4 * G4x * G5ppx - 2 * G5p * G5ppx + 2 * G3x * G5px - 8 * G4px * G5px + 2 * G5pp * G5px - G3px * G5x + 2 * G4ppx * G5x + 2 * G4xx * G5x * pow(H, 2) - 3 * G5px * G5x * pow(H, 2) + 10 * G4x * G5xx * pow(H, 2) - 2 * G5p * G5xx * pow(H, 2) - 2 * G4 * G5xxx * pow(H, 2))) +
                  pow(a, 5) * pow(dphi, 4) * (pow(dphi, 2) * H * (2 * G3x * G4px - 4 * pow(G4px, 2) - 2 * G3px * G4x + 4 * G4ppx * G4x + G3px * G5p - 2 * G4ppx * G5p - G3x * G5pp + 2 * G4px * G5pp + 24 * G4x * G4xx * pow(H, 2) - 8 * G4 * G4xxx * pow(H, 2) - 4 * G4xx * G5p * pow(H, 2) - 18 * G4x * G5px * pow(H, 2) + 3 * G5p * G5px * pow(H, 2) + 6 * G4 * G5pxx * pow(H, 2) + G5pp * G5x * pow(H, 2)) + dH * dphi * (G3x * (2 * G4x - G5p) + G4px * (-4 * G4x + 2 * G5p) - (2 * G4x * G5x + G5p * G5x - 2 * G4 * G5xx) * pow(H, 2)) + ddphi * H * (2 * G3xx * G4 - 4 * G3x * G4x + 12 * G4px * G4x - 4 * G4p * G4xx + G3x * G5p - 4 * G4px * G5p + 2 * G4p * G5px - 8 * G4x * G5x * pow(H, 2) - G5p * G5x * pow(H, 2) - 4 * G4 * (G4pxx - 4 * G5xx * pow(H, 2)))) -
                  pow(a, 6) * pow(dphi, 3) * H * (4 * dH * dphi * G4p * G5x + 4 * ddphi * (2 * pow(G4x, 2) - 10 * G4 * G4xx - 2 * G4x * G5p + 6 * G4 * G5px + G4p * G5x) * H + pow(dphi, 2) * H * (2 * G3xx * G4 - 4 * G3x * G4x + 12 * G4px * G4x + 4 * G4p * G4xx + G3x * G5p - 4 * G4x * G5pp - 2 * G4pp * G5x - 8 * G4x * G5x * pow(H, 2) - G5p * G5x * pow(H, 2) + 4 * G4 * (-3 * G4pxx + G5ppx + 4 * G5xx * pow(H, 2)))) +
                  pow(a, 3) * pow(dphi, 6) * H * (ddphi * (4 * G4pxx * G4x + 2 * G3x * G4xx - 4 * G4px * G4xx - 2 * G4pxx * G5p + G3xx * (-2 * G4x + G5p) - G3x * G5px + 2 * G4px * G5px - 2 * G4xx * G5x * pow(H, 2) + 3 * G5px * G5x * pow(H, 2) - 10 * G4x * G5xx * pow(H, 2) + 2 * G5p * G5xx * pow(H, 2) + 2 * G4 * G5xxx * pow(H, 2)) + dphi * H * (dH * (4 * G4xx * G5x - 2 * G5px * G5x - 2 * G4x * G5xx + G5p * G5xx) + dphi * H * (-8 * pow(G4xx, 2) + 8 * G4x * G4xxx - 4 * G4xxx * G5p + 12 * G4xx * G5px - 4 * pow(G5px, 2) - 6 * G4x * G5pxx + 3 * G5p * G5pxx + G3xx * G5x - 6 * G4pxx * G5x + 2 * G5ppx * G5x - G3x * G5xx + 4 * G4px * G5xx - G5pp * G5xx + 2 * G5x * G5xx * pow(H, 2)))) +
                  4 * pow(a, 9) * (dH * dphi * G4 * G4p - H * (ddphi * G4 * G4p + pow(dphi, 2) * (-pow(G4p, 2) + G4 * (G4pp + 4 * (G4x - G5p) * pow(H, 2)))))) /
                 (pow(a, 4) * pow(H, 2) * pow(2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H, 2));

    double aK = (pow(dphi, 2) * (3 * a * pow(dphi, 4) * (2 * G4xxx - G5pxx) * pow(H, 2) + pow(dphi, 5) * G5xxx * pow(H, 3) +
                                 pow(a, 3) * pow(dphi, 2) * (G2xx - G3px + 3 * (8 * G4xx - 5 * G5px) * pow(H, 2)) +
                                 6 * pow(a, 4) * dphi * H * (G3x - 3 * G4px + G5x * pow(H, 2)) + pow(a, 2) * pow(dphi, 3) * H * (3 * G3xx - 6 * G4pxx + 7 * G5xx * pow(H, 2)) +
                                 pow(a, 5) * (G2x - 2 * (G3p + 3 * (-G4x + G5p) * pow(H, 2))))) /
                (pow(a, 4) * pow(H, 2) * (2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H));

    double daK = (dphi * (-2 * pow(a, 2) * dH * dphi * (2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H) *
                              (3 * a * pow(dphi, 4) * (2 * G4xxx - G5pxx) * pow(H, 2) + pow(dphi, 5) * G5xxx * pow(H, 3) +
                               pow(a, 3) * pow(dphi, 2) * (G2xx - G3px + 3 * (8 * G4xx - 5 * G5px) * pow(H, 2)) +
                               6 * pow(a, 4) * dphi * H * (G3x - 3 * G4px + G5x * pow(H, 2)) +
                               pow(a, 2) * pow(dphi, 3) * H * (3 * G3xx - 6 * G4pxx + 7 * G5xx * pow(H, 2)) +
                               pow(a, 5) * (G2x - 2 * (G3p + 3 * (-G4x + G5p) * pow(H, 2)))) -
                          pow(dphi, 2) * H * (3 * a * pow(dphi, 4) * (2 * G4xxx - G5pxx) * pow(H, 2) + pow(dphi, 5) * G5xxx * pow(H, 3) + pow(a, 3) * pow(dphi, 2) * (G2xx - G3px + 3 * (8 * G4xx - 5 * G5px) * pow(H, 2)) + 6 * pow(a, 4) * dphi * H * (G3x - 3 * G4px + G5x * pow(H, 2)) + pow(a, 2) * pow(dphi, 3) * H * (3 * G3xx - 6 * G4pxx + 7 * G5xx * pow(H, 2)) + pow(a, 5) * (G2x - 2 * (G3p + 3 * (-G4x + G5p) * pow(H, 2)))) *
                              (2 * pow(a, 5) * G4p + 2 * pow(a, 4) * dphi * (G4x - G5p) * H - ddphi * pow(dphi, 3) * G5xx * H -
                               pow(a, 2) * dphi * (dH * dphi * G5x - 2 * pow(dphi, 2) * (G4xx - G5px) * H + 3 * ddphi * G5x * H) +
                               a * (ddphi * pow(dphi, 2) * (-2 * G4xx + G5px) + pow(dphi, 4) * G5xx * pow(H, 2)) +
                               pow(a, 3) * (2 * ddphi * (-G4x + G5p) + pow(dphi, 2) * (-2 * G4px + G5pp + 3 * G5x * pow(H, 2)))) +
                          H * (-2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (2 * G4x - G5p) + pow(dphi, 3) * G5x * H) *
                              (-(ddphi * pow(dphi, 7) * G5xxxx * pow(H, 3)) + a * pow(dphi, 6) * pow(H, 2) * (-6 * ddphi * G4xxxx + 3 * ddphi * G5pxxx + pow(dphi, 2) * G5xxxx * pow(H, 2)) +
                               2 * pow(a, 8) * dphi * H * (G2x - 2 * (G3p + 3 * (-G4x + G5p) * pow(H, 2))) +
                               pow(a, 4) * pow(dphi, 3) * (pow(dphi, 2) * H * (G2xxx - 4 * G3pxx + 6 * G4ppxx + 60 * G4xxx * pow(H, 2) - 40 * G5pxx * pow(H, 2)) + ddphi * H * (-21 * G3xx + 48 * G4pxx - 41 * G5xx * pow(H, 2)) - 3 * dH * dphi * (G3xx - 2 * G4pxx + 7 * G5xx * pow(H, 2))) -
                               pow(a, 2) * pow(dphi, 5) * H * (dphi * H * (3 * dH * G5xxx - 6 * dphi * G4xxxx * H + 4 * dphi * G5pxxx * H) + ddphi * (3 * G3xxx - 6 * G4pxxx + 14 * G5xxx * pow(H, 2))) +
                               pow(a, 6) * dphi * (-18 * ddphi * H * (G3x - 3 * G4px + G5x * pow(H, 2)) - 6 * dH * dphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2)) + pow(dphi, 2) * H * (5 * G2xx + 6 * (-2 * G3px + 3 * G4ppx + 17 * G4xx * pow(H, 2) - 12 * G5px * pow(H, 2)))) +
                               pow(a, 7) * (-12 * dH * dphi * (G4x - G5p) * H - 2 * ddphi * (G2x - 2 * (G3p + 3 * (-G4x + G5p) * pow(H, 2))) +
                                            pow(dphi, 2) * (-G2px + 2 * G3pp + 6 * pow(H, 2) * (3 * G3x - 10 * G4px + G5pp + 3 * G5x * pow(H, 2)))) +
                               pow(a, 5) * pow(dphi, 2) * (-6 * dH * dphi * (8 * G4xx - 5 * G5px) * H + ddphi * (-5 * G2xx + 6 * (G3px + (-17 * G4xx + 11 * G5px) * pow(H, 2))) + pow(dphi, 2) * (-G2pxx + G3ppx + pow(H, 2) * (21 * G3xx - 72 * G4pxx + 15 * G5ppx + 41 * G5xx * pow(H, 2)))) +
                               pow(a, 3) * pow(dphi, 4) * (ddphi * (-G2xxx + G3pxx - 60 * G4xxx * pow(H, 2) + 33 * G5pxx * pow(H, 2)) + dphi * H * (6 * dH * (-2 * G4xxx + G5pxx) + dphi * H * (3 * G3xxx - 12 * G4pxxx + 3 * G5ppxx + 14 * G5xxx * pow(H, 2))))))) /
                 (4. * pow(a, 12) * pow(H, 3) * pow(G4 - (pow(dphi, 2) * (2 * a * G4x - a * G5p + dphi * G5x * H)) / (2. * pow(a, 3)), 2));

    double aT = (pow(dphi, 2) * (2 * pow(a, 2) * (G4x - G5p) - ddphi * G5x + 2 * a * dphi * G5x * H)) /
                (a * (2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H));

    double daT = (dphi * (pow(ddphi, 2) * pow(dphi, 4) * (-2 * G4xx * G5x + G5px * G5x + 2 * G4x * G5xx - G5p * G5xx) + 8 * pow(a, 7) * dphi * G4 * (-G4x + G5p) * H +
                          4 * pow(a, 5) * dphi * (dH * dphi * G4 * G5x + 5 * ddphi * G4 * G5x * H + pow(dphi, 2) * (pow(G4x, 2) - G4 * G4xx - G4x * G5p + 2 * G4 * G5px - G4p * G5x) * H) +
                          2 * pow(a, 3) * pow(dphi, 3) * (-(dH * dphi * G4x * G5x) + pow(dphi, 2) * (G4xx * G5p - 2 * G4x * G5px + G4px * G5x) * H + ddphi * (-6 * G4x * G5x + G5p * G5x + 3 * G4 * G5xx) * H) -
                          a * pow(dphi, 3) * (pow(ddphi, 2) * pow(G5x, 2) * H - dddphi * dphi * pow(G5x, 2) * H + ddphi * dphi * (dH * pow(G5x, 2) + dphi * (-4 * G4xx * G5x + G5px * G5x + 4 * G4x * G5xx - G5p * G5xx) * H)) -
                          pow(a, 2) * pow(dphi, 2) * (dddphi * dphi * (-2 * G4x + G5p) * G5x + pow(ddphi, 2) * (-2 * G4x * G5x + 2 * G4 * G5xx) + 2 * pow(dphi, 4) * (G4xx * G5x - G4x * G5xx) * pow(H, 2) + ddphi * pow(dphi, 2) * (2 * G4xx * G5p - 4 * G4x * G5px + G5p * G5px + 2 * G4px * G5x - G5pp * G5x + pow(G5x, 2) * pow(H, 2))) +
                          4 * pow(a, 6) * (2 * ddphi * G4 * (G4x - G5p) + pow(dphi, 2) * (G4p * (-G4x + G5p) + G4 * (G4px - G5pp - 3 * G5x * pow(H, 2)))) -
                          2 * pow(a, 4) * (2 * pow(ddphi, 2) * G4 * G5x + dddphi * dphi * G4 * G5x + ddphi * pow(dphi, 2) * (2 * pow(G4x, 2) - 2 * G4 * G4xx - 2 * G4x * G5p + 3 * G4 * G5px - G4p * G5x) + pow(dphi, 4) * (G4px * G5p + 2 * G4 * G5xx * pow(H, 2) - G4x * (G5pp + 3 * G5x * pow(H, 2)))))) /
                 (pow(a, 2) * pow(2 * pow(a, 3) * G4 + a * pow(dphi, 2) * (-2 * G4x + G5p) - pow(dphi, 3) * G5x * H, 2));

    *M = M2;
    *dM = dM2;
    *ddM = ddM2;
    *alphaM = aM;
    *dalphaM = daM;
    *alphaB = aB;
    *dalphaB = daB;
    *alphaK = aK;
    *dalphaK = daK;
    *alphaT = aT;
    *dalphaT = daT;
}

/** function to compute coefficients of background equations*/
void horndeski_backeq_coefficients(
    double a,
    double H,
    double dphi,
    double rptot,
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
    double *cp1,
    double *ch1,
    double *cs1,
    double *cp2,
    double *ch2,
    double *cs2)
{
    double X = 0.5 * dphi * dphi / a / a;

    *cp1 = (-2 * a * G4p + 4 * dphi * H * (G4x - G5p + 2 * G4xx * X - G5px * X) +
            2 * a * X * (G3x - 2 * G4px + pow(H, 2) * (3 * G5x + 2 * G5xx * X))) /
           (3. * pow(a, 2));

    *ch1 = (4 * dphi * G5x * H * X) / (3. * a) - (4 * (G4 + (-2 * G4x + G5p) * X)) / 3.;

    *cs1 = (4 * dphi * H * (G4p - X * (2 * G3x - 6 * G4px + G5pp + 3 * G5x * pow(H, 2) + 2 * G5xx * pow(H, 2) * X)) -
            a * (rptot + 2 * X * (G2x - 2 * (G3p - G4pp + pow(H, 2) * (-5 * G4x + 5 * G5p - 10 * G4xx * X + 6 * G5px * X))))) /
           3.;

    *cp2 = (6 * dphi * (G3x - 3 * G4px + (G3xx - 2 * G4pxx) * X) + (a * (G2x - 2 * (G3p + (-G2xx + G3px) * X))) / H +
            6 * a * H * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X)) +
            2 * dphi * pow(H, 2) * (3 * G5x + X * (7 * G5xx + 2 * G5xxx * X))) /
           pow(a, 3);

    *ch2 = (-6 * G4p + 6 * G3x * X - 12 * G4px * X) / (a * H) +
           (dphi * (12 * (G4x - G5p) + 12 * (2 * G4xx - G5px) * X)) / pow(a, 2) +
           (H * (18 * G5x * X + 12 * G5xx * pow(X, 2))) / a;

    *cs2 = (2 * dphi * (G2x - 2 * G3p + 6 * G4x * pow(H, 2) - 6 * G5p * pow(H, 2) - G2xx * X + 4 * G3px * X - 6 * G4ppx * X - 6 * G4xx * pow(H, 2) * X + 8 * G5px * pow(H, 2) * X - 12 * G4xxx * pow(H, 2) * pow(X, 2) + 8 * G5pxx * pow(H, 2) * pow(X, 2))) / a -
           (G2p + 12 * G4p * pow(H, 2) + 2 * X * (-G2px + G3pp + pow(H, 2) * (-3 * G3x + 3 * G5pp - 3 * G5x * pow(H, 2) + 6 * G3xx * X - 24 * G4pxx * X + 6 * G5ppx * X + 8 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2)))) / H;
}

/** function to compute conformal time derivative of coefficients of background equations*/
void horndeski_backeq_coefficients_p(
    double a,
    double H,
    double dH,
    double dphi,
    double ddphi,
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
    double *cp1p,
    double *ch1p,
    double *cs1p,
    double *cp2p,
    double *ch2p,
    double *cs2p)
{
    double X = 0.5 * dphi * dphi / a / a;

    *cp1p = (2 * ddphi * dphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2) + G3xx * X - 2 * G4pxx * X + 7 * G5xx * pow(H, 2) * X + 2 * G5xxx * pow(H, 2) * pow(X, 2)) -
             2 * pow(a, 2) * (-2 * dH * H * X * (3 * G5x + 2 * G5xx * X) + dphi * (G4pp + 4 * G4x * pow(H, 2) - 4 * G5p * pow(H, 2) - G3px * X + 2 * G4ppx * X + 20 * G4xx * pow(H, 2) * X - 15 * G5px * pow(H, 2) * X + 8 * G4xxx * pow(H, 2) * pow(X, 2) - 6 * G5pxx * pow(H, 2) * pow(X, 2))) -
             2 * pow(a, 3) * H * (-G4p + X * (3 * G3x - 12 * G4px + 4 * G5pp + 9 * G5x * pow(H, 2) + 2 * G3xx * X - 12 * G4pxx * X + 4 * G5ppx * X + 16 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2))) + 4 * a * (dH * dphi * (G4x - G5p + 2 * G4xx * X - G5px * X) + ddphi * H * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X)))) /
            (3. * pow(a, 3));

    *ch1p = (-4 * (2 * pow(a, 3) * H * X * (G4x - G5p + 2 * (G4xx - G5px) * X) +
                   ddphi * dphi * (-G4x + G5p + (-2 * G4xx + G5px) * X) -
                   a * X * (dH * dphi * G5x + ddphi * H * (3 * G5x + 2 * G5xx * X)) +
                   pow(a, 2) * dphi * (G4p + X * (-2 * G4px + G5pp + pow(H, 2) * (3 * G5x + 2 * G5xx * X))))) /
            (3. * pow(a, 2));

    *cs1p = (-2 * ddphi * dphi * (G2x - 2 * G3p + 2 * G4pp + 10 * G4x * pow(H, 2) - 10 * G5p * pow(H, 2) + G2xx * X - 2 * G3px * X + 2 * G4ppx * X + 50 * G4xx * pow(H, 2) * X - 34 * G5px * pow(H, 2) * X + 20 * G4xxx * pow(H, 2) * pow(X, 2) - 12 * G5pxx * pow(H, 2) * pow(X, 2)) +
             pow(a, 3) * H * (-rptot + 2 * X * (G2x + 2 * (-G3p + 3 * G4pp + 5 * G4x * pow(H, 2) - 5 * G5p * pow(H, 2) + G2xx * X - 6 * G3px * X + 14 * G4ppx * X - 2 * G5ppp * X + 40 * G4xx * pow(H, 2) * X - 34 * G5px * pow(H, 2) * X + 20 * G4xxx * pow(H, 2) * pow(X, 2) - 16 * G5pxx * pow(H, 2) * pow(X, 2)))) +
             pow(a, 2) * (-drptot - 2 * dphi * G2px * X -
                          8 * dH * H * X * (5 * G4x - 5 * G5p + 10 * G4xx * X - 6 * G5px * X) +
                          4 * dphi * X * (G3pp - G4pxx * (1 + 22 * pow(H, 2) * X) + pow(H, 2) * (4 * G3x - 19 * G4px + 7 * G5pp + 6 * G5x * pow(H, 2) + 4 * G3xx * X + 8 * G5ppx * X + 14 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2)))) +
             4 * a * (dH * dphi * (G4p - X * (2 * G3x - 6 * G4px + G5pp + 9 * G5x * pow(H, 2) + 6 * G5xx * pow(H, 2) * X)) + ddphi * H * (G4p - X * (6 * G3x - 20 * G4px + 3 * G5pp + 9 * G5x * pow(H, 2) + 4 * G3xx * X - 12 * G4pxx * X + 2 * G5ppx * X + 16 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2))))) /
            (3. * a);

    *cp2p = (ddphi * dphi * H * (3 * G2xx - 4 * G3px + 54 * G4xx * pow(H, 2) - 36 * G5px * pow(H, 2) + 2 * G2xxx * X - 2 * G3pxx * X + 96 * G4xxx * pow(H, 2) * X - 54 * G5pxx * pow(H, 2) * X + 24 * G4xxxx * pow(H, 2) * pow(X, 2) - 12 * G5pxxx * pow(H, 2) * pow(X, 2)) -
             2 * pow(a, 3) * pow(H, 2) * (G2x - 2 * G3p + 6 * G4x * pow(H, 2) - 6 * G5p * pow(H, 2) + 5 * G2xx * X - 12 * G3px * X + 18 * G4ppx * X + 102 * G4xx * pow(H, 2) * X - 72 * G5px * pow(H, 2) * X + 2 * G2xxx * pow(X, 2) - 8 * G3pxx * pow(X, 2) + 12 * G4ppxx * pow(X, 2) + 120 * G4xxx * pow(H, 2) * pow(X, 2) - 80 * G5pxx * pow(H, 2) * pow(X, 2) + 24 * G4xxxx * pow(H, 2) * pow(X, 3) - 16 * G5pxxx * pow(H, 2) * pow(X, 3)) +
             2 * a * pow(H, 2) * (2 * dH * dphi * H * (3 * G5x + 7 * G5xx * X + 2 * G5xxx * pow(X, 2)) + ddphi * (3 * G3x - 9 * G4px + 3 * G5x * pow(H, 2) + 15 * G3xx * X - 36 * G4pxx * X + 27 * G5xx * pow(H, 2) * X + 6 * G3xxx * pow(X, 2) - 12 * G4pxxx * pow(X, 2) + 24 * G5xxx * pow(H, 2) * pow(X, 2) + 4 * G5xxxx * pow(H, 2) * pow(X, 3))) -
             pow(a, 2) * (dH * (G2x - 2 * (G3p + 3 * G4x * pow(H, 2) - 3 * G5p * pow(H, 2) - G2xx * X +
                                           G3px * X + 24 * G4xx * pow(H, 2) * X - 15 * G5px * pow(H, 2) * X +
                                           12 * G4xxx * pow(H, 2) * pow(X, 2) - 6 * G5pxx * pow(H, 2) * pow(X, 2))) +
                          dphi * H * (-G2px + 2 * (G3pp + 9 * G3x * pow(H, 2) - 30 * G4px * pow(H, 2) + 3 * G5pp * pow(H, 2) + 9 * G5x * pow(H, 4) - G2pxx * X + G3ppx * X + 21 * G3xx * pow(H, 2) * X - 72 * G4pxx * pow(H, 2) * X + 15 * G5ppx * pow(H, 2) * X + 41 * G5xx * pow(H, 4) * X + 6 * G3xxx * pow(H, 2) * pow(X, 2) - 24 * G4pxxx * pow(H, 2) * pow(X, 2) + 6 * G5ppxx * pow(H, 2) * pow(X, 2) + 28 * G5xxx * pow(H, 4) * pow(X, 2) + 4 * G5xxxx * pow(H, 4) * pow(X, 3))))) /
            (pow(a, 4) * pow(H, 2));

    *ch2p = 6 * (G4p + (ddphi * dphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2) + G3xx * X - 2 * G4pxx * X + 7 * G5xx * pow(H, 2) * X + 2 * G5xxx * pow(H, 2) * pow(X, 2))) / (pow(a, 3) * H) -
                 X * (3 * G3x - 12 * G4px + 4 * G5pp + 9 * G5x * pow(H, 2) + 2 * G3xx * X - 12 * G4pxx * X +
                      4 * G5ppx * X + 16 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2)) +
                 (2 * ddphi * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X))) / pow(a, 2) +
                 (-(dphi * H * (G4pp + 4 * G4x * pow(H, 2) - 4 * G5p * pow(H, 2) - G3px * X + 2 * G4ppx * X + 20 * G4xx * pow(H, 2) * X - 15 * G5px * pow(H, 2) * X + 8 * G4xxx * pow(H, 2) * pow(X, 2) - 6 * G5pxx * pow(H, 2) * pow(X, 2))) +
                  dH * (G4p + X * (-G3x + 2 * G4px + pow(H, 2) * (3 * G5x + 2 * G5xx * X)))) /
                     (a * pow(H, 2)));

    *cs2p = 2 * a * X * (G2px - 2 * (G3pp + 3 * G3x * pow(H, 2) - 12 * G4px * pow(H, 2) + 3 * G5pp * pow(H, 2) + 3 * G5x * pow(H, 4) + 2 * G2pxx * X - 5 * G3ppx * X + 6 * G4pppx * X - 9 * G3xx * pow(H, 2) * X + 54 * G4pxx * pow(H, 2) * X - 23 * G5ppx * pow(H, 2) * X - 13 * G5xx * pow(H, 4) * X - 6 * G3xxx * pow(H, 2) * pow(X, 2) + 36 * G4pxxx * pow(H, 2) * pow(X, 2) - 14 * G5ppxx * pow(H, 2) * pow(X, 2) - 20 * G5xxx * pow(H, 4) * pow(X, 2) - 4 * G5xxxx * pow(H, 4) * pow(X, 3))) +
            (ddphi * dphi * (G2px - 2 * (G3pp - 3 * G3x * pow(H, 2) + 6 * G4px * pow(H, 2) + 3 * G5pp * pow(H, 2) - 3 * G5x * pow(H, 4) - G2pxx * X + G3ppx * X + 9 * G3xx * pow(H, 2) * X - 48 * G4pxx * pow(H, 2) * X + 15 * G5ppx * pow(H, 2) * X + 13 * G5xx * pow(H, 4) * X + 6 * G3xxx * pow(H, 2) * pow(X, 2) - 24 * G4pxxx * pow(H, 2) * pow(X, 2) + 6 * G5ppxx * pow(H, 2) * pow(X, 2) + 20 * G5xxx * pow(H, 4) * pow(X, 2) + 4 * G5xxxx * pow(H, 4) * pow(X, 3)))) /
                (pow(a, 2) * H) +
            (2 * (ddphi * (G2x - 2 * G3p + 6 * G4x * pow(H, 2) - 6 * G5p * pow(H, 2) -
                           G2xx * X + 8 * G3px * X - 18 * G4ppx * X - 6 * G4xx * pow(H, 2) * X + 12 * G5px * pow(H, 2) * X -
                           2 * G2xxx * pow(X, 2) + 8 * G3pxx * pow(X, 2) - 12 * G4ppxx * pow(X, 2) -
                           72 * G4xxx * pow(H, 2) * pow(X, 2) + 56 * G5pxx * pow(H, 2) * pow(X, 2) -
                           24 * G4xxxx * pow(H, 2) * pow(X, 3) + 16 * G5pxxx * pow(H, 2) * pow(X, 3)) +
                  4 * dH * dphi * H * (3 * G4x - 3 * G5p + X * (-3 * G4xx + 4 * G5px - 6 * G4xxx * X + 4 * G5pxx * X)))) /
                a +
            (-(dphi * H * (G2pp + 2 * (G2x * pow(H, 2) - 2 * G3p * pow(H, 2) + 6 * G4pp * pow(H, 2) + 6 * G4x * pow(H, 4) - 6 * G5p * pow(H, 4) - G2ppx * X + G3ppp * X - G2xx * pow(H, 2) * X + 5 * G3px * pow(H, 2) * X - 18 * G4ppx * pow(H, 2) * X + 3 * G5ppp * pow(H, 2) * X - 6 * G4xx * pow(H, 4) * X + 9 * G5px * pow(H, 4) * X - 2 * G2xxx * pow(H, 2) * pow(X, 2) + 14 * G3pxx * pow(H, 2) * pow(X, 2) - 36 * G4ppxx * pow(H, 2) * pow(X, 2) + 6 * G5pppx * pow(H, 2) * pow(X, 2) - 72 * G4xxx * pow(H, 4) * pow(X, 2) + 64 * G5pxx * pow(H, 4) * pow(X, 2) - 24 * G4xxxx * pow(H, 4) * pow(X, 3) + 20 * G5pxxx * pow(H, 4) * pow(X, 3)))) +
             dH * (G2p - 2 * (6 * G4p * pow(H, 2) +
                              X * (G2px - G3pp + 3 * pow(H, 2) * (-G3x + G5pp - 3 * G5x * pow(H, 2) + 2 * G3xx * X - 8 * G4pxx * X + 2 * G5ppx * X + 8 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2)))))) /
                pow(H, 2);
}

/** function to compute coefficients of designer background equations */
void horndeski_designereq_coefficients(
    double a,
    double H,
    double dH,
    double ddH,
    double dphi,
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
    double *d0,
    double *d1,
    double *d2)
{
    double X = 0.5 * dphi * dphi / a / a;

    *d0 = (4 * dphi * G5x * (pow(dH, 2) + ddH * H) * X -
           4 * a * (ddH * (G4 + (-2 * G4x + G5p) * X) + 2 * dH * dphi * X * (G3x - 4 * G4px + G5pp + 6 * G5x * pow(H, 2) + 4 * G5xx * pow(H, 2) * X)) + pow(a, 3) * H * (-rptot + 2 * X * (G2x + 2 * (-G3p + 3 * G4pp + 5 * G4x * pow(H, 2) - 5 * G5p * pow(H, 2) + G2xx * X - 6 * G3px * X + 14 * G4ppx * X - 2 * G5ppp * X + 40 * G4xx * pow(H, 2) * X - 34 * G5px * pow(H, 2) * X + 20 * G4xxx * pow(H, 2) * pow(X, 2) - 16 * G5pxx * pow(H, 2) * pow(X, 2)))) +
           pow(a, 2) * (-drptot + 2 * X * (-(dphi * G2px) - 8 * dH * H * (3 * G4x - 3 * G5p + 6 * G4xx * X - 4 * G5px * X) + 2 * dphi * (G3pp - G4pxx * (1 + 22 * pow(H, 2) * X) + pow(H, 2) * (4 * G3x - 19 * G4px + 7 * G5pp + 6 * G5x * pow(H, 2) + 4 * G3xx * X + 8 * G5ppx * X + 14 * G5xx * pow(H, 2) * X + 4 * G5xxx * pow(H, 2) * pow(X, 2)))))) /
          (3. * a);

    *d1 = (-2 * (4 * dH * dphi * (-G4x + G5p - 2 * G4xx * X + G5px * X) +
                 a * (-4 * dH * H * X * (3 * G5x + 2 * G5xx * X) +
                      dphi * (G2x - 2 * G3p + 3 * G4pp + 14 * G4x * pow(H, 2) - 14 * G5p * pow(H, 2) +
                              G2xx * X - 3 * G3px * X + 4 * G4ppx * X + 70 * G4xx * pow(H, 2) * X -
                              49 * G5px * pow(H, 2) * X + 28 * G4xxx * pow(H, 2) * pow(X, 2) -
                              18 * G5pxx * pow(H, 2) * pow(X, 2))) +
                 pow(a, 2) * H * (-3 * G4p + X * (15 * G3x - 52 * G4px + 10 * G5pp + 27 * G5x * pow(H, 2) + 10 * G3xx * X - 36 * G4pxx * X + 8 * G5ppx * X + 48 * G5xx * pow(H, 2) * X + 12 * G5xxx * pow(H, 2) * pow(X, 2))))) /
          (3. * pow(a, 2));

    *d2 = (2 * dphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2) + G3xx * X - 2 * G4pxx * X + 7 * G5xx * pow(H, 2) * X + 2 * G5xxx * pow(H, 2) * pow(X, 2)) +
           4 * a * H * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X))) /
          (3. * pow(a, 3));
}

/** function to compute coefficients of designer background equations */
void horndeski_designereq_coefficients_p(
    double a,
    double H,
    double dH,
    double ddH,
    double dddH,
    double dphi,
    double ddphi,
    double rptot,
    double drptot,
    double ddrptot,
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
    double *dp0,
    double *dp1)
{
    double X = 0.5 * dphi * dphi / a / a;

    *dp0 = (-2 * pow(dphi, 9) * G5xxxx * pow(H, 5) +
            4 * pow(a, 2) * pow(ddphi, 3) * dphi * H *
                (9 * G4xx - 6 * G5px + X * (16 * G4xxx - 9 * G5pxx + 4 * G4xxxx * X - 2 * G5pxxx * X)) +
            4 * pow(a, 4) * ddphi * (3 * ddH * dphi * (G4x - G5p + 2 * G4xx * X - G5px * X) - ddphi * dphi * H * (9 * G3x - 33 * G4px + 6 * G5pp + 18 * G5x * pow(H, 2) + 21 * G3xx * X - 78 * G4pxx * X + 18 * G5ppx * X + 82 * G5xx * pow(H, 2) * X + 6 * G3xxx * pow(X, 2) - 24 * G4pxxx * pow(X, 2) + 6 * G5ppxx * pow(X, 2) + 56 * G5xxx * pow(H, 2) * pow(X, 2) + 8 * G5xxxx * pow(H, 2) * pow(X, 3)) + 3 * ddphi * dH * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X))) -
            2 * pow(a, 9) * (2 * (G2ppx - 2 * G3ppp + 2 * G4pppp) * pow(X, 2) + pow(H, 2) * (rptot + 2 * pow(X, 2) * (3 * G2xx - 22 * G3px + 72 * G4ppx - 18 * G5ppp + 2 * G2xxx * X - 20 * G3pxx * X + 72 * G4ppxx * X - 20 * G5pppx * X)) + 20 * pow(H, 4) * pow(X, 2) * (9 * G4xx - 9 * G5px + 4 * X * (4 * G4xxx - 4 * G5pxx + G4xxxx * X - G5pxxx * X))) +
            2 * pow(a, 3) * pow(ddphi, 2) * (ddphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2) + 5 * G3xx * X - 12 * G4pxx * X + 27 * G5xx * pow(H, 2) * X + 2 * G3xxx * pow(X, 2) - 4 * G4pxxx * pow(X, 2) + 24 * G5xxx * pow(H, 2) * pow(X, 2) + 4 * G5xxxx * pow(H, 2) * pow(X, 3)) + 6 * dH * dphi * H * (3 * G5x + X * (7 * G5xx + 2 * G5xxx * X))) -
            pow(a, 8) * (2 * drptot * H + dH * rptot - 4 * dphi * G2px * H * X +
                         8 * dphi * H * X * (G3pp - 2 * G4pxx + 2 * G3x * pow(H, 2) - 12 * G4px * pow(H, 2) + 6 * G5pp * pow(H, 2) + 3 * G5x * pow(H, 4) - G2pxx * X + 4 * G3ppx * X - 8 * G4pppx * X + G5pppp * X + 10 * G3xx * pow(H, 2) * X - 72 * G4pxx * pow(H, 2) * X + 36 * G5ppx * pow(H, 2) * X + 27 * G5xx * pow(H, 4) * X + 4 * G3xxx * pow(H, 2) * pow(X, 2) - 32 * G4pxxx * pow(H, 2) * pow(X, 2) + 16 * G5ppxx * pow(H, 2) * pow(X, 2) + 24 * G5xxx * pow(H, 4) * pow(X, 2)) -
                         2 * dH * X * (G2x + 2 * (-G3p + 3 * G4pp + 27 * G4x * pow(H, 2) - 27 * G5p * pow(H, 2) + G2xx * X - 10 * G3px * X + 30 * G4ppx * X - 6 * G5ppp * X + 216 * G4xx * pow(H, 2) * X - 198 * G5px * pow(H, 2) * X + 108 * G4xxx * pow(H, 2) * pow(X, 2) - 96 * G5pxx * pow(H, 2) * pow(X, 2)))) -
            2 * pow(a, 5) * (-2 * dphi * G5x * (3 * ddH * dH + dddH * H) * X + pow(ddphi, 2) * (G2x - 2 * G3p + 3 * G4pp + 18 * G4x * pow(H, 2) - 18 * G5p * pow(H, 2) + 5 * G2xx * X - 15 * G3px * X + 24 * G4ppx * X + 306 * G4xx * pow(H, 2) * X - 225 * G5px * pow(H, 2) * X + 2 * G2xxx * pow(X, 2) - 8 * G3pxx * pow(X, 2) + 12 * G4ppxx * pow(X, 2) + 360 * G4xxx * pow(H, 2) * pow(X, 2) - 246 * G5pxx * pow(H, 2) * pow(X, 2) + 72 * G4xxxx * pow(H, 2) * pow(X, 3) - 48 * G5pxxx * pow(H, 2) * pow(X, 3)) - 6 * ddphi * (pow(dH, 2) * X * (3 * G5x + 2 * G5xx * X) + ddH * H * X * (3 * G5x + 2 * G5xx * X) + dH * dphi * H * (-10 * G4x + 10 * G5p + X * (-50 * G4xx + 37 * G5px - 20 * G4xxx * X + 14 * G5pxx * X)))) +
            pow(a, 7) * (-ddrptot + 2 * X *
                                        (-5 * ddphi * G2px - 4 * ddH * H * (7 * G4x - 7 * G5p + 2 * (7 * G4xx - 5 * G5px) * X) -
                                         12 * pow(dH, 2) * (2 * G4x - 2 * G5p + 4 * G4xx * X - 3 * G5px * X) +
                                         4 * dH * dphi * H * (6 * G3x - 33 * G4px + 15 * G5pp + 24 * G5x * pow(H, 2) + 6 * G3xx * X - 42 * G4pxx * X + 18 * G5ppx * X + 56 * G5xx * pow(H, 2) * X + 16 * G5xxx * pow(H, 2) * pow(X, 2)) +
                                         2 * ddphi * (5 * G3pp + 27 * G3x * pow(H, 2) - 126 * G4px * pow(H, 2) + 45 * G5pp * pow(H, 2) + 45 * G5x * pow(H, 4) - 2 * G2pxx * X + 5 * G3ppx * X - 6 * G4pppx * X + 63 * G3xx * pow(H, 2) * X + 129 * G5ppx * pow(H, 2) * X + 205 * G5xx * pow(H, 4) * X + 18 * G3xxx * pow(H, 2) * pow(X, 2) - 108 * G4pxxx * pow(H, 2) * pow(X, 2) + 42 * G5ppxx * pow(H, 2) * pow(X, 2) + 140 * G5xxx * pow(H, 4) * pow(X, 2) + 20 * G5xxxx * pow(H, 4) * pow(X, 3) - 6 * G4pxx * (1 + 57 * pow(H, 2) * X)))) -
            2 * pow(a, 6) * (2 * dddH * (G4 + (-2 * G4x + G5p) * X) + ddphi * (-3 * dH * G4p + 3 * dH * X * (9 * G3x - 36 * G4px + 10 * G5pp + 63 * G5x * pow(H, 2) + 6 * G3xx * X - 28 * G4pxx * X + 8 * G5ppx * X + 112 * G5xx * pow(H, 2) * X + 28 * G5xxx * pow(H, 2) * pow(X, 2)) - 2 * dphi * H * (G2x - 2 * G3p + 6 * G4pp + 12 * G4x * pow(H, 2) - 12 * G5p * pow(H, 2) + 5 * G2xx * X - 27 * G3px * X + 66 * G4ppx * X - 9 * G5ppp * X + 204 * G4xx * pow(H, 2) * X - 174 * G5px * pow(H, 2) * X + 2 * G2xxx * pow(X, 2) - 14 * G3pxx * pow(X, 2) + 36 * G4ppxx * pow(X, 2) - 6 * G5pppx * pow(X, 2) + 240 * G4xxx * pow(H, 2) * pow(X, 2) - 200 * G5pxx * pow(H, 2) * pow(X, 2) + 48 * G4xxxx * pow(H, 2) * pow(X, 3) - 40 * G5pxxx * pow(H, 2) * pow(X, 3))) + 2 * dphi * (9 * pow(dH, 2) * H * X * (3 * G5x + 2 * G5xx * X) + ddH * (G4p + X * (2 * G3x - 10 * G4px + 3 * G5pp + 15 * G5x * pow(H, 2) + 10 * G5xx * pow(H, 2) * X))))) /
           (3. * pow(a, 6));

    *dp1 = (6 * ddphi * dphi * (G3x - 3 * G4px + 3 * G5x * pow(H, 2) + G3xx * X - 2 * G4pxx * X + 7 * G5xx * pow(H, 2) * X + 2 * G5xxx * pow(H, 2) * pow(X, 2)) -
            2 * pow(a, 2) * (-6 * dH * H * X * (3 * G5x + 2 * G5xx * X) + dphi * (G2x - 2 * G3p + 4 * G4pp + 18 * G4x * pow(H, 2) - 18 * G5p * pow(H, 2) + G2xx * X - 4 * G3px * X + 6 * G4ppx * X + 90 * G4xx * pow(H, 2) * X - 64 * G5px * pow(H, 2) * X + 36 * G4xxx * pow(H, 2) * pow(X, 2) - 24 * G5pxx * pow(H, 2) * pow(X, 2))) -
            4 * pow(a, 3) * H * (-2 * G4p + X * (9 * G3x - 32 * G4px + 7 * G5pp + 18 * G5x * pow(H, 2) + 6 * G3xx * X - 24 * G4pxx * X + 6 * G5ppx * X + 32 * G5xx * pow(H, 2) * X + 8 * G5xxx * pow(H, 2) * pow(X, 2))) +
            12 * a * (dH * dphi * (G4x - G5p + 2 * G4xx * X - G5px * X) + ddphi * H * (G4x - G5p + X * (8 * G4xx - 5 * G5px + 4 * G4xxx * X - 2 * G5pxx * X)))) /
           (3. * pow(a, 3));
}
