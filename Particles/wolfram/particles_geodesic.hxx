/* particles_geodesic.hxx */
/* Produced with Generato */

#ifndef PARTICLES_GEODESIC_HXX
#define PARTICLES_GEODESIC_HXX

#include "particles_powerinline.hxx"
#include "../src/Particles.hxx"

namespace Particles {

CCTK_HOST CCTK_DEVICE inline void
calc_rhs_geodesic(VectR &dtxpos, VectR &dtpmom, const VectR &pmom,
                  const ScalR alpha, const VectR &beta, const SmatR &gam,
                  const dScalR &dalpha, const dVectR &dbeta, const dSmatR &dgam) {

auto &dtxpos1 = dtxpos[0];
auto &dtxpos2 = dtxpos[1];
auto &dtxpos3 = dtxpos[2];
auto &dtpmom1 = dtpmom[0];
auto &dtpmom2 = dtpmom[1];
auto &dtpmom3 = dtpmom[2];

const auto &pmom1 = pmom[0];
const auto &pmom2 = pmom[1];
const auto &pmom3 = pmom[2];

const auto &beta1 = beta[0];
const auto &beta2 = beta[1];
const auto &beta3 = beta[2];
const auto &gam11 = gam[0];
const auto &gam12 = gam[1];
const auto &gam13 = gam[2];
const auto &gam22 = gam[3];
const auto &gam23 = gam[4];
const auto &gam33 = gam[5];

const auto &dalpha1 = dalpha[0];
const auto &dalpha2 = dalpha[1];
const auto &dalpha3 = dalpha[2];
const auto &dbeta11 = dbeta[0][0];
const auto &dbeta12 = dbeta[0][1];
const auto &dbeta13 = dbeta[0][2];
const auto &dbeta21 = dbeta[1][0];
const auto &dbeta22 = dbeta[1][1];
const auto &dbeta23 = dbeta[1][2];
const auto &dbeta31 = dbeta[2][0];
const auto &dbeta32 = dbeta[2][1];
const auto &dbeta33 = dbeta[2][2];
const auto &dgam111 = dgam[0][0];
const auto &dgam112 = dgam[0][1];
const auto &dgam113 = dgam[0][2];
const auto &dgam121 = dgam[1][0];
const auto &dgam122 = dgam[1][1];
const auto &dgam123 = dgam[1][2];
const auto &dgam131 = dgam[2][0];
const auto &dgam132 = dgam[2][1];
const auto &dgam133 = dgam[2][2];
const auto &dgam221 = dgam[3][0];
const auto &dgam222 = dgam[3][1];
const auto &dgam223 = dgam[3][2];
const auto &dgam231 = dgam[4][0];
const auto &dgam232 = dgam[4][1];
const auto &dgam233 = dgam[4][2];
const auto &dgam331 = dgam[5][0];
const auto &dgam332 = dgam[5][1];
const auto &dgam333 = dgam[5][2];

const auto
detinvgam
=
1/(-(Power(gam13,2)*gam22) + 2*gam12*gam13*gam23 - gam11*Power(gam23,2) -
    Power(gam12,2)*gam33 + gam11*gam22*gam33)
;

const auto
invgam11
=
detinvgam*(-Power(gam23,2) + gam22*gam33)
;

const auto
invgam12
=
detinvgam*(gam13*gam23 - gam12*gam33)
;

const auto
invgam13
=
detinvgam*(-(gam13*gam22) + gam12*gam23)
;

const auto
invgam22
=
detinvgam*(-Power(gam13,2) + gam11*gam33)
;

const auto
invgam23
=
detinvgam*(gam12*gam13 - gam11*gam23)
;

const auto
invgam33
=
detinvgam*(-Power(gam12,2) + gam11*gam22)
;

const auto
dinvgam111
=
-(dgam111*Power(invgam11,2)) - 2*dgam121*invgam11*invgam12 -
  dgam221*Power(invgam12,2) - 2*dgam131*invgam11*invgam13 -
  2*dgam231*invgam12*invgam13 - dgam331*Power(invgam13,2)
;

const auto
dinvgam112
=
-(dgam112*Power(invgam11,2)) - 2*dgam122*invgam11*invgam12 -
  dgam222*Power(invgam12,2) - 2*dgam132*invgam11*invgam13 -
  2*dgam232*invgam12*invgam13 - dgam332*Power(invgam13,2)
;

const auto
dinvgam113
=
-(dgam113*Power(invgam11,2)) - 2*dgam123*invgam11*invgam12 -
  dgam223*Power(invgam12,2) - 2*dgam133*invgam11*invgam13 -
  2*dgam233*invgam12*invgam13 - dgam333*Power(invgam13,2)
;

const auto
dinvgam121
=
-(dgam111*invgam11*invgam12) - dgam131*invgam12*invgam13 -
  dgam221*invgam12*invgam22 - dgam231*invgam13*invgam22 -
  dgam121*(Power(invgam12,2) + invgam11*invgam22) -
  dgam131*invgam11*invgam23 - dgam231*invgam12*invgam23 -
  dgam331*invgam13*invgam23
;

const auto
dinvgam122
=
-(dgam112*invgam11*invgam12) - dgam132*invgam12*invgam13 -
  dgam222*invgam12*invgam22 - dgam232*invgam13*invgam22 -
  dgam122*(Power(invgam12,2) + invgam11*invgam22) -
  dgam132*invgam11*invgam23 - dgam232*invgam12*invgam23 -
  dgam332*invgam13*invgam23
;

const auto
dinvgam123
=
-(dgam113*invgam11*invgam12) - dgam133*invgam12*invgam13 -
  dgam223*invgam12*invgam22 - dgam233*invgam13*invgam22 -
  dgam123*(Power(invgam12,2) + invgam11*invgam22) -
  dgam133*invgam11*invgam23 - dgam233*invgam12*invgam23 -
  dgam333*invgam13*invgam23
;

const auto
dinvgam131
=
-(dgam111*invgam11*invgam13) - dgam131*Power(invgam13,2) -
  dgam221*invgam12*invgam23 - dgam231*invgam13*invgam23 -
  dgam121*(invgam12*invgam13 + invgam11*invgam23) -
  dgam131*invgam11*invgam33 - dgam231*invgam12*invgam33 -
  dgam331*invgam13*invgam33
;

const auto
dinvgam132
=
-(dgam112*invgam11*invgam13) - dgam132*Power(invgam13,2) -
  dgam222*invgam12*invgam23 - dgam232*invgam13*invgam23 -
  dgam122*(invgam12*invgam13 + invgam11*invgam23) -
  dgam132*invgam11*invgam33 - dgam232*invgam12*invgam33 -
  dgam332*invgam13*invgam33
;

const auto
dinvgam133
=
-(dgam113*invgam11*invgam13) - dgam133*Power(invgam13,2) -
  dgam223*invgam12*invgam23 - dgam233*invgam13*invgam23 -
  dgam123*(invgam12*invgam13 + invgam11*invgam23) -
  dgam133*invgam11*invgam33 - dgam233*invgam12*invgam33 -
  dgam333*invgam13*invgam33
;

const auto
dinvgam221
=
-(dgam111*Power(invgam12,2)) - 2*dgam121*invgam12*invgam22 -
  dgam221*Power(invgam22,2) - 2*dgam131*invgam12*invgam23 -
  2*dgam231*invgam22*invgam23 - dgam331*Power(invgam23,2)
;

const auto
dinvgam222
=
-(dgam112*Power(invgam12,2)) - 2*dgam122*invgam12*invgam22 -
  dgam222*Power(invgam22,2) - 2*dgam132*invgam12*invgam23 -
  2*dgam232*invgam22*invgam23 - dgam332*Power(invgam23,2)
;

const auto
dinvgam223
=
-(dgam113*Power(invgam12,2)) - 2*dgam123*invgam12*invgam22 -
  dgam223*Power(invgam22,2) - 2*dgam133*invgam12*invgam23 -
  2*dgam233*invgam22*invgam23 - dgam333*Power(invgam23,2)
;

const auto
dinvgam231
=
-(dgam111*invgam12*invgam13) - dgam131*invgam13*invgam23 -
  dgam221*invgam22*invgam23 - dgam231*Power(invgam23,2) -
  dgam121*(invgam13*invgam22 + invgam12*invgam23) -
  dgam131*invgam12*invgam33 - dgam231*invgam22*invgam33 -
  dgam331*invgam23*invgam33
;

const auto
dinvgam232
=
-(dgam112*invgam12*invgam13) - dgam132*invgam13*invgam23 -
  dgam222*invgam22*invgam23 - dgam232*Power(invgam23,2) -
  dgam122*(invgam13*invgam22 + invgam12*invgam23) -
  dgam132*invgam12*invgam33 - dgam232*invgam22*invgam33 -
  dgam332*invgam23*invgam33
;

const auto
dinvgam233
=
-(dgam113*invgam12*invgam13) - dgam133*invgam13*invgam23 -
  dgam223*invgam22*invgam23 - dgam233*Power(invgam23,2) -
  dgam123*(invgam13*invgam22 + invgam12*invgam23) -
  dgam133*invgam12*invgam33 - dgam233*invgam22*invgam33 -
  dgam333*invgam23*invgam33
;

const auto
dinvgam331
=
-(dgam111*Power(invgam13,2)) - 2*dgam121*invgam13*invgam23 -
  dgam221*Power(invgam23,2) - 2*dgam131*invgam13*invgam33 -
  2*dgam231*invgam23*invgam33 - dgam331*Power(invgam33,2)
;

const auto
dinvgam332
=
-(dgam112*Power(invgam13,2)) - 2*dgam122*invgam13*invgam23 -
  dgam222*Power(invgam23,2) - 2*dgam132*invgam13*invgam33 -
  2*dgam232*invgam23*invgam33 - dgam332*Power(invgam33,2)
;

const auto
dinvgam333
=
-(dgam113*Power(invgam13,2)) - 2*dgam123*invgam13*invgam23 -
  dgam223*Power(invgam23,2) - 2*dgam133*invgam13*invgam33 -
  2*dgam233*invgam23*invgam33 - dgam333*Power(invgam33,2)
;

const auto
p2
=
invgam11*Power(pmom1,2) + 2*invgam12*pmom1*pmom2 + invgam22*Power(pmom2,2) +
  2*invgam13*pmom1*pmom3 + 2*invgam23*pmom2*pmom3 + invgam33*Power(pmom3,2)
;

const auto
pmomt
=
sqrt(p2)/alpha
;


dtxpos1
=
(invgam11*pmom1 + invgam12*pmom2 + invgam13*pmom3 - beta1*pmomt)/pmomt
;

dtxpos2
=
(invgam12*pmom1 + invgam22*pmom2 + invgam23*pmom3 - beta2*pmomt)/pmomt
;

dtxpos3
=
(invgam13*pmom1 + invgam23*pmom2 + invgam33*pmom3 - beta3*pmomt)/pmomt
;

dtpmom1
=
dbeta11*pmom1 - (dinvgam111*Power(pmom1,2))/2. + dbeta21*pmom2 -
  dinvgam121*pmom1*pmom2 - (dinvgam221*Power(pmom2,2))/2. + dbeta31*pmom3 -
  dinvgam131*pmom1*pmom3 - dinvgam231*pmom2*pmom3 -
  (dinvgam331*Power(pmom3,2))/2. - alpha*dalpha1*pmomt
;

dtpmom2
=
dbeta12*pmom1 - (dinvgam112*Power(pmom1,2))/2. + dbeta22*pmom2 -
  dinvgam122*pmom1*pmom2 - (dinvgam222*Power(pmom2,2))/2. + dbeta32*pmom3 -
  dinvgam132*pmom1*pmom3 - dinvgam232*pmom2*pmom3 -
  (dinvgam332*Power(pmom3,2))/2. - alpha*dalpha2*pmomt
;

dtpmom3
=
dbeta13*pmom1 - (dinvgam113*Power(pmom1,2))/2. + dbeta23*pmom2 -
  dinvgam123*pmom1*pmom2 - (dinvgam223*Power(pmom2,2))/2. + dbeta33*pmom3 -
  dinvgam133*pmom1*pmom3 - dinvgam233*pmom2*pmom3 -
  (dinvgam333*Power(pmom3,2))/2. - alpha*dalpha3*pmomt
;


}

} // namespace Particles

#endif // #ifndef PARTICLES_GEODESIC_HXX

/* particles_geodesic.hxx */
