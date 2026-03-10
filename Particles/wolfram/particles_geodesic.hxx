/* particles_geodesic.hxx */
/* Produced with Generato */

#ifndef PARTICLES_GEODESIC_HXX
#define PARTICLES_GEODESIC_HXX

#include "util_powerinline.hxx"
#include "../src/Particles.hxx"

namespace Particles {
using namespace UtilForge;

CCTK_HOST CCTK_DEVICE inline void
calc_rhs_geodesic(VectR &dtpmom, VectR &dtxpos, const VectR &pmom,
                  const ScalR ADMalpha, const VectR &ADMbeta, const SmatR &ADMgam,
                  const dScalR &ADMdalpha, const dVectR &ADMdbeta, const dSmatR &ADMdgam) {

auto &dtxpos1 = dtxpos[0];
auto &dtxpos2 = dtxpos[1];
auto &dtxpos3 = dtxpos[2];
auto &dtpmom1 = dtpmom[0];
auto &dtpmom2 = dtpmom[1];
auto &dtpmom3 = dtpmom[2];

const auto &pmom1 = pmom[0];
const auto &pmom2 = pmom[1];
const auto &pmom3 = pmom[2];

const auto &ADMgam11 = ADMgam[0];
const auto &ADMgam12 = ADMgam[1];
const auto &ADMgam13 = ADMgam[2];
const auto &ADMgam22 = ADMgam[3];
const auto &ADMgam23 = ADMgam[4];
const auto &ADMgam33 = ADMgam[5];
const auto &ADMbeta1 = ADMbeta[0];
const auto &ADMbeta2 = ADMbeta[1];
const auto &ADMbeta3 = ADMbeta[2];

const auto &ADMdgam111 = ADMdgam[0][0];
const auto &ADMdgam112 = ADMdgam[0][1];
const auto &ADMdgam113 = ADMdgam[0][2];
const auto &ADMdgam121 = ADMdgam[1][0];
const auto &ADMdgam122 = ADMdgam[1][1];
const auto &ADMdgam123 = ADMdgam[1][2];
const auto &ADMdgam131 = ADMdgam[2][0];
const auto &ADMdgam132 = ADMdgam[2][1];
const auto &ADMdgam133 = ADMdgam[2][2];
const auto &ADMdgam221 = ADMdgam[3][0];
const auto &ADMdgam222 = ADMdgam[3][1];
const auto &ADMdgam223 = ADMdgam[3][2];
const auto &ADMdgam231 = ADMdgam[4][0];
const auto &ADMdgam232 = ADMdgam[4][1];
const auto &ADMdgam233 = ADMdgam[4][2];
const auto &ADMdgam331 = ADMdgam[5][0];
const auto &ADMdgam332 = ADMdgam[5][1];
const auto &ADMdgam333 = ADMdgam[5][2];
const auto &ADMdalpha1 = ADMdalpha[0];
const auto &ADMdalpha2 = ADMdalpha[1];
const auto &ADMdalpha3 = ADMdalpha[2];
const auto &ADMdbeta11 = ADMdbeta[0][0];
const auto &ADMdbeta12 = ADMdbeta[0][1];
const auto &ADMdbeta13 = ADMdbeta[0][2];
const auto &ADMdbeta21 = ADMdbeta[1][0];
const auto &ADMdbeta22 = ADMdbeta[1][1];
const auto &ADMdbeta23 = ADMdbeta[1][2];
const auto &ADMdbeta31 = ADMdbeta[2][0];
const auto &ADMdbeta32 = ADMdbeta[2][1];
const auto &ADMdbeta33 = ADMdbeta[2][2];

const auto
detinvgam
=
1/(-(Power(ADMgam13,2)*ADMgam22) + 2*ADMgam12*ADMgam13*ADMgam23 -
    ADMgam11*Power(ADMgam23,2) - Power(ADMgam12,2)*ADMgam33 +
    ADMgam11*ADMgam22*ADMgam33)
;

const auto
invgam11
=
(-Power(ADMgam23,2) + ADMgam22*ADMgam33)*detinvgam
;

const auto
invgam12
=
(ADMgam13*ADMgam23 - ADMgam12*ADMgam33)*detinvgam
;

const auto
invgam13
=
(-(ADMgam13*ADMgam22) + ADMgam12*ADMgam23)*detinvgam
;

const auto
invgam22
=
(-Power(ADMgam13,2) + ADMgam11*ADMgam33)*detinvgam
;

const auto
invgam23
=
(ADMgam12*ADMgam13 - ADMgam11*ADMgam23)*detinvgam
;

const auto
invgam33
=
(-Power(ADMgam12,2) + ADMgam11*ADMgam22)*detinvgam
;

const auto
dinvgam111
=
-(ADMdgam111*Power(invgam11,2)) - 2*ADMdgam121*invgam11*invgam12 -
  ADMdgam221*Power(invgam12,2) - 2*ADMdgam131*invgam11*invgam13 -
  2*ADMdgam231*invgam12*invgam13 - ADMdgam331*Power(invgam13,2)
;

const auto
dinvgam112
=
-(ADMdgam112*Power(invgam11,2)) - 2*ADMdgam122*invgam11*invgam12 -
  ADMdgam222*Power(invgam12,2) - 2*ADMdgam132*invgam11*invgam13 -
  2*ADMdgam232*invgam12*invgam13 - ADMdgam332*Power(invgam13,2)
;

const auto
dinvgam113
=
-(ADMdgam113*Power(invgam11,2)) - 2*ADMdgam123*invgam11*invgam12 -
  ADMdgam223*Power(invgam12,2) - 2*ADMdgam133*invgam11*invgam13 -
  2*ADMdgam233*invgam12*invgam13 - ADMdgam333*Power(invgam13,2)
;

const auto
dinvgam121
=
-(ADMdgam111*invgam11*invgam12) - ADMdgam131*invgam12*invgam13 -
  ADMdgam221*invgam12*invgam22 - ADMdgam231*invgam13*invgam22 -
  ADMdgam121*(Power(invgam12,2) + invgam11*invgam22) -
  ADMdgam131*invgam11*invgam23 - ADMdgam231*invgam12*invgam23 -
  ADMdgam331*invgam13*invgam23
;

const auto
dinvgam122
=
-(ADMdgam112*invgam11*invgam12) - ADMdgam132*invgam12*invgam13 -
  ADMdgam222*invgam12*invgam22 - ADMdgam232*invgam13*invgam22 -
  ADMdgam122*(Power(invgam12,2) + invgam11*invgam22) -
  ADMdgam132*invgam11*invgam23 - ADMdgam232*invgam12*invgam23 -
  ADMdgam332*invgam13*invgam23
;

const auto
dinvgam123
=
-(ADMdgam113*invgam11*invgam12) - ADMdgam133*invgam12*invgam13 -
  ADMdgam223*invgam12*invgam22 - ADMdgam233*invgam13*invgam22 -
  ADMdgam123*(Power(invgam12,2) + invgam11*invgam22) -
  ADMdgam133*invgam11*invgam23 - ADMdgam233*invgam12*invgam23 -
  ADMdgam333*invgam13*invgam23
;

const auto
dinvgam131
=
-(ADMdgam111*invgam11*invgam13) - ADMdgam131*Power(invgam13,2) -
  ADMdgam221*invgam12*invgam23 - ADMdgam231*invgam13*invgam23 -
  ADMdgam121*(invgam12*invgam13 + invgam11*invgam23) -
  ADMdgam131*invgam11*invgam33 - ADMdgam231*invgam12*invgam33 -
  ADMdgam331*invgam13*invgam33
;

const auto
dinvgam132
=
-(ADMdgam112*invgam11*invgam13) - ADMdgam132*Power(invgam13,2) -
  ADMdgam222*invgam12*invgam23 - ADMdgam232*invgam13*invgam23 -
  ADMdgam122*(invgam12*invgam13 + invgam11*invgam23) -
  ADMdgam132*invgam11*invgam33 - ADMdgam232*invgam12*invgam33 -
  ADMdgam332*invgam13*invgam33
;

const auto
dinvgam133
=
-(ADMdgam113*invgam11*invgam13) - ADMdgam133*Power(invgam13,2) -
  ADMdgam223*invgam12*invgam23 - ADMdgam233*invgam13*invgam23 -
  ADMdgam123*(invgam12*invgam13 + invgam11*invgam23) -
  ADMdgam133*invgam11*invgam33 - ADMdgam233*invgam12*invgam33 -
  ADMdgam333*invgam13*invgam33
;

const auto
dinvgam221
=
-(ADMdgam111*Power(invgam12,2)) - 2*ADMdgam121*invgam12*invgam22 -
  ADMdgam221*Power(invgam22,2) - 2*ADMdgam131*invgam12*invgam23 -
  2*ADMdgam231*invgam22*invgam23 - ADMdgam331*Power(invgam23,2)
;

const auto
dinvgam222
=
-(ADMdgam112*Power(invgam12,2)) - 2*ADMdgam122*invgam12*invgam22 -
  ADMdgam222*Power(invgam22,2) - 2*ADMdgam132*invgam12*invgam23 -
  2*ADMdgam232*invgam22*invgam23 - ADMdgam332*Power(invgam23,2)
;

const auto
dinvgam223
=
-(ADMdgam113*Power(invgam12,2)) - 2*ADMdgam123*invgam12*invgam22 -
  ADMdgam223*Power(invgam22,2) - 2*ADMdgam133*invgam12*invgam23 -
  2*ADMdgam233*invgam22*invgam23 - ADMdgam333*Power(invgam23,2)
;

const auto
dinvgam231
=
-(ADMdgam111*invgam12*invgam13) - ADMdgam131*invgam13*invgam23 -
  ADMdgam221*invgam22*invgam23 - ADMdgam231*Power(invgam23,2) -
  ADMdgam121*(invgam13*invgam22 + invgam12*invgam23) -
  ADMdgam131*invgam12*invgam33 - ADMdgam231*invgam22*invgam33 -
  ADMdgam331*invgam23*invgam33
;

const auto
dinvgam232
=
-(ADMdgam112*invgam12*invgam13) - ADMdgam132*invgam13*invgam23 -
  ADMdgam222*invgam22*invgam23 - ADMdgam232*Power(invgam23,2) -
  ADMdgam122*(invgam13*invgam22 + invgam12*invgam23) -
  ADMdgam132*invgam12*invgam33 - ADMdgam232*invgam22*invgam33 -
  ADMdgam332*invgam23*invgam33
;

const auto
dinvgam233
=
-(ADMdgam113*invgam12*invgam13) - ADMdgam133*invgam13*invgam23 -
  ADMdgam223*invgam22*invgam23 - ADMdgam233*Power(invgam23,2) -
  ADMdgam123*(invgam13*invgam22 + invgam12*invgam23) -
  ADMdgam133*invgam12*invgam33 - ADMdgam233*invgam22*invgam33 -
  ADMdgam333*invgam23*invgam33
;

const auto
dinvgam331
=
-(ADMdgam111*Power(invgam13,2)) - 2*ADMdgam121*invgam13*invgam23 -
  ADMdgam221*Power(invgam23,2) - 2*ADMdgam131*invgam13*invgam33 -
  2*ADMdgam231*invgam23*invgam33 - ADMdgam331*Power(invgam33,2)
;

const auto
dinvgam332
=
-(ADMdgam112*Power(invgam13,2)) - 2*ADMdgam122*invgam13*invgam23 -
  ADMdgam222*Power(invgam23,2) - 2*ADMdgam132*invgam13*invgam33 -
  2*ADMdgam232*invgam23*invgam33 - ADMdgam332*Power(invgam33,2)
;

const auto
dinvgam333
=
-(ADMdgam113*Power(invgam13,2)) - 2*ADMdgam123*invgam13*invgam23 -
  ADMdgam223*Power(invgam23,2) - 2*ADMdgam133*invgam13*invgam33 -
  2*ADMdgam233*invgam23*invgam33 - ADMdgam333*Power(invgam33,2)
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
sqrt(p2)/ADMalpha
;


dtxpos1
=
(invgam11*pmom1 + invgam12*pmom2 + invgam13*pmom3 - ADMbeta1*pmomt)/pmomt
;

dtxpos2
=
(invgam12*pmom1 + invgam22*pmom2 + invgam23*pmom3 - ADMbeta2*pmomt)/pmomt
;

dtxpos3
=
(invgam13*pmom1 + invgam23*pmom2 + invgam33*pmom3 - ADMbeta3*pmomt)/pmomt
;

dtpmom1
=
-0.5*(dinvgam111*Power(pmom1,2) + 2*dinvgam121*pmom1*pmom2 +
     dinvgam221*Power(pmom2,2) + 2*dinvgam131*pmom1*pmom3 +
     2*dinvgam231*pmom2*pmom3 + dinvgam331*Power(pmom3,2) -
     2*ADMdbeta11*pmom1*pmomt - 2*ADMdbeta21*pmom2*pmomt -
     2*ADMdbeta31*pmom3*pmomt + 2*ADMalpha*ADMdalpha1*Power(pmomt,2))/pmomt
;

dtpmom2
=
-0.5*(dinvgam112*Power(pmom1,2) + 2*dinvgam122*pmom1*pmom2 +
     dinvgam222*Power(pmom2,2) + 2*dinvgam132*pmom1*pmom3 +
     2*dinvgam232*pmom2*pmom3 + dinvgam332*Power(pmom3,2) -
     2*ADMdbeta12*pmom1*pmomt - 2*ADMdbeta22*pmom2*pmomt -
     2*ADMdbeta32*pmom3*pmomt + 2*ADMalpha*ADMdalpha2*Power(pmomt,2))/pmomt
;

dtpmom3
=
-0.5*(dinvgam113*Power(pmom1,2) + 2*dinvgam123*pmom1*pmom2 +
     dinvgam223*Power(pmom2,2) + 2*dinvgam133*pmom1*pmom3 +
     2*dinvgam233*pmom2*pmom3 + dinvgam333*Power(pmom3,2) -
     2*ADMdbeta13*pmom1*pmomt - 2*ADMdbeta23*pmom2*pmomt -
     2*ADMdbeta33*pmom3*pmomt + 2*ADMalpha*ADMdalpha3*Power(pmomt,2))/pmomt
;


}

} // namespace Particles

#endif // #ifndef PARTICLES_GEODESIC_HXX

/* particles_geodesic.hxx */
