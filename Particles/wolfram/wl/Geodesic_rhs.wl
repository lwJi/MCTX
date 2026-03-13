(* ::Package:: *)

(* Geodesic_rhs.wl *)

(* (c) Liwei Ji, 03/2025 *)

Module[{Mat, invMat},
  Mat = Table[ADMgam[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> detinvgam};
  SetEQNDelayed[detinvgam[], 1 / Det[Mat] // Simplify];
  SetEQNDelayed[invgam[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify];
];

SetEQN[dinvgam[i_, j_, k_], -invgam[i, l] invgam[j, m] ADMdgam[-l, -m, k]];

SetEQN[p2[], invgam[i, j] pmom[-i] pmom[-j]];

SetEQN[pt[], (p2[]) ^ (1/2) / ADMalpha[]];

SetEQN[dtxpos[i_], invgam[i, j] pmom[-j] / pt[] - ADMbeta[i]];

SetEQN[dtpmom[i_], -ADMalpha[] pt[] ADMdalpha[i] + pmom[-j] ADMdbeta[j, i] - 1/(2 pt[]) pmom[-j] pmom[-k] dinvgam[j, k, i]];

SetEQN[pmomt[], pmom[-i] ADMbeta[i] - ADMalpha[] (p2[]) ^ (1/2)];

