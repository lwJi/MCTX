(* ::Package:: *)

(* Geod_rhs.wl *)

(* (c) Liwei Ji, 08/2025 *)

SetEQN[dinvgam[i_, j_, k_], -invgam[i, l] invgam[j, m] dgam[-l, -m, k]];

SetEQN[p2[], invgam[i, j] pmom[-i] pmom[-j]];

SetEQN[pmomt[], (p2[]) ^ (1/2) / alpha[]];

SetEQN[dtxpos[i_], invgam[i, j] pmom[-j] / pmomt[] - beta[i]];

SetEQN[dtpmom[i_], -alpha[] pmomt[] dalpha[i] + pmom[-j] dbeta[j, i] - 1/2 pmom[-j] pmom[-k] dinvgam[j, k, i]];
