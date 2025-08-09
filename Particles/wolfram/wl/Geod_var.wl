(* ::Package:: *)

(* Geod_vars.wl *)

(* (c) Liwei Ji, 08/2025 *)


ADMVarlist =
  TempTensors[
    {alpha[],                          PrintAs -> "\[Alpha]"},
    {beta[i],                          PrintAs -> "\[Beta]"},
    {gam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Gamma]"}
  ];

dADMVarlist =
  TempTensors[
    {dalpha[-k],                            PrintAs -> "\[PartialD]\[Alpha]"},
    {dbeta[i, -k],                          PrintAs -> "\[PartialD]\[Beta]"},
    {dgam[-i, -j, -k], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[Gamma]"}
  ];

TempVarlist =
  TempTensors[
    {detinvgam[],                          PrintAs -> "1/\[Gamma]"},
    {invgam[i, j], Symmetric[{i, j}],      PrintAs -> "\[Gamma]"},
    {dinvgam[i, j, -k], Symmetric[{i, j}], PrintAs -> "\[PartialD]\[Gamma]"},
    {p2[],    PrintAs -> "\!\(\*SuperscriptBox[\(p\), \(2\)]\)"},
    {pmomt[], PrintAs -> "\!\(\*SuperscriptBox[\(p\), \(t\)]\)"}
  ];

ParticleVarlist =
  GridTensors[
    {xpos[i],  PrintAs -> "x"},
    {pmom[-i], PrintAs -> "p"}
  ];

dtParticleVarlist =
  GridTensors[
    {dtxpos[i],  PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)x"},
    {dtpmom[-i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)p"}
  ];
