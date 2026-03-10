(* ::Package:: *)

(* Geodesic_var.wl *)

(* (c) Liwei Ji, 03/2026 *)

(**********************)

(* Particle Variables *)

(**********************)

TempVarlist =
  TempTensors[
    {detinvgam[],                          PrintAs -> "1/\[Gamma]"},
    {invgam [i, j],     Symmetric[{i, j}], PrintAs -> "\[Gamma]"},
    {dinvgam[i, j, -k], Symmetric[{i, j}], PrintAs -> "\[PartialD]\[Gamma]"},
    {p2[],    PrintAs -> "\!\(\*SuperscriptBox[\(p\), \(2\)]\)"},
    {pmomt[], PrintAs -> "\!\(\*SuperscriptBox[\(p\), \(t\)]\)"}
  ];

ParticleVarlist =
  GridTensors[
    {xpos [i], PrintAs -> "x"},
    {pmom[-i], PrintAs -> "p"}
  ];

dtParticleVarlist =
  GridTensors[
    {dtxpos [i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)x"},
    {dtpmom[-i], PrintAs -> "\!\(\*SubscriptBox[\(\[PartialD]\), \(t\)]\)p"}
  ];
