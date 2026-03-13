(* ::Package:: *)

(* particles_geodesic.wl *)

(* (c) Liwei Ji, 07/2024 *)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"}]]

SetPVerbose[False];

SetPrintDate[False];

(*SetPrintHeaderMacro[False];*)

SetGridPointIndex[""];

SetTempVariableType["auto"];

DefManifold[M3, 3, IndexRange[a, z]];

DefChart[cart, M3, {1, 2, 3}, {X[], Y[], Z[]}, ChartColor -> Blue];

(**********************************)
(* Define Variables and Equations *)
(**********************************)

Get[FileNameJoin[{Environment["GENERATO"], "module/ADM/ADM_var.wl"}]];

Get["wl/Geodesic_var.wl"];

Get["wl/Geodesic_rhs.wl"];

(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "particles_geodesic.hxx"}]];

SetMainPrint[

  pr["#include \"util_powerinline.hxx\""];
  pr["#include \"../src/Particles.hxx\""];
  pr[];

  pr["namespace Particles {"];
  pr["using namespace UtilForge;"];
  pr[];

  pr["CCTK_HOST CCTK_DEVICE inline void"];
  pr["calc_rhs_geodesic(VectR &dtpmom, VectR &dtxpos, const VectR &pmom,"];
  pr["                  const ScalR ADMalpha, const VectR &ADMbeta, const SmatR &ADMgam,"];
  pr["                  const dScalR &ADMdalpha, const dVectR &ADMdbeta, const dSmatR &ADMdgam) {"];
  pr[];

  PrintInitializations[{Mode -> "MainOut"}, dtParticleVarlist];
  pr[];

  PrintInitializations[{Mode -> "Temp"}, ParticleVarlist[[2;;2]]];
  pr[];

  PrintInitializations[{Mode -> "Temp"}, Drop[ADMVarlist, {2, 3}]];
  pr[];

  PrintInitializations[{Mode -> "Temp", TensorType -> "Vect"}, ADMdVarlist];
  pr[];

  PrintEquations[{Mode -> "Temp", ExtraReplaceRules -> {Sqrt[p2] -> sqrt[p2]}}, TempVarlist];
  pr[];

  PrintEquations[{Mode -> "MainOut"}, dtParticleVarlist];
  pr[];

  pr["}"];
  pr[];

  pr["CCTK_HOST CCTK_DEVICE inline void"];
  pr["calc_pt(ScalR &pmomt, const VectR &pmom,"];
  pr["        const ScalR ADMalpha, const VectR &ADMbeta, const SmatR &ADMgam) {"];
  pr[];

  PrintInitializations[{Mode -> "Temp"}, Drop[ADMVarlist, {2, 3}]];
  pr[];

  PrintInitializations[{Mode -> "Temp"}, ParticleVarlist[[2;;2]]];
  pr[];

  PrintEquations[{Mode -> "Temp"}, Drop[Drop[TempVarlist, {3}], {-1}]];
  pr[];

  PrintEquations[{Mode -> "MainOut", ExtraReplaceRules -> {Sqrt[p2] -> sqrt[p2]}}, ParticleVarlist[[3;;3]]];
  pr[];

  pr["}"];
  pr[];

  pr["} // namespace Particles"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/AMReX.wl"}]];
