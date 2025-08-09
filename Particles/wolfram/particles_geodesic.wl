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

(* Define Variables *)

<<wl/Geod_var.wl

<<wl/Geod_rhs.wl

Module[{Mat, invMat},
  Mat = Table[gam[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> detinvgam};
  SetEQNDelayed[detinvgam[], 1 / Det[Mat] // Simplify];
  SetEQNDelayed[invgam[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];

(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "particles_geodesic.hxx"}]];

SetMainPrint[

  pr["#include \"particles_powerinline.hxx\""];
  pr[];

  pr["namespace Particles {"];
  pr[];

  pr["template <typename T>"];
  pr["CCTK_HOST CCTK_DEVICE inline void"];
  pr["calc_rhs_geodesic(std::array<T, 3> &dtxpos, std::array<T, 3> &dtpmom,"];
  pr["                  const std::array<T, 3> &pmom,"];
  pr["                  const T alpha,"];
  pr["                  const std::array<T, 3> &beta,"];
  pr["                  const std::array<T, 6> &gam,"];
  pr["                  const std::array<T, 3> &dalpha,"];
  pr["                  const std::array<std::array<T, 3>, 3> &dbeta,"];
  pr["                  const std::array<std::array<T, 3>, 6> &dgam) {"];
  pr[];

  PrintInitializations[{Mode -> "MainOut"}, dtParticleVarlist];
  pr[];

  PrintInitializations[{Mode -> "MainIn"}, Drop[ParticleVarlist, 1]];
  pr[];

  PrintInitializations[{Mode -> "Temp"}, Drop[ADMVarlist, 1]];
  pr[];

  PrintInitializations[{Mode -> "Temp", TensorType -> "Vect"}, dADMVarlist];
  pr[];

  PrintEquations[{Mode -> "Temp", ExtraReplaceRules -> {Sqrt[p2] -> sqrt[p2]}}, TempVarlist];
  pr[];

  PrintEquations[{Mode -> "Main"}, dtParticleVarlist];
  pr[];

  pr["}"];
  pr[];

  pr["} // namespace Particles"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/AMReX.wl"}]];
