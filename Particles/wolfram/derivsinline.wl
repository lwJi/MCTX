(* ::Package:: *)

(* derivsinline.wl *)

(* (c) Liwei Ji, 02/2025 *)

(******************)
(* Configurations *)
(******************)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"}]]

SetPVerbose[False];

SetPrintDate[False];

PrintIndexes3DAMReX[accuracyOrd_?IntegerQ, fdOrd_?IntegerQ, strDir_?StringQ] :=
  Module[{stencils, solution},
    (* Get stencils and finite difference coefficients *)
    stencils = GetCenteringStencils[accuracyOrd];
    solution = GetFiniteDifferenceCoefficients[stencils, fdOrd];
    (* Construct the finite difference index expression *)
    Do[
      index = stencils[[i]];
      If[(Subscript[c, index] /. solution) == 0, Continue[]];
      buf = "  const T " <> ToString[GetGFIndexName[index]] <>
        If[index == 0,
          " = gf(i, j, k);"
          ,
          " = gf("
            <> "i + (D == 0 ? " <> ToString[index] <> " : 0), "
            <> "j + (D == 1 ? " <> ToString[index] <> " : 0), "
            <> "k + (D == 2 ? " <> ToString[index] <> " : 0));"
        ];
      pr[buf]
      ,
      {i, 1, Length[stencils]}
    ];
  ];

PrintFDExpressionAMReX[accuracyOrd_?IntegerQ, fdOrd_?IntegerQ, strIdx_?StringQ] :=
  Module[{stencils, solution, buf, rule},
    (* Rules for string replacements *)
    rule = {
      "invdx" -> strIdx <> "[D]"
    };
    (* Get stencils and finite difference coefficients *)
    stencils = GetCenteringStencils[accuracyOrd];
    solution = GetFiniteDifferenceCoefficients[stencils, fdOrd];
    (* Construct the finite difference expression *)
    buf = "    " <> ToString[CForm[
      (Sum[
        index = stencils[[i]];
        (Subscript[c, index] /. solution) GetGFIndexName[index],
        {i, 1, Length[stencils]}] // Simplify)
      Product[invdx, {i, 1, fdOrd}]
    ]] <> ";";
    pr[StringReplace[buf, rule]];
  ];

(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "derivsinline.hxx"}]];

SetMainPrint[
  pr["#include <loop_device.hxx>"];
  pr[];
  pr["#include <array>"];
  pr["#include <cmath>"];
  pr[];
  pr["#include \"powerinline.hxx\""];
  pr[];

  pr["namespace Particles {"];
  pr["using namespace Loop;"];
  pr[];

  (****************************************)
  (* First-order Finite Difference Scheme *)
  (****************************************)

  Do[
    pr["template <int D, typename T>"];
    pr["CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T"];
    pr["fd_1_o" <> ToString[aOrd]
                <> "(amrex::Array4<T const> const &gf, "
                <> "int i, int j, int k, int comp, "
                <> "amrex::GpuArray<T, 3> const &dxi) {"];
    PrintIndexes3DAMReX[aOrd, 1, "D"];
    pr["  return"];
    PrintFDExpressionAMReX[aOrd, 1, "dxi"];
    pr["}"];
    pr[]
    ,
    {aOrd, 2, 8, 2}
  ];

  pr["} // namespace Particles"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
