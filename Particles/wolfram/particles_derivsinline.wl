(* ::Package:: *)

(* particles_derivsinline.wl *)

(* (c) Liwei Ji, 02/2025 *)

(******************)
(* Configurations *)
(******************)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"}]]

SetPVerbose[False];

SetPrintDate[False];

(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "particles_derivsinline.hxx"}]];

SetMainPrint[
  pr["#include <loop_device.hxx>"];
  pr[];
  pr["#include <array>"];
  pr["#include <cmath>"];
  pr[];
  pr["#include \"particles_powerinline.hxx\""];
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
    PrintIndexes3D[aOrd, 1, "D"];
    pr["  return"];
    PrintFDExpression[aOrd, 1, "dxi"];
    pr["}"];
    pr[]
    ,
    {aOrd, 2, 8, 2}
  ];

  pr["} // namespace Particles"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/AMReX.wl"}]];
