(* ::Package:: *)
(* Performance comparison benchmarks for RTNI (Mathematica).
   Computes the same unitary integrals as bench_intu.jl.

   Usage:
     1. Place this file alongside the RTNI package (RTNI.wl + precomputedWG/).
     2. math -script bench_rtni.wl

   Prerequisites: Wolfram Mathematica with RTNI package.
*)

(* --- Setup --- *)
If[$InputFileName =!= "",
  SetDirectory[DirectoryName[$InputFileName]],
  SetDirectory[NotebookDirectory[]]
];

Needs["RTNI`"];

nWarmup = 2;
nSamples = 10;
slowThreshold = 2.0;
slowSamples = 3;

(* Benchmark helper: returns {medianTime, result, nSamples} *)
benchmark[func_, nwarm_:nWarmup, nsamp_:nSamples] := Module[
  {times, result, probeTime, actualSamples},
  Do[func[], {nwarm}];
  {probeTime, result} = AbsoluteTiming[func[]];
  actualSamples = If[probeTime > slowThreshold, slowSamples, nsamp];
  times = Table[First[AbsoluteTiming[func[]]], {actualSamples}];
  {Median[times], result, actualSamples}
];

results = <||>;

runAndReport[name_String, func_] := Module[
  {med, res, n, ms},
  Print["  Running: ", name, " ..."];
  {med, res, n} = benchmark[func];
  ms = med * 1000;
  Print["    ", NumberForm[ms, {6, 2}], " ms  (N=", n, ", result: ", res, ")"];
  results[name] = <|"median_ms" -> ms, "result" -> ToString[res], "samples" -> n|>;
];


(* ============================================================================ *)
(* Element integrals: E[|U_{1,1}|^{2k}] via MultinomialexpectationvalueHaar    *)
(*                                                                              *)
(* |U_{1,1}|^{2k} = U_{1,1}^k * Conjugate[U_{1,1}]^k                          *)
(*   = tr(E_{11} U E_{11} U ... E_{11} U\[Dagger] E_{11} U\[Dagger] ...)       *)
(* where E_{11} is the rank-1 projector |1><1|.                                *)
(*                                                                              *)
(* We use eps = {1,...,1, 2,...,2} (k ones then k twos)                         *)
(* with Xvars = {x, x, ..., x} (symbolic projector).                           *)
(* Since E_{11} is a rank-1 projector, Tr[E_{11}^n] = 1 for all n,             *)
(* so we substitute all Tr terms to 1 in the result.                            *)
(* ============================================================================ *)

(* Helper: compute E[|U_{1,1}|^{2k}] with symbolic dimension dd *)
u11moment[k_] := Module[{eps, xvars, raw},
  eps = Join[ConstantArray[1, k], ConstantArray[2, k]];
  xvars = Table[x, {2 k}];
  raw = MultinomialexpectationvalueHaar[dd, eps, xvars, usetrace -> True];
  (* E_{11} is a rank-1 projector: all traces of products of x evaluate to 1 *)
  raw /. Tr[__] -> 1
];

(* Helper: compute E[|U_{1,1}|^{2k}] with numeric dimension dval *)
u11momentNumeric[k_, dval_] := Module[{eps, xvars, e11},
  eps = Join[ConstantArray[1, k], ConstantArray[2, k]];
  e11 = SparseArray[{1, 1} -> 1, {dval, dval}];
  xvars = Table[e11, {2 k}];
  MultinomialexpectationvalueHaar[dval, eps, xvars, usetrace -> True]
];


Print["\n=== Easy: Unitary |U_11|^{2k}, symbolic d ==="];

runAndReport["U_|U11|^2_sym", Function[{}, u11moment[1]]];
runAndReport["U_|U11|^4_sym", Function[{}, u11moment[2]]];
runAndReport["U_|U11|^6_sym", Function[{}, u11moment[3]]];


Print["\n=== Harder: Unitary |U_11|^{2k}, symbolic d ==="];

runAndReport["U_|U11|^8_sym", Function[{}, u11moment[4]]];
runAndReport["U_|U11|^10_sym", Function[{}, u11moment[5]]];


Print["\n=== Harder: Unitary |U_11|^{2k}, numeric d ==="];

runAndReport["U_|U11|^10_d=10", Function[{}, u11momentNumeric[5, 10]]];
runAndReport["U_|U11|^10_d=50", Function[{}, u11momentNumeric[5, 50]]];


(* ============================================================================ *)
(* Trace moments: E[|tr(U)|^{2k}] via integrateHaarUnitary graph interface      *)
(*                                                                              *)
(* |tr(U)|^{2k} = tr(U)^k tr(U\[Dagger])^k                                    *)
(* Graph: k copies of U with self-loops (trace), k copies of Ub with self-loops *)
(* ============================================================================ *)

trMomentGraph[k_] := Module[{g},
  g = Join[
    Table[{{"U", i, 0, 1}, {"U", i, 1, 1}}, {i, k}],
    Table[{{"Ub", i, 0, 1}, {"Ub", i, 1, 1}}, {i, k}]
  ];
  integrateHaarUnitary[g, "U", {dd}, {dd}, dd]
];


Print["\n=== Trace moments: |tr(U)|^{2k}, symbolic d ==="];

runAndReport["U_|trU|^4_sym", Function[{}, trMomentGraph[2]]];
runAndReport["U_|trU|^6_sym", Function[{}, trMomentGraph[3]]];
runAndReport["U_|trU|^8_sym", Function[{}, trMomentGraph[4]]];


(* ============================================================================ *)
(* Trace polynomials via MultinomialexpectationvalueHaar                        *)
(*                                                                              *)
(* E[tr(U A U\[Dagger] B)] : eps = {1, 2}, Xvars = {A, B}                      *)
(* E[tr(U A U\[Dagger] B U C U\[Dagger] D)] : eps = {1,2,1,2}, Xvars={A,B,C,D} *)
(* ============================================================================ *)

Print["\n=== Trace polynomials: symbolic d ==="];

runAndReport["U_trUAUdB_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2},
      {Subscript[X, 1], Subscript[X, 2]}, usetrace -> True]
  ]
];

runAndReport["U_tr(UAUdB)^2_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2, 1, 2},
      {Subscript[X, 1], Subscript[X, 2], Subscript[X, 3], Subscript[X, 4]},
      usetrace -> True]
  ]
];


(* ============================================================================ *)
(* Save results                                                                 *)
(* ============================================================================ *)
Print["\n", StringJoin@ConstantArray["=", 72]];
Print["Summary (median times in ms)"];
Print[StringJoin@ConstantArray["=", 72]];

Do[
  Print[StringPadRight[name, 30], "  ",
    NumberForm[results[name]["median_ms"], {8, 2}]
  ],
  {name, Keys[results]}
];

Export["results_rtni.json", Normal[results]];
Print["\nResults saved to results_rtni.json"];
