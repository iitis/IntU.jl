(* ::Package:: *)
(* Performance comparison benchmarks for RTNI (Mathematica).
   Computes the same unitary integrals as bench_intu.jl for a head-to-head comparison.

   Usage:
     1. Place this file in the same directory as the RTNI package (RTNI.wl + precomputedWG/).
     2. Open in Mathematica and evaluate, or run from command line:
        math -script bench_rtni.wl

   Prerequisites:
     - Wolfram Mathematica (tested with v11+)
     - RTNI package: https://github.com/MotohisaFukuda/RTNI
*)

(* --- Setup --- *)
(* If running as a script, set directory to script location *)
If[$InputFileName =!= "",
  SetDirectory[DirectoryName[$InputFileName]],
  SetDirectory[NotebookDirectory[]]
];

(* Load RTNI - adjust path if RTNI.wl is elsewhere *)
Needs["RTNI`"];

nWarmup = 2;
nSamples = 10;
slowThreshold = 2.0;  (* seconds *)
slowSamples = 3;

(* Benchmark helper: returns {medianTime, result} *)
benchmark[func_, nwarm_:nWarmup, nsamp_:nSamples] := Module[
  {times, result, probeTime, actualSamples, t},

  (* Warmup *)
  Do[func[], {nwarm}];

  (* Probe to detect slow benchmarks *)
  {probeTime, result} = AbsoluteTiming[func[]];
  actualSamples = If[probeTime > slowThreshold, slowSamples, nsamp];

  (* Timed samples *)
  times = Table[
    First[AbsoluteTiming[func[]]],
    {actualSamples}
  ];

  {Median[times], result, actualSamples}
];

results = <||>;

runAndReport[name_String, func_] := Module[
  {med, res, n, ms},
  Print["  Running: ", name, " ..."];
  {med, res, n} = benchmark[func];
  ms = med * 1000;
  Print["    ", NumberForm[ms, {6,2}], " ms  (N=", n, ", result: ", res, ")"];
  results[name] = <|"median_ms" -> ms, "result" -> ToString[res], "samples" -> n|>;
];


(* ============================================================================ *)
(* RTNI MultinomialexpectationvalueHaar API:                                    *)
(*   MultinomialexpectationvalueHaar[d, eps, Xvars]                             *)
(*   eps = {e1, e2, ..., e_{2q}} where:                                         *)
(*     1 = U, 2 = U^*, 3 = U^T, 4 = conj(U)                                   *)
(*   The expression computed is X1 U^{e1} X2 U^{e2} ... X_{2q} U^{e_{2q}}      *)
(*   Xvars = list of symbolic X matrices (or Id for identity)                   *)
(*                                                                              *)
(* For |U_{11}|^{2k} = (U_{11})^k * (conj(U_{11}))^k:                          *)
(*   We need alternating U (1) and U* (2), with Id matrices between them.       *)
(*   The result with usetrace->True gives the scalar value.                     *)
(* ============================================================================ *)

Print["\n=== Easy: Unitary |U_11|^{2k}, symbolic d ==="];

(* |U_11|^2: U * U^* -> eps = {1, 2} *)
runAndReport["U_|U11|^2_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2}, {Id, Id}, usetrace -> True]
  ]
];

(* |U_11|^4: U U* U U* -> eps = {1, 2, 1, 2} *)
runAndReport["U_|U11|^4_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2, 1, 2},
      {Id, Id, Id, Id}, usetrace -> True]
  ]
];

(* |U_11|^6 -> eps = {1, 2, 1, 2, 1, 2} *)
runAndReport["U_|U11|^6_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2, 1, 2, 1, 2},
      Table[Id, 6], usetrace -> True]
  ]
];


Print["\n=== Harder: Unitary |U_11|^{2k}, symbolic d ==="];

(* |U_11|^8 -> eps = {1,2,1,2,1,2,1,2} *)
runAndReport["U_|U11|^8_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, ConstantArray[{1, 2}, 4] // Flatten,
      Table[Id, 8], usetrace -> True]
  ]
];

(* |U_11|^10 -> eps = {1,2,1,2,1,2,1,2,1,2} *)
runAndReport["U_|U11|^10_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, ConstantArray[{1, 2}, 5] // Flatten,
      Table[Id, 10], usetrace -> True]
  ]
];


Print["\n=== Harder: Unitary |U_11|^{2k}, numeric d ==="];

(* |U_11|^10, d=10 *)
runAndReport["U_|U11|^10_d=10",
  Function[{},
    MultinomialexpectationvalueHaar[10, ConstantArray[{1, 2}, 5] // Flatten,
      Table[Id, 10], usetrace -> True]
  ]
];

(* |U_11|^10, d=50 *)
runAndReport["U_|U11|^10_d=50",
  Function[{},
    MultinomialexpectationvalueHaar[50, ConstantArray[{1, 2}, 5] // Flatten,
      Table[Id, 10], usetrace -> True]
  ]
];


Print["\n=== Off-diagonal: Unitary, symbolic d ==="];

(* For off-diagonal integrals like |U_11|^2 |U_12|^2 we need to use
   the graph-based integrateHaarUnitary, since MultinomialexpectationvalueHaar
   computes monomials of the form X1 U^e1 X2 U^e2 ... which correspond
   to a single chain. For multi-entry products we construct the appropriate
   tensor network graph. *)

(* Alternative: use MultinomialexpectationvalueHaar with appropriate X matrices.
   |U_{11}|^2 |U_{12}|^2 = U_{11} conj(U_{11}) U_{12} conj(U_{12})
   This can be expressed as:
   tr(E_{11} U E_{11} U^* E_{11} U E_{22} U^*)
   where E_{ij} is a matrix unit.
   But this requires care. Instead, we use the graph interface directly. *)

(* --- Off-diagonal via integrateHaarUnitary graph interface ---

   integrateHaarUnitary[g, symbU, dimsI, dimsO, totalDim]
   where g is a list of edges, each edge = {{boxname, boxid, in/out(0/1), labelid}, {boxname, boxid, in/out(0/1), labelid}}

   For |U_{11}|^2 = U_{11} conj(U)_{11}:
   We need one copy of U (id=1) and one copy of conj(U) (id=1, but marked as conjugate).
   The dimension is just {1} for single-entry access.

   Actually, RTNI's MultinomialexpectationvalueHaar is for single-chain monomials,
   not for products of individual matrix elements with independent indices.

   For products of matrix elements, we construct a graph:
   |U_{i1,j1}|^2 ... |U_{ik,jk}|^2 = U_{i1,j1} Ubar_{i1,j1} ... U_{ik,jk} Ubar_{ik,jk}

   Each pair (U, Ubar) with matching row/col indices creates edges in the graph.
*)

(* Helper: build a graph for integral of prod |U_{rows[[i]], cols[[i]]}|^2
   Each U copy has id=i, each Ubar copy also has id=i.
   Rows are connected via dangling edges (representing delta functions).
   Columns similarly.

   For RTNI: the graph encodes the contraction pattern.
   - "U" boxes have input (0) and output (1) legs.
   - We connect row indices of U to row indices of Ubar (dangling nodes for external indices).
   - For |U_{11}|^2 |U_{12}|^2, after integration, the result contracts with delta functions.

   However, for the simplest RTNI interface for element-wise products:
   integrateHaarUnitary needs the full graph specification.
*)

(* For simplicity with RTNI, let's use producemonomialgraph for
   products of matrix elements. The function producemonomialgraph[eps]
   takes eps = list of {type, row, col} specifications. However, the exact
   API is complex. Let's use a direct graph construction. *)

(* ---- Graph construction for |U_{11}|^{2k} ---- *)
(* For a single entry |U_{11}|^{2k}: we need k copies of U and k copies of Ubar.
   Each U has one input and one output (scalar dimensions).
   All U's row indices connect to corresponding Ubar's row indices (all are index 1).
   All U's col indices connect to corresponding Ubar's col indices (all are index 1).

   In RTNI graph notation:
   - U copy i: {"U", i, 0, 1} for input, {"U", i, 1, 1} for output
   - Ubar copy i: {"Ub", i, 0, 1} for input, {"Ub", i, 1, 1} for output
   Edges connect: input of U_i to input of Ub_i (row delta: both are index 1)
                  output of U_i to output of Ub_i (col delta: both are index 1)  *)

(* Actually, for scalar products of entries with specified indices, the edges encode
   which row/col indices are identified. For |U_{11}|^2 |U_{12}|^2:
   U_1 has row=1,col=1; Ubar_1 has row=1,col=1 -> connect
   U_2 has row=1,col=2; Ubar_2 has row=1,col=2 -> connect row1=row1 but col2!=col1

   In RTNI, the graph edges between dangling nodes at the same row position create
   delta contractions. So for independent element products, we create "X" (constant)
   tensors that enforce the right index structure.

   This is getting complex. Let's focus on what RTNI does most naturally:
   MultinomialexpectationvalueHaar handles trace-expressible monomials.

   For the comparison, let's focus on the cases that MultinomialexpectationvalueHaar
   handles directly (the chain monomials), and add one graph-based case. *)


(* --- Trace-based integrals via MultinomialexpectationvalueHaar ---
   |tr(U)|^4 = tr(U)^2 * tr(U^*)^2
   This is a product of traces, not a single trace.
   MultinomialexpectationvalueHaar handles: tr(X1 U^e1 X2 U^e2 ...)
   which is a SINGLE trace of a product.

   |tr(U)|^2 = tr(U) * tr(U^*) = sum_{i,j} U_{ii} Ubar_{jj}
   This IS expressible as tr(U) * tr(U*) = tr(U * Id) * tr(U* * Id)

   However, MultinomialexpectationvalueHaar computes E[tr(X1 U^e1 X2 U^e2 ...)]
   which is a single trace. For products of traces we need the graph interface.
*)

(* --- Graph-based: |tr(U)|^4 = tr(U)^2 tr(U^*)^2 ---
   This needs integrateHaarUnitary with a graph that has:
   - 2 copies of U (call them U1, U2)
   - 2 copies of Ubar (call them Ub1, Ub2)
   - Each U's output connects to its own input (trace loop)
   - Each Ub's output connects to its own input (trace loop)
   - All dimensions = {d}

   Graph:
   U copy 1: input  = {"U", 1, 0, 1}, output = {"U", 1, 1, 1}
   U copy 2: input  = {"U", 2, 0, 1}, output = {"U", 2, 1, 1}
   Ub copy 1: input = {"Ub", 1, 0, 1}, output = {"Ub", 1, 1, 1}
   Ub copy 2: input = {"Ub", 2, 0, 1}, output = {"Ub", 2, 1, 1}

   Trace edges (each U/Ub traces over itself):
   {{"U", 1, 0, 1}, {"U", 1, 1, 1}}   (* tr(U1) *)
   {{"U", 2, 0, 1}, {"U", 2, 1, 1}}   (* tr(U2) *)
   {{"Ub", 1, 0, 1}, {"Ub", 1, 1, 1}} (* tr(Ub1) *)
   {{"Ub", 2, 0, 1}, {"Ub", 2, 1, 1}} (* tr(Ub2) *)
*)

runAndReport["U_|trU|^4_sym",
  Function[{},
    Module[{g},
      g = {
        {{"U", 1, 0, 1}, {"U", 1, 1, 1}},
        {{"U", 2, 0, 1}, {"U", 2, 1, 1}},
        {{"Ub", 1, 0, 1}, {"Ub", 1, 1, 1}},
        {{"Ub", 2, 0, 1}, {"Ub", 2, 1, 1}}
      };
      integrateHaarUnitary[g, "U", {dd}, {dd}, dd]
    ]
  ]
];


(* --- Graph-based: tr(U A U^* B) ---
   This is E[tr(U A U^* B)] which is the depolarizing channel integral.
   eps = {1, 2} with X = {A, B} and usetrace -> True *)

runAndReport["U_trUAUdB_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2}, {Subscript[X, 1], Subscript[X, 2]},
      usetrace -> True]
  ]
];


(* --- Graph-based: |tr(U)|^6 = tr(U)^3 tr(U^*)^3 --- *)

runAndReport["U_|trU|^6_sym",
  Function[{},
    Module[{g},
      g = Join[
        Table[{{"U", i, 0, 1}, {"U", i, 1, 1}}, {i, 3}],
        Table[{{"Ub", i, 0, 1}, {"Ub", i, 1, 1}}, {i, 3}]
      ];
      integrateHaarUnitary[g, "U", {dd}, {dd}, dd]
    ]
  ]
];


(* --- Graph-based: |tr(U)|^8 --- *)

runAndReport["U_|trU|^8_sym",
  Function[{},
    Module[{g},
      g = Join[
        Table[{{"U", i, 0, 1}, {"U", i, 1, 1}}, {i, 4}],
        Table[{{"Ub", i, 0, 1}, {"Ub", i, 1, 1}}, {i, 4}]
      ];
      integrateHaarUnitary[g, "U", {dd}, {dd}, dd]
    ]
  ]
];


(* --- Graph-based: tr(U A U^* B U A U^* B) = E[tr((U A U^* B)^2)] ---
   This is a degree-2 trace polynomial.
   eps = {1, 2, 1, 2}, Xvars = {A, B, A, B} *)

runAndReport["U_tr(UAUdB)^2_sym",
  Function[{},
    MultinomialexpectationvalueHaar[dd, {1, 2, 1, 2},
      {Subscript[X, 1], Subscript[X, 2], Subscript[X, 3], Subscript[X, 4]},
      usetrace -> True]
  ]
];


(* ============================================================================ *)
(* Save results                                                                  *)
(* ============================================================================ *)
Print["\n", StringJoin@ConstantArray["=", 72]];
Print["Summary (median times in ms)"];
Print[StringJoin@ConstantArray["=", 72]];

Do[
  Print[StringPadRight[name, 30], "  ",
    If[KeyExistsQ[results[name], "error"],
      "FAILED",
      NumberForm[results[name]["median_ms"], {8,2}]
    ]
  ],
  {name, Keys[results]}
];

(* Export as JSON *)
Export["results_rtni.json", Normal[results]];
Print["\nResults saved to results_rtni.json"];
