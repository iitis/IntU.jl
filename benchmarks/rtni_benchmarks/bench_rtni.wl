(* ::Package:: *)
(* Performance comparison benchmarks for RTNI (Mathematica).
   Computes the same unitary integrals as bench_intu.jl.

   Usage:
     1. Place this file alongside the RTNI package (RTNI.wl + precomputedWG/).
     2. math -script bench_rtni.wl

   Prerequisites: Wolfram Mathematica with RTNI package.
*)

(* --- Setup --- *)
(* When running as a script (math -script), NotebookDirectory[] is unavailable.
   RTNI.wl internally uses NotebookDirectory[] to locate precomputedWG/.
   We override it so RTNI can find its files. *)
If[$InputFileName =!= "",
  scriptDir = DirectoryName[$InputFileName],
  scriptDir = NotebookDirectory[]
];
SetDirectory[scriptDir];

(* Override NotebookDirectory so RTNI can find precomputedWG/ *)
Unprotect[NotebookDirectory];
NotebookDirectory[] := scriptDir;
Protect[NotebookDirectory];

Get["RTNI.wl"];

rtniPath = ExpandFileName["RTNI.wl"];
rtniHash = If[
  FileExistsQ[rtniPath],
  IntegerString[FileHash[rtniPath, "SHA256"], 16, 64],
  "missing"
];
scriptPath = If[$InputFileName =!= "", ExpandFileName[$InputFileName], "interactive"];
timestampUTC = DateString[
  TimeZoneConvert[Now, 0],
  {"Year", "-", "Month", "-", "Day", "T", "Hour", ":", "Minute", ":", "Second", "Z"}
];

nWarmup = 2;
nSamples = 10;
slowThreshold = 2.0;
slowSamples = 3;
rowTimeoutSec = 120;

projectorSymbols = {p11, p22};

buildProjector[dval_, i_] := SparseArray[{i, i} -> 1, {dval, dval}];

buildAlternatingProjectors[rowProjectors_List, colProjectors_List] :=
  Flatten @ Transpose[{rowProjectors, colProjectors}];

projectorTraceValue[args__] := Module[{lst = {args}},
  If[Length[lst] == 0, 0, If[SameQ @@ lst, 1, 0]]
];

simplifyProjectorTraces[expr_] := expr //. {
  Tr[Dot[args__]] :> projectorTraceValue[args],
  Tr[s_Symbol] /; MemberQ[projectorSymbols, s] :> 1,
  Tr[Power[s_Symbol, n_Integer?Positive]] /; MemberQ[projectorSymbols, s] :> 1
};

uMomentFromProjectors[dim_, rowProjectors_List, colProjectors_List] := Module[
  {k, eps, xvars, raw},
  k = Length[rowProjectors];
  eps = Join[ConstantArray[1, k], ConstantArray[2, k]];
  xvars = buildAlternatingProjectors[rowProjectors, colProjectors];
  raw = MultinomialexpectationvalueHaar[dim, eps, xvars, usetrace -> True];
  raw
];

uMomentSymbolic[rowProjectors_List, colProjectors_List] :=
  simplifyProjectorTraces[uMomentFromProjectors[dd, rowProjectors, colProjectors]];

uMomentNumeric[rowProjectors_List, colProjectors_List, dval_] := Module[
  {rowMats, colMats},
  rowMats = rowProjectors /. {
    p11 -> buildProjector[dval, 1],
    p22 -> buildProjector[dval, 2]
  };
  colMats = colProjectors /. {
    p11 -> buildProjector[dval, 1],
    p22 -> buildProjector[dval, 2]
  };
  uMomentFromProjectors[dval, rowMats, colMats]
];

u11moment[k_] := uMomentSymbolic[ConstantArray[p11, k], ConstantArray[p11, k]];

u11momentNumeric[k_, dval_] :=
  uMomentNumeric[ConstantArray[p11, k], ConstantArray[p11, k], dval];

mixedCaseDefinitions = <|
  "U_|U11|^2|U12|^2" -> <|
    "rows" -> {p11, p11},
    "cols" -> {p11, p22}
  |>,
  "U_|U11|^2|U12|^4" -> <|
    "rows" -> {p11, p11, p11},
    "cols" -> {p11, p22, p22}
  |>,
  "U_|U11|^2|U22|^2" -> <|
    "rows" -> {p11, p22},
    "cols" -> {p11, p22}
  |>
|>;

withTimeLimit[func_] := Module[{timed},
  timed = AbsoluteTiming[TimeConstrained[func[], rowTimeoutSec, $Aborted]];
  If[
    timed[[2]] === $Aborted,
    <|"ok" -> False, "status" -> "timeout", "time_s" -> Null, "result" -> "TIMEOUT"|>,
    <|"ok" -> True, "status" -> "ok", "time_s" -> timed[[1]], "result" -> timed[[2]]|>
  ]
];

benchmark[func_, nwarm_:nWarmup, nsamp_:nSamples] := Module[
  {run, probe, actualSamples, times = {}, lastResult = "N/A", i},
  Do[
    run = withTimeLimit[func];
    If[
      !TrueQ[run["ok"]],
      Return[<|"status" -> run["status"], "median_s" -> Null, "result" -> run["result"], "samples" -> 0|>]
    ],
    {nwarm}
  ];

  probe = withTimeLimit[func];
  If[
    !TrueQ[probe["ok"]],
    Return[<|"status" -> probe["status"], "median_s" -> Null, "result" -> probe["result"], "samples" -> 0|>]
  ];

  actualSamples = If[probe["time_s"] > slowThreshold, slowSamples, nsamp];
  AppendTo[times, probe["time_s"]];
  lastResult = probe["result"];

  For[i = 2, i <= actualSamples, i++,
    run = withTimeLimit[func];
    If[
      !TrueQ[run["ok"]],
      Return[<|
        "status" -> "timeout_partial",
        "median_s" -> Median[times],
        "result" -> lastResult,
        "samples" -> Length[times]
      |>]
    ];
    AppendTo[times, run["time_s"]];
    lastResult = run["result"];
  ];

  <|"status" -> "ok", "median_s" -> Median[times], "result" -> lastResult, "samples" -> Length[times]|>
];

scalarizeTraceResult[result_] := Module[{attempt, s},
  attempt = Quiet@Check[FullSimplify[result], result];
  s = StringTrim[ToString[attempt, InputForm]];
  If[
    s === "{}",
    <|"value" -> result, "scalarized" -> False|>,
    <|"value" -> attempt, "scalarized" -> True|>
  ]
];

results = <||>;

Options[runAndReport] = {traceRow -> False};

runAndReport[name_String, func_, OptionsPattern[]] := Module[
  {bench, status, res, n, ms, scalarized = True, traceMeta},
  Print["  Running: ", name, " ..."];
  bench = benchmark[func];
  status = bench["status"];
  res = bench["result"];
  n = bench["samples"];
  ms = If[NumberQ[bench["median_s"]], bench["median_s"] * 1000, Null];

  If[TrueQ[OptionValue[traceRow]],
    traceMeta = scalarizeTraceResult[res];
    res = traceMeta["value"];
    scalarized = traceMeta["scalarized"];
    If[status === "ok" && !scalarized, status = "non_scalar"];
    If[status === "timeout_partial" && !scalarized, status = "non_scalar_timeout_partial"];
  ];

  If[
    ms === Null,
    Print["    ", status, "  (N=", n, ", result: ", res, ")"],
    Print[
      "    ",
      NumberForm[ms, {8, 2}],
      " ms  (N=",
      n,
      ", status=",
      status,
      ", result: ",
      res,
      ")"
    ]
  ];

  results[name] = <|
    "median_ms" -> ms,
    "result" -> ToString[res, InputForm],
    "samples" -> n,
    "status" -> status,
    "scalarized" -> scalarized
  |>;
];

(* ============================================================================ *)
(* Element integrals: diagonal moments                                           *)
(* ============================================================================ *)

Print["\n=== Easy: Unitary |U_11|^{2k}, symbolic d ==="];
runAndReport["U_|U11|^2_sym", Function[{}, u11moment[1]]];
runAndReport["U_|U11|^4_sym", Function[{}, u11moment[2]]];
runAndReport["U_|U11|^6_sym", Function[{}, u11moment[3]]];

Print["\n=== Harder: Unitary |U_11|^{2k}, symbolic d ==="];
runAndReport["U_|U11|^8_sym", Function[{}, u11moment[4]]];

Print["\n=== Practical numeric coverage (no 10th-power rows) ==="];
runAndReport["U_|U11|^8_d=10", Function[{}, u11momentNumeric[4, 10]]];

(* ============================================================================ *)
(* Mixed element moments                                                         *)
(* ============================================================================ *)

Print["\n=== Mixed element moments, symbolic d ==="];
runAndReport[
  "U_|U11|^2|U12|^2_sym",
  Function[{}, uMomentSymbolic[mixedCaseDefinitions["U_|U11|^2|U12|^2"]["rows"], mixedCaseDefinitions["U_|U11|^2|U12|^2"]["cols"]]]
];
runAndReport[
  "U_|U11|^2|U12|^4_sym",
  Function[{}, uMomentSymbolic[mixedCaseDefinitions["U_|U11|^2|U12|^4"]["rows"], mixedCaseDefinitions["U_|U11|^2|U12|^4"]["cols"]]]
];
runAndReport[
  "U_|U11|^2|U22|^2_sym",
  Function[{}, uMomentSymbolic[mixedCaseDefinitions["U_|U11|^2|U22|^2"]["rows"], mixedCaseDefinitions["U_|U11|^2|U22|^2"]["cols"]]]
];

Print["\n=== Mixed element moments, numeric d ==="];
runAndReport[
  "U_|U11|^2|U12|^2_d=10",
  Function[{}, uMomentNumeric[mixedCaseDefinitions["U_|U11|^2|U12|^2"]["rows"], mixedCaseDefinitions["U_|U11|^2|U12|^2"]["cols"], 10]]
];
runAndReport[
  "U_|U11|^2|U12|^4_d=10",
  Function[{}, uMomentNumeric[mixedCaseDefinitions["U_|U11|^2|U12|^4"]["rows"], mixedCaseDefinitions["U_|U11|^2|U12|^4"]["cols"], 10]]
];
runAndReport[
  "U_|U11|^2|U22|^2_d=10",
  Function[{}, uMomentNumeric[mixedCaseDefinitions["U_|U11|^2|U22|^2"]["rows"], mixedCaseDefinitions["U_|U11|^2|U22|^2"]["cols"], 10]]
];

(* ============================================================================ *)
(* Trace moments: graph API                                                      *)
(* ============================================================================ *)

trMomentGraph[k_] := Module[{g},
  g = Join[
    Table[{{"U", i, 0, 1}, {"U", i, 1, 1}}, {i, k}],
    Table[{{"Ub", i, 0, 1}, {"Ub", i, 1, 1}}, {i, k}]
  ];
  integrateHaarUnitary[g, "U", {dd}, {dd}, dd]
];

Print["\n=== Trace moments: |tr(U)|^{2k}, symbolic d ==="];
runAndReport["U_|trU|^4_sym", Function[{}, trMomentGraph[2]], traceRow -> True];
runAndReport["U_|trU|^6_sym", Function[{}, trMomentGraph[3]], traceRow -> True];
runAndReport["U_|trU|^8_sym", Function[{}, trMomentGraph[4]], traceRow -> True];

(* ============================================================================ *)
(* Trace polynomials                                                             *)
(* ============================================================================ *)

Print["\n=== Trace polynomials: symbolic d ==="];
runAndReport[
  "U_trUAUdB_sym",
  Function[{},
    MultinomialexpectationvalueHaar[
      dd,
      {1, 2},
      {Subscript[X, 1], Subscript[X, 2]},
      usetrace -> True
    ]
  ]
];
runAndReport[
  "U_tr(UAUdB)^2_sym",
  Function[{},
    MultinomialexpectationvalueHaar[
      dd,
      {1, 2, 1, 2},
      {Subscript[X, 1], Subscript[X, 2], Subscript[X, 3], Subscript[X, 4]},
      usetrace -> True
    ]
  ]
];

results["_meta"] = <|
  "timestamp_utc" -> timestampUTC,
  "host" -> <|
    "hostname" -> $MachineName,
    "os" -> $OperatingSystem,
    "arch" -> $SystemID,
    "machine" -> $MachineName
  |>,
  "runtime" -> <|"name" -> "Mathematica", "version" -> $Version|>,
  "packages" -> <|"RTNI" -> "unknown"|>,
  "sources" -> <|
    "RTNI.wl" -> <|
      "path" -> rtniPath,
      "sha256" -> rtniHash
    |>
  |>,
  "script" -> scriptPath
|>;

(* ============================================================================ *)
(* Save results                                                                  *)
(* ============================================================================ *)

Print["\n", StringJoin@ConstantArray["=", 72]];
Print["Summary (median times in ms)"];
Print[StringJoin@ConstantArray["=", 72]];
Do[
  Print[
    StringPadRight[name, 30],
    "  ",
    If[
      NumberQ[results[name]["median_ms"]],
      NumberForm[results[name]["median_ms"], {8, 2}],
      results[name]["status"]
    ]
  ],
  {name, Keys[results]}
];

Export["results_rtni.json", Normal[results]];
Print["\nResults saved to results_rtni.json"];
