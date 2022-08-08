(* ::Package:: *)

(* ::Title:: *)
(*Package for obtaining the transformation matrix to canonical basis*)


(* ::Input::Initialization:: *)
BeginPackage["INITIAL`",{"FiniteFlow`"}];


(* ::Input::Initialization:: *)
(*Print["Package for obtaining the transformation matrix to canonical basis"]
Print["Author: Christoph Dlapa and Kai Yan"];
Print["Version: 1.2dev"];
Print["Last changes: 03.11.2020"];*)


(* ::Input:: *)
(*SetOptions[EvaluationNotebook[],ShowGroupOpener->True]*)


(* ::Subsection::Closed:: *)
(*Descriptions*)


(* ::Input::Initialization:: *)
functionNames={"psiStep","FFpsiStep","psiCalc","phiStep","FFphiStep","phiCalc","makePsiInvGraph","checkDegrees","eps0Calc","equStep","solCalc","FFInvMatMul","basisChange","denominators","BCalc","TCalc"};


(* ::Input::Initialization:: *)
Unprotect@@functionNames//Quiet;

Clear@@functionNames;

Remove@@functionNames//Quiet;


(* ::Input::Initialization:: *)
psiStep::usage="\
psiStep[A[1],v[1].A[i-1],x] calculates the first row of A[i], which is v[1].A[i].";


(* ::Input::Initialization:: *)
FFpsiStep::usage="\
FFpsiStep[A[1],v[1].A[i-1],x] calculates the first row of A[i], which is v[1].A[i]. For this calculation, FiniteFlow is used.";


(* ::Input::Initialization:: *)
psiCalc::usage="\
psiCalc[AMatrix,x] calculates the psi-matrix.";


(* ::Input::Initialization:: *)
phiStep::usage="\
phiStep[B[1],{v[1].B[i-1](j-1),v[1].B[i-1](j)},{TVariables[i-1](j-1),TVariables[i-1](j),sol,x] calculates v[1].B[i](j) in the variables \"TVariables\" while taking into account the solution \"sol\" frome the previous order. Output is a list of length two. The first entry gives the coefficients of the variables given in the second entry.";


(* ::Input::Initialization:: *)
FFphiStep::usage="\
FFphiStep[B[1],{v[1].B[i-1](j-1),v[1].B[i-1](j)},{TVariables[i-1](j-1),TVariables[i-1](j),sol,x] calculates v[1].B[i](j) in the variables \"TVariables\" while taking into account the solution \"sol\" frome the previous order. Output is a list of length two. The first entry gives the coefficients of the variables given in the second entry. For this calculation, FiniteFlow is used.";


(* ::Input::Initialization:: *)
phiCalc::usage="\
phiCalc[ansatzFunctions,phi(j-1),TVariables,sol,x,eps] calculates the next order eps^j of the phi-matrix from an ansatz based on \"ansatzFunctions\" and the previous order \"phi(j-1)\". A previous found solution is taken into account through the substitution rules in \"sol\".";


(* ::Input::Initialization:: *)
makePsiInvGraph::usage="\
makePsiInvGraph[psi,graphName] creates a graph named \"graphName\" that calculates the matrix-inverse of \"psi\". Relevant nodes for further calculations are \"input\", \"psinode\", \"psiInvnode\" and \"bsnode\".";


(* ::Input::Initialization:: *)
checkDegrees::usage="\
checkDegrees[psi,eps] computes the degrees in epsilon of the first row of Inverse[psi]. Output is the following list: {{dmin},{dmax,b1,...,bn}}, where dmin is the minimum power of epsilon of the GCD, dmax is the maximal degree of the GCD and b1 to bn are the maximal degrees of the bs.";


(* ::Input::Initialization:: *)
eps0Calc::usage="\
eps0Calc[psi,eps] computes the eps^0 part of the bs (not normalized)."


(* ::Input::Initialization:: *)
equStep::usage="\
equStep[psi,phi,epsOrder,TVariables,degrees,eps] computes the equation at order \"epsOrder\". \"TVariables\" are the Dot products used to describe the list \"phi\", and \"degrees\" is the output of checkDegrees.";


(* ::Input::Initialization:: *)
solCalc::usage="\
solCalc[psi,ansatzFunctions,degrees,x,eps] recursively tries to calculate a solution to the equation from an ansatz based on \"ansatzFunctions\". If no solution is found, a partial solution for the T-variables may be returned.";


(* ::Input::Initialization:: *)
FFInvMatMul::usage="\
FFInvMatMul[mat1,mat2,...] uses FiniteFlow to evaluate and reconstruct mat1.mat2...\n\
FFInvMatMul[mat1,mat2,...,\"InvertInput\"->n] inverts \"matn\" before the multiplication."


(* ::Input::Initialization:: *)
basisChange::usage="\
basisChange[AMatrix,TMatrix,x] uses FiniteFlow to compute the basis change given through \"TMatrix\" of the differential equation given through \"AMatrix\"."


(* ::Input::Initialization:: *)
denominators::usage="\
denominators[AMatrix] returns all unique denominator factors in \"AMatrix\".
denominators[AMatrix,x] returns all unique denominator factors in \"AMatrix\" that depend on \"x\".
denominators[AMatrix,x,eps] returns all unique denominator factors in \"AMatrix\" that depend on \"x\" and not on \"eps\"."


(* ::Input::Initialization:: *)
BCalc::usage="\
BCalc[sol,ansatzFunctions,n,x,eps] computes the BMatrix of size \"n x n\" from the \"ansatzFunctions\" and the solution \"sol\" given by solCalc."


(* ::Input::Initialization:: *)
TCalc::usage="\
TCalc[AMatrix,ansatzFunctions,x,eps] uses FiniteFlow to compute the basis change of the differential equation given through \"AMatrix\" to a canonical basis. The ansatz for the canonical form is given through the functions in the List \"ansatzFunctions\", e.g. ansatzFunctions=D[Log[#],x]&/@letters."


(* ::Input::Initialization:: *)
nthO::badmatrix1 = "The first argument must be a square matrix.";
nthO::badmatrix2= "The second argument must be a square matrix.";
nthO::badlist="The first argument must be a list of length n.";
nthO::badvars="Not enough variables.";
nthO::badtestnums="\"TestNumbers\" should be a list of 64-bit integers with length Length[Variables]-1.";
nthO::singmatrix="The matrix is singular, suggesting that the derivatives do not couple to all master integrals.";


(* ::Chapter:: *)
(*1 Functions*)


(* ::Input::Initialization:: *)
If[!ValueQ[Global`$UsedSymbols],Global`$UsedSymbols={Global`T,Global`v,Global`m}];


(* ::Input::Initialization:: *)
Begin["`Private`"]


(* ::Input::Initialization:: *)
{T,v,m}=Global`$UsedSymbols;
Set@@@Transpose[{Global`$UsedSymbols,Global`$UsedSymbols}];


(* ::Section:: *)
(*1.1 Computing psi*)


(* ::Subsection::Closed:: *)
(*1.1.1 psiStep*)


(* ::Input::Initialization:: *)
Options[psiStep]:={"Variables"->Automatic,"DerivativeRules"->{},Method->Automatic};
psiStep[startA_,stepA_,xv_,OptionsPattern[]]:=Module[{opt,sz,graphvars,outA,Dreps},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[psiStep];
(*testing some things*)

If[!SquareMatrixQ[startA],Message[nthO::badmatrix1];Return[$Failed]];

sz=Length[startA];

If[Dimensions[stepA]=!={sz},Message[nthO::badlist];Return[$Failed]];
Dreps=Flatten[{OptionValue["DerivativeRules"]}];

graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[startA]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*xv=graphvars[[2]];*)

outA=D[stepA,xv]/.Dreps;
outA=outA+stepA . startA;

If[$KernelCount>0,
outA=ParallelMap[Together,outA,Sequence@@FilterRules[opt,Options[ParallelMap]]];
,
outA=Together[outA];
];

outA];


(* ::Subsection::Closed:: *)
(*1.1.2 FFpsiStep*)


(* ::Input::Initialization:: *)
Options[FFpsiStep]:=Join[{"Variables"->Automatic,"MaxDegree"->1000,"MaxPrimes"->150,"DerivativeRules"->{}},DeleteCases[Options[FFReconstructFunction],("MaxDegree"->_)|("MaxPrimes"->_)]];
FFpsiStep[startAIn_,stepAIn_,xv_,OptionsPattern[]]:=Module[{outA,pA,qA,psiGraph,input,starANode,stepANode,derNode1,derNode2,graphvars,mulNode,addNode,sz,startA,stepA,
der1,der2,time,derNode21,derNode22,derNode23,opt,Dreps},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFpsiStep];
(*testing some things*)

If[!SquareMatrixQ[startAIn],Message[nthO::badmatrix1];Return[$Failed]];

sz=Length[startAIn];

If[Dimensions[stepAIn]=!={sz},Message[nthO::badlist];Return[$Failed]];

Dreps=Flatten[{OptionValue["DerivativeRules"]}];

graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[startAIn]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*xv=graphvars[[2]];*)
graphvars=Union[graphvars,Dreps[[All,1]]];

startA=Flatten[startAIn];
stepA=Flatten[stepAIn];

(*Computing the derivative*)

{pA,qA}=Transpose[NumeratorDenominator[stepA]];
time=SessionTime[];
der1=D[pA,xv]/qA;
(*der2=-((pA D[qA,xv])/qA^2);*)
(*Second derivative not computed explicitely because faster*);

FFDeleteGraph[psiGraph];
FFNewGraph[psiGraph,input,graphvars];
FFAlgRatFunEval[psiGraph,starANode,{input},graphvars,startA];
FFAlgRatFunEval[psiGraph,stepANode,{input},graphvars,stepA];
FFAlgRatFunEval[psiGraph,derNode1,{input},graphvars,der1];
(*time=SessionTime[];*)
FFAlgRatFunEval[psiGraph,derNode21,{input},graphvars,1/qA];
FFAlgRatFunEval[psiGraph,derNode22,{input},graphvars,-pA];
FFAlgRatFunEval[psiGraph,derNode23,{input},graphvars,D[qA,xv]];
FFAlgMul[psiGraph,derNode2,{derNode21,derNode21,derNode22,derNode23}];

(*Now do the matrix multiplication and add with the derivative*)
FFAlgMatMul[psiGraph,mulNode,{stepANode,starANode},1,sz,sz];
FFAlgAdd[psiGraph,addNode,{derNode1,derNode2,mulNode}];
FFGraphOutput[psiGraph,addNode];

(*Print["Reconstructing: "];*)
outA=FFReconstructFunction[psiGraph,graphvars,
Sequence@@FilterRules[opt,Options[FFReconstructFunction]]];

FFDeleteGraph[psiGraph];
If[outA===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];Return[$Failed]];
If[outA===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];Return[$Failed]];

outA/.Dreps];


(* ::Subsection::Closed:: *)
(*1.1.3 psiCalc*)


(* ::Input::Initialization:: *)
Options[psiCalc]:=Join[{"Silent"->False,"FiniteFlow"->False,"MaxDerivatives"->Automatic,"TakeRow"->1},Options[FFpsiStep]];
psiCalc[Am_,xv_,OptionsPattern[]]:=Module[{Aco,sz,psi1,time,graphvars,i,opt,printQ,FFQ,nthreads,tr},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[psiCalc];

tr=OptionValue["TakeRow"];
If[Head[tr]===List,
psi1=Join@@Table[
psiCalc[Am,xv,Sequence@@DeleteCases[opt,("TakeRow"->_)|("MaxDerivatives"->_)],"TakeRow"->tr[[i,1]],"MaxDerivatives"->tr[[i,2]]]
,{i,Length[tr]}];
Return[psi1];
];

(*Checking things*)
printQ=!OptionValue["Silent"];
If[!SquareMatrixQ[Am],Message[nthO::badmatrix1];Return[$Failed]];

sz=OptionValue["MaxDerivatives"];
If[sz===Automatic,
sz=Length[Am]];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[Am]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*xv=graphvars[[2]];*)

nthreads=OptionValue["NThreads"];
If[nthreads===Automatic,nthreads=FFAutomaticNThreads[]];
FFQ=OptionValue["FiniteFlow"];
If[FFQ===Automatic,
If[((sz>20)&&(nthreads>20))||(Length[graphvars]>2),
FFQ=True,FFQ=False];
];

(*Computing psi step-by-step*)
time=SessionTime[];
Aco[1]=Am[[tr]]//Together;
Catch[
If[printQ,
Monitor[
Do[
If[FFQ,
Aco[i]=FFpsiStep[Am,Aco[i-1],xv,Sequence@@FilterRules[opt,Options[FFpsiStep]]],
Aco[i]=psiStep[Am,Aco[i-1],xv,Sequence@@FilterRules[opt,Options[psiStep]]];
];
If[Aco[i]===$Failed,Throw[Return[$Failed]]];
,{i,2,sz}];
,{i,sz}];
,
Do[
If[FFQ,
Aco[i]=FFpsiStep[Am,Aco[i-1],xv,Sequence@@FilterRules[opt,Options[FFpsiStep]]],
Aco[i]=psiStep[Am,Aco[i-1],xv,Sequence@@FilterRules[opt,Options[psiStep]]];
];
If[Aco[i]===$Failed,Throw[Return[$Failed]]];
,{i,2,sz}];
];
];
(*Assembling the result*)
psi1=Table[Aco[i],{i,sz}];
If[printQ,Print[ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];

psi1];


(* ::Section:: *)
(*1.2 Computing phi*)


(* ::Subsection::Closed:: *)
(*1.2.1 phiStep*)


(* ::Input::Initialization:: *)
Options[phiStep]:=Join[{"Last"->False},Options[psiStep]];
phiStep[startB_,stepBIn_,varsIn_,sol1_,xv_,OptionsPattern[]]:=Module[{opt,graphvars,stepBmul,stepBder,varsTd,lsz,stepout,vars,oldPhi,stepBmulM,zpos,Dreps},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[phiStep];

(*testing input*)
Dreps=Flatten[{OptionValue["DerivativeRules"]}];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[startB]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*xv=graphvars[[2]];*)
vars=varsIn;

{stepBmul,stepBder}=DeleteCases[#,0,{1}]&/@stepBIn;

lsz=Length[startB];

If[!OptionValue["Last"],

stepBder=(D[stepBder,xv]/.Dreps) . vars[[2]];

varsTd[1]=Table[Append[#,m[l]]&/@vars[[1]],{l,lsz}]//Flatten;
varsTd[1]=varsTd[1]//.sol1;
stepBmulM=Table[stepBmul startB[[l]],{l,lsz}]//Flatten;

stepout=stepBder+stepBmulM . varsTd[1];
varsTd[1]=Cases[varsTd[1],T[__],\[Infinity]]//Union;
varsTd[1]=Join[varsTd[1],vars[[2]]]//Union;
stepout=Coefficient[stepout,varsTd[1]];
,
varsTd[1]=stepout={};
];

(*Now also apply sol to the previous order in epsilon*)
varsTd[0]=vars[[1]]//.sol1;
oldPhi=stepBmul . varsTd[0];
varsTd[0]=Cases[varsTd[0],T[__],\[Infinity]]//Union;
oldPhi=Coefficient[oldPhi,varsTd[0]];

If[$KernelCount>0,
stepout=ParallelMap[Together,Join[oldPhi,stepout],Sequence@@FilterRules[opt,Options[ParallelMap]]];
oldPhi=stepout[[1;;Length[oldPhi]]];
stepout=stepout[[(Length[oldPhi]+1);;-1]];
,
{oldPhi,stepout}=Together[{oldPhi,stepout}];
];

zpos=Position[oldPhi,0,{1}];
oldPhi=Delete[oldPhi,zpos];
varsTd[0]=Delete[varsTd[0],zpos];
zpos=Position[stepout,0,{1}];
stepout=Delete[stepout,zpos];
varsTd[1]=Delete[varsTd[1],zpos];

{{oldPhi,varsTd[0]},{stepout,varsTd[1]}}];


(* ::Subsection::Closed:: *)
(*1.2.2 FFphiStep*)


(* ::Input::Initialization:: *)
Options[FFphiStep]:=Join[{"Last"->False},Options[FFpsiStep]];
FFphiStep[startB_,stepBIn_,varsIn_,sol1_,xv_,OptionsPattern[]]:=Module[{graphvars,pB,qB,phiGraph,input,derNode1,derNode2,mul,solRedLs,solRedLsVars,
addpattrn,inputnodes,takepattern,varsTd,deletepos,
nummul,tp,stepout,num,take,muloutNode,stepBoutNode,lsz,opt,nzlearn,nznode,
stepBder,stepBmul,phim1Node,chainNode,vars,oldPhi,varsTdOld,Dreps},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFphiStep];

If[OptionValue["Last"],Return[FFphiStepLast[startB,stepBIn,varsIn,sol1,Sequence@@opt]]];

(*testing input*)
Dreps=Flatten[{OptionValue["DerivativeRules"]}];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[startB]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*xv=graphvars[[2]];*)
graphvars=Union[graphvars,Dreps[[All,1]]];

{stepBmul,stepBder}=stepBIn;

{pB,qB}=Transpose[NumeratorDenominator[stepBder]];
(*sol1=Join[sol1In,{T[a___]:>0/;(Length[{a}]>(epsor+1))}];*)
varsTd[0]=varsIn[[1]];
lsz=Length[startB];

FFDeleteGraph[phiGraph];
FFNewGraph[phiGraph,input,graphvars];

(*build input nodes*)
Do[FFAlgRatFunEval[phiGraph,m[i],{input},graphvars,{startB[[i]]}],{i,Length[startB]}];
(*Now we also apply the solution of previous order to the phi of previous order*)
FFAlgRatFunEval[phiGraph,phim1Node,{input},graphvars,stepBmul];
varsTd[0]=varsTdOld[0]=makeSolRedNode[phiGraph,mul[lsz+1],{phim1Node},{varsTd[0]},sol1];

(*Then take the derivative*)
FFAlgRatFunEval[phiGraph,derNode1,{input},graphvars,D[pB,xv]/qB];
FFAlgRatFunEval[phiGraph,derNode2,{input},graphvars,-((pB D[qB,xv])/qB^2)];

(*Then do the multiplication*)
If[varsTd[0]=!={},
(*Do the multiplication with auxiliary nodes*)
Do[FFAlgMatMul[phiGraph,mul[l],{m[l],mul[lsz+1]},1,1,Length[varsTd[0]]],{l,lsz}];
varsTd[1]=Table[Append[#,m[l]]&/@varsTd[0],{l,lsz}];

(*put red. here*)
varsTd[1]=makeSolRedNode[phiGraph,muloutNode,Array[mul,lsz],varsTd[1],sol1];
,
varsTd[1]={}];

(*Then add everything*)
varsTd[0]=varsIn[[2]];
inputnodes={derNode1,derNode2,muloutNode};
addpattrn=varsTd/@{0,0,1};
deletepos={#}&/@(Position[addpattrn,{},{1}][[All,1]]);
inputnodes=Delete[inputnodes,deletepos];
addpattrn=Delete[addpattrn,deletepos];

If[inputnodes=!={},
varsTd[1]=Union[Flatten[addpattrn]];
takepattern=Position[addpattrn,#,\[Infinity]][[All,1;;2]]&/@varsTd[1];
takepattern=Apply[tp,takepattern,{2}];
takepattern=Table[tp[i,j],{i,Length[addpattrn]},{j,Length[addpattrn[[i]]]}]->takepattern;
FFAlgTakeAndAdd[phiGraph,stepBoutNode,inputnodes,takepattern];
,
Return[{{0},{}}]];

FFAlgChain[phiGraph,chainNode,{mul[lsz+1],stepBoutNode}];

FFGraphOutput[phiGraph,chainNode];
FFAlgNonZeroes[phiGraph,nznode,{chainNode}];
FFGraphOutput[phiGraph,nznode];
nzlearn=FFNonZeroesLearn[phiGraph];

stepout=FFReconstructFunction[phiGraph,graphvars,Sequence@@FilterRules[opt,Options[FFReconstructFunction]]];
FFDeleteGraph[phiGraph];

If[stepout===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];Return[{{$Failed,$Failed},{$Failed,$Failed}}]];
If[stepout===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];Return[{{$Failed,$Failed},{$Failed,$Failed}}]];

varsTd[0]=varsTdOld[0];
vars=Join[varsTd[0],varsTd[1]];
varsTd[1]=vars[[Select["NonZero"/.nzlearn,#>Length[varsTd[0]]&]]];
varsTd[0]=vars[[Select["NonZero"/.nzlearn,#<=Length[varsTd[0]]&]]];

oldPhi=stepout[[1;;Length[varsTd[0]]]]/.Dreps;
stepout=stepout[[(Length[varsTd[0]]+1);;-1]]/.Dreps;

{{oldPhi,varsTd[0]},{stepout,varsTd[1]}}];


(* ::Input::Initialization:: *)
Options[FFphiStepLast]:=Options[FFphiStep];
FFphiStepLast[startB_,stepBIn_,varsIn_,sol1_,OptionsPattern[]]:=Module[{opt,graphvars,stepBmul,varsTd,lsz,phiGraph,input,phim1Node,mul,nznode,nzlearn,stepout,oldPhi},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFphiStepLast];

(*testing input*)
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[startB]];

stepBmul=stepBIn[[1]];

varsTd[0]=varsIn[[1]];
lsz=Length[startB];

FFDeleteGraph[phiGraph];
FFNewGraph[phiGraph,input,graphvars];

(*build input nodes*)
Do[FFAlgRatFunEval[phiGraph,m[i],{input},graphvars,{startB[[i]]}],{i,Length[startB]}];
(*Now we also apply the solution of previous order to the phi of previous order*)
FFAlgRatFunEval[phiGraph,phim1Node,{input},graphvars,stepBmul];
varsTd[0]=makeSolRedNode[phiGraph,mul[lsz+1],{phim1Node},{varsTd[0]},sol1];

FFGraphOutput[phiGraph,mul[lsz+1]];
FFAlgNonZeroes[phiGraph,nznode,{mul[lsz+1]}];
FFGraphOutput[phiGraph,nznode];
nzlearn=FFNonZeroesLearn[phiGraph];

stepout=FFReconstructFunction[phiGraph,graphvars,Sequence@@FilterRules[opt,Options[FFReconstructFunction]]];
FFDeleteGraph[phiGraph];

If[stepout===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];Return[{{$Failed,$Failed},{$Failed,$Failed}}]];
If[stepout===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];Return[{{$Failed,$Failed},{$Failed,$Failed}}]];

varsTd[0]=varsTd[0][["NonZero"/.nzlearn]];
oldPhi=stepout;

{{oldPhi,varsTd[0]},{{},{}}}];


(* ::Input::Initialization:: *)
Options[makeSolRedNode]:={};
makeSolRedNode[phiGraph_,muloutNode_,innodes_,Tvars_,sol1_]:=Module[{varsTd,solRedLs,solRedLsVars,takepattern,take,lsz,mul,num,nummul,tp,addpattrn,inputnodes,opt},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[makeSolRedNode];
varsTd[1]=Tvars;
lsz=Length[innodes];
Do[mul[i]=innodes[[i]],{i,lsz}];

(*Check the reduction of the Td variables*)
solRedLs=varsTd[1]//.sol1//Expand;
solRedLs=Map[If[Head[#]===Plus,List@@#,{#}]&,solRedLs,{2}];
solRedLsVars=Cases[solRedLs,T[___],\[Infinity]]//Union;

(*Take individual elements*)
If[solRedLsVars=!={},

Do[
takepattern={varsTd[1][[i]]}->{varsTd[1][[i,j]]};
FFAlgTake[phiGraph,take[i,j],{mul[i]},takepattern];
,{i,lsz},{j,Length[varsTd[1][[i]]]}];

(*Do all Times in the rules with auxiliary nodes*)
Do[
If[Head[solRedLs[[i,j,k]]]===Times,
FFAlgRatNumEval[phiGraph,num[i,j,k],{solRedLs[[i,j,k,1]]}];
FFAlgMatMul[phiGraph,nummul[i,j,k],{num[i,j,k],take[i,j]},
1,1,1];
,
FFAlgTake[phiGraph,nummul[i,j,k],{take[i,j]},{{tp[1]}}->{tp[1]}];
]
,{i,lsz},{j,Length[solRedLs[[i]]]},{k,Length[solRedLs[[i,j]]]}];

(*Assemble the result*)
addpattrn=Position[solRedLs,#,\[Infinity]][[All,1;;3]]&/@solRedLsVars;
inputnodes=Table[nummul[i,j,k],{i,Length[solRedLs]},{j,Length[solRedLs[[i]]]},{k,Length[solRedLs[[i,j]]]}]//Flatten;
takepattern=Table[tp[i,j,k],{i,Length[solRedLs]},{j,Length[solRedLs[[i]]]},{k,Length[solRedLs[[i,j]]]}];
takepattern={#}&/@Flatten[takepattern];
takepattern=takepattern->Apply[tp,addpattrn,{2}];
FFAlgTakeAndAdd[phiGraph,muloutNode,inputnodes,takepattern];

varsTd[1]=solRedLsVars;
,
varsTd[1]={};

FFAlgRatNumEval[phiGraph,muloutNode,{}]];

varsTd[1]];


(* ::Subsection::Closed:: *)
(*1.2.3 phiCalc*)


(* ::Input::Initialization:: *)
Options[phiCalc]:=Join[Options[psiCalc],{"TakeRow"->1}];
phiCalc[ansatzF_,phim1_,varsIn_,sol1_,xv_,eps_,OptionsPattern[]]:=Module[{phi1,varsTd,time,startB,Bco,graphvars,allVars,i,opt,print,epsor,sz,vars,nthreads,FFQ,oldPhi,oldVars,tr},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[phiCalc];

time=SessionTime[];
graphvars=OptionValue["Variables"];
(*ansatzF=OptionValue["AnsatzFunctions"];
If[ansatzF===Automatic,
ansatzF=letters];*)
If[graphvars===Automatic,graphvars=Variables[ansatzF]];
graphvars=Join[{eps,xv},DeleteCases[graphvars,eps|xv]];
If[Length[graphvars]<2,Message[nthO::badvars];Return[$Failed]];
(*eps=graphvars[[1]];
xv=graphvars[[2]];*)
print=!OptionValue["Silent"];
(*If[OptionValue["AnsatzFunctions"]===Automatic,
ansatzF=D[Log[#],xv]&/@ansatzF];*)

epsor=Dimensions[phim1][[1;;2]];

sz=OptionValue["MaxDerivatives"];
If[sz===Automatic,
sz=epsor[[1]]];
epsor=epsor[[2]];

tr=OptionValue["TakeRow"];

nthreads=OptionValue["NThreads"];
If[nthreads===Automatic,nthreads=FFAutomaticNThreads[]];
FFQ=OptionValue["FiniteFlow"];
If[FFQ===Automatic,
If[((sz>20)&&(nthreads>20))||(Length[graphvars]>2),
FFQ=True,FFQ=False];
];

If[epsor=!=0,
phi1=phim1;
vars=varsIn;
,
epsor=1;
sz=sz+1;
phi1=Transpose[{{IdentityMatrix[sz][[1]]}},{3,2,1}];
vars=Table[{{}},sz];
vars[[1]]={{T[v[tr]]}};
];

Do[Bco[i]={0};varsTd[i]={},{i,epsor}];

startB=Table[eps Together[ansatzF[[k]]],{k,Length[ansatzF]}];
Catch[
If[print,
Monitor[
Do[
If[FFQ,
{{oldPhi,oldVars},{Bco[i],varsTd[i]}}=FFphiStep[startB,{phi1[[i-1,epsor]],Bco[i-1]},{vars[[i-1,epsor]],varsTd[i-1]},sol1,xv,
Sequence@@FilterRules[Join[opt,{"Last"->(i===(sz+1))}],Options[FFphiStep]]];
,
{{oldPhi,oldVars},{Bco[i],varsTd[i]}}=phiStep[startB,{phi1[[i-1,epsor]],Bco[i-1]},{vars[[i-1,epsor]],varsTd[i-1]},sol1,xv,
Sequence@@FilterRules[Join[opt,{"Last"->(i===(sz+1))}],Options[phiStep]]];
];
If[Bco[i]===$Failed,Throw[Return[$Failed]]];
phi1[[i-1,epsor]]=oldPhi;
vars[[i-1,epsor]]=oldVars;
,{i,epsor+1,sz+1}];
,{i,sz}];
,
Do[
If[FFQ,
{{oldPhi,oldVars},{Bco[i],varsTd[i]}}=FFphiStep[startB,{phi1[[i-1,epsor]],Bco[i-1]},{vars[[i-1,epsor]],varsTd[i-1]},sol1,xv,
Sequence@@FilterRules[Join[opt,{"Last"->(i===(sz+1))}],Options[FFphiStep]]];
,
{{oldPhi,oldVars},{Bco[i],varsTd[i]}}=phiStep[startB,{phi1[[i-1,epsor]],Bco[i-1]},{vars[[i-1,epsor]],varsTd[i-1]},sol1,xv,
Sequence@@FilterRules[Join[opt,{"Last"->(i===(sz+1))}],Options[phiStep]]];
];
If[Bco[i]===$Failed,Throw[Return[$Failed]]];
phi1[[i-1,epsor]]=oldPhi;
vars[[i-1,epsor]]=oldVars;
,{i,epsor+1,sz+1}];
];
];

If[Dimensions[phim1][[2]]===0,
phi1=Table[{Bco[i]},{i,2,sz}];
allVars=Table[{varsTd[i]},{i,2,sz}];
,
phi1=Table[Join[phi1[[i]],{Bco[i]}],{i,sz}];
allVars=Table[Join[vars[[i]],{varsTd[i]}],{i,sz}];
];
If[print,Print[ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];
{phi1,allVars}];


(* ::Section:: *)
(*1.3 psiInv*)


(* ::Subsection::Closed:: *)
(*1.3.1 makePsiInvGraph*)


(* ::Input::Initialization:: *)
Options[makePsiInvGraph]={"Variables"->Automatic,"TakeRow"->1};
makePsiInvGraph[psi1_,psiInvGraph_,OptionsPattern[]]:=Module[{graphvars,t1,xt,yt,m,t,id,sz,solverlearn,singularTest,tr},
If[!SquareMatrixQ[psi1],Message[nthO::badmatrix1];Return[$Failed]];
sz=Length[psi1];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[psi1]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
tr=OptionValue["TakeRow"];
(*Print["Taking row number ",tr];*)

FFDeleteGraph[psiInvGraph];
FFNewGraph[psiInvGraph,input,graphvars];

(*Build equation system for inverting psi*)
FFAlgRatFunEval[psiInvGraph,psinode,{input},graphvars,Flatten[psi1]];
FFAlgRatNumEval[psiInvGraph,id,Flatten[Append[#,0]&/@(-IdentityMatrix[sz])]];
FFAlgTake[psiInvGraph,t1,{psinode,id},{Flatten[Table[m[i,j],{i,sz},{j,sz}]],Flatten[Table[t[i,j],{i,sz},{j,sz+1}]]}->Flatten[Table[Join[Table[m[i,j],{j,sz}],Table[t[i,j],{j,sz+1}]],{i,sz}]]];
(*Create and learn solver*)
FFAlgNodeDenseSolver[psiInvGraph,psiInvnode,{t1},sz,
Join[Array[xt,sz],Array[yt,sz]]];
FFSolverOnlyHomogeneous[psiInvGraph,psiInvnode];
FFGraphOutput[psiInvGraph,psiInvnode];
solverlearn = FFDenseSolverLearn[psiInvGraph,Join[Array[xt,sz],Array[yt,sz]]];
singularTest={"DepVars","IndepVars","ZeroVars"}/.solverlearn;
If[singularTest=!={Array[xt,sz],Array[yt,sz],{}},Return[FFImpossible]];

(*We only need the first row*)
FFAlgTake[psiInvGraph,bsnode,{psiInvnode},{Flatten[Table[t[i,j],{i,sz},{j,sz}]]}->Table[t[tr,j],{j,sz}]];
FFGraphOutput[psiInvGraph,bsnode];

solverlearn
];


(* ::Subsection::Closed:: *)
(*1.3.2 Check Degrees of Denominator and bs*)


(* ::Input::Initialization:: *)
Options[checkDegrees]:=Join[Options[makePsiInvGraph],{"TestNumbers"->Automatic,"MaxDegree"->1000,"MaxPrimes"->150},DeleteCases[Options[FFReconstructFunctionMod],("MaxDegree"->_)|("MaxPrimes"->_)]];
checkDegrees[psi1_,eps_,OptionsPattern[]]:=Module[{t,detNode2,detNode3,denomIn,denomgraph,bnode,epsnode,subnode,denom,out,bs,
sz,graphvars,psiInvGraph,bsnode,comp,opt,testnums,learn,tr},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[checkDegrees];
(*testing some things*)
If[!SquareMatrixQ[psi1],Message[nthO::badmatrix1];Return[$Failed]];

graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[psi1]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
graphvars=Join[{eps},DeleteCases[graphvars,eps]];
(*eps=graphvars[[1]];*)

tr=OptionValue["TakeRow"];

testnums=OptionValue["TestNumbers"];
If[testnums===Automatic,testnums=RandomInteger[{10,10000},Length[graphvars]-1]];
If[(Length[testnums]=!=(Length[graphvars]-1))||(!AllTrue[testnums,Head[#]===Integer&]),Message[nthO::badtestnums];Return[$Failed]];
sz=Length[psi1];

FFDeleteGraph[psiInvGraph];

(*calculate bs and denom*)
learn=makePsiInvGraph[psi1,psiInvGraph,"Variables"->graphvars,"TakeRow"->tr];
If[learn===FFImpossible,Message[nthO::singmatrix];Return[$Failed]];
(*build super-graph to check the degrees of the sub-graph for specific values of x*)
(*build graph for denominator*)
FFNewGraph[denomgraph,denomIn,{eps}];
FFAlgRatNumEval[denomgraph,epsnode,testnums];
FFAlgChain[denomgraph,comp,{denomIn,epsnode}];
FFAlgSimpleSubgraph[denomgraph,subnode,{comp},psiInvGraph];
FFAlgTakeAndAdd[denomgraph,detNode2,{subnode},{Table[t[j],{j,sz}]}->{Table[t[j],{j,sz}]}];
FFGraphOutput[denomgraph,detNode2];

denom=FFReconstructFunctionMod[denomgraph,{eps},Sequence@@FilterRules[opt,Options[FFReconstructFunctionMod]]];
If[denom===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];
FFDeleteGraph[denomgraph];Return[$Failed]];
If[denom===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];
FFDeleteGraph[denomgraph];Return[$Failed]];
denom=Denominator[denom][[1]];
out={Coefficient[denom,eps,0],Min[Cases[Expand[{denom} eps^2],eps^a_:>a,\[Infinity]]]-2,Exponent[denom,eps]};

(*multiply bs with denominator*)
FFAlgRatFunEval[denomgraph,detNode3,{denomIn},{eps},{denom}];
FFAlgMatMul[denomgraph,bnode,{detNode3,subnode},1,1,sz];
FFGraphOutput[denomgraph,bnode];
(*no explicit reconstruction is needed, just the degrees*)
bs=FFAllDegrees[denomgraph,Sequence@@FilterRules[opt,Options[FFAllDegrees]]][[All,1]];
If[bs===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];FFDeleteGraph[denomgraph];Return[$Failed]];
If[bs===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];FFDeleteGraph[denomgraph];Return[$Failed]];

FFDeleteGraph[denomgraph];
FFDeleteGraph[psiInvGraph];

(*{out[[1;;2]],Join[{out[[3]]},bs]}*)
{out[[2]],Join[{out[[3]]},bs]}
];


(* ::Subsection::Closed:: *)
(*1.3.3 Compute eps^0 of bs analytically*)


(* ::Input::Initialization:: *)
Options[eps0Calc]:=Join[Options[FFpsiStep],Options[makePsiInvGraph]]//Union;
eps0Calc[psi1_,eps_,OptionsPattern[]]:=Module[{opt,graphvars,expandvars,sz,LaurentGraph,psiInvGraph,learn,denomnode,t,laurent,laurentlearn,stepout,tr},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[eps0Calc];
(*testing some things*)
If[!SquareMatrixQ[psi1],Message[nthO::badmatrix1];Return[$Failed]];

graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[psi1]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
(*eps=graphvars[[1]];*)
(*expandvars=graphvars[[2;;-1]];*)
graphvars=Join[{eps},DeleteCases[graphvars,eps]];
expandvars=DeleteCases[graphvars,eps];

tr=OptionValue["TakeRow"];

sz=Length[psi1];

FFDeleteGraph[LaurentGraph];
FFDeleteGraph[psiInvGraph];

(*calculate bs and denom*)
learn=makePsiInvGraph[psi1,psiInvGraph,"Variables"->graphvars,"TakeRow"->tr];
If[learn===FFImpossible,Message[nthO::singmatrix];Return[$Failed]];
(*add the bs*)
(*FFAlgTakeAndAdd[psiInvGraph,denomnode,{bsnode},{Table[t[j],{j,sz}]}->{Table[t[j],{j,sz}]}];*)
(*now compute the eps^0 part*)
FFGraphOutput[psiInvGraph,bsnode];

FFNewGraph[LaurentGraph,input,expandvars];
FFAlgLaurent[LaurentGraph,laurent,{input},psiInvGraph,0,Sequence@@FilterRules[opt,Options[FFAlgLaurent]]];
FFGraphOutput[LaurentGraph,laurent];

laurentlearn=FFLaurentLearn[LaurentGraph];
If[laurentlearn===$Failed,Print["FFLaurentLearn failed. Try increasing MaxDegree."];FFDeleteGraph[psiInvGraph];FFDeleteGraph[LaurentGraph];Return[$Failed]];

stepout=FFReconstructFunction[LaurentGraph,expandvars,Sequence@@FilterRules[opt,Options[FFReconstructFunction]]];

FFDeleteGraph[LaurentGraph];
FFDeleteGraph[psiInvGraph];

If[stepout===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];Return[{$Failed,$Failed}]];
If[stepout===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];Return[{$Failed,$Failed}]];

stepout=FFLaurentSol[stepout,eps,laurentlearn];

stepout
];


(* ::Section:: *)
(*1.4 The equation*)


(* ::Subsection::Closed:: *)
(*1.4.1 Calculate the Equation*)


(* ::Input::Initialization:: *)
Options[equStep]=Join[{"Silent"->False,"MaxDegree"->1000,"MaxPrimes"->150},Options[makePsiInvGraph],DeleteCases[Options[FFGraphEvaluateMany],"PrimeNo"->_]];
equStep[psi1_,phi1In_,epsor_,varsTdIn_,degsIn_,eps_,OptionsPattern[]]:=Module[{graphvars,sz,varsTd,LaurentGraph,equGraph,
bst,phit,am,t,tp,equNode,expandvars,laurent,laurentlearn,reconstructed,equs,normNode,detNode1,amNorm,
inputnodes,addpattrn,deletepos,takepattern,minusnormNode,
testpoints,testpointsnew,matrix,rank,opt,learn,nthreadsm,nthreads,print,time,prime,
tpnode,zeronode,SolveGraph,evalNode,nonzerovarsTd,takenode,solveNode,z,i,solverlearn,nulltest,primetest,primetestpoint,matrixrec,matrixnew,deg,
leadingOrder,phi1,varsTdFlat,tr},
time=SessionTime[];
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[equStep];
(*test input*)
If[!SquareMatrixQ[psi1],Message[nthO::badmatrix1];Return[$Failed]];

deg=degsIn[[1]];

print=!OptionValue["Silent"];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[{psi1,phi1In}]];
(*If[Length[graphvars]<1,Print["here"];Message[nthO::badvars];Return[$Failed]];*)
(*eps=graphvars[[1]];*)
graphvars=Join[{eps},DeleteCases[graphvars,eps]];
expandvars=DeleteCases[graphvars,eps];

tr=OptionValue["TakeRow"];

phi1=DeleteCases[Flatten[#],0,{1}]&/@phi1In;
varsTdFlat=Flatten/@varsTdIn;

sz=Length[psi1];
varsTd[0]={T[v[tr]]};
Do[varsTd[i]=varsTdFlat[[i]],{i,1,sz}];

FFDeleteGraph[equGraph];
FFDeleteGraph[LaurentGraph];

(*calculate bs and denom*)
learn=makePsiInvGraph[psi1,equGraph,Sequence@@DeleteCases[FilterRules[opt,Options[makePsiInvGraph]],("Variables"->_)],"Variables"->graphvars];
If[learn===FFImpossible,Message[nthO::singmatrix];Return[$Failed]];
FFAlgRatFunEval[equGraph,normNode,{input},graphvars,{eps^deg}];

FFAlgRatFunEval[equGraph,minusnormNode,{input},graphvars,{-eps^deg}];

(*multiply phi with psi and don't forget to normalize in eps*)
Do[
FFAlgTake[equGraph,bst[i],{bsnode},{Table[t[j],{j,sz}]}->{t[i]}];
FFAlgRatFunEval[equGraph,phit[i],{input},graphvars,phi1[[i]]],
{i,sz}];

Do[
FFAlgMul[equGraph,amNorm[i],{normNode,bst[i]}];
FFAlgMatMul[equGraph,am[i],{amNorm[i],phit[i]},1,1,Length[varsTd[i]]];
,
{i,sz}];

(*add with the derivative*)

inputnodes=Join[{minusnormNode},Table[am[i],{i,sz}]];
addpattrn=Table[varsTd[i],{i,0,sz}];
deletepos={#}&/@(Position[addpattrn,{},{1}][[All,1]]);
inputnodes=Delete[inputnodes,deletepos];
addpattrn=Delete[addpattrn,deletepos];

If[inputnodes=!={},
varsTd[sz+1]=Union[Flatten[addpattrn]];
takepattern=Position[addpattrn,#,\[Infinity]][[All,1;;2]]&/@varsTd[sz+1];
takepattern=Apply[tp,takepattern,{2}];
takepattern=Table[tp[i,j],{i,Length[addpattrn]},{j,Length[addpattrn[[i]]]}]->takepattern;
FFAlgTakeAndAdd[equGraph,equNode,inputnodes,takepattern];
,
varsTd[sz+1]={}];

(*expandvars=graphvars[[2;;-1]];*)

FFGraphOutput[equGraph,equNode];
(*expand in epsilon*)

(*testres=FFGraphEvaluate[equGraph,{13,23}];
Return[testres];*)

FFNewGraph[LaurentGraph,input,expandvars];
FFAlgLaurent[LaurentGraph,laurent,{input},equGraph,epsor,Sequence@@FilterRules[opt,Options[FFAlgLaurent]]];
FFGraphOutput[LaurentGraph,laurent];

laurentlearn=FFLaurentLearn[LaurentGraph];
(*Print[laurentlearn];*)
If[laurentlearn===$Failed,Print["FFLaurentLearn failed. Try increasing MaxDegree."];FFDeleteGraph[equGraph];FFDeleteGraph[LaurentGraph];Return[$Failed]];
If[Union[laurentlearn[[2]]]=!={epsor},Print["Error in epsor 1"];Return[$Failed]];

(*If[!AllTrue[laurentlearn[[1]],(#>=epsor)&],Print["Error in epsor 2"];Return[$Failed]];*)

If[AllTrue[laurentlearn[[1]],(#>epsor)&],Return[{}]];
nonzerovarsTd=Pick[varsTd[sz+1],LessEqual@@@Transpose[laurentlearn]];
takepattern=laurentlearn[[2]]-laurentlearn[[1]]+1;
takepattern=Table[t[i,j],{i,Length[takepattern]},{j,takepattern[[i]],1,-1}]//Flatten;
takepattern={takepattern}->Cases[takepattern,t[_,1]];
FFAlgTake[LaurentGraph,leadingOrder,{laurent},takepattern];
FFGraphOutput[LaurentGraph,leadingOrder];
(*Print[laurentlearn,FFNParsOut[LaurentGraph]];*)
(*nonzerovarsTd=Table[#[[1]],#[[2]]]&/@Transpose[{varsTd[sz+1],laurentlearn[[2]]-laurentlearn[[1]]+1}]//Flatten;*)

(*evaluate multiple points to build equation system and find its rank*)
nthreads=OptionValue["NThreads"];
If[nthreads===Automatic,nthreads=FFAutomaticNThreads[]];

prime=FFPrimeNo[0];
(*preparing first with testprime*)
testpoints=RandomInteger[{10,10000},{nthreads,Length[graphvars]-1}];
testpoints=DeleteDuplicates[testpoints];
testpoints=Transpose[Join[Transpose[testpoints],{Join[{1},ConstantArray[0,Length[testpoints]-1]]}]];
matrix=FFGraphEvaluateMany[LaurentGraph,testpoints,
Sequence@@FilterRules[opt,DeleteCases[Options[FFGraphEvaluateMany],"PrimeNo"->_]]];
primetest=matrix[[1]];
primetestpoint=testpoints[[1,1;;-2]];
matrix=RowReduce[matrix[[2;;-1]],Modulus->prime];
rank=MatrixRank[matrix,Modulus->prime];
testpoints=testpoints[[2;;-1,1;;-2]];

While[rank>=Length[matrix],
testpointsnew=RandomInteger[{10,10000},{nthreads,Length[graphvars]-1}];
testpointsnew=Complement[DeleteDuplicates[testpointsnew],testpoints];
testpoints=Join[testpoints,testpointsnew];
matrix=Join[FFGraphEvaluateMany[LaurentGraph,testpointsnew,"PrimeNo"->0,Sequence@@FilterRules[opt,DeleteCases[Options[FFGraphEvaluateMany],"PrimeNo"->_]]],matrix];
matrix=RowReduce[matrix,Modulus->prime];
rank=MatrixRank[matrix,Modulus->prime];
];
If[rank===0,Return[matrix . nonzerovarsTd]];
testpoints=Join[{primetestpoint},testpoints[[1;;rank]]];
matrix=matrix[[1;;rank]];
If[print,Print["Rank of system: ",rank]];

(*Now reconstruct the full matrix*)
(*FFNewGraph[SolveGraph];
Do[FFAlgRatNumEval[SolveGraph,tpnode[i],testpoints[[i]]],{i,rank+1}];
FFAlgRatNumEval[SolveGraph,zeronode,ConstantArray[0,(rank+1)]];

FFAlgSubgraphMap[SolveGraph,evalNode,Table[tpnode[i],{i,rank+1}],LaurentGraph];
takepattern={Flatten[Table[t[i,j],{i,(rank+1)},{j,Length[nonzerovarsTd]}]],Flatten[Table[z[i],{i,(rank+1)}]]}\[Rule]Flatten[Table[Join[Table[t[i,j],{j,Length[nonzerovarsTd]}],{z[i]}],{i,(rank+1)}]];
FFAlgTake[SolveGraph,takenode,{evalNode,zeronode},takepattern];
FFAlgNodeDenseSolver[SolveGraph,solveNode,{takenode},(rank+1),nonzerovarsTd];

FFGraphOutput[SolveGraph,solveNode];
solverlearn=FFDenseSolverLearn[SolveGraph,nonzerovarsTd];
equs=FFReconstructNumeric[SolveGraph];
equs=(#[[1]]-#[[2]])&/@FFDenseSolverSol[equs,solverlearn];*)

(*testprime=FFGraphEvaluateMany[LaurentGraph,{{13,1}},Sequence@@FilterRules[opt,Options[FFGraphEvaluateMany]]];
Print[testprime];*)

matrixrec=Map[FFRatRec[#,prime]&,matrix,{2}];
(*matrix=DeleteCases[matrix,{0..}];*)

(*Testing if first evaluation was enough*)
If[matrixrec=!={},
nulltest=NullSpace[matrixrec],
nulltest={}];
If[nulltest==={},
nulltest=False;
,
nulltest=nulltest . primetest;
nulltest=FFRatMod[#,FFPrimeNo[1]]&/@nulltest;
nulltest=!AllTrue[nulltest,(#===0&)];
];

For[i=1,nulltest&&(i<OptionValue["MaxPrimes"]),i++,
If[print,Print["Calculating prime number ",i+1]];

testpointsnew=Transpose[Join[Transpose[testpoints],{Join[{i+1},ConstantArray[i,rank]]}]];
matrixnew=FFGraphEvaluateMany[LaurentGraph,testpointsnew,
Sequence@@FilterRules[opt,DeleteCases[Options[FFGraphEvaluateMany],"PrimeNo"->_]]];
primetest=matrixnew[[1]];
matrixnew=RowReduce[matrixnew[[2;;-1]],Modulus->FFPrimeNo[i]];

matrix=Transpose[{matrix,matrixnew},{3,1,2}];
matrix=Map[ChineseRemainder[#,{prime,FFPrimeNo[i]}]&,matrix,{2}];
(*matrix=ArrayReshape[matrix,{rank,Length[nonzerovarsTd]}];*)
prime=prime FFPrimeNo[i];
matrixrec=Map[FFRatRec[#,prime]&,matrix,{2}];

nulltest=NullSpace[matrixrec] . primetest;
nulltest=FFRatMod[#,FFPrimeNo[i+1]]&/@nulltest;
nulltest=!AllTrue[nulltest,(#===0&)];
];

FFDeleteGraph[SolveGraph];
FFDeleteGraph[LaurentGraph];
FFDeleteGraph[equGraph];
(*equs=CoefficientList[varsTd[sz+1].Normal[FFLaurentSol[#,eps,laurentlearn]],eps]&/@matrix;*)

equs=matrixrec . nonzerovarsTd;
(*outvars=varsTd[sz+1];*)

If[print,Print["Total time for building equation system: ",ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];

equs];




(* ::Subsection:: *)
(*1.4.2 Solve the Equation*)


(* ::Input::Initialization:: *)
Options[solCalc]:=DeleteDuplicates[Join[{"MaxIterations"->20,"StartIterations"->1,"ExtraCheck"->False,"TakeRow"->1},Options[phiCalc],DeleteCases[Options[FFDenseSolve],("MaxDegree"->_)|("MaxPrimes"->_)]]];
solCalc[psi1_,ansatzF_,degsIn_,xv_,eps_,OptionsPattern[]]:=Module[{phi1,sol1,oldsol,equ,equs,sz,graphvars,allVars,startit,maxit,
phi2,opt,print,time,a,ts,Qtest,Tvars,
extracheck,checksol,oldsolequs,tr,phi1temp,allVarstemp,i,equtemp,equstemp,degs,phiold,allVarsold},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[solCalc];
(*checking input*)
(*ansatzF=OptionValue["AnsatzFunctions"];
If[ansatzF===Automatic,
ansatzF=letters];*)
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[psi1]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
graphvars=Join[{eps,xv},DeleteCases[graphvars,eps|xv]];

If[!SquareMatrixQ[psi1],Message[nthO::badmatrix1];Return[$Failed]];

print=!OptionValue["Silent"];
extracheck=OptionValue["ExtraCheck"];
time=SessionTime[];
sz=Length[psi1];
startit=OptionValue["StartIterations"];
maxit=OptionValue["MaxIterations"];
(*initialize loop*)
oldsol={};
ts={};

tr=OptionValue["TakeRow"];
If[Head[tr]=!=List,
tr={{tr,sz}}];
degs=degsIn;
If[Head[degs[[1]]]===Integer,
degs={degs}];

allVarsold=phiold=Table[{},Length[tr],sz];

Catch[
Do[If[print,Print["Starting epsilon order ",epsor]];

(*calculate step-by-step the phi-matrix, build the equation system and solve it*)

If[print,Print["Calculating phi."]];
phi1={};
allVars={};

Do[
phi1temp=phiCalc[ansatzF,phiold[[i]],allVarsold[[i]],oldsol,xv,eps,Sequence@@DeleteCases[FilterRules[opt,Options[phiCalc]],("TakeRow"->_)|("MaxDerivatives"->_)],"TakeRow"->tr[[i,1]],"MaxDerivatives"->tr[[i,2]]];
If[phi1temp===$Failed,Throw[Return[$Failed]]];
{phi1temp,allVarstemp}=phi1temp;
phi1=Append[phi1,phi1temp];
allVars=Append[allVars,allVarstemp];
,{i,Length[tr]}];
(*Print[{allVars,oldsol}];*)
phiold=phi1;
allVarsold=allVars;
phi1=Join@@phi1;
allVars=Join@@allVars;

If[print,Print["Calculating the equation."]];
equ={};
equs={};
Do[
equtemp=equStep[psi1,phi1,epsor,allVars,degs[[i]],eps,Sequence@@DeleteCases[FilterRules[opt,Options[equStep]],("TakeRow"->_)|("MaxDerivatives"->_)],"TakeRow"->tr[[i,1]]];
If[equtemp===$Failed,Throw[Return[$Failed]]];
equstemp=CoefficientList[equtemp,graphvars]//Flatten//Union;
(*new, need to test: apply old solution to equations*)
equstemp=Collect[equstemp//.oldsol,T[__],Together];
equ=Join[equ,equtemp];
equs=Join[equs,equstemp];
,{i,Length[tr]}];

equs=DeleteCases[equs,0];
If[equs==={},Continue[]];

oldsolequs=(#[[1]]-#[[2]])&/@(oldsol/.{a_T:>a[[1;;-2]]});
oldsolequs=Collect[oldsolequs,T[__],Together];
equs=Join[equs,oldsolequs];

Tvars=Union[Cases[equs,T[__],\[Infinity]]];
(*Tvars=Join[ts,Tvars]//DeleteDuplicates//Reverse;*)
Tvars=DeleteCases[Tvars,T[_]];
Tvars=Tvars//Reverse;

sol1=FFDenseSolve[equs==0//Thread,Tvars,Sequence@@FilterRules[opt,Options[FFDenseSolve]]];
If[sol1===$Failed,Print["Solver failed. Try increasing MaxDegree or MaxPrimes."];Throw[Return[$Failed]]];
If[sol1===FFMissingPrimes,Print["Solver failed. Try increasing MaxPrimes."];Throw[Return[$Failed]]];

(*If[MemberQ[sol1,T[_]->0],Break[]];*)
If[(sol1===FFImpossible)||(Flatten[sol1]==={}),Print["No solution found."];Throw[Return[$Failed]]];
(*Now format the solution*)
sol1=List@@@sol1;
sol1[[All,1]]=sol1[[All,1]]/.T[b__]->T[b,\!\(\*
TagBox[
StyleBox[
RowBox[{"Pattern", "[", 
RowBox[{"a", ",", 
RowBox[{"BlankNullSequence", "[", "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)];
sol1[[All,2]]=sol1[[All,2]]/.T[b__]->T[b,a];
sol1=Rule@@@sol1;
(*sol1=Join[oldsol,sol1];*)
checksol=sol1;
Quiet[Check[
sol1[[All,2]]=sol1[[All,2]]//.sol1;
,Throw[Return[$Failed]]]];
sol1=sol1//Union;
oldsol=sol1;

ts=Union[Cases[sol1[[All,2]],T[__],\[Infinity]]][[All,1;;-2]];

ts=Table[Append[#,m[i]]&/@ts//.sol1,{i,Length[ansatzF]}];
ts=Union[Cases[ts,T[__],\[Infinity]]];

If[print,Print["Number of independent variables found so far: ",Length[ts]]];

If[Length[ts]===sz,
Qtest=Table[Append[#,m[i]],{i,Length[ansatzF]}]&/@ts//.sol1;
Qtest=Union[Cases[Qtest,T[__],\[Infinity]]];
If[Length[Qtest]===sz,
If[Qtest===ts,
If[(extracheck===False)||(extracheck===(epsor-1)),
Break[],
extracheck=epsor],
Print["Error in variables."];Throw[Return[$Failed]]]
,
If[print,Print["Number of relations missing: ",Length[Qtest]-sz]];
];
];
If[epsor===maxit,
Print["No solution found after ",maxit, " iterations, returning partial solution"];
Throw[Return[sol1]];
];
,
{epsor,startit,maxit}];
];

If[print,Print["Total time for finding the solution: ",ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];
(*If a sensible solution has been found, build the full phi-matrix, else only the solution*)
(*If[Length[ts]==sz,
Print[allVars];
phi2=Table[phi1[[i]].allVars[[i]],{i,sz}];
ts=ts->IdentityMatrix[sz]//Thread;
phi2=phi2/.ts];*)
(*Return[phi2],
Return[sol1]];*)
sol1];


(* ::Section:: *)
(*1.5 Functions for analyzing the result*)


(* ::Subsection::Closed:: *)
(*1.5.1 FFInvMatMul*)


(* ::Input::Initialization:: *)
Options[FFInvMatMul] := FilterRules[Join[DeleteCases[Options[FFDenseSolve],("MaxDegree"->_)|("MaxPrimes"->_)],{"Sparse"->False,"InvertInput"->0,"MaxDegree"->1000,"MaxPrimes"->150}], Except["NeededVars"|"IndepVarsOnly"]];
FFInvMatMul[matricesIn__, OptionsPattern[]]:=Module[
    {eqs, len,leninv, varx, vary, varsx, varsy, sol, params,res,graph,in,learn,sparse,opt,dim,m,matmul,inv,matrices,i,test,minv},
    opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFInvMatMul];
inv=OptionValue["InvertInput"];
sparse=OptionValue["Sparse"];

matrices={matricesIn};
    len=Length[matrices];
If[(inv<0)||(inv>len),Print["Wrong \"InvertInput\" "];Return[$Failed]];
(*If[len===1,Return[matricesIn]];*)

Do[dim[i]=Dimensions[matrices[[i]]],{i,len}];

test=Table[dim[i][[2]]===dim[i+1][[1]],{i,len-1}];
    If[!(And@@test), Print["Matrices have wrong dimensions for multiplication"]; Return[$Failed]];

FFNewGraph[graph];

    res = Catch[
      params = OptionValue["Parameters"];
      params = If[params===Automatic,
                  Variables[matrices],
                  params];
      If[params =!= {}, FiniteFlow`Private`CheckVariables[params]];
      FFGraphInputVars[graph,in,params];

(*Inverting the first matrix*)
If[inv=!=0,

minv=matrices[[inv]];
If[dim[inv][[1]]=!=dim[inv][[2]],Print["Matrix for inverse not a square matrix"];Return[$Failed]];

    leninv= Length[minv];
    
    varsx = varx/@Range[leninv];
    varsy = vary/@Range[leninv];
    
    eqs = Table[minv[[ii]] . varsx-vary[ii]==0,{ii,leninv}];
    
      res = If[sparse,
               FFAlgSparseSolver[graph,m[inv],{in},params,eqs,Join[varsx,varsy],
                                    "VarsPattern"->(DeleteDuplicates[Cases[{#},(_varx|_vary),Infinity]]&),
                                     Sequence@@FilterRules[{opt}, Options[FFAlgSparseSolver]]],
               FFAlgDenseSolver[graph,m[inv],{in},params,eqs,Join[varsx,varsy],
                                   Sequence@@FilterRules[{opt}, Options[FFAlgDenseSolver]]]];
      If[res==$Failed,Throw[$Failed]];
      FFSolverOnlyHomogeneous[graph,m[inv]];
     
      FFGraphOutput[graph,m[inv]];
      learn = If[sparse,
                 FFSparseSolverLearn[graph,Join[varsx,varsy]],
                 FFDenseSolverLearn[graph,Join[varsx,varsy]]];
      If[TrueQ[learn==FFImpossible],Throw[FFSingularMatrix]];
      If[!TrueQ[learn[[0]]==List],Throw[learn]];
      If[!TrueQ[("DepVars"/.learn) == varsx && ("IndepVars"/.learn) == varsy],Throw[FFSingularMatrix]];
];
(*Do the matrix multiplication *)

Do[If[params==={},
FFAlgRatNumEval[graph,m[i], Join@@(matrices[[i]])],
FFAlgRatFunEval[graph,m[i],{in},params, Join@@(matrices[[i]])]
],
{i,DeleteCases[Range[len],inv]}];

If[len===1,
matmul[1]=m[1],
FFAlgMatMul[graph,matmul[2],{m[1],m[2]},dim[1][[1]],dim[1][[2]],dim[2][[2]]]
];

Do[
FFAlgMatMul[graph,matmul[i+1],{matmul[i],m[i+1]},dim[1][[1]],dim[i][[2]],dim[i+1][[2]]],
{i,2,len-1}];
      
      FFGraphOutput[graph,matmul[len]];
      
       res = If[TrueQ[params == {}],
               FFReconstructNumeric[graph, Sequence@@FilterRules[{opt}, Options[FFReconstructNumeric]]],
               If[TrueQ[Length[params]==1],
FFParallelReconstructUnivariate[graph,params, Sequence@@FilterRules[{opt}, Options[FFParallelReconstructUnivariate]]],
FFReconstructFunction[graph,params, Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
]
             ];
       If[!TrueQ[res[[0]]==List],Throw[res]];
       ArrayReshape[res,{dim[1][[1]],dim[len][[2]]}]
    ];
    FFDeleteGraph[graph];
    res
];


(* ::Subsection::Closed:: *)
(*1.5.2 basisChange*)


(* ::Input::Initialization:: *)
Options[basisChange]:=Join[{"DerivativeRules"->{},"MaxDegree"->1000,"MaxPrimes"->150},Options[makePsiInvGraph],DeleteCases[Options[FFReconstructFunction],("MaxDegree"->_)|("MaxPrimes"->_)]];
basisChange[amatrix_,Tmatrix_,xv_,OptionsPattern[]]:=Module[{x,graphvars,sz,Tgraph,Anode,TAnode,AddNode,outnode,opt,learn,
pA,qA,der1,derNode1,derNode21,derNode22,derNode23,derNode2,atest,Dreps},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[basisChange];
(*testing input*)
If[!SquareMatrixQ[amatrix],Message[nthO::badmatrix1];Return[$Failed]];
If[!SquareMatrixQ[Tmatrix],Message[nthO::badmatrix2];Return[$Failed]];

Dreps=Flatten[{OptionValue["DerivativeRules"]}];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[{amatrix,Tmatrix}]];
(*If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];*)
(*xv=graphvars[[2]];*)
sz=Length[amatrix];
graphvars=Union[graphvars,Dreps[[All,1]]];

(*make graph for computing Tinv*)
learn=makePsiInvGraph[Tmatrix,Tgraph,Sequence@@DeleteCases[FilterRules[opt,Options[makePsiInvGraph]],("Variables"->_)],"Variables"->graphvars];
If[learn===FFImpossible,Message[nthO::singmatrix];Return[$Failed]];

(*Now compute the derivatives*)
{pA,qA}=Transpose[NumeratorDenominator[Flatten[Tmatrix]]];
der1=D[pA,xv]/qA;

FFAlgRatFunEval[Tgraph,Anode,{input},graphvars,Flatten[amatrix]];
FFAlgRatFunEval[Tgraph,derNode1,{input},graphvars,der1];
(*time=SessionTime[];*)
FFAlgRatFunEval[Tgraph,derNode21,{input},graphvars,1/qA];
FFAlgRatFunEval[Tgraph,derNode22,{input},graphvars,-pA];
FFAlgRatFunEval[Tgraph,derNode23,{input},graphvars,D[qA,xv]];
FFAlgMul[Tgraph,derNode2,{derNode21,derNode21,derNode22,derNode23}];

(*Now do the matrix multiplication and add with the derivative*)
FFAlgMatMul[Tgraph,TAnode,{psinode,Anode},sz,sz,sz];
FFAlgAdd[Tgraph,AddNode,{derNode1,derNode2,TAnode}];
FFAlgMatMul[Tgraph,outnode,{AddNode,psiInvnode},sz,sz,sz];

FFGraphOutput[Tgraph,outnode];
atest=FFReconstructFunction[Tgraph,graphvars,Sequence@@FilterRules[opt,Options[FFReconstructFunction]]];
If[atest===$Failed,Print["Reconstruction failed. Try increasing MaxDegree."];FFDeleteGraph[Tgraph];Return[$Failed]];
If[atest===FFMissingPrimes,Print["Reconstruction failed. Try increasing MaxPrimes."];FFDeleteGraph[Tgraph];Return[$Failed]];
atest=ArrayReshape[atest,{sz,sz}];

FFDeleteGraph[Tgraph];

atest/.Dreps];


(* ::Subsection::Closed:: *)
(*1.5.3 denominators*)


(* ::Input::Initialization:: *)
getDenomFacs[in_List]:=getDenomFacs/@in;
getDenomFacs[in_Plus]:=getDenomFacs/@(List@@in);
getDenomFacs[in:Except[_List|_Plus]]:=Union[FactorList[Denominator[in]][[All,1]]];


(* ::Input::Initialization:: *)
denominators[amatrix_]:=Select[Union[Flatten[getDenomFacs[amatrix]]],Variables[#]=!={}&];
denominators[amatrix_,xv_]:=Select[denominators[amatrix],MemberQ[{#},xv,Infinity]&];
denominators[amatrix_,xv_,eps_]:=Select[denominators[amatrix,xv],FreeQ[{#},eps,Infinity]&];


(* ::Section:: *)
(*1.6 Combining everything into one function*)


(* ::Subsection::Closed:: *)
(*1.6.1 BCalc*)


(* ::Input::Initialization:: *)
Options[BCalc]:={"Silent"->False,"Variables"->Automatic};
BCalc[sol1_,ansatzF_,sz_,x_,eps_,OptionsPattern[]]:=Module[{ts,ms,Qtest,QtestTs,bmatrix,print,temp,graphvars},

(*ansatzF=OptionValue["AnsatzFunctions"];
If[ansatzF===Automatic,
ansatzF=letters];*)
print=!OptionValue["Silent"];
graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[ansatzF]];
graphvars=Join[{eps,x},DeleteCases[graphvars,eps|x]];
(*eps=graphvars[[1]];
x=graphvars[[2]];*)
(*If[OptionValue["AnsatzFunctions"]===Automatic,
ansatzF=D[Log[#],x]&/@ansatzF];*)

ts=Cases[sol1[[All,2]],T[__],\[Infinity]]//Union;

If[Length[ts]===0,
Print["Relations missing."];Return[$Failed]];

ts=ts[[All,1;;-2]];
ts=Table[Append[#,m[i]]&/@ts//.sol1,{i,Length[ansatzF]}];
ts=Union[Cases[ts,T[__],\[Infinity]]];
(*sz=Length[ts];*)
If[(!(Length[ts]===sz)),
Print["Relations missing."];Return[$Failed]];

Qtest=Table[Append[#,m[i]]&/@ts//.sol1,{i,Length[ansatzF]}];
QtestTs=Union[Cases[Qtest,T[__],\[Infinity]]];

If[(!(Length[QtestTs]===sz))||(Length[QtestTs]===0),
Print["Relations missing."];Return[$Failed]];
If[!(QtestTs===ts),
Print["Relations not closed.";Return[$Failed]]
];

If[print,Print["Calculating m's."]];
ms=Table[temp=CoefficientArrays[Qtest[[i]],ts]//Normal;
If[Union[temp[[1]]]=!={0},Print["Error in Qtest."]];
If[Length[temp]===1,
If[print,Print["Result is independent of letter ",i]];
ConstantArray[0,{sz,sz}],
temp[[2]]],
{i,Length[ansatzF]}];

If[print,Print["Building B-matrix."]];
bmatrix=Sum[Together[eps ansatzF[[i]]]ms[[i]],{i,Length[ansatzF]}]//Together;

bmatrix];


(* ::Subsection::Closed:: *)
(*1.6.2 TCalc*)


(* ::Input::Initialization:: *)
Options[TCalc]:=DeleteCases[Join[Options[checkDegrees],Options[solCalc],Options[FFInvMatMul]]//DeleteDuplicates,"InvertInput"->_];
TCalc[Am_,ansatzF_,xv_,eps_,OptionsPattern[]]:=Module[{graphvars,opt,sz,psi1,time,print,degs,phi2,Tm,sol1,bmatrix,timetot,pmax,tr,i},
opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[TCalc];
print=!OptionValue["Silent"];

timetot=SessionTime[];

If[!SquareMatrixQ[Am],Message[nthO::badmatrix1];Return[$Failed]];

graphvars=OptionValue["Variables"];
If[graphvars===Automatic,graphvars=Variables[Am]];
If[Length[graphvars]<1,Message[nthO::badvars];Return[$Failed]];
graphvars=Join[{eps,xv},DeleteCases[graphvars,eps|xv]];
If[print,Print["Variables: ",graphvars]];

(*letters=OptionValue["Letters"];
ansatzF=OptionValue["AnsatzFunctions"];
If[ansatzF===Automatic,
If[letters===Automatic,
If[print,Print["No letters given. Using the following letters:"]];
letters=FactorList/@(Denominator/@Flatten[Am]);
letters=Flatten[letters,1][[All,1]]//Union;
letters=Select[letters,FreeQ[{#},graphvars[[1]],\[Infinity]]&&MemberQ[{#},Alternatives@@graphvars[[2;;-1]],\[Infinity]]&];
If[print,Print[letters]];
,
If[print,Print["Letters: ",letters]];
];
,
If[print,Print["AnsatzFunctions: ",ansatzF]];
];*)

sz=Length[Am];
If[print,Print["Size of the matrix: ",sz]];

tr=OptionValue["TakeRow"];
If[Head[tr]=!=List,tr={{tr,sz}}];

If[print,Print["Calculating psi."]];
psi1=psiCalc[Am,xv,Sequence@@FilterRules[opt,Options[psiCalc]]];
If[psi1===$Failed,Return[$Failed]];

If[print,Print["Calculating degrees."]];
time=SessionTime[];
degs=Table[checkDegrees[psi1,eps,Sequence@@DeleteCases[FilterRules[opt,Options[checkDegrees]],("TakeRow"->_)],"TakeRow"->tr[[i,1]]],{i,Length[tr]}];
If[print,Print[ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];
If[MemberQ[degs,$Failed],Return[$Failed]];
(*Print[degs];*)
(*checking the degrees*)
pmax=sz (sz+1)/2;
pmax=And@@Flatten[Table[degs[[i,2,l+1]]<=(pmax-l),{i,Length[degs]},{l,0,sz}]];
If[!pmax,Print["Degree check failed, the integral cannot be a UT integral."];Return[$Failed]];
Catch[
Do[
If[degs[[i,1]]<1,Print["b-tilde[0](0) not vanishing. Try a different normalization for integral ",i];Throw[Return[$Failed]]];
,{i,Length[degs]}];
];

If[print,Print["Calculating solution."]];
sol1=solCalc[psi1,ansatzF,degs,xv,eps,Sequence@@DeleteCases[FilterRules[opt,Options[solCalc]],("TakeRow"->_)],"TakeRow"->tr];
If[sol1===$Failed,Return[$Failed]];
(*If[Length[Union[Cases[sol1[[All,2]],T[__],\[Infinity]]]]=!=sz,
Print["Full solution not found, returning only partial solution"];Return[sol1]];*)

If[print,Print["Calculating full phi-matrix."]];
bmatrix=BCalc[sol1,ansatzF,sz,xv,eps,Sequence@@FilterRules[opt,Options[BCalc]]];
If[bmatrix===$Failed,
Print["Full solution not found, returning only partial solution"];Return[sol1]];
phi2=psiCalc[bmatrix,xv,Sequence@@FilterRules[opt,Options[psiCalc]]];
If[bmatrix===$Failed,Return[$Failed]];

If[print,Print["Calculating T-matrix."]];
time=SessionTime[];
Tm=FFInvMatMul[phi2,psi1,Sequence@@FilterRules[Join[{"InvertInput"->1},opt],Options[FFInvMatMul]]];
If[print,Print[ToString[SetAccuracy[SessionTime[]-time,3]]," s"]];
If[Tm===$Failed,Return[$Failed]];
If[print,Print["Done"]];
If[print,Print["Total time for all operations: ",ToString[SetAccuracy[SessionTime[]-timetot,3]]," s"]];
(*If[print,Print["Maximum memory used: ",ToString[SetAccuracy[N[MaxMemoryUsed[]/10^6],3]]," MB"]];*)

Tm];


(* ::Chapter:: *)
(*End*)


(* ::Input::Initialization:: *)
End[];


(* ::Input::Initialization:: *)
Protect@@functionNames;
Clear[functionNames];
Remove[functionNames];


(* ::Input::Initialization:: *)
EndPackage[];
