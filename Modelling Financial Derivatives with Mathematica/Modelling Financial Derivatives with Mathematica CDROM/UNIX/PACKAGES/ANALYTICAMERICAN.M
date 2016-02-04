(* :Title: 
	Approximate American Options *)

(* :Context: 
	Derivatives`AnalyticAmerican` *)

(* :Author: 
	Dr William T. Shaw *)

(* :Summary: 
	*)

(* :Package Version: 
	1.0 *)

(* :Mathematica Version: 
	3.0 *)

(* :History:
	V.1.0, January 1997, by William T. Shaw *)


BeginPackage["Derivatives`AnalyticAmerican`", "Derivatives`BlackScholes`"]

(* Usage Information *)

Front::usage = 
"Front[strike, interest, vol, t] returns the critical price frontier for a vanilla zero-dividend American Put according to an approximation of the Carr-Jarrow-Myneni integral equation."

AmericanPutCJM::usage = 
"AmericanPutCJM[price,strike,vol,r,t]  returns an APPROXIMATE value for a vanilla zero-dividend American Put according to an approximation of the Carr-Jarrow-Myneni integral equation."


AmericanPutDeltaCJM::usage = 
"AmericanPutDeltaCJM[price,strike,vol,r,t]  returns an APPROXIMATE value for the delta of a vanilla zero-dividend American Put according to an approximation of the Carr-Jarrow-Myneni integral equation."


AmericanPutGammaCJM::usage = 
"AmericanPutGammaCJM[price,strike,vol,r,t]  returns an APPROXIMATE value for the gamma of a vanilla zero-dividend American Put according to an approximation of the Carr-Jarrow-Myneni integral equation."


(* END OF USAGE INFORMATION *)

Begin["`Private`"]

Front[strike_, r_, sigma_, t_] :=
If[t==0,strike,Module[{l,m},
l= m/. FindRoot[Exp[-m*sigma*Sqrt[t] - (r+sigma^2/2)*t]*(1+Erf[-m/Sqrt[2]]) \
==
 r*Sqrt[2/(r+(sigma - (r/(2*sigma)+sigma/4) +m/(2*Sqrt[t]))^2/2)]*Erf[Sqrt[(
r+(sigma-(r/(2*sigma)+sigma/4)+m/(2*Sqrt[t]))^2/2)*t]]/sigma,
{m,-(sigma/2+r/sigma)*Sqrt[t]},MaxIterations -> 30];
strike*Exp[-l*sigma*Sqrt[t] - (r+sigma^2/2)*t]]]

AmericanPutCJM[price_,strike_,sigma_,r_,t_] :=
If[price <= Front[strike, r, sigma, t], strike - price,
BlackScholesPut[price,strike,sigma,r,0,t]+r*strike*NIntegrate[Exp[-r*u]*Norm[(\

Log[Front[strike,r,sigma,t-u]/price]-(r-sigma^2/2)*u)/(sigma*Sqrt[u])],
{u,0,t}]]

AmericanPutDeltaCJM[price_,strike_,sigma_,r_,t_] :=
If[price <= Front[strike, r, sigma, t],  -1,
N[BlackScholesPutDelta[price,strike,sigma,r,0,t] - (r*strike/(price*sigma*
Sqrt[2.*Pi]))*NIntegrate[Exp[-r*u]*Exp[-((Log[Front[strike,r,sigma,t-u]/
price]-(r-sigma^2/2)*u)/(sigma*Sqrt[u]))^2/2]/Sqrt[u],
{u, 0, t}]]]

AmericanPutGammaCJM[price_,strike_,sigma_,r_,t_] :=
If[price <= Front[strike, r, sigma, t], 0,
N[
BlackScholesPutGamma[price,strike,sigma,r,0,t]+(r*strike/(price^2*sigma*Sqrt[
2.*Pi]))*NIntegrate[Exp[-r*u]*Exp[-((Log[Front[strike,r,sigma,t-u]/price]-(r-
sigma^2/2)*u)/(sigma*Sqrt[u]))^2/2]*(1/Sqrt[u] - (Log[Front[strike,r,sigma,t-
u]/price]-(r-sigma^2/2)*u)/(sigma*Sqrt[u])/(sigma*u)),
{u, 0, t}]]]


End[]

EndPackage[]
