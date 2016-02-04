(* :Title: 
	Basic Options *)

(* :Context: 
	Derivatives`BlackScholes` *)

(* :Author: 
	Dr William T. Shaw *)

(* :Summary: 
	*)

(* :Package Version: 
	1.1 *)

(* :Mathematica Version: 
	3.0, but compatible with 2.X *)

(* :History:
	V.1.0, January 1996, by William T. Shaw;
    V.1.1, December 1997, by William T. Shaw (spelling corrected) *)

BeginPackage["Derivatives`BlackScholes`"]

(* Usage Information *)

CallPayoff::usage = 
"CallPayoff[price, strike] returns the payoff value of a vanilla 
European Call."

PutPayoff::usage = 
"PutPayoff[price, strike] returns the payoff value of a vanilla 
European Put."

BlackScholesCall::usage = 
"BlackScholesCall[price, strike, vol, riskfree, divyield, expiry] 
returns the Black-Scholes value of a vanilla European Call."

BlackScholesPut::usage = 
"BlackScholesPut[price, strike, vol, riskfree, divyield, expiry] 
returns the Black-Scholes value of a vanilla European Put."

BlackScholesCallDelta::usage = 
"BSCallDelta[price, strike, vol, riskfree, divyield, expiry] 
returns the delta of a vanilla European Call."

BlackScholesPutDelta::usage = 
"BlackScholesPutDelta[price, strike, vol, riskfree, divyield, expiry] 
returns the delta of a vanilla European Put."

BlackScholesCallTheta::usage = 
"BlackScholesCallTheta[price, strike, vol, riskfree, divyield, expiry] 
returns the theta of a vanilla European Call."

BlackScholesPutTheta::usage = 
"BlackScholesPutTheta[price, strike, vol, riskfree, divyield, expiry] 
returns the theta of a vanilla European Put."

BlackScholesCallRho::usage = 
"BlackScholesCallRho[price, strike, vol, riskfree, divyield, expiry] 
returns the rho of a vanilla European Call."

BlackScholesPutRho::usage = 
"BlackScholesPutRho[price, strike, vol, riskfree, divyield, expiry] 
returns the rho of a vanilla European Put."

BlackScholesCallGamma::usage = 
"BlackScholesCallGamma[price, strike, vol, riskfree, divyield, expiry] 
returns the gamma of a vanilla European Call."

BlackScholesPutGamma::usage = 
"BlackScholesPutGamma[price, strike, vol, riskfree, divyield, expiry] 
returns the gamma of a vanilla European Call."

BlackScholesCallVega::usage = 
"BlackScholesCallVega[price, strike, vol, riskfree, divyield, expiry] 
returns the vega (lambda) of a vanilla European Call."

BlackScholesPutVega::usage = 
"BlackScholesPutVega[price, strike, vol, riskfree, divyield, expiry] 
returns the vega (lambda) of a vanilla European Put."

Norm::usage = 
"Norm[x] returns the cumulative normal distribution function evaluated at x."

BlackScholesCallImpVol::usage = 
"BlackScholesCallImpVol[price, strike, riskfree, divyield, expiry, \
optionprice] 
returns the implied volatility of a vanilla European Call."

BlackScholesPutImpVol::usage = 
"BlackScholesPutImpVol[price, strike, riskfree, divyield, expiry, \
optionprice] 
returns the implied volatility of a vanilla European Call."



(* END OF USAGE INFORMATION *)

Begin["`Private`"]

Norm[z_?NumberQ]:= N[0.5 + 0.5 Erf[z/Sqrt[2]]];

Norm[x_] := (1 + Erf[x / Sqrt[2]]) / 2;

dtwo[s_, v_, k_, t_, r_, q_] :=(Log[s / k] + (r-q) t) / (v Sqrt[t]) -1/2 v \
Sqrt[t];
   
done[s_, v_, k_, t_, r_, q_] :=(Log[s / k] + (r-q) t) / (v Sqrt[t]) + 1/2 v \
Sqrt[t];
   
BlackScholesCall[s_,k_, v_, r_, q_, t_] := 
s Exp[-q*t] Norm[done[s, v, k, t, r, q]] - k*Exp[-r*t] Norm[dtwo[s, v, k, t, \
r, q]];

BlackScholesPut[s_,k_, v_, r_, q_, t_]:= 
k*Exp[-r*t] Norm[-dtwo[s, v, k, t, r, q]]-s Exp[-q*t] Norm[-done[s, v, k, t, \
r, q]];

BlackScholesCallDelta[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesCall[s,k, v, r, q, t], s]];

BlackScholesPutDelta[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesPut[s,k, v, r, q, t], s]];

BlackScholesCallGamma[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesCall[s,k, v, r, q, t], {s, 2}]];

BlackScholesPutGamma[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesPut[s,k, v, r, q, t], {s, 2}]];

BlackScholesCallTheta[s_,k_, v_, r_, q_, t_]=
-Evaluate[D[BlackScholesCall[s,k, v, r, q, t], t]];

BlackScholesPutTheta[s_,k_, v_, r_, q_, t_]=
-Evaluate[D[BlackScholesPut[s,k, v, r, q, t], t]];

BlackScholesCallRho[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesCall[s,k, v, r, q, t], r]];

BlackScholesPutRho[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesPut[s,k, v, r, q, t], r]];

BlackScholesCallVega[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesCall[s,k, v, r, q, t], v]];

BlackScholesPutVega[s_,k_, v_, r_, q_, t_]= 
Evaluate[D[BlackScholesPut[s,k, v, r, q, t], v]];

CallPayoff[s_, k_] := Max[0, s-k];

PutPayoff[s_, k_] := Max[0, k-s];

BlackScholesCallImpVol[s_,k_,r_,q_,t_,optionprice_] := 
sd /. FindRoot[BlackScholesCall[s,k, sd, r, q, t]== optionprice,{sd,0.2}];

BlackScholesPutImpVol[s_,k_,r_,q_,t_,optionprice_] := 
sd /. FindRoot[BlackScholesPut[s,k, sd, r, q, t]== optionprice,{sd,0.2}];

End[]

EndPackage[]
