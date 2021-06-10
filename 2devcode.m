(* ::Package:: *)

T=77; \[Alpha]0=0.6;\[Alpha]1=1.6;\[Omega]c=53.3;hmcons=2\[Pi]/53.3;

\[Mu]A=\[Mu]B=1;\[Mu]vibAe=\[Mu]vibBe=\[Mu]vibAg=\[Mu]vibBg=1;\[Epsilon]g=0;

H=Table[0,{i,1,6},{j,1,6}]
\[Omega]A=1030;\[Omega]\[Omega]A=950;\[Omega]B=1030;\[Omega]\[Omega]B=950;\[Epsilon]A=12000;\[Epsilon]B=12900;\[CapitalSigma]=0.0025;
S00=Exp[-\[CapitalSigma]/2];Sg0e1=Sqrt[\[CapitalSigma]]Exp[-\[CapitalSigma]/2];Sg1e0=-Sqrt[\[CapitalSigma]]Exp[-\[CapitalSigma]/2];S11=Exp[-\[CapitalSigma]/2](1-\[CapitalSigma]);J= hmcons 250;

H[[1,1]]=hmcons \[Epsilon]A;H[[4,4]]=hmcons \[Epsilon]B;H[[2,2]]=hmcons (\[Epsilon]A+\[Omega]\[Omega]A); H[[5,5]]=hmcons (\[Epsilon]B+\[Omega]\[Omega]B);H[[3,3]]=hmcons (\[Epsilon]A+\[Omega]B);H[[6,6]]=hmcons (\[Epsilon]B+\[Omega]A);H[[1,2]]=H[[2,1]]=H[[1,3]]=H[[3,1]]=H[[2,3]]=H[[3,2]]=H[[4,5]]=H[[5,4]]=H[[4,6]]=H[[6,4]]=H[[5,6]]=H[[6,5]]=0;H[[1,4]]=H[[4,1]]=J S00^2;H[[1,5]]=H[[5,1]]= J S00 Sg0e1;H[[1,6]]=H[[6,1]]=J Sg1e0 S00;H[[2,4]]=H[[4,2]]=J Sg0e1 S00;H[[2,5]]=H[[5,2]]=J Sg0e1^2;H[[2,6]]=H[[6,2]]=J S11 S00;H[[3,4]]=H[[4,3]]=J S00 Sg1e0;H[[3,5]]=H[[5,3]]=J S00 S11;H[[3,6]]=H[[6,3]]=J Sg1e0 Sg1e0;
\[Mu]el[1]=\[Mu]A S00;\[Mu]el[2] = \[Mu]A Sg0e1; \[Mu]el[3] = 0;\[Mu]el[4] = \[Mu]B S00; \[Mu]el[5] = \[Mu]B Sg0e1; \[Mu]el[6] = 0;
\[Mu]vib = Table[0,{i,1,6},{j,1,6}]
\[Mu]vib[[1,2]]=\[Mu]vib[[2,1]]=\[Mu]vibAe; \[Mu]vib[[1,3]]=\[Mu]vib[[3,1]]=\[Mu]vibBg;\[Mu]vib[[4,5]]=\[Mu]vib[[5,4]]=\[Mu]vibBe;\[Mu]vib[[4,6]]=\[Mu]vib[[6,4]]=\[Mu]vibAg;
vector = Chop[Eigenvectors[H]];conjvec=ConjugateTranspose[vector];energy = Eigenvalues[H];
\[Mu]elg[m_]:=Sum[vector[[m,i]]\[Mu]el[i],{i,1,6}]
Table[\[Mu]elg[m],{m,1,6}]
\[Mu]vibex[m_,n_]:=Sum[conjvec[[i,m]]vector[[n,j]]\[Mu]vib[[i,j]],{i,1,6},{j,1,6}]
SS=Table[Chop[\[Mu]vibex[m,n]],{m,1,6},{n,1,6}]


\[Lambda]Ae=50;\[Lambda]Be=50;\[Lambda]Av=5;\[Lambda]Bv=5;\[Omega]cAe=\[Omega]cBe=\[Omega]cAv=\[Omega]cBv=\[Omega]c=53.3;T = 77;
\[Lambda][1]=\[Lambda]Ae + \[Alpha]0^2 \[Lambda]Av;\[Lambda][2]=\[Lambda]Ae + \[Alpha]1^2 \[Lambda]Av;\[Lambda][3]=\[Lambda]Ae + \[Alpha]0^2 \[Lambda]Av + \[Lambda]Bv;\[Lambda][4]=\[Lambda]Be + \[Alpha]0^2 \[Lambda]Bv;\[Lambda][5]=\[Lambda]Be + \[Alpha]1^2 \[Lambda]Bv;\[Lambda][6]=\[Lambda]Be + \[Alpha]0^2 \[Lambda]Bv + \[Lambda]Av;
\[Tau]c = N[1/(53.3 3 10^10)];kB=1.3806 10^-23 ; hc= 1.05457266/10^34 ; 
coeff1 = hc 2\[Pi]/(2 kB T \[Tau]c)
coeff2[i_] := 8 \[Pi] (2\[Pi] \[Lambda][i]3 10^10\[Tau]c)(kB T \[Tau]c/hc)
\[Nu] = 2\[Pi] kB T \[Tau]c/hc
pp[t_]:=Sum[(Exp[-(n*\[Nu])*t]+(n*\[Nu])*t-1)/((n*\[Nu])*((\[Nu]*n)^2-4 Pi^2)),{n,1,2000}];
\[Alpha][i_]:=\[Lambda][i]/\[Omega]c
\[Phi]re[i_,t_]:=(\[Alpha][i] Cot[coeff1](Exp[-2 Pi t]+2 Pi t-1)+coeff2[i]*pp[t]);
\[Phi]im[i_,t_]:=- (\[Alpha][i]*(Exp[-2 Pi t]+2 Pi t-1))
\[Phi][i_]:=Table[{t,\[Phi]re[i,t]+I \[Phi]im[i,t]},{t,0,4,0.003125}]
\[Phi]tbl=Table[\[Phi][i],{i,1,6}];
\[Phi]total[i_,t_]:=\[Phi]tbl[[i]][[IntegerPart[t/0.003125]+1,2]]
Haml[k_]:=Table[If[i==j&&j==k,H[[i,j]]+Q,H[[i,j]]],{i,1,6},{j,1,6}]
egv[j_,k_]:=Table[{i,Eigenvalues[Haml[j]/.Q->i][[k]]},{i,0,100}]
quadfit[j_,k_]:= Fit[egv[j,k], {1,x,x^2,x^3,x^4,x^5,x^6,x^7,x^8,x^9,x^10}, x]
f[x_,i_,j_]:=quadfit[i,j]
derv[lc_,adb_]:=Chop[D[f[x,lc,adb],x]/.x->0]
XX[m_,n_,t_]:=If[t >=0,Sum[derv[j,m]derv[j,n](\[Phi]re[j,t]+I \[Phi]im[j,t]),{j,1,6}],Conjugate[Sum[derv[j,m]derv[j,n](\[Phi]re[j,t]+I \[Phi]im[j,t]),{j,1,6}]]]
X[m_,n_,t_]:=If[t >=0,Sum[derv[j,m]derv[j,n](\[Phi]total[j,t]),{j,1,6}],Conjugate[Sum[derv[j,m]derv[j,n](\[Phi]total[j,Abs[t]]),{j,1,6}]]]
Y[m_,n_,t_,\[Tau]_]:=If[t > \[Tau],X[m,n,t]- X[m,n,t-\[Tau]]+Conjugate[X[m,n,\[Tau]]],X[m,n,t]- Conjugate[X[m,n,Abs[t-\[Tau]]]]+Conjugate[X[m,n,\[Tau]]]]
YY[m_,n_,t_,\[Tau]_]:=If[t > \[Tau],XX[m,n,t]- XX[m,n,t-\[Tau]]+Conjugate[XX[m,n,\[Tau]]],XX[m,n,t]- Conjugate[XX[m,n,Abs[t-\[Tau]]]]+Conjugate[XX[m,n,\[Tau]]]]
Z[m_,n_,t_]:=X[m,n,t]+Conjugate[X[m,n,t]]
ZZ[m_,n_,t_]:=XX[m,n,t]+Conjugate[XX[m,n,t]]
dec[m5_,m2_,m3_,m4_,t3_,t2_,t1_]:=Conjugate[X[m2,m2,t2]+X[m2,m2,t3]]+X[m3,m3,t3]+X[m4,m4,t2]+X[m5,m5,t1]+Y[m2,m2,t2,t3]-Y[m2,m3,t2,t3]-Y[m2,m5,t2,t1]-Y[m2,m4,t3,t2]-Y[m2,m5,t3,t1]+Y[m3,m4,t3,t2]-Z[m2,m4,t2]-Z[m2,m3,t3]+Y[m3,m5,t3,t1]+Y[m4,m5,t2,t1]


specdens[\[Omega]_,i_]:=2 \[Lambda][i]/\[Omega]c 2\[Pi] \[Omega]/(\[Omega]^2+4 \[Pi]^2)
popln[b_,c_]:=Table[{vector[[m,c]]*conjvec[[b,m]]},{m,1,6}]
coh[m_,n_,b_,c_]:={vector[[n,c]]*conjvec[[b,m]]}
jump[m_,n_]:=Sum[(vector[[m,i]]*conjvec[[i,n]])^2,{i,1,6}]
w[m_,n_]:=energy[[m]]-energy[[n]]
wT=(6.626 10^-34 /(1.3806 10^-23 T 2 N[Pi]\[Tau]c))
poplnrate[m_,n_]:=If[m==n,0,If[energy[[m]]>energy[[n]],(Sum[(vector[[m,j]]*conjvec[[j,n]])^2*2 N[Pi]specdens[w[m,n],j],{j,1,6}])/(Exp[ w[m,n]/wT]-1),(Sum[(vector[[m,j]]*conjvec[[j,n]])^2*2 N[Pi]specdens[w[n,m],j],{j,1,6}])/(1-Exp[-( w[n,m]/wT)])]]
self[n_]:=Sum[If[m==n,0,poplnrate[m,n]],{m,1,6}]
rate=Table[If[m==n,-Chop[self[n]],Chop[poplnrate[m,n]]],{m,1,6},{n,1,6}]
q[t_,b_,c_]:=MatrixExp[ rate t]. popln[b,c]
coh[t_,m_,n_,b_,c_]:=Exp[-(poplnrate[m,n]+poplnrate[n,m])/2 t]coh[m,n,b,c]
relax[t_,m_,n_,b_,c_]:=If[m==n,q[t,b,c][[m,1]],coh[t,m,n,b,c][[1]]]
qq[t_,m1_,m2_,m3_,m4_]:=Chop[Sum[vector[[m2,b]]conjvec[[c,m3]]relax[t,m1,m4,b,c],{b,1,6},{c,1,6}]]

response1[t1_,t2_,t3_]:=Sum[Chop[Chop[\[Mu]elg[3]\[Mu]vibex[4,m]\[Mu]vibex[m,4]\[Mu]elg[3]]Exp[-I (energy[[3]]-\[Epsilon]g)t1] Exp[I (energy[[4]]-energy[[4]])t2]Exp[I (energy[[4]]-energy[[m]])t3]Exp[-dec[3,4,m,4,t3,t2,t1]]qq[t2,4,3,3,4]],{m,1,1}]



\[Beta]x = 0.003125; mx = 256; \[Beta]y = 0.003125; my = 256; basx =
 Table[i \[Beta]x, {i, -mx/2, mx/2 - 1}];
basy = Table[i \[Beta]y, {i, -my/2, my/2 - 1}];

ss1 = ParallelTable[Chop[response1[basx[[k + 1]], N[1], basy[[j+ 1]]]],{j, 0, my - 1}, {k, 0, mx - 1}]



FFT2DTransformtlsR[F_, \[Beta]x_, \[Beta]y_, mx_, my_] :=
 Module[{basx, basy, a1, a2, a3,
   a4, \[Gamma]x = 2 \[Pi]/(mx \[Beta]x), \[Gamma]y =
    2 \[Pi]/(my \[Beta]y)},
  basx = Table[i \[Beta]x, {i, -mx/2, mx/2 - 1}];
  basy = Table[i \[Beta]y, {i, -my/2, my/2 - 1}];
  a1 = Table[(-1)^j F[[k + 1, j + 1]], {k, 0, my - 1}, {j, 0,
     mx - 1}]; a2 = Table[0, {j, 0, my - 1}, {i, 0, mx - 1}];
  a3 = Transpose[a2]; a4 = Transpose[a2];
  For[l = 1, l <= my, l++,
   a2[[l]] = Fourier[a1[[l]], FourierParameters -> {-1,  1}];];
  For[l = 1, l <= mx, l++,
   a3[[l]] = Table[Transpose[a2][[l, i + 1]] (-1)^i, {i, 0, my - 1}];
   a4[[l]] = Fourier[a3[[l]], FourierParameters -> { -1, 1}];];
  data = Flatten[
    Table[{(kx - mx/2) \[Gamma]x, (ky -
         my/2) \[Gamma]y, (-1)^(kx + ky) \[Beta]x \[Beta]y a4[[kx + 1,
         ky + 1]]}, {kx, 0, mx - 1}, {ky, 0, my - 1}], 1]]


FTdataR1 =
  FFT2DTransformtlsR[ss1, \[Beta]x, \[Beta]y, mx, my];

