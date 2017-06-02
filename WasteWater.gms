$ontext

Model for waste water treatment design problem of the MINO challenge.

Param and set definitions in *.gmsdat files, passed with --instance=*.gmsdat

Specify an output file with --out=*

Sven Wiese, Unibo, April 2016.

$offtext

$if not set out $set out results.txt

***************
* START *.dat *
***************
$include %instance%
***************
* END *.dat   *
***************

alias(NetElements,N,V)
alias(Splitters,S)
alias(Mixers,M)
alias(T,TI);

Set Range / 1*1000 /;

Set Units(N);
Units(N) = TU(N) + PU(N);

Scalars AR / 0.1 /
        alpha / 0.7 /
        gamma / 0.6 /
        pfix / 6 /
        pvar / 100 /
        pop / 0.006 /
        mill / 1000000 /
;

Variable    cost;

Positive Variables
            fin(N,T)    inflow
            fout(N,T)   outflow
            f(N,V,T)    flow between two elements
            cin(N,C,T)  inflow concentration
            cout(N,C,T) outflow concentration
            con(N,V,C,T)  concentration of flow between two nodes
            h(N,V)      helper variable for objective function
            win(N,C,T)    product for inflows
            wout(N,C,T)   product for outflows
            w(N,V,C,T)    product for flows
;

*bounds for flows when recycle == 0
f.lo(Units,Units,T)$(recycle eq 0) = 0;
f.up(Units,Units,T)$(recycle eq 0) = 0;

*bounds for sources
fout.lo(WS,T) = SourceMinFlow(WS,T);
fout.up(WS,T) = SourceMaxFlow(WS,T);

* bounds fo PUs
fin.lo(PU,T) = PUminFlow(PU);
fin.up(PU,T) = PUmaxFlow(PU);
fout.lo(PU,T) = PUminFlow(PU) + PUbalance(PU);
fout.up(PU,T) = PUmaxFlow(PU) + PUbalance(PU);
cin.up(PU,C,T) = PUcontInflow(PU,C);
cout.up(PU,C,T) = PUcontOutflow(PU,C);
* reset some to zero due to start and end times
fin.lo(PU,T)$(PUstart(PU) > ord(T) or PUend(PU) < ord(T)) = 0;
fin.up(PU,T)$(PUstart(PU) > ord(T) or PUend(PU) < ord(T)) = 0;
fout.lo(PU,T)$(PUstart(PU) > ord(T) or PUend(PU) < ord(T)) = 0;
fout.up(PU,T)$(PUstart(PU) > ord(T) or PUend(PU) < ord(T)) = 0;

*bounds for TUs
fin.lo(TU,T) = TUminFlow(TU);
fin.up(TU,T) = TUmaxFlow(TU);
fout.lo(TU,T) = TUminFlow(TU);
fout.up(TU,T) = TUmaxFlow(TU);
cin.up(TU,C,T) = TUcontin(TU,C);
cout.up(TU,C,T)$(TUtype(TU) eq 1) = TUcontin(TU,C);
cout.up(TU,C,T)$(TUtype(TU) eq 2) = min(TUcontin(TU,C),TUcontTypeChar(TU,C));

*bounds for WDs
fin.lo(WD,T) = DestMinFlow(WD,T);
fin.up(WD,T) = DestMaxFlow(WD,T);
cin.up(WD,C,T) = DestCont(WD,T,C);

*bounds for products
win.lo(M,C,T) = fin.lo(M,T) * cin.lo(M,C,T);
win.up(M,C,T) = fin.up(M,T) * cin.up(M,C,T);
wout.lo(Units,C,T) = fout.lo(Units,T) * cout.lo(Units,C,T);
wout.up(Units,C,T) = fout.up(Units,T) * cout.up(Units,C,T);

*some params for easier conditioning of equations
Parameter in(N,T);
Parameter out(N,T);
Parameter flow(N,V,T);
Loop ((M,T),
    If (fin.lo(M,T) eq 0 and fin.up(M,T) eq 0, in(M,T) = 0; else in(M,T) = 1);
) ;
Loop ((S,T),
    If (fout.lo(S,T) eq 0 and fout.up(S,T) eq 0, out(S,T) = 0; else out(S,T) = 1);
) ;
Loop ((S,M,T),
    If ((f.lo(S,M,T) eq 0 and f.up(S,M,T) eq 0) or (out(S,T) eq 0) or (in(M,T) eq 0), flow(S,M,T) = 0; else flow(S,M,T) = 1);
) ;


Binary Variables
            b(N,V)      existence of a pipe
            z(TU)       existence of a treatment unit
            xi(TU,C,T)  auxiliary binary for treatment units of type 2
            zeta(N,V,Range) helper binary variable for objective function
;

Integer Variable
            y(N,V)      pipe size
;

y.lo(N,V) = 0;
y.up(N,V) = PipeCapInt;

Equations
        obj
        SplitterFlow(N,T)       flow balance for splitters
        SplitterCon(N,V,C,T)    contaminant conservation for splitters
        MixerFlow(N,T)          flow balance for mixers
        MixerCon(N,C,T)         contaminant conservation for mixers
        PUaddFlow(N,T)          additional flow entering at PUs
        PUdisc(N,C,T)           discharge balance at PUs
        TUbalance(N,T)          flow balance at TUs
        TUconcOutType1(N,C,T)   concentration reduction in TUs of type 1
        TUconcOutType2_1(N,C,T)   concentration reduction in TUs of type 2
        TUconcOutType2_2(N,C,T)   concentration reduction in TUs of type 2
        TUconcOutType2_3(N,C,T)   concentration reduction in TUs of type 2
        TUconcOutType2_4(N,C,T)   concentration reduction in TUs of type 2
        TUconcOutType2_5(N,C,T)   concentration reduction in TUs of type 2
        PipeCapacity(N,V,T)     flow bounded by pipe capacity
        PipeNumber              maximum number of pipes
        PipeExistence(N,V)    existence of a pipe
        TUnumber                maximum number of TUs
        TUexistence(TU,T)       existence of a TU
        ConBal(C)             overall contamination balance
        objGUB1(N,V)             helper expression for objective function
        objGUB2(N,V)
        quadinflow(N,C,T)
        quadoutflow(N,C,T)
        quadflow(N,V,C,T)
;

obj..   cost =e= AR * sum(TU, TUmaxFlow(TU)**alpha * TUinvcost(TU) * z(TU)) + AR * pfix * sum((S,M)$(smax(T, flow(S,M,T)) eq 1), b(S,M)) + AR * pvar * sum((S,M)$(smax(T, flow(S,M,T)) eq 1), h(S,M)) + cycles * sum((WS,T),fout(WS,T) * WaterCost(WS,T)) + cycles * sum((TU,T)$(in(TU,T) eq 1), TUopCost(TU) * fin(TU,T)) + cycles * pop * sum((S,M,T)$(flow(S,M,T) eq 1), f(S,M,T));

SplitterFlow(S,T)$(out(S,T) eq 1).. fout(S,T) =e= sum(M$(flow(S,M,T) eq 1), f(S,M,T));

SplitterCon(Units,M,C,T)$(flow(Units,M,T) eq 1)..  con(Units,M,C,T) =e= cout(Units,C,T);

MixerFlow(M,T)$(in(M,T) eq 1)..    fin(M,T) =e= sum(S$(flow(S,M,T) eq 1),f(S,M,T));

MixerCon(M,C,T)$(in(M,T) eq 1)..   win(M,C,T) =e= sum(Units$(flow(Units,M,T) eq 1), w(Units,M,C,T)) + sum(WS$(flow(WS,M,T) eq 1), f(WS,M,T) *  SourceCont(WS,T,C));

PUaddFlow(PU,T)$(out(PU,T) eq 1 and in(PU,T) eq 1)..   fout(PU,T) =e= fin(PU,T) + PUbalance(PU);

PUdisc(PU,C,T)$(out(PU,T) eq 1 and in(PU,T) eq 1)..    wout(PU,C,T) =e= win(PU,C,T) + mill * ContDischarge(PU,C);

TUbalance(TU,T)$(in(TU,T) eq 1 and out(TU,T) eq 1)..   fin(TU,T) =e= fout(TU,T);

TUconcOutType1(TU,C,T)$((TUtype(TU) eq 1) and in(TU,T) eq 1 and out(TU,T) eq 1).. cout(TU,C,T) =e= (1 - TUcontTypeChar(TU,C)/100) * cin(TU,C,T);

TUconcOutType2_1(TU,C,T)$((TUtype(TU) eq 2) and in(TU,T) eq 1 and out(TU,T) eq 1)..    cin(TU,C,T) - TUcontTypeChar(TU,C) =l= (cin.up(TU,C,T) - TUcontTypeChar(TU,C)) * xi(TU,C,T);
TUconcOutType2_2(TU,C,T)$((TUtype(TU) eq 2) and in(TU,T) eq 1 and out(TU,T) eq 1)..    cout(TU,C,T) - cin(TU,C,T) =l= (cin.up(TU,C,T) - cin.lo(TU,C,T)) * xi(TU,C,T);
TUconcOutType2_3(TU,C,T)$((TUtype(TU) eq 2) and in(TU,T) eq 1 and out(TU,T) eq 1)..    cout(TU,C,T) - cin(TU,C,T) =g= (cin.lo(TU,C,T) - cin.up(TU,C,T)) * xi(TU,C,T);
TUconcOutType2_4(TU,C,T)$((TUtype(TU) eq 2) and in(TU,T) eq 1 and out(TU,T) eq 1)..    cout(TU,C,T) - TUcontTypeChar(TU,C) =l= (cin.up(TU,C,T) - TUcontTypeChar(TU,C)) * (1 - xi(TU,C,T));
TUconcOutType2_5(TU,C,T)$((TUtype(TU) eq 2) and in(TU,T) eq 1 and out(TU,T) eq 1)..    cout(TU,C,T) - TUcontTypeChar(TU,C) =g= (cin.lo(TU,C,T) - TUcontTypeChar(TU,C)) * (1 - xi(TU,C,T));

PipeCapacity(S,M,T)$(flow(S,M,T) eq 1)..   f(S,M,T) =l= maxPipeCap * y(S,M) / PipeCapInt;

PipeNumber..    sum((S,M)$(smax(T, flow(S,M,T)) eq 1), b(S,M)) =l= maxPipes;

PipeExistence(S,M)$(smax(T, flow(S,M,T)) eq 1)..    y(S,M) =l= PipeCapInt * b(S,M);

TUnumber..      sum(TU, z(TU)) =l= maxTU;

TUexistence(TU,T)$(in(TU,T) eq 1).. fin(TU,T) =l= fin.up(TU,T) * z(TU);

*CHECK WHAT IS THE FIRST PERIOD IN BATCH MODE ??? 0 OR 1?
ConBal(C).. sum((T,WS)$(out(WS,T) eq 1), wout(WS,C,T)) + mill * sum(PU$(PUduration(PU) eq 0), ContDischarge(PU,C) * (PUend(PU) - PUstart(PU) + 1)) =e= sum((WD,T)$(in(WD,T) eq 1), win(WD,C,T)) + sum((TU,T)$(TUduration(TU) eq 0 and in(TU,T) eq 1), win(TU,C,T) - wout(TU,C,T)) ;

objGUB1(S,M)$(smax(T,flow(S,M,T)) eq 1)..  h(S,M) =e= sum(Range$(ord(Range) - 1 <= PipeCapInt), zeta(S,M,Range) * (maxPipeCap * (ord(Range) - 1) / PipeCapInt)**gamma);

objGUB2(S,M)$(smax(T,flow(S,M,T)) eq 1).. 1 =e= sum(Range$(ord(Range) - 1 <= PipeCapInt), zeta(S,M,Range));

quadoutflow(Units,C,T)$(out(Units,T) eq 1)..    wout(Units,C,T) =e= fout(Units,T) * cout(Units,C,T);

quadinflow(M,C,T)$(in(M,T) eq 1)..    win(M,C,T) =e= fin(M,T) * cin(M,C,T);

quadflow(Units,M,C,T)$(flow(Units,M,T) eq 1)..    w(Units,M,C,T) =e= f(Units,M,T) * con(Units,M,C,T);

model wastewater / all /;

option miqcp = scip;

wastewater.reslim=3600;
wastewater.optFile=1;

solve wastewater using miqcp min cost;

