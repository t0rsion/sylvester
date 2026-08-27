R = ZZ/1073741827[x1,x2,x3,x4, MonomialOrder => GRevLex];
I = ideal(
  x1+x2+x3+x4,
  x1*x2+x2*x3+x1*x4+x3*x4,
  x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4,
  x1*x2*x3*x4-1
);
t = elapsedTiming groebnerBasis I;
G = t#1;
<< "TIME_S " << t#0 << endl;
<< "SIZE " << numColumns G << endl;
scan(first entries leadTerm G, m -> << "LM " << toString m << endl);
scan(first entries G, m -> << "POLY " << toString m << endl);
exit 0;
