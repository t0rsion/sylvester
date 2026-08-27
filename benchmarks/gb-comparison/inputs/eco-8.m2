R = ZZ/1073741827[x1,x2,x3,x4,x5,x6,x7,x8, MonomialOrder => GRevLex];
I = ideal(
  x1*x2*x8+x2*x3*x8+x3*x4*x8+x4*x5*x8+x5*x6*x8+x6*x7*x8+x1*x8-1,
  x1*x3*x8+x2*x4*x8+x3*x5*x8+x4*x6*x8+x5*x7*x8+x2*x8-2,
  x1*x4*x8+x2*x5*x8+x3*x6*x8+x4*x7*x8+x3*x8-3,
  x1*x5*x8+x2*x6*x8+x3*x7*x8+x4*x8-4,
  x1*x6*x8+x2*x7*x8+x5*x8-5,
  x1*x7*x8+x6*x8-6,
  x7*x8-7,
  x1+x2+x3+x4+x5+x6+x7+1
);
t = elapsedTiming groebnerBasis I;
G = t#1;
<< "TIME_S " << t#0 << endl;
<< "SIZE " << numColumns G << endl;
scan(first entries leadTerm G, m -> << "LM " << toString m << endl);
scan(first entries G, m -> << "POLY " << toString m << endl);
exit 0;
