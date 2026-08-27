R = ZZ/1073741827[x1,x2,x3,x4,x5, MonomialOrder => GRevLex];
I = ideal(
  10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2-11*x1+10,
  10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2-11*x2+10,
  10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2-11*x3+10,
  10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2-11*x4+10,
  10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5-11*x5+10
);
t = elapsedTiming groebnerBasis I;
G = t#1;
<< "TIME_S " << t#0 << endl;
<< "SIZE " << numColumns G << endl;
scan(first entries leadTerm G, m -> << "LM " << toString m << endl);
scan(first entries G, m -> << "POLY " << toString m << endl);
exit 0;
