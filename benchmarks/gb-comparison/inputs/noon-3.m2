R = ZZ/1073741827[x1,x2,x3, MonomialOrder => GRevLex];
I = ideal(
  10*x1*x2^2+10*x1*x3^2-11*x1+10,
  10*x1^2*x2+10*x2*x3^2-11*x2+10,
  10*x1^2*x3+10*x2^2*x3-11*x3+10
);
t = elapsedTiming groebnerBasis I;
G = t#1;
<< "TIME_S " << t#0 << endl;
<< "SIZE " << numColumns G << endl;
scan(first entries leadTerm G, m -> << "LM " << toString m << endl);
scan(first entries G, m -> << "POLY " << toString m << endl);
exit 0;
