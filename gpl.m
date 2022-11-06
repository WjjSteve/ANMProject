function g = gpl(t,a,C,epsillon)

g = a^2*C/(2*epsillon)*(1-tanh(a*(-1-C*t)/(2*epsillon))^2);

end