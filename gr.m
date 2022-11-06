function g = gr(t,a,C,epsillon)

g = C-a*tanh(a*(1-C*t)/(2*epsillon));

end