function R = ReactMat2D(p,t)

np = size(p,2); 
nt = size(t,2);
R = sparse(np,np);

for K = 1:nt 
    loc2glb = t(1:3,K); 
    x = p(1,loc2glb);
    y = p(2,loc2glb); 
    area = polyarea(x,y); 
    RK = [2 1 1;
          1 2 1;
          1 1 2]/12*area; 
    R(loc2glb,loc2glb) = R(loc2glb,loc2glb)+ RK; 
end

end