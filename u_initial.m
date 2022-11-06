function u0 = u_initial(x1,x2)

if ((x1)^2+(x2)^2<=1)
    u0 = 14/4/pi;
else
    u0 = pi/4;
end

end