% This function is used to compute the best spectacle correction 

[sphere cylinder axis] = spectacleOSA(c);


cm3=c(3)/1e6;
cm4=c(4)/1e6;
cm5=(c(5)/1e6);

if (cm5==0)
   phi=(-1)*sign(cm3)*pi/4;% special case when cm5 is equal to zero
else
   phi = (-1)*0.5*atan(cm3/cm5);
end

if (abs(cm5)<abs(cm3))
   Acoeff=cm3*2*sqrt(6)/sin(2*phi);
else
   Acoeff=-cm5*2*sqrt(6)/cos(2*phi);
end

axis = (-1)*180*phi/pi;

if (Acoeff<=0)
   Acoeff=(-1)*Acoeff;
   axis=axis-90;
end

if (axis<=0)
   axis=axis+180;
end

Dcoeff=cm4*2*sqrt(3)-Acoeff/2;

cylinder = 2*Acoeff/((PARAMS(3)/2000)^2);
defocus = 2*Dcoeff/((PARAMS(3)/2000)^2);
if (abs(cylinder)<0.01)
   cylinder =0;axis=0;
end
if (abs(defocus)<0.01)
   defocus=0;
end



