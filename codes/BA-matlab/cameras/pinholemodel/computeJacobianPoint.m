function grad_pt = computeJacobianPoint(point, cameraparams)
%COMPUTEJACOBIANPOINT
%    GRAD_PT = COMPUTEJACOBIANPOINT(X,Y,Z,FX,FY,TX,TY,TZ,WX,WY,WZ,X0,Y0)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    28-Nov-2018 13:44:12

X = point(1);
Y = point(2);
Z = point(3);
fx = cameraparams(1);
fy = cameraparams(2);
x0 = cameraparams(3);
y0 = cameraparams(4);
wx = cameraparams(5);
wy = cameraparams(6);
wz = cameraparams(7);
tx = cameraparams(8);
ty = cameraparams(9);
tz = cameraparams(10);

if (0.0 < wx.^2+wy.^2+wz.^2)
    t0 = (fx.*((wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-x0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fx.*tx-Y.*(fx.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-x0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*x0+X.*(fx.*((wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-x0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(x0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+fx.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t0 = fx./(Z+tz);
end
if (0.0 < wx.^2+wy.^2+wz.^2)
    t1 = (fy.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-y0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fy.*ty+X.*(fy.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-y0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*y0+Y.*(fy.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+y0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(y0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-fy.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t1 = 0.0;
end
if (0.0 < wx.^2+wy.^2+wz.^2)
    t2 = -(fx.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-x0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))-(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fx.*tx-Y.*(fx.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-x0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*x0+X.*(fx.*((wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-x0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(x0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+fx.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t2 = 0.0;
end
if (0.0 < wx.^2+wy.^2+wz.^2)
    t3 = (fy.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+y0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))-(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fy.*ty+X.*(fy.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-y0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*y0+Y.*(fy.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+y0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(y0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-fy.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t3 = fy./(Z+tz);
end
if (0.0 < wx.^2+wy.^2+wz.^2)
    t4 = (x0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+fx.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))-((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fx.*tx-Y.*(fx.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-x0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*x0+X.*(fx.*((wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-x0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(x0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+fx.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t4 = -1.0./(Z+tz).^2.*(X.*fx+Z.*x0+fx.*tx+tz.*x0)+x0./(Z+tz);
end
if (0.0 < wx.^2+wy.^2+wz.^2)
    t5 = (y0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-fy.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))-((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0).*1.0./(tz+Z.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-X.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))+Y.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))).^2.*(fy.*ty+X.*(fy.*(wz.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wx.*wy.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))-y0.*(wy.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wx.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+tz.*y0+Y.*(fy.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wz.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)+y0.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)-(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)))+Z.*(y0.*((wx.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+(wy.^2.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2)+1.0)-fy.*(wx.*sin(sqrt(wx.^2+wy.^2+wz.^2)).*1.0./sqrt(wx.^2+wy.^2+wz.^2)+(wy.*wz.*(cos(sqrt(wx.^2+wy.^2+wz.^2))-1.0))./(wx.^2+wy.^2+wz.^2))));
else
    t5 = -1.0./(Z+tz).^2.*(Y.*fy+Z.*y0+fy.*ty+tz.*y0)+y0./(Z+tz);
end
grad_pt = reshape([t0,t1,t2,t3,t4,t5],[2,3]);
