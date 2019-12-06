function [x, xPhys, change] = OC(x, dc, dv, H, Hs, volfrac, nele, ft, move, beta)
l1 = 0; l2 = 1e9;
while (l2-l1)/(l1+l2) > 1e-4
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    elseif ft == 3
        xTilde(:) = (H*xnew(:))./Hs;
        xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    end
    if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
end
change = max(abs(xnew(:)-x(:)));
x = xnew;
end