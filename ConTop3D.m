function ConTop3D(Macro_struct, Micro_struct, penal, rmin)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200; displayflag = 1; E0 = 1; Emin = 1e-9; nu = 0.3;
Macro.length = Macro_struct(1); Macro.width = Macro_struct(2); Macro.Height = Macro_struct(3);
Micro.length = Micro_struct(1); Micro.width = Micro_struct(2); Micro.Height = Micro_struct(3);
Macro.nelx   = Macro_struct(4); Macro.nely  = Macro_struct(5); Macro.nelz   = Macro_struct(6);
Micro.nelx   = Micro_struct(4); Micro.nely  = Micro_struct(5); Micro.nelz   = Micro_struct(6);
Macro.Vol = Macro_struct(7); Micro.Vol = Micro_struct(7);
Macro.Elex   = Macro.length/Macro.nelx; Macro.Eley = Macro.width/Macro.nely; Macro.Elez = Macro.Height/Macro.nelz;
Macro.nele = Macro.nelx*Macro.nely*Macro.nelz; Micro.nele = Micro.nelx*Micro.nely*Micro.nelz;
Macro.ndof = 3*(Macro.nelx+1)*(Macro.nely+1)*(Macro.nelz+1);
% PREPARE FINITE ELEMENT ANALYSIS
[load_x,load_y,load_z] = meshgrid(Macro.nelx/2, Macro.nely, Macro.nelz/2);
loadnid = load_z*(Macro.nelx+1)*(Macro.nely+1)+load_x*(Macro.nely+1)+(Macro.nely+1-load_y);
F = sparse(3*loadnid(:) - 1,1,-1,Macro.ndof,1);
U = zeros(Macro.ndof,1);
[fixed_x,fixed_y,fixed_z] = meshgrid([0 Macro.nelx],0,[0 Macro.nelz]);
fixednid = fixed_z*(Macro.nelx+1)*(Macro.nely+1)+fixed_x*(Macro.nely+1)+(Macro.nely+1-fixed_y);
fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
freedofs = setdiff(1:Macro.ndof,fixeddofs);
nodegrd = reshape(1:(Macro.nely+1)*(Macro.nelx+1),Macro.nely+1,Macro.nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),Macro.nely*Macro.nelx,1);
nodeidz = 0:(Macro.nely+1)*(Macro.nelx+1):(Macro.nelz-1)*(Macro.nely+1)*(Macro.nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofMat = repmat(3*nodeids(:)+1,1,24)+ repmat([0 1 2 3*Macro.nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(Macro.nely+1)*(Macro.nelx+1)+[0 1 2 3*Macro.nely + [3 4 5 0 1 2] -3 -2 -1]],Macro.nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*Macro.nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*Macro.nele,1);
% PREPARE FILTER
[Macro.H,Macro.Hs] = filtering3d(Macro.nelx, Macro.nely, Macro.nelz, Macro.nele, rmin);
[Micro.H,Micro.Hs] = filtering3d(Micro.nelx, Micro.nely, Micro.nelz, Micro.nele, rmin);
% INITIALIZE ITERATION
Macro.x = repmat(Macro.Vol,[Macro.nely,Macro.nelx,Macro.nelz]);
Micro.x = ones(Micro.nely,Micro.nelx,Micro.nelz);
Micro.x(Micro.nely/2:Micro.nely/2+1,Micro.nelx/2:Micro.nelx/2+1,Micro.nelz/2:Micro.nelz/2+1) = 0;
beta = 1;
Macro.xTilde = Macro.x; Micro.xTilde = Micro.x;
Macro.xPhys = 1-exp(-beta*Macro.xTilde)+Macro.xTilde*exp(-beta);
Micro.xPhys = 1-exp(-beta*Micro.xTilde)+Micro.xTilde*exp(-beta);
loopbeta = 0; loop = 0; Macro.change = 1; Micro.change = 1;
while loop < maxloop || Macro.change > 0.01 || Micro.change > 0.01
    loop = loop+1; loopbeta = loopbeta+1;
    % FE-ANALYSIS AT TWO SCALES
    [DH, dDH] = EBHM3D(Micro.xPhys, Micro.length, Micro.width, Micro.Height, E0, Emin, nu, penal);
    Ke = elementMatVec3D(Macro.Elex/2, Macro.Eley/2, Macro.Elez/2, DH);
    sK = reshape(Ke(:)*(Emin+Macro.xPhys(:)'.^penal*(1-Emin)),24*24*Macro.nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*Ke).*U(edofMat),2),[Macro.nely,Macro.nelx,Macro.nelz]);
    c = sum(sum(sum((Emin+Macro.xPhys.^penal*(1-Emin)).*ce)));
    Macro.dc = -penal*(1-Emin)*Macro.xPhys.^(penal-1).*ce;
    Macro.dv = ones(Macro.nely, Macro.nelx, Macro.nelz);
    Micro.dc = zeros(Micro.nely, Micro.nelx, Micro.nelz);
    for i = 1:Micro.nele
        dDHe = [dDH{1,1}(i) dDH{1,2}(i) dDH{1,3}(i) dDH{1,4}(i) dDH{1,5}(i) dDH{1,6}(i);
                dDH{2,1}(i) dDH{2,2}(i) dDH{2,3}(i) dDH{2,4}(i) dDH{2,5}(i) dDH{2,6}(i);
                dDH{3,1}(i) dDH{3,2}(i) dDH{3,3}(i) dDH{3,4}(i) dDH{3,5}(i) dDH{3,6}(i);
                dDH{4,1}(i) dDH{4,2}(i) dDH{4,3}(i) dDH{4,4}(i) dDH{4,5}(i) dDH{4,6}(i);
                dDH{5,1}(i) dDH{5,2}(i) dDH{5,3}(i) dDH{5,4}(i) dDH{5,5}(i) dDH{5,6}(i);
                dDH{6,1}(i) dDH{6,2}(i) dDH{6,3}(i) dDH{6,4}(i) dDH{6,5}(i) dDH{6,6}(i)];
        [dKE] = elementMatVec3D(Macro.Elex, Macro.Eley, Macro.Elez, dDHe);
        dce = reshape(sum((U(edofMat)*dKE).*U(edofMat),2),[Macro.nely,Macro.nelx,Macro.nelz]);
        Micro.dc(i) = -sum(sum(sum((Emin+Macro.xPhys.^penal*(1-Emin)).*dce)));
    end
    Micro.dv = ones(Micro.nely, Micro.nelx, Micro.nelz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    Macro.dx = beta*exp(-beta*Macro.xTilde)+exp(-beta); Micro.dx = beta*exp(-beta*Micro.xTilde)+exp(-beta);
    Macro.dc(:) = Macro.H*(Macro.dc(:).*Macro.dx(:)./Macro.Hs); Macro.dv(:) = Macro.H*(Macro.dv(:).*Macro.dx(:)./Macro.Hs);
    Micro.dc(:) = Micro.H*(Micro.dc(:).*Micro.dx(:)./Micro.Hs); Micro.dv(:) = Micro.H*(Micro.dv(:).*Micro.dx(:)./Micro.Hs);
    % OPTIMALITY CRITERIA UPDATE MACRO AND MICRO ELELMENT DENSITIES
    [Macro.x, Macro.xPhys, Macro.change] = OC(Macro.x, Macro.dc, Macro.dv, Macro.H, Macro.Hs, Macro.Vol, Macro.nele, 0.02, beta);
    [Micro.x, Micro.xPhys, Micro.change] = OC(Micro.x, Micro.dc, Micro.dv, Micro.H, Micro.Hs, Micro.Vol, Micro.nele, 0.02, beta);
    Macro.xPhys = reshape(Macro.xPhys, Macro.nely, Macro.nelx, Macro.nelz);
    Micro.xPhys = reshape(Micro.xPhys, Micro.nely, Micro.nelx, Micro.nelz);
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Macro_Vol.:%7.3f Micro_Vol.:%7.3f Macro_ch.:%7.3f Micro_ch.:%7.3f\n',...
        loop,c,mean(Macro.xPhys(:)),mean(Micro.xPhys(:)), Macro.change, Micro.change);
    if displayflag, clf; display_3D(Macro.xPhys); end
    if displayflag, clf; display_3D(Micro.xPhys); end
    % UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if beta < 512 && (loopbeta >= 50 || Macro.change <= 0.01 || Micro.change <= 0.01 )
        beta = 2*beta; loopbeta = 0; Macro.change = 1; Micro.change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
end
%% SUB FUNCTION:filtering3D
function [H,Hs] = filtering3d(nelx, nely, nelz, nele, rmin)
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);
end
%% SUB FUNCTION: EBHM3D
function [DH, dDH] = EBHM3D(den, lx, ly, lz, E0, Emin, nu, penal)
% the initial definitions of the PUC
D0 = E0/(1+nu)/(1-2*nu)*...
    [ 1-nu   nu   nu     0          0          0     ;
        nu 1-nu   nu     0          0          0     ;
        nu   nu 1-nu     0          0          0     ;
         0    0    0 (1-2*nu)/2     0          0     ;
         0    0    0     0      (1-2*nu)/2     0     ;
         0    0    0     0          0      (1-2*nu)/2];
[nely, nelx, nelz] = size(den); nele = nelx*nely*nelz;
dx = lx/nelx; dy = ly/nely; dz = lz/nelz;
Ke = elementMatVec3D(dx/2, dy/2, dz/2, D0);
Num_node = (1+nely)*(1+nelx)*(1+nelz);
nodenrs = reshape(1:Num_node,1+nely,1+nelx,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nelx*nely*nelz,1);
edofMat = repmat(edofVec,1,24)+repmat([0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 ...
    3*(nelx+1)*(nely+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]], nele, 1);
% 3D periodic boundary formulation
% the nodes classification
n1 = [nodenrs(end, [1 end], 1) nodenrs(1, [end 1], 1) nodenrs(end, [1 end], end) nodenrs(1, [end 1], end)];
d1 = reshape([3*n1-2; 3*n1-1; 3*n1],3*numel(n1),1);
n3 = [reshape(squeeze(nodenrs(end,1,2:end-1)),1,numel(squeeze(nodenrs(end,1,2:end-1))))...              % AE
      reshape(squeeze(nodenrs(1, 1, 2:end-1)),1,numel(squeeze(nodenrs(1, 1, 2:end-1))))...              % DH
      reshape(squeeze(nodenrs(end,2:end-1,1)),1,numel(squeeze(nodenrs(end,2:end-1,1))))...              % AB
      reshape(squeeze(nodenrs(1, 2:end-1, 1)),1,numel(squeeze(nodenrs(1, 2:end-1, 1))))...              % DC
      reshape(squeeze(nodenrs(2:end-1, 1, 1)),1,numel(squeeze(nodenrs(2:end-1, 1, 1))))...              % AD
      reshape(squeeze(nodenrs(2:end-1,1,end)),1,numel(squeeze(nodenrs(2:end-1,1,end))))...              % EH
      reshape(squeeze(nodenrs(2:end-1, 2:end-1, 1)),1,numel(squeeze(nodenrs(2:end-1, 2:end-1, 1))))...  % ABCD
      reshape(squeeze(nodenrs(2:end-1, 1, 2:end-1)),1,numel(squeeze(nodenrs(2:end-1, 1, 2:end-1))))...  % ADHE
      reshape(squeeze(nodenrs(end,2:end-1,2:end-1)),1,numel(squeeze(nodenrs(end,2:end-1,2:end-1))))];   % ABFE                   
d3 = reshape([3*n3-2; 3*n3-1; 3*n3],3*numel(n3),1);
n4 = [reshape(squeeze(nodenrs(1, end, 2:end-1)),1,numel(squeeze(nodenrs(1, end, 2:end-1))))...          % CG
      reshape(squeeze(nodenrs(end,end,2:end-1)),1,numel(squeeze(nodenrs(end,end,2:end-1))))...          % BF
      reshape(squeeze(nodenrs(1, 2:end-1, end)),1,numel(squeeze(nodenrs(1, 2:end-1, end))))...          % HG
      reshape(squeeze(nodenrs(end,2:end-1,end)),1,numel(squeeze(nodenrs(end,2:end-1,end))))...          % EF
      reshape(squeeze(nodenrs(2:end-1,end,end)),1,numel(squeeze(nodenrs(2:end-1,end,end))))...          % FG
      reshape(squeeze(nodenrs(2:end-1, end, 1)),1,numel(squeeze(nodenrs(2:end-1, end, 1))))...          % BC
      reshape(squeeze(nodenrs(2:end-1,2:end-1,end)),1,numel(squeeze(nodenrs(2:end-1,2:end-1,end))))...  % EFGH
      reshape(squeeze(nodenrs(2:end-1,end,2:end-1)),1,numel(squeeze(nodenrs(2:end-1,end,2:end-1))))...  % BCGF
      reshape(squeeze(nodenrs(1, 2:end-1, 2:end-1)),1,numel(squeeze(nodenrs(1, 2:end-1, 2:end-1))))];   % DCGH
d4 = reshape([3*n4-2; 3*n4-1; 3*n4],3*numel(n4),1);
n2 = setdiff(nodenrs(:),[n1(:);n3(:);n4(:)]); d2 = reshape([3*n2-2; 3*n2-1; 3*n2],3*numel(n2),1);
% the imposing of six linearly independent unit test strains
e = eye(6); ufixed = zeros(24,6);
vert_cor = [0  lx, lx,  0,  0, lx, lx,  0;
            0   0, ly, ly,  0,  0, ly, ly;
            0,  0,  0,  0, lz, lz, lz, lz];
for i = 1:6
    epsilon = [  e(i,1), e(i,4)/2, e(i,6)/2;
               e(i,4)/2,   e(i,2), e(i,5)/2;
               e(i,6)/2, e(i,5)/2,   e(i,3)];
    ufixed(:,i) = reshape(epsilon*vert_cor,24,1);
end
% 3D boundary constraint equations
wfixed = [repmat(ufixed(  7:9,:),numel(squeeze(nodenrs(end,1,2:end-1))),1);                    % C
          repmat(ufixed(  4:6,:)-ufixed(10:12,:),numel(squeeze(nodenrs(1, 1, 2:end-1))),1);    % B-D
          repmat(ufixed(22:24,:),numel(squeeze(nodenrs(end,2:end-1,1))),1);                    % H
          repmat(ufixed(13:15,:)-ufixed(10:12,:),numel(squeeze(nodenrs(1, 2:end-1, 1))),1);    % E-D
          repmat(ufixed(16:18,:),numel(squeeze(nodenrs(2:end-1, 1, 1))),1);                    % F
          repmat(ufixed(  4:6,:)-ufixed(13:15,:),numel(squeeze(nodenrs(2:end-1,1,end))),1);    % B-E
          repmat(ufixed(13:15,:),numel(squeeze(nodenrs(2:end-1, 2:end-1, 1))),1);              % E
          repmat(ufixed(  4:6,:),numel(squeeze(nodenrs(2:end-1, 1, 2:end-1))),1);              % B
          repmat(ufixed(10:12,:),numel(squeeze(nodenrs(end,2:end-1,2:end-1))),1)];             % D
% the reduced elastic equilibrium equation to compute the induced displacement field
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
sK = reshape(Ke(:)*(Emin+den(:)'.^penal*(1-Emin)),24*24*nele,1);
K = sparse(iK(:), jK(:), sK(:)); K = (K+K')/2;
Kr = [K(d2,d2), K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2),K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:)= ufixed;
U([d2;d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4);K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:) = U(d3,:) + wfixed;
% homogenization to evaluate macroscopic effective properties
qe = cell(6,6);DH = zeros(6,6); dDH = cell(6,6);
cellVolume = lx*ly*lz;
for i = 1:6
    for j = 1:6
        U1 = U(:,i); U2 = U(:,j);
        qe{i,j} = reshape(sum((U1(edofMat)*Ke).*U2(edofMat),2),nely,nelx,nelz);
        DH(i,j) = 1/cellVolume*sum(sum(sum((Emin+den.^penal*(1-Emin)).*qe{i,j})));
        dDH{i,j} = 1/cellVolume*(penal*(1-Emin)*den.^(penal-1).*qe{i,j});
    end
end
disp('--- Homogenized elasticity tensor ---'); disp(DH)
end
%% SUB FUNCTION: elementMatVec3D
function Ke = elementMatVec3D(a, b, c, DH)
GN_x=[-1/sqrt(3),1/sqrt(3)]; GN_y=GN_x; GN_z=GN_x; GaussWeigh=[1,1];
Ke = zeros(24,24); L = zeros(6,9);
L(1,1) = 1; L(2,5) = 1; L(3,9) = 1;
L(4,2) = 1; L(4,4) = 1; L(5,6) = 1;
L(5,8) = 1; L(6,3) = 1; L(6,7) = 1;
for ii=1:length(GN_x)
    for jj=1:length(GN_y)
        for kk=1:length(GN_z)
            x = GN_x(ii);y = GN_y(jj);z = GN_z(kk);
            dNx = 1/8*[-(1-y)*(1-z)  (1-y)*(1-z)  (1+y)*(1-z) -(1+y)*(1-z) -(1-y)*(1+z)  (1-y)*(1+z)  (1+y)*(1+z) -(1+y)*(1+z)];
            dNy = 1/8*[-(1-x)*(1-z) -(1+x)*(1-z)  (1+x)*(1-z)  (1-x)*(1-z) -(1-x)*(1+z) -(1+x)*(1+z)  (1+x)*(1+z)  (1-x)*(1+z)];
            dNz = 1/8*[-(1-x)*(1-y) -(1+x)*(1-y) -(1+x)*(1+y) -(1-x)*(1+y)  (1-x)*(1-y)  (1+x)*(1-y)  (1+x)*(1+y)  (1-x)*(1+y)];
            J = [dNx;dNy;dNz]*[ -a  a  a  -a  -a  a  a  -a ;  -b  -b  b  b  -b  -b  b  b; -c -c -c -c  c  c  c  c]';
            G = [inv(J) zeros(3) zeros(3);zeros(3) inv(J) zeros(3);zeros(3) zeros(3) inv(J)];
            dN(1,1:3:24) = dNx; dN(2,1:3:24) = dNy; dN(3,1:3:24) = dNz;
            dN(4,2:3:24) = dNx; dN(5,2:3:24) = dNy; dN(6,2:3:24) = dNz;
            dN(7,3:3:24) = dNx; dN(8,3:3:24) = dNy; dN(9,3:3:24) = dNz;
            Be = L*G*dN;
            Ke = Ke + GaussWeigh(ii)*GaussWeigh(jj)*GaussWeigh(kk)*det(J)*(Be'*DH*Be);
        end
    end
end
end
%% SUB FUNCTION: OC
function [x, xPhys, change] = OC(x, dc, dv, H, Hs, volfrac, nele, move, beta)
l1 = 0; l2 = 1e9;
while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    xTilde(:) = (H*xnew(:))./Hs; xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else, l2 = lmid; end
end
change = max(abs(xnew(:)-x(:))); x = xnew;
end
%% SUB FUNCTION: display_3D
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0)
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
%======================================================================================================================%
% Function ConTop3D:                                                                                                   %
% A compact and efficient MATLAB code for Concurrent topology optimization of multiscale composite structures          %
% in Matlab.                                                                                                           %
%                                                                                                                      %
% Developed by: Jie Gao, Zhen Luo, Liang Xia and Liang Gao*                                                            %
% Email: gaoliang@mail.hust.edu.cn (GabrielJie_Tian@163.com)                                                           %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Jie Gao, Zhen Luo, Liang Xia, Liang Gao. Concurrent topology optimization of multiscale composite structures     %
% in Matlab. Accepted in Structural and multidisciplinary optimization.                                                %
%                                                                                                                      %
% (2) Xia L, Breitkopf P. Design of materials using topology optimization and energy-based homogenization approach in  %
% Matlab. % Structural and multidisciplinary optimization, 2015, 52(6): 1229-1241.                                     %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%
