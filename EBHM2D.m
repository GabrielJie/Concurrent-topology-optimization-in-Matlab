function [DH, dDH] = EBHM2D(den, lx, ly, E0, Emin, nu, penal)
% the initial definitions of the PUC
D0=E0/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];  % the elastic tensor
[nely, nelx] = size(den);
nele = nelx*nely;
dx = lx/nelx; dy = ly/nely;
Ke = elementMatVec2D(dx/2, dy/2, D0);
Num_node = (1+nely)*(1+nelx);
nodenrs = reshape(1:Num_node,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nele,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
% 3D periodic boundary formulation
alldofs = (1:2*(nely+1)*(nelx+1));
n1 = [nodenrs(end,[1,end]),nodenrs(1,[end,1])];
d1 = reshape([(2*n1-1);2*n1],1,8);
n3 = [nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
d2 = setdiff(alldofs,[d1,d3,d4]);
e0 = eye(3);
ufixed = zeros(8,3);
for j = 1:3
    ufixed(3:4,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[lx;0];
    ufixed(7:8,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;ly];
    ufixed(5:6,j) = ufixed(3:4,j)+ufixed(7:8,j);
end
wfixed = [repmat(ufixed(3:4,:),nely-1,1);repmat(ufixed(7:8,:),nelx-1,1)];
% the reduced elastic equilibrium equation to compute the induced displacement field
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK = reshape(Ke(:)*(Emin+den(:)'.^penal*(1-Emin)),64*nelx*nely,1);
K  = sparse(iK,jK,sK); K = (K + K')/2;
Kr = [K(d2,d2),K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2),K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:)= ufixed;
U([d2,d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4);K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:) = U(d3,:) + wfixed;
% homogenization to evaluate macroscopic effective properties
DH = zeros(3); qe = cell(3,3); dDH = cell(3,3);
cellVolume = lx*ly;
for i = 1:3
    for j = 1:3
        U1 = U(:,i); U2 = U(:,j);
        qe{i,j} = reshape(sum((U1(edofMat)*Ke).*U2(edofMat),2),nely,nelx)/cellVolume;
        DH(i,j) = sum(sum((Emin+den.^penal*(1-Emin)).*qe{i,j}));
        dDH{i,j} = penal*(1-Emin)*den.^(penal-1).*qe{i,j};
    end
end
disp('--- Homogenized elasticity tensor ---'); disp(DH)
end
%======================================================================================================================%
% Subfunction EBHM2D:                                                                                                  %
% A compact and efficient MATLAB code to evaluate the 2D material effective property                                   %
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