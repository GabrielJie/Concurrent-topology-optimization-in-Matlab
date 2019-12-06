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
vert_cor = [0, lx, lx,  0,  0, lx, lx,  0;
            0,  0, ly, ly,  0,  0, ly, ly;
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
%======================================================================================================================%
% Subfunction EBHM3D:                                                                                                  %
% A compact and efficient MATLAB code to evaluate the 3D material effective property                                   %
%                                                                                                                      %
% Developed by: Jie Gao, Zhen Luo, Liang Xia and Liang Gao*                                                            %
% Email: gaoliang@mail.hust.edu.cn (GabrielJie_Tian@163.com)                                                           %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Jie Gao, Zhen Luo, Liang Xia, Liang Gao. Concurrent topology optimization of multiscale composite structures     %
% in Matlab. Accepted in Structural and multidisciplinary optimization.                                                %
%                                                                                                                      %
% (2) Gao J, Li H, Gao L, et al. Topological shape optimization of 3D micro-structured materials using energy-based    %
% homogenization method. Advances in Engineering Software, 2018, 116: 89-102.                                          %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%