function [S] = VSSI_L1N(B,L,VertConn,sparseweight,varargin)
%% Description: Extended sources reconstruction based on sparse regularization and nuclear-norm
% Input:
%         B(d_b x T):               M/EEG Measurement
%         L(d_b x d_s):             Leadfiled Matrix
%         VertConn:                 Cortex Connectivity Condition
%         sparseweight:             sparseweight for the sparse variation
% Output:
%         S:                        Estimated Sources

[Nsensor,Nsource] = size(L);
T = size(B,2);
V = VariationEdge(VertConn);
tol = 1e-3;
QUIET = 1;
rou_update = 1;
cost = 0;
rou1 = 1e10;
rou2 = rou1/100;

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'transform'
                transform = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'lam'
                lam = varargin{i+1};
            case 'alp'
                alp = varargin{i+1};
            case 'p'
                p = varargin{i+1};
            case 'dic'
                Dic = varargin{i+1};
        end
    end
end       
  Edge = VariationEdge(VertConn);      
if strcmp(transform, 'Variation')
    V = VariationEdge(VertConn);
elseif strcmp(transform, 'Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
elseif strcmp(transform, 'Laplacian+Variation')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [opts.laplacian*V;VariationEdge(VertConn)];
elseif strcmp(transform,'Sparse+Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);V];
elseif strcmp(transform,'Sparse+Variation')
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);VariationEdge(VertConn)];
elseif strcmp(transform, 'Sparse')
    V = sparse(1:Nsource,1:Nsource,1);
end

ADMM_iter = 600;

% Initial values
TMNE = MNE(B,[],L,[],'MNE');
S_MNE = TMNE*B;
S = S_MNE;
% S = zeros(Nsource,T);

% U = V*S;  U_old = U;
% Z = S; Z_old = Z;
% S = zeros(Nsource,T);
U = V*S;  U_old = U;
Z = S;  Z_old = Z;
W = zeros(size(V,1),T);
C = zeros(size(V,2),T);
rou1_old = rou1;
rou2_old = rou2;
S_old = S;
tmp = sqrt(sum((V*L'*B).^2,2));
lam = (1/(max(tmp)*0.01))^(-1);
alp = 0.02;

% approximation of inverse (Woodbury inversion lemma)
Lambda_MAX = eigs(V'*V+rou2/rou1*speye(Nsource),1,'lm');
LLt = L*L'; LtB = L'*B;
mu = 0.9/(rou1*Lambda_MAX);
Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;

tic 
%% ADMM
alpha = 0.6;
w = 1;
for iter_ADMM = 1:ADMM_iter
%----------------------------S update----------------------------%
    S = Precision*(LtB + (1/mu)*S - rou1*V'*(V*S - U + W) - rou2*(S - Z + C));
%----------------------------U update----------------------------%
    VS = V * S;
    VS_hat = alpha*VS + (1-alpha)*U_old;
    
    U = proxl21ARD(VS_hat + W,w,lam,rou1); 
%----------------------------Z update----------------------------%
    SZ_hat = alpha*S + (1-alpha)*Z_old;

    Z = (N_prox((SZ_hat + C) * Dic,lam/alp,rou2))*Dic';
%----------------------------W update----------------------------%
    W = W + (VS_hat - U);
%----------------------------C update----------------------------%
    C = C + (SZ_hat - Z);
%-------------------------stop criterion-------------------------%
    primerror1 = norm(VS - U,'fro');
    dualerror1 = norm(rou1*V'*(U - U_old),'fro');
    U_old = U;

    primerror2 = norm(S - Z,'fro');
    dualerror2 = norm(rou2*(Z - Z_old),'fro');
    Z_old = Z;

    tol_prim1 = 1e-6*max(norm(U,'fro'),norm(VS,'fro'));
    tol_dual1 = 1e-6*rou1*norm(V'*W,'fro');
    tol_prim2 = 1e-6*max(norm(Z,'fro'),norm(S,'fro'));
    tol_dual2 = 1e-6*rou2*norm(C,'fro');

    Rprimerror1 = primerror1/max(norm(U,'fro'),norm(VS,'fro'));
    Rdualerror1 = dualerror1/norm(rou1*V'*W,'fro');
    Rprimerror2 = primerror2/max(norm(Z,'fro'),norm(S,'fro'));
    Rdualerror2 = dualerror2/norm(rou2*C,'fro');

    if ~QUIET && mod(iter_ADMM,10) == 0
        fprintf('ADMM : iteration: %g, Data-Fit : %g, PrimError: %g, DualError: %g\n', iter_ADMM, norm(B - L*S,'fro')/norm(B,'fro'), primerror1, dualerror1);
    end

    if (primerror1 < tol_prim1 && dualerror1 < tol_dual1) || (primerror2 < tol_prim2 && dualerror2 < tol_dual2)
        break;
    end


%---------------------------yita update--------------------------%
%     if max(W) / max(C) > 1e10
%         yita = yita * 10;
%     elseif max(C) / max(W) > 1e10
%         yita = yita * 0.1;
%     end
%---------------------------rou update---------------------------%
    if rou_update && mod(iter_ADMM,10) == 0
        ratio1 = -1;
        ratio2 = -1;
        if Rdualerror1~=0
            ratio1 = sqrt(Rprimerror1/Rdualerror1);
        end
        if Rdualerror2~=0
            ratio2 = sqrt(Rprimerror2/Rdualerror2);
        end

        tau_max = 2;
        if ratio1>=1 && ratio1<tau_max, tau1 = ratio1;
        elseif ratio1>1/tau_max && ratio1<1, tau1 = 1/ratio1;
        else tau1 = tau_max;
        end
        if ratio2>=1 && ratio2<tau_max, tau2 = ratio2;
        elseif ratio2>1/tau_max && ratio2<1, tau2 = 1/ratio2;
        else tau2 = tau_max;
        end

        if Rprimerror1 > 10 * Rdualerror1
            rou1 = tau1*rou1; W = W./tau1;
        elseif Rdualerror1 > 10 * Rprimerror1
            rou1 = rou1/tau1; W = W.*tau1;
        end
        if Rprimerror2 > 10 * Rdualerror2
            rou2 = tau2*rou2; C = C./tau2;
        elseif Rdualerror2 > 10 * Rprimerror2
            rou2 = rou2/tau2; C = C.*tau2;
        end
        if ~QUIET
            fprintf('rou = %g, Rprimerror = %g, Rdualerror = %g\n',rou, Rprimerror1, Rdualerror1);
        end
        if rou1 ~= rou1_old || rou2 ~= rou2_old
            Lambda_MAX = eigs(V'*V+rou2/rou1*speye(Nsource),1,'lm');
            mu = 0.9/(rou1*Lambda_MAX);
            Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;
        end
        rou1_old = rou1;
        rou2_old = rou2;
    end
    if ~mod(iter_ADMM,100)
        SS{iter_ADMM/100} = S;
    end
end
toc
S = SS;

function V = VariationEdge(VertConn)
Nsource = size(VertConn,1);
Nedge = numel(find(VertConn(:)~=0))/2;
V = sparse(Nedge,Nsource);
edge = 0;
for i = 1:Nsource
    idx = find(VertConn(i,:)~=0);
    idx = idx(idx<i);
    for j = 1:numel(idx)
        V(j+edge,i) = 1;
        V(j+edge,idx(j)) = -1;
    end
    edge = edge + numel(idx);
end
end


function A = N_prox(Y,lam,rou)
[Kkk, Sss, Vvv] = svd(Y, 'econ');
thr = lam/rou;
Sss_trunc = Sss;
Sss_trunc(Sss_trunc<thr) = 0;
ind = find(Sss_trunc > 0);
Sss_trunc(ind) = Sss(ind) - thr;
A = Kkk * Sss_trunc * Vvv';
end

% function x = L21_prox(y,G,lam,rou)
% temp = lam/rou;
% numcolumns = size(y,2);
% x = zeros(size(y));
% w_s = sum(G.^2,1);
% for t = 1:numcolumns
%     threshold = temp * sqrt(w_s(t));
%     x(:,t) = y(:,t).*(1-min(1,threshold./norm(y(:,t),2))).^max(0,1);
%      o = max(threshold./norm(y(:,t),2));
%      i = min(threshold./norm(y(:,t),2));
% end
% end

function Z = proxl21ARD(Y,w,lam,rou)
% Z = arg min_Y 0.5*rou*|| Y-Z ||_2^2 + lam * (sum_k sqrt(sum_t w_k * Z(k,t)^2))
[m,n] = size(Y);
temp = lam*sqrt(w)./sqrt(sum(Y.^2,2))/rou;
scale = pos(ones(m,1) - temp);
Z = Y.*repmat(scale,1,n);
end

% function x = pos(y)
% x = max(0, y);
% end

end


    
