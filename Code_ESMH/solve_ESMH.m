function [F,W1_final, W2_final,B_final,G_final,mu1, mu2] = solve_ESMH(phi_I, phi_T, S,Y, param)
%% random initialization
[col, row] = size(phi_I);
[~, rowt] = size(phi_T);
[Nlab,~] = size(Y);
threshold = 0.1;
lastF = 100000000;
iter = 1;

alpha = param.alpha;
beta = param.beta;
lambda = param.lambda;
eta = param.eta;
bits = param.bits;

phi_I1 = phi_I(1:Nlab,:);
phi_I2 = phi_I(Nlab+1:end,:);
phi_T1 = phi_T(1:Nlab,:);
phi_T2 = phi_T(Nlab+1:end,:);

SL = param.SL;
SU = param.SU;
mm = SL*ones(Nlab,1);
E1 = diag(mm);
nn = SU*ones(col-Nlab,1);
E2 = diag(nn);
mn = [mm;nn];
E = diag(mn);
mu1 = 0.8;
mu2 = 0.2;
%% Initialize
B = randn(col,bits)>0;B=B*2-1;
B1 = B(1:Nlab,:);
Z = randn(Nlab,bits)>0;Z=Z*2-1;
H = B1-Z;
W2 = randn(rowt,bits);
A = S*Y;

%% Iterative algorithm
while (iter<50)    
    % Update G
    G = pinv(Y'*Y)*(alpha*A'*B1+beta*Y'*B1)/(alpha*(B1'*B1)+beta*eye(bits));
    
    % Update Z H
    Z = sign(-alpha*B1*G'*(Y'*Y)*G+eta*B1+H);
    H = H+eta*(B1-Z); 
    
    % update W1
    W1 = (mu1.*mu1*((E*phi_I)'*(E*phi_I))+lambda*eye(row))\(mu1*(E*phi_I)'*E*B-mu1.*mu2*(E*phi_I)'*(E*phi_T)*W2);
    W2 = (mu2.*mu2*((E*phi_T)'*(E*phi_T))+lambda*eye(rowt))\(mu2*(E*phi_T)'*E*B-mu1.*mu2*(E*phi_T)'*(E*phi_I)*W1);
    
    % update B
    P1 = E1'*E1*mu1*phi_I1*W1+E1'*E1*mu2*phi_T1*W2;
    P2 = E2'*E2*mu1*phi_I2*W1+E2'*E2*mu2*phi_T2*W2;
    P = [P1;P2];       
    B1 = sign(2*P1+2*alpha*A*G-alpha*Z*G'*(Y'*Y)*G+2*beta*Y*G+eta*Z-H);   
    B2 = sign(2*P2);
    B = [B1;B2];
    
    %%  compute objective function
    norm1 = sum(sum((B-P).^2));
    norm2 = alpha*sum(sum((S - B1*(Y*G)').^2));
    norm3 = beta*sum(sum((B1-Y*G).^2));
    norm4 = lambda*(sum(sum((W1).^2))+sum(sum((W2).^2)));
    currentF = norm1 + norm2+norm3+norm4;
    F(iter) = currentF;
    W1_final = W1;
    W2_final = W2;
    B_final = B;
    G_final = G;
    
%     fprintf('\ncurrentF at iteration %d: %.2f; obj: %.4f\n',iter, currentF, lastF - currentF);
    if (lastF - currentF) < threshold
        if(iter > 3)
            return;
        end
    end
    iter = iter + 1;
    lastF = currentF;
end
return;

