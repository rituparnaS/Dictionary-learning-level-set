function [phi,p1_in,p2_out,c1_vec,c2_vec,cn1,cn2] = ChanVeseDL(opt)
%CHANVESE Segmentation by Chan-Vese's method

Dirac_global            = @(x,e) ((e/pi)./(e^2.+ x.^2));
Heaviside               = @(y,e) (0.5*(1+(2/pi)*atan(y/e)));

its         = 0;
max_iter    = opt.num_iter;
u           = opt.Img;
phi         = opt.init_phi;
convg_err   = opt.convg_error;
reg_term    = opt.length_term;
l1          = opt.lambda1;
l2          = opt.lambda2;
display     = opt.display_intrvl;
mag         = opt.img_magnify;
color       = opt.contour_color;
count_lim   = opt.convg_count;
B           = opt.basis_vect;             % vectorized basis, each col is the basis vector
lambda_l2   = opt.lambda_l2;              % constant for L2 regularization 

stop = 0;
count = 0;
figure;

II = eye(size(B,2),size(B,2));
cn1 = [];
cn2=[];
while (its < max_iter && stop == 0)
    
    h_phi = Heaviside(phi,2);
    inside_mask = h_phi;
    outside_mask = 1-h_phi;
    u_in = u.*inside_mask;
    u_out = u.*outside_mask;
    
    u_in = u_in(:);
    u_out = u_out(:);
    
    inside_indicator = inside_mask(:);
    outside_indicator = outside_mask(:);
    
    A1 = B';    % ( each row contains a basis vector)
    A2 = A1.*(repmat(inside_indicator',size(A1,1),1));   % A1, with each row multiplied with hphi
    B2 = A1.*(repmat(outside_indicator',size(A1,1),1));   % A1, with each row multiplied with hphi
    
    cn1 = [cn1;cond(A1*A2')];
    cn2 = [cn2;cond(A1*B2')];
    
    c1_vec = (A1*A2' + lambda_l2*II)\(A1*u_in);
    c2_vec = (A1*B2' + lambda_l2*II)\(A1*u_out);
    
    p1_vec = B*c1_vec;
    p2_vec = B*c2_vec;

    p1 = reshape(p1_vec,size(u));
    p2 = reshape(p2_vec,size(u));

    p1_in = p1.*h_phi;
    p2_out = p2.*(1-h_phi);
    
    curvature   = curvature_central(phi);
    delta_phi   = Dirac_global(phi,2);
    
    evolve_force = delta_phi.*(-l1*(u-p1).^2 + l2*(u-p2).^2);
    
    reg_force    = reg_term*curvature;
    
    dphi_dt = evolve_force./(max(abs(evolve_force(:)))+eps) + reg_force;
    delta_t = 1/(max(abs(dphi_dt(:)))+eps);          % Step size using CFL
%     delta_t = 2;
    
    prev_mask = phi >=0;
    
    phi = phi + delta_t*dphi_dt;
    phi = SussmanReinitLS(phi,0.5);
    phi = NeumannBoundCond(phi);
    
    if display > 0
        if mod(its,display) == 0
              displayContour(phi,u,mag,color); drawnow; title(num2str(its));     
        
        if opt.evolution 
            frm = getframe();
            fname = strcat(opt.fnameevolve,num2str(its),'.png');
            imwrite(frm.cdata,fname);
        end
        end
    end
    
    curr_mask = phi >=0 ;
    
    count = convergence(prev_mask,curr_mask,convg_err,count);
    % count how many succesive times we have attained convergence, reduce local minima
    if count <= count_lim
        its = its + 1;
        disp(its);
    else
        stop = 1;
        disp('converged');
    end
    
end
end




% Convergence Test
function c = convergence(p_mask,n_mask,thresh,c)
    diff = p_mask - n_mask;
    n_diff = sum(abs(diff(:)));
    if n_diff < thresh
        c = c + 1;
    else
        c = 0;
    end
end


% Compute curvature    
function k = curvature_central(u)                       

    [ux,uy] = gradient(u);                                  
    normDu = sqrt(ux.^2+uy.^2+1e-10);	% the norm of the gradient plus a small possitive number 
                                        % to avoid division by zero in the following computation.
    Nx = ux./normDu;                                       
    Ny = uy./normDu;
    nxx = gradient(Nx);                              
    [~,nyy] = gradient(Ny);                              
    k = nxx+nyy;                        % compute divergence
end


% Check boundary condition
function g = NeumannBoundCond(f)
    
    [nrow,ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
end