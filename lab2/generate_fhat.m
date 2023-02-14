function [c,f_hat] = generate_fhat(N,F,T,t_dense,A,alpha,beta)
% Generates a Nx1 column vector c with the coefficients for the N interpolating sinc functions parameterized by N 
% (determines the scaling factor delta and shifting factor k, in order for scanning the entire interval I = [0,1]), returning as well a  
% Lx1 column vector simulating the continuous interpolated function taken in the time samples t_dense. 
% This function is extended to the other optimization problems that take into account the Energy Function defined in the equation (9) - A and alfa are additional inputs - 
% and (10) - A, alpha and beta are aditional inputs. 

%Inputs:
%N - number of interpolating functions 
%F - column vector with the set of observations 
%T - column vector with random time samples 
%t_dense - column vector with the time samples for simulating a continuous representation
%A - (N-1)xN matrix for performing the sum(c_k - c_(k-1)), for k = 1 to k = N - 1
%alpha - hyperparameter that regulates the convergency of the interpolation 
%beta - (N-1)x(N-1) matrix, having in its diagonal the inverse of the absolute difference between contiguous coefficients

    delta = 1/(N - 1); %scaling factor taking into account the interval I = [0,1]
    k = 0:N-1; % shifting factor for covering the interval
    phi_discrite = sinc(T/delta - k); %MxN matrix, being in each column the estimate of the Nth interpolating sinc function for the T random time samples. 
		
	%Taking into consideration that when N>M the obtained linear system of equations is underdetermined, i.e, involves more unknowns than equations, 
	%infinitely many solutions are available. The matrix left division operation in MATLAB finds a basic least-squares solution, being of particular interest finding 
	%the solution with minimum norm. Thus, the lsqminnorm built-in function was used, computing the minimum-norm least-squares solution. 
	%for AX = B, X = lsqminnorm(A,B)
	
	%estimated coefficients of the N interpolating sinc functions, derived by minimizing the respective Energy Function
    if nargin == 4 % optimization problem regarding the Energy function for the equation (2) 
        c = lsqminnorm(phi_discrite'*phi_discrite,phi_discrite'*F); 
        %c = ((phi_discrite'*phi_discrite)\phi_discrite')*F;
		%Expression derived in annex II.1) 
		
    elseif nargin == 6 % optimization problem regarding the Energy function for the equation (9) 
        c = lsqminnorm(phi_discrite'*phi_discrite + alpha*(A'*A),phi_discrite'*F);
        %c = ((phi_discrite'*phi_discrite + alpha*(A'*A))\phi_discrite')*F;
		%Expression derived in annex II.5)
		
    else % optimization problem regarding the Energy function for the equation (10) 
        c = lsqminnorm(phi_discrite'*phi_discrite + alpha*(A'*beta*A),phi_discrite'*F);
        %c = ((phi_discrite'*phi_discrite + alpha*(A'*beta*A))\phi_discrite')*F;
		%Expression derived in annex II.8) 
		
    end
	
    f_hat = sinc(t_dense/delta - k)*c; % Lx1 column vector containing a simulation of a continuous representation of the interpolated function f_hat for the time samples 
	%t_dense, obtained by linear combination (weighted by the estimated coefficients) of the interpolating sinc functions in those same time samples 
	
end