function y = chirpTone(T,f1,f2,fs)
    %returns a sinusoidal function modulated by f(t), being digitalized 
    %by the sampling frequency fs. The time duration of the signal is T, 
    %with the function f(t) discritization being normalized to the f2-f1 scale, having
    %N samples that scale. 

    %number of samples: optimization of the computation to be performed for the 
    %frequency discritization, since it will not be needed to compute the 
    %length of a vector (the time complexety would be linear and proportional 
    %to the N length. Instead, by making this computation, the time complexety is constant) 
    N = (T*fs) + 1;
    
    % Ts = 1/fs, corresponding to the time between adjacent samples. Given the
    % desired time duration T of the considered signal, it can be calculated
    % how many samples are used to perform the signal digitalization/sampling. 
    % The added sample is due to the fact that T*fs gives the number of 
    % pairs of adjacent samples, which corresponds to N - 1 pairs.
    
    %time discretization scale vector 
    t = ((0:N-1)/fs)'; % tk = nk/fs = nk*Ts, with Ts = 1/fs. tk gives the time point 
    %corresponding to the sumatory of time intervals between the previous
    %adjacent samples, until reaching k. This time discritization ends at
    %N-1/fs, because it must contain only N samples in total,  which account 
    % for the first one being automaticaly counted (t = 0 seconds)
    
    % higher time complexity computation, O(N):
    %f = (f1:(f2-f1)/(length(t)-1):f2)';
    
    %lower time complexity computation, O(c), with c being constant. 
    %linear discritized function, for increasing frequency, being stored
    %into a column vector:
    f = (f1:(f2-f1)/(N-1):f2)'; 
    % function f(t) contains values from f1 to f2, being its
    % discritization normalized to the f2-f1 scale, having N discretized 
    % samples within that scale, accounting for the first one being 
    % automaticaly counted - reason for dividing by N-1 
    
    %output vector
    y = sin(2*pi.*f.*t); 
end