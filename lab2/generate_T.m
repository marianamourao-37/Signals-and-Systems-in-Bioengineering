function T = generate_T(M)
%Given a integer M higher than 0, it generates a Mx1 column vector T that contains M random time samples uniformly distributed in the interval I = [0,1]. 
%Accounts for the repetition of time samples, in which is repeated the random process for generating T, given a new seed 
%(avoids repeatability of the obtained results)

    i = 1; %initializing the seed, in case repetition of time samples has occured 
    T = rand(M,1); %generates a Mx1 column vector that contains M random time samples uniformly distributed in the interval I = [0,1]

	% check if there are any repetitions of time samples (not intended) 
    statement = true;
    while statement 
        T_unique = unique(T); % removes duplicates 
        if length(T) ~= length(T_unique) %if there was any duplicates 
            i = i + 1; % creating new seed for avoiding repeatability of the obtained results
            rng(i); % new seed 
            T = rand(M,1); 
        else
            statement = false; % there aren't any duplicates 
        end
    end
end