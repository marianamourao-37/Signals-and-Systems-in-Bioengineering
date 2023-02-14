function [error_matrix, error_min, x_min, y_min] = Error(x,y,g,p1,p2,N)

%receives as inputs the x-y plane axes vectors of the error surface (x and y), 
%a given vector g contained in an unknown space, the basis vectors 
%(p1 and p2) and the number of discretized points (N).
%As outputs, it retrieves the error matrix (error_matrix), the minimum error (error_min) 
%and the value of both x and y where the error is minimum (x_min and y_min)

    error_matrix = zeros(N); %initializes the error matrix, being square (NxN)
    
    for i=1:N  %iterates over x
        for j=1:N  %iterates over y
            v = (g-x(i)*p1-y(j)*p2); % for the coefficients x(i) and y(j), computes the 
            % difference between g and g_approx = x(i)*p1-y(j)*p2
            
            error_matrix(i,j)=v*v'; % error corresponding to the coefficients 
            %x(i) and y(j) 
        end
    end
    
    %indexes corresponding to the minimum value in the error matrix 
    [i_min,j_min] = find(error_matrix == min(error_matrix(:)));
    error_min = error_matrix(i_min, j_min); % minimum error 
    x_min=x(1,i_min); % x corresponding to the minimum error, i.e, coefficient 
    %of p1 which minimizes the norm of the error 
    y_min=y(1,j_min); % y corresponding to the minimum error, i.e, coefficient 
    %of p2 which minimizes the norm of the error  
end