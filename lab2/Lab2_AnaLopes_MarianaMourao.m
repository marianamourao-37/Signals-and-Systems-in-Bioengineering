%Group 22: 
%-Ana Rita Lopes nº98587
%-Mariana Mourão nº98473

%LAB#2 - Interpolation

%% I - Introduction

%% II - 1D INTERPOLATION

    %% 1)  
% Derivation of the expression in order to obtain the vector c that minimizes the equation 2 for 
% the Energy function E(c) (cost function), given 
% a set of random samples T and respective observations F.   
% It's used the matrix in expression 4 to re-write the equation 2 in matrice notation, being the 
% minimization process performed by applying the gradient operator, in which grad(E(c))=0 gives the stationary solutions.
% The final expression is given by: c*=((phi'*phi)^-1*phi')*F
% The complete derivation can be accessed on the paper sent in pdf format 


    %% 2) 

M1 = 10; % number of random points for sampling f (t) to F1

%Calls an auxiliary function which abstracts the method for preventing the repetition of time samples
T1 = generate_T(M1); % M1x1 column vector with M1 random samples within the interval I = [0,1]

F1 = sin(2*pi*T1).*exp(-20*(T1-0.5).^2); % M1x1 column vector with the set of observations F1 taken in T1 


%Simulates a continuous representation of f(t) 
L = 500; % total number of densed samples
t_dense = (0:1/(L-1):1)'; % Lx1 column vector with L time samples equally spaced in I = [0,1]
f_dense = sin(2*pi*t_dense).*exp(-20*(t_dense-0.5).^2); %Lx1 column vector with the values of f(t) 
%taken in the time points represented in t_dense 

figure(1);
plot(T1,F1,'rx',t_dense,f_dense,'b');
title('Sampling of f(t) for M = 10 random samples');
legend('Continuous function', 'Discrete Samples'); 
xlabel('t');
ylabel('f(t)');

    %% 3) 
N1 = 10; % number of interpolating functions 

%Calls an auxiliary function which abstracts the method for generating f^(t), given a set of random samples T1 and 
%respective observations F1 previously computed for M1 = 10, the number of interpolating functions N1 and the vector 
%t_dense with high density time points for simulating a continuous representation
[~,f_hat] = generate_fhat(N1,F1,T1,t_dense);

% Display graphically the function f(t) and the interpolated one f^(t), for the same density time samples
figure(2);
hold on;
plot(t_dense,f_dense,'LineWidth',2);
plot(t_dense, f_hat);
plot(T1,F1,'rx');
legend('Continuous function','Interpolated function','Discrete samples');
xlabel('t');
ylabel('f(t)');
title('Interpolation for N=10 interpolation functions and M=10 random samples');

%The signal to noise ratio (dB) is given by the equation (8), being simplefied by taking into account the euclidian norm and the cancelation of the (1/L)^2 factor.

SNR1 = 10*log10((f_dense'*f_dense)/((f_dense - f_hat)'*(f_dense - f_hat)));

%Prints the value of the SNR into the command window
fprintf('The SNR between the original function and the estimated is : %4.2f  dB\n', SNR1);

%Due to the implicite randomness involved in the selection of samples from which the f_hat is derived, this value will vary depending on how the samples are distributed. 
%However, for the consider N and M, the respective SNR are positive and generaly higher than 30 dB, demonstrating that the deployed reconstruction of f(t) was successful

    %% 4)

N2 = 50; % number of interpolating functions 

%Calls an auxiliary function which abstracts the method for generating f^(t), given a set of random samples T1 and 
%respective observations F1 previously computed for M1 = 10, the number of interpolating functions N2 and the vector 
%t_dense with high density time points for simulating a continuous representation
[~,f_hat_1] = generate_fhat(N2,F1,T1,t_dense);

% Display graphically the function f(t) and the interpolated one f^(t), for the same density time samples
figure(3);
hold on;
plot(t_dense,f_dense,'LineWidth',2);
plot(t_dense, f_hat_1);
plot(T1,F1,'rx');
legend('Continuous function','Interpolated function','Discrete samples');
xlabel('t');
ylabel('f(t)');
title('Interpolation for N=50 interpolating functions and M=10 random samples');

SNR2 = 10*log10((f_dense'*f_dense)/((f_dense - f_hat_1)'*(f_dense - f_hat_1)));
%Prints the value of the SNR into the command window
fprintf('The SNR between the original function and the estimated is : %4.2f  dB\n', SNR2);

M2 = 100; % number of random points for sampling f (t) to F

%Calls an auxiliary function which abstracts the method for preventing the repetition of time samples
T2 = generate_T(M2); % M2x1 column vector with M random samples within the interval I = [0,1]

F2 = sin(2*pi*T2).*exp(-20*(T2-0.5).^2); % M2x1 column vector with the set of observations F2 taken in T2 

%Calls an auxiliary function which abstracts the method for generating f^(t), given a set of random samples T2 and 
%respective observations F2 previously computed for M2 = 100, the number of interpolating functions N2 and the vector 
%t_dense with high density time points for simulating a continuous representation
[~,f_hat_2] = generate_fhat(N2,F2,T2,t_dense);

% Display graphically the function f(t) and the interpolated one f^(t), for the same density time samples
figure(4);
hold on;
plot(t_dense,f_dense,'LineWidth',2);
plot(t_dense, f_hat_2);
plot(T2,F2,'rx');
legend('Continuous function','Interpolated function','Discrete samples');
xlabel('t');
ylabel('f(t)');
title('Interpolation for N=50 interpolating functions and M=100 random samples');

SNR3 = 10*log10((f_dense'*f_dense)/((f_dense - f_hat_2)'*(f_dense - f_hat_2)));
%Prints the value of the SNR into the command window
fprintf('The SNR between the original function and the estimated is : %4.2f  dB\n', SNR3);

%For N>M, the linear system of equations for calculating the coefficients c is underdetermined, 
%since it has more parameters to estimate than equations, existing an infinite number of 
%solutions (althougth the method applied retrieves the minimum-norm least-squares solution). 
%In order to be determined, the number of samples would have to be M>=50, case where the 
%considered model would try to infer patterns just based on existing samples, having a higher 
%generalization capability (doesn’t overfitts the data). As a consequence of overfitting the data, 
%the noise generated in the interpolation process dominates the percentage of the 
%reconstructed function, thus obtaining SNR < 0, being the error significantly large for this 
%reconstruction.
%For M>N, the linear system of equations for computing c is determined, resulting as well in a better 
%sampling of the signal (average distancing between samples is lower), making the model generalize the 
%function with a higher precision. This translates into a higher SNR, reveling that the reconstruction of the function
%was successful

    %% 5)  

% Derivation of the expression in order to obtain the vector c that minimizes the equation (9) for 
% the Energy function E(c) (cost function), given a set of random samples T and respective observations F, as well as alpha. 
% The final expression is given by: c* =((?'*?+a*A'A)^-1*?')*F, with A being a (N-1)xN matrix:
% where A = [-1 1 0 0 ... 0
    %        0 -1 1 0 ... 0
	%        ..............
	%		 ..............
	%		 .............. 
    %        0 ... 0 -1 1 0
    %        0 0 ... 0 -1 1]
% The complete derivation can be accessed on the paper sent in pdf format 


    %% 6) 
% log scale in order to search for alphas closer to 0
list_alpha = [0, logspace(-10,log10(10),15)]; % logspace generates 15 points between decades 10^-10 and 10^(log10(10)) = 10. 
%In total, it's a row with 16 elements (accounting for alpha = 0, giving equation (2), without averaging the error (1/M factor)) 

num_alpha = length(list_alpha); % number of elements = 16 
A = [zeros(N2-1,1), eye(N2-1)] + [-eye(N2-1), zeros(N2-1,1)]; % (N2-1)xN2 matrix 

fig5 = figure(5);
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text( 0.5, 0, 'Sinc Interpolation for N = 50 basis functions and M = 10 samples',...
     'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'Bottom' ) ;

%test the model for different alphas (hyperparameter):
SNR_v = zeros(1,num_alpha); %allocating memory for storing the corresponding SNR for a given alpha 
for i = 1:num_alpha
    alpha = list_alpha(i);

	%Calls an auxiliary function which abstracts the method for generating f^(t), taking into account the minimization of the Energy Function given by 
	%equation (9)
    [~,f_hat] = generate_fhat(N2,F1,T1,t_dense,A,alpha);

    SNR_v(i) = 10*log10((f_dense'*f_dense)/((f_dense - f_hat)'*(f_dense - f_hat)));
    
	% Display graphically the function f(t) and the interpolated one f^(t) (for the ith alpha), for the same density time samples
    figure(5);
    subplot(4,4,i);
    hold on;
    plot(t_dense,f_dense,'LineWidth',2);
    plot(t_dense, f_hat);
    plot(T1,F1,'rx');
    hold off;
    xlabel('t');
    ylabel('f(t)');
    title(['alpha = ', num2str(alpha)]);
    
end
legend('Continuous function','Interpolated Function','Discrete samples');

% Retrieve the alpha from list_alpha that led to the maximum SNR computed 
[~,i_max] = max(SNR_v);
alpha_opt = list_alpha(i_max);

%Calls an auxiliary function which abstracts the method for generating f^(t), taking into account the minimization of the Energy Function given by equation (9), 
%considering the optimal alpha that led to better results in the previous case 
[~,f_hat] = generate_fhat(N2,F2,T2,t_dense,A,alpha_opt);

% Display graphically the function f(t) and the interpolated one f^(t), for the same density time samples
figure(6)
hold on;
plot(t_dense,f_dense,'LineWidth',2);
plot(t_dense, f_hat);
plot(T2,F2,'rx');
legend('Continuous function','Ideal Interpolation','Discrete samples');
xlabel('t');
ylabel('f(t)');
title(['Interpolation for N = 50, M= 100 and alpha = ', num2str(alpha_opt)]);

% Display graphically the computed values of SNR for the different alphas considered 
figure(7)
plot(list_alpha, SNR_v, '*-');
xlabel('alpha');
ylabel('SNR (dB)');
xlim([-0.5 10]) % sets the x-axis limits for the current axes or chart.
title('SNR after Interpolation for N = 50 and M= 10');

%INCLUDES THE ANSWER FOR II.7)
%As previously explained in exercise II.4), when there are more coefficients c to compute (N) then time samples (M), i.e, N>M, the considered model overfits the data, 
%since it has more parameters to estimate than equations, making the linear system underdetermined. In order to prevent this, an additional dependency constraint must be 
%introduced in order to effectively diminish the pseudo-order of the model, i.e, the number of free parameters. This can be achieved by taking into account the distance 
%between contiguous interpolators, being its minimization regulated by the hyperparameter alpha. If alpha =  0, the expression translates into equation (2), without 
%averaging the error (1/M factor), although this hasn’t any effect on deriving c (the result was the same as for exercise II.4). When alpha is considered, a dependency 
%between interpolators is introduced, which forces the interpolators to have more similar coefficients for deriving f_hat, making the system converge into a solution. For the 
%alphas being tested, the optimal value was 10^-10 (leads to an higher SNR), suggesting that the model didn't need much guidance to converge. As alpha increased, it started 
%to be observed a deterioration of the reconstructed function, since it is converging for the limit situation where the coefficients of the interpolators are equal, leading 
%to f_hat converging to a constant function. 


    %% 8) 
% Considering the Energy function E(c) (cost function) given by equation (10), the coefficients c that minimizes E(c) was derived by reformulating 
% |c_k - c_(k-1)| = (c_k - c_(k-1))^2/|c_k - c_(k-1)|, establishing beta_k = |c_k - c_(k-1)|^(-1)
% The final expression is given by: c*=((phi'*phi + alpha*(A'*beta*A))\phi')*F
% The complete derivation can be accessed on the paper sent in pdf format 


%It is of interest to test the designed algorithm for the critical situation where N>M

%initial alpha being the optimal one obtained in the previous exercise 
alpha8 = alpha_opt; 

%initial estimates for the coefficients, taking into consideration the minimization of the Energy Function given by equation (9), for the optimal alpha
[c_i,f_hat_init] = generate_fhat(N2,F1,T1,t_dense,A,alpha8);
SNR_prev =10*log10((f_dense'*f_dense)/((f_dense - f_hat_init)'*(f_dense - f_hat_init))); % initial estimate of SNR

%Computation of beta given c 
c_diff = abs(c_i(2:N2) - c_i(1:N2-1)); % (N2-1)x1 vector containing the absolute difference between contiguous coefficients
beta = c_diff'.^(-1).*eye(N2-1); % (N2-1)x(N2-1) matrix, having in its diagonal the inverse of the absolute difference between contiguous coefficients

max_iter = 500; %maximum number of iterations 
i = 1; % initialize counter of iterations 

SNR_best = SNR3; % higher SNR obtained previously 
incr_alpha = 10^(-11); % increment of alpha to be performed in the iterative process 

figure(8)
while SNR_prev < 0.8*SNR_best && i <= max_iter
    %stop when there's convergence
	
	%Calls an auxiliary function which abstracts the method for generating f^(t), taking into account the minimization of the Energy Function given by 
	%equation (10)
    [c,f_hat8] = generate_fhat(N2,F1,T1,t_dense,A,alpha8,beta);
        
    SNR_after = 10*log10((f_dense'*f_dense)/((f_dense - f_hat8)'*(f_dense - f_hat8)));

    dif_SNR = SNR_after - SNR_prev; % determine the difference between consecutive iterations 
    
    if dif_SNR > 0 %if the SNR is increasing 
        alpha8 = alpha8 - incr_alpha; % increases the alpha to accelerate convergency 
    else %if the SNR is decreasing 
        alpha8 = alpha8 + incr_alpha; % diminishes the alpha to approach the direction of convergence
    end 
    
     hold on;
     plot(i,SNR_prev,'*-k');
     ylabel('SNR(dB)');xlabel('No. of iterations')
    
	%These lines can be uncommented if desired to re-estimate the beta parameter given the new estimates of c. Being commented, implies that the 
	%iterative process fixed the initial estimate of beta
%     c_diff = abs(c(2:N2) - c(1:N2-1));
%     c_diff(c_diff == 0) = 0.00001;
%     beta = c_diff'.^(-1).*eye(N2-1);
    
    SNR_prev = SNR_after; 
    i = i + 1; %new iteration 
end

hold off;
  
% Display graphically the function f(t) and the interpolated one f^(t) obtained in the final iteration, for the same density time samples
figure(9)
hold on;
plot(t_dense,f_dense,'g','LineWidth',2);
plot(t_dense,f_hat8,'r');
plot(T1,F1,'*k');
legend('Continuous function','Ideal Interpolation','Discrete samples');
title(['M=10, N=50, alpha = ', num2str(alpha8)]);
hold off
 
 
 %% III - 2D INTERPOLATION   
 
    %% A)
    
    %% 1)

%Load the image to be analysed
x=imread('len_gray.jpg');
figure(10)
subplot(1,2,1)
colormap gray;  %Gray colormap array
imagesc(x)
title('Original Lena Image')

    %% 2)

y=x(1:10:end, 1:10:end);
subplot (1,2,2)
colormap gray
imagesc(y)
title('Image of Lena reduced 10 times - Downsampling')

%The performed operation of removing 9 pixels every 10 results in the reduction by a factor of 10 of the total number of pixels representing the image, 
%corresponding to resampling the image taking into consideration fs_y=1/10 fs_x, being this process called downsampling. 
%Instead of resizing the dimensions of the image, this process results in the increasing of the individual size of each 
%pixel (thus the decrease of fs_y), leading to degradation of the image resolution due to pixelization (spatial aliasing)

    %% 3)
B=fir1(100,0.1); %computes the N+1 (100 + 1 = 101) coefficients of a 100-order lowpass FIR filter with 0.1 cut-off frequency 
%(0.1 = normalized angular frequency/pi, assuming 1 if being the Nyquist frequency), considering a Hamming window

%convolving the filter with the image in each dimension 
xLP=conv2(B,B',x,'same');
%By specifying 'same', The conv2 function returnes the central part of the convolution

%reduction by a factor of 10 of the number of pixels of the filtered image 
w=xLP(1:10:end, 1:10:end);

figure(11)
subplot(1,2,1)
colormap gray;
imagesc(xLP)
title('Low pass Filtered image of Lena')
subplot(1,2,2)
colormap gray
imagesc(w)
title('Image filtered of Lena reduced 10 times - Decimation')

%In this exercise, the original image is low pass filtered before downsampling, in order to mitigate the distortion caused by aliasing, corresponding to a  
%process called decimation. By firstly filtering the image and consequently attenuating frequencies higher than pi/10, the high frequencies 
%associated to the intensity transitions between pixels are eliminated, being associated with lower frequencies (more smooth transitions).
%Given this, by performing downsampling on the filtered image the resulting image will have a higher resolution due to less pixelization. 
%However, aliasing isn't completely avoided, since the applied filter isn't ideal.

    %% B)

    %% 1)

%
xDirac = zeros(size(x));
xDirac(1:10:end,1:10:end)= y;

%Removes 9 pixels every 10 by setting them to 0, with the kept pixels remaining in their original spatial locations. This results 
%in the segmentation of the image into blocks, each one containing in the upper left corner the pixel value of the 
%original image and the remaining pixels of the block equal to zero. 

figure(12)
subplot(1,2,1)
colormap gray;
imagesc(x)
title('Original Lena Image')
subplot(1,2,2)
colormap gray;
imagesc(xDirac)
title('Image of Lena split in 10*10 pixel block')
 
%Once obtained the reconstructed image, we have high difficulty in
%recognizing the original image an its details, since we only obtain
%one color pixel in every 10*10pixel, so the distance between them
%makes only visible a set of points (the only values different from zero).

    %% 2)
pixelling_effect = conv2(xDirac,ones(10),'same');
%zero order interpolation by putting all pixels equal to the upper left corner one in each block.
%By specifying 'same', the conv2 function returns the central part of the convolution, which is the same size as xDirac

figure(13)
subplot(1,3,1)
colormap gray;
imagesc(x)
title('Lena')
subplot(1,3,2)
imagesc(y)
title('Image of Lena reduced 10 times - Downsampling')
subplot(1,3,3)
colormap gray;
imagesc(pixelling_effect)
title('Reconstruction of xDirac')

%In this exercise, we obtained an image in which the detail is
%distinguished, as the blocks of 10 by 10 pixels present the same intensity value, obtained by applyging piecewise constant interpolation
%We observed a similar result to the one visualize in group A in question A2), since the only diference is that 
%the total number of pixels are equal to the orignal image, but having blocks of 10x10 pixels with the same intensity value 
%as the one perserved in the downsampling image of exercise A2). By comparing the images, it's concluded that 
%each pixel in the downsampling image of exercise A2) is equivalent to the 10x10 pixel blocks in the pixelling_effect image. 

    %% 3) and 5) 

P=10;
Q=5*P;   

%computes the Q+1 coefficients of the Q-order lowpass FIR filter with 0.1 cut-off frequency 
%(0.1 = normalized angular frequency/pi, assuming 1 if being the Nyquist frequency), considering a Hamming window
h=fir1(Q,1/P,hamming(Q+1));  
%fir1 uses a least-squares approximation to compute the filter coefficients and then smooths the impulse response 
%by applying the hamming window

%convolving the filter with the image in each dimension 
xSinc=conv2(h,h',xDirac,'same');
%By specifying 'same', the conv2 function returnes the central part of the convolution

figure(14)
subplot(1,2,1)
colormap gray;
imagesc(pixelling_effect)
title('Reconstruction of xDirac')
colormap gray
subplot(1,2,2)
imagesc(xSinc)
title('Image xDirac with low-pass filter')

%The xDirac image can be interpreted as upsampling the image obtained in A.2) by a factor of 10, consisting of inserting 9 zeros between every pixel of the 
%y image. In order to infer the additional pixel intensities, it's necessary to interpolate those values given the known pixel intensities. In the previous 
%exercise, piecewise constant interpolation was performed, giving rise to pixelization due to abrupt intensity transitions between pixels. In this exercise, 
%it's applied a low pass filter in which its coefficients varies according to the defined hamming window, reducing the contrast and abrupt changes 
%of the intensity transitions between pixels, like it has been seen in exercise A.3), giving an apperance of cohesion and softening to the reconstructed image, 
%diminishing the pixeliing effect that was verified in the question before.

    
    %% 4)

grid_scale = -1:0.1:1; 
[X,Y] = meshgrid(grid_scale); 
vectors_2D = sqrt(X.^2+Y.^2); % computes the corresponding vectors in 2D
sinc_reconst = conv2(xDirac,sinc(vectors_2D)); %Ideal interpolation
    
figure(15)
subplot(1,3,1); colormap gray; imagesc(x);
title('Original Lena Image')
subplot(1,3,2); colormap gray; imagesc(xSinc);
title('Reconstructed with 50-order filter')
subplot(1,3,3); colormap gray; imagesc(sinc_reconst);
title('Ideally Reconstructed Image - Sinc Interpolation')

figure(16)
mesh(X,Y,sinc(vectors_2D));
title('2D Sinc Interpolator')
%projection in 2D plane gives rise to circular patterns

%As confirmed in group II of this Lab, the sinc function has great capability for interpolating functions. Generalizing to 2D, 
%an ideal image reconstruction can be obtained by considering a sinc interpolator 2D kernel, by performing the convolution operation between this kernel 
%and the xDirac image. For each pixel, the kernel is centered, making the weighted (given by the coefficients) average of the intensities of the neighboring 
%pixels within an area of 21x21 pixels (dimension of the kernel), with the weights of the neighboring pixels being lower as more distant they are from the
%central pixel. For the pixels that contain information (different from 0), this weighted average gives its value, since the central pixel of the kernel equals 1.
%The image interpolated with the filter of order 50 seems better, but an analysis of the corresponding SNR would discern better. 


 