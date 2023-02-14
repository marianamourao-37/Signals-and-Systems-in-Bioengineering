%Group 22: 
%-Ana Rita Lopes nº98587
%-Mariana Mourão nº98473

%LAB#1 - Space of Signals
%% 1)
    %% a) function define in the file named 'ChirpTone'

    %% b) 
T = input("Enter the time duration (in seconds) of the signal, for the 1. exercise: ");
f1 = 100;
f2 = 2000;
fs = 4000;
        
x = chirpTone(T,f1,f2,fs);
soundsc(x,fs);

%The considered signal has its frequencies modulated by the time-varying phase, tetha(t) = 2*pi*f(t)*t. For each time t, the frequency 
%of the sinusoidal is called instantaneous frequency, being the time derivative d/dt of the signal's phase, theta(t), giving rise to 
%the following expression for the angular instantaneous frequency w_inst:
% w_inst = (d/dt)*theta(t) = 2*pi*[f'(t)*t + f(t)]
% f_inst = w_inst/(2*pi) = f(t) + f'(t)*t

% f(t) is the growing frequency defined as: f(t) = f1 + (f2 - f1)*t/T
% Thus, the time derivative of f(t) is: f'(t) = (f2 - f1)/T_bs

%Substituing f(t) and f'(t) into the expression of f_inst, arises the following:
%f_inst = f1 + 2*(f2-f1)*t/T

%Thus, for t = T, the heared frequency is f_inst = 2*f2 - f1

%After hearing the digitalized signal with the mentioned specifications, it would be expected that its frequency would increase with time. 
%Although, considering the above derivation of f_inst, it's concluded that the created signal contains frequencies from 
%f1 (=100 Hz) up to 2*f2 - f1 (=3900 Hz). Due to an inadequate digital sampling rate, the contained frequencies that aren't lower than fs/2 
%(Nyquist frequency) will be underestimated, resulting in the observed distortion of the digitalized signal. Modifying fs = 3900*2 = 7800 Hz would
%produce a signal with growing frequencies over time.


    %% c) 
%function to write an audio file, given a specific sampling frequency (fs).
%it is important to preserve the sampling frequency between processing and audition, 
%in order to ensure that we're audio-accessing the same digitalized signal, i.e, the same 
%dicretized frequencies and time samples. If desired re-sampling a digitalized signal, 
%appropriate interpolation is needed. 
audiowrite('chirp.wav',x,fs);

    %% d)
	
%Problem formulation:
%In order to simplify the analysis, it will be considered an ambulance
%moving restritively in the horizontal direction, from left to right.
%Additionaly, the ambulance movement in relation to the stationary person
%will be symmetricaly formulated, being restrited to an equal initial and final 
%relative positions:
d_horizontal = input("Enter the relative horizontal distance between stationary person and initial and final positions of the ambulance: ");

%Ambulance's horizontal velocity
amb_v = input("Enter a high velocity (in m/s) experimented by an ambulance: ");

t_half = d_horizontal/amb_v;  % time (in seconds) that the ambulance takes to travel to 
%the stationary person

%The signal of a siren can be decomposed into a basis signal, consisting in
%2 cycles (one in which the frequencies increase, and other where decrease)
Time_siren = input("Enter a time duration for the basis signal of a siren: ");

%frequencies to define the 2 distinguish cycles on the siren sound.
%Attention must be taken in order to not surpass fs/2, again taking into account that f(t = T) = 2*f2 - f1 
f1 = 960; 
f2 = 770;

%Basis signal, with each cycle having T_bs/2 duration.
cycle_1 = chirpTone(Time_siren/2, f1, f2, fs);
cycle_2 = chirpTone(Time_siren/2, f2, f1, fs);
sirene_bs = [cycle_1; cycle_2(2:end)]; 

%Looping of the basis signal, by repeating n times within the total time of 
%the recording (2*t_half)
sirene = repmat(sirene_bs,[floor(2*t_half*fs/length(sirene_bs)),1]);

%repmat considers complete repetitions of a basis signal. However, the
%relation between the duration of the basis signal (T_bs) and the total duration 
%of the looped signal (t_half*2) may require truncation of the basis signal
nppt = floor(2*t_half*fs)-length(sirene); %sample where the basis signal is truncated 
sirene = [sirene; sirene_bs(1:nppt)]; % siren sound

%The auditory system has underlying physiological mechanisms that allow stereo perception (i.e, if a sound is 
%located on the left or on the right side), including the detection of intensity-differences between both ears. 
%Given that it was considered that the ambulance is moving from left to right, the intensity captured for the left 
%channel will be higher than for the right channel until it crosses the person, position where intensities equalize 
%and from thereafter the intensity for the right channel becomes progressively higher than for the left channel.
%This intensity-differences can be modulated by the following attenuation functions:
right_ch = (0:1/(length(sirene) - 1):1); % attenuation for the right channel (attenuation decreases)
left_ch = 1-right_ch; % attenuation for the left channel (attenuation increases)
vector_attenuation = [left_ch' right_ch'];

%element-wise multiplication in order to create a stereo audio
sirene_stereo = [sirene sirene].*vector_attenuation;

soundsc(sirene_stereo, fs);

audiowrite('Ambulance_Sirene.wav',sirene_stereo,fs);

%The moving ambulance has a different pitch as it approaches, when it is closest to us 
%and as it passes us and drives away. This phenomena is called doppler
%effect. A better formulation of the problem would take this into account

%% 2. 
    %% a) 
[sd, fs]=audioread('Let It Be.mp3');
%reads data from the file named 'Let It Be', returning sampled data, sd, and a sample rate 
%(in hertz) for that data, fs.
    %% b) 
N = length(sd);
T = (N-1)/fs; 

% Ts = 1/fs, being the time interval between adjacent samples on
% the discritized signal. Thereby, the total time duration of the sampling
% signal can be computed by summing all time intervals (Ts) between
% adjacent samples. This can be obtained by multiplying N-1 samples by Ts, because
% for N samples there are N-1 pairs of adjacent samples. 

    %% c)
sd_backwards_l = sd(end:-1:1,1); %% left channel, being the 1st column on the matrix  
sd_backwards_r = sd(end:-1:1,2); %% right channel, being the 2nd column on the matrix  
sd_backwards = [sd_backwards_l sd_backwards_r]; 
% Concatenation of the two vectors, in order to form a matrix with two columns, 
%which corresponds to the left and right channels.

%Note that the audio inversion can be computed by using a built-in matlab function (flipud):
%sd_backwards=flipud(sd); 

%% Demonstration of audio inversion:
    %% 1st: Signal Plotting 
%we see in the plots the amplitudes of the signal in both channels. In figure 1 it´s the 
%amplitude for  the original song and in figure 2 the amplitude of the song when it is played backwards

figure(1)
subplot(1,1,1);
subplot(2,1,1), plot(1:N, sd(:,1));
title('Original Audio Signal - Left Channel');
xlabel('Samples');
ylabel('Signal Amplitude');

subplot(2,1,2), plot(1:N, sd(:,2));
title('Original Audio Signal - Right Channel');
xlabel('Samples');
ylabel('Signal Amplitude');

figure(2)
subplot(2,1,1), plot(1:N, sd_backwards(:,1));
title('Inverted Audio Signal - Left Channel');
xlabel('Samples');
ylabel('Signal Amplitude');

subplot(2,1,2), plot(1:N, sd_backwards(:,2));
title('Inverted Audio Signal - Right Channel');
xlabel('Samples');
ylabel('Signal Amplitude');

    %% 2nd Hearing the inverted audio vs original audio
soundsc(sd_backwards,fs);
soundsc(sd,fs);
    %% d) 
%Generate a stereo audio signal where the volume decreases to zero in the left channel and increases from zero 
%in the right one. The input values from audioread() are dimensionless, scaled to -1<=x<1, typically scaled 
%from a voltage reading that might be positive and negative or might be positive only.

mag_incr = 0:1/(N-1):1; %vector of increasing amplitudes                  
mag_descr = 1-mag_incr; %vetor of decrescent amplitudes 

%the initial sound magnification is totaly stored in the left channel, 
%being progressevely transfered to the right channel
total_mag = sd(:,1) + sd(:,2);
ch_l = total_mag.*mag_descr';
ch_r = total_mag.*mag_incr';

sd_movement= [ch_l ch_r]; % vetor concatenation, creating the resulting stereo audio

soundsc(sd_movement,fs);
 

%% 3.
    %% a)
T = 2;
fs = 4000; 
t = (0:1/fs:T); %time discretization scale vector 

    %% b) 
        %% i, ii and iii)
f = 440; % signal frequency

% sinusoidal signals with pure frequency f are created, being discretized 
% with a sample rate fs. The sinusoidal amplitudes, for each time sample,
% are obtain by applying element-wise multiplication 
p1 = sin(2*pi*f.*t); 
p2 = sin(2*pi*f.*t + pi/4);
p3 = sin(2*pi*f.*t)+ t;

    %% c) 
% The norm of a vector can be defined as the positive square root of the inner product of the vector with itself: 
%(Note: it could be used the norm built-in function)
p1_norm = sqrt(p1*p1');
p2_norm = sqrt(p2*p2');
p3_norm = sqrt(p3*p3');

fprintf('||p1|| = %f', p1_norm)
fprintf('\n||p2|| = %f', p2_norm)
fprintf('\n||p3|| = %f', p3_norm)

% Given the inner product mathematical formulation, cos(ang(a,b))=<a,b>/(||a||*||b||), 
% resulting in ang(a,b) = arccos(<a,b>/(||a||*||b||))
% Note: for computing the inner product it could have been used the dot built-in function 

ang_p1p2 = acosd(p1*p2'/(p1_norm*p2_norm));
ang_p1p3 = acosd(p1*p3'/(p1_norm*p3_norm));
ang_p2p3 = acosd(p2*p3'/(p2_norm*p3_norm));

fprintf('\np1^p2 = %f', ang_p1p2)
fprintf('\np1^p3 = %f', ang_p1p3)
fprintf('\np2^p3 = %f', ang_p2p3)

%In order to verify orthogonality, the angle between two vectors must be
%90º, being equivalent to the inner product of said vectors equal to
% zero, that is, the projection of each vector onto the other is zero. 

if ang_p1p2 == 90
    fprintf('\np1 and p2 are ortogonal');
else
    fprintf('\np1 e p2 are not ortogonal');
end

if ang_p1p3 == 90
    fprintf('\np1 and p3 are ortogonal');
else
    fprintf('\np1 and p3 are not ortogonal');
end

if ang_p2p3 == 90
    fprintf('\np2 and p3 are ortogonal');
else
    fprintf('\np2 and p3 are not ortogonal');
end

% The angle between p1 and p2 traduces the phase shift of 45º.
% The angle between p1 and p3 is 68º, and between p2 and p3 is 58º, being explained by the 
% linear growth (rate of 1/fs) of p3' amplitudes. 

    %% d)
% linear coefficients 
alpha = 0.5;
beta = 1;
gamma = 1.5;

% composed signal obtained by linear combination of base signals (p1, p2 and p3)
z=(alpha*p1 + beta*p2 + gamma*p3)';

%since every basis component has been discritized with a sampling frequency 
%fs = 4000 Hz, and there weren't added any more frequencies superior to the 
%Nyquist frequency (given by fs/2 = 2000 Hz), the recommended sampling frequency 
%for auditing is equal or higher than 4000 Hz. 

soundsc(z,fs);

    %% e) 
% The Graminian R is computed as R = A^H*A, with A being the matrix which 
% contains the basis vectors in the columns (A = [p1 p2 p3]), and A^H the hermitian, 
% i.e, transposed conjugated matrix. Since A is real, A is transposable, thus being 
% appliable R = A'*A
% At each position (i,j), the R matrix gives the inner product of the vectors 
%p(i) and p(j): R(i,j) = <p(i),p(j)>

A = [p1' p2' p3']; 
R = A'*A;

% If the vectors compose an orthonormated basis, then <pi,pj> = 1
% if i=j, being 0 otherwise, generating a graminian matrix R that corresponds to the
% identity matrix. 
% In the context of this problem, p1, p2 and p3 are not orthogonal vectors, 
% as seen in the previous exercise, so R does not correspond to the diagonal matrix.
[n, m] = size(R);

if R == eye(n)
    fprintf('\nThe Graminian is an identity matrix, implying that the basis vectors are ortogonal');
else
    fprintf('\nThe Graminian is not an identity matrix, implying that the basis vectors are not ortogonal');
end

%% 4.
    %% a) i) 

%According to the orthogonality principle, the norm of the error vector e is 
%minimum when the error vector is orthogonal to each basis vector (p1,p2 and p3), 
%i.e, <e,p_j> = 0, j = 1,2,3. Given that e(a) = z-z_aprox(a), where 
%z_aprox(n) = Sum{k = 1, k = L}(a_k(n)*p_k(n)), with L being the number of 
%basis vectors, the following can be derived:
% <e, p_j> , j =1,2,3
% <z - Sum{k = 1, k = L}(a_k*p_k), p_j> = 0
% <z, p_j> - <Sum{k = 1, k = L}(a_k*p_k), p_j> = 0
% Taking the matrice notation, A^H*z - R*a = 0, thus a = R^(-1)*A^H*z 
%(possible because the basis vectors are linearly independent, making R invertable)

a = R\(A'*z);

%Backslash or matrix left division. If A is a square matrix, A\B is roughly the 
%same as inv(A)*B, except it is computed in a different way. If A is an n-by-n 
%matrix and B is a column vector with n components, or a matrix with several 
%such columns, then X = A\B is the solution to the equation AX = B computed 
%by Gaussian elimination.


    %% a) ii) 

% Considering the gradient method, e(a) = z-z_aprox can be simplified into 
% e(a) = z - Phi*a, where z_aprox(n) = Sum{k = 1, k = L}(a_k(n)*p_k(n)), 
% with L being the number of basis vectors. Traducing z_aprox into matrice
% notation, we get z_aprox = p(NxL) * a(Lx1), with N being the number of
% samples. By analogy, it's concluded that phi = p(NxL) = [p1 p2 p3] 

phi = [p1' p2' p3'];

% Aiming the minimization of the cost function J(a) = ||e||^2, 
% it can be considered J(a)= ||(z-phi*a)||^2 = (z-phi*a)'*(z-phi*a). 
% Applying the gradient operator, Grad(J(a)), the stationary points 
% (maximum and minimum) are where it is equal to 0. Deriving the expression
% Grad(J(a)) = 0, it's obtained the solution a = (phi'*phi)^(-1)*phi'*z. 

a = inv((phi')*phi)*(phi')*z;

 
    %% a) iii)

% As expected, both methods applied in a)i) and a)ii) converge to the 
% coefficients (alfa,beta e gamma) used to define z, thus having minimize 
% the norm of the error vector to 0, i.e, z vector belongs to the space 
% defined by the basis vectors, where z_aprox belongs as well

    %% b)
    
%Mesh cells are used as discrete local approximations of the continuous space domain.
g = alpha*p1 + beta*p2;
xi = 0;
xf = 1;
yi = 0.5;
yf = 1.5;
N = 100;
x = xi:(xf - xi)/(N - 1):xf;
y = yi:(yf - yi)/(N - 1):yf;

[error_matrix, error_min, x_min, y_min] = Error(x,y,g,p1,p2,N);
%Function defined in the file named 'Error'

figure(1)
mesh(x,y,error_matrix)
hold on
xlabel('X');
ylabel('Y');
zlabel('Error');
plot3(x_min,y_min,error_min,'r.','MarkerSize',15)
hold off

%The obtained error surface is shaped as an elliptic paraboloid, converging unsymmetrically to a minimum , 
%depending on the direction of the descent, resulting in a differential speed of convergence to the
%minimum (if taken an iterative approach). The stationary point of the surface theoretically corresponds to 
%X(i) = 0.5 and Y(j) = 1.0, where the expected error = 0, since g = g_aprox.
%Although, the obtained numerical error on the minimum was 0.0598, due to the grid detail, i.e, the 
%discretization of the x-y plane defined by N. If N (grid detail) increases, we obtain a minimum error closer to 0.



    %% c) 
    
p4=sin(2*pi*f.*t);
p5=sin(2*pi*(2*f).*t);

g2=alpha*p4 + beta*p5;

[error_matrix, error_min, x_min, y_min] = Error(x,y,g2,p4,p5,N);

figure(1)
mesh(x,y,error_matrix)
hold on
xlabel('X');
ylabel('Y');
zlabel('Error');
plot3(x_min,y_min,error_min,'r.','MarkerSize',15)
hold off

%In contrast to b), the considered basis vectors are orthogonal, with the defined space having 
% a diagonal Gramian matrix, resulting in a surface error being shaped as a symmetrical elliptic paraboloid 
%(can be confirmed by a contour graph, displaying the equipotential lines, being circles in c) and ellipses in b)). 
%As a consequence of orthogonality, the concavity is more closed (larger curvature), and if adopted an iterative approach to 
%find the minimum, it would converge faster (less computational power required) to the stationary point, independently of the 
%direction of the descent path. Moreover, considering again an iterative approach, for a given defined epsilon, 
%deltaE < epsilon gives the condition to stop the iterative process, being deltaE the difference between the computed errors for  
%consecutive iterations. The flatter the curve, the worse the estimate of the minimum error, and consequently of the coefficients. 








    