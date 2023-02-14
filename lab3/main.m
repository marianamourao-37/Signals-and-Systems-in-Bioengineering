%Grupo 22: 
%-Ana Rita Lopes nº98587
%-Mariana Mourão nº98473
 
%LAB#3 - Space of Signals
 
clear all;
close all;
clc;
 
%% Parte 1)
    %% 1) Função definida o script designado 'DFT_mat'
    
%% 2)    
 %% i. 
 
f=440;   %frequência em Hz do sinal sinusóide puro
fs=4000; %frequência de amostragem em Hz
T=1;     %duração do sinal em segundos 

t = linspace(0,T,fs*T)'; %vetor de discretização do tempo (em segundos) 

%Implementacao da funcao de tom puro de frequência f
A = sin(2*pi*f*t); %vetor coluna  

%soundsc(A); % opcional de reproduzir o som

figure(1);
plot(A); 
xlabel('Samples N');
ylabel('Amplitude');
title(['Part1-2.i - Representation of a pure Tone with ', num2str(f), ' Hz']);
axis tight
%savefig('Part1-2.i.fig')

%% ii.
N=length(A);

% Admitiu-se uma DTF com comprimento igual a N
M = 2*N;
%M = input(['Enter the M-length of the DFT to be computed, in order to M >= ', num2str(N), ':']);

DFT_A = DFT_mat(A, M);  %computa a DFT de comprimento M usando a função definida em I.1)

k=(0:M-1)'; %fks = k*fs/M;
f = k*fs/M; %como fs = M, f = k 

%reorganiza a DFT movendo a componente DC de frequência zero para o centro do array, sendo 
%útil para termos de visualização da DFT
f_shift = (-length(DFT_A)/2:length(DFT_A)/2-1)*(fs/length(DFT_A)); % zero-centered frequency range

mag_dft = abs(DFT_A); %modulo da DFT 

figure(2);
subplot(121);plot(f,mag_dft);
xlabel('Frequencies (Hz)');
ylabel('Magnitude of DFT');
title('DFT of A Tone');
subplot(122);plot(f_shift,fftshift(mag_dft));
xlabel('Frequencies (Hz)');
ylabel('Magnitude of DFT');
title('zero-centered DFT of A Tone - fftshift command');

%máximo da DFT, devolvendo com 4 casas decimais 
XX=round(abs(DFT_A),4);
%indices dos coeficientes correspondentes ao máximo da DFT
Indxs = find(XX == max(XX));

ks = Indxs - 1; %devido à indexação do matlab, em que k=0 corresponde ao indice 1

%A correspondência entre indice k da DFT de compriemnto M e a frequência no sinal original, x, 
%é dada por f_k = f_s*k/M
fks = ks*fs/M;

fprintf('\n I.2.iii) O maior coeficiente ocorre para %d e %d\n',fks)

%% iii
%A DFT duma função seno corresponde a dois deltas
%dirac para a frequência fundamental (440) e a sua simétrica (-440). Como a dft assume uma 
%extensão periódica do sinal ao infinito em ambas as direções (isto é, o espectro é repetido a
%cada M amostras, ou seja, DFT(k)=DFT(k+M)), e devido ao facto de só se calcular a DFT para 
%frequências positivas, obtém-se para além do pico em 440 Hz um 2ºpico 
%em 4000-440 Hz = 3560 Hz que é a primeira harmônica do pico de -440 Hz no espectro original.
%A correspondência entre indice k da DFT de compriemnto M e a frequência no sinal original, x, 
%é dada por f_k = f_s*k/M

%Note-se que ao contrário do previsto através da transformada de fourier
%contínua, os picos obtidos não correspondem a deltas dirac. Tal deve-se ao
%facto de o suporte finito computacional apenas permitir uma representação
%finita do sinal x, correspondendo no dominio temporal a uma multiplicação do sinal x original
%por uma janela retangular. No dominio espectral tal é equivalente a convolver o espectro 
%do sinal original com o espectro da janela retangular. Uma vez que o espectro da janela retangular 
%é a função sinc(x), observam-se lóbulos sin(x)/x em torno de cada componente 
%espectral original conforme eles são convolvidos com o espectro da janela.
%Este fenomeno designa-se spectral leakage.

%% 3.
T=1; % duracao do sinal em segundos 
f = 440; %frequência em Hz do sinal sinusóide puro

%A implementação otimizada do algoritmo FFT implica que o sinal tenha dimensão 
%N = 2^u (com u inteiro)
N_pow2=2^(nextpow2(N));
M=N_pow2;
fs_1 = (N_pow2-1)/T; %frequencia de amostragem de forma a obter um sinal de tamanho N_pow2

t = (0:1/fs_1:T)'; %vetor de discretização do tempo (em segundos) 

A=sin(2*pi*f*t); % Implementacao da funcao de tom puro de frequência f

%usando a função tic e toc, é possivel obter uma medição precisa do tempo de 
%processamento de um determinado bloco.
%tempo de execução do algoritmo da DFT implementado pela função DFT_mat
tic;
[DFT_A]=DFT_mat(A,M);
tempo_dft=toc; 

fprintf('\n I.3. O tempo de computação da DFT usando a função matDFT é: %d\n',tempo_dft); 

k=(0:M-1)'; %= fk = k*fs/M, com fs = M
subplot(121)
plot(k,abs(DFT_A));
xlabel('Frequencies (Hz)');
ylabel('Magnitude of DFT');
title('I.3: DFT')

%tempo de computação da dft utilizando o algoritmo da FFT
tic;
FFT_A=fft(A); 
tempo_fft=toc; 

fprintf('\n I.3. O tempo de computação da DFT usando a FFT é: %d\n',tempo_fft); 

subplot(122)
plot(k,abs(FFT_A));
title('I.3: FFT')

    %% GENERAL DEMONSTRATION
%pow2_list = 2:12;
length_signal_list = 2:2:2048;
n_points = length(length_signal_list);
times = zeros(n_points,2);
T = 1;
f=440;

for i = 1:n_points
    length_signal = length_signal_list(i);
    %length_signal = 2^pow2_list(i);
    fs_1 = (length_signal-1)/T; %frequencia de amostragem de forma a obter um sinal de tamanho N_pow2
    t = (0:1/fs_1:T)';
    signal=sin(2*pi*f*t); 

    % Tempo de execução da FFT 
    tic;
    X_FFT = fft(signal);
    times(i,1) = toc;
    
    % Tempo de execução da DFT
    tic;
    X_DFT = DFT_mat(signal,length_signal);
    times(i,2) = toc;
end

% length_signal_list = 2.^pow2_list;

% Representação gráfica das complexidades temporais para sinais de diferentes complrimentos 
figure(4);
subplot(121);
hold on;
plot(length_signal_list,log(times(:,1)'));
plot(length_signal_list,log(times(:,2)'));
xlabel('Length of Signal');
ylabel('log(Computation time)');
hold off;
legend('FFT', 'DFT');
title('Time Complexity: DFT vs FFT');
axis tight;

subplot(122);
hold on;
plot(log2(length_signal_list),log(times(:,1)'));
plot(log2(length_signal_list),log(times(:,2)'));
xlabel('Power of 2');
ylabel('log(Computation time)');
hold off;
legend('FFT', 'DFT');
title('Time Complexity: DFT vs FFT');
axis tight;

% Comprova-se que a complexidade temporal da DFT é muito superior à da FFT,
% resultando do facto de para cada um dos N pontos multiplicam-se N
% exponenciais complexas por forma a obter o coeficiente da DFT correspondente, 
% envolvendo um custo total de N^2 (O(N^2)). 
% Já relativamente à complexidade temporal da FFT, o seu método de cálculo intrinseco 
% segue uma abordagem divide and conquer, eliminando as redudâncias associadas a
% cálculos repetidos na computação da DFT, reduzindo a complexidade quadrática 
% da DFT para uma complexidade O(NlogN). A otimização da DFT é então
% conseguida através do algoritmo FFT, sendo mais notório para N elevados. 

    %% Part II
    
    %% 1)   
f1=200; 
f2=1000;
fs=4000; %frequencia de amostragem 
T=1; %duração do sinal 
Chirp=chirpTone(T,f1,f2,fs); %vetor coluna contendo os valores da função chirp, 
%considerando a função chirpTone.m desenvolvida para o lab1.  

%soundsc(Chirp); %opcional, caso se pretenda ouvir o sinal 

t = (0:1/fs:T)'; %vetor de discretização do tempo (em segundos)

%Representação gráfica do sinal chirp 
figure(5)
plot(t,Chirp);
xlabel('Time(s)');
ylabel('Amplitude');
title('Chirp');
 
DFT_chirp = fft(Chirp); %Otimização da computação da DFT usando o algoritmo FFT
 
nfft=length(DFT_chirp);
%k = (0:nfft-1)';
f = (0:nfft-1)*fs/nfft;
f_shift = (-nfft/2:nfft/2-1)*(fs/nfft); % gama de frequências centrada em zero 

figure(6)
subplot(121);plot(f,abs(DFT_chirp));
xlabel('Frequency (Hz)');
ylabel('Magnitude of FFT');
title('FFT of Chirp Signal');
subplot(122);plot(f_shift,fftshift(abs(DFT_chirp))); 
xlabel('Frequency (Hz)');
ylabel('Magnitude of FFT');
title('Zero-Centered FFT of Chirp Signal - fftshift command');

%% 2 

% Como anteriormente constatado, uma função seno pura é representada no domínio espectral por duas
% funções Dirac, +f e -f. No entanto, neste caso tem-se um sinal chirp, em que a frequência do 
% sinal não depende somente do valor pelo qual se multiplica o tempo t no argumento do sin, 
% dependendo tambem da taxa de variação da dita frequência. 
% Para cada tempo t, a frequência da sinusoidal designa-se de frequência instantânea, sendo 
% analiticamente obtida pela derivada temporal da fase do sinal, theta (t), dando origem à seguinte 
% expressão para a frequência instantânea angular w_inst:
% w_inst = (d / dt) * theta (t) = 2 * pi * [f '(t) * t + f (t)]
% f_inst = w_inst / (2 * pi) = f (t) + f '(t) * t

% f (t) é a frequência crescente definida como: f (t) = f1 + (f2 - f1) * t / T
% Assim, a derivada de tempo de f (t) é: f '(t) = (f2 - f1) / T_bs

% Substituindo f (t) e f '(t) na expressão de f_inst, surge o seguinte:
% f_inst = f1 + 2 * (f2-f1) * t / T

% Assim, para t = T, a frequência instântanea é de f_inst = 2 * f2 - f1 = 1800 Hz
% Para t = 0, tem-se que a frequência instantanea é f_inst = f1 = 200

%O estudo analitico da frequência deste sinal chirp pode ser também consultado em anexo. 

%Deste modo, demonstra-se que as frequências instantaneas do sinal chirp
%estão contidas entre 200-1800 Hz, ao invés de 200-1000 Hz. Ao analisar o
%gráfico da FFT, denotam-se bandas para k=[200,1800] e k=[2200,3800], sendo
%tal mais uma vez explicado pela periodicidade da DFT e a sua simetria
%relativamente a d*(fs/2), com d inteiro e fs = 4000 Hz, sendo o espectro em
%k=[2200,3800] o simétrico do espectro em k=[200,1800], correspondendo a segundo banda ao 
%conjunto das primeiras harmónicas de [-1800;-200]. 

%Note-se que a DFT apenas permite obter o conteúdo em frequências do sinal,
%perdendo-se toda a informação temporal, sendo todas as frequências
%assumidas na duração do sinal incluidas no espectro do sinal.

%Mais uma vez, importa mencionar que seria expectável obter duas bandas
%retangulares compostas pelos deltas dirac positivos e negativos que
%compoem a DFT do sin para as diferentes frequências puras contidas
%no sinal, considerando o sinal chirp e a tranformada de fourier no dominio continuo.
%No entando, constata-se que a amplitude das bandas nao é
%uniforme, observando-se ondulações, designadas ondulações de Fresnel,
%mais predominantes nas extremidades das duas bandas. Tal deve-se ao facto de o suporte finito 
%computacional apenas permitir uma representação finita do sinal chirp, correspondendo no dominio 
%temporal a uma multiplicação do sinal chirp continuo por uma janela
%retangular, introduzindo descontinuidades repentinas no inicio e fim do
%sinal chirp.  Uma vez que o espectro da janela retangular 
%é a função sinc(x), observam-se lóbulos sin (x)/x em torno de cada componente 
%espectral original conforme eles são convolvidos com o espectro da janela.
%Este fenomeno designa-se spectral leakage. 

%De forma a reduzir os ripples observados, o sinal chirp pode ser multiplicado no dominio 
%temporal por uma outra janela que possibilite uma transição mais suave no inicio e fim do sinal. 

Chirp_window=Chirp.*window(@tukeywin,length(Chirp)); %multiplicação no dominio temporal pela janela
%tukeywin

DFT_chirp_window = fft(Chirp_window); 

nfft_window=length(DFT_chirp_window);
f = (0:nfft_window-1)*fs/nfft_window;
f_shift = (-nfft_window/2:nfft_window/2-1)*(fs/nfft_window); % gama de frequências centrada em zero

figure(7)
subplot(121);plot(f,abs(DFT_chirp_window));
xlabel('Frequency (Hz)');
ylabel('Magnitude of FFT');
title('FFT of Chirp Signal Convolved with Tukey Window');
subplot(122);plot(f_shift,fftshift(abs(DFT_chirp_window))); 
xlabel('Frequency (Hz)');
ylabel('Magnitude of FFT');
title('Zero-Centered FFT of Chirp Signal Convolved with Tukey Window');

%% 3)
%Calcular o espectograma para os seguintes parametros:
n_window =256; % comprimento da janela hamming 
n_overlap = 250; %número de amostras sobrepostas
nfft = 256; % número de amostras usadas para calcular a fft (nfft)
fs = 4000;

figure(8)
spectrogram(Chirp,n_window,n_overlap,nfft,fs,'yaxis');
title('Spectrogram of Chirp');

resolucao_temporal = (n_window-n_overlap)/fs; %o sinal chirp é dividido em segmentos 
%do mesmo comprimento da janela, partilhando n_overlap amostras com os
%segmentos adjacentes. Deste modo, os segmentos considerados distam
%n_window-n_overlap = 6 amostras, ou seja, 1.5 ms, definindo a resolução
%temporal 

resolucao_espectral = fs/n_window; % inverso do espaçamento temporal entre amostras 
%no segmento considerado, o qual tem o tamanho da janela.

fprintf('\n II.3. A resolução temporal é %d segundos e a resolução espectral é %d Hz\n',...
    resolucao_temporal,resolucao_espectral); 

%A variação tempo-frequência do spectrogram segue a variação imposta pelo
%chirp, verificando-se que o sinal contém frequencias contidas em
%[200,1800] Hz e as quais variam linearmente ao longo do tempo. 
    
%% 4) 
f1=1000;
f2=1500;
fs = 4000; %frequencia de amostragem 
T = 1; %duração do sinal em segundo 
t=linspace(0,T,fs*T)'; %vetor de discritização do tempo (em segundos) 

f=sin(2*pi*f1*t)+sin(2*pi*f2*t); % sinal f
g=sin(2*pi*f1*t).*sin(2*pi*f2*t); % sinal g

F=fft(f); % DFT do sinal f computada pelo algoritmo FFT 
G=fft(g); % DFT do sinal g computada pelo algoritmo FFT 

f_F = (0:length(F)-1)*fs/length(F);
f_G = (0:length(G)-1)*fs/length(G);

%plot signals
figure(9)
plot(t,f,t,g)
xlabel('Time(s)');
ylabel('Amplitude');
title('II.4) Temporal representation of f(t) and g(t)');
legend('f(t)', 'g(t)');
    
figure(10)
subplot(121)
plot(f_F,abs(F));
xlabel('Frequency(Hz)');
ylabel('Magnitude of the FFT');
title('II.3: FFT do sinal f')
axis tight

subplot(122)
plot(f_G,abs(G));
xlabel('Frequency(Hz)');
ylabel('Magnitude of the FFT');
title('II.3: FFT do sinal G')
axis tight

% O espetro do sinal f reflte a propriedade linear da DFT, 
% em que o espetro da soma de dois sinais resulta na soma dos
% espetros de cada sinal. Como se pode ver pela dedução feita em anexo, a DFT de f vai ter diracs
% em f1=1000 Hz e f2=1500 Hz. Como é periódica e simétrica relativamente a d*(fs/2), 
% com d inteiro e fs = 4000 Hz, os diracs repetem-se em
% 4000-f1=3000 Hz e 4000-f2=2500 Hz, correspondendo às primeiras harmónicas
% de -1000 Hz e -1500 Hz do espectro original.

% Em relacao ao sinal g, o produto dos dois sinais no dominio temporal corresponde à
% convolução no dominio espectral, pelo que o espectro resultante é a convolução dos
% espetros de cada um dos senos. Analiticamente, derivou-se que g equivale 
%à soma de 2 sinusóides, mas de frequências que diferem de f1 e f2:
%g(t) = (1/2)*cos(2*pi*500*t) - (1/2)*cos(2*pi*2500*t)
%Visto a frequência da segunda sinusoide ser superior à frequência de
%Nyquist (2500>fs/2 = 2000), ocorre aliasing. Para a sinusóide de 500 Hz, 
%obtém-se no espectro um delta dirac em 500 Hz, obtendo-se a 1ºharmónica de
%-500Hz em -500+4000 = 3500 Hz. Devido ao aliasing, é representado
%um delta dirac em 2000-(2500-2000) = 1500 Hz, repetindo-se em 4000-1500 =
%2500 Hz. 

%Novamnete refletindo a lineariedade da DFT, quanto à magnitude da FFT constata-se que as
%amplitudes dos picos de G correspondem aproximadamente a metade da
%amplitude dos picos de F. 

    %% Part III
    
    %% 1)
    
f=100; % Frequência do sinal x(t) em Hz
f0=1000; % Frequência da onda portadora, c(t), em Hz
fs=2400; % Frequencia de amostragem em Hz
 
t=linspace(0,T,fs*T)';
x=sin(2*pi*f*t); % Sinal x(t)
X=fft(x); % DFT do sinal f computada pelo algoritmo FFT 

nfft = length(X);
f = (0:nfft-1)*fs/nfft;
f_shift = (-nfft/2:nfft/2-1)*(fs/nfft); % gama de frequências centrada em zero


figure(11);
subplot(121);plot(f,abs(X));
xlabel('Frequencies (Hz)');
ylabel('Magnitude of FFT');
title('FFT of  x(t)');
subplot(122);plot(f_shift,fftshift(abs(X)));
xlabel('Frequencies (Hz)');
ylabel('Magnitude of FFT');
title('zero-centered FFT of x(t) - fftshift command');

%Como referido anteriormente, a DFT duma função seno corresponde a dois deltas
%dirac para a frequência fundamental (100) e a sua simétrica (-100).
%Como a dft assume uma extensão periódica do sinal ao infinito em ambas as
%direções, considerando apenas frequências positivas, o 2º delta dirac que
%se obtém corresponde à primeira harmónica de -100Hz no espectro original, correspondente a
%2400-100 = 2300 Hz. 

%% 2)

% A Modulação em Amplitude é um sistema onde a amplitude máxima da onda portadora varia de acordo 
% com o valor instantâneo (amplitude) do sinal modulante (mensagem ou banda base).

list_alpha = logspace(0,-3,4); % logspace gera 4 pontos entre 10^0=1 and 10^-3

num_alpha = length(list_alpha);

c=cos(2*pi*f0*t); % Onda portadora, c(t)

figure(12)
for i = 1:num_alpha
    alpha = list_alpha(i);
    xam=(1+alpha*x).*c; % sinal de amplitude modulada
    subplot(4,2,2*i-1);plot(t, xam);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['alpha = ', num2str(alpha)]);
    
    Xam=fft(xam);  % DFT do sinal f computada pelo algoritmo FFT 
    nfft = length(Xam);
    f = (-nfft/2:nfft/2-1)*(fs/nfft); % gama de frequências centrada em zero
    subplot(4,2,2*i);plot(f,fftshift(abs(Xam)));
    xlabel('Frequencies (Hz)');
    ylabel('Magnitude of FFT');
    title(['alpha = ', num2str(alpha)]);
    axis tight
end

% Desenvolvendo a expressão xam(t,alpha) = (1 + alfa*x(t))*c(t), tem-se
% que: xam(t,alfa) = cos(2*Pi*f0*t) + alfa*x(t)*cos(2*Pi*f0*t), ou seja
% equivale à soma da onda portadora com o produto do sinal
% alfa*x*portadora.
% Desta forma, tendo em conta as proprieddaes da DFT, o espectro obtido
% será dado pela soma do espetro da onda portadora com o espectro
% resultante da convolução entre os espectros X e C. 

% Em anexo demonstra-se que x_am equivale a:  
% x_am(t,alpha) = cos(2*pi*f0*t) + (alpha/2)*(sin(2*pi*(f0+f)*t) + sin(2*pi*(f0-f)*t)

%Considerando a lineariedade da DFT, o espectro vai resultar da
%soma dos espectros do cosseno (dois deltas dirac em +-f0 = +-1000Hz), 
%com os espectros dos senos (dois deltas dirac em +-(f0+f)=+-1100Hz e em 
% +-(f0-f)=+-900Hz), apresentando os espectros dos senos uma magnitude de
%alpha/2 da sua original, pelo que quanto maior o alfa, maior a magnitude 
%associada aos deltas dirac dos espectros dos dois sinais seno. Desta
%forma, quanto maior o valor de alfa, uma maior quantidade de informação do
%sinal x é transmitida pela onda portadora. 

% Comparando o espectro de xam com o espectro de x, tem-se que o o espectro de xam corresponde 
% ao delta dirac na frequência f0 da onda portadora sobreposto ao espectro
% do sinal x centrado na frequência f0 da onda portadora, devido à
% convolução.

%% 3)

% FM em vez de modular a magnitude como antes, modula a frequência. f0 é a frequência 
% da onda portadora que mudará de acordo com o valor de x(t).

f0=1000; % Frequência da onda portadora, c(t), em Hz
fs=2400; % Frequencia de amostragem em Hz

figure(13)
for i = 1:num_alpha
    alpha = list_alpha(i);
    c = cos(2*pi*f0*(1+alpha*x).*t);
    subplot(4,2,2*i-1);plot(t, c);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['alpha = ', num2str(alpha)]);

    C=fft(c); % DFT do sinal f computada pelo algoritmo FFT
    nfft = length(C);
    f = (-nfft/2:nfft/2-1)*(fs/nfft); % gama de frequências centrada em zero
    subplot(4,2,2*i);plot(f,fftshift(abs(C)));
    xlabel('Frequencies (Hz)');
    ylabel('Magnitude of FFT');
    title(['alpha = ', num2str(alpha)]);
    axis tight
end

% A onda portadora c modulada em frequência, à semelhança do sinal chirp, 
% tem frequências que variam a uma dada taxa, sendo que no dominio
% continuo, para um dado instante t associa-se uma frequência instantanea
% dada por (derivado em anexo):
% f_inst = f0 + f0*alpha*(2*pi*f*cos(2*pi*f*t)*t + sin(2*pi*f*t))

% Analisando o intervalo de frequências instantaneas contidas no sinal, 
% obtidas para os instantes temporais extremos, obtem-se 
% para alpha=0.1 f_inst=[1000;63835.85], e para para alpha=1
% f_inst=[1000;629318.53]. Deste modo, a banda de frequências inclui valores que 
% estão significativamente acima da frequência de Nyquist (fs/2 = 1200 Hz),
% sendo que para essas frequências o sinal sofre aliasing. Assim, o espectro
% desses sinais é a superposição do espectro puro com as suas harmónicas,
% devido à periodicidade da DFT. 

% Enquanto na modulação em amplitude (AM) obtêm-se deltas dirac em frequências definidas, 
% na modulação em frequência (FM), visto a frequência da onda portadora variar em função do 
% sinal a transmitir, obtem-se uma maior ocupação da banda de frequência.
% Constata-se também que a onda modulada em frequência assemalha-se à
% modulada em amplitude quando alfa é diminuido, contudo ter-se-á ruido, tendo-se a situação
% limite quando alfa = 0 e nao se transmite informação relativa a x(t). De
% forma a aumentar a qualidade da transmissao, aumenta-se o alfa,
% resultando numa lagura de banda muito elevada.  

% Realizados os espectrogramas da onda portadora modulada em frequência por x(t), para os 
% alphas considerados, os mesmos confirmam a ampla faixa de frequências presente a cada instante.
for i = 1:num_alpha
    alpha = list_alpha(i);
    c = cos(2*pi*f0*(1+alpha*x).*t);
  
    %Calcular o espectograma para os seguintes parametros:
    n_window =256; % comprimento da janela hamming 
    n_overlap = 250; %número de amostras sobrepostas
    nfft = 256; % número de amostras usadas para calcular a fft (nfft)
    fs = 2400;

    figure(13+i)
    spectrogram(c,n_window,n_overlap,nfft,fs,'yaxis');
    title(['Spectrogram of C - FM Modulation for alpha = ', num2str(alpha)]);
end

%spectrogram(cos(2*pi*f0*(1+0*x).*t),n_window,n_overlap,nfft,fs,'yaxis');

%% 4) 

% Para a modulação em amplitude (AM), a informação a transmitir x(t)
% encontra-se codificada na amplitude da onda portadora, com a fase e a
% frequência a permanecerem constantes. Como verificado, quanto menor o
% alfa, maior a dependência da amplitude do sinal modulado com a informação 
% que se pretende transmitir, possibilitando uma melhor reconstrução do sinal
% x(t).
% Relativamnete à modulação em frequência (FM), a frequência da onda
% portadora varia com a informação a transmitir x(t), permanecendo a
% amplitude e a fase constantes. 

%Visto o ruido afetar a amplitude da onda modulada, os sinais modulados em FM são
%menos afetados, pois o ruido nao afeta a frequência da onda modulada, 
%ao contrário da modulação em AM, em que a informação
%transmitida na presença de ruido tem pouca qualidade, sendo por isso um
%método de modulação menos eficaz que FM. 
 
% Como tal, infere-se que a FM é melhor para transmitir informação
% devido ao facto de possibilitar uma melhor reconstrução do sinal x(t). Contudo, 
% como a qualidade da transmissão em FM implica que a largura de banda da
% onda modulada seja grande, a distância de transmissão é condicionada, isto é, em FM não se 
% pode transmitir a longas distância, sendo que por isso é necessário existirem mais antenas 
% menos espaçadas de forma a que os obstaculos fisicos nao afetem a transmissao. Em contrapartida, 
%visto AM possibilitar uma menor largura de banda, permite ter mais estações disponíveis em 
%qualquer gama de frequência.
