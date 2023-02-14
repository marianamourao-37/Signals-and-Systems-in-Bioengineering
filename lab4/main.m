%Grupo 22: 
%-Ana Rita Lopes nº98587
%-Mariana Mourão nº98473
 
%LAB#4 - Filtering
 
clear all;
close all;
clc;
 
%% Parte 1 - 1D FILTERING
    %% 1) Função definida no script designado 'filtro.m'
    
% Em anexo, demonstra-se a transformada Z inversa da função de 
% transferência do filtro H(z), obtendo-se a resposta impulsional h(n).
    
    %% 2)

% Implementação de um filtro Butterworth digital lowpass de 10ª ordem, com frequência de
% corte w = pi/2. Para tal, recorreu-se à função built-in do matlab
% butter, determinando-se os coeficientes a e b do filtro. w é dada
% em radianos/amostra (rad/sample) e é normalizada entre 0 a pi. Como tal, 
% o filtro low-pass requerido tem uma frequência de corte normalizada wn = 1/2

[b,a] = butter(10, 1/2);

% Reposta em frequência do filtro, apresentando a sua magnitude e fase 
figure(1)
freqz(b,a);
title('Frequency Response of a 10th-order Butterworth filter (wn = 1/2)');

% Observando-se o diagrama de bode para a magnitude da resposta em
% frequência do filtro, podem-se caracterizar três regiões:
% - passband: faixa de frequência em que o filtro passa a energia do sinal. Normalmente 
% definido como a faixa de frequência em que a resposta em frequência do filtro é igual ou 
% superior a -3 dB. Verifica-se que a resposta em frequência é constante (0 dB) na banda passante, 
% correspondendo às frequências inferiores à frequência de corte, confirmando-se
% de que se trata de um filtro passa baixo. Constata-se também que a banda passante nao apresenta ondulações 
% ou ripples, as quais, se existentes, iriam conduzir à distorção dos sinais que ocorrem na banda passante. 
% - stopband: banda de frequências, entre limites especificados, em que não é permitida a passagem de sinais ou 
% a atenuação está acima do nível de atenuação da banda de paragem necessário. 
% - banda de transição: limite entre a passband e a stopband é a banda de transição, que ocorre em torno da 
% frequência de corte do filtro. A largura (declive da resposta entre passband e stopband) desta banda de transição é o
% roll-off do filtro, verificando-se que a mesma é pouco abrupta. Para que a zona de 
% transição seja mais abrupta, temos de aumentar a ordem do filtro, contudo tal implica 
% um maior nº de coeficientes e, como tal, maior tempo de computação.

% Relativamente ao diagrama de bode para a fase da resposta em frequência do filtro, conclui-se que 
% a variação da fase com a frequência é aproximadamente linear longe da 
% frequencia de corte, contudo é altamente não linear na zona da frequência
% de corte. O observado sugere que ocorrerá distorção de fase do sinal à
% entrada, atrasando diferentemente as diversas componentes em frequência do sinal, sendo que o 
% maior atraso temporal ocorre sensivelmente para a frequência de corte do filtro. Para compensar a
% distorção de fase introduzida pelo filtro, o matlab tem uma função
% built-in designada filtfilt, a qual permite realizar zero-phase
% filtering, procedendo à filtragem sequencial do sinal de entrada em ambas
% as direções, direta e reversa. O resultado é uma sequência filtrada que tem distorção de fase zero e 
% uma ordem de filtro que é duplicada, sendo por isso a atenuação do sinal de entrada também
% duplicada. 

polos = roots(a); %localização dos polos
zeros = roots(b); %localização dos zeros

figure(2)
zplane(zeros,polos);
grid
title('zero-pole plot of 10th-order Butterworth filter (wn = 1/2)');
legend('Zeros','Polos');

% Relativamente à distribuição dos polos e zeros no plano Z, tem-se que existem 
% 10 pólos (ordem 10) distribuídos ao longo do eixo imaginário igualmente espaçados, 
% correspondendo ao esperado para um filtro Butterworth passa-baixo 
% discreto, sendo que para um filtro Butterworth contínuo espera-se que os 
% polos ocorram na circunferência de raio unitário igualmente espaçados. 
% Quanto aos zeros, constata-se que se distribuem no semi-plano complexo
% esquerdo, em redor de z=-1. 
% Uma vez que todos os pólos se encontram localizados dentro do círculo unitário, 
% tal constitui uma condição suficiente para o filtro ser considerado
% estável. 


% os filtros IIR, por serem recursivos, são altamente sensíveis à 
% quantificação dos seus coeficientes, podendo originar elevados erros 
% de arredondamento e, consequentemente, tornar-se instável. Tal é 
% particularmente critico em caso de filtros de elevada ordem, visto que 
% para um filtro de ordem n existem n coeficientes a ser quantificados, e 
% os quais, para os filtros IIR, dependem da quantificação dos anteriores. 
% Nesse sentido, ao invés de implementar um filtro de elevada ordem, implementa-se 
% uma cascata de filtros de ordens inferiores, nomeadamente os designados Second 
% Order Sections (SOS), que permitem conjugar raízes (polos ou zeros) complexas, 
% ao contrário de filtros de 1ªordem
[z,p,k] = butter(10,1/2); 
sos = zp2sos(z,p,k);
figure(3)
freqz(sos) 

figure(4)
zplane(z,p);
grid
title('zero-pole plot of 10th-order Butterworth filter (wn = 1/2) - Second Order Sections');
legend('Zeros','Polos');

% Analisando a distribuição dos polos e zeros no plano Z para o SOS, tem-se
% que verifica-se a mesma distribuição dos polos no eixo imaginário, com os
% zeros a estarem sobrepostos e lozalizados na circunferencia de raio unitario, em z =
% -1. 

    %% 3. 
        %% a) 

%ler os ficheiros de audio fornecidos, retornando sampling data (sd_awgn e sd_saltpeper) 
%e uma taxa de amostragem (em hertz) para os ficheiros correspondentes.

[sd_awgn, fs_awgn]=audioread('fuggesAWGN.mp3');
[sd_spn, fs_spn]=audioread('fuggesSaltPeper.mp3');

%ouvir cada ficheiro de audio:
%soundsc(sd_awgn, fs_awgn);
%soundsc(sd_saltpeper, fs_saltpeper);

% numero de amonstras de cada sinal digitalizado 
N_awgn = length(sd_awgn); 
N_spn = length(sd_spn);  

% vetores de tempo discretizados, segundo a respetiva frequência de
% amostragem e número de amostras
t_awgn = ((0:N_awgn-1)/fs_awgn)';
t_spn = ((0:N_spn-1)/fs_spn)'; 

% Representação gráfica dos sinais de aúdio digitalizados, para uma melhor comparação visual entre
% os mesmos
figure(5)
subplot(3,1,1), plot(t_awgn, sd_awgn);
title('Audio corrupted by Additive white Gaussian Noise');
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
axis tight

subplot(3,1,2), plot(t_spn, sd_spn);
title('Audio Signal corrupted by Salt Pepper Noise');
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
axis tight

        %% b) 

%Para uma melhor comparação entre os ruídos que corrompem a musica, 
%decidiu-se avaliar o sinal original, nao corrompido, fugges.wav

[sd_original, fs_original] = audioread('fugges.wav');
%soundsc(sd_original, fs_original);

N_original = length(sd_original);

t_original = ((0:N_original-1)/fs_original)'; 

subplot(3,1,3), plot(t_original, sd_original);
title('Uncorrupted Original Audio Signal');
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
axis tight  

% Ao comparar-se o sinal original com o corrompido pelo ruido branco
% gaussiano aditivo (AWGN), constata-se uma mancha densa central que 
% afeta todo o sinal original. Ouvindo-se os áudios respetivos, denota-se 
% um ruído de fundo a corromper o sinal original, distorcendo todas as amostras 
% do sinal, tornando-o menos percetivel. Ora, o constatado é
% devido às caracteristicas do AWGN: aditivo, pelo que é adicionado ao sinal original, 
% podendo ser relativo a ruído intrínseco ao sistema; branco, pelo que todas as amostras de ruído 
% são independentes e não correlacionadas entre si, apresentando uma densidade 
% espectral de potência constante em todas as frequências; gaussiano, pelo que as 
% amostras de ruído seguem uma distribuição normal (média zero e variância finita) 
% ao longo do tempo. Desta forma, se se subtrair o sinal original ao
% sinal corrompido, idealmente obter-se-ia uma quantificação do AWGN adicionado.

% Quanto ao sinal corrompido pelo ruido salt and peper (SPN), constata-se
% que o mesmo apresenta amostras corrompidas e outras que permanecem 
% conservadas (não corrompidas), distinguido-se picos aleatórios no dominio temporal.
% Ora, o constatado é devido às caracteristicas 
% do SPN, designado ruido "impulsivo", o qual perturba o sinal em amostras
% aleatórias segundo a mesma magnitude em módulo, originando saturações positivas 
% (magnitude +1) e/ou negativas (magnitude -1). Se subtraido o sinal original ao sinal 
% corrompido, obter-se-iam impulsos aleatoriamente distribuidos de magnitude -1 e +1. Ouvindo-se o áudio
% respetivo, constata-se  que o mesmo é perturbado por "estalos" ou
% "cliques", os quais correspondem às saturações mencionadas. Se a
% densidade do ruído fosse superior, a musica tornar-se-ia progressivamente
% mais irreconhecivel.

        %% c) 

% Filtrar os sinais corrompidos com o filtro concebido em 2.
filter_sd_awgn = filtro(sd_awgn, a', b');
%soundsc(filter_sd_awgn,fs);
filter_sd_spn = filtro(sd_spn, a', b');
%soundsc(filter_sd_spn,fs);

% Representação gráfica do sinal corrompido por AWGN, sinal corrompido filtrado e 
% do sinal correspondente à subtração entre os dois 
figure(6);
subplot(3,1,1);
plot(t_awgn,sd_awgn);
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal AWGN');

subplot(3,1,2);
axis tight
plot(t_awgn,filter_sd_awgn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Linearly filtered signal AWGN');

subplot(3,1,3);
plot(t_awgn,sd_awgn - filter_sd_awgn)
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Difference between signals');

% Representação gráfica do sinal corrompido por SPN, do sinal filtrado, e
% do sinal correspondente à subtração entre os dois 
figure(7);
subplot(3,1,1);
plot(t_spn,sd_spn);
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal SNP');

subplot(3,1,2);
plot(t_spn,filter_sd_spn);
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Linearly filtered signal SNP');

subplot(3,1,3);
plot(t_spn,sd_spn - filter_sd_spn)
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Difference between signals');

% Enquanto análise geral, constata-se que o filtro linear conduziu à
% redução da amplitude do ruído de ambos os sinais corrompidos, sendo
% audível essa redução. Por outro lado, também se verifica que os sinais
% filtrados apresentam um tom mais grave, uma vez terem sido
% submetidos a uma filtragem passa-baixo, levando a que a informação da
% música correspodente a altas frequências (nomeadamente superiores à
% frequência de corte do filtro) - som mais agudo - seja atenuada e, em
% parte, perdida.

% Analisando mais especificamente as representações gráficas no dominio 
% temporal para o sinal corrompido por AWGN e o sinal filtrado, conclui-se 
% que a filtragem linear efetuada nao foi particularmente eficiente,
% visualizando-se ainda a macha estreita em torno de 0 que afeta o sinal 
% em toda a sua extensão, contudo agora com uma menor magnitude. Isto
% deve-se ao facto de o filtro calcular cada amostra do sinal filtrado como
% uma combinação linear, ponderada por diferentes pesos, do sinal à saida y
% previamente reproduzido e do sinal à entrada x até ao momento filtrado. 
% Uma vez que o ruido já está distribuido em toda a extensão do sinal 
% (ruido aditivo), a combinação linear efetuada pelo filtro conduz à
% atenuação da magnitude do ruido. Possivelmente, o ruido poderia ser mais 
% eficazmente removido caso se implementasse um filtro Butterworth mais
% adequado para o sinal corrompido em causa. Contudo, é de notar que visto 
% a largura de banda do AWGN ser teoricamente infinita, é de esperar que 
% o mesmo nunca seja efetivamente removido. 

% Quanto às representações gráficas no dominio temporal para o sinal 
% corrompido por SPN e o sinal filtrado, constata-se que pouca diferença
% existe entre os mesmos, pelo que se conclui que a remoção de SPN foi
% ineficaz com o filtro linear Butterworth concebido. Visto a filtragem 
% passa-baixo sensivelmente corresponder a uma média do sinal numa dada 
% janela deslizante, tem-se que o mesmo é sensivel e altamente influenciado por 
% valores extremos (outliers), os quais são produzidos pelo SPN na
% intensidade do sinal. Como o AWGN afeta a intensidade do sinal segundo
% uma média de 0 (distribuição normal gaussiana), conclui-se que o processo
% de filtragem é mais eficiente para este ultimo caso em comparação com a
% filtragem de SPN. No entanto, visto o filtro aplicado ser um passa-baixo,
% verifica-se que as altas frequências (superiroes à frequência de corte)
% são atenuadas e, em parte, removidas, sendo mais notória a atenuação dos
% picos de saturação. Contudo, visto o filtro ser linear e, como tal, 
% cada amostra do sinal filtrado resultar numa combinação linear tanto 
% de entradas como de saídas anteriores, tem-se que o ruido dispersa-se por 
% todas as amostras, espalhando-o por instantes que originalmente não tinham 
% ruido.

        %% d) 
        
N=3; % dimensão da janela
% Filtragem dos sinais corrompidos com o filtro de mediana 
filtmed_sd_awgn = medfilt1(sd_awgn, N);
filtmed_sd_spn = medfilt1(sd_spn, N);

% Representação gráfica do sinal original, do sinal corrompido por AWGN, e do sinal filtrado 
figure(8);
subplot(3,1,1), plot(t_original, sd_original);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Uncorrupted Original Audio Signal');
axis tight

subplot(3,1,2);
plot(t_awgn,sd_awgn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal AWGN');
axis tight

subplot(3,1,3);
plot(t_awgn,filtmed_sd_awgn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original Corrupted signal AWGN with Median filter, window 3');
axis tight

% Representação gráfica do sinal original, do sinal corrompido por SPN, e do sinal filtrado 
figure(9);
subplot(3,1,1), plot(t_original, sd_original);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Uncorrupted Origonal Audio Signal');
axis tight

subplot(3,1,2);
plot(t_spn,sd_spn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal SPN');
axis tight

subplot(3,1,3);
plot(t_spn,filtmed_sd_spn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original Corrupted signal SPN with Median filter, window 3');
axis tight

% O filtro de mediana é um filtro não linear que computa cada amostra n do 
% sinal filtrado de saida como a mediana das amostras do sinal 
% à entrada dentro de uma janela de tamanho N centrada na amostra n 
% (mesma posição) do sinal de entrada, correspondendo a uma abordagem 
% computacionalmente pesada. Contudo, visto a mediana não ser tão sensivel 
% a outliers como a média, tem-se que esta filtragem é mais eficiente para
% a remoção de SPN, ao suprimir-se os valores extremos nas janelas
% consideradas, substituindo-os por amostras não corrompidas, sendo de
% notar que não correspondem à amostra original, mas a uma amostra não
% corrompida vizinha (vizinhança de 1 amostra, visto a dimensão da janela 
% ser N = 3) da original, pelo que se pode considerar que estejam
% correlacionadas. Contudo, analisando a representação gráfica do sinal 
% corrompido por SPN filtrado, constata-se que a filtragem não foi tão
% eficiente como seria de esperar, podendo dever-se ao facto de o desempenho do
% filtro ser altamente dependente da densidade do ruído do SPN
% e do tamanho da janela: se a densidade for muito alta ou a janela muito 
% pequena, a remoção de ruído não é eficaz. 

% Quanto à filtragem do ruido AWGN, a mediana de um sinal com ruido aditivo
% branco (ruido afeta toda a extensão do sinal) corresponde a uma amostra
% ruidosa, ao passo que no SPN apenas algumas amostras aleatórias sao
% corrompidas, pelo que podem ser suprimidas. Desta forma, a remoção do 
% AWGN nao é eficiente, permanecendo o ruido de fundo, podendo a magnitude
% do ruido ser diferente. Deste modo, conclui-se que a filtragem linear passa-baixo 
% do AWGN é mais eficiente que a filtragem não linear dada pela mediana.

        %% e) 

% Representação gráfica do sinal corrompido por AWGN
figure(10)
subplot(4,1,1)
plot(t_awgn,sd_awgn);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal AWGN');
axis tight

% Representação gráfica do sinal corrompido por SPN
figure(11)
subplot(4,1,1)
plot(t_spn,sd_spn);
axis tight
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Original corrupted signal SPN');

list_N = 10:10:30; % vetor com as dimensões da janela a serem testadas 

for i = 1:length(list_N)
    % Filtragem dos sinais corrompidos com o filtro de mediana 
    filtmed_sd_awgn = medfilt1(sd_awgn, list_N(i));
    filtmed_sd_spn = medfilt1(sd_spn, list_N(i));
    
    % Representação gráfica do sinal corrompido por AWGN filtrado
    figure(10)
    subplot(4,1,1+i)
    plot(t_awgn,filtmed_sd_awgn);
    axis tight
    xlabel('Time (seconds)');
    ylabel('Signal Amplitude');
    title(['Median filtered signal AWGN for a window size of ',...
        num2str(list_N(i))]);
    
    % Representação gráfica do sinal corrompido por SPN filtrado
    figure(11)
    subplot(4,1,1+i)
    plot(t_spn,filtmed_sd_spn);
    axis tight
    xlabel('Time (seconds)');
    ylabel('Signal Amplitude');
    title(['Median filtered signal SNP for a window size of ',...
        num2str(list_N(i))]);
end

% Representação gráfica do sinal original
figure(12)
plot(t_original, sd_original);
xlabel('Time (seconds)');
ylabel('Signal Amplitude');
title('Uncorrupted Original Audio Signal');
axis tight

% Após testar várias dimensões da janela do filtro de mediana aplicado aos
% sinais corrompidos, verifica-se que a partir duma dada dimensão da janela
% a qualidade do sinal é deteriorada, sendo particularmente visivel através
% da redução da amplitude das amostras iniciais dos sinais, ainda que o 
% mesmo efeito seja reproduzido em toda a extensão dos sinais 
% (menos percetivel). Auditivamnete, o referido traduz-se em sinais com tons
% mais graves (remoção das frequências mais elevadas) à medida que a dimensão 
% da janela aumenta. Tal deve-se ao facto de janelas de dimensões 
% progressivamente maiores conduzirem ao aumento da probabilidade de a amostra n ser substituida 
% por uma amostra menos correlacionada (mais distante temporalmente) com a 
% original. 

% Analise do processo de filtragem dos dois tipos de ruido consoante a 
% dimensão da janela, para valores inferiores à dimensão critica em que 
% ocorre um deterioramento progressivo da qualidade do sinal devido à
% substituição por amostras não correlacionadas: 
% -> Analisando mais especificamente a filtragem do AWGN consoante a dimensão
% da janela, tem-se que o processo de filtragem de AWGN produz sensivelmente 
% o mesmo sinal filtrado independemente da dimensão da janela, não 
% correspondendo a uma filtragem eficiente visto a mediana de um sinal com 
% ruído branco aditivo (ruido distribuido em toda a extensão do sinal) ser 
% uma amostra com ruído. 
% -> Relativamente à filtragem do SPN consoante a dimensão da janela, obtem-se
% uma filtragem mais eficiente à medida que a dimensão da janela aumenta, 
% com a maioria dos picos de elevada intensidade a serem suprimidos. Para
% janelas demasiado estreitas, pode existir a possibilidade de duas
% amostras consecutivas estarem corrompidas, não permitindo uma eficaz
% remoção do ruido aleatório SPN. 

% Como tal, identifica-se um trade-off entre a dimensão da janela e a
% qualidade do sinal. 

        %% f) 

% Sumariando todas as conclusões retiradas dos exercicios anteriores,
% tem-se que o filtragem linear pela aplicação do filtro de Butterworth 
% apenas foi considerado razoável para remover AWGN, visto computar cada
% amostra do sinal filtrado como uma combinação linear das saidas e
% entradas anteriores, conduzindo à atenuação do ruido. Contudo, para a
% remoção de SPN verificou-se que a filtragem linear conduz à dispersão do
% ruido por toda a extensão do sinal, adicionando ruido a amostras que 
% previamente nao apresentavam, pelo que não é eficaz na sua remoção.

% Já relativamente à filtragem não-linear aplicada pelo filtro de mediana, 
% conclui-se que é mais eficiente na remoção de SNP, visto computar cada 
% amostra do sinal filtrado como a mediana das amostras do sinal à entrada 
% compreendidas numa janela centrada na mesma posição no sinal de entrada, 
% conduzindo à supressão de outliers, os quais pouco influenciam a
% mediana. Contudo, para a remoção de AWGN verificou-se que a filtragem
% não-linear oferece poucas melhorias, visto a mediana de um sinal
% corrompido por ruido branco aditivo ser uma amostra ruidosa. 


%% Parte II - 2D FILTERING

    %% 1) 
    
I_original = imread('len_gray.jpg'); % Carrega a imagem a ser analisada

I_original = mat2gray(I_original); % conversão da imagem para 
% uma matriz de intensidades, contendo valores no intervalo de 0 (preto) a 
% 1 (branco).

figure(13);
subplot(1,2,1); 
colormap gray 
imagesc(I_original);
title('Original image')

subplot(1,2,2);
mesh(I_original); 
colorbar
title('Original image (mesh)')

list_noise = logspace(-2,0,5);
% vetor com a quantidade de noise a ser aplicado. logspace gera 
% 5 pontos entre 10^a e 10^b, em que a=-2 e b=0.

for i = 1:length(list_noise)
    %Imagem corrompida para uma determinada densidade de ruido
    I_corrupted = imnoise(I_original,'salt & pepper',list_noise(i)); 
    
    figure(13+i)
    subplot(1,2,1); colormap gray; imagesc(I_corrupted);
    title(['Noisy Image (Salt and Pepper) corrupted by noise = ', ...
        num2str(list_noise(i))]);
    
    subplot(1,2,2); 
    mesh(I_corrupted);
    colorbar
    title(['Noisy Image (mesh) corrupted by noise = ', ...
        num2str(list_noise(i))]);
end

% O SPN caracteriza-se por substituir aleatoriamente, consoante a densidade 
% do ruido, um determinado numero de pixeis da imagem original por pixeis 
% brancos (saturação positiva, com máxima intensidade) ou por pixeis 
% pretos (saturação negativa, com minima intensidade). Para as densidades
% de ruido testadas, verificou-se que para 0.01 obteve-se a imagem menos
% corrompida, apenas com alguns pixeis a serem saturados, permitindo ainda
% percecionar facilmente a imagem. Contudo, para a densidade de ruido de 1,
% verifiou-se que a imagem resultante é constituida por um conjunto de
% pixeis brancos e pretos, nao permitindo percecionar a informação da imagem 
% original devido ao elevado ruido (afeta 100% das amostras).  
% Como tal, quanto maior a densidade de ruido, um maior numero de pixeis serão
% saturados e menos percetiveis são as imagens. O comando mesh, ao
% representar tridimensionalmente as imagens, isto é, para cada pixel
% apresenta a sua magnitude/intensidade enquanto a altura, permite analisar
% a saturação dos pixeis na imagem com o ruido adicionado

    %% 2) 

% Imagem corrompida para uma densidade de ruido de 0.05, afetando
% aproximadamente 5% dos pixeis
I_corrupted = imnoise(I_original,'salt & pepper');

pow2 = 6:9;  % vetor de potências
N = 2.^(pow2); % vetor de potências de base 2 

% Correspondência entre o fator de largura α e o desvio-padrão σ da função 
% de densidade de probabilidade gaussiana é dada por:
% σ = (N – 1)/(2α), com N como sendo o comprimento do filtro 

figure(19)
% Diferentes comprimentos da janela, para um mesmo desvio-padrao de ~2.5, pelo 
% que α = 0.2*(N-1) 
for i = 1:length(N)
    alpha = (1/5)*(N(i)-1);
    w_gaussian = gausswin(N(i), alpha); % janela gaussiana 
    n = -(N(i)-1)/2:(N(i)-1)/2;
    
    % representação grafica da janela gaussiana 1D no dominio discreto 
    subplot(2,4,i)
    plot(n,w_gaussian)
    title(['Square Gaussian Window [', num2str(N(i)), ', ', ...
        num2str(N(i)), ']']);
    axis tight

    % filtragem 2D por convolução no dominio discreto, aplicando-se o filtro 
    % em ambas as dimensões da imagem 
    I_gaufilt = conv2(w_gaussian,w_gaussian',I_corrupted,'same');
    subplot(2,4,i+4)
    colormap gray; imagesc(I_gaufilt);
    title('Gaussian Filtered SPN Image')
end

figure(20)
n = -(N(1)-1)/2:(N(1)-1)/2;

% Diferentes desvios-padrao, para um mesmo comprimento da janela N(1) = 64
for i = 1:length(N)
    alpha = (1/5)*(N(i)-1);
    w_gaussian = gausswin(N(1),alpha); % janela gaussiana 

    subplot(2,4,i);
    plot(n,w_gaussian);
    title(['Square Gaussian Window with stdev = ', ...
        num2str(round((N(1)-1)/(2*alpha),2))]);
    axis tight

    % filtragem 2D por convolução no dominio discreto, aplicando-se o filtro 
    % em ambas as dimensões da imagem 
    I_gaufilt = conv2(w_gaussian,w_gaussian',I_corrupted,'same');
    subplot(2,4,i+4)
    colormap gray; imagesc(I_gaufilt);
    title('Filtered image')
end

% resposta impulsional do filtro gaussiano 2D 
h_g = fspecial('gaussian',N(1),2.5); % N = 64, σ = 2.5

figure(21);
mesh(h_g);
colormap(hot); caxis([0 2]);
title('Impulse Responde of a 2D gaussian  - Discrete Domain')
xlabel('N','FontSize',10); ylabel('M','FontSize',10);
zlabel('Amplitude','FontSize',10);

figure(22)
mesh(fftshift(abs(fft2(h_g)))); title('Surface Plot da resposta em frequencia do filtro gaussiano');
colormap(hot); caxis([0 2]);
xlabel('Horizontal Frequency','FontSize',10); ylabel('Vertical Frequency','FontSize',10);
zlabel('Magnitude','FontSize',10);
% No dominio espectral, o filtro gaussiano aplicado também corresponde a
% uma gaussiana, com um desvio padrão proporcional. 

% Matematicamente, a convolução da imagem com a função gausiana no dominio 
% discreto corresponde a uma filtragem passa-baixo, com a imagem filtrada
% traduzida por y(n,m)=sum(i,j=1)h(i,j)*x(n-i,m-j), com a intensidade do
% pixel (n,m) a corresponder a uma média ponderada dos pixeis vizinhos, 
% com pesos que diminuem monotonicamente com a distância ao pixel central 
% (n,m) do kernel/mascara. 
% Como tal, a atenuação do ruido no pixel central da janela dá mais peso
% aos pixeis mais proximos e menos peso aqueles que estão mais distantes,
% sendo esta uma propriedade importante de forma a que não ocorra distorção
% da imagem por se dar mais relevância a pixeis menos correlacionados (mais
% distantes). Deste modo, a filtragem gaussiana 2D aplicada a píxeis com 
% ruído resulta na sua atenuação devido à contribuição dos pixeis vizinhos 
% não corrompidos. Por outro lado, pixeis sem ruído são afetados pelos
% pixeis vizinhos que tenham ruído.

% Analisando o impacto dos parâmetros da janela gaussiana, tem-se que a sua 
% dimensão, para um mesmo desvio-padrão, não parece produzir diferenças
% significativas observáveis na remoção de SPN. 

% Contrariamente, para uma mesma dimensão da janela, mas diferentes desvios padrão, verifica-se uma
% performance diferenciada da filtragem, visto o desvio-padrão definir a
% frequência de corte. Para os valores do desvio padrão testados (definidos
% a partir do parâmetro alpha da função gausswin), tem-se que para o
% desvio padrão de 2.5 obteve-se um bom equilibrio entre remoção do ruido e
% a conservação do detalhe da imagem. Considerando desvios-padrão
% superiores, obtem-se uma gaussiana mais ampla no dominio discreto, o que
% tem como correspondente no dominio espectral uma gaussiana mais estreita,
% tendo associada uma frequência de corte mais reduzida e, como tal, uma
% filtragem mais extensa de frequências, removendo-se o ruido mas
% sacrificando-se o detalhe da imagem de forma mais pronunciada. Por outro
% lado, se o desvio-padrão é demasiadamente reduzido, obtem-se uma
% gaussiana estreita no dominio discreto, sendo o sua equivalente no
% dominio espectral uma gaussiana mais ampla, tendo associado uma
% frequência de corte mais elevada e, como tal, uma menor eficiência na
% remoção do ruido, contudo melhor preservando o detalhe da imagem. 

% Como tal, existe um trade-off entre a remoção do ruido (melhor para σ 
% maiores, janelas mais amplas no dominio temporal), e o detalhe da imagem 
% (melhor para σ menores, janelas menos amplas no dominio temporal). Como
% tal, conclui-se que a filtragem gaussiana 2D nao é ideal para a remoção
% de SPN, visto estar associada a perdas do sinal original. 

    %% 3)
% Imagem corrompida para uma densidade de ruido de 0.05, afetando
% aproximadamente 5% dos pixeis
I_corrupted = imnoise(I_original,'salt & pepper');

size_janela = [2,5,10,15]; % vetor contendo as dimensões da janela a serem testadas 

for i = 1:length(size_janela)
    % filtragem por convolução no dominio discreto, aplicando-se o filtro
    % de mediana 2D com uma janela quadrada 
    I_medfilt = medfilt2(I_corrupted,[size_janela(i) size_janela(i)]);
    
    figure(23+i); 
    colormap gray; imagesc(I_medfilt);
    title(['Median Filtering with square window size [', ...
        num2str(size_janela(i)), ', ', num2str(size_janela(i)), ']']);
end

%filtragem por convolução no dominio discreto, aplicando-se o filtro
%de mediana 2D com uma janela retangular
I_medfilt = medfilt2(I_corrupted,[size_janela(1) size_janela(3)]);

figure(28);
subplot(1,2,1);
colormap gray; imagesc(I_medfilt);
title(['Median Filtering with [', num2str(size_janela(1)), ', ', ...
    num2str(size_janela(3)), '] rectangular window size']);

%filtragem por convolução no dominio discreto, aplicando-se o filtro
%de mediana 2D com uma janela retangular
I_medfilt = medfilt2(I_corrupted,[size_janela(3) size_janela(1)]);

subplot(1,2,2);
colormap gray; imagesc(I_medfilt);
title(['Median Filtering with [', num2str(size_janela(3)), ', ', ...
    num2str(size_janela(1)), '] rectangular window size']);

% Como já anteriormente mencionado, o filtro de mediana corresponde a um
% tipo de filtragem não-linear que computa a intensidade de um pixel como 
% o valor mediano de intensidade de um conjunto de pixels vizinhos dentro 
% de uma janela pré-definida centrada no pixel em questão, sendo que para
% os pixeis na borda da imagem o cálculo é realizado por zero-padding.

% Analisando o impacto dos parâmetros da janela do filtro de mediana, tem-se 
% que o aumento da sua dimensão até a um determinado valor critico [5,5] conduz 
% à progressiva remoção dos valores de intensidade extremos introduzidos pela
% saturação produzida pelo SPN, perservando-se o detalhe da imagem, ao
% contrário do verificado para a filtragem passa-baixo gaussiana, em que
% uma eficaz remoção de ruido está associada à desfocagem da imagem (perda
% de detalhe). Contudo, para janelas de dimensão superior ao valor
% critico, verifica-se que o ruido é removido mas ocorre desfoque (perda de
% detalhe) da imagem, atuando de forma similar a um filtro passa-baixo. Tal
% deve-se ao facto de para janelas de dimensão consideravel, aumentar-se 
% a probabilidade de ser considerado para um dado pixel um valor de intensidade 
% relativa a um pixel distante, de reduzida correlação face ao pixel 
% correspondente da imagem original, conduzindo ao desfoque da imagem. 

% Por outro lado, destaca-se que a performance do filtro é também afetada
% pela densidade de ruido que corrompe a imagem, sendo que para densidades
% consideravelmente elevadas a supressão do SPN é insuficiente para
% janelas de pequena dimensão, em que a mediana tenderá igualmente para valores
% extremos (corrompidos). Considerando janelas de maior dimensão, a remoção
% do ruido é atingida mas ocorre desfoque (perda de detalhe) da imagem, 
% atuando de forma similar a um filtro passa-baixo. Como tal, conclui-se
% que para imagens corrompidas por SPN de elevada densidade, não se
% consegue atingir um trade-off razoável entre remoção de ruido e detalhe
% da imagem (preservação dos edges). 

% De forma geral, conclui-se que o filtro de mediana é um método mais
% adequado, face à filtragem gaussiana, para remover ruido impulsional 
% como o SPN 

    %% 4) 

% Este operador utiliza duas matrizes 3 por 3 que sao convolvidas com a 
% imagem original de modo a calcular valores aproximados das derivadas 
% para variaçoes horizontais e verticais. Obtêm-se três estimativas do gradiente, 
% realizando-se depois a média, reduzindo-se o ruido associado à estimativa
% do gradiente 

%mascaras 2D
Gx = [-1 0 1; -2 0 2; -1 0 1]; 
Gy = [1  2  1; 0  0  0; -1 -2 -1];

%convolve as mascaras com  a imagem 
I_gx = conv2(I_original,Gx,'same'); % primeira derivada na direção horizontal
I_gy = conv2(I_original,Gy,'same'); % primeira derivada na direção vertical
I_g_magnitude = sqrt(I_gx.^2 + I_gy.^2); % magnitude do gradiente para cada pixel 

figure(29)
subplot(2,2,1); colormap gray; imagesc(I_original);
title('Original image')
subplot(2,2,2); colormap gray; imagesc(I_gx);
title('Horizontal Gradient - Vertical Component of Boundaries')
subplot(2,2,3); colormap gray; imagesc(I_gy);
title('Vertical Gradient - Horizontal Component of Boundaries')
subplot(2,2,4); colormap gray; imagesc(I_g_magnitude);
title('Gradient Magnitude - Bondaries of the Filtered Image')

    %% 5) 

% As fronteiras numa imagem podem ser detetadas através da análise de 
% descontinuidades locais da intensidade da imagem que, à partida, 
% serão máximas ou mínimas nas fronteiras. Para a sua deteção, utilizam-se 
% filtros 2D representados por uma matriz de convolução (kernel) que tanto 
% pode estimar, em cada pixel, as derivadas discretas de primeira ordem
% (detetar os máximos locais ou minimos da primeira derivada) como de 
% segunda ordem (detetar zero-crossings da segunda derivada). De acordo com 
% a sua dimensão, simetria e coeficientes, estes serão mais adequados 
% consoante o tipo de fronteiras (orientação, ruído e estrutura – graduais 
% ou abruptas).

% Considerando as máscaras aplicadas no exercicio anterior, tem-se que 
% Gx representa as derivadas parciais na direção horizontal x ([-1 0 1],
% derivada com respeito ao pixel central), e Gy as derivadas na direção
% vertical y ([1 0 -1]'). O kernel do operador Sobel aplica em conjunto um
% filtro passa-baixo na direção perpendicular à da derivada, ao realizar uma
% média ponderada das estimativas da derivada, dando mais peso à
% estimativa da derivada na linha central. 

% As máscaras Gx e Gy são concebidas de forma a detetar fronteiras numa
% direção unica, com a direção do gradiente a ser perpendicular às
% fronteiras a detetar: fronteiras verticais são apenas detetadas pela
% convolução da imagem com a mascara Gx, visto detetar variações na 
% horizontal, sendo que a mascara Gy nesta situação devolve um valor de 0; 
% fronteiras horizontais são detetadas pela convolução da imagem com a
% máscara Gy, visto detetar variações na vertical. Por outro lado, ainda se
% distinguem dois valores diferentes das estimativas das derivadas
% parciais: se a variação da intensidade diminuir da esquerda para a
% direita, resulta numa estimativa negativa da derivada horizontal x; se a
% variação da intensidade aumentar da esquerda para a direita, resulta uma
% estimativa positiva da derivada horizontal x. Por outro lado, se não for
% detetada nenhuma variação numa dada direção, a estimativa da derivada
% parcial é 0. Um raciocinio semelhante pode ser aplicado às estimativas 
% das derivadas verticais y. 

% Desta forma, quando as máscaras são aplicadas separadamente à imagem, 
% as mesmas representam a componente do gradiente em cada direção, detetando 
% fronteiras horizontais e verticais. Contudo, fronteiras inclinidas podem
% tambem ser detetadas ao combinar-se I_gx e I_gy de modo a computar a 
% magnitude absoluta e direção do gradiente em cada pixel. A direção do 
% gradiente, isto é, a direção da máxima variação, é dada por: 
% theta = atan(I_gy/I_gx) (em graus) = atan(I_gy/I_gx)*180/pi (em radianos) 

% Quanto à magnitude absoluta do gradiente, tem-se que o mesmo é calculado 
% através da norma do vetor que aponta no sentido da máxima variação, sendo 
% a norma máxima quando as variações (componentes em ambas as direções x e 
% y) são máximas. Se as variações forem nulas, a imagem filtrada para esse 
% pixel apresenta intensidade nula. 

% Visto as imagens serem representadas na escala de cinzentos (imagens de 
% intensidade), o matlab atribui a cor preta ao valor de intensidade mais 
% baixo, sendo para o caso da magnitude absoluta o valor 0, ao passo que para 
% as imagens que representam os gradientes nas direções x e y ao valor 0 
% era atribuida uma cor cinzenta. Deste modo, os pixeis a branco correspondem
% às fronteiras, que dizem respeito às regiões de maior variação da 
% intensidade. 

% Por fim, destaca-se que filtros que aproximam derivadas atuam como
% filtros passa-alto, os quais tendem a amplificar o ruido, existindo o
% trade-off de que filtros maiores conduzem à redução da amplificação do
% ruido, contudo pioram a deteção da correta localização das fronteiras. 

        %% 6) 
    
%mascaras 2D
L4 = [0 -1  0; -1  4 -1; 0 -1  0];
L8 = [-1 -1 -1;-1  8 -1; -1 -1 -1];

%convolve as mascaras com  a imagem, considerando duas representações do 
%operador Laplaciano
I_L4 = conv2 (I_original,L4,'same');
I_L8 = conv2 (I_original,L8,'same');

figure(31)
subplot(1,3,1); colormap gray; imagesc(I_original);
title('Original image');
subplot(1,3,2); colormap gray; imagesc(I_L4);
title('Image Convolved with mask L4');
subplot(1,3,3); colormap gray; imagesc(I_L8);
title('Image Convolved with mask L8');

% O Laplaciano corresponde a uma medida isotrópica 2D (isto é, não depende 
% da direção) da segunda derivada espacial de uma imagem, calculando a divergência do 
% gradiente. Sendo a imagem de entrada representada como uma matriz 2D de 
% pixeis discretos, podem ser concebidos kernels de convolução discretos que
% aproximem as segundas derivadas segundo a definição do Laplaciano, sendo 
% L4 e L8 dois tipos de kernels aplicados para este fim, correspodendo a 
% operadores laplacianos negativos, isto é, detetam fronteiras internas. De destacar que 
% a máscara L8 possibilita uma deteção mais eficaz por ser mais sensivel a detetar 
% fronteiras diagonais ao pixel central, uma vez o kernel considerar os 
% pixeis na diagonal principal. 

% Como anteriromente mencionado, fronteiras podem ser identificadas ao
% detetarem-se zero-crossings da segunda derivada da intensidade,
% correspondendo às variações máximas da primeira derivada (computada por
% exemplo pelo filtro Sobel aplicado em 4), realçando 
% regiões de rápida variação de intensidade. Assim, tal como para o filtro
% sobel, convolvendo a imagem com estas mascaras é possivel detetar as
% fronteiras da imagem, isto é, as variações bruscas nos valores de
% intensidade dos pixeis, correspondendo a um filtro passa-alto. 

% De notar que, contudo, o laplaciano é frequentemente aplicado a uma
% imagem após filtragem passa-baixo (filtro gaussiano 2D, por exemplo), 
% por forma a reduzir a sensibilidade ao ruído, ignorando os zero-crossings
% produzidos por pequenas mudanças na intensidade da imagem. Dado o
% discutido previamente, decidiu-se aplicar uma janela gaussiana de N = 64
% e σ = 2.5

h_g = fspecial('gaussian',64,2.5);

% Visto a operação de convolução ser associativa, pode-se primeiramente 
% convolver o filtro passa-baixo gaussiano com o filtro Laplaciano e,
% depois, convolver este filtro passa-banda hibrido com a imagem: 
% I_filtrada = (I*G)*L = I*(L*G)

filtro_bandpass = conv2(h_g, L8, 'same');
I_filtrada_LoG = conv2(I_original,filtro_bandpass,'same');

figure(32);
subplot(1,3,1); colormap gray; imagesc(I_original);
title('Original image');
subplot(1,3,2); colormap gray; imagesc(I_L8);
title('Image Convolved with mask L8');
subplot(1,3,3);colormap gray; imagesc(I_filtrada_LoG);
title('Image Convolved with mask LoG');

% De notar que a localização de zero-crossing na imagem filtrada gaussiana 
% tende para a sua correta localização quando σ → 0, ou seja, quanto mais
% estreito o filtro gaussiano for mais o kernel LoG aproxima-se do kernel
% Laplaciano. 

    %% 7) 

% obtenção da imagem binária das fronteiras, obtendo-se os thresholds 
% determinados por default para ambos os métodos de deteção de fronteiras
[I_sobel,thresh_default_sobel] = edge(I_original, 'sobel');
[I_canny,thresh_default_canny] = edge(I_original, 'canny');

figure(33);
subplot(1,4,1); colormap gray; imagesc(I_g_magnitude);
title('Sobel Edge');

subplot(1,4,2); colormap gray; imagesc(I_sobel);
title(['Sobel Edge (MATLAB), threshold = ',...
    num2str(round(thresh_default_sobel,2))]);

subplot(1,4,3); colormap gray; imagesc(I_filtrada_LoG);
title('L8oG edge - Convolved with L8 and Gaussian');

subplot(1,4,4); colormap gray; imagesc(I_canny);
title(['Canny Edge, threshold = [', ...
    num2str(round(thresh_default_canny(1),2)), ', ', ....
    num2str(round(thresh_default_canny(2),2)), ']']);

% Filtros passa-alto para deteção de fronteira implementados pela função
% built-in do matlab edge, a qual devolve como output uma imagem binária
% das fronteiras detetadas. Esta imagem binária pode ser convertida para
% uma máscara complementar (recorrendo-se ao comando imcomplement), a qual, 
% por multiplicação ponto a ponto com a matriz de de intensidade da imagem, 
% atribui intensidade nula aos pixeis das fronteiras detetadas, produzindo 
% uma imagem em que as fronteiras sao removidas. 
figure(34);
I_rem_sobel2 = I_original .* imcomplement(I_sobel);
subplot(1,2,1); colormap gray; imagesc(I_rem_sobel2);
title('Sobel (MATLAB) Edge Removed Image');

I_rem_canny = I_original .* imcomplement(I_canny);
subplot(1,2,2); colormap gray; imagesc(I_rem_canny);
title('Canny Edge Removed Image');


[I_LoG,thresh_default_LoG] = edge(I_original, 'LoG');

figure(35);
subplot(1,2,1);colormap gray; imagesc(I_LoG);
title(['LoG Edge, threshold = ',...
    num2str(thresh_default_LoG)]);

I_rem_LoG = I_original .* imcomplement(I_LoG);
subplot(1,2,2); colormap gray; imagesc(I_rem_LoG);
title('LoG Edge Removed Image');

% Ao contrário do procedimento previamente aplicado de filtragem passa-alto
% para detetar as fronteiras da imagem, tem-se que a função built-in edge
% do matlab procede a dois passos adicionais automaticos: filtragem passa-baixo 
% realizada a priori, considerando a posteriori da filtragem passa-alto 
% um processo de thresholding, sendo os thresholds relativos à magnitude 
% do gradiente (filtros que computam derivadas de 1ªordem) ou à magnitude 
% dos zero-crossings (filtros que computam derivadas de 2ªordem), 
% considerando apenas pixeis como pertencentes a fronteiras aqueles que
% apresentem uma maginutude da 1ª ou 2ªderivada superior aos thresholds 
% definidos. 

% Deste modo, a diferente implementação dos algoritmos explica 
% as diferenças obtidas: para o filtro de sobel, o comando edge realiza o
% threshold da imagem filtrada, sendo eliminadas alguns edges, visto o
% criterio ser mais exigente. 

% Já relativamente ao filtro Canny, verifica-se o melhor desempenho na 
% deteção das fronteiras, devendo-se à maior robustez do algoritmo, o qual
% primeiramente convolve um filtro gaussiano com a imagem. Após isto,
% procura maximos locais do gradiente (calculado através da
% derivada de um filtro gaussiano), realizando uma operação de thresholding,
% classificando weak edges como edges que tenham uma magnitude do 
% gradiente superior ao low threshold e encontrem-se conetados a edges 
% com uma magnitude do gradiente superior ao high threshold (strong edges),
% permitindo suprimir maximos locais que correspondecem a ruido. O algoritmo 
% descrito permite uma maior conservação da integridade das fronteiras
% detetadas. 

% Contrariamente, o filtro Sobel deteta menos fronteiras e, quando detetadas, 
% com mais descontinuidades. Tal deve-se não só à preferência na deteção 
% de fronteiras horizontais e verticais, com consequente inferior sensibilidade 
% a outras orientações que representem uma continuidade das anteriormente 
% detetadas, como também à eliminação por threshold de possíveis weak edges,
% apresentando por isso maior sensibilidade a transições de intensidade 
% mais abruptas.
% Considerando thresholds superiores, prevê-se que os os filtros diminuam 
% a sua sensibilidade a detetar fronteiras. Nessas situações, o filtro LoG 
% (para threshold adequado) pode constituir uma alternativa para detetar 
% fronteiras entre regiões de menor contraste (reduzida magnitude do 
% gradiente), por efetivamente corresponderem a zero-crossings. 

    %% 8) 

% As equações às diferenças que representam os operadores Laplacianos L4 e L8 
% encontram-se derivadas em anexo. 

% Esquematicamente, as equações às diferenças obtêm-se sobrepondo a máscara (L4 ou L8) à
% imagem original, realizando uma soma ponderada das entradas consoante os coeficientes do kernel. 
% Para kernels de 3x3, ponderam-se unicamente os vizinhos que distam 1
% pixel do pixel central, obtendo-se uma imagem y(n,m) final: 

% L4
% y(n,m)=4*x(n,m)-x(n-1,m)-x(n,m-1)-x(n,m+1)-x(n+1,m)

% L8
% y(n,m)=8*x(n,m)-x(n-1,m-1)-x(n-1,m)-x(n-1,m+1)-x(n,m-1)-x(n,m+1)-x(n+1,m-1)-x(n+1,m)-x(n+1,m+1)

figure(36);
subplot(1,2,1); colormap gray; imagesc(conv2(I_original,-L4,'same'));
title('Convolved with -L4');
subplot(1,2,2); colormap gray; imagesc(conv2(I_original,L4,'same'));
title('Convolved with L4');

% Note-se que L4 representa o operador Laplaciano negativo, pelo que deteta
% fronteiras internas (direção do gradiente da fronteira para o centro da
% imagem), podendo-se constatar na figura 36 a imagem original convolvida
% com o laplaciano positivo (deteta fronteiras externas) e o laplaciano negativo. 