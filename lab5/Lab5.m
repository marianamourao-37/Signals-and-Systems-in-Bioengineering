%Grupo 22: 
%-Ana Rita Lopes nº98587
%-Mariana Mourão nº98473
 
%LAB#5 - Adaptive Filtering

clear all;
close all;
clc;

%% PART I - SYNTHETIC DATA

%% 1)

% Implementaçao da funcao pedida, adaptFilter.m, assim considerando 
% tambem o critério de otimização relativo à minimização da norma 
% do vetor de erros obtidos para uma dada dimensão de janela 

%% 2) Gerar o sinal

%Gera um sinal designado s(n), sendo a funcao degrau deslocada 0.5 segundos.

fs = 4000; %frequencia de amostragem 
T = 1;  %duração 
t = linspace(0,T,fs)';  %vetor coluna contendo o de tempo discretizado 

t_descontinuidade = 0.5; % tempo (em segundos) em que ocorre a descontinuidade da função step 
n_05 = t_descontinuidade*fs; % amostra em que ocorre a descontinuidade da função step

s = zeros(T*fs,1); % inicialização da função step 
s(n_05:end,1)=1; % atualização do degrau 

% Testar a implementação do filtro adaptativo, considerando o seguinte:
% d(n)=s(n); x(n) = 1

d = s;
x = ones(length(s),1);

% testar para ambos os critérios de otimização, isto é, minimização 
% do erro em cada instante n (1º critério) e minimização da norma 
% do vetor de erros obtidos para uma dada dimensão de janela

%% - variação da ordem do filtro para os dois critérios:
p_vetor = [2,5,10,25,50]; % vetor das ordens do filtro a serem testadas 
w = 200; % fixar o valor da janela 

figure;
for i = 1:length(p_vetor)
    [e1,y1,~] = adaptFilter(d,x,p_vetor(i));

    [e2,y2,~] = adaptFilter(d,x,p_vetor(i),w);

    subplot(5,2,2*i-1)
    
    plot(t,e1,t,y1);
    title(['Filtro Adaptativo (1ºcriterio), p = ', num2str(p_vetor(i))]);
    xlabel('Tempo (s)');
    legend ('e(n)','y(n)');
    axis tight
    
    subplot(5,2,2*i)
    plot(t,e2,t,y2);
    title(['Filtro Adaptativo (2ºcriterio), p = ',num2str(p_vetor(i))]);
    xlabel('Tempo (s)');
    legend ('e(n)','y(n)');
    axis tight
end
suptitle('Teste para diferentes ordens do filtro');

%% - variação da dimensão da janela:
w_vector = [50,200,500,1000]; % vetor das dimensões da janela a serem testadas
p = 25; % fixar a ordem do filtro

figure
for i = 1:length(w_vector)
  
    [e,y,~] = adaptFilter(d,x,p,w_vector(i));
    
    subplot(2,2,i)
    plot(t,e,t,y);
    title(['Filtro adaptativo, w = ',num2str(w_vector(i))]);
    xlabel('Tempo (s)');
    legend ('e(n)','y(n)');
    axis tight
end
suptitle('Teste para diferentes dimensões da janela')

%% resposta impulsiva hn(tau)
p = 50; % ordem do filtro  
[~,~,m_coeffs] = adaptFilter(d,x,50);

figure;
j = 400:400:4000; % estimativas de hn(tau) a aceder
for i = 1:length(j)
    subplot(5,2,i);
    plot(0:p,m_coeffs(:,j(i)));

    xlabel('Coeficientes');
    ylabel('Pesos');
    title(sprintf('Resposta impulsiva h_{%i}', j(i)));
    axis tight
end

% 0.019607843137255*51 = 1
% x(3000:-1:3000-50)'*ones(p+1,1)*0.019607843137255;
    
% Conclusões:
% Considerando o primeiro critério de otimização (sem janela), 
% o filtro reproduz de forma exata o sinal original d(n). 
% Relativamente ao 2º critério de otimização (com janela), 
% a reconstrução do sinal original d(n) é afetada devido ao 
% facto dos coeficientes do filtro ajustarem-se de forma mais
% gradual, introduzindo-se um maior delay, neste caso durante a 
% adaptação à descontinuidade do sinal d(n), resultando num erro 
% maximo para t = 0.5 segundos, o qual decresce à medida que 
% o ajuste é conseguido. Quanto maior o tamanho da janela, maior 
%será o erro na região de transição, t= 0.5 segundos. 

% Devido ao facto do filtro adaptativo dizer respeito a um sistema 
% linear variante no tempo, para cada amostra n computa-se uma 
% nova resposta impulsiva hn(tau) do filtro, com os seus coeficientes 
% ajustados de forma a minimizar o erro d(n) - y(n), de acordo com a
% referência x(n). Se d(n) = s(n) + miu(n), com x(n) e s(n) 
% não correlacionados, mas miu(n) e x(n) correlacionados (mas não iguais), 
% então o filtro irá tentar aproximar d(n) e y(n) através 
% da correlação entre x(n) e miu(n), atingindo-se 
% performance máxima de cancelamento ativo de ruido, com o erro 
% minimo a corresponder a s(n). Contudo, se existir uma correlação entre
% s(n) e miu(n), o processo de cancelamento ativo de ruido pode 
% resultar na perda de informação. 

% Para o caso especifico deste problema, tem-se que
% para t<0.5 segundos, x(n) = 1 e d(n) = 0. Pra t >= 0.5 
% segundos, x(n) = d(n) = 1. Para a segunda situação, 
% os coeficientes devem convergir para o delta dirac, 
% ainda que uma equipartição dos pesos para todos 
% os p+1 coeficientes também atinge y(n) = 1, e, como tal, a minimização da
% a energia do erro e(n), visto y(n) = sum(i=0,i=p) hi*x_(n-i), 
% com x(n) = 1 para todos os n, pelo que h(n) = 1/(p+1) para todos os n 
% tambem resulta na minimização da energia do erro. 


%% 3) Gerar o ruido

%% 3a) 
%gerar um sinal de ruido branco gaussiano (gaussian white noise signal), n(n), 
%com as mesmas dimensoes do sinal s(n) e filtrar esse sinal com um filtro
%Butterworth passa baixo de 10ºordem 

length_n=length(s);
miu = wgn(length_n,1,-30); % ruido gaussiano branco aditivo
%gera um vetor coluna de comprimento igual ao sinal de s(n), 
% com amostras de ruido gaussiano branco aditivo (em volts), 
% especificando-se a potência de ruído em -30 dBW.

[b,a] = butter(10,1/3); 
% Implementação de um filtro Butterworth digital passa baixo de 10ª ordem, 
%com frequência de corte w = pi/3. Para tal, recorreu-se à 
%função built-in do matlab butter, determinando-se os coeficientes
% a e b do filtro. w é dada em radianos/amostra (rad/sample) e 
%é normalizada entre 0 a pi. Como tal, o filtro low-pass 
%requerido tem uma frequência de corte normalizada wn=1/3.

%epsilon(n)=g(n)*miu(n), sendo uma versao filtrada do ruido branco 
% obtido por convolução linear com  g(n), resposta impulsiva do filtro IIR
% Butterworth. Disto resulta um sinal designado ruido colorido 
epsilon=filter(b,a,miu);

%d(n) é a soma da função degrau, s(n), com o ruído colorido epsilon(n)
d = s + epsilon;

%Comparar o sinal miu(n) e epsilon(n)
figure(3);
plot(t,miu,t,epsilon);
title('Sinais Ruidosos');
xlabel('Time (s)');
legend ('Ruido Branco Gaussiano','Ruido Colorido');    
axis tight

% ruido colorido resultante da convolução temporal entre o filtro 
% IIR passa-baixo e ruído branco gaussiano, para uma frequências de corte 
% normalizada wn = 1/3 relativamente à frequencia de amostragem fs, 
% pelo que fc = 4000/3 = 1333,(3) Hz. O ruido colorido de largura de banda 
% inferior varia mais lentamente em amplitude, por conter frequências menos
% elevadas, em consequência da aplicação do filtro passa baixo, 
% sendo que quanto menor fc, mais notorio será o efeito da média 
% ponderada sobre o ruido branco, passando as amostras a assumir correlação
% com outras adjacentes. Isto pode ser mais aprofundado ao estudarem-se 
% as funções de autocorrelação de miu(n) e epsilon(n). 

%% EXTRA

% FWHM do pico máximo da função de autocorrelação enquanto critério 
% de medição da sua largura, estando o pico definido em cor=1, pelo 
% que a meia altura é cor=0.5. Contudo, devido à discretização da 
% função de autocorrelação, o mesmo poderá não conter pontos em 
% cor=0.5, pelo que se recorreu à função find(cor>0.5, 1, ‘last’) 
% que retorna o índice do último elemento do vetor cor que satisfaz 
% cor>0.5. Com base no ponto obtido com a função find e o 
% imediatamente depois, procedeu-se à interpolação 
% por regressão linear do lags correspondente a cor=0.5, através da 
% função regression que devolve o declive m e ordenada na origem b 
% da reta. Por fim, recordando a simetria característica da função 
% de autocorrelação, calculou-se a FWHM como o dobro do lag 
% obtido por interpolação,

N = length(miu);
text_vector = ["Ruido Branco Gaussiano miu","Ruido Colorido epsilon"];

[cor, lags] = xcorr(miu,'coeff'); % cor normalizado pelo 
    % parâmetro ‘coeff’ (-1≤cor≤1, com cor=1 e cor=-1 para sinais 
    %em fase e anti-fase, respetivamente, permitindo análise 
    % independente da escala das variáveis).

maximo = max(cor); % É igual a 1 devido à normalização realizada pelo 
        %parâmetro 'coeff', correspondendo ao momento em que os sinais estão em 
        %fase.

abovehalfmax = cor>maximo/2; % Find where it's more than half the max. 
%(=0.5 em consequência da linha anterior)

% Pontos imediatamente antes e depois de cor=0.5: 
lastindex = find(abovehalfmax, 1, 'last'); % Índice do último valor 
%cor>0.5, ou seja, imediatamente antes de cor=0.5
index2 = lastindex+1; % Índice imediatamente depois do último índice 
%correspondente a cor>0.5

% Obtenção da reta de regressão linear tendo em conta os pontos
%[cor(lastindex),lags(lastindex)] e [cor(index2), lags(index2)]:
[~,m,b_r]=regression([lags(lastindex) lags(index2)],[cor(lastindex) cor(index2)]);
% Obtém-se o declive m e a ordenada na origem b, para o posterior 
% cálculo dos índices interpolados que correspondem a cor=0.5. 

% Interpolação por regressão linear tendo em conta que a equação da
%reta obtida é cor=m*lags + b. Pretende-se saber o valor 
%interpolado por regressão linear do índice para o qual cor=0,5:  
lags_interpolacao = (0.5 - b_r)/m;

largura_rxx(i,1)= 2*lags_interpolacao; % Devido à simetria da função 
%de autocorrelação
figure;
plot(lags,cor);
title("Ruido Branco Gaussiano miu");
grid on;
axis([-N/2 N/2 -.5 1.1]);  
xlabel('Lags (amostras)');
ylabel('Rxx');

hold on;
% Linhas que evidenciam graficamente a largura da função de
%autocorrelação a partir do ponto interpolado (linha azul) e a
%partir dos pontos usados na interpolação (linhas tracejadas):
line([-lags_interpolacao lags_interpolacao -lags_interpolacao;  ...
   -lags_interpolacao lags_interpolacao lags_interpolacao], ...
   [-maximo/2 -maximo/2 maximo/2; maximo/2 maximo/2 maximo/2], ...
   'Color', 'b', 'LineWidth', 1,'LineStyle', '--');

%texto a adicionar aos gráficos da função de autocorrelação,
%explicitando-se o valor da FWHM com e sem interpolação:
text(-48,0.9,['FWHM=', num2str(2*lags_interpolacao)],'Color',...
   'b','FontSize',7); %FWHM interpolada

wn = (2:10).^(-1); % Vetor wn que contém frequências de corte 
% normalizadas para [0,pi].

for k = 1:length(wn) 
    
    [b_w,a_w] = butter(10,wn(k)); 
    epsilon=filter(b_w,a_w,miu);
        
    [cor, lags] = xcorr(epsilon,'coeff'); % cor normalizado pelo 
    % parâmetro ‘coeff’ (-1≤cor≤1, com cor=1 e cor=-1 para sinais 
    %em fase e anti-fase, respetivamente, permitindo análise 
    % independente da escala das variáveis).
    
    maximo = max(cor); % É igual a 1 devido à normalização realizada pelo 
    %parâmetro 'coeff', correspondendo ao momento em que os sinais estão em 
    %fase.

    abovehalfmax = cor>maximo/2; % Find where it's more than half the max. 
    %(=0.5 em consequência da linha anterior)

    % Pontos imediatamente antes e depois de cor=0.5: 
    lastindex = find(abovehalfmax, 1, 'last'); % Índice do último valor 
    %cor>0.5, ou seja, imediatamente antes de cor=0.5
    index2 = lastindex+1; % Índice imediatamente depois do último índice 
    %correspondente a cor>0.5

    % Obtenção da reta de regressão linear tendo em conta os pontos
    %[cor(lastindex),lags(lastindex)] e [cor(index2), lags(index2)]:
    [~,m,b_r2]=regression([lags(lastindex) lags(index2)],[cor(lastindex) cor(index2)]);
    % Obtém-se o declive m e a ordenada na origem b, para o posterior 
    % cálculo dos índices interpolados que correspondem a cor=0.5. 

    % Interpolação por regressão linear tendo em conta que a equação da
    %reta obtida é cor=m*lags + b. Pretende-se saber o valor 
    %interpolado por regressão linear do índice para o qual cor=0,5:  
    lags_interpolacao = (0.5 - b_r2)/m;

    largura_rxx(i,1)= 2*lags_interpolacao; % Devido à simetria da função 
    %de autocorrelação
    
    figure;
    plot(lags,cor);
    title(['Ruido Colorido epsilon; Largura de banda = ', num2str(wn(k))]); % Describe bandwidth
    grid on;
    axis([-N/2 N/2 -.5 1.1]);  
    xlabel('Lags (amostras)');
    ylabel('Rxx');

    hold on;
   % Linhas que evidenciam graficamente a largura da função de
   %autocorrelação a partir do ponto interpolado (linha azul) e a
   %partir dos pontos usados na interpolação (linhas tracejadas):
   line([-lags_interpolacao lags_interpolacao -lags_interpolacao;  ...
       -lags_interpolacao lags_interpolacao lags_interpolacao], ...
       [-maximo/2 -maximo/2 maximo/2; maximo/2 maximo/2 maximo/2], ...
       'Color', 'b', 'LineWidth', 1,'LineStyle', '--');

   %texto a adicionar aos gráficos da função de autocorrelação,
   %explicitando-se o valor da FWHM com e sem interpolação:
   text(-48,0.9,['FWHM=', num2str(2*lags_interpolacao)],'Color',...
       'b','FontSize',7); %FWHM interpolada
   
end
   
% após filtragem passa baixo, as amostras passam a assumir correlação com outras adjacentes, 
% sendo que a função da autocorrelação deixa de ser o delta dirac 
% em lag=0 (sinal sem memória), passando a assumir uma largura 
% dependente da largura de banda do sinal, sendo o grau de 
% predictabilidade dos valores futuros dum sinal com base nos 
% passados (memória do sinal) inversamente proporcional à sua 
% largura de banda. Contudo, para lags sucessivamente maiores a, 
% função de autocorrelação tende para valores nulos

%% 3b)

%Comparar o sinal s(n) e d(n)
figure(4);
plot(t,s,t,d);
title('Adição de Ruido Colorido à Função Step');
xlabel('Time (s)');
legend ('Original','Corrompida');    
axis tight

%% 4) Testar o filtro desenvolvido

%% 4a)

x = miu; % sinal de referencia 
p = 50; % ordem do filtro fixa 

figure;
w_vector = [0 100, 750, 1500]; % vetor das dimensões da janela a serem testadas
for i = 1:length(w_vector)
    [e,~,m_coeff] = adaptFilter(d,x,p,w_vector(i));
    subplot(2,2,i)
    plot(t,s,t,e);
    title(['Filtro adaptativo, w = ',num2str(w_vector(i))]);
    xlabel('Tempo (s)');
    legend ('s(n)','e(n)');
    axis tight
end

suptitle('Filtragem Adaptativa do Ruido Colorido');

figure;
j = 400:400:4000;  % estimativas de hn(tau) a aceder
for i = 1:length(j)
    subplot(5,2,i);
    plot(0:p,m_coeff(:,j(i)));
    xlabel('Coeficientes');
    ylabel('Pesos');
    title(sprintf('Resposta impulsiva h_{%i}', j(i)));    
end

% Aplicando-se o primeiro critério de otimização, o qual minimiza o
% erro em cada instante n, constata-se que o ruido branco 
% ainda é possivel de ser correlacionado com o sinal original d(n). 
% Quanto à filtragem adaptativa que considera o segundo critério 
% de otimização (minimização da norma do vetor de erros obtidos 
% para uma dada dimensão de janela), constata-se que para t>0.5 
% não existe uma boa adapatação do sistema, traudzindo-se em erros 
% superiores. Tal deve-se ao facto dos coeficientes do filtro tornarem.-se
% instanveis, nao se conseguindo adaptar após a descontinuidade em 
% t=0.5 segundos. 
% Ainda assim, para janelas de tamanho superior, o filtro considera 
% mais informação para correlacionar o sinal de referência x(n) com o sinal
% desejado d(n), obtendo-se melhores performances do cancelamento 
% ativo de ruido do que aquelas atingidas para janela sinferiores. Contudo,
% apresenta o trade-off de ser computacioanlmente mais pesado, por computar
% em cada n a pseudo-inverse, isto é, inversão de uma matrix 
% autocorrelação de dimensao consideravel. 

%% 4b)

imp = zeros(4000,1); imp(1) = 1; % delta dirac 
ir_g = filter(b,a,imp); % resposta impulsiva do filtro IIR Butterworth passa-baixo
t_k = round((0:0.1:T)*fs); % vetor com as amostras correspodentes aos 
% t tempos uniformemente distribuídos separados por 0,1 segundos
[~,~,h_coeffs] = adaptFilter(d,x,p,1000);

figure;
for i = 2:length(t_k) 
    subplot(2,5,i-1)
    plot(0:p,ir_g(1:p+1),'-r');
    hold on;
    plot(0:p,h_coeffs(:,t_k(i)),'b');
    
    title(sprintf('Resposta impulsiva h_{%i}; t = %i', ...
        t_k(i),t_k(i)/fs));
    
    title(['Coeficientes = ' num2str(t_k(i))]);
    xlabel('Coeficientes');

    axis tight  
end    

hold off;

J = zeros(length(s),1);

figure;
for i=1:length(s)
    J(i,1) = norm(ir_g(1:p+1,1) - h_coeffs(:,i));
end
plot(J)
xlabel('Amostras');
ylabel('Função de Custo J(n)')

%% 4c)

% Para n<2000 (t<0,5s), h(n) apresenta uma boa adaptação ao sistema, obtendo-se uma
% funcao de custo J(n) perto de 0. Isto ocorre enquanto a conversão do ruído 
%miu para o sinal observado é feita por um sistema linear. 
%A partir de n=2000 (t=0,5s), o adaptive filter não consegue representar a não linearidade do sistema, pelo 
% que este diverge do sistema e apresenta oscilações. 

%% 5) Repete o exercicio 4 

%Repete o exercicio 4 para diferentes miu(n) com o algoritmo LMS (Least mean
%square filter), um algoritmo adaptativo, em que mu é o parametro 
% step-size que controla as caracteristicas de convergência do algoritmo 

%vetor com diferentes step sized a serem testados
mu_p = [0.01,0.06,0.16,0.26,0.5,1];
for i = 1:length(mu_p)

    [~,e,w] = LMSfunction(x,d,mu_p(i),p);
    
    figure(5+i)
    for j = 2:length(t_k)
        subplot(2,5,j-1)
        plot(ir_g(1:p+1),'-r');
        hold on;
        plot(w(:,t_k(j)),'b');
        
        title(sprintf('h_{%i}, Step Size mu = %0.2f', t_k(j),mu_p(i)));
        xlabel('Amostras');
        axis tight 
    end
    
    J = zeros(size(s));
    for q=1:length(s)
        J(q) = norm(ir_g(1:p) - w(:,q));
    end
    
    figure(5);
    subplot(3,2,i)
    plot(0:length(x)-1,J)
    title(['J(n) = ||g(\tau)-h_n(\tau)||, Step size mu= ',num2str(mu_p(i))]);
    xlabel('Amostras')
    
    figure(12)
    t = linspace(0,T,fs)'; 
    subplot(3,2,i)
    plot(t,e)
    hold on
    plot(t,s)
    title(['Step size Mu=',num2str(mu_p(i))]);
    xlabel('Tempo, s')
    legend('e(n)','s(n)')
end

% Em cada iteração, não se obtêm a solucoes otimas, mas sim uma solucao que 
%tenta convergir para a designada solucao otimizada. Invés de minimizar a função 
%de custo J(n) de forma imediata e obter a solução, propôs-se um processo iterativo.  

% Os pontos do gradiente dão nos a variação maxima da função, pelo que o inverso dá a 
%variação minima da funcao.
%De modo a convergir para o mínimo tem de se considerar o inverso do gradiente da função
%de custo. Calcula-se o gradiente da função de custo, e considera-se a direção inversa  

% Inves de se tentar minimizar a funcao de custo para cada amostra,
% realiza-se essa minimização de forma iterativa.

%A solucao otimizada no intante n, hn, pode ser obtida pela variaçao da solucao 
%otimizada na iteracao anterior hn-1, atualizando pelo termo negativo, sendo que mu
%caracteriza o parametro que define a dimensao do salto.Para um salto pequeno, 
%alcançamos a solucao de forma muito lenta, pelo que pode não se atingir a mesma em tempo útil. 
% Contrariamente, com um salto grando, atinge-se o valor mínimo de forma mais rápida, 
%mas com o risto de condicionar a coonvergencia para o mínimo global  

% A parte mais interessante é a descontinuidade que ocorre a metade do intervalo, 
%tendo se rápida convergenecia no inicio para a correta resposta g.  
 
% Quanto maior mu, mais rápida a adaptação, mas tambem maior 
% instabolidade na convergência (maiores oscilações)   

% Neste caso, é observável que para valores de mu mais baixos, h demora
% mais tempo a aproximar-se de g (convergência mais lenta), mas quando o
% faz (para t>0,5s) a oscilação é menor quando comparado com valores de
% mu mais elevados. Tal acontece porque a ponderação dada aos coeficientes do
% filtro é maior. Pelo contrário, para valores de mu mais elevados, apesar
% da oscilação final ser maior, a convergência para o sistema real é mais 
% rápida. Conclui-se ainda que este sistema apresenta um comportamento
% semelhante ao do adaptive filter independentemente do tamanho da window.

%% PART II - REAL DATA

%% 1)
% Separar voz e musica

%sinais têm a mesma frequência de amostragem (8000)
[Noise,fs1] = audioread('musica.wav');
[NoiseVoice,fs2] = audioread('musica-voz.wav');

% Dowsampling ratio
% Decimate: reduz a taxa de amostragem do sinal por um fator (3)
% Filtro + Down Sample dos audiofiles y e w 
% Noise=decimate(xh,3);  %apenas musica
% NoiseVoice=decimate(dh,3); %musica mais voz
% fs_decimate = fs2/3; 

%testing for different step sizes
mu = [1, 0.1, 0.01, 0.001];

figure
for i = 1:(length(mu)) 
    [y,e,w] = LMSfunction(Noise, NoiseVoice,mu(i),32); 
    voice = e; %A voz corresponde ao erro que retiramos a partir da filtragem adaptativa
    Nv=length(voice);
    subplot(2,2,i)
    plot(0:Nv-1,NoiseVoice,0:Nv-1,voice);
    title(['Remoção do ruído para mu= ', num2str(mu(i))]);
    legend('Musica','Voz'); 
end
% O LMS utiliza um algoritmo no qual há um ajuste contínuo dos coeficientes
% do filtro designados como pesos de derivacao. Estes tem como o objetivo o de 
%minimizar a diferençaa média quadrática entre os dois sinais de entrada, sendo estes
%atualizados pela expressão: W(n+1)=w(n)+,mu*e(n)*input, onde mu é o step size 
%e determina a estabilidade e a taxa de convergência do algoritmo.

% Ao estudar o efeito da variação do step size (mu) na separação da voz, sendo esta
%correspondente ao erro da filtragem da música, podemos observar que quando
%avaliamos para pequenos step sizes, o erro entre os sinais da música+voz com música é 
%pequeno pelo que a convergência é mais lenta e é possível ouvir música por detrás da voz. 
%Ao aumentar o valor de mu, a velocidade da convergência aumenta,verificando-se uma 
%eliminação da música da voz mais eficiente. Contudo quando mu toma 
% valores elevados, pode ocorrer divergência, podendo fazer com que o
% filtro salte os valores minimos. Neste caso não se consegue obter a voz apenas um sinal muito curto. 

% Para os valores testados, observou-se que o teste que permitiu obter uma voz mais clara,
% ou seja, o teste que permitiu uma melhor extração, foi para o valor de mu=0.1. 
%Tal é observado nos gráficos onde se vê uma diferença maior entre a voz e o sinal música+voz para 
%mu=0.1. Mesmo assim neste teste ainda é possível ouvir um pouco da música de fundo. 

%Os resultados obtidos podem ser reproduzidos usando a função built-in do Matlab (dsp.LMSFilter), 
%obtendo-se para mu=1 o melhor resultado. Ora, as diferenças na performance podem dever-se à propria 
%implemnetação de ambos os algoritmos, podendo o filtro do matlab nao considerar o fator 2 no calculo 
%do jacobiano (-2*e(n)*input), ou seja, no nosso algoritmo equivale a considerar mu = mu/2. 
%Na realidade, para mu = 1/2, obtem-se uma separação dos audios satisfatoria. 

%% Função built-in do matlab 
figure
for i = 1:(length(mu)) 
    %Utilizar o algoritmo do MATLAB LMS filter
    lms1 = dsp.LMSFilter(32,'StepSize',mu(i)); 
    
    %filtragem
    [y1,e1,w1] = lms1(Noise,NoiseVoice);
    voice_1=e1;
    N= length(NoiseVoice);
    subplot(1,length(mu),i)
    plot(0:N-1,NoiseVoice,0:N-1,voice_1);
    title(['Remoção do ruído com LMS para mu= ' num2str(mu(i))]);
    xlabel('Amostras');
    legend('Música','Voz');
    axis tight
end

%Decidimos testar a separação da ruido da musica+voz a partir do algoritmo do Matlab LMS 
%(dsp.LMSFilter), que consideramos anteriormente como sendo mais eficiente.

%% Reproduzir resultados - Função implementada 

%De modo a reproduzir o melhar extracao possivel, utilizamos o mu=1

[y,e,w] = LMSfunction(Noise,NoiseVoice,0.1,32); 
voice = e;
soundsc(voice,fs1);

%% - Função Built-in do Matlab 
mu=1;
lms1 = dsp.LMSFilter(32,'StepSize',mu); 
[y1,e1,w1] = lms1(Noise,NoiseVoice);
voice_1 = e1;
soundsc(voice_1,fs1);

%% 2)

Tempo=15; %duração de 15 segundos

%Leitura dos ficheiros de audio
[d_h,fs_H] = audioread('HouseCuddy.mp3');
[x_h,fs_H] = audioread('RSHouse.mp4');
%A frequencia de amostragem é igual para ambos os ficheiros (fs=44100)

% Utilizamos a funçao Decimate de forma a reduzir a taxa de amostragem do sinal 
%por um fator de 3
Noise_1=decimate(d_h,3); 
Noise_2=decimate(x_h,3);
fs_decimate = fs_H/3; 

%Os seguintes valores foram escolhidos como tentativa e erro
%Compr = 275;
%mu_1 = 0.03;
%mu = [1, 0.1, 0.01, 0.001];

HC= Noise_1(1:Tempo*fs_H);
RSH = Noise_2(1:Tempo*fs_H);

lms_2= dsp.LMSFilter(10,'StepSize', 0.01);  %com o valor 10, diminui o eco
[y_H,e_H] = lms_2(RSH', HC');

soundsc(e_H,fs_H);
t_House=(0:(length(HC)-1))/fs_H;

%Ver o gráfico para o valor escolhido
figure
plot(t_House, RSH)
hold on
plot(t_House, e_H)
hold off
legend('Musica - RSHouse', 'HouseCuddy com música')
soundsc(e_H,fs_decimate)

%% PART III - Encryption

clear all
close all

IHouse = imread('House8000.png');

figure(1)
imagesc(IHouse); axis square
title('Imagem Original');

% Componentes RGB (Red, Green & Blue) da imagem:
 
R = double(IHouse(:,:,1)); % Componente Vermelha 
G = double(IHouse(:,:,2)); % Componente Verde 
B = double(IHouse(:,:,3)); % Componente Azul

% Imagem da componente Red 
figure
colormap gray
subplot (1,3,1);
imagesc(R); axis square
title('Componente Vermelha');
% Imagem da componentes Green 
subplot (1,3,2);
imagesc(G); axis square
title('Componente Verde');
% Imagem da componentes Blue 
subplot (1,3,3);
imagesc(B); axis square
title('Componente Azul');
% Ao separarmos a imagem pelas componentes RGB, verificamos que na
% componente vermelha, artefactos, sendo estes inexistentes para as componentes
%azul e verde.


%% Histograms

figure;
subplot(3,1,1);
histogram(R);
title('Componente Vermelha')
subplot(3,1,2);
histogram(G);
title('Componente Verde')
subplot(3,1,3);
histogram(B);
title('Componente Azul')

 %% Descobrir a música incriptada na imagem

% Transfoma-se a matriz da componente vermelha de 681*1024=697344 elementos 
% em uma matriz de 697344 linhas por uma coluna   
NewR = reshape(R,697344,1);

fs = 8000; % Este valor foi escolhido por ser o nome da imagem dada
soundsc(NewR,fs);


% Ao reproduzir-se o som da matriz R à frequência de amostragem de 8000Hz, 
% descobre-se que a musica "you cant always get what you want" está 
% incriptada na imagem.  