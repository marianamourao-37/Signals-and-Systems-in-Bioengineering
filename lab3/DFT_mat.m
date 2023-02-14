function X = DFT_mat(x,M)
%Função que recebe como input um vetor coluna x com os valores do sinal discretizado, e
%M como o número de coeficientes da DFT a computar. Devolve como output
%vetor X da DFT de comprimento M 

%comprimento do sinal discretizado (nº de amostras) 
N=length(x);

k = (0:(M-1)); % vetor k de comprimento M
n = (0:(N-1)); % vetor n de comprimento N 

matrix_exp = exp(-1i*2*pi*k'.*n/M); %matriz MxN, que contem, na entrada (k,n), exp(-j*2*pi*k*n/N)
X = matrix_exp*x; %vetor coluna Mx1 com a DFT do sinal x 
end