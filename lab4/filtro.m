function y = filtro(x, A, B)
    % A função do filtro determina a saída y para uma entrada x, 
    % tendo em conta um sistema linear com uma função de transferência 
    % digital H(z) (Eq 1) conforme especificado pelos coeficientes A e B. 

    % (Eq 1) H(z) = (sum(k=0, K) b(k)*z^(-n))/(sum(l=0,L) a(l)*z^(-n))

    % A transformada Z é uma operação digital análoga à transformada de Laplace
    % no domínio analógico, pelo que a resposta de um filtro quando aplicado um
    % input X(z) é dada por:
    % (Eq 2) Y(z) = H(z)X(z) 

    % Aplicando a traformada Z inversa a Y(z), tem-se que:
    % (Eq 3) y(n) = sum(k=0, K) b(k)x(n-k) - sum(l=0,L)a(l)y(n-l)

    % Enquanto os coeficientes b(n) operam apenas nos valores de entrada x(n), 
    % os coeficientes a(n) operam em valores passados da saída y(n) e são, portanto, 
    % às vezes referidos como coeficientes recursivos.

    L = length(A); % nº de coeficientes a(l) -> ordem do denominador 
    K = length(B); % nº de coeficientes b(k) -> ordem do numerador  
	
	% Se a(1) não for igual a 1, o filtro normaliza os coeficientes do filtro 
    % por a(1). Portanto, a(1) deve ser diferente de zero.
    if(A(1) ~= 1)
        A = A ./ A(1);
        B = B ./ A(1);
    end
    
    % inicialização da saída y(n), tendo a mesma dimensão da entrada x(n) 
    y = zeros(size(x));
    
    % Os dados de saída filtrados não começam até que o número n de pontos 
    % de dados passados pelo filtro seja maior que K, criando-se um delay
    % da resposta de saída y face à entrada x, isto é, o sinal de saida é
    % shifted no tempo em relação ao sinal de entrada. 

    % De notar que a devido à indexação no matlab (começa no indice 1 e não 0), a Eq 3
    % é alterada nesse sentido. 
    for n = K:length(x)
		y(n)=B'*x(n:-1:n-K+1)-A(2:end)'*y(n-1:-1:n-L+1);
    end
   
    % observe que se os coeficientes a(l) na Eq 3 são zero (com exceção 
    % de a(0) = 1), Eq 3 reduz-se a convolução (filtros FIR, sendo os outros 
    % os filtros IIR)
end
