function [y,e,m_coeff] = LMSfunction(x,d,mu,p)

% Inputs: 
% x = sinal de referencia, normalmente ruido 
% d = sinal desejado, normalmente corrompido por ruido 
% mu = parametro step-size que controla as caracteristicas de 
% convergência do algoritmo
% p = comprimento ou ordem do filtro FIR

% Outputs:
% m_coeff = matriz contendo em cada coluna a resposta impulsiva estimada 
% y = ALE output
% e = erro residual 

N = length(x); % comprimento dos sinais 
% inicialização de outputs 
h = zeros(p,1);
y = zeros(N,1);
m_coeff = zeros(p,N);

for n = p:N
    a = x(n:-1:n-p+1); % Isolate reversed signal segment
    y(n) = h'*a; % Convolution
    e(n)=d(n)-y(n); % Calculate error signal
    
    jacb = -2*e(n)*a;
    
    h = h - mu*jacb; % Adjust filter coefficients. Eq. 8.17
    m_coeff(:,n)=h;
    
end
e = e';
end