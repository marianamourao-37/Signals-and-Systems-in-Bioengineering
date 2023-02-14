function [e,y,M_coeff] = adaptFilter(d,x,p,w)

% p - ordem do filtro
% x - sinal de referencia 
% d - sinal desejado (normalmente corrompido com ruido) 
% w - dimensoa da janela 

if nargin < 4 % se forem entregues menos de 4 inputs, 
    % aplica-se o primeiro criteiro de otimização 
    w = 0; % dimensao da janela nula 
end

N = length(x); % comprimento dos sinais definidos a partir do 
% sinal de referencia 

y = zeros(size(x));
e = zeros(size(x));
M_coeff = zeros(p+1,length(x));
A = zeros(w+1,p+1);

for n = p+1+w:N
    a = x(n:-1:n-p)'; 
    
    for i=0:w
        A(i+1,:)=x(n-i:-1:n-i-p)';
    end
    
    D = d(n:-1:n-w);
    h=pinv(A)*D;
    
    M_coeff(:,n) = h;
    
    y(n) = a*h;
    e(n) = d(n)-y(n);
end
end