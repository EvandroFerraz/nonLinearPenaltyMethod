%Método de Penalidades
clear all
clc

%beta = 3, mi = 100, x0 = [5, 5] e E = 0.0001

%Dados de entrada
fprintf('---Método das Penalidades com Newton--- \n');
a = input('Digite o ponto inicial x1: ');
b = input('Digite o ponto inicial x2: ');
eps = input('Digite o valor limite de tolerância: ');
mi = input('Digite o parâmetro de penalidade de mi: ');
beta = input('Digite o valor de beta: ');
k = 1;

%Variáveis e Função objetivo
syms x1 x2;
x = [x1,x2];

f = (x1-1)^2 +(x2-2)^2;
v(k,:) = [a,b];
 
 while 1
    penalidade = (mi/2)*(x2-x1-1)^2 + (mi/2)*(x1+x2-2)^2; %calculo da penalidade
    fpen = f + penalidade;  %função irrestrita
    
    %método de newton
    %define o gradiente de fpen
    for m = 1:2
        gradf(m) = diff(fpen,x(m));
    end
    
    %define a matriz hessiana de fpen
    for i = 1:2
        for j = 1:2
            H(i,j) = diff(gradf(i),x(j));
        end
    end
    
    fpenEval(k) = subs(fpen,x,v(k,:));
    
    %valida o gradiente 
    gradfEval(k,:) = subs(gradf,x,v(k,:));
    
    %matriz hessiana evaluada
    h(:,:,k) = subs(H,x,v(k,:)); 
    
    %inversa da hessiana
    invH(:,:,k) = inv(h(:,:,k));
    
    %atualização de xk
    v(k+1,:) = v(k,:) - (invH(:,:,k)*gradfEval(k,:)')';
    
    penEval = subs(penalidade,x,v(k,:)); %validade a função de penalidade
    
    if(penEval < eps) %critério de parada
        break;
    end
    
    k = k+1; 
    mi = mi * beta; %incrementa o valor de penalidade
 end
 clc
 %Resultados
 fprintf('Pontos [x1,x2]: \n');
 v(k+1,:)
 fprintf(1, 'Valor da Função = %6.4f \n', subs(f,x,v(k+1,:)));
 fprintf(1, 'Valor de k = %6.4f \n', k);
