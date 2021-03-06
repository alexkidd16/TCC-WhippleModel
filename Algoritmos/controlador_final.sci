funcprot(0)
clear; clc; clearglobal;
global phi
global delta
global r

function [degrees] = radians2degrees(radians)  //Utilizado para imprimir a saida em graus
    degrees = radians*(180/%pi); 
endfunction 

//Auto-valores do ciclo fechado
s = [ -1, -1, -1+(%i*1),  -1+(-%i*1) ]
//s = [ -1, -1, -1+(%i*20),  -1+(-%i*20) ]
//s = [ -20, -20, -1+(%i*1),  -1+(-%i*1) ]
//s = [ -1, -1, -20+(%i*1),  -20+(-%i*1) ]
//s = [ 1, 1, -1+(%i*1),  -1+(-%i*1) ]

//Variaveis já utilizadas e calculadas no outro script
v = 5
g = 9.81
c_base = 0.08
lambda = %pi/10
w = 1.02

M = [80.81722, 2.3194133; 2.3194133, 0.29784189]
c = [0, 33.866414; -0.8503564, 1.685404]
K0 = [-80.95, -2.5995169; -2.5995169, -0.8032949]
K2 = [0, 76.597346; 0, 2.6543152]


//Criando a representação dos estados de espaço
invM = inv(M)
KV = - invM * (g * K0 + K2 * v^2)
CV = - invM * c * v

A = eye(4,4)
A(1,:) = [0,0,1,0]
A(2,:) = [0,0,0,1]
A(3,:) = [KV(1,1), KV(1,2), CV(1,1), CV(1,2)]
A(4,:) = [KV(2,1), KV(2,2), CV(2,1), CV(2,2)]

B = [0; 0; invM(2,1); invM(2,2)]
C = [0,1,0,0]

Wr = [B A*B A^2*B A^3*B]
if det(Wr) == 0 then
    disp("O sistema não é controlável pois a matriz Wr tem determinante igual 0")
else
    disp("O sistema é controlável.")
end

//Calculo dos ganhos do controlador, dados ´pr det(sI-A-BK)
K = ppol( A, B, s)

//variaveis da simulação
t0 = 0; tmax=20;
h = 0.1
t = t0:h:tmax
x0 = [0;0;0;0;0;0]
Kr = -1 / ( C * inv(A-B*K) * B )
r = 0

//Criação do caminho a ser percorrido
target(1)=0
for i=2:length(t)
    if i <= 10 then
        target(i) = 0
    elseif i < 40 then
        target(i) = target(i-1) + 0.02
    else
        target(i) = target(i-1)
    end
end


function ydot=f(t,y)
    global phi; global delta; global r;
    
    i=t*10+1
    if i > 201 then
        i=201
    end
    
    // olhando para a frente
    if( i < (tmax*10-20)) then    
        targ = target(i+20)
    else
        targ = target(tmax*10)
    end
    
    // atribuindo valor para o r, que é o valor que queremos que o delta atinja
    if t < 1 then
        r = targ
    elseif modulo(int(i), 2) == 0 then
        erro = targ-y(6)
        r = y(2)+erro
    end
    
    //controlador
    yTemp = [y(1) y(2) y(3) y(4)]'
    u=-K*yTemp+Kr*r
    //integral
    ydot(1:4) = A*yTemp + B*u
    ydot(5) = (v*ydot(2)+c_base*ydot(4))*(cos(lambda)/w)
    ydot(6) = ydot(5)*v*h
    
    //Guarda o valor de phi e delta, pois a saida y(1) e y(2) é o resultado
    //da integral, e portanto, é uma somatoria de todos os valores ate aquele momento
    phi(i) = ydot(1)
    delta(i) = ydot(2)
endfunction
y = ode(x0, t0, t, f)


//Desenhando o grafico.
clf();
plot2d(t,y(6,:), style=color("red"))
plot2d(t,target, style=color("blue"))
xlabel("Tempo (s)")
ylabel("Distância (m)")
title("s = [ -1, -1, -1+i,  -1-i ]")

b=newaxes()
b.filled="off"
xset("line style", 3)
plot2d(t,radians2degrees(delta), axesflag=3, style=color("green"))
plot2d(t,radians2degrees(phi), axesflag=3, style=color("magenta"))
xset("line style", 1)
ylabel("Graus (º)")
legends(["$Posição\;da\;bicicleta$";"$Caminho\;traçado$";'$\delta\;delta\;(^\circ)$';'$\varphi\;phi\;(^\circ)$'], [[5;1],[2;1],[3;3],[6;3]] ,opt="lr")
