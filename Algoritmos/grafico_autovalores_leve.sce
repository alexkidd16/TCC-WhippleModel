w = 1.02                            // base da roda (m)
c = 0.08                            // trilha (m)
lambda = %pi/10                     // steer axis tilt (rad)
g = 9.81                            // gravidade (N/kg_1)

//Roda traseira R
rR = 0.3                            // raio da roda (m)
mR = 2                              // massa da roda(kg)
IRxx = 0.0603
IRyy = 0.12                         // momento de inercia (kg/m²)

//Quadro da bicicleta B
[xB,zB] = (0.3, -0.9)               // Centro de massa (m)
mB = 50                             // massa do quadro (kg)
IBxx = 9.2
IBxz = 2.4
IByy = 11
IBzz = 2.8                          // momento de inercia (kg/m²)

//Garfo da bicicleta H
[xH,zH] = (0.9, -0.7)               // Centro de massa (m)
mH = 4                              // massa do garfo (kg)
IHxx = 0.05892
IHxz = -0.00756
IHyy = 0.06
IHzz = 0.00708                      // momento de inercia (kg/m²)

// Roda dianteira F
rF = 0.35                           // raio da roda (m)
mF = 3                              // massa da roda(kg)
IFxx = 0.1405
IFyy = 0.28                         // momento de inercia (kg/m²)

// Coeficientes das equações lineares
mT = mR + mB + mH + mF
xT = (xB * mB + xH * mH + w * mF) / mT
zT = (-rR*mR + zB*mB + zH*mH - rF*mF) / mT
ITxx = IRxx + IBxx + IHxx + IFxx + mR*rR^2 + mB*zB^2 + mH*zH^2 + mF*rF^2
ITxz = IBxz + IHxz - mB*xB*zB - mH*xH*zH + mF*w*rF
IRzz = IRxx
IFzz = IFxx
ITzz = IRzz + IBzz + IHzz + IFzz + mB*xB^2 + mH*xH^2 + mF*w^2
mA = mH + mF
xA = (xH*mH+w*mF)/mA
zA = (zH*mH-rF*mF)/mA
IAxx = IHxx + IFxx + mH*(zH-zA)^2 + mF*(rF+zA)^2
IAxz = IHxz - mH*(xH-xA)*(zH-zA)+mF*(w-xA)*(rF+zA)
IAzz = IHzz + IFzz + mH*(xH-xA)^2 + mF*(w-xA)^2
uA = (xA-w-c)*cos(lambda) - zA * sin(lambda)
IAll = mA*uA^2 + IAxx*sin(lambda)^2 + 2*IAxz * sin(lambda) * cos(lambda) + IAzz*cos(lambda)^2
IAlx = -mA*uA*zA + IAxx*sin(lambda) + IAxz*cos(lambda)
IAlz = mA*uA*xA + IAxz*sin(lambda) + IAzz*cos(lambda)
mi = (c/w) * cos(lambda)
SR = IRyy/rR
SF = IFyy/rF
ST = SR+SF
SA = mA*uA + mi*mT*xT
M = [ITxx, IAlx+mi*ITxz; IAlx+mi*ITxz, IAll+2*mi*IAlz+mi^2*ITzz]
K0 = [mT*zT, -SA; -SA, -SA*sin(lambda)]
K2 = [0, ((ST - mT*zT)/w)*cos(lambda); 0, ((SA+SF*sin(lambda))/w)*cos(lambda)]
C1 = 0
C2 = mi*ST+SF*cos(lambda)+(ITxz/w)*cos(lambda)-mi*mT*zT
C3 = -(mi*ST + SF*cos(lambda))
C4 = (IAlz/w)*cos(lambda) + mi*( SA+ (ITzz/w)*cos(lambda))
c = [C1, C2; C3, C4]

invM = inv(M)

v = 0                                               // velocidade (m/s)
x = 1
while v <= 10
    KV = - invM * (g * K0 + K2 * v^2)
    CV = - invM * c * v

    A = eye(4,4)
    A(1,:) = [0,0,1,0]
    A(2,:) = [0,0,0,1]
    A(3,:) = [KV(1,1), KV(1,2), CV(1,1), CV(1,2)]
    A(4,:) = [KV(2,1), KV(2,2), CV(2,1), CV(2,2)]
    
    B = [0; 0; invM(2,1); invM(2,2)]
    C = [0,1,0,0]
    
    //Verificado os autovalores do sistema.
    a = spec(A)
    if v == 0
        auto(1,x)= a(3) 
        auto(2,x)= a(4)
        auto(3,x)= a(2)
        auto(4,x)= a(1)
    elseif x < 5
        auto(1,x)= a(1) 
        auto(2,x)= a(2)
        auto(3,x)= a(4)
        auto(4,x)= a(3)        
    elseif x < 40 then
        auto(:,x) = a
    else
        auto(1,x)= a(2) 
        auto(2,x)= a(3)
        auto(3,x)= a(1)
        auto(4,x)= a(4)
    end
    
    v=v+0.1
    x=x+1
end

vel = 0:0.1:10
hor = -16:0.1:6
clf
plot(vel,auto(1,:),'r')
xrects([3.75,6,1.31,22]', 7 )
plot(vel,auto(1,:),'r')
plot(vel,auto(2,:),'r')
plot(vel,imag(auto(1,:)),'-.r')
plot(vel,imag(auto(2,:)),'-.r')
plot(vel,auto(3,:),'b')
plot(vel,auto(4,:),'g')
plot(vel,0, 'black')
xrect(3.75,6,1.31,22)
xlabel('Velocidade(m/s)')
ylabel('Auto-valores')
title('Bicicleta leve')
mtlb_axis([0, 10, -10, 5])
