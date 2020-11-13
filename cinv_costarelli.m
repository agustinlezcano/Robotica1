clc
clear
close all

dh = [0 0.520 .160 pi/2 0;
    0 0 0.980 0 0;
    0 0 0 pi/2 0;
    0 0.8545 0 -pi/2 0;
    0 0 0 pi/2 0;
    0 0.1605 0 0 0];

%matriz de prueba
T=[1 0 0 .9;
    0 1 0 .5;
    0 0 1 1.43;
    0 0 0 1];

R=SerialLink(dh);
q0 = [0 pi/4 pi/3 0 0 0]
T = R.fkine(q0);
T = T.double

Pm = T(1:3,4) - dh(6,2) * T(1:3,3); %posicion de la muñeca resto parametro d6

%8 soluciones
q01=[0,0,0,0,0,0];
q02=[0,0,0,0,0,0];
q03=[0,0,0,0,0,0];
q04=[0,0,0,0,0,0];
q05=[0,0,0,0,0,0];
q06=[0,0,0,0,0,0];
q07=[0,0,0,0,0,0];
q08=[0,0,0,0,0,0];

q1a=atan2(Pm(2),Pm(1)); %q1 normal
q1b=q1a+pi;              %q1 opuesto

a2=dh(2,3);
a3=dh(4,2);

%% geometrico
%primera parte problema pieper

r0 = sqrt(Pm(1)^2 + Pm(2)^2);

r01 = r0 - dh(1,3);
%r02=r0+dh(1,3); %cuando q1 es sentido opuesto el parametro dh lo aleja mas
 
z01 = Pm(3) - dh(1,2); %altura triangulo medido desde el sistema 1

hip1 = sqrt(r01^2 + z01^2); %hip del triangulo desde sistema 1 hasta comienzo de muñeca

%hip2=sqrt(r02^2+z01^2);

beta1 =  atan2(z01, r01);

%beta2=atan2(z01,r02);


alfa1=real(acos((a2^2+hip1^2-a3^2)/(2*a2*hip1))); %q1 normal
%alfa2=real(acos((a2^2+hip1^2-a3^2)/(2*a2*r02))); %q1 opuesto

q2a1=beta1+alfa1; %codo arriba para q1 normal
q2b1=beta1-alfa1;%codo abajo para q1 normal

%q2a2=beta2+alfa2; %codo arriba para q1 opuesto
%q2b2=beta2-alfa2; %codo abajo para q1 opuesto

w=real(acos((a2^2+a3^2-hip1^2)/(2*a2*a3)));
q3=pi-w;


q3a=pi/2-q3;%codo arriba
q3b=pi/2+q3;%codo abajo

%codo arriba q1+
q01(1)=q1a;
q01(2)=q2a1;
q01(3)=q3a;

%codo abajo q1+
q02(1)=q1a;
q02(2)=q2b1;
q02(3)=q3b;

R.plot(q0)
%R.teach;
hold on
trplot(T);

rad2deg(q01)
R.fkine(q01)

rad2deg(q02)
R.fkine(q02)

%% Segunda parte: soluci�n por Pieper
syms q1 q2 q3 q4 q5 q6 real

fprintf('Matriz 1 respecto a 0:\n')
T1 = R.links(1).A(q01(1)).double;
T1
fprintf('\nMatriz 2 respecto a 1:\n')
T2 = R.links(2).A(q01(2)).double;
T2
fprintf('\nMatriz 3 respecto a 2:\n')
T3 = R.links(3).A(q01(3)).double;
T3

fprintf('\nMatriz 4 respecto a 3:\n')
T4 = R.links(4).A(q4).double;
T4
fprintf('\nMatriz 5 respecto a 4:\n')
T5 = R.links(5).A(q5).double;
T5
fprintf('\nMatriz 6 respecto a 5:\n')
T6 = R.links(6).A(q6).double;
T6


%T06 = T03 * T36 -> T36 = inv(T03) * T06
T03 = T1*T2 *T3;
T36_a = inv(T03) * T;
T36_b = simplify(T4 * T5 *T6);

[(1:16)',T36_a(:), T36_b(:)]

%De la ecuacion 15 se puede despejar q5 y luego tratar de sacar q4 con la
%ecuacion 10. SI no ver de despejar (isolate o parecido) un seno o coseno e igualar y ver si
%anda o no.

solve(T36_a == T36_b, q5)