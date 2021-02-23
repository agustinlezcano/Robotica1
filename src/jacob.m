
function [path, Pos] = jacob(R)
%Jacobiano
q = zeros(6,10000);
q(1,:) = (randi([-185,  185], [1,10000]));
q(2,:) = (randi([-155,  95], [1,10000]));
q(3,:) = (randi([-85, 128], [1,10000]));
q(4,:) = (randi([-165,  165], [1,10000]));
q(5,:) = (randi([-115,  140], [1,10000]));
q(6,:) = (randi([-350,  350], [1,10000]));

q = q * pi/180;
h=1;
Ja=[];
for i=1:10000
    i
    %Jacobiano y determinante de 10000 posiciones articuares aleatorias
    J = (R.jacob0(q(:,i)));
    Ja(:,:,i)=J;
    DJ(i) = (det(J));
    if abs(DJ(i))<0.0009
        %Se almacena el indice de las posiciones con det menor al umbral.
        Pos(h)=i;
        h = h + 1;
    end
end

fprintf('Cantidad de posiciones singulares encontradas:');
h
fprintf('Valor minimo de los determinates analizados');
min(abs(DJ))
fprintf('Valor maximo de los determinates analizados');
max(abs(DJ))

path = zeros(6,5);
for i=1:5
    n = randi([1 length(Pos)]);
    path(:,i) = q(:,Pos(n));             %Poner Pos(n) porque sino toma posiciones que son validas
end

end