function y = parameterfun(C,X,Y,Z,En)

A=zeros(1,numel(X));

for i = 1:numel(X)
A(i)= 1/sqrt((C(1)-X(i))^2+(C(2)-Y(i))^2+(C(3)-Z(i))^2);
end

y = (A*En*En'*A')/(A*A');

