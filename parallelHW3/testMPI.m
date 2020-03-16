clear;
clc


a = zeros(10,10);

x = [1;2;3;4;5;1;2;3;4;5];


for j = 1:10
  for i = 1:5
    a(i,j) = 55;
  endfor
endfor


for j = 1:10
  for i = 6:10
    a(i,j) = 99;
  endfor
endfor

y = a*x