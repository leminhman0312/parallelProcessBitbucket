clear;
clc;
s = zeros(10,10);


for i = 1:9
  s(i,i+1) = 0.5;
endfor

for i = 2:10
  s(i,i-1) = 0.5;
endfor

s(10,1) = 0.5;
s(9,10) = 1.0;








q = 0.15;
b = zeros(10,10);
numpage = 10;
b(:) = (1.00/numpage);


g = (1-q)*s + q*b;

x = zeros(10,1);
x(:) = 1.00/numpage;



for i = 1:3
  y = g*x
  x = y;
endfor

  





















