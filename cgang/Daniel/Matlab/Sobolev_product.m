% Sobolev_product.m
% by Daniel Borrero 08/3/2012
% ----------------------------------------------------------------
% Calculates the inner product of two vectors for the 2-mode
% system of Porter & Knobloch using the Sobolev H^1 norm whose 
% metric tensor is given by
%                   [[1 0 0 0];
%               g =  [0 1 0 0];
%                    [0 0 4 0];
%                    [0 0 0 4]]

function n = Sobolev_product(x,y)

% Calculate inner product of two 4D vectors x and y using the H^1 Sobolev
% norm
n = x(1)*y(1) + x(2)*y(2) + 4*x(3)*y(3) + 4*x(4)*y(4);