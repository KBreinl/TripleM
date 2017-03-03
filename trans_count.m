function [matx, rw, st, ed] = trans_count(d, Order)
% TRANS_COUNT
% Copyright (c) 2014, Shawn Pethel
% https://se.mathworks.com/matlabcentral/fileexchange/40188-whittle-surrogate/content/trans_count.m
%[matx rw st ed] = trans_count(d,Order)
%
% Given a sequence of integers trans_count 
% generates a Markov transition count matrix of
% order Order, a list of words, and a beginning
% and start word. The outputs are intended for use
% with whittle_surrogate.m
%
% INPUTS: 
% d:        a symbol sequence
% Order:    the Markov order of the sequence
%
% OUPUTS:
% matx:     a sparse matrix of transitions counts
% rw:       a list of words
% st:       index of the first word in d -> rw(st)
% ed:       index of the last word in d -> rw(ed)
%
% EXAMPLE:
% >> d = [0 1 1 0 1 0 1 1 1 0 0 1];
% >> [f w u v] = trans_count(d,1);
% >> full(f)
% ans =
%    1     4
%    3     3
% w
% w =
%       0
%       1
%
% There is one '0' to '0' word transition, four '01's, three '10's, 
% and three '11's in d
%
% Shawn Pethel, 2012

d = d(:);
L = Order + 1;
[rw, c] = unique_count(d);
if L == 1
    matx = c;
    [~, st] = ismember(d(1), rw);
    [~, ed] = ismember(d(end), rw);
    return
end
%Find all L length sequences
a = embed(d, L);
[rwL, bL] = unique_count(a);
%Find all L-1 length sequences
a = embed(d, L-1);
[rw, ~] = unique_count(a);
if L == 2
    [~, st] = ismember(d(1), rw);
    [~, ed] = ismember(d(end), rw);
else
    [~, st] = ismember(d(1:L - 1)', rw, 'rows');
    [~, ed] = ismember(d(end - L + 2:end)', rw, 'rows');
end
indx = (1:size(rwL, 1));
e1 = zeros(size(rwL, 1), 1);
e2 = zeros(size(rwL, 1), 1);
% Store the transistion table as a sparse matrix
for i = 1:size(rwL, 1)
    y = rwL(i * ones(1, size(rw,1)), 1:L - 1);
    c1 = all(rw == y, 2);
    y = rwL(i * ones(1, size(rw, 1)), 2:L);
    c2 = all(rw == y, 2);
    e1(i) = indx(c1);
    e2(i) = indx(c2);
end
matx = sparse(e1, e2, bL, size(rw, 1), size(rw, 1), length(bL));


function [y, c] = unique_count(x)
% Unique_count sorts and counts the rows of x
% y lists the unique rows and c their counts
if size(x, 1) == 1
    x = x(:);
end
y = sortrows(x);
d = any(y(1:end - 1 , :) - y(2:end, :), 2);
d = [d; true];
ind = (1:size(y, 1));
g = ind(d);
n = length(g);
c(1) = g(1);
c(2:n) = g(2:n) - g(1:n - 1);
c = c';
y = y(g, :);


function z = embed(x, d)
% rearranges a one dimensional array x into
% overlapping rows of length d
x = x(:);
N = length(x);
i = (1:d);
j = (0:N - d)';
a = i(ones(N - d + 1, 1), :);
b = j(:,ones(1, d));
z = x(a + b);
