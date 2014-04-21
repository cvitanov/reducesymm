function [ func, gradfunc] = linkfunc( symbolseq )
w = 0.3;
r = 1;

% symbolseq can be cell structure

% symbols and number
% a b c d e f  A B C D E F
% 0 2 4 6 8 10 1 3 5 7 9 11

nx = length(symbolseq);

Rh = zeros(2,11);

% nearest neighbor
for jj = 0:2:10
    ang = jj*pi/6;
    Rh(:, jj+1) = (2*r+w)*[cos(ang);sin(ang)];
end

% next nearest neight
for jj = 1:2:11
    ang = (jj-1)*pi/6 + pi/6;
    Rh(:, jj+1) = sqrt(3)*(2*r+w)*[cos(ang);sin(ang)];
end

tmpstr = '';
for ii = 1:nx
    number = symbolseq(ii)+1;
    rrvec1 = strcat('[cos(th(',num2str(ii),'));sin(th(',num2str(ii),'))]*', num2str(r));
    rrvec2 = strcat('[cos(th(',num2str(ii+1),'));sin(th(',num2str(ii+1),'))]*', num2str(r));
    tmpstr = strcat(tmpstr,'veclen(Rh(:,',num2str(number),')-', rrvec1, '+', rrvec2, ');');
end
tmpstr = strcat('[',tmpstr,']');

funcv = [];
eval(strcat('funcv = @(th) ',tmpstr,';'));

func = @(th) sum(funcv([th,th(1)]));
% func = @(th) funcv([th,th(1)]);


%% now we will explicitly evaluate the gradient function and compare with numerical results

% first row [sin(th1), -cos(th1)]*(tmpl1xy/l1-tmplnxy/ln);
% the ith row is given by [sin(thi),-cos(thi)]*(tmplixy/li-tmpl(i-1)xy/l(i-1))

gradstr = '';
lastlxy = '';
lastlen = '';
for ii = 1:nx
    number = symbolseq(ii)+1;
    thvecn = strcat('[sin(th(',num2str(ii),')),-cos(th(',num2str(ii),'))]');
    rrvec1 = strcat('[cos(th(',num2str(ii),'));sin(th(',num2str(ii),'))]*', num2str(r));
    rrvec2 = strcat('[cos(th(',num2str(ii+1),'));sin(th(',num2str(ii+1),'))]*', num2str(r));
    if ii == nx
        rrvec2 = strcat('[cos(th(',num2str(1),'));sin(th(',num2str(1),'))]*', num2str(r));
    end
    tmplxy = strcat('Rh(:,',num2str(number),')-', rrvec1, '+', rrvec2);
    tmplen = strcat('veclen(',tmplxy,')');
    tmpval = strcat(thvecn,'*((',tmplxy,')/',tmplen,'-(',lastlxy,')/',lastlen,');');
    if ii >=2
        gradstr = strcat(gradstr,tmpval);
    end
    lastlxy = tmplxy;
    lastlen = tmplen;
end

number = symbolseq(ii)+1;
thvecn = strcat('[sin(th(1)),-cos(th(1))]');
rrvec1 = strcat('[cos(th(1));sin(th(1))]*', num2str(r));
rrvec2 = strcat('[cos(th(2));sin(th(2))]*', num2str(r));
tmplxy = strcat('Rh(:,',num2str(number),')-', rrvec1, '+', rrvec2);
tmplen = strcat('veclen(',tmplxy,')');
tmpval = strcat(thvecn,'*((',tmplxy,')/',tmplen,'-(',lastlxy,')/',lastlen,');');
gradstr = strcat('[',tmpval,gradstr,']');

gradfunc = [];
eval(strcat('gradfunc = @(th)', gradstr, ';'));

    function [ len ] = veclen(tmpvec)
        len = sqrt(dot(tmpvec, tmpvec));
    end

end

