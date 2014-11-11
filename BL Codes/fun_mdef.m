function mdefj = fun_mdef(uej,dstarj,xj)

mdefj = ((uej(2:end)-uej(1:end-1)).*(dstarj(2:end)-dstarj(1:end-1)))./(xj(2:end)-xj(1:end-1));

