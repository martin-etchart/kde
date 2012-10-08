function y = kde_dct_with_matlab( x )

N=length(x);
x_n=x*sqrt(2*N);
y=dct(x_n);
y(1)=y(1)/sqrt(2);

end

