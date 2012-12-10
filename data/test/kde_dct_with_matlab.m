function y = kde_dct_with_matlab( x )

n=length(x);

x_n=x*sqrt(2*n);
y=dct(x_n);
y(1)=y(1)/sqrt(2);

end

