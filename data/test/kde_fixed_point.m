%################################################################
function  out=kde_fixed_point(t,N,I,a2)

% this implements the function t-zeta*gamma^[l](t)
l=7;
sum_f=sum(I.^l.*a2.*exp(-I*pi^2*t));
f=2*pi^(2*l)*sum_f;
%fprintf('f: %g\n',f);
%fprintf('sum_f: %g\n',sum_f);
for s=l-1:-1:2
    K0=prod([1:2:2*s-1])/sqrt(2*pi);  
    const=(1+(1/2)^(s+1/2))/3;
    time=(2*const*K0/N/f)^(2/(3+2*s));
    sum_f=sum(I.^s.*a2.*exp(-I*pi^2*time));
    f=2*pi^(2*s)*sum_f;
    %fprintf('s: %d c: %g k0: %g tim: %g f: %g \n',s,const,K0,time,f);
    %fprintf('sum_f: %g\n',sum_f);
end
out=t-(2*N*sqrt(pi)*f)^(-2/5);

end

