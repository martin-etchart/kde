%################################################################
function  out=kde_fixed_point(t,N,I,a2)
verbose=0;

if verbose, bones_xml('in'); end

% this implements the function t-zeta*gamma^[l](t)
l=7;
sum_f=sum(I.^l.*a2.*exp(-I*pi^2*t));
sum_f_2=sum(I.^l.*exp(-I*pi^2*t));
f=2*pi^(2*l)*sum_f;

if verbose, 
    bones_log('t: %g\n',t);
    bones_log('sum_f_2: %g sum_f: %g f: %g\n',sum_f_2,sum_f,f);
end

for s=l-1:-1:2
    K0=prod([1:2:2*s-1])/sqrt(2*pi);  
    const=(1+(1/2)^(s+1/2))/3;
    tim=(2*const*K0/N/f)^(2/(3+2*s));
    sum_f=sum(I.^s.*a2.*exp(-I*pi^2*tim));
    f=2*pi^(2*s)*sum_f;
    if verbose,
        bones_log('s: %d c: %g k0: %g tim: %g f: %g \n',s,const,K0,tim,f);
        bones_log('sum_f: %g\n',sum_f);
    end
end
out=t-(2*N*sqrt(pi)*f)^(-2/5);


if verbose, bones_xml('out'); end
end

