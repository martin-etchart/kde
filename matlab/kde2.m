function [bandwidth,density,xmesh,cdf]=kde2(data,n,MIN,MAX)
% Reliable and extremely fast kernel density estimator for one-dimensional data;
%        Gaussian kernel is assumed and the bandwidth is chosen automatically;
%        Unlike many other implementations, this one is immune to problems
%        caused by multimodal densities with widely separated modes (see example). The
%        estimation does not deteriorate for multimodal densities, because we never assume
%        a parametric model for the data.
% INPUTS:
%     data    - a vector of data from which the density estimate is constructed;
%          n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12;
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where Range=max(data)-min(data);
% OUTPUTS:
%   bandwidth - the optimal bandwidth (Gaussian kernel assumed);
%     density - column vector of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%             - If no output is requested, then the code automatically plots a graph of
%               the density estimate.
%        cdf  - column vector of length 'n' with the values of the cdf
%  Reference: 
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957. 

%
%  Example:
%           data=[randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%              kde(data,2^14,min(data)-5,max(data)+5);


%  Notes:   If you have a more reliable and accurate one-dimensional kernel density
%           estimation software, please email me at botev@maths.uq.edu.au

opt.show_plots=1;

data=data(:); %make data a column vector
if nargin<2 % if n is not supplied switch to the default
    n=2^14;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;
if nargin<4 %define the default  interval [MIN,MAX]
    minimum=min(data); maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/10; MAX=maximum+Range/10;
end
% set up the grid over which the density estimate is computed;
R=MAX-MIN; dx=R/(n-1); xmesh=MIN+[0:dx:R]; N=length(unique(data));
%bin the data uniformly using the grid defined above;
initial_data=histc(data,xmesh)/N;  initial_data=initial_data/sum(initial_data);
a=kde_dct1d(initial_data); % discrete cosine transform of initial data

%jfc_vector_save_simple(initial_data,'initial_data.txt');

if opt.show_plots
    figure(1)
    plot(data,'-b')
    hold on
    if(exist('data_fftw'))
    plot(data_fftw,'-r')
    norm(data-data_fftw)
    end
    hold off
end

if opt.show_plots
    figure(13)
    plot(xmesh,'-b')
    hold on
    if(exist('xmesh_fftw'))
    plot(xmesh_fftw,'-r')
    norm(xmesh-xmesh_fftw)
    end
    hold off
end

if opt.show_plots
    figure(11)
    plot(initial_data,'-b')
    hold on
    if(exist('initial_data_fftw'))
    plot(initial_data_fftw,'-r')
    norm(initial_data-initial_data_fftw)
    end
    hold off
end

if opt.show_plots
    figure(12)
    plot(a,'-b')
    hold on
    if(exist('a_fftw'))
    plot(a_fftw,'-r')
    norm(a-a_fftw')
    end
    hold off
end

% now compute the optimal bandwidth^2 using the referenced method
I=[1:n-1]'.^2; a2=(a(2:end)/2).^2;

tt=kde_fixed_point(0.01,N,I,a2)
tt=kde_fixed_point(0.0,N,I,a2)
tt=kde_fixed_point(0.1,N,I,a2)



%plot fixed point
tin=[-0.1:1e-5:0.1];
lt=length(tin);
tout=zeros(1,lt);
figure(15)
for i=1:lt
    ti=tin(i);
    tout(i)=kde_fixed_point(ti,N,I,a2);
end
plot(tin,tin, '-b')
hold on
plot(tin,tout, '-r')
hold off
figure(16)
plot(tin,tout-tin)



% use  fzero to solve the equation t=zeta*gamma^[5](t)
try
    options=optimset('fzero');
    t_star=fzero(@(t)kde_fixed_point(t,N,I,a2),[0,.1],options);
catch
    t_star=.28*N^(-2/5);
end
%just for debug with fixed values
%t_star=.28*N^(-2/5);
%t_star=1.04904e-06;


% smooth the discrete cosine transform of initial data using t_star
a_t=a.*exp(-[0:n-1]'.^2*pi^2*t_star/2);


if opt.show_plots
   figure(2)
    plot(a2,'-b')
    hold on
    if(exist('a2_fftw'))
    plot(a2_fftw,'-r')
    norm(a2-a2_fftw)
    end
    hold off
end

if opt.show_plots
   figure(21)
    plot(a_t,'-b')
    hold on
    if(exist('a_t_fftw'))
    plot(a_t_fftw,'-r')
    norm(a_t-a_t_fftw)
    end
    hold off
end


% now apply the inverse discrete cosine transform
if (nargout>1)|(nargout==0)
    density=kde_idct1d(a_t)/R;
end
% take the rescaling of the data into account
bandwidth=sqrt(t_star)*R;
if nargout==0
    figure(1), plot(xmesh,density)
end
% for cdf estimation
if nargout>3
    f=2*pi^2*sum(I.*a2.*exp(-I*pi^2*t_star));
    t_cdf=(sqrt(pi)*f*N)^(-2/3);
    % now get values of cdf on grid points using IDCT and cumsum function
    a_cdf=a.*exp(-[0:n-1]'.^2*pi^2*t_cdf/2);
    cdf=cumsum(kde_idct1d(a_cdf))*(dx/R);
    % take the rescaling into account if the bandwidth value is required
    bandwidth_cdf=sqrt(t_cdf)*R;
end

if opt.show_plots
   figure(22)
    plot(xmesh,density,'-b')
    hold on
    if(exist('density_fftw'))
    plot(xmesh,density_fftw,'-r')
    norm(density-density_fftw)
    end
    hold off
end

keyboard
end









