function [Iseg,sep,thr] = otsu(I,n)
%OTSU Gray-level image segmentation using Otsu's method.
%   Iseg = OTSU(I,n) computes a segmented image (Iseg) containing n classes
%   by means of Otsu's n-thresholding method (Otsu N, A Threshold Selection
%   Method from Gray-Level Histograms, IEEE Trans. Syst. Man Cybern.
%   9:62-66;1979). Thresholds are computed to maximize a separability
%   criterion of the resultant classes in gray levels.
%
%   OTSU(I) is equivalent to OTSU(I,2). By default, n=2 and the
%   corresponding Iseg is therefore a binary image. The pixel values for
%   Iseg are [0 1] if n=2, [0 0.5 1] if n=3, [0 0.333 0.666 1] if n=4, ...
%
%   [Iseg,sep] = OTSU(I,n) returns the value (sep) of the separability
%   criterion within the range [0 1]. Zero is obtained only with images
%   having less than n gray level, whereas one (optimal value) is obtained
%   only with n-valued images.
%
%   Notes:
%   -----
%   It should be noticed that the thresholds generally become less credible
%   as the number of classes (n) to be separated increases (see Otsu's
%   paper for more details).
%
%   If n=2 or 3, the separability criterion is directly maximized by simply
%   using the MAX function. For n values >= 4, a minimization method is
%   used by means of the FMINSEARCH function.
%
%   The OTSU function works with I of any size. An RGB image I will thus be
%   considered as a gray-level 3D-array!
%   
%   Example:
%   -------
%   load clown
%   subplot(221)
%   imshow(X/max(X(:)))
%   title('Original','FontWeight','bold')
%   for n = 2:4
%     Iseg = otsu(X,n);
%     subplot(2,2,n), colormap(gray)
%     imshow(Iseg)
%     title(['n = ' int2str(n)],'FontWeight','bold')
%   end
%
%   -- Damien Garcia -- 2007/08

%% Checking n (number of classes)
if nargin==1
    n = 2;
elseif n==1;
    Iseg = NaN(size(I),'single');
    sep = 0;
    return
elseif n~=abs(round(n)) || n==0
    error('n must be a strictly positive integer!')
end

%% Probability distribution
I = single(I); % the HIST function only accepts single or double
unI = sort(unique(I));
nbins = min(length(unI),256);
if nbins==n
    Iseg = single(ones(size(I)));
    graycol = linspace(0,1,n);
    for i = 1:n-1
        Iseg(I==unI(i)) = graycol(i);
    end
    sep = 1;
    return
elseif nbins<n
    Iseg = NaN(size(I),'single');
    sep = 0;
    return
elseif nbins<256
    [histo,pixval] = hist(I(:),unI);
else
    [histo,pixval] = hist(I(:),256);
end
P = histo/sum(histo);
clear unI

%% Zeroth- and first-order cumulative moments
w = cumsum(P);
mu = cumsum((1:nbins).*P);

%% Maximal sigmaB^2 and Segmented image
if n==2
    sigma2B =...
        (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
    
    [maxsig,k] = max(sigma2B);
        
    % segmented image
    Iseg = single(ones(size(I)));
    Iseg(I<=pixval(k+1)) = 0;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
    thr=pixval(k+1);
    
elseif n==3
    w0 = w;
    w2 = fliplr(cumsum(fliplr(P)));
    [w0,w2] = ndgrid(w0,w2);
   
    mu0 = mu./w;
    mu2 = fliplr(cumsum(fliplr((1:nbins).*P))./cumsum(fliplr(P)));
    [mu0,mu2] = ndgrid(mu0,mu2);
    
    w1 = 1-w0-w2;
    w1(w1<=0) = NaN;

    sigma2B =...
        w0.*(mu0-mu(end)).^2 + w2.*(mu2-mu(end)).^2 +...
        (w0.*(mu0-mu(end)) + w2.*(mu2-mu(end))).^2./w1;
    sigma2B(isnan(sigma2B)) = 0; % zeroing if k1>= k2;
    
    [maxsig,k] = max(sigma2B(:));
    [k1,k2] = ind2sub([nbins nbins],k);
    
    % segmented image
    Iseg = single(ones(size(I)));
    Iseg(I<=pixval(k1)) = 0;
    Iseg(I>pixval(k1) & I<=pixval(k2)) = 0.5;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
    thr = [pixval(k1); pixval(k2)];

else
    % Threshold positions are adjusted using a horizontal mass-spring
    % system. Threshold positions can thus be adjusted by modulating the
    % force exerted on each mass.
    opt = [];
    K0 = 0.1;
    [F,y] = fminsearch(@sig_func,ones(n-1,1),opt,[n nbins K0 P]);
    
    K = K0*(diag(-2*ones(n-1,1))+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1));
    x0 = zeros(n-1,1); x0(n-1) = -K0;
    k = K\(x0+F-1);
    k = round(k*(nbins-1)+1);
    k = [0;k;nbins];
    
    % segmented image
    Iseg = single(ones(size(I)));
    graycol = linspace(0,1,n);
    for i = 1:n-1
        Iseg(I>=pixval(k(i)+1) & I<=pixval(k(i+1))) = graycol(i);
    end

    % separability criterion
    sep = 1-y; 
    
    thr=0;
    
end


%% Function to be minimized if n>=4
function y = sig_func(F,par)

n = par(1);
nbins = par(2);
K0 = par(3);
P = par(4:end);

K = K0*(diag(-2*ones(n-1,1))+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1));
x0 = zeros(n-1,1); x0(n-1) = -K0;
k = K\(x0+F-1);
k = round(k*(nbins-1)+1);

if any(k<1 | k>nbins) || any(diff(k)<1)
    y = 1;
    return
end

muT = sum((1:nbins).*P);
sigma2T = sum(((1:nbins)-muT).^2.*P);

k = [0;k;nbins];
sigma2B = 0;
for i = 1:n
    wi = sum(P(k(i)+1:k(i+1)));
    if wi==0
        y = 1;
        return
    end
    mui = sum((k(i)+1:k(i+1)).*P(k(i)+1:k(i+1)))/wi;
    sigma2B = sigma2B + wi*(mui-muT)^2;
end

y = 1-sigma2B/sigma2T; % within the range [0 1]

