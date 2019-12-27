function [mu, ksmu , uwpks, uwpkf] = fwdbwdsmoother(Y, nsettg)

switch nargin
case 2
    isub = nsettg.sub;
    ilsp = nsettg.lsp;
    itr = nsettg.iter;
    ielectrode = nsettg.e;
    param.fc = nsettg.f;
    saveDirec = nsettg.saveDirec;
    flag_type = 1;
    otherwise 
    flag_type = 0;
end

    fc = Y.fc;
    fs = Y.dwnfs;
    dt = Y.hy;
    R = Y.R;
    factor = Y.cfactor;

    w0 = 2 * pi * fc / fs;  
    T = length(dt);
    N = size(Y,1);
    mu = zeros(1,T);
    ksmu = zeros(1,T);
    p  = zeros(1,T);
    % define if the variance of noise is constant or varying over time
    %[vr, vq] = estimate_rq(Y);
    vr = R;
    A = mean(dt);
    B = [A(1) (exp(1i * w0) * A(1:end-1))];
    q = factor*var(A - B);
    %q = 0.000017; %0.000005;
    N = 1;

    p(1) = var(dt(:,1));
    mu(1) = mean(dt(:,1));

        % forward pass
        for k = 2:T 
            % if we assume observations y are independent, r and q can be estimated
            % as the variance and mean of the data. Else, they are initialized as
            % constant.
           r = vr(k);
           %q = vq(k);
           mu(k) = exp(1i * w0) * mu(k-1);
           p(k) = p(k-1) + q;
           % compute Kalman Gain Factor
           KG(k) = p(k)/(p(k) + (r/N));
           % corrections based on the previous observation
           tilday = mean(dt(:,k));
           mu(k) = mu(k) + KG(k) * (tilday - mu(k));
           p(k) = (r/N) * KG(k);
        end


        % plot the phase before and after applying the KF
        %h(1) = figure;
        %subplot(2,1,1), plot(abs(mean(dt))),  % before applying KF
        %hold on, plot(abs(mu),'g'), xlim([1 length(abs(mu))]), grid on,  % after applying KF

        %subplot(2,1,2), plot(unwrap(angle(mean(dt))) - (2 * pi * fc / fs) * (1:length(dt)),'b'),
        %hold on, plot(unwrap(angle(mu)) - (2 * pi * fc / fs) * (1:length(mu)),'g'),
        %xlim([1 length(dt)]); title('IP w KF'), grid on;
        
        
        uwpkf = unwrap(angle(mu)) - (2 * pi * fc / fs) * (1:length(mu));
        

        % backward pass
        v = zeros(1,T);
        v(T) = p(T);
        ri = exp(-1i * w0);
        ksmu(T) = mu(T);

        for k = (T-1):-1:1
            lambda = p(k)/(p(k) + q);
            ksmu(k)  = mu(k) + lambda * (ri * ksmu(k+1) - mu(k));
            v(k)   = p(k) + lambda^2 * (v(k+1) - p(k) - q);
        end

        % plot the phase after applying the KS        
        %subplot(2,1,1), plot(abs(ksmu),'k'), xlim([1 length(abs(dt))]),
        %subplot(2,1,2), plot(unwrap(angle(ksmu)) - (2 * pi * fc / fs) * (1:length(ksmu)),'k');
        
        uwpks = unwrap(angle(ksmu)) - (2 * pi * fc / fs) * (1:length(ksmu));
        
        if isequal(flag_type,1) 
                fname = ...
            strcat(saveDirec,'fsub_',num2str(isub),'_lsp',num2str(ilsp),...
            '_itr',num2str(itr),'_e',num2str(ielectrode),'_f',num2str(param.fc),'.fig');
                savefig(h,fname);
                close all
        end
        
        close all
        
end
