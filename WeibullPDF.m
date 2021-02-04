function pdf = WeibullPDF(Vw, vav, k, Type, PLOT)
Vw = Vw(:);
% A = vav/(sqrt(pi)/k);   % weibull scale factor ( from Burton Chap.2)
A = vav/exp(gammaln(1+1/k)); % Correct for variable shape factor k

[Vw_sort, ~, ind_sort] = unique(Vw);

BinSize = Vw_sort(2:end)-Vw_sort(1:end-1);
BinSize = [BinSize(1); (BinSize(1:end-1)+BinSize(2:end))/2; BinSize(end)];
BinSize = BinSize(ind_sort);

if strcmp(Type,'pdf')
    %% Point approach
    pdf = k/A*(Vw/A).^(k-1).*exp(-(Vw/A).^k); % weibull probability density function
    % pdf = (pi/(2*vav^2)*Vw).*exp(-pi/4*(Vw./vav).^2); % Rayleigh pdf - same as weibull but with k=2

elseif strcmp(Type,'cdf')
    %% Binned approach
    Vw_edges = round(Vw)+BinSize.*[-1,1]/2;
    pdf = (1 - exp(-(Vw_edges(:,2)/A).^k)) - (1 - exp(-(Vw_edges(:,1)/A).^k)); % Weibull cumulative density function
end

if PLOT
    Vfine = 0:0.01:40;
%     switch Type
%         case 'pdf'
%             pdf_fine = k/A*(Vfine/A).^(k-1).*exp(-(Vfine/A).^k);
%         case 'cdf'
%             Vw_edges = round(Vfine')+BinSize.*[-1,1]/2;
%             pdf_fine = (1 - exp(-(Vw_edges(:,2)/A).^k)) - (1 - exp(-(Vw_edges(:,1)/A).^k)); % Weibull cumulative density function
%     end
    pdf_fine = k/A*(Vfine/A).^(k-1).*exp(-(Vfine/A).^k);
    figure, hold on,
    plot(Vfine, pdf_fine,'-')
    plot(Vw,pdf,'*')
    legend('Full PDF','PDF at user-defined points')
end


end % end function