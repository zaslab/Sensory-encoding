function [c_pvalues, c_alpha, h, extra] = fwer_bonf(pvalues, alpha, plotting)
    % Error detection
    if nargin < 3, plotting = false; end
    if nargin < 2, error('Not enough parameters.'); end
    if ~isnumeric(pvalues) && ~isnumeric(alpha)
        error('Parameters pvalues and alpha must be numeric.');
    end
    pvalues = pvalues(:);
    if length(pvalues) < 2, error('Not enough tests to perform the correction.'); end
    
    % Parameters
    m = length(pvalues);    % No. tests
    
    % Corrected pvalues
    c_pvalues = min(pvalues.*m,1);
    
    % Corrected significance levels
    c_alpha = (alpha/m).*ones(length(pvalues),1);
    
    % Rejected H0
    h = pvalues(:) < c_alpha(:);
    
    % Extra information
    extra.s_pvalues = sort(pvalues,'ascend');
    extra.s_c_pvalues = sort(c_pvalues,'ascend');
    extra.alpha = alpha;
    extra.pvalues = pvalues;
    
    % Plotting
    if plotting
        figure;
        subplot(2,2,1:2);
        plot(extra.s_pvalues, extra.s_c_pvalues, 'b', 'linewidth',2);
        ylabel('Adj. p-values'); xlabel('p-values');
        title('Bonferroni');
        
        subplot(2,2,3);
        hist(pvalues); xlabel('p-values'); ylabel('Histogram');
        
        subplot(2,2,4);
        hist(c_pvalues); xlabel('Adj. p-values'); ylabel('Histogram');
    end
end