function [ image ] = Fn_k2x( data,dims,dont_shift )
%% ==================================================================
%FN_K2X do ifft on selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-07-23
%
%   [INPUTS]
%   ---- Required ----
%   data                    signal in freq/k domain
%   dims                    dimensions along which to do the fft
%
%   ---- Optional ----
%   dont_shift              if true, don't do any fft shifts [false]
%
%   [OUTPUTS]
%   image                   signal in time/spatial domain
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/07/23
%
%   Example:
%       F3_k2x(data) is equivalent to Fn_k2x(data,[1,2,3])
%       F_k2x(data)  is equivalent to Fn_k2x(data,[1,2])
%       F1_k2x(data) is equivalent to Fn_k2x(data,1)
%
%   See also F1_K2X, F_K2X, F3_K2X, F1_X2K, F_X2K, F3_X2K
%--------------------------------------------------------------------------

    if ~exist('dont_shift','var')||isempty(dont_shift)
        dont_shift = false;
    end
    scale_fctr = sqrt(numel_dims(data,dims));
    if dont_shift
        for ind_dim = 1:length(dims)
            data = ifft(data,[],dims(ind_dim));
        end
    else
        for ind_dim = 1:length(dims)
            data = ifftshift(ifft(fftshift(data,dims(ind_dim)),[],dims(ind_dim)),dims(ind_dim));
        end
    end
    image = data*scale_fctr;
    return;

end

