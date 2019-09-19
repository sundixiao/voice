function [ data ] = Fn_x2k( image,dims,dont_shift )
%% ==================================================================
%FN_X2K do fft on selected dimensions
% ===================================================================
%   Author: Yibo Zhao @ UIUC
%   Created: 2018-07-23
%
%   [INPUTS]
%   ---- Required ----
%   image                   signal in time/spatial domain
%   dims                    dimensions along which to do the fft
%
%   ---- Optional ----
%   dont_shift              if true, don't do any fft shifts [false]
%
%   [OUTPUTS]
%   data                    signal in freq/k domain
%
%   Change log:
%       Created by  Yibo Zhao @ UIUC, 2018/07/23
%
%   Note:
%       F3_x2k(image) is equivalent to Fn_x2k(image,[1,2,3])
%       F_x2k(image)  is equivalent to Fn_x2k(image,[1,2])
%       F1_x2k(image) is equivalent to Fn_x2k(image,1)
%
%   See also F1_X2K, F_X2K, F3_X2K, F1_K2X, F_K2X, F3_K2X
%--------------------------------------------------------------------------

    if ~exist('dont_shift','var')||isempty(dont_shift)
        dont_shift = false;
    end
    scale_fctr = 1/sqrt(numel_dims(image,dims));
    if dont_shift
        for ind_dim = 1:length(dims)
            image = fft(image,[],dims(ind_dim));
        end
    else
        for ind_dim = 1:length(dims)
            image = ifftshift(fft(fftshift(image,dims(ind_dim)),[],dims(ind_dim)),dims(ind_dim));
        end
    end
    data = image*scale_fctr;
    return;

end

