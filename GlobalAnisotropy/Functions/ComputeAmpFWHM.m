function [ currAmplitude, FWHM_degree ] = ComputeAmpFWHM( normPDF, OffSet_thres, rads )
%COMPUTEAMPFWHM Computes Amp and FWHM from Norm PDF of angles
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3

% Description by Zuhui: See SupMaterial of CTCF NatChemBio paper for more detail
% Briefly, 72 degress is considered to be background of isotropy and
% substracted from the angle PDF and to see the anisotropy effect.

% find the range in rads for the offset:
[~,idx] = min(abs(rads-OffSet_thres));
currOffSet = mean(normPDF(1:idx));
% calculate PDF without off-set:
normPDF_OffSet = normPDF - currOffSet;
currAmplitude = sum(normPDF_OffSet);

% make a fine-interpolated version for FWHM calculations:
fine_rads = 0:0.01:2*pi;
rads_for_interp1 = 0.5*(rads(2)-rads(1))+ rads(1:end-1);
interp1PDF = interp1(rads_for_interp1, normPDF_OffSet, fine_rads);

% find the FWHM: basically find the x-val closest to max/2:
HalfMax = max(interp1PDF)/2;
% function is symmetric, so consider only the first half
Half_interp1PDF = interp1PDF(1:round(length(interp1PDF)/2));
Half_fine_rads = fine_rads(1:round(length(interp1PDF)/2));
[~, FWHM_idx] = min(abs(Half_interp1PDF-HalfMax));
%FWHM_val = Half_interp1PDF(FWHM_idx);
FWHM_rad = Half_fine_rads(FWHM_idx);
FWHM_degree = 2*radtodeg(pi - FWHM_rad);

end

