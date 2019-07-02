function [H_V1] = HV(XH_magfilt,XV_magfilt)

H_V1 = XH_magfilt./XV_magfilt;
end