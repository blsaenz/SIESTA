% coversion of land_mask
% ======================================================================
% SIESTA extecting text file with projected order, inverted from lat/lon
% below.  Don't ask why, just do as you are told.  1 = water, 0 = not water

fid = fopen('Nl.LOCImask_land50_coast0km.721x721.bin','rb');
mask_n = fread(fid,[721,721],'uint8');
mask_n(mask_n ~= 255)=0;
mask_n(mask_n == 255)=1;
mask_ni = mask_n';
save EASE_25km_land_mask_721x721_NORTH.txt mask_ni -ascii

% Then using text editor strip all the x0000000 so that file has one space and one digit,
% like this: ' 0 0 0 0 1 1 0 0 0' etc in a 721x721 text matrix.

% conversion of lat/lon grids - FORTRAN/SIESTA expecting same order as delivered
% ======================================================================
% not-in-grid lat/lons are set to -200 

fid = fopen('NLLATLSB','rb');
lat_n = fread(fid, [721, 721], 'int32')/100000;
lat_n(lat_n==1.431655765000000e+04)=-200;
fid2=fopen('EASE_25km_lat_invert_NORTH.double.le','wb');
fwrite(fid2,lat_n,'double');
fclose(fid)
fclose(fid2)