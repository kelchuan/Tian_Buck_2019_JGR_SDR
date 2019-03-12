%load FLAC data
function [particles_4]=SDR_load_data()
if exist('platform.s','file') == 0
  plat = 'native';
elseif exist('platform.s','file') == 2
  plat = 'ieee-le';
end

fid_t = fopen('time.0','r',plat);

load nxnz.0
nx=nxnz(1);
ny=nxnz(2);

t = fread(fid_t,'single');
%nt=length(t);
nt=length(t);

num_particles = 360;
%{
if exist('particles_02.0','file')>0
    fid_part = fopen('particles_02.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_2 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
%}
if exist('particles_04.0','file')>0
    fid_part = fopen('particles_04.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_4 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
