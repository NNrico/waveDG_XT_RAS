function  [Solutions]  = Save_Solution(UH,femregion)


np = femregion.ndof_p;
ne = femregion.ndof_e;
na = femregion.ndof_a;

up_h  = [];
wp_h  = [];
phi_h = [];
ue_h  = [];
dot_up_h  = [];
dot_wp_h  = [];
dot_phi_h = [];
dot_ue_h  = [];


if (np > 0 && na > 0 && ne > 0)
    % poro - acoustic - elastic 
    
    up_h   = UH(1         : 2*np);
    wp_h   = UH(2*np+1    : 4*np);
    phi_h  = UH(4*np+1    : 4*np+na);
    ue_h   = UH(4*np+na+1 : 4*np+na+2*ne);
    
    ndof_vel = 4*np+na+2*ne;
    
    dot_up_h   = UH(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = UH(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_phi_h  = UH(ndof_vel+4*np+1     : ndof_vel+4*np+na);
    dot_ue_h   = UH(ndof_vel+4*np+na+1  : ndof_vel+4*np+na+2*ne);
    
              
    
elseif (np > 0 && na > 0 && ne == 0)
    % poro - acoustic
    
    up_h   = UH(1         : 2*np);
    wp_h   = UH(2*np+1    : 4*np);
    phi_h  = UH(4*np+1    : 4*np+na);
    
    ndof_vel = 4*np+na;
    
    dot_up_h   = UH(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = UH(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_phi_h  = UH(ndof_vel+4*np+1     : ndof_vel+4*np+na);
    
              
elseif (np > 0 && na == 0 && ne > 0)
    % poro  - elastic 
    
    up_h   = UH(1         : 2*np);
    wp_h   = UH(2*np+1    : 4*np);
    ue_h   = UH(4*np+1    : 4*np+2*ne);
    
    ndof_vel = 4*np+2*ne;
    
    dot_up_h   = UH(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = UH(ndof_vel+2*np+1     : ndof_vel+4*np);
    dot_ue_h   = UH(ndof_vel+4*np+1     : ndof_vel+4*np+2*ne);
    
              
              
elseif (np == 0 && na > 0 && ne > 0)
    % acoustic - elastic 
    
    phi_h  = UH(1    : na);
    ue_h   = UH(na+1 : na+2*ne);
    
    ndof_vel = na+2*ne;
    
    dot_phi_h  = UH(ndof_vel+1     : ndof_vel+na);
    dot_ue_h   = UH(ndof_vel+na+1  : ndof_vel+na+2*ne);
    
              
elseif (np > 0 && na == 0 && ne == 0)
    % poro
    
    up_h   = UH(1         : 2*np);
    wp_h   = UH(2*np+1    : 4*np);
   
    ndof_vel = 4*np;
    
    dot_up_h   = UH(ndof_vel + 1        : ndof_vel+2*np);
    dot_wp_h   = UH(ndof_vel+2*np+1     : ndof_vel+4*np);
        
 
elseif (np == 0 && na > 0 && ne == 0)
    % acoustic 
    
    phi_h  = UH(1    : na);
 
    ndof_vel = na;
    
    dot_phi_h  = UH(ndof_vel+1     : ndof_vel+na);
    

elseif (np == 0 && na == 0 && ne > 0)
    % elastic 
    
    ue_h   = UH(1 : 2*ne);
    
    ndof_vel = 2*ne;
    
    dot_ue_h   = UH(ndof_vel+1  : ndof_vel+2*ne);
    
              
end


Solutions.up_h      = up_h;
Solutions.wp_h      = wp_h;
Solutions.phi_h     = phi_h;
Solutions.ue_h      = ue_h;
Solutions.dot_up_h  = dot_up_h;
Solutions.dot_wp_h  = dot_wp_h;
Solutions.dot_phi_h = dot_phi_h;
Solutions.dot_ue_h  = dot_ue_h;

