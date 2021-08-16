% initialize the tests

n = 16;

% Markers
msize = 5;

marker_cosm     =  'v';
marker_cosm_adv  =  's';
marker_cosm_exp  =  'o';
marker_cosm_mp  =  '^';
marker_cosm_mp_s   =  'p';
marker_cosm_tay  =  'x';
marker_cosm_pol  =  '+';

marker_fre_blk = 'v';
marker_fre_exp = 'o';
marker_fre_mp = '^';
marker_fre_mp_s = 'p';

color_cosm    = [0, 0, 1];
color_cosm_adv = [0.0 0.9 0.9];
color_cosm_exp = [1.0 0.5 0.0]; 
color_cosm_mp = [0.8 0.0 0.8];
color_cosm_mp_s  = [0.25, 0.25, 0.25];
color_cosm_tay = [0.0 0.5 0.0];
color_cosm_pol = [0.75, 0.75, 0]; 

color_fre_blk    = [0, 0, 1];
color_fre_exp = [1.0 0.5 0.0]; 
color_fre_mp = [0.8 0.0 0.8];
color_fre_mp_s  = [0.25, 0.25, 0.25];

color_cond    = [0.25, 0.25, 0.25];

ls_cond = '-';

ls_cosm = '-';
ls_cosm_adv = '-.';
ls_cosm_exp = '-.';
ls_cosm_mp = '--';
ls_cosm_mp_s = '--';
ls_cosm_tay = '-.';
ls_cosm_pol = '-.';

ls_fre_blk = '-';
ls_fre_exp = '-.';
ls_fre_mp = '--';
ls_fre_mp_s = '--';

lw = 1.0; % linewidth
lw_cond = 1.5; %linewidth of condu