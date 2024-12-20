function [mpc,n_ren] = get_g(data_file,res_bus,flag_ren)
%get_g   Get the mpc MATPOWER case with added RES with PST machine parameters.
%
%   Inputs:                          
%       
%       data_file - Name of the MATPOWER testcase. 
%       res_bus - Vector that specifies the buses where RESs are connected.
%       flag_ren - A boolean flag that says if there are RESs in the system.
%
%
%   Outputs:
%
%       mpc - Updated MATPOWER testcase. 
%       n_ren - Number of RESs in the system.


    % 'case118'
    if length(strsplit(data_file,'.mat')) == 2

        mpc = load(fullfile(pwd,'data\Synthethic Grids\',data_file));
        mpc = mpc.mpc;
    else
        mpc = loadcase(data_file);
    end

   % THIS COMMENT WAS TAKEN FROM PST documentation:
   % https://sites.ecse.rpi.edu/~chowj/PSTMan.pdf  as well as the
   % parameters
   % tg_con matrix format 
   %column        data         unit 
   %  1    turbine model number (=1)    
   %  2    machine number   
   %  3    speed set point   wf        pu 
   %  4    steady state gain 1/R       pu 
   %  5    maximum power order  Tmax   pu on generator base 
   %  6    servo time constant   Ts    sec 
   %  7    governor time constant  Tc  sec 
   %  8    transient gain time constant T3 sec 
   %  9    HP section time constant   T4   sec 
   % 10    reheater time constant    T5    sec     

    mpc.tg_con = [... 
        1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  2  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  3  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  4  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  5  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  6  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  7  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  8  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  9  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  10  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  11  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  12  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  13  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  14  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  15  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
        1  16  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  17  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  18  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  19  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  20  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  21  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  22  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  23  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  24  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  25  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  26  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  27  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  28  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  29  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  30  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  31  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  32  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  33  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  34  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  35  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  36  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  37  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  38  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  39  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  40  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  41  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  42  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  43  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  44  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  45  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  46  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  47  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  48  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  49  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  50  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  51  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  52  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  53  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
        1  54  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;];

    % THIS COMMENT WAS TAKEN FROM PST documentation:
    % https://sites.ecse.rpi.edu/~chowj/PSTMan.pdf  as well as the
    % parameters
    % Machine data format 
    %       1. machine number, 
    %       2. bus number, 
    %       3. base mva, 
    %       4. leakage reactance x_l(pu), 
    %       5. resistance r_a(pu), 
    %       6. d-axis sychronous reactance x_d(pu), 
    %       7. d-axis transient reactance x'_d(pu), 
    %       8. d-axis subtransient reactance x"_d(pu), 
    %       9. d-axis open-circuit time constant T'_do(sec), 
    %      10. d-axis open-circuit subtransient time constant  T"_do(sec),             
    %      11. q-axis sychronous reactance x_q(pu), 
    %      12. q-axis transient reactance x'_q(pu), 
    %      13. q-axis subtransient reactance x"_q(pu), 
    %      14. q-axis open-circuit time constant T'_qo(sec), 
    %      15. q-axis open circuit subtransient time constant T"_qo(sec),
    %      16. inertia constant H(sec), 
    %      17. damping coefficient d_o(pu), 
    %      18. dampling coefficient d_1(pu), 
    %      19. bus number 
    %                 

     mpc.mac_con =  [...
          1  53     300  0.003  	0  	0.969   0.248     0.147  	12.6     0.045  ...
      								    0.600   0.250     0      	0.035    0  ...
                                        3.4     0   	  0  		53  	  0.0654  0.5743;% hydro unit
          2  54     800  0.035  	0  	1.8     0.42529   0.30508   6.56  0.05 ...
         							    1.7207  0.3661    0.30508   1.5   0.035  ...
                                        4.9494  0   	  0  		54  	  0.0654  0.5743;
          3  55     800	 0.0304  0  	1.8  	0.38309   0.32465   5.7   0.05 ...
      								    1.7098  0.36072   0.32465   1.5   0.035  ...
      								    4.9623  0   	  0  		55  	  0.0654  0.5743;
   	      4  56  	800  0.0295  0  	1.8  	0.29954   0.24046   5.69  0.05 ...
      								    1.7725  0.27481   0.24046   1.5   0.035  ...
      								    4.1629  0   	  0  		56  	  0.0654  0.5743;
   	      5  57  	700   0.027  	0  	1.8     0.36  	  0.27273   5.4   0.05 ...
      								    1.6909  0.32727   0.27273   0.44  0.035  ...
      								    4.7667  0  	      0  		57  	  0.0654  0.5743;
   	      6  58  	900  0.0224  	0  	1.8  	0.35433   0.28346   7.3   0.05 ...
      								    1.7079  0.3189    0.28346   0.4   0.035  ...
      								    4.9107  0   	  0  		58 	  0.0654  0.5743;
   	     7  59  	800  0.0322  	0  	1.8  	0.29898   0.24407   5.66  0.05 ...
      								    1.7817  0.27458   0.24407   1.5   0.035  ...
      								    4.3267  0   	  0  		59  	  0.0654  0.5743;
   	     8  60  	800  0.028      0  	1.8  	0.35379   0.27931   6.7   0.05 ...
      								    1.7379  0.31034   0.27931   0.41  0.035   ...
      								    3.915  	0   	  0  		60  	  0.0654  0.5743;
   	    9   61     1000  0.0298  	0  	1.8  	0.48718   0.38462   4.79  0.05 ...
      								    1.7521  0.42735   0.38462   1.96  0.035  ...
      								    4.0365  0   	  0  		61  	  0.0654  0.5743;
   	    10  62     1200  0.0199  	0  	1.8  	0.48675   0.42604   9.37  0.05 ...
      								    1.2249  0.47929   0.42604   1.5   0.035  ...
      								    2.9106  0   	  0  		62  	  0.0654  0.5743;
   	    11  63     1600  0.0103  	0  	1.8  	0.25312   0.16875   4.1   0.05 ...
      								    1.7297  0.21094   0.16875   1.5   0.035  ...
      								    2.0053  0   	  0  		63  	  0.0654  0.5743;
   	    12  64     1900  0.022  	0  	1.8  	0.55248   0.44554   7.4   0.05 ...
      								    1.6931  0.49901   0.44554   1.5  	0.035  ...
      								    5.1791  0   	  0  		64  	  0.0654  0.5743;
   	    13  65    12000  0.003  	0  	1.8  	0.33446   0.24324   5.9   0.05 ...
      								    1.7392  0.30405   0.24324   1.5     0.035  ...
      								    4.0782  4.0782    0  		  65  	0.0654  0.5743;
   	    14  66    10000  0.0017  	0  	1.8    	0.285     0.23       4.1    0.05   ...
      								    1.73    0.25      0.23       1.5    0.035       ...
      								    3  		3   	  0  		  66  	0.0654  0.5743;
   	    15  67    10000  0.0017  	0  	1.8    	0.285     0.23       4.1    0.05   ...
      								    1.73    0.25      0.23   	 1.5    0.035       ...
      								    3  		3   	  0  		 67  	0.0654  0.5743;
   	    16  68    11000  0.0041  	0  	1.8  	0.35899   0.27809    7.8    0.05 ...
      								    1.6888  0.30337   0.27809    1.5    0.035    ...
                                        4.45  	4.45      0  		 68  	0.0654  0.5743;
    ];

    % Extend these parameters to the number of machines of MATPOWER testcase and assign it to the right buses 
    n_machines = size(mpc.gen,1);

    if n_machines > size(mpc.tg_con,1)
        mpc.tg_con = repmat(mpc.tg_con ,floor(n_machines/size(mpc.tg_con,1)),1);
        machines_left = mod(n_machines,size(mpc.tg_con,1));
        mpc.tg_con(end:end+machines_left,:) = mpc.tg_con(end-machines_left:end,:);
    else
        mpc.tg_con(n_machines+1:end,:) = [];
    end

    if n_machines > size(mpc.mac_con,1)
        mpc.mac_con = repmat(mpc.mac_con ,floor(n_machines/size(mpc.mac_con,1)),1);
        machines_left = mod(n_machines,size(mpc.mac_con,1));
        mpc.mac_con(end:end+machines_left,:) = mpc.mac_con(end-machines_left:end,:);
    else
        mpc.mac_con(n_machines+1:end,:) = [];
    end
    mpc.mac_con(:,1) = 1:n_machines;
    mpc.mac_con(:,2) = mpc.gen(:,1);
    

    % Assign RESs buses 
    
    idx = [];
    mpc.isgen = true(size(mpc.gen,1),1);
    mpc.isolar = [];
    mpc.isolar_mask = false(size(mpc.gen,1),1);



    n_ren = size(res_bus,1);

    if (n_ren ~= 0 || (size(mpc.bus,1) == 118 && flag_ren))
        if size(mpc.bus,1) == 118
            load('data\res_profile_118.mat');
            res_bus = res.bus;
            n_ren = size(res.bus,1);
        end
        gen_cost_to_add = repmat(mpc.gencost(1,:),size(res_bus,1),1);
        
        gen_to_add = repmat(mpc.gen(1,:),size(res_bus,1),1);
        
        gen_to_add(:,1) = res_bus;

        %Initial Active and Reactive Power
        gen_to_add(:,2) = 0 ;
        gen_to_add(:,3) = 0 ;

        gen_to_add(:,6) = 1 ;
    
        %Active power limits
        gen_to_add(:,9) = 0 ;
        gen_to_add(:,10) = 0;
        
        %Reactive power limits
        gen_to_add(:,4) = 0;
        gen_to_add(:,5) = 0;
        
        [mpc,idx] = addgen2mpc(mpc,gen_to_add, gen_cost_to_add, 'solar');
        
        mask = false(length(mpc.gen(:,1)),1);
        mask(mpc.isolar) = 1;
        mpc.isolar_mask = mask;
        mpc.isgen = ~mpc.isolar_mask;
        %n_res = size(data.bus,1);

    end



end

