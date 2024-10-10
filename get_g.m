function [mpc,n_ren,idx] = get_g(data_file,res_bus)

    % 'case118'
    if length(strsplit(data_file,'.mat')) == 2

        mpc = load(fullfile(pwd,'data\Synthethic Grids\',data_file));
        mpc = mpc.mpc;
    else
        mpc = loadcase(data_file);
    end
    % 
    % g.mac.mac_con=mac_con;
    % %g.tg.tg_con=tg_con;
    % g.lmod.lmod_con=lmod_con;
    % 
    % g.line = line;
    % 
    % %g.bus = line;
    % g.bus.bus_int = bus(:,1);


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


     mpc.mac_con =  [
         1  10 991 0.046  0 0.295 0.081  0  5.70  0 ...
                            0.286 0.081  0  0.035 0 ...
                            30.3  30.3   0  36    0  0   1.00 1.00;
         2  21 991 0.029  0 0.254 0.048  0  5.70  0 ...
                            0.241 0.048  0  0.035 0 ...
                            34.8  34.80  0  21    0  0   1.00 1.00;
         3  22 991 0.0295 0 0.262 0.0436 0 5.69   0 ...
                            0.258 0.0436 0 0.035  0 ...
                            28.6  28.6   0  22    0  0   1.00 1.00;
         4  23 991 0.0333 0 0.58  0.0815 0 5.60   0 ...
                            0.56  0.0815 0  0.035 0 ...
                            7.34  7.34   0  23    0  0   0.55 0.55;
         5  23 991 0.0745 0 0.705 0.182  0 5.20   0 ...
                            0.680 0.182  0  0.035 0 ...
                            18.6  18.6   0  23    0  0   0.45 0.45;
         6  24 991 0.029  0 0.254 0.048  0 7.30   0 ...
                            0.241 0.048  0 0.035  0 ...
                            34.8  34.8   0  24    0  0   1.00 1.00;
         7  25 991 0.0322 0 0.295 0.048  0 5.70   0 ...
                            0.292 0.048  0 0.035  0 ...
                            26.4  26.4   0  25    0  0   1.00 1.00;
         8  26 991 0.313 0 2.10 0.58  0 4.79   0 ...
                            2.05 0.58  0 0.035 0 ...
                            3.42  3.42   0  26    0  0   1.00 1.00;
         9  27 991 0.030  0 0.290 0.055  0 6.70   0 ...
                            0.280 0.055  0  0.035 0 ...
                            24.30 24.30  0  27    0  0   1.00 1.00;
        10  42 991 0.0445 0 0.3162 0.0595 0 5.00  0 ...
                            0.3027 0.0595 0 0.035 0 ...
                            18.86  18.86  0 42    0  0   1.00 1.00;
        11  47 991 0.0782 0 0.298  0.110  0 5.00  0 ...
                            0.171  0.110  0 0.035 0 ...
                            15.17  15.17  0 47    0  0   1.00 1.00;
        12  48 991 0.0782 0 0.298  0.110  0 5.00  0 ...
                            0.171  0.110  0 0.035 0 ...
                            15.17  15.17  0 48    0  0   1.00 1.00;
        13  50 991 0.059  0 0.275  0.078  0 6.70  0 ...
                            0.261  0.078  0 0.035 0 ...
                            34.44  34.44  0 50    0  0   1.00 1.00;
        14  51 991 0.0413 0 0.3103 0.0625 0 5.50  0 ...
                            0.2941 0.0625 0 0.035 0 ...
                            21.46  21.46  0 51    0  0   1.00 1.00;
        15  53 991 0.0000 0 0.     0.070  0 0.00  0 ...
                            0      0      0  0    0 ...
                           37.0    37.0   0 53    0  0   1.00 1.00;
        16  54 991 0.050  0 0.1322 0.0566 0 7.00 0 ...
                            0.0900 0.0566 0  0.035 0 ...
                            51.24  51.24  0 54    0  0   0.50 0.50;
        17  54 991 0.050 0  0.1322 0.0566  0 7.00 0 ...
                            0.0900 0.0566  0 0.035 0 ...
                            51.24   51.24  0 54    0  0   0.50 0.50;
        18  55 991 0.0302 0 0.1016 0.0364 0 7.00  0 ...
                            0.0643 0.0364 0 0.035 0 ...
                            59.78  59.78  0 55    0  0   1.00 1.00;
        19  56 991 0.0515 0 0.3847 0.0715 0 5.50  0 ...
                            0.3697 0.0715 0 0.035 0 ...
                            21.940 21.94  0 56    0  0   1.00 1.00;
        20  57 991 0.0464 0 0.3373 0.0699 0 5.60  0 ...
                            0.3227 0.0699 0  0.035 0 ...
                            17.71  17.71  0 57    0  0   1.00 1.00;
        21  60 991 0.0523 0 0.3852 0.0722 0 5.50  0 ...
                            0.3705 0.0722 0  0.035 0 ...
                            21.74  21.74  0 60    0  0   1.00 1.00;
        22  61 991 0.0947 0 0.6384 0.1248 0 5.00  0 ...
                            0.6114 0.1248 0 0.035 0 ...
                            8.320  8.32   0 61    0  0   1.00 1.00;
        23  65 991 0.0000 0 0      0.1220 0 0  0 ...
                            0      0      0 0  0 ...
                            10.71  10.71  0 65    0  0   1.00 1.00;
        24  68 991 0.0000 0 0      0.2686 0 0.00 0 ...
                            0      0 0 0 0 ...
                            3.77   3.77   0 68    0  0   1.00 1.00;
        25  71 991 0.0000 0 0      0.1900 0 0.00 0 ...
                            0      0      0 0    0 ...
                            8.50   8.500  0 71    0  0   1.00 1.00;
        26  72 991 0.0000 0 0      0.1185 0 0.00 0 ...
                            0      0 0 0 0 ...
                            11.77  11.77  0 72    0  0   1.00 1.00;
        27  78 991 0.0000 0 0      0.0001 0 0.00 0 ...
                            0      0 0 0 0 ...
                            1000.0 1000.0 0 78    0  0   1.00 1.00;
        28  79 991 0.0280 0 0.1740 0.032  0 8.00  0  ...
                            0.1700 0.032  0 0.035 0 ...
                            48.0   48.0   0 79    0  0   1.00 1.00;
        29  80 991 0.0250 0 0.213  0.047  0 4.00  0 ...
                            0.193  0.047  0 0.035 0 ...
                            23.8   23.8   0 80    0  0   1.00 1.00;
        30  82 991 0.0500 0 0.3090 0.0620 0 5.90 0 ...
                            0.3020 0.0620 0 0.035 0 ...
                            19.600 19.6   0 82    0  0   1.00 1.00;
        31 101 991 0.0170 0 0.0750 0.0260 0 5.50 0 ...
                            0.0480 0.0260 0 0.035 0 ...
                            55.000 55.00  0 101   0  0   1.00 1.00;
        32  86 991 0.0100 0 0.100  0.020  0 6.00 0 ...
                            0.0840 0.020  0 0.035 0 ...
                            79.000 79.000 0 86    0  0   1.00 1.00;
        33  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        34  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        35  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;                   
        36  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        37  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        38  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        39  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        40  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        41  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        42  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        43  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        44  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        45  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        46  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        47  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
        48  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00; 
                            
        49  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;                    
                            
        50  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;                    
                            
        51  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;                    
                            
        52  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;                   
        53  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;
                            
        54  91 991 0.0000 0 0      0.0150 0 0.00 0 ...
                            0      0 0 0 0 ...
                            98.7   98.7   0 91    0  0   1.00 1.00;];

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
    
    idx = [];
    mpc.isgen = true(size(mpc.gen));
    mpc.isolar = [];



    n_ren = size(res_bus,1);

    if (n_ren ~= 0 || size(mpc.bus,1) == 118)
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

