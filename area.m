classdef area
    %AREAS Summary of this class goes here
    %   Detailed explanation goes here
    

    properties
        inertia
        damping
        bias_factor
        machines
        tg_con
        tg_sig
        mac_base
        bus
        to_bus
        mac_bus;
        lmod
        mac_nr
        res_bus
        res_nr
        res
        A
        B
        C
        C_mech
        C_res
        W_res
        W
        freq_feedback
        mec_index
        u_index
        res_index
        disturbance_res_index
    end

    methods
        function obj = area(bus)
            obj.inertia = 0;
            % 
            obj.damping = 0;
            % 
            obj.machines = 0;

            obj.bias_factor = 0;
            % 
            % 
            obj.tg_con = [];
            
            
            obj.tg_sig = [];
            
            obj.mac_base = [];

            obj.to_bus = [];

            obj.mac_bus = [];

            obj.bus = bus;
            obj.mac_nr = [];

            obj.res_nr = [];
            obj.res_bus = [];
            obj.res = 0;
            obj.freq_feedback = [];



            mec_index = [];
            u_index = [];
            res_index = [];
            disturbance_res_index = [];
            
        end
    end
end

