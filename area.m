classdef area
    %AREAS Summary of this class goes here
    %   Detailed explanation goes here
    

    properties
        inertia
        damping
        machines
        tg_con
        tg_sig
        mac_base
        bus
        to_bus
        lmod
    end

    methods
        function obj = area(bus)
            obj.inertia = 0;
            % 
            obj.damping = 0;
            % 
            obj.machines = 0;
            % 
            % 
            obj.tg_con = [];
            
            
            obj.tg_sig = [];
            
            obj.mac_base = [];

            obj.to_bus = [];

            obj.bus = bus;
        end
    end
end

