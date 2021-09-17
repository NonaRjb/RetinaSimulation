classdef Horizontal_2 < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        % constants Horizontal
		gA = 15.0;  % [nS]
        gKv = 4.5;  % [nS]
        gl = 0.5;   % [nS]
        gKa = 4.5;  % [nS]
        gCa = 9.0;  % [nS]
        gNa = 2.4;  % [nS]
        Ek = -80;   % [mV]
        El = -80;   % [mV]
        ENa = 55;   % [mV]
        C = 0.106;  % [nF]
        Cao = 2.0;  % [mM]
        Cad = 30;   % [uM]
        % variables
        Y
    end
    
    methods
        function obj = Horizontal_2(Y0, buffer_size, dt)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==1
                buffer_size = 10000;
                dt = 1e-06;
            elseif nargin == 2
                dt = 1e-06;
            end
            
            obj.t = 1;
            obj.dt = dt;
            obj.buffer_size = buffer_size;
            obj.flag = 0;
            
            obj.t_vec = zeros(buffer_size, 1); 
            obj.Y = zeros(8, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
            % obj.V0 = Y0(1);
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, curr_t, c]  = solve(obj,Vpre)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
			
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                k = obj.t+1;
            else
                k = obj.t;
            end
			% explicit functions
            %%% glu %%%
            glu = 100;
			%%% A %%%
            Vk = obj.Y(1, k) - obj.Ek;
            amA = 2400/(1+exp(-(obj.Y(1, k)-50)/28));
            bmA = 80*exp(-obj.Y(1, k)/36);
            ahA = exp(-obj.Y(1, k)/60);
            bhA = 20/(exp(-(obj.Y(1, k)+40)/5)+1);
            iA = obj.gA*obj.Y(7, k)^3*obj.Y(8, k)*Vk;
            %%% Kv %%%
            amKv = 0.40*(65-obj.Y(1, k))/(exp((65-obj.Y(1, k))/50)-1);
            bmKv = 4.8*exp((45-obj.Y(1, k))/85);
            ahKv = 1500/(exp((obj.Y(1, k)+92)/7)+1);
            bhKv = 80/(exp((obj.Y(1, k)+100)/15)+1)+0.02;
            iKv = obj.gKv*obj.Y(5, k)^4*obj.Y(6, k)*Vk;
            %%% Ka %%%
            mKa = 1/(1+exp((obj.Y(1, k)+60)/12));
            iKa = obj.gKa*mKa^5*Vk;
            %%% Ca %%%            
            amCa = 240*(68-obj.Y(1, k))/(exp((68-obj.Y(1, k))/21)-1);
            bmCa = 800/(exp((55+obj.Y(1, k))/55)+1);
            ECa = 12.9*log(obj.Cao/obj.Cad);
            iCa = obj.gCa*obj.Y(4, k)^4*(obj.Y(1, k)-ECa);
            %%% Na %%%
            amNa = 200*(38-obj.Y(1, k))/(exp((38-obj.Y(1, k))/25)-1);
            bmNa = 2000*exp(-(55+obj.Y(1, k))/18);
            ahNa = 1000*exp(-(80+obj.Y(1, k))/8);
            bhNa = 800/(exp((80-obj.Y(1, k))/75)+1);
            iNa = obj.gNa*obj.Y(2, k)^3*obj.Y(3, k)*(obj.Y(1, k)-obj.ENa);
            %%% glu %%%
%             iGlu = 800*(exp((obj.Y(1, k)-obj.Eglu)/125)-1)/...
%                 (exp(-(obj.Y(1, k)+obj.Eglu)/125)+1)*glu^3/(glu^3+obj.Kglu^3);
            %%% leak %%%
            iL = obj.gl*(obj.Y(1, k)-obj.El);
			
            consts = [obj.C];
            
            vars = [iCa;iA;iKa;iKv;iL;iNa;amA;bmA;ahA;bhA;amKv;bmKv;ahKv;...
                bhKv;amCa;bmCa;amNa;bmNa;ahNa;bhNa];
            
            D = f(obj.Y(:,k), vars, consts);
            
            %%%% specify dt %%%%
%             
%             if abs(max(D, [], 'all')) > 300
%                 obj.dt = 1e-06;
%             elseif abs(max(D, [], 'all')) <= 300 && abs(max(D, [], 'all')) > 200
%                 obj.dt = 2*1e-06;
%             elseif abs(max(D, [], 'all')) <= 200 && abs(max(D, [], 'all')) > 100
%                 obj.dt = 5*1e-06;
%             elseif abs(max(D, [], 'all')) <= 100 && abs(max(D, [], 'all')) > 10
%                 obj.dt = 1e-05;
%             elseif abs(max(D, [], 'all')) <= 100 && abs(max(D, [], 'all')) > 10
%                 obj.dt = 5*1e-05;
%             else
%                 obj.dt = 1e-04;
%             end
            
%             maxD = abs(max(D, [], 'all'));
%             if maxD > 1e-03
%                 obj.dt = 1/maxD*1e-03;
%             else
%                 obj.dt = 0.1;
%             end
            %objdt = obj.dt;
            
            %%%%% runge-kutta 4
%             k1 = obj.dt * D;
%             k2 = obj.dt * f(obj.Y(:, k)+k1/2, vars, consts);
%             k3 = obj.dt * f(obj.Y(:, k)+k2/2, vars, consts);
%             k4 = obj.dt * f(obj.Y(:, k)+k3, vars, consts);
%             k_tot = 1/6*(k1+2*k2+2*k3+k4);
            %%%%% forward euler
            k_tot = obj.dt*D;
            
            
            % values of variables at the next time step
            
            % we have already reached maximum buffer size so we should discard
            % primary values of the variables
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                obj.Y = [obj.Y(:, 2:end), obj.Y(:, k)+k_tot];
                obj.t_vec = [obj.t_vec(2:end); obj.t_vec(end)+obj.dt];
                y = obj.Y(:, end);
                curr_t = obj.t_vec(end);
            % we have just reached maximum buffer size so e should change
            % flag value to 1
            elseif (obj.t+1 == obj.buffer_size) && (obj.flag == 0)
                obj.Y(:, k+1) = obj.Y(:, k)+k_tot;
                obj.t_vec(k+1) = obj.t_vec(k)+obj.dt;
                obj.flag = 1;
                y = obj.Y(:, k+1);
                curr_t = obj.t_vec(k+1);
            % we have not yet reached the buffer size
            else
                obj.Y(:, k+1) = obj.Y(:, k)+k_tot;
                obj.t_vec(k+1) = obj.t_vec(k)+obj.dt;
                y = obj.Y(:, k+1);
                curr_t = obj.t_vec(k+1);
            end
                        
            c = [iCa; iNa; iKv; iA; iKa; iL];
           
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function tVector = get_tvec(obj)
            tVector = obj.t_vec;
        end
    end
end


%%%% F %%%%
function D = f(Y, vars, consts)

% consts
C = consts(1);

% vars
iCa = vars(1);
iA = vars(2);
iKa = vars(3);
iKv = vars(4);
iL = vars(5);
iNa = vars(6);
amA = vars(7);
bmA = vars(8);
ahA = vars(9);
bhA = vars(10);
amKv = vars(11);
bmKv = vars(12);
ahKv = vars(13);
bhKv = vars(14);
amCa = vars(14);
bmCa = vars(15);
amNa = vars(16);
bmNa = vars(17);
ahNa = vars(18);
bhNa = vars(19);


D = zeros(8, 1);
D(1) = -(iCa+iA+iKa+iKv+iL+iNa-10)/C;
D(2) = amNa*(1-Y(2))-bmNa*Y(2);
D(3) = ahNa*(1-Y(3))-bhNa*Y(3);
D(4) = amCa*(1-Y(4))-bmCa*Y(4);
D(5) = amKv*(1-Y(5))-bmKv*Y(5);
D(6) = ahKv*(1-Y(6))-bhKv*Y(6);
D(7) = amA*(1-Y(7))-bmA*Y(7);
D(8) = ahA*(1-Y(8))-bhA*Y(8);

% if D(8) < 0
%     D(8)
% end
end
