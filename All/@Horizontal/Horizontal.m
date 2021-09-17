classdef Horizontal < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        method
        % constants Horizontal
		gA = 2.2;   % [nS]
        gKv = 0.4;  % [nS]
        gl = 0.2;   % [nS]
        gKa = 24;   % [nS]
        gNa = 200;  % [nS]
        Ek = -96;   % [mV]
        El = -70;   % [mV]
        Eglu = -10; % [mV]
        ENa = 81;   % [mV]
        Kglu = 20;  % [uM]
        Kka = 1;    % [mM]
        nu = 6.5e-14;   % [m^3]
        F = 9.646e04;   % [C/mol]
        R = 8.31;       % [j/molK]
        T = 293;        % [K]
        Cao = 2.5;      % [mM]
        Btot = 20;      % [uM]
        KCa = 0.5;      % [uM]
        PCa = 2;        % [nm/s]
        Iex = 5.0;      % [pA]
        Ipump = 5.0;    % [pA]
        alpha = 0.15;
        fB1 = 4.0;      % [1/uMs]
        fB2 = 100.0;    % [1/uMs]
        bB1 = 2.0;      % [1/s]
        bB2 = 100;      % [1/s]
        B1t = 10.0;     % [uM]
        B2t = 20.0;     % [uM]
        kb = 1.5;       % [uM]
        kf = 0.1;       % [uM]
        bSB = 0.8;      % [1/uMs]
        fSB = 0.2;      % [1/uMs]
        gSA = 1.5;      % [nS]
        C_S = 0.0125;   % [nF]
        iAT = 6.48;     % [pA]
        tfb = 20;       % [sec^(-1)]
        gs = 42.0;      % [nS]
        gAT = 420;      % [nS]
        C_AT = 0.00625; % [nF]
        hAT = 0.17;     % [nS]
        V0 = 0;
        % variables
        Y
    end
    
    methods
        function obj = Horizontal(Y0, buffer_size, dt, method)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==1
                buffer_size = 10000;
                dt = 1e-06;
                method = 'euler';
            elseif nargin == 2
                dt = 1e-06;
                method = 'euler';
            elseif nargin == 3
                method = 'euler';
            end
            
            obj.t = 1;
            obj.dt = dt;
            obj.buffer_size = buffer_size;
            obj.flag = 0;
            obj.method = method;
            
            obj.t_vec = zeros(buffer_size, 1); 
            obj.Y = zeros(10, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
            obj.V0 = Y0(1);
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
            %%% obj.Y(11, k) %%%
            % glu = 1/(1+exp(-Vpre+obj.alpha*obj.Y(10,k)));
            glu = exp((Vpre+obj.alpha*obj.Y(10,k))/30);
%             if (Vpre-obj.alpha*obj.Y(10, k)-(-36.185963)) >= 0
%                 glu = tanh((Vpre-obj.alpha*obj.Y(10, k)-(-36.185963)));
%             else
%                 glu = 0;
%             end
%             if obj.t_vec(k) < 1
%                 glu = 0;
%             else
%                 glu = 100;
%             end
            %glu = tanh((Vpre-obj.alpha*obj.Y(10, k))/20);
			%%% A %%%
            Vk = obj.Y(1, k) - obj.Ek;
            amA = 80*(75-obj.Y(1, k))/(exp((75-obj.Y(1, k))/14)-1);
            bmA = 8*exp(-(obj.Y(1, k)-40)/60);
            ahA = 0.02*exp(-obj.Y(1, k)/20);
            bhA = 4/(exp(-(obj.Y(1, k)+30)/8)+1);
            iA = obj.gA*obj.Y(2, k)^3*obj.Y(3, k)*Vk;
            %%% Kv %%%
            amKv = 10*(70-obj.Y(1, k))/(exp((70-obj.Y(1, k))/25)-1);
            bmKv = 2.5*exp(-obj.Y(1, k)/20);
            iKv = obj.gKv*obj.Y(4, k)^3*Vk;
            %%% Ka %%%
            mKa = 1/(1+exp(Vk/10));
            iKa = obj.gKa*(1-glu^3/(glu+obj.Kka^3))*mKa^4*Vk;
            %%% Ca %%%
            kb1 = obj.fB1/obj.bB1*obj.B1t/(1+obj.Y(1, 9)*obj.fB1/obj.bB1)^2;
            kb2 = obj.fB2/obj.bB2*obj.B2t/(1+obj.Y(1, 9)*obj.fB2/obj.bB2)^2;
            
            amCa = 10*(40-obj.Y(1, k))/(exp((40-obj.Y(1, k))/20)-1);
            bmCa = 13.2*exp(-obj.Y(1, k)/15);
            hCa = obj.KCa/(obj.KCa+obj.Y(9, k));
            eps = 2*obj.Y(1, k)*obj.F/(obj.R*obj.T)*10^(-3);
            iCa = obj.PCa*obj.Y(5, k)^2*hCa*eps*...
                (obj.Y(9, k)-obj.Cao*exp(-eps))/(1-exp(-eps));
            dCain = -iCa/(2*obj.F*obj.nu)*10^(-9);
            iex = obj.Iex*exp((obj.Y(1, k)+14)/70)*obj.Y(9, k)/(obj.Y(9, k)+3);
            ipump = obj.Ipump*obj.Y(9, k)/(obj.Y(9, k)+0.1);
            dCaef = (iex + ipump)/(2*obj.F*obj.nu)*10^(-9);
            %%% Na %%%
            amNa = 15*(35-obj.Y(1, k))/(exp((35-obj.Y(1, k))/32)-1);
            bmNa = 150*exp(-(30+obj.Y(1, k))/10);
            ahNa = 60480*exp(-(100+obj.Y(1, k))/5);
            bhNa = 2304/(exp((100-obj.Y(1, k))/75)+1);
            iNa = obj.gNa*obj.Y(6, k)^3*obj.Y(7, k)*(obj.Y(1, k)-obj.ENa);
            %%% glu %%%
            iGlu = 800*(exp((obj.Y(1, k)-obj.Eglu)/125)-1)/...
                (exp(-(obj.Y(1, k)+obj.Eglu)/125)+1)*glu^3/(glu^3+obj.Kglu^3);
            % iGlu = 2.56*glu*(obj.Y(1, k)-0);
            %%% leak %%%
            iL = obj.gl*(obj.Y(1, k)-obj.El);
			
            consts = [obj.C_S; obj.Btot; obj.kb; obj.kf; obj.bSB;...
                obj.fSB; obj.tfb; obj.V0];
            
            vars = [iCa; iA; iKa; iKv; iGlu; iL; iNa; amA; bmA; ahA; bhA; amKv;...
                bmKv; amCa; bmCa; amNa; bmNa; ahNa; bhNa; kb1; kb2; dCain; dCaef];
            
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
            
            if strcmp(obj.method, 'rk4')
            %%%%% runge-kutta 4
                k1 = obj.dt * D;
                k2 = obj.dt * f(obj.Y(:, k)+k1/2, vars, consts);
                k3 = obj.dt * f(obj.Y(:, k)+k2/2, vars, consts);
                k4 = obj.dt * f(obj.Y(:, k)+k3, vars, consts);
                k_tot = 1/6*(k1+2*k2+2*k3+k4);
            elseif strcmp(obj.method, 'euler')
            %%%%% forward euler
                k_tot = obj.dt*D;
            else
                k_tot = obj.dt*D;
            end
            
            
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
                        
            c = [iCa; iNa; iKv; iA; iKa; iGlu; iL; glu];
           
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function gaba = get_gaba(obj)
            gaba = obj.Y(10, :);
        end
        
        function tVector = get_tvec(obj)
            tVector = obj.t_vec;
        end
    end
end


%%%% F %%%%
function D = f(Y, vars, consts)

% consts
C_S = consts(1);
Btot = consts(2);
kb = consts(3);
kf = consts(4);
bSB = consts(5);
fSB = consts(6);
tfb = consts(7);
V0 = consts(8);

% vars
iCa = vars(1);
iA = vars(2);
iKa = vars(3);
iKv = vars(4);
iGlu = vars(5);
iL = vars(6);
iNa = vars(7);
amA = vars(8);
bmA = vars(9);
ahA = vars(10);
bhA = vars(11);
amKv = vars(12);
bmKv = vars(13);
amCa = vars(14);
bmCa = vars(15);
amNa = vars(16);
bmNa = vars(17);
ahNa = vars(18);
bhNa = vars(19);
kb1 = vars(20);
kb2 = vars(21);
dCain = vars(22);
dCaef = vars(23);
%Vth = vars(24);
%Vpre = vars(25);
%glu_inf = vars(26);

D = zeros(10, 1);
D(1) = -(iCa+iA+iKa+iKv+iGlu+iL+iNa)/C_S;
D(2) = amA*(1-Y(2))-bmA*Y(2);
D(3) = ahA*(1-Y(3))-bhA*Y(3);
D(4) = amKv*(1-Y(4))-bmKv*Y(4);
D(5) = amCa*(1-Y(5))-bmCa*Y(5);
D(6) = amNa*(1-Y(6))-bmNa*Y(6);
D(7) = ahNa*(1-Y(7))-bhNa*Y(7);
D(8) = Y(9)^3/(Y(9)^3+kb^3)*bSB*(Btot-Y(8))-Y(9)^3/(Y(9)^3+kf^3)*fSB*Y(8)*Y(9);
D(9) = 1/(1+kb1+kb2)*(dCain+dCaef+D(8));
D(10) = tfb*(Y(1)-V0-Y(10));
% if Vpre <= Vth
%     D(11) = 0;
% elseif abs(Vpre-Vth)>= 0
%     D(11) = (glu_inf-Y(11))/((1-glu_inf)*0.01);
% else
%     D(11) = 0;
% end



end
