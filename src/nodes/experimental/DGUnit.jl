@doc doc"""
```Julia
DGUnit(;I_r)
```


"""
@DynamicNode DGUnit(Pref, Qref) begin
    MassMatrix()
end  begin
end [[Pmeas, dPmeas], [Qmeas, dQmeas]] begin
    # Pmeas and Qmeas from wrapper!
    # Add limiter of Johannis
    # what about dQref and dPref?
    # neglected: delay and u_hold_in

# Local controller

# global inputs (for wrapper!)
control_mode =
power_limits =
power_functions =
pdirect =
qdirect =
EDI =
voltage_vector =
DG_number = # check with Paul if PowerDynamics node number can be extracted

# local inputs
PQ_local = # @Holm
power_vector =
S_max =
G_adj =

# define: control_mode,power_limits,power_functions,pdirect,qdirect, EDI,voltage_vector,DG_number, PQ_local,power_vector,S_max,G_adj
global_error =  GenerateControlErrors(control_mode,power_limits,power_functions,pdirect,qdirect, EDI,voltage_vector,DG_number, PQ_local,power_vector,S_max,G_adj)

e_P = global_error[1]
e_Q = global_error[2]

e_filter_P = DynamicRateLimiter(e_P) # also denoted as power_err
e_filter_Q = DynamicRateLimiter(e_Q) # @Holm: same limiter for P and Q?

################# I-controller block #####################
# e,power_lim_reach,quantyn,dPQ

power_lim_reach_P = # @Holm
power_lim_reach_Q = # @Holm
quantyn_P =
DPQ_P =
quantyn_Q =
DPQ_Q =


e_r_P = DGIntegralController(e_filter_P,power_lim_reach_P,quantyn_P,DPQ_P)
e_r_Q = DGIntegralController(e_filter_Q,power_lim_reach_Q,quantyn_Q,DPQ_Q)

e_out_P = e_r_P[1]
reset_P = e_r_P[2] # @Holm: what for?

e_out_Q = e_r_Q[1]
reset_Q = e_r_Q[2] # @Holm: what for?

############################# finish I-controller- function block ############

# Integrator ---> ADD Johannis limiter here??!
K_P = # @Holm
de_out_P = K_P * y_P  # @Holm: Initial condition y_P from step before?

K_Q = # @Holm
de_out_Q = K_Q * y_Q # @Holm: IC y_Q from step before?

########## Quantizer + PT1 (Holm) #################
 # check equivalence in Daniels model! Only limiter?

 ymax = 10e12
 ymin = -10e12

y_P = Limiter(y_P,ymax,ymin)
y_Q = Limiter(y_Q,ymax,ymin)

###############################################################
# DG Unit
###############################################################
    # y_P and y_Q from local controller above
    P_LC_ref = -0.5 * y_P
    Q_LC_ref = -0.5 * y_Q
    Pref_min = 0; # @Holm: WERT?
    Pref_max = 1; # @Holm: WERT?
    Qref_min = 0; # @Holm: WERT?
    Qref_max = 1; # @Holm: WERT?

P_ref = PT1Limiter(Pref,Pref_max,Pref_min,10.)[1] # delay of 0.01 missing
Q_ref = PT1Limiter(Qref,Qref_max,Qref_min,10.)[1] # delay of 0.01 missing

P_LC_ref = PT1Limiter(Pref,Pref_max,Pref_min,10.)[2]
Q_LC_ref = PT1Limiter(Qref,Qref_max,Qref_min,10.)[2]

    Va_amp = abs(u)
    Va_phase = angle(u)
    id = 2//3 * Pref / Va_amp
    iq = - 2//3 * im * Qref/V_amp
    # PLL
    K_pll = 0.1
    vq= imag(u)*cos(theta)-real(u)*sin(theta) # in simulink imag(u) and real(u) are calculted from V_amp and V_phase
    dvq = K_pll*theta
    ix=id*cos(theta)-iq*sin(theta);
    iy=id*sin(theta)+iq*cos(theta);
    I_r_xy = ix + im * iy
    # xy to abc missing
    du = i - I_r
    dPmeas = 0.
    dQmeas = 0.
end


function Limiter(var,max,min)
if var >= max
        var = max
    elseif var <= ymin
        var = min
end

function PT1Limiter(var,max,min,gain)
    if max <= min
        error("upper limit needs to be larger than lower limit")
    end

    if var >= max
        var = max
    elseif var <= var_min
        var = min
    else
        dvar = inp - gain * var # delay of 0.01 missing
    end
    (var, inp)
end


function DynamicRateLimiter(u)#,u_hold_in)

    ##dynamic Q rate limiter

    #Calcuting the rate of the limiter, only works for FixedStep = 0.005
    dy = 12500/8.0;

    # #to save the reference change
    # if u ~= u_last
    #     #change of the reference
    #     u_hold_out = u;
    # else
    #     #no change
    #     u_hold_out = u_hold_in;
    # end

    ####################################
    #If no step is detected
    if y_last == u
        y = u;
    #to avoid oscillation around zero
    elseif (u < dy && u > 0) || (u > -dy && u < 0)
        y = 0;
    #Regulation Value limitation
    elseif y_last> 1e6 || y_last < -1e6
        y = 1e6*u/abs(u);
    #step detection from pos to neg
    elseif u < 0 && u_last >=0
        y = -dy;
    #step detection from neg to pos
    elseif  u > 0 && u_last <=0
        y=dy;
    #step detection with no zero crossing
    elseif y_last < u_hold_out
        y = y_last +dy;
    #step detection with no zero crossing
    elseif y_last > u_hold_out
        y = y_last -dy;
    #should be without use, but matlab need this
    else
        y  = u;

    (y,u_hold_out)
end


# [e_P,e_Q,logfile] = fcn(control_mode,power_limits,power_functions,pdirect,qdirect, ...
#                         EDI,voltage_vector,DG_number, ...
#                         PQ_local,power_vector,S_max,G_adj)

function GenerateControlErrors(control_mode,power_limits,power_functions,pdirect,qdirect, EDI,voltage_vector,DG_number, PQ_local,power_vector,S_max,G_adj)
    logfile = zeros(1,4);
    # Generate local values
    voltage = voltage_vector[DG_number];
    EDI_local = EDI[DG_number,:];

    EDI_local = 1./EDI_local;
    EDI_local = EDI_local./sum(EDI_local);

    # Beginning and critical area of control action
    v_critical = 0.06; # 6
    v_pref = 0.03; # 3

    scale_max = 100;

    if abs(voltage-1) >= v_critical
        KUU = scale_max;
    elseif abs(voltage-1) >= v_pref
        KUU = (abs(voltage-1)-v_pref)/(v_critical-v_pref)*scale_max;
    else
        KUU = 0;
    end

    ref_weight = [EDI_local KUU];
    local_voltage = voltage;

##################################################
################### Controller ###################
##################################################


    e_P = 0;
    e_Q = 0;

    e = power_functions[1:length(power_functions)/2] - power_functions[length(power_functions)/2+1:length(power_functions)];

    # Generate control deviation depending on control mode
    # M = 1, P/Q/U controlled
    # M = 2, P/U controlled
    # M = 5, DG off
    # M = 8, setpoints (open loop control)
    # M = 11, P/Q
    # M = 51, Multinode PQU

    # ref_weight = [i_1 i_2 ... bus_voltage]

    v_scale = max(abs(e[[2 4]]));

    switch control_mode
        case 1
            # P/Q/U controlled
            #{
            eps = 1e6*0.05;
            if (abs(ref_weight[[1 2]]*e[[1 3]]) <= eps) ...
                    && (sum(abs(e[[1 3]])) > eps)
                ref_weight[ref_weight==min(ref_weight)]=0;
            elseif (abs(ref_weight[[1 2]]*e[[2 4]]) <= eps) ...
                    && (sum(abs(e[[2 4]])) > eps)
                ref_weight[ref_weight==min(ref_weight)]=0;
            end
            #}

            e_P = ref_weight[[1 2]]*e[[1 3]];
            e_Q = ref_weight[[1 2]]*e[[2 4]]+ ...
                  (-1)*(1-local_voltage)*v_scale*ref_weight[3];
            logfile = [ref_weight[[1 2]] e([1 3])'];
        case 2
            # P controlled, voltage controlled
            e_P = ref_weight[[1 2]]*e[[1 3]];
            e_Q = 0*ref_weight[[1 2]]*e[[2 4]]+(-1)*(1-local_voltage)*v_scale*ref_weight[3];
        case 5
            # DG off
            e_P = 0;
            e_Q = 0;
        case 8
            # direct ref
            e_P = -pdirect + PQ_local(1);
            e_Q = -qdirect + PQ_local(2);
        case 11
            # P/Q controlled, no voltage control
            e_P = ref_weight[[1 2]]*e[[1 3]];
            e_Q = ref_weight[[1 2]]*e[[2 4]]+ 0*(1-local_voltage)*v_scale*ref_weight[3];
        case 51
            # Multinode control

            # Construct refweight vector
            EDI = 1./EDI;
            ref_vec = [EDI[1:14,1]./sum(EDI[1:14,:],2) EDI[1:14,2]./sum(EDI[1:14,:],2) zeros(14,1)];
            #{
            sEDI = sum(EDI[1:14,:],2);
            ref_vec = zeros(14,3);
            for in = 1:14
                for ii = 1:2
                    ref_vec[in,ii] = EDI[in,ii]/sEDI[in];
                end
            end
            #}

            for i = 1:size(ref_vec,1)
                if abs(voltage_vector[i]-1) >= v_critical
                    KUU = scale_max;
                elseif abs(voltage_vector[i]-1) >= v_pref
                    KUU = (abs(voltage_vector[i]-1)-v_pref)/(v_critical-v_pref)*scale_max;
                else
                    KUU = 0;
                end
                ref_vec[i,3] = KUU;
            end

            # Generate control error vector
            Ve_P = ref_vec[:,[1 2]]*e[[1 3]];
            Ve_Q = ref_vec[]:,[1 2]]*e[[2 4])+(-1)*(1-voltage_vector).*v_scale.*ref_vec[:,3];

            s_crit = 0.05;
            # Truth vector
            S = (power_vector[1,:].^2+power_vector[2,:].^2).^0.5;
            T = 1-1/s_crit*abs(S_max-S')./S_max;
            T[T<0] = 0; T[DG_number] = 1; G_adj[DG_number,DG_number] = 1;
            TGAdj = G_adj[:,DG_number].*T;
            TGAdj = TGAdj./sum(TGAdj);

            e_P = TGAdj'*Ve_P;
            e_Q = TGAdj'*Ve_Q;
            logfile = TGAdj[[2 3 4 8]]';



    (e_P, e_Q, logfile)
end


function DGIntegralController(e,power_lim_reach,quantyn,dPQ)

    e_out = 0;
    reset = 0;

    # P_ref < P_DG_min => P_lim -1 (positive e_P_out still allowed)
    # P_ref > P_DG_max => P_lim 1  (negative e_P_out still allowed)
    if power_lim_reach == -1
        if e <= 0; e_out = 0; end
    elseif power_lim_reach == 1
        if e >= 0; e_out = 0; end
    else
        if quantyn == 0
            e_out = e;
        else
            if abs(e) < dPQ
                e_out = 0;
            else
                e_out = e;
            end
        end
    (e_out,reset)
end


export DGUnit
