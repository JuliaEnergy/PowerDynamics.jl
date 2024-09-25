@mtkmodel AVRTypeI begin
    @components begin
        vh = RealInput() # generator terminal voltage
        vref = RealInput() # reference voltage
        vf = RealOutput(guess=1) # field voltage
    end
    @parameters begin
        Ka, [description="amplifier gain"]
        Ke, [description="field circuit integral deviation"]
        Kf, [description="stabilizer gain"]
        Ta, [description="amplifier time constant"]
        Tf, [description="stabilizer time constant"]
        Te, [description="field circuit time constant"]
        Tr, [description="measurement time constant"]
        vr_min, [description="minimum regulator voltage"]
        vr_max, [description="maximum regulator voltage"]
        Ae, [description="1st ceiling coeff"]
        Be, [description="2st ceiling coeff"]
    end
    @variables begin
        vr1(t), [guess=1,description="intenral amplifier state TODO: implement AVR limiter"]
        vr2(t), [guess=0,description="internal feedback state"]
        vm(t), [guess=1, description="terminal voltage measurement (lagged)"]
        vfceil(t), [description="ceiled field voltage"]
    end
    @equations begin
        vfceil ~ Ae * exp(Be * abs(vf.u))
        Dt(vf.u) ~ -(vf.u  * (Ke+vfceil) - vr1)/Te
        Dt(vr1) ~ (Ka*(vref.u - vm - vr2 - Kf/Tf*vf.u) - vr1)/Ta
        Dt(vr2) ~ -(Kf/Tf*vf.u + vr2)/Tf
        Dt(vm) ~ (vh.u - vm)/Tr
    end
end
