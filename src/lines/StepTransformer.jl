"""
```Julia
StepTransformer(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl)
```
A transformer representation uses the Π model,
assuming an ideal transformer in series with an admittance.
The admittance is here taken to be on the high-voltage side.

#Sr      = 25.; #Bemessungsscheinleistung in MVA
#Uos     = 110.; # Bemessungsspannung Oberspannungsseite in kV
#Uus     = 20.; # Bemessungsspannung Unterspannungsseite in kV
#uk      = 12.; # Kurzschlussspannung in %
#Pvk     = 25.; # Kupferverluste in kW
#Pvl     = 0.; # Eisenverluste in kW
#iLeer   = 0.; # Leerlaufstrom in %
#k     = 0. ; # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
# stufe_us = 0.625 ; # Stufung pro Stufenschalterstellung in %
#ssp = 10. # Stufenschalterposition

"""
@Line StepTransformer(Uos, Uus, k, ssp, stufe_us, Sr, uk, Pvk, Pvl) begin
    ue=(Uos/Uus)*exp(im*k*(π/6));
    # Übersetzungsverhältnis
    ZL = ((uk/100)*(Uos*1000)^2)/(Sr*10^6);
    # Berechnung der Längsimpedanz (Betrag)
    Rl = ((Pvk*10^3)*(Uos*10^3)^2)/((Sr*10^6)^2);
    # Berechnung der Längswiderstände
    Xl = sqrt((ZL^2)-(Rl^2));
    # Berechnung der Längsreaktanz
    LaengsImp = (Rl+im*Xl);
    # Berechnung der Längsimpedanz (als komplexe Zahl)

    if iLeer == 0
        Zh=0;
        RFE=0;
        Xh=0;
        QuerImp=0;
        Ys=0;
    else
        Zh =(Uos*10^3)^2/(Sr*10^6*(iLeer/100));
        # Berechnung der Querimpedanz (Betrag)
        RFE=((Uos*10^3)^2)/(Pvl*10^3);
        # Berechnung des Querwiderstandes (über Eisenverluste)
        Xh =Zh/sqrt(1-((Zh/RFE)^2));
        # Berechnung der Querreaktanz
        QuerImp =1/(1/RFE + 1/(Xh*im));
        # Berechnung der Querimpedanz (als komplexe Zahl)
        Ys=(1/QuerImp)/2;
    end

    # Transformation T-Ersatzschaltbild in Vierpoldarstellung
    YB = Ys;
    YC = 1/LaengsImp;

    # %%%%% YAA und YBB vertauschen bei Änderungen des Spannungsbezugs beim
    # %%%%% Slackknoten %%%%%
    YAA = (YB + YC); # A entspricht OS-Seite
    YAB = -YC * ue;
    YBA = -YC * conj(ue);
    YBB = (YB + YC)*abs(ue^2); # B entspricht der US-Seite

    # %%%%% Speichern der Vierpol-Darstellung in einer Betriebsmittelmatrix
    # %%%%% YTT (Diagonal) %%%%%
    YTT = ones(Complex, 2, 2)
    YTT[1  , 1]= YAA; # Die Spannungsseite wo Slack
    YTT[2,  1]= YAB;
    YTT[1  ,2]= YBA;
    YTT[2, 2]= YBB; # Die Spannungsseite wo kein Slack

    #### Stufung unterspannungsseitig (Längsregler = Phasenverschiebung konstant)

    trafo_step = 1. / (1. + ssp * stufe_us / 100.)
    Y = YTT .* [1 trafo_step; trafo_step trafo_step^2]

    current_vector = Y * [source_voltage, destination_voltage]
end

export StepTransformer
