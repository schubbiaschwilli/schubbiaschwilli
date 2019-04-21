PROC FCMP OUTLIB=sasuser.funktionen.paket;

	FUNCTION BSMOption(PutCallFlag $, S, T, K, r, sigma, Greek $);

		IF PutCallFlag="Call" AND Greek="Price" AND T>0 THEN DO;
			d1=(log(S/K) + (r + sigma**2  / 2) * T) / (sigma * sqrt(T));
			d2 = d1 - sigma * sqrt(T);
			value = (S * probnorm(d1) - K * exp(-r*T) * probnorm(d2));
		END;		

		IF PutCallFlag="Call" AND Greek="Price" AND T=0 THEN DO;
			value = max(S-K, 0);
		END;	

		IF PutCallFlag="Call" AND Greek="Delta" AND T>0 THEN DO;
			d1=(log(S/K) + (r + sigma**2  / 2) * T) / (sigma * sqrt(T));
			value = probnorm(d1);
		END; 

		IF PutCallFlag="Call" AND Greek="Delta" AND T=0 THEN DO;
			value = 0;
		END;	

		IF PutCallFlag="Put" AND Greek="Price" AND T>0 THEN DO;
			value = -1 * BSMOption("Call", S, T, K, r, -1 * sigma, "Price");
		END;

		IF PutCallFlag="Put" AND Greek="Price" AND T=0 THEN DO;
			value = max(K-S, 0);
		END;

		IF PutCallFlag="Put" AND Greek="Delta" THEN DO;
			value = -1 * BSMOption("Call", S, T, K, r, -1 * sigma, "Delta");
		END;

		return(value);

	ENDSUB;
RUN; 

OPTIONS CMPLIB=sasuser.funktionen;

data Results(keep=S Hedge PayOff Hedgeerror);

	PutCallFlag = "Call";
	S_0 = 100;
	T = 0.5;
	K = S_0 * 1.15;
	r = 0.05;
	mu = 0.1;
	sigma = 0.2;
	Npaths = 5000;
	Nhedgepoints = 52;

	call streaminit(123456);

	dt = T/Nhedgepoints;

	DO j=1 TO NPaths; ** Loop over NPaths;
		
		** Init ; 
		S= S_0;
		V= BSMOption(PutCallFlag, S, T, K, r, sigma, "Price"); * initial investment ;
		a= BSMOption(PutCallFlag, S, T, K, r, sigma, "Delta"); * stock position = delta;
		b= V - a * S; * rest in bank self - fin. Cond. ;

		DO i=1 TO Nhedgepoints; * Loop over Nhedgepoints ;
			eps = rand('NORMAL'); * eps= random.normal () simulate outcome of N(0 ,1);
			S=S* exp ((mu -0.5 * sigma**2) *dt + sigma * sqrt (dt)* eps);
			V= a* S+b* exp (r*dt);
 			a= BSMOption(PutCallFlag, S, T-i*dt, K, r, sigma ,"Delta");
			b=V-a*S;
		end;
		
		PayOff=BSMOption(PutCallFlag, S, 0, K, r, sigma, "Price");
		Hedge=V;
		Hedgeerror=PayOff-Hedge;
		output;
	end;
run;

Title "BSM: Delta-Hedging; Overview PayOfff vs Hedge";
Symbol1 V=Point C=Red; 
Symbol2 V=Point C=Green;

* Horizontal axis;
axis1  label=('S(t)');
***vertical axis;
axis2 label=('Payoff / Hedge');

Proc Gplot Data=Results;
	* <y-Variable>*<x-Variable> ;
	Plot Hedge*S Payoff*S / overlay haxis=axis1 vaxis=axis2;
Run; Quit; 

Title "Overview Hedgeerror";
PROC UNIVARIATE Data=Results PLOT NORMAL;
     VAR Hedgeerror;
RUN;

