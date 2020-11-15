using Distributions

function pchisq_ding(;x, F, theta)
	# Computes cumulative non-central chi-squared with f degrees of freedom and non-central parameter theta,
	chi2nc = x
  
	if x <= 0
		return missing
	end
	
	lam = theta / 2
	n = 1
	u = exp(-lam)
	v = u
	X2 = x / 2
	f2 = F / 2

	tmp1 = f2*log(X2) - X2 - gammln(xx=(f2 + 1))
  
	if tmp1 > 700
		return chilarge(z2=x, v2=F, k2=theta)
	end
  
	t = exp(tmp1)
  
	if t <= 1e-15
		return chilarge(z2=x, v2=F, k2=theta) 
	end

	term = v*t
	chi2nc = term
    
	# https://docs.julialang.org/en/v1/manual/control-flow/#man-loops
	while true
		if F + 2*n - x > 0
			bound = t*x / (F + 2*n - x)

			if bound < 1e-15
				return chi2nc
			end
		end
		
		u = u*lam / n
		v = v + u
		t = t*x / (F + 2*n)
		term = v*t
		chi2nc = chi2nc + term
		n = n + 1   
    end
end

function chilarge(;z2, v2, k2)
	z = 0.5*z2
	v = 0.5*v2
	k = 0.5*k2
  
	h = 1 - 2*(2*v + 2*k)*(2*v + 6*k) / ((2*v + 4*k)*(2*v + 4*k)*3)
	p = (2*v + 4*k) / ((2*v + 2*k)*(2*v + 2*k))
	m = (h - 1)*(1 - 3*h)
	A = 1 - h*p*(1 - h + 0.5*(2 - h)*m*p) - (2*z / (2*v + 2*k)) ^ h
	b = h*sqrt(2*p*(1 + m*p))
	rr = A / b
  
	x = cdf.(Normal(), rr)
	return (1 - x)
end

function gammln(;xx)
	# Define constants
	c1 = 76.1800917294715
	c2 = -86.5053203294168
	c3 = 24.0140982408309
	c4 = -1.23173957245015
	c5 = 1.20865097386618E-03
	c6 = -5.395239384953E-06
    
	if xx <= 0
		return 0
	end
	
	Y = xx
	x = xx
	tmp = x + 5.5
	tmp = tmp - (x + 0.5) * log(tmp)
	ser = 1.00000000019001

	Y = Y + 1
	ser = ser + c1 / Y
	Y = Y + 1
	ser = ser + c2 / Y
	Y = Y + 1
	ser = ser + c3 / Y
	Y = Y + 1
	ser = ser + c4 / Y
	Y = Y + 1
	ser = ser + c5 / Y
	Y = Y + 1
	ser = ser + c6 / Y
  
	return (-tmp + log(2.506628274631 * ser / x))
end

function BSMOption(;PutCallFlag,S,T,K,r,sigma,Greek)

	if PutCallFlag == "Call" && Greek == "Price" && T>0
		d1 = (log.(S/K) .+ (r .+ 0.5 .* sigma.^2).* T)./(sigma .* sqrt.(T))
		d2 = d1 .- sigma * sqrt.(T)
		return S .* cdf.(Normal(0,1), d1) .- K .* exp.(-r * T) .* cdf.(Normal(0,1), d2)
	end

	if PutCallFlag=="Call" && Greek=="Price" && T==0
		return max(S - K,0)
	end

	if PutCallFlag=="Put" && Greek=="Price" && T>0
		return -1 * BSMOption(PutCallFlag="Call",S=S,T=T,K=K,r=r,sigma=(-1 * sigma),Greek="Price")
	end

	if PutCallFlag=="Put" && Greek=="Price" && T==0
		return max(K - S,0)
	end
end

function CEVOption(;PutCallFlag,S,T,K,r,sigma,alpha)
	
	if alpha==1
		return BSMOption(PutCallFlag=PutCallFlag,S=S,T=T,K=K,r=r,sigma=sigma,Greek="Price")
	end
  
	if alpha!=1
		q = 0
		ny = (sigma^2)/(2*(r-q)*(alpha-1))*(exp(2*(r-q)*(alpha-1)*T)-1)
    
		a = (K*exp(-(r-q)*T))^(2*(1-alpha))/(((1-alpha)^2)*ny)
		b = 1/(1-alpha)
		c = (S^(2*(1-alpha)))/((1-alpha)^2*ny)
	end
  
	if alpha>0 && alpha<1 && PutCallFlag=="Call"
		return S*exp(-q*T)*(1-pchisq_ding(x=a, F=b+2, theta=c)) - K*exp(-r*T)*pchisq_ding(x=c, F=b, theta=a)
	end
  
	if alpha>0 && alpha<1 && PutCallFlag=="Put"
		return K*exp(-r*T)*(1-pchisq_ding(x=c, F=b, theta=a)) - S*exp(-q*T)*pchisq_ding(x=a, F=b+2, theta=c)
	end
  
	if alpha>1 && PutCallFlag=="Call"

		return S*exp(-q*T)*(1-pchisq_ding(x=c, F=-b, theta=a)) - K*exp(-r*T)*pchisq_ding(x=a, F=2-b, theta=c)
	end
  
	if alpha>1 && PutCallFlag=="Put"
		return K*exp(-r*T)*(1-pchisq_ding(x=a, F=2-b, theta=c)) - S*exp(-q*T)*pchisq_ding(x=c, F=-b, theta=a)
	end
end

# Load dataframe with test cases
using CSV
data = CSV.read(".../CEVTestData.csv")

using Dates
data[:TSStart] = DateTime.(data[:TSStart], Dates.DateFormat("yyyy-mm-dd HH:MM:SS"))
data[:TSEnd] = DateTime.(data[:TSEnd], Dates.DateFormat("yyyy-mm-dd HH:MM:SS"))
data[:Price] = float(data[:Price])

using DataFrames

for i = 1:nrow(data)
	if data[i,8] == -1
		data[i,9] = Dates.DateTime(Dates.now())
		data[i,8] = CEVOption(PutCallFlag=data[i,1],S=data[i,3],T=data[i,2],K=data[i,4],r=data[i,5],sigma=data[i,7],alpha=data[i,6])
		data[i,10] = Dates.DateTime(Dates.now())
	end
end
	
# Save Results
CSV.write(".../CEVTestData_Julia_With_Ding.csv", data)
