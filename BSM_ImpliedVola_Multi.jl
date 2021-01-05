# Creates C:\_Data\h_da\Derivate\Kurse\Modells\CEV\Julia\BSM_ImpliedVola\Blog\Plots\ImpliedVola_*Handelstag*.png

using Distributions
using Plots
using StatsPlots
using DataFrames
using Dates
using RData
using CSV
using HypothesisTests

cd("C:/_Data/h_da/Derivate/Kurse/Modells/CEV/Julia/BSM_ImpliedVola/Blog/Output_Multi")

# Load data
data = load("C:/_Data/h_da/_Blog/Themen/Options/BSM_ImpliedVola_AllSingleOptions/CreataDataFile/EurexOptionsDaxPlus.RData", convert=true)

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

	if (PutCallFlag=="Call" || PutCallFlag=="Put") && Greek=="Vega" && T>0
		d1 = (log.(S/K) .+ (r .+ 0.5 .* sigma.^2).* T)./(sigma .* sqrt.(T))
		return S * exp.(-r * T) * pdf.(Normal(0,1), d1) * sqrt.(T)
	end

	if PutCallFlag=="Put" && Greek=="Price" && T==0
		return max(K - S, 0)
	end
end

function BSMOptionImpliedVolatility(;PutCallFlag,S,T,K,r,Price)
	epsilon = 10^-6
	sigma_lower = epsilon
	sigma_upper = 10
	Price_lower = BSMOption(PutCallFlag=PutCallFlag,S=S,T=T,K=K,r=r,sigma=sigma_lower,Greek="Price")
	Price_upper = BSMOption(PutCallFlag=PutCallFlag,S=S,T=T,K=K,r=r,sigma=sigma_upper,Greek="Price")
	
	if Price < BSMOption(PutCallFlag=PutCallFlag,S=S,T=T,K=K,r=r,sigma=sigma_lower,Greek="Price")
			return 0
		else
			while abs(sigma_upper-sigma_lower) >= epsilon
				Price_new = BSMOption(PutCallFlag=PutCallFlag,S=S,T=T,K=K,r=r,sigma=(sigma_lower+sigma_upper)/2,Greek="Price")

				if  (Price_lower-Price)*(Price_new-Price) < 0
    					sigma_upper = (sigma_lower+sigma_upper)/2
    					Price_upper = Price_new
    				else
    					sigma_lower = (sigma_lower+sigma_upper)/2
    					Price_lower = Price_new
				end
			end

			return (sigma_lower+sigma_upper)/2
	end
end

# Create dataframe
df = DataFrame(data["EurexOptionsDax"])
df.Handelstag_mod = repeat([Date(df[1,1], UTC)], nrow(df))

for i = 1:nrow(df)
	df.Handelstag_mod[i] = Date(df[i,1], UTC)
end

df[!, :Handelstag] = df.Handelstag_mod
df.ImpliedVola_old = df.ImpliedVola
df.ImpliedVola = repeat([-.999], nrow(df))

Handelstage = df.Handelstag
Handelstage = unique(Handelstage)
Handelstage = DataFrame(Handelstag = Handelstage)
Handelstage.Mean_Calls = repeat([-.999], nrow(Handelstage))
Handelstage.Mean_Puts = repeat([-.999], nrow(Handelstage))
Handelstage.StdDev_Calls = repeat([-.999], nrow(Handelstage))
Handelstage.StdDev_Puts = repeat([-.999], nrow(Handelstage))
Handelstage.RelStdDev_Calls = repeat([-.999], nrow(Handelstage))
Handelstage.RelStdDev_Puts = repeat([-.999], nrow(Handelstage))
Handelstage.NumberOfCalls = repeat([-1], nrow(Handelstage))
Handelstage.NumberOfPuts = repeat([-1], nrow(Handelstage))

for i = 1:nrow(Handelstage)
	println(Dates.Time(Dates.now()), " : ", i)
	subset = df[(df.Handelstag .== Handelstage[i,:Handelstag]) .& (0.8 .< df.Moneyness .< 1.2) .& (1/12 .< df.t_delta .< 1),:]

	# Calc Implied Vola	
	for j = 1:nrow(subset)
		subset[j,:ImpliedVola] = BSMOptionImpliedVolatility(PutCallFlag=subset[j,:OptionType],S=subset[j,:SchlusspreisBasiswert],T=subset[j,:t_delta],K=subset[j,:StrikePrice],r=log(1+subset[j,:EONIA]),Price=subset[j,:TaeglicherAbrechnungspreis])
	end
	
	# subset = 
	subset_Calls = subset[(subset.OptionType .=="Call"),:]
	subset_Puts = subset[(subset.OptionType .=="Put"),:]

	Handelstage[i,:Mean_Calls] = mean(subset_Calls.ImpliedVola)
	Handelstage[i,:Mean_Puts] = mean(subset_Puts.ImpliedVola)
	Handelstage[i,:StdDev_Calls] = std(subset_Calls.ImpliedVola)
	Handelstage[i,:StdDev_Puts] = std(subset_Puts.ImpliedVola)
	Handelstage[i,:NumberOfCalls] = nrow(subset_Calls)
	Handelstage[i,:NumberOfPuts] = nrow(subset_Puts)
	Handelstage[i,:RelStdDev_Calls] = Handelstage[i,:StdDev_Calls] / Handelstage[i,:NumberOfCalls]
	Handelstage[i,:RelStdDev_Puts] = Handelstage[i,:StdDev_Puts] / Handelstage[i,:NumberOfPuts]

	# Create plot
	filename = string("ImpliedVola_", Handelstage[i,:Handelstag], ".png") 

	if !isfile(filename)
		# Single Plots
		Scatter_Calls = scatter(subset_Calls.StrikePrice, subset_Calls.t_delta, subset_Calls.ImpliedVola, title=string("Handelstag ",  Handelstage[i,:Handelstag], " : Calls"), xlabel="Strike", ylabel="T", zlabel="Implied Volatility", zcolor=-subset_Calls.ImpliedVola, legend=false)
		Scatter_Puts = scatter(subset_Puts.StrikePrice, subset_Puts.t_delta, subset_Puts.ImpliedVola, title="Puts", xlabel="Strike", ylabel="T", zlabel="Implied Volatility", zcolor=-subset_Puts.ImpliedVola, legend=false)
		Histogram_Calls = histogram(subset_Calls.ImpliedVola, title=string("Mean : ", round(mean(subset_Calls.ImpliedVola), digits=4)), xlabel="Implied Volatility", ylabel="Counts", legend=false)
		Histogram_Puts = histogram(subset_Puts.ImpliedVola, title=string("Mean : ", round(mean(subset_Puts.ImpliedVola), digits=4)), xlabel="Implied Volatility", ylabel="Counts", legend=false)
	
		# Combined 2x2 Plot
		plot(Scatter_Calls, Scatter_Puts, Histogram_Calls, Histogram_Puts, layout=(2,2), size=(800,600))
	
		savefig(filename)
	end 
end

# Save Results Handelstage
CSV.write("MeanImpliedVolas.csv", Handelstage)

# Plot Overview
	# Plot Means
	TSMeans = plot(Handelstage.Handelstag, Handelstage.Mean_Calls, label="Calls", ylabel="Mean(sigma)", size=(800, 600))
	TSMeans = plot!(Handelstage.Handelstag, Handelstage.Mean_Puts, label="Puts")

	# Plot Scatter
	t = string("Correlation: ", round(cor(Handelstage.Mean_Calls, Handelstage.Mean_Puts), digits=4))
	Scatter = scatter(Handelstage.Mean_Calls, Handelstage.Mean_Puts, xlabel="Calls Mean(sigma)", ylabel="Puts Mean(sigma)", title=t, legend=false)

	# Histograms
	Histogram_Calls = histogram(Handelstage.Mean_Calls, title=string("Calls; Mean : ", round(mean(Handelstage.Mean_Calls), digits=4)), xlabel="Mean(sigma)", ylabel="Counts", legend=false)
	Histogram_Puts = histogram(Handelstage.Mean_Puts, title=string("Puts; Mean : ", round(mean(Handelstage.Mean_Puts), digits=4)), xlabel="Mean(sigma)", ylabel="Counts", legend=false)

	plot(TSMeans, Histogram_Calls, Scatter, Histogram_Puts, layout=(2,2), size=(1200,800))

	savefig("Overview_Sigmas.png")

# Plot TS StdDev
	plot(Handelstage.Handelstag, Handelstage.StdDev_Calls, label="Calls", ylabel="StdDev(sigma)", size=(800, 600))
	plot!(Handelstage.Handelstag, Handelstage.StdDev_Puts, label="Puts")

	savefig("StdDev_Sigmas.png")

# Plot TS RelStdDev / relative standard deviation (RSD)
	plot(Handelstage.Handelstag, Handelstage.RelStdDev_Calls, label="Calls", ylabel="Relative StdDev(sigma)", size=(800, 600))
	plot!(Handelstage.Handelstag, Handelstage.RelStdDev_Puts, label="Puts")

	savefig("RelStdDev_Sigmas.png")


# Histogram and QQPlot Deltas
Mean_Deltas = Handelstage.Mean_Calls - Handelstage.Mean_Puts
Histogram_Deltas = histogram(Mean_Deltas, title=string("Mean : ", round(mean(Mean_Deltas), digits=4)), xlabel="Mean(delta(sigmas))", ylabel="Counts", legend=false, size=(800,600))
QQPlot_Deltas = plot(qqnorm(Mean_Deltas))

plot(Histogram_Deltas, QQPlot_Deltas, layout=(1,2), size=(800,400))
savefig("Deltas.png")

OneSampleTTest(Mean_Deltas)
JarqueBeraTest(Mean_Deltas)

# Test Histograms - Didn't work
#histogram([Handelstage.Mean_Calls,Handelstage.Mean_Puts])

#histogram(Handelstage.Mean_Calls, bar_width=1/250, bar_edges=false)
#histogram!(Handelstage.Mean_Puts, bar_width=1/250) # , bar_edges=false)