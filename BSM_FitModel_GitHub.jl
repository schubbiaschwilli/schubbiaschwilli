using Distributions
using Optim
using Plots
using DataFrames
using Dates
using RData
using CSV
using HypothesisTests
using StatsPlots
using GLM

# Load data
data = load(".../EurexOptionsDaxPlus.RData", convert=true)

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
		return max(K-S, 0)
	end
end

function BSMLossFunction(sigma, data)
	RSS = 0
	for i = 1:nrow(data)
		RSS += (data[i,:TaeglicherAbrechnungspreis] - BSMOption(PutCallFlag=data[i,:OptionType],S=data[i,:SchlusspreisBasiswert],T=data[i,:t_delta],K=data[i,:StrikePrice],r=log(1+data[i,:EONIA]),sigma=sigma[1],Greek="Price"))^2
	end
	return RSS
end

# Create dataframe
df = DataFrame(data["EurexOptionsDax"])
df.Handelstag_mod = repeat([Date(df[1,1], UTC)], nrow(df))

for i = 1:nrow(df)
	df.Handelstag_mod[i] = Date(df[i,1], UTC)
end

df[!, :Handelstag] = df.Handelstag_mod

Handelstage = df.Handelstag 
Handelstage = unique(Handelstage)
Handelstage = DataFrame(Handelstag = Handelstage)

Handelstage.SchlusspreisBasiswert = repeat([-.999], nrow(Handelstage))
Handelstage.EONIA = repeat([-.999], nrow(Handelstage))

Handelstage.All_sigma_mean = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_sigma_mean = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_sigma_mean = repeat([-.999], nrow(Handelstage))
Handelstage.All_sigma_fit = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_sigma_fit = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_sigma_fit = repeat([-.999], nrow(Handelstage))
Handelstage.All_RSS = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RSS = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RSS = repeat([-.999], nrow(Handelstage))
Handelstage.All_RMSE = -repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RMSE = -repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RMSE = -repeat([-.999], nrow(Handelstage))
Handelstage.All_RMSE_S = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RMSE_S = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RMSE_S = repeat([-.999], nrow(Handelstage))
Handelstage.All_NumberOfDatasets = repeat([-1], nrow(Handelstage))
Handelstage.Calls_NumberOfDatasets = repeat([-1], nrow(Handelstage))
Handelstage.Puts_NumberOfDatasets = repeat([-1], nrow(Handelstage))
Handelstage.All_RSS_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RSS_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RSS_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.All_RMSE_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RMSE_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RMSE_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.All_RMSE_S_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Calls_RMSE_S_SteadyState = repeat([-.999], nrow(Handelstage))
Handelstage.Puts_RMSE_S_SteadyState = repeat([-.999], nrow(Handelstage))

Handelstage[1,:All_RSS_SteadyState] = NaN
Handelstage[1,:Calls_RSS_SteadyState] = NaN
Handelstage[1,:Puts_RSS_SteadyState] = NaN
Handelstage[1,:All_RMSE_SteadyState] = NaN
Handelstage[1,:Calls_RMSE_SteadyState] = NaN
Handelstage[1,:Puts_RMSE_SteadyState] = NaN
Handelstage[1,:All_RMSE_S_SteadyState] = NaN
Handelstage[1,:Calls_RMSE_S_SteadyState] = NaN
Handelstage[1,:Puts_RMSE_S_SteadyState] = NaN

for i = 1:nrow(Handelstage)

	println(Dates.Time(Dates.now()), " : ", i)

	subset_All = df[(df.Handelstag .== Handelstage[i,:Handelstag]) .& (0.8 .< df.Moneyness .< 1.2) .& (1/12 .< df.t_delta .< 1),:]
	subset_Calls = df[(df.Handelstag .== Handelstage[i,:Handelstag]) .& (df.OptionType .=="Call") .& (0.8 .< df.Moneyness .< 1.2) .& (1/12 .< df.t_delta .< 1),:]
	subset_Puts = df[(df.Handelstag .== Handelstage[i,:Handelstag]) .& (df.OptionType .=="Put") .& (0.8 .< df.Moneyness .< 1.2) .& (1/12 .< df.t_delta .< 1),:]	

	Handelstage[i,:SchlusspreisBasiswert] = subset_Calls[1,:SchlusspreisBasiswert]
	Handelstage[i,:EONIA] = subset_Calls[1,:EONIA]
	
	Handelstage[i,:All_NumberOfDatasets] = nrow(subset_All)
	Handelstage[i,:Calls_NumberOfDatasets] = nrow(subset_Calls)
	Handelstage[i,:Puts_NumberOfDatasets] = nrow(subset_Puts)
	Handelstage[i,:All_sigma_mean] = mean(subset_All.ImpliedVola)
	Handelstage[i,:Calls_sigma_mean] = mean(subset_Calls.ImpliedVola)
	Handelstage[i,:Puts_sigma_mean] = mean(subset_Puts.ImpliedVola)

	# Model calibration
	res_All = optimize(sigma -> BSMLossFunction(sigma, subset_All), [Handelstage[i,:All_sigma_mean]], Newton())
	Handelstage[i,:All_sigma_fit] = Optim.minimizer(res_All)[1]
	Handelstage[i,:All_RSS] = BSMLossFunction(Handelstage[i,:All_sigma_fit], subset_All)
	Handelstage[i,:All_RMSE] = sqrt(Handelstage[i,:All_RSS] / nrow(subset_All))
	Handelstage[i,:All_RMSE_S] = Handelstage[i,:All_RMSE] / subset_All[1,:SchlusspreisBasiswert]

	res_Calls = optimize(sigma -> BSMLossFunction(sigma, subset_Calls), [Handelstage[i,:Calls_sigma_mean]], Newton())
	Handelstage[i,:Calls_sigma_fit] = Optim.minimizer(res_Calls)[1]
	Handelstage[i,:Calls_RSS] = BSMLossFunction(Handelstage[i,:Calls_sigma_fit], subset_Calls)
	Handelstage[i,:Calls_RMSE] = sqrt(Handelstage[i,:Calls_RSS] / nrow(subset_Calls))
	Handelstage[i,:Calls_RMSE_S] = Handelstage[i,:Calls_RMSE] / subset_Calls[1,:SchlusspreisBasiswert]
	
	res_Puts = optimize(sigma -> BSMLossFunction(sigma, subset_Puts), [Handelstage[i,:Puts_sigma_mean]], Newton())
	Handelstage[i,:Puts_sigma_fit] = Optim.minimizer(res_Puts)[1]
	Handelstage[i,:Puts_RSS] = BSMLossFunction(Handelstage[i,:Puts_sigma_fit], subset_Puts)
	Handelstage[i,:Puts_RMSE] = sqrt(Handelstage[i,:Puts_RSS] / nrow(subset_Puts))
	Handelstage[i,:Puts_RMSE_S] = Handelstage[i,:Puts_RMSE] / subset_Puts[1,:SchlusspreisBasiswert]

	if i>1
		Handelstage[i,:All_RSS_SteadyState] = BSMLossFunction(Handelstage[i-1,:All_sigma_mean], subset_All)
		Handelstage[i,:All_RMSE_SteadyState] = sqrt(Handelstage[i,:All_RSS_SteadyState] / nrow(subset_All))
		Handelstage[i,:All_RMSE_S_SteadyState] = Handelstage[i,:All_RMSE_SteadyState] / subset_All[1,:SchlusspreisBasiswert]

		Handelstage[i,:Calls_RSS_SteadyState] = BSMLossFunction(Handelstage[i-1,:Calls_sigma_mean], subset_Calls)
		Handelstage[i,:Calls_RMSE_SteadyState] = sqrt(Handelstage[i,:Calls_RSS_SteadyState] / nrow(subset_Calls))
		Handelstage[i,:Calls_RMSE_S_SteadyState] = Handelstage[i,:Calls_RMSE_SteadyState] / subset_Calls[1,:SchlusspreisBasiswert]

		Handelstage[i,:Puts_RSS_SteadyState] = BSMLossFunction(Handelstage[i-1,:Puts_sigma_mean], subset_Puts)
		Handelstage[i,:Puts_RMSE_SteadyState] = sqrt(Handelstage[i,:Puts_RSS_SteadyState] / nrow(subset_Puts))
		Handelstage[i,:Puts_RMSE_S_SteadyState] = Handelstage[i,:Puts_RMSE_SteadyState] / subset_Puts[1,:SchlusspreisBasiswert]
	end
end

# Export Data/Resuluts
CSV.write("BSM_FitModel_Results.csv", Handelstage, header=true)

# Plot Fitted sigmas Calls&Puts
	plot(Handelstage.Handelstag, Handelstage.Calls_sigma_fit, label="Calls", ylabel="Volatility/sigma", title="BSM-Fitted", size=(800,600))
	plot!(Handelstage.Handelstag, Handelstage.Puts_sigma_fit, label="Puts")
	plot!(Handelstage.Handelstag, Handelstage.All_sigma_fit, label="All")
	savefig("BSM_FitModel_TS_All.png")

# Plot Fitted sigmas Calls&Puts
	plot(Handelstage.Handelstag, Handelstage.Calls_sigma_fit, label="Calls", ylabel="Volatility/sigma", title="BSM-Fitted", size=(800,600))
	plot!(Handelstage.Handelstag, Handelstage.Puts_sigma_fit, label="Puts")
	savefig("BSM_FitModel_TS.png")

# Plot Scatter Calls vs Puts
	t = string("Correlation: ", round(cor(Handelstage.Calls_sigma_fit, Handelstage.Puts_sigma_fit), digits=4))
	Scatter = scatter(Handelstage.Calls_sigma_fit, Handelstage.Puts_sigma_fit, xlabel="Calls Fitted sigma", ylabel="Puts Fitted sigma", title=t, legend=false, size=(800,600))
	savefig("BSM_FitModel_Scatter.png")

	cor(Handelstage.All_sigma_fit, Handelstage.Calls_sigma_fit)
	cor(Handelstage.All_sigma_fit, Handelstage.Puts_sigma_fit)

### Plot Fit Calls minus Puts
	Fit_Calls_minus_Puts = Handelstage.Calls_sigma_fit - Handelstage.Puts_sigma_fit
	Histogram_Deltas = histogram(Fit_Calls_minus_Puts, title=string("Mean : ", round(mean(Fit_Calls_minus_Puts), digits=4)), xlabel="delta(Fitted Calls - Fitted Puts))", ylabel="Counts", legend=false)
	QQPlot_Deltas = plot(qqnorm(Fit_Calls_minus_Puts))

	plot(Histogram_Deltas, QQPlot_Deltas, layout=(1,2), size=(800,400))
	savefig("Fit_Calls_minus_Puts.png")

### Tests Fit Calls minus Puts
	OneSampleTTest(vec(Fit_Calls_minus_Puts))
	JarqueBeraTest(vec(Fit_Calls_minus_Puts))

# Plot sigmas Mean vs Fit
	Calls_Mean_vs_Fit = plot(Handelstage.Handelstag, Handelstage.Calls_sigma_mean, label="Mean", ylabel="Volatility/sigma", title="Calls")
	Calls_Mean_vs_Fit = plot!(Handelstage.Handelstag, Handelstage.Calls_sigma_fit, label="BSM-Fitted")

	Puts_Mean_vs_Fit = plot(Handelstage.Handelstag, Handelstage.Puts_sigma_mean, label="Mean", ylabel="Volatility/sigma", title="Puts")
	Puts_Mean_vs_Fit = plot!(Handelstage.Handelstag, Handelstage.Puts_sigma_fit, label="BSM-Fitted")

	plot(Calls_Mean_vs_Fit, Puts_Mean_vs_Fit, layout=(2,1), size=(800,600))
	savefig("BSM_Mean_vs_FitModel.png")

# Scatterplot Fit vs Mean
	t = string("Calls; Correlation: ", round(cor(Handelstage.Calls_sigma_mean, Handelstage.Calls_sigma_fit), digits=4))
	Calls_Scatter = scatter(Handelstage.Calls_sigma_mean, Handelstage.Calls_sigma_fit, xlabel="Mean(sigma)", ylabel="sigma fitted", title=t, legend=false)

	t = string("Puts; Correlation: ", round(cor(Handelstage.Puts_sigma_mean, Handelstage.Puts_sigma_fit), digits=4))
	Puts_Scatter = scatter(Handelstage.Puts_sigma_mean, Handelstage.Puts_sigma_fit, xlabel="Mean(sigma)", ylabel="sigma fitted", title=t, legend=false)

	plot(Calls_Scatter, Puts_Scatter, layout=(1,2), size=(800,400))
	savefig("BSM_Mean_vs_FitModel_Correl.png")

# Plot TS RMSE/S Current vs SteadyState
	Calls = plot(Handelstage.Handelstag, Handelstage.Calls_RMSE_S, title="Calls", label="Current", ylabel="RMSE/S")
	Calls = plot!(Handelstage.Handelstag, Handelstage.Calls_RMSE_S_SteadyState, label="Steady State Prediction")

	Puts = plot(Handelstage.Handelstag, Handelstage.Puts_RMSE_S, title="Puts", label="Current", ylabel="RMSE/S")
	Puts = plot!(Handelstage.Handelstag, Handelstage.Puts_RMSE_S_SteadyState, label="Steady State Prediction")

	plot(Calls, Puts, layout=(2,1), size=(800,600))
	savefig("RMSE_S_TS.png")

# Histogram Deltas Current vs Previous
	Fit_Calls_Pred = Handelstage.Calls_sigma_fit[2:end,:] - Handelstage.Calls_sigma_fit[1:end-1,:]
	Histogram_Deltas_Calls = histogram(Fit_Calls_Pred, title=string("Mean : ", round(mean(Fit_Calls_Pred), digits=4)), xlabel="Delta Sigma Calls Current vs Previous", ylabel="Counts", legend=false)
	QQPlot_Deltas_Calls = plot(qqnorm(vec(Fit_Calls_Pred)))

	Fit_Puts_Pred = Handelstage.Puts_sigma_fit[2:end,:] - Handelstage.Puts_sigma_fit[1:end-1,:]
	Histogram_Deltas_Puts = histogram(Fit_Puts_Pred, title=string("Mean : ", round(mean(Fit_Puts_Pred), digits=4)), xlabel="Delta Sigma Puts Current vs Previous", ylabel="Counts", legend=false)
	QQPlot_Deltas_Puts = plot(qqnorm(vec(Fit_Puts_Pred)))

	plot(Histogram_Deltas_Calls, QQPlot_Deltas_Calls, Histogram_Deltas_Puts, QQPlot_Deltas_Puts, layout=(2,2), size=(800,600))
	savefig("Histogram_QQPlot_Deltas_CurrentVsPrevious.png")
	
	### Tests Fit Deltas Current vs Previous
	OneSampleTTest(vec(Fit_Calls_Pred))
	JarqueBeraTest(vec(Fit_Calls_Pred))

	OneSampleTTest(vec(Fit_Puts_Pred))
	JarqueBeraTest(vec(Fit_Puts_Pred))

# Plot Scatter Current Vs Pred / https://juliastats.org/GLM.jl/latest/examples/
	# Calls
	OLSData = DataFrame(Current=vec(Handelstage.Calls_sigma_fit[1:end-1,:]), Pred=vec(Handelstage.Calls_sigma_fit[2:end,:]))
	reg = lm(@formula(Pred ~ Current), OLSData)
	alpha = GLM.coef(reg)[1]
	beta = GLM.coef(reg)[2]
	R² = cor(OLSData.Current, OLSData.Pred)^2

	t = string("Calls : y = ", round(beta, digits=4),  "x + ", round(alpha, digits=4), "; R² = ", round(R², digits=4))
	MinMax = [minimum(OLSData.Current), maximum(OLSData.Current)]
	f_MinMax = (beta * MinMax) + repeat([alpha], 2) # [alpha,alpha]
	Calls = scatter(OLSData.Current, OLSData.Pred, xlabel="Current", ylabel="Predicted", title=t, legend=false)
	Calls = plot!(MinMax,f_MinMax, linewidth=2)
	# Puts
	OLSData = DataFrame(Current=vec(Handelstage.Puts_sigma_fit[1:end-1,:]), Pred=vec(Handelstage.Puts_sigma_fit[2:end,:]))
	reg = lm(@formula(Pred ~ Current), OLSData)
	alpha = GLM.coef(reg)[1]
	beta = GLM.coef(reg)[2]
	R² = cor(OLSData.Current, OLSData.Pred)^2

	t = string("Puts : y = ", round(beta, digits=4),  "x + ", round(alpha, digits=4), "; R² = ", round(R², digits=4))
	MinMax = [minimum(OLSData.Current), maximum(OLSData.Current)]
	f_MinMax = (beta * MinMax) + repeat([alpha], 2) # [alpha,alpha]
	Puts = scatter(OLSData.Current, OLSData.Pred, xlabel="Current", ylabel="Predicted", title=t, legend=false)
	Puts = plot!(MinMax,f_MinMax, linewidth=2)
	# Plot
	plot(Calls, Puts, layout=(1,2), size=(1000,500))
	savefig("Scatter_CurrentVsPred.png")
	
# ScatterPlot RMSE/S Current vs SteadyState
	t = string("Calls; Correlation: ", round(cor(Handelstage.Calls_RMSE_S[2:end,:], Handelstage.Calls_RMSE_S_SteadyState[2:end,:])[1], digits=4))
	Calls_RMSE_S_Scatter = scatter(Handelstage.Calls_RMSE_S[2:end,:], Handelstage.Calls_RMSE_S_SteadyState[2:end,:], xlabel="Current RMSE/S", ylabel="Steady State RMSE/S", title=t, legend=false)

	t = string("Puts; Correlation: ", round(cor(Handelstage.Puts_RMSE_S[2:end,:], Handelstage.Puts_RMSE_S_SteadyState[2:end,:])[1], digits=4))
	Puts_RMSE_S_Scatter = scatter(Handelstage.Puts_RMSE_S[2:end,:], Handelstage.Puts_RMSE_S_SteadyState[2:end,:], xlabel="Current RMSE/S", ylabel="Steady State RMSE/S", title=t, legend=false)

	plot(Calls_RMSE_S_Scatter, Puts_RMSE_S_Scatter, layout=(1,2), size=(800,400))
	savefig("ScatterPlot_RMSE_S_Current_vs_SteadyState.png")
