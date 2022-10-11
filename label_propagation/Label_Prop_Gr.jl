#This code reproduces the fit of the label propagation model for granulocyte data
#Date: 10.10.2022; julia v1.6.1
using Plots #v1.23.5
using DataFrames #v1.2.2
using DifferentialEquations #v6.18.0
using Turing #v0.17.4

begin
	#select the number of steps and read data
	nsteps = 10
	timepoint1 = [0.0,3.0,6.0,10.0,11.0]
	timepoint2 = [0.0,3.0,6.0,10.0]
	timepoint3 = [0.0,3.0,6.0,10.0,11.0]
	data1=[0.0,0.72,1.58,8.62,12.2]
	data2=[0.0,0.45,3.1,12.8]
	data3=[0.0,0.0,0.02,2.23,4.44]
end

begin
	#label propagation models for three mice
	#compartment size ratios used where available
	function model1(df,f,p,t)
		n_ratio=t -> (190.0/(1.0+exp(-0.8*(t-3.0)))/44.0)^(1/2)
		df[1] = p[1]*(1/n_ratio(t))*(p[2]-f[1])
		df[2] = p[1]*(1/n_ratio(t))*(f[1]-f[2])
		df[3] = p[1]*(1/22.0)*(f[2]-f[3])
		df[4] = p[1]*(1/2.0)*(f[3]-f[4])
		df[5] = p[1]*(f[4]-f[5])
		df[6] = p[1]*(f[5]-f[6])
		df[7] = p[1]*(f[6]-f[7])
		df[8] = p[1]*(f[7]-f[8])
		df[9] = p[1]*(f[8]-f[9])
		df[10] = p[1]*(f[9]-f[10])
		df[11] = p[1]*(f[10]-f[11])
		df[12] = p[1]*(f[11]-f[12])
	end
	function model2(df,f,p,t)
		n_ratio=t -> (190.0/(1.0+exp(-0.8*(t-3.0)))/44.0)^(1/2)
		df[1] = p[1]*(1/n_ratio(t))*(p[2]-f[1])
		df[2] = p[1]*(1/n_ratio(t))*(f[1]-f[2])
		df[3] = p[1]*(1/22.0)*(f[2]-f[3])
		df[4] = p[1]*(1/2.0)*(f[3]-f[4])
		df[5] = p[1]*(f[4]-f[5])
		df[6] = p[1]*(f[5]-f[6])
		df[7] = p[1]*(f[6]-f[7])
		df[8] = p[1]*(f[7]-f[8])
		df[9] = p[1]*(f[8]-f[9])
		df[10] = p[1]*(f[9]-f[10])
		df[11] = p[1]*(f[10]-f[11])
		df[12] = p[1]*(f[11]-f[12])
	end
	function model3(df,f,p,t)
		df[1] = p[1]*(1/2.8)*(p[2]-f[1])
		df[2] = p[1]*(1/1.6)*(f[1]-f[2])
		df[3] = p[1]*(1/22.0)*(f[2]-f[3])
		df[4] = p[1]*(1/2.0)*(f[3]-f[4])
		df[5] = p[1]*(f[4]-f[5])
		df[6] = p[1]*(f[5]-f[6])
		df[7] = p[1]*(f[6]-f[7])
		df[8] = p[1]*(f[7]-f[8])
		df[9] = p[1]*(f[8]-f[9])
		df[10] = p[1]*(f[9]-f[10])
		df[11] = p[1]*(f[10]-f[11])
		df[12] = p[1]*(f[11]-f[12])
	end
	#specify parameter vector [differentiation rate, initial label frequency in MPPs]
	pars = [0.01,34.0]
	#specify initial condition
	f0s = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	t_span = (.0,15.0)

	prob1 = ODEProblem(model1,f0s,t_span,pars)
	prob2 = ODEProblem(model2,f0s,t_span,pars)
	prob3 = ODEProblem(model3,f0s,t_span,pars)
	sol1 = solve(prob1,saveat=timepoint1)
	sol2 = solve(prob2,saveat=timepoint2)
	sol3 = solve(prob3,saveat=timepoint3)
end

@model function fit_alpha(data1, prob1, data2, prob2, data3, prob3)
	#specify priors for the fitted parameters
	a1 ~ Uniform(0.0,4.0)
	a2 ~ Uniform(0.0,4.0)
	a3 ~ Uniform(0.0,4.0)
	#
	fittedError ~ Uniform(0.1,5.0)
	#
	par1 = [a1, 34.0]
	par2 = [a2, 31.4]
	par3 = [a3, 24.0]
	#
	prob_model1 = remake(prob1,u0=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], p=par1)
	prob_model2 = remake(prob2,u0=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], p=par2)
	prob_model3 = remake(prob3,u0=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], p=par3)
	#
    predicted1 = solve(prob_model1,Tsit5(),saveat=timepoint1)
	predicted2 = solve(prob_model2,Tsit5(),saveat=timepoint2)
	predicted3 = solve(prob_model3,Tsit5(),saveat=timepoint3)
	#compute distance between model and data
	for i = 1:5
		data1[i] ~ Normal(predicted1[nsteps,i], fittedError)
	end
	for i = 1:4
		data2[i] ~ Normal(predicted2[nsteps,i], fittedError)
	end
	for i = 1:5
		data3[i] ~ Normal(predicted3[nsteps,i], fittedError)
	end
end

#fit label propagation model to the data
model = fit_alpha(data1, prob1, data2, prob2, data3, prob3)
#run sampler
nruns = 1000
chain = sample(model,NUTS(0.65), nruns, burnin = 20)

begin
	#simulate model with parameters from the posterior and save results
	nruns = 1000
	sol_mat1=zeros(151,nruns)
	sol_mat2=zeros(151,nruns)
	sol_mat3=zeros(151,nruns)
	parmatrix = DataFrame(chain)
	for i = 1:nruns
		a1 = parmatrix.a1[i]
		a2 = parmatrix.a2[i]
		a3 = parmatrix.a3[i]
		#
		fittedError = median(parmatrix.fittedError)
		#
		solt = collect(0.0:0.1:15.0)
		prob1 = ODEProblem(model1,f0s,t_span,[a1, 34.0])
		prob2 = ODEProblem(model2,f0s,t_span,[a2, 31.4])
		prob3 = ODEProblem(model3,f0s,t_span,[a3, 24.0])
		sol1 = solve(prob1,saveat=solt)
		sol2 = solve(prob2,saveat=solt)
		sol3 = solve(prob3,saveat=solt)
		#
		sol_mat1[:,i] = sol1[nsteps,:]
		sol_mat2[:,i] = sol2[nsteps,:]
		sol_mat3[:,i] = sol3[nsteps,:]
	end
	#select significance level and extract prediction bands
	alevel=0.95
	bands=zeros(151,9)
	for j=1:151
		bands[j,1]=quantile(sol_mat1[j,:],(1-alevel)/2)
		bands[j,2]=median(sol_mat1[j,:])
		bands[j,3]=quantile(sol_mat1[j,:],1-(1-alevel)/2)
		#
		bands[j,4]=quantile(sol_mat2[j,:],(1-alevel)/2)
		bands[j,5]=median(sol_mat2[j,:])
		bands[j,6]=quantile(sol_mat2[j,:],1-(1-alevel)/2)
		#
		bands[j,7]=quantile(sol_mat3[j,:],(1-alevel)/2)
		bands[j,8]=median(sol_mat3[j,:])
		bands[j,9]=quantile(sol_mat3[j,:],1-(1-alevel)/2)
	end
end

begin
	#plot results
	solt = collect(0.0:0.1:15.0)
	plot(solt,bands[:,2],legend=false,ylim=[0.0,15.0],seriescolor=:firebrick,linewidth=2,
		ribbon=(bands[:,2].-bands[:,1],bands[:,3].-bands[:,2]),
		ylabel=string("Label frequency in granulocytes"),xlabel="Time after treatment (days)",
		grid=false,tick_dir=:out)
	plot!(solt,bands[:,5],seriescolor=:firebrick,linewidth=2,
		ribbon=(bands[:,5].-bands[:,4],bands[:,6].-bands[:,5]))
	plot!(solt,bands[:,8],seriescolor=:steelblue,linewidth=2,
		ribbon=(bands[:,8].-bands[:,7],bands[:,9].-bands[:,8]))
	#data
	plot!(timepoint1,data1,seriestype=:scatter,seriescolor=:firebrick,markersize=4.5)
	plot!(timepoint2,data2,seriestype=:scatter,seriescolor=:firebrick,markersize=4.5)
	plot!(timepoint3,data3,seriestype=:scatter,seriescolor=:steelblue,markersize=4.5)
end
