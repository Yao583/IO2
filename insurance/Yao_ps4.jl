using CSV, DataFrames, GLM, Plots, Statistics

data = CSV.read("insdata.csv", DataFrame);
# extract data and assign it to variables
cl = 0.90.*data.m;
ch = 0.96.*data.m;
c  = ch - cl;
data.c = c;

# estimate the demand function D(p) = α + βp
demand = lm(@formula(d ~ 1 + p), data);
# extract the coefficients
coefficients = coef(demand);
α = coefficients[1]
β = coefficients[2]
se1 = stderror(demand)

# estimate the AICh curve AIC_H(p) = γh + δh*p

AICh = lm(@formula(c ~ 1 + p), data[data.d .== 1,:]);
coefficients = coef(AICh);
γh = coefficients[1]
δh = coefficients[2]
se2 = stderror(AICh)
# estimate the AICl curve AIC_L(p) = γl + δl*p
AICl = lm(@formula(c ~ 1 + p), data[data.d .== 0,:]);
coefficients = coef(AICl);
γl = coefficients[1]
δl = coefficients[2]
se3 = stderror(AICl)

# estimate the TC curve
tc = data.d.*ch + (1 .- data.d).*cl;
data.tc = tc;
TC = lm(@formula(tc ~ 1 + p + p^2), data);
coefficients = coef(TC);
ϕ1 = coefficients[1]
ϕ2 = coefficients[2]
ϕ3 = coefficients[3]
se4 = stderror(TC)

# calculate the intercept and slope of MCh, MCl and MCs
μh = γh + α*δh/β
νh = 2*δh
μl = γl - (1-α)*δl/β
νl = 2*δl
μs = ϕ2/β
νs = 2*ϕ3/β
# calculate the corresponding standard errors
seμh = sqrt(se2[1]^2 + se2[2]^2*(α/β)^2 + se1[1]^2*(δh/β)^2 + (α*δh/(β^2))^2*se1[2]^2)
seνh = 2*se2[2]
# plot AIC_H, MC_H as a function of D(p) and the data points
q = 0:1;
p = (q.-α)./β;
AICh = γh .+ δh.*(q.-α)./β;
MCh  = μh .+ νh.*(q.-α)./β;
plot(q, p, label = "Demand")
plot!(q, AICh, label = "AIC_H")
plot!(q, MCh, label = "MC_H")
grouped_data = groupby(data, :p);
sums = combine(grouped_data, :d => sum => :sum_d);
means = combine(grouped_data,:d=>mean=>:mean_d);
sizes = sums.sum_d/200;
scatter!(means.mean_d, means.p, markersize = sizes, label = "Data points")
savefig("plot1.png");
# add MC_L
MCl = μl .+ νl.*(q.-α)./β;
plot!(q, MCl, label = "MC_L")

savefig("plot2.png");
# leave out MC_L but add MC_S
MCs = μs .+ νs.*(q.-α)./β;
plot!(q, MCs, label = "MC_S")
savefig("plot3.png");






# calculating standard errors???
