

lines=split(read(ARGS[1], String),"\n")[2:end-1]
lines=map(l->split(l,"\t"),lines)
pressures=[parse(Float64,l[2]) for l in lines]
rhos=[parse(Float64,l[3]) for l in lines]
pref=split(ARGS[1],".")[1]
f=open(pref*".dat","a")
for (r,p) in zip(rhos,pressures)
println(f,r," ",p)
end
println(f)
print(f,"T=",lines[1][1],"K ,format density(kg*m^-3) pressure(MPa)")
close(f)
