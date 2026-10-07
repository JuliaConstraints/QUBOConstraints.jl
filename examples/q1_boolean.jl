using QUBOConstraints
x = codebook(:x, 0:1; encoding=:domain_wall)
y = codebook(:y, 0:1; encoding=:domain_wall)
q = QUBOComponent(vcat(x.bits,y.bits); offset=1,
    linear=[(1,-1),(2,-1)],quadratic=[(1,2,2)],codebooks=[x,y],
    provenance="Boolean allDifferent equality indicator")
report = exhaustive_check(q, values -> values[1] != values[2];
    oracle_id="allDifferent/boolean/v1",profile=:indicator_exact)
@assert report.status === :pass && report.proof === :exhaustive
@assert report.evaluations == 4
if isempty(ARGS)
    write_component(stdout,q;report)
else
    open(ARGS[1],"w") do io
        write_component(io,q;report)
    end
end
