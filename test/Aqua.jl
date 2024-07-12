@testset "Aqua.jl" begin
    import Aqua
    import QUBOConstraints

    # TODO: Fix the broken tests and remove the `broken = true` flag
    Aqua.test_all(
        QUBOConstraints;
        ambiguities = (broken = true,),
        deps_compat = false,
        piracies = (broken = false,)
    )

    @testset "Ambiguities: QUBOConstraints" begin
        # Aqua.test_ambiguities(QUBOConstraints;)
    end

    @testset "Piracies: QUBOConstraints" begin
        Aqua.test_piracies(QUBOConstraints;)
    end

    @testset "Dependencies compatibility (no extras)" begin
        Aqua.test_deps_compat(
            QUBOConstraints;
            check_extras = false,
            ignore = [:LinearAlgebra]
        )
    end
end
