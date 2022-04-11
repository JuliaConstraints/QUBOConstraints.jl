@testset "Encoding: One-Hot" begin
    x = [0,1,2,3,4]
    b = [1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1]

    @test collect(binarize(x)) == b
    @test integerize(b) == x

    d = domain(x)
    @test integerize(collect(binarize(x, d)), d) == x
end

@testset "Encoding: Domain Wall" begin
    x = [0,1,2,3,4,5]

    d = domain([0,1,2,3,4])
    y = get_domain(d)

    a = [0,0,0,0,0, 1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]
    b = [1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]

    @test collect(binarize(x; encoding = :domain_wall)) == a
    @test integerize(a; encoding = :domain_wall) == x

    # @info collect(binarize(y, d; encoding = :domain_wall))
    @info d integerize(b, d; encoding = :domain_wall)

    # @test collect(binarize(y, domain([0,1,2,3,4]); encoding = :domain_wall)) == b
    # @test integerize(b, domain([0,1,2,3,4]); encoding = :domain_wall) == y
end
