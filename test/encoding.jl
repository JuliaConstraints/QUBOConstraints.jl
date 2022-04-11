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
    y = [0,1,2,3,4]

    a = [0,0,0,0,0, 1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]

    @test collect(binarize(x; encoding = :domain_wall)) == a
    @test integerize(a; encoding = :domain_wall) == x

    @test integerize(collect(binarize(y, d; encoding = :domain_wall)), d; encoding = :domain_wall) == y
end
