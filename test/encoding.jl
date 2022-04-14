@testset "Encoding: One-Hot" begin
    x = [0,1,2,3,4]
    b = [1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1]

    @test collect(binarize(x)) == b
    @test debinarize(b) == x

    d = domain(x)
    @test debinarize(collect(binarize(x, d)), d) == x
end

@testset "Encoding: Domain Wall" begin
    x = [0,1,2,3,4,5]

    d = domain([0,1,2,3,4])
    y = [0,1,2,3,4]

    a = [0,0,0,0,0, 1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]

    @test collect(binarize(x; binarization = :domain_wall)) == a
    @test debinarize(a; binarization = :domain_wall) == x

    @test debinarize(collect(binarize(y, d; binarization = :domain_wall)), d; binarization = :domain_wall) == y
end
