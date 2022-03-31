@testset "Encoding: One-Hot" begin
    x = [0,1,2,3,4]
    b = [1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1]

    @test collect(binarize(x)) == b
    @test integerize(b) == x
end

@testset "Encoding: Domain Wall" begin
    x = [0,1,2,3,4,5]
    b = [0,0,0,0,0, 1,0,0,0,0, 1,1,0,0,0, 1,1,1,0,0, 1,1,1,1,0, 1,1,1,1,1]

    @test collect(binarize(x; encoding = :domain_wall)) == b
    @test integerize(b; encoding = :domain_wall) == x
end